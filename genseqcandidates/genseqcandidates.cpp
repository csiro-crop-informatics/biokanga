// genseqcandidates.cpp : Defines the entry point for the console application.
// generate sequencing candidates for reads of interest
// Primarily directed towards nucleosome detection in C&M
// Input is a csv file containing sequence loci of interest
// These loci are then spacially with the objective function of maximising the
// number of sequence loci within a user specified block size - defaults to 1Knt
// Contained subsequences of a user specified length - defaults to 25 - are then
// aligned against the whole genome, if these contained subsequences are non-unique then
// the block will be penalised with additional penalty if the subsequence was contained within
// an original sequence loci of interest.
// The generated results file contains a list of blocks, sequences of interest contained within each block plus block scores
// Note that sequences of interest are not unique to a single block, multiple blocks may contain the same sequence.
#include "stdafx.h"

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#if _WIN32
#include "../libbiokanga/commhdrs.h"
#else
#include "../libbiokanga/commhdrs.h"
#endif
#include "IncExclChroms.h"

const unsigned int cProgVer = 109;		// increment with each release

//const int cMaxIncludeChroms = 20;		// max number of include chromosomes regular expressions
//const int cMaxExcludeChroms = 20;		// max number of exclude chromosomes regular expressions

const int cMinMinRegionLen = 10;
const int cDfltMinRegionLen = 147;
const int cMaxMinRegionLen = 300;

const int cMinDeltaLen = -2000;
const int cDfltDeltaLen = 0;
const int cMaxDeltaLen = 2000;

const int cMinOfsLoci = -2000;
const int cDfltOfsLoci = 0;
const int cMaxOfsLoci = 2000;

const int cMinTruncLen = 10;
const int cDfltTruncLen = 147;
const int cMaxTruncLen = 300;



const int cMinSeqBlockLen = 300;	// minimum block sequence length
const int cDfltSeqBlockLen = 1000;	// default block sequence length
const int cMaxSeqBlockLen = 2000;	// blocks can hold upto this length sequence

const int cMinSubSeqLen = 15;
const int cDfltSubSeqLen = 25;
const int cMaxSubSeqLen = 100;



// processing modes
typedef enum TAG_ePMode {		
	ePMdefault,					// default processing mode
	ePMplaceholder				// used to set the enumeration range
	} etPMode;

// output format modes
typedef enum TAG_eFMode {		
	eFMdefault,					// default is for csv
	eFMplaceholder				// used to set the enumeration range
	} etFMode;

#pragma pack(1)
typedef struct TAG_sSeqBlock {
	int BlockID;			// uniquely identifies this block
	int Score;				// score associated with this block
	int NumRegions;			// number of regions mapping to this block
	int ChromID;			// block is on this chromosome
	int StartLoci;			// starts at this loci
	int EndLoci;			// ends at this loci
	int SeqLen;				// and is of this length
	UINT8 Seq[cMaxSeqBlockLen];		// block sequence, low order bits contain bases, bit 5 marks if this base is from a sequence of interest, bit 6 marks if containing window is none-unique
} tsSeqBlock;


typedef struct TAG_sRegion {
	int RegionID;					// uniquely identifies this region
	int SrcID;						// original region source identifier
	int BlockID;					// region is contained within this sequence block
	int NumInstances;				// number of instances of this region
	int ChromID;					// region is on this chromosome
	int StartLoci;					// loci at which region starts
	int EndLoci;					// loci at which region ends
	char Strand;					// is on '+' or '-' strand
} tsRegion;
#pragma pack()

int
Process(etPMode PMode,				// processing mode
		int BlockSeqLen,			// block sequences to be maximally of this length
		int SubSeqLen,				// with subsequence matches of this length
		int OfsLoci,				// offset region start loci by this many nt
		int DeltaLen,				// change region length by this many nt
		int TruncLength,			// truncate regions to be a maximum of this length
		int MinLength,				// regions to be of this minimum length	
		char *pszSfxFile,			// file containing genome suffix array
		char *pszRegionsFile,		// file containing region locii of interest
		char *pszOutFile,			// results to this file
		int	NumIncludeChroms,			// number of chromosome regular expressions to include
		char **ppszIncludeChroms,		// array of include chromosome regular expressions
		int	NumExcludeChroms,			// number of chromosome expressions to exclude
		char **ppszExcludeChroms);		// array of exclude chromosome regular expressions

int LocateBlockMatches(int SubSeqLen);	// subsequences are of this length

int LoadRegions(char *pszRegionFile,	// csv file containing region loci
		int OfsLoci,				// offset region start loci by this many nt
		int DeltaLen,				// change region length by this many nt
		int TruncLength,			// truncate regions to be a maximum of this length
		int MinLength,				// regions to be of this minimum length	
		int	NumIncludeChroms,			// number of chromosome regular expressions to include
		char **ppszIncludeChroms,		// array of include chromosome regular expressions
		int	NumExcludeChroms,			// number of chromosome expressions to exclude
		char **ppszExcludeChroms);		// array of exclude chromosome regular expressions

int AddBlock(int ChromID,			// block is on this chromosome
			 int StartLoci,			// with this start loci
			 int EndLoci);			// ending at this loci

void Init(void);
void Reset(void);
int SortRegions( const void *arg1, const void *arg2);

CStopWatch gStopWatch;
CDiagnostics gDiagnostics;				// for writing diagnostics messages to log file
char gszProcName[_MAX_FNAME];			// process name
bool gbActivity;						// used to determine if activity messages vi printf's can be used - output activity if eDLInfo or less

#ifdef _WIN32
// required by str library
#if !defined(__AFX_H__)  ||  defined(STR_NO_WINSTUFF)
HANDLE STR_get_stringres()
{
	return NULL;	//Works for EXEs; in a DLL, return the instance handle
}
#endif

const STRCHAR* STR_get_debugname()
{
	return _T("alignreads");
}
// end of str library required code
#endif

#ifdef _WIN32
int _tmain(int argc, char* argv[])
{
// determine my process name
_splitpath(argv[0],NULL,NULL,gszProcName,NULL);
#else
int 
main(int argc, const char** argv)
{
// determine my process name
CUtility::splitpath((char *)argv[0],NULL,gszProcName);
#endif
int iScreenLogLevel;		// level of screen diagnostics
int iFileLogLevel;			// level of file diagnostics
char szLogFile[_MAX_PATH];	// write diagnostics to this file
int Rslt = 0;   			// function result code >= 0 represents success, < 0 on failure
int Idx;
int LenChromList;

etPMode PMode;				// processing mode
int MinLength;				// regions must be of at least this length
int OfsLoci;				// offset region loci by this many bases
int DeltaLen;				// change region lengths by this many bases
int TruncLength;			// and truncate to be no longer than this length

int BlockSeqLen;			// blocks to be of this length
int SubSeqLen;				// align block subsequences of this length
char szInFile[_MAX_PATH];	// input element loci from this file
char szInSfxFile[_MAX_PATH];	// input bioseq file containing suffix array
char szRsltsFile[_MAX_PATH];	// output block loci to this file

int NumIncludeChroms;
char *pszIncludeChroms[cMaxIncludeChroms];
int NumExcludeChroms;
char *pszExcludeChroms[cMaxExcludeChroms];

// command line args
struct arg_lit  *help    = arg_lit0("hH","help",                "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_int *ScreenLogLevel=arg_int0("S", "ScreenLogLevel",	"<int>","Level of diagnostics written to screen 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"diagnostics log file");

struct arg_int *pmode = arg_int0("m","mode","<int>",		    "processing mode: 0 - generate blocks (default = 0)");
struct arg_int *subseqlen = arg_int0("s","subseqlen","<int>",	"subsequence alignment length (default: 25)");
struct arg_int *blockseqlen = arg_int0("b","blockseqlen","<int>", "generated block lenth (default: 1000)");
struct arg_int  *minlength = arg_int0("l","minlength","<int>",	"minimum element length (default 147)");
struct arg_int  *trunclength = arg_int0("T","truncatelength","<int>","truncate regions to be no  longer than this length (default 147)");
struct arg_int  *ofsloci   = arg_int0("u","offset","<int>",	    "offset element loci by this many bases, -2048 to +2048 (default 0)");
struct arg_int  *deltalen  = arg_int0("U","deltalen","<int>",	"delta element lengths by this many bases, -2048 to +2048 (default 0)");

struct arg_file *infile = arg_file1("i","in","<file>",			"input regions of interest from this csv loci file");
struct arg_file *sfxfile = arg_file1("I","sfx","<file>",		"align generated block subsequences against this suffix array file");
struct arg_file *outfile = arg_file1("o","out","<file>",		"output blocks containing regions to this file");
struct arg_str  *ExcludeChroms = arg_strn("Z","chromexclude","<string>",0,cMaxExcludeChroms,"high priority - regular expressions defining chromosomes to exclude");
struct arg_str  *IncludeChroms = arg_strn("z","chromeinclude","<string>",0,cMaxIncludeChroms,"low priority - regular expressions defining chromosomes to include");

struct arg_end *end = arg_end(20);

void *argtable[] = {help,version,FileLogLevel,ScreenLogLevel,LogFile,
					pmode,subseqlen,blockseqlen,minlength,trunclength,ofsloci,deltalen,infile,sfxfile,outfile,ExcludeChroms,IncludeChroms,
					end};

char **pAllArgs;
int argerrors;
argerrors = CUtility::arg_parsefromfile(argc,(char **)argv,&pAllArgs);
if(argerrors >= 0)
	argerrors = arg_parse(argerrors,pAllArgs,argtable);

/* special case: '--help' takes precedence over error reporting */
if (help->count > 0)
        {
		printf("\n%s ", gszProcName);
        arg_print_syntax(stdout,argtable,"\n");
        arg_print_glossary(stdout,argtable,"  %-25s %s\n");
		printf("\nNote: Parameters can be entered into a parameter file, one parameter per line.");
		printf("\n      To invoke this parameter file then precede it's name with '@'");
		printf("\n      e.g. %s @myparams.txt\n\n",gszProcName);
		exit(1);
        }

    /* special case: '--version' takes precedence error reporting */
if (version->count > 0)
        {
		printf("\n%s Version %d.%2.2d",gszProcName,cProgVer/100,cProgVer%100);
		exit(1);
        }

if (!argerrors)
	{
	iScreenLogLevel = ScreenLogLevel->count ? ScreenLogLevel->ival[0] : eDLInfo;
	if(iScreenLogLevel < eDLNone || iScreenLogLevel > eDLDebug)
		{
		printf("\nError: ScreenLogLevel '-S%d' specified outside of range %d..%d",iScreenLogLevel,eDLNone,eDLDebug);
		exit(1);
		}
	gbActivity = iScreenLogLevel >= eDLInfo ? true : false;
	if(FileLogLevel->count && !LogFile->count)
		{
		printf("\nError: FileLogLevel '-f%d' specified but no logfile '-F<logfile>'",FileLogLevel->ival[0]);
		exit(1);
		}

	iFileLogLevel = FileLogLevel->count ? FileLogLevel->ival[0] : eDLInfo;
	if(iFileLogLevel < eDLNone || iFileLogLevel > eDLDebug)
		{
		printf("\nError: FileLogLevel '-l%d' specified outside of range %d..%d",iFileLogLevel,eDLNone,eDLDebug);
		exit(1);
		}
	if(LogFile->count)
		{
		strncpy(szLogFile,LogFile->filename[0],_MAX_PATH);
		szLogFile[_MAX_PATH-1] = '\0';
		}
	else
		{
		iFileLogLevel = eDLNone;
		szLogFile[0] = '\0';
		}

	PMode = (etPMode)(pmode->count ? pmode->ival[0] : 0);
	if(PMode < 0 || PMode >= ePMplaceholder)
		{
		printf("\nError: Processing mode '-m%d' specified outside of range %d..%d",PMode,0,(int)ePMplaceholder-1);
		exit(1);
		}

	MinLength = minlength->count ? minlength->ival[0] : cDfltMinRegionLen;
	if(MinLength < cMinMinRegionLen || MinLength > cMaxMinRegionLen)
		{
		printf("\nError: minimum length '-l%d' specified outside of range %d..%d",MinLength,cMinMinRegionLen,cMaxMinRegionLen);
		exit(1);
		}

	OfsLoci = ofsloci->count ? ofsloci->ival[0] : cDfltOfsLoci;
	if(OfsLoci < cMinOfsLoci || OfsLoci > cMaxOfsLoci)
		{
		printf("\nError: Offset region loci '-u%d' specified outside of range %d..%d",OfsLoci,cMinOfsLoci,cMaxOfsLoci);
		exit(1);
		}


	DeltaLen = deltalen->count ? deltalen->ival[0] : cDfltDeltaLen;
	if(DeltaLen < cMinDeltaLen || DeltaLen > cMaxDeltaLen)
		{
		printf("\nError: Delta region length '-U%d' specified outside of range %d..%d",DeltaLen,cMinDeltaLen,cMaxDeltaLen);
		exit(1);
		}

	TruncLength = trunclength->count ? trunclength->ival[0] : cDfltTruncLen;
	if(TruncLength < cMinTruncLen || TruncLength > cMaxTruncLen)
		{
		printf("\nError: Truncate region length '-T%d' specified outside of range %d..%d",TruncLength,cMinTruncLen,cMaxTruncLen);
		exit(1);
		}

	BlockSeqLen = blockseqlen->count ? blockseqlen->ival[0] : cDfltSeqBlockLen;
	if(BlockSeqLen < cMinSeqBlockLen || BlockSeqLen > cMaxSeqBlockLen)
		{
		printf("\nError: Block sequence length '-b%d' specified outside of range %d..%d",BlockSeqLen,cMinSeqBlockLen,cMaxSeqBlockLen);
		exit(1);
		}

	SubSeqLen = subseqlen->count ? subseqlen->ival[0] : cDfltSubSeqLen;
	if(SubSeqLen < cMinSubSeqLen || SubSeqLen > cMaxMinRegionLen)
		{
		printf("\nError: Block subsequence align length '-s%d' specified outside of range %d..%d",SubSeqLen,cMinSubSeqLen,cMaxSubSeqLen);
		exit(1);
		}

	NumIncludeChroms = IncludeChroms->count;
	for(Idx=0;Idx < IncludeChroms->count; Idx++)
		{
		LenChromList = (int)strlen(IncludeChroms->sval[Idx]);
		pszIncludeChroms[Idx] = new char [LenChromList+1];
		strcpy(pszIncludeChroms[Idx],IncludeChroms->sval[Idx]);
		CUtility::TrimQuotes(pszIncludeChroms[Idx]);
		}

	NumExcludeChroms = ExcludeChroms->count;
	for(Idx=0;Idx < ExcludeChroms->count; Idx++)
		{
		LenChromList = (int)strlen(ExcludeChroms->sval[Idx]);
		pszExcludeChroms[Idx] = new char [LenChromList+1];
		strcpy(pszExcludeChroms[Idx],ExcludeChroms->sval[Idx]);
		CUtility::TrimQuotes(pszExcludeChroms[Idx]);
		}

	strcpy(szInFile,infile->filename[0]);
	strcpy(szInSfxFile,sfxfile->filename[0]);
	strcpy(szRsltsFile,outfile->filename[0]);

			// now that command parameters have been parsed then initialise diagnostics log system
	if(!gDiagnostics.Open(szLogFile,(etDiagLevel)iScreenLogLevel,(etDiagLevel)iFileLogLevel,true))
		{
		printf("\nError: Unable to start diagnostics subsystem.");
		if(szLogFile[0] != '\0')
			printf(" Most likely cause is that logfile '%s' can't be opened/created",szLogFile);
		exit(1);
		}

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Version: %d.%2.2d Processing parameters:",cProgVer/100,cProgVer%100);
	const char *pszDescr;
	switch(PMode) {
		case ePMdefault:				
			pszDescr = "Generate blocks containing regions";
			break;
		}
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"processing mode is : '%s'",pszDescr);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"generate blocks of length : %d",BlockSeqLen);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"align block subsequences of length : %d",SubSeqLen);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"offset region loci by: %d",OfsLoci);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"change region length by: %d",DeltaLen);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"truncate regions to be no longer than: %d",TruncLength);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"slough regions of less than: %d",MinLength);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"genome suffix array file : '%s'",szInSfxFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"region loci file : '%s'",szInFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"write block loci to : '%s'",szRsltsFile);

	for(Idx = 0; Idx < NumIncludeChroms; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"reg expressions defining chroms to include: '%s'",pszIncludeChroms[Idx]);
	for(Idx = 0; Idx < NumExcludeChroms; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"reg expressions defining chroms to exclude: '%s'",pszExcludeChroms[Idx]); 


#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	gStopWatch.Start();
	Rslt = Process(PMode,BlockSeqLen,SubSeqLen,OfsLoci,DeltaLen,TruncLength,MinLength,szInSfxFile,szInFile,szRsltsFile,NumIncludeChroms,pszIncludeChroms,NumExcludeChroms,pszExcludeChroms);
	gStopWatch.Stop();
	Rslt = Rslt >=0 ? 0 : 1;
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Exit code: %d Total processing time: %s",Rslt,gStopWatch.Read());
	exit(Rslt);
	}
else
	{
	arg_print_errors(stdout,end,gszProcName);
	arg_print_syntax(stdout,argtable,"\nend of help\n");
	exit(1);
	}
return 0;
}


CSfxArrayV3 *m_pSfxArray;				// suffix array over genome of interest
etPMode m_PMode;					// processing mode
etFMode m_FMode;					// output format mode

int m_hOutFile;						// file handle for output results
char *m_pszOutBuff;					// used for buffering output results
int m_AllocdOutBuff;				// how many chars have been allocated to m_pszOutBuff
int m_UsedOutBuff;					// how many chars of results currently in m_pszOutBuff

CCSVFile *m_pCSV;					// used whilst loading regions from CSV file
const int cRegionAlloc = 100000;	// allocate for this many tsRegion
tsRegion *m_pRegions;				// pts to loaded regions of interest
int m_AllocdRegions;				// memory allocated for this many regions
int m_NumRegions;					// this many regions have been loaded


const int cBlocksAlloc = 10000;		// allocation for this number of tsSeqBlock's
tsSeqBlock *m_pSeqBlocks;			// pts to memory allocated for holding linked list of sequence blocks
int m_AllocdBlocks;					// allocated for this many blocks
int m_NumBlocks;					// this many blocks used

int
Process(etPMode PMode,				// processing mode
		int BlockSeqLen,			// block sequences to be maximally of this length
		int SubSeqLen,				// with subsequence matches of this length
		int OfsLoci,				// offset region start loci by this many nt
		int DeltaLen,				// change region length by this many nt
		int TruncLength,			// truncate regions to be a maximum of this length
		int MinLength,				// regions to be of this minimum length	
		char *pszSfxFile,			// file containing genome suffix array
		char *pszRegionsFile,		// file containing region locii of interest
		char *pszOutFile,			// results to this file
		int	NumIncludeChroms,			// number of chromosome regular expressions to include
		char **ppszIncludeChroms,		// array of include chromosome regular expressions
		int	NumExcludeChroms,			// number of chromosome expressions to exclude
		char **ppszExcludeChroms)		// array of exclude chromosome regular expressions
{
int Rslt;
Init();

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Loading suffix array file '%s'", pszSfxFile);
if((m_pSfxArray = new CSfxArrayV3()) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CSfxArray");
	Reset();
	return(eBSFerrObj);
	}

if((Rslt=m_pSfxArray->Open(pszSfxFile,false))!=eBSFSuccess)
	{
	while(m_pSfxArray->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,m_pSfxArray->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open input bioseq suffix array file '%s'",pszSfxFile);
	Reset();
	return(Rslt);
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Loading regions of interest from file '%s'", pszRegionsFile);
if((Rslt = LoadRegions(pszRegionsFile,OfsLoci,DeltaLen,TruncLength,MinLength,NumIncludeChroms,ppszIncludeChroms,NumExcludeChroms,ppszExcludeChroms)) <= 0)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to load regions file '%s'",pszRegionsFile);
	Reset();
	return(Rslt);
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing %d regions into contigous blocks of size %d...",m_NumRegions,BlockSeqLen);

tsSeqBlock *pBlock;
UINT8 *pSeq;
tsRegion *pRegion;
tsRegion *pRegion1;
int RegionIdx;
int RegionIdx2;
int StartBlockRegionIdx;
int EndBlockRegionIdx;
int ChromID;
int ChromLen;
int BlockID;
int StartLoci;
int EndLoci;
int DLoci;
int SeqLen;
int Score;
int SeqIdx;
int Matches;
bool bNonCannonical;
bool bInRegion;
int LowScoreBlocks;

pRegion = m_pRegions;
ChromID = pRegion->ChromID;
ChromLen = m_pSfxArray->GetSeqLen(ChromID);
StartLoci = pRegion->StartLoci;
StartBlockRegionIdx = 0;
EndBlockRegionIdx = 0;
for(RegionIdx = 0; RegionIdx < m_NumRegions; RegionIdx++,pRegion++)
	{
	if((RegionIdx + 1) == m_NumRegions || // if no more regions or
		pRegion->ChromID != ChromID ||    // no longer on same chromosome or
		pRegion->EndLoci >= (StartLoci + BlockSeqLen)) // including this region would make block too large
		{
LastRegion:			// hate goto's but need to handle last region if region needs it's own separate block
		EndLoci = StartLoci + BlockSeqLen - 1;
		// centralise regions within the block
		int leftloci = m_pRegions[StartBlockRegionIdx].StartLoci;
		int rightloci = m_pRegions[EndBlockRegionIdx].EndLoci;
		DLoci = ((m_pRegions[StartBlockRegionIdx].StartLoci - StartLoci) +
				(EndLoci - m_pRegions[EndBlockRegionIdx].EndLoci))/2;
		StartLoci = m_pRegions[StartBlockRegionIdx].StartLoci - DLoci;
		if(StartLoci < 0)
			StartLoci = 0;
		EndLoci = StartLoci + BlockSeqLen - 1;
		if(EndLoci >= ChromLen)
			{
			EndLoci = ChromLen - 1;
			StartLoci = 1 + EndLoci - BlockSeqLen;
			}
		if((BlockID = AddBlock(ChromID,StartLoci,EndLoci)) < 1)
			{
			Reset();
			return(BlockID);
			}
		pBlock = &m_pSeqBlocks[BlockID-1];
		pBlock->NumRegions = 0;
		pRegion1 = &m_pRegions[StartBlockRegionIdx];
		for(RegionIdx2 = StartBlockRegionIdx; RegionIdx2 <= EndBlockRegionIdx; RegionIdx2++,pRegion1++)
			{
			pRegion1->BlockID = BlockID;
			pBlock->NumRegions += 1;
			// mark those block bases which are from the regions of interest
			pSeq = &pBlock->Seq[pRegion1->StartLoci - StartLoci];
			SeqLen = 1 + pRegion1->EndLoci - pRegion1->StartLoci;
			for(SeqIdx = 0; SeqIdx < SeqLen; SeqIdx++,pSeq++)
				*pSeq |= 0x080;
			}

		if(pRegion->BlockID > 0 && (RegionIdx + 1) == m_NumRegions)
			continue;
		if(pRegion->ChromID != ChromID)
			{
			ChromID = pRegion->ChromID;
			ChromLen = m_pSfxArray->GetSeqLen(ChromID);
			}
		StartBlockRegionIdx = RegionIdx;
		StartLoci = pRegion->StartLoci;
		if((RegionIdx + 1) == m_NumRegions) // hate goto's but need to handle last region if region needs it's own separate block
			{
			EndBlockRegionIdx = RegionIdx;
			goto LastRegion;
			}
		}
	EndBlockRegionIdx = RegionIdx;
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Generated %d sequence blocks from %d regions, processing for aligned subsequences of length %d uniqueness...",m_NumBlocks,m_NumRegions,SubSeqLen);
Rslt = LocateBlockMatches(SubSeqLen);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Aligned all subsequences, now scoring...");

// iterate over blocks and associate a score
// aggregate score by -
// a) score non-unique subsequences x 1 unless subsequence part of a region in which case score by x 3 
// b) any non-acgt's then score as though that subsequence had multiple matches, x SubSeqLen
LowScoreBlocks = 0;
pBlock = m_pSeqBlocks;
for(BlockID = 0; BlockID < m_NumBlocks; BlockID++, pBlock++)
	{
	Score = 0;
	pSeq = pBlock->Seq;
	for(SeqIdx = 0; SeqIdx < pBlock->SeqLen; SeqIdx++,pSeq++)
		{
		Matches = (*pSeq & 0x060) >> 5;
		bNonCannonical = *pSeq & 0x010 ? 1 : 0;
		bInRegion = *pSeq & 0x080 ? true : false;
		if(Matches > 1)
			Score += bInRegion ? 3 : 1;
		if(bNonCannonical)
			Score += bInRegion ? SubSeqLen * 2 : SubSeqLen;
		*pSeq &= 0x07;
		}
	if(Score <= SubSeqLen)
		LowScoreBlocks++;
	pBlock->Score = Score;
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Completed scoring, from %d blocks there are %d blocks with scores of %d or less, now generating output file...",m_NumBlocks,LowScoreBlocks,SubSeqLen);

char szChrom[cMaxDatasetSpeciesChrom];
char szBuff[cMaxSeqBlockLen * 4];
int BuffIdx;
char szSpecies[cMaxDatasetSpeciesChrom];

ChromID = -1;
strcpy(szSpecies,m_pSfxArray->GetDatasetName());

// output as the usual csv loci file -
#ifdef _WIN32
	if((m_hOutFile = open(pszOutFile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
	if((m_hOutFile = open(pszOutFile, O_RDWR | O_CREAT | O_TRUNC, S_IREAD | S_IWRITE ))==-1)
#endif
		{
		gDiagnostics.DiagOut(eDLFatal,"Process","Unable to create or truncate results file %s error: %s",pszOutFile,strerror(errno));
		Reset();
		return(eBSFerrCreateFile);
		}

BuffIdx = 0;
pBlock = m_pSeqBlocks;
for(BlockID = 0; BlockID < m_NumBlocks; BlockID++, pBlock++)
	{
	if(ChromID != pBlock->ChromID)
		{
		m_pSfxArray->GetIdentName(pBlock->ChromID,sizeof(szChrom),szChrom);
		ChromID = pBlock->ChromID;
		}
	BuffIdx += sprintf(&szBuff[BuffIdx],"%d,\"bk\",\"%s\",\"%s\",%d,%d,%d,\"+\",%d,",
			pBlock->BlockID,szSpecies,szChrom,pBlock->StartLoci,pBlock->EndLoci,pBlock->SeqLen,pBlock->Score);
	BuffIdx+=sprintf(&szBuff[BuffIdx],"\"%s\"\n",CSeqTrans::MapSeq2Ascii(pBlock->Seq,pBlock->SeqLen));
	if((BuffIdx + 200 + BlockSeqLen) > sizeof(szBuff))
		{
		if(write(m_hOutFile,szBuff,BuffIdx)!=BuffIdx)
			{
			gDiagnostics.DiagOut(eDLFatal,"Process","Unable to write results to file %s - error %s",pszOutFile,strerror(errno));
			Reset();
			return(eBSFerrFileAccess);
			}
		BuffIdx = 0;
		}
	}
if(BuffIdx && write(m_hOutFile,szBuff,BuffIdx)!=BuffIdx)
	{
	gDiagnostics.DiagOut(eDLFatal,"Process","Unable to write results to file %s - error %s",pszOutFile,strerror(errno));
	Reset();
	return(eBSFerrFileAccess);
	}

BuffIdx = sprintf(szBuff,"\n\n");
// after the blocks, write out regions contained within each block
BuffIdx = 0;
ChromID = -1;
pRegion = m_pRegions;
for(RegionIdx = 0; RegionIdx < m_NumRegions; RegionIdx++, pRegion++)
	{
	if(ChromID != pRegion->ChromID)
		{
		m_pSfxArray->GetIdentName(pRegion->ChromID,sizeof(szChrom),szChrom);
		ChromID = pRegion->ChromID;
		}

	BuffIdx += sprintf(&szBuff[BuffIdx],"%d,\"region\",\"Block_%d\",\"%s\",%d,%d,%d,\"%c\",%d\n",
		pRegion->SrcID,pRegion->BlockID,szChrom,pRegion->StartLoci,pRegion->EndLoci,1 + pRegion->EndLoci-pRegion->StartLoci,pRegion->Strand,pRegion->NumInstances);
	if((BuffIdx + 200) > sizeof(szBuff))
		{
		if(write(m_hOutFile,szBuff,BuffIdx)!=BuffIdx)
			{
			gDiagnostics.DiagOut(eDLFatal,"Process","Unable to write results to file %s - error %s",pszOutFile,strerror(errno));
			Reset();
			return(eBSFerrFileAccess);
			}
		BuffIdx = 0;
		}
	}
if(BuffIdx && write(m_hOutFile,szBuff,BuffIdx)!=BuffIdx)
	{
	gDiagnostics.DiagOut(eDLFatal,"Process","Unable to write results to file %s - error %s",pszOutFile,strerror(errno));
	Reset();
	return(eBSFerrFileAccess);
	}

Reset();
return(Rslt);
}

int
AddBlock(int ChromID,int StartLoci,int EndLoci)
{
int Rslt;
tsSeqBlock *pSeqBlock;
int SeqLen = 1 + EndLoci - StartLoci;

if(m_pSeqBlocks == NULL || m_NumBlocks == m_AllocdBlocks)
	{
	if(m_pSeqBlocks == NULL)
		{
		if((pSeqBlock = (tsSeqBlock *)malloc(sizeof(tsSeqBlock) * cBlocksAlloc))==NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to malloc memory (%d bytes requested) for region filtering",sizeof(tsSeqBlock) * cBlocksAlloc);
			Reset();
			return(eBSFerrMem);
			}
		m_AllocdBlocks = 0;
		m_NumBlocks = 0;
		}
	else
		{
		if((pSeqBlock = (tsSeqBlock *)realloc(m_pSeqBlocks,sizeof(tsSeqBlock) * (m_AllocdBlocks+cBlocksAlloc)))==NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to realloc memory (%d bytes requested) for region filtering",sizeof(tsSeqBlock) * (m_AllocdBlocks+cBlocksAlloc));
			Reset();
			return(eBSFerrMem);
			}
		}
	m_pSeqBlocks = pSeqBlock;
	m_AllocdBlocks += cBlocksAlloc;
	}
pSeqBlock = &m_pSeqBlocks[m_NumBlocks++];
pSeqBlock->BlockID = m_NumBlocks;
pSeqBlock->ChromID = ChromID;
pSeqBlock->StartLoci= StartLoci;
pSeqBlock->EndLoci = EndLoci;
pSeqBlock->SeqLen = 1 + EndLoci - StartLoci;

if((Rslt=(int)m_pSfxArray->GetSeq(ChromID,StartLoci,pSeqBlock->Seq,pSeqBlock->SeqLen)) != pSeqBlock->SeqLen)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to get expected sequence from suffix array (ChromID: %d StartLoci: %d SeqLen: %d)!",ChromID,StartLoci,pSeqBlock->SeqLen);
	Reset();
	return(eBSFerrInternal);
	}

return(m_NumBlocks);
}

void
Init(void)
{
m_pCSV = NULL;
m_pSfxArray = NULL;
m_pszOutBuff = NULL;
m_pRegions = NULL;
m_pSeqBlocks = NULL;
m_hOutFile = -1;
Reset();
}

void
Reset(void)
{
if(m_hOutFile != -1)
	{
	close(m_hOutFile);
	m_hOutFile = -1;
	}
if(m_pCSV != NULL)
	{
	delete m_pCSV;
	m_pCSV = NULL;
	}
if(m_pszOutBuff != NULL)
	{
	delete m_pszOutBuff;
	m_pszOutBuff = NULL;
	}
if(m_pSfxArray != NULL)
	{
	delete m_pSfxArray;
	m_pSfxArray = NULL;
	}
if(m_pSeqBlocks != NULL)
	{
	delete m_pSeqBlocks;
	m_pSeqBlocks = NULL;
	}
if(m_pRegions != NULL)
	{
	free(m_pRegions);				// NOTE - use free() as memory malloc'd/realloc'd
	m_pRegions = NULL;
	}
m_AllocdOutBuff = 0;				// how many chars have been allocated to m_pszOutBuff
m_UsedOutBuff = 0;					// how many chars of results currently in m_pszOutBuff
m_AllocdRegions = 0;				// memory allocated for this many regions
m_NumRegions = 0;					// this many regions have been loaded
m_AllocdBlocks = 0;					// memory allocated for this many blocks
m_NumBlocks = 0;					// this many blocks
}

int
AddRegion(int SrcID,				// region source identifier
		  int Instances,			// number of instances to associate with this region
		  char *pszChrom,			// is on this chrom
		  int StartLoci,			// starting at this loci
		  int EndLoci,				// and ending at this loci
		  bool bOnMinStrand)		// true if on minus strand
{
int ChromID;
tsRegion *pRegion;
if((ChromID = m_pSfxArray->GetIdent(pszChrom)) < 1)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to locate region chromosome '%s'",pszChrom);
	return(eBSFerrChrom);
	}

if(m_pRegions == NULL || m_NumRegions == m_AllocdRegions)
	{
	if(m_pRegions == NULL)
		{
		if((pRegion = (tsRegion *)malloc(sizeof(tsRegion) * cRegionAlloc))==NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to malloc memory (%d bytes requested) for region filtering",sizeof(tsRegion) * cRegionAlloc);
			Reset();
			return(eBSFerrMem);
			}
		m_AllocdRegions = 0;
		m_NumRegions = 0;
		}
	else
		{
		if((pRegion = (tsRegion *)realloc(m_pRegions,sizeof(tsRegion) * (m_AllocdRegions+cRegionAlloc)))==NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to realloc memory (%d bytes requested) for region filtering",sizeof(tsRegion) * (m_AllocdRegions+cRegionAlloc));
			Reset();
			return(eBSFerrMem);
			}
		}
	m_pRegions = pRegion;
	m_AllocdRegions += cRegionAlloc;
	}
pRegion = &m_pRegions[m_NumRegions++];
pRegion->RegionID = m_NumRegions;
pRegion->SrcID = SrcID;
pRegion->BlockID = 0;
pRegion->NumInstances = Instances;
pRegion->ChromID = ChromID;
pRegion->StartLoci = StartLoci;
pRegion->EndLoci = EndLoci;
pRegion->Strand = bOnMinStrand ? '-' : '+';
return(m_NumRegions);
}


int			// actual length, or if under MinLen then -1, or if start would be < 0 or end < start then returns -2
AdjLoci(bool bOnMinStrand,		// true if element is on '-' strand
		int *pStartLoci,		// starting loci on '+' strand
		int *pEndLoci,			// ending loci on '+' strand
		int OfsLoci,			// offset loci by
		int DeltaLen,			// adjust length by
		int MinLen,				// must be at least this length
		int TruncLen)			// truncate to this length
{
int Ofs;
int DLen;

if(*pStartLoci < 0 || *pStartLoci > *pEndLoci)
	return(-2);

if(OfsLoci == 0 && DeltaLen == 0 && TruncLen == 0)
	return(1 + *pEndLoci -  *pStartLoci);

if(OfsLoci != 0 || DeltaLen != 0)
	{
	if(bOnMinStrand)
		{
		Ofs = OfsLoci * -1;
		DLen = DeltaLen * -1;
		}
	else
		{
		Ofs = OfsLoci;
		DLen = DeltaLen;
		}

	*pStartLoci += Ofs;
	*pEndLoci += Ofs;
	if(bOnMinStrand)
		*pStartLoci += DLen;
	else
		*pEndLoci += DLen;
	}

if(*pStartLoci > *pEndLoci || (!bOnMinStrand && *pStartLoci < 0))
	return(-2);
DLen = 1 + *pEndLoci - *pStartLoci;
if(DLen < MinLen)
	return(-1);
if(DLen > TruncLen)
	{
	if(bOnMinStrand)
		{
		*pStartLoci = 1 + *pEndLoci - TruncLen;
		if(*pStartLoci < 0)
			return(-2);
		}
	else
		*pEndLoci = *pStartLoci + TruncLen - 1;
	DLen= TruncLen;
	}
return(DLen);
}


int
LoadRegions(char *pszRegionFile,	// csv file containing region loci
		int OfsLoci,				// offset region start loci by this many nt
		int DeltaLen,				// change region length by this many nt
		int TruncLength,			// truncate regions to be a maximum of this length
		int MinLength,				// regions to be of this minimum length	
		int	NumIncludeChroms,			// number of chromosome regular expressions to include
		char **ppszIncludeChroms,		// array of include chromosome regular expressions
		int	NumExcludeChroms,			// number of chromosome expressions to exclude
		char **ppszExcludeChroms)		// array of exclude chromosome regular expressions
{
int Rslt;
int Idx;
bool bProcChrom;
int IncChromIdx;
int ExclChromIdx;
int NumFields;
int NumProcessed;

int Len;
int SrcID;
char *pszElType;
char *pszRefSpecies;
char *pszChrom;
char *pszStrand;
int StartLoci;
int EndLoci;
int Instances;
bool bOnMinStrand;
int NumUnderlen;
int NumStartLess0;
int NumFiltChrom;

if(pszRegionFile == NULL || pszRegionFile[0] == '\0')
	return(eBSFSuccess);

#ifdef _WIN32
RegexpMatch mc;
#else
regmatch_t mc;
int RegErr;				// regular expression parsing error
char szRegErr[128];			// to hold RegErr as textual representation ==> regerror();
#endif

#ifdef _WIN32
	Regexp *IncludeChromsRE[cMaxIncludeChroms];	// compiled regular expressions
	Regexp *ExcludeChromsRE[cMaxExcludeChroms];
#else
	regex_t IncludeChromsRE[cMaxIncludeChroms];	// compiled regular expressions
	regex_t ExcludeChromsRE[cMaxExcludeChroms];
#endif

#ifdef _WIN32
try {
	for(Idx=0;Idx < NumIncludeChroms;Idx++)
		{
		IncludeChromsRE[Idx] = new Regexp();
		IncludeChromsRE[Idx]->Parse((const STRCHAR *)ppszIncludeChroms[Idx],false);	// note case insensitive
		}
	}
catch(...)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to process include regexpr chrom '%s'",ppszIncludeChroms[Idx]);
	return(eBSFerrMem);
	}
try {
	for(Idx=0;Idx < NumExcludeChroms;Idx++)
		{
		ExcludeChromsRE[Idx] = new Regexp();
		ExcludeChromsRE[Idx]->Parse((const STRCHAR *)ppszExcludeChroms[Idx],false);	// note case insensitive
		}
	}
catch(...)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to process exclude regexpr chrom '%s'",ppszExcludeChroms[Idx]);
	return(eBSFerrMem);
	}

#else
for(Idx=0;Idx < NumIncludeChroms;Idx++)
	{

	RegErr=regcomp(&IncludeChromsRE[Idx],ppszIncludeChroms[Idx],REG_EXTENDED | REG_ICASE);	// note case insensitive
	if(RegErr)
		{
		regerror(RegErr,&IncludeChromsRE[Idx],szRegErr,sizeof(szRegErr));
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to process include chrom '%s' error: %s",ppszIncludeChroms[Idx],szRegErr);
		return(eBSFerrMem);
		}
	}
for(Idx=0;Idx < NumExcludeChroms;Idx++)
	{
	RegErr = regcomp(&ExcludeChromsRE[Idx],ppszExcludeChroms[Idx],REG_EXTENDED | REG_ICASE);	// note case insensitive
	if(RegErr)
		{
		regerror(RegErr,&ExcludeChromsRE[Idx],szRegErr,sizeof(szRegErr));
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to process exclude chrom '%s' error: %s",ppszExcludeChroms[Idx],szRegErr);
		return(eBSFerrMem);
		}
	}
#endif

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Loading regions from CSV loci file '%s'",pszRegionFile);

m_pCSV = new CCSVFile;
if(m_pCSV == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CCSVfile");
	return(eBSFerrObj);
	}

if((Rslt=m_pCSV->Open(pszRegionFile))!=eBSFSuccess)
	{
	while(m_pCSV->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,m_pCSV->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open file: %s",pszRegionFile);
	Reset();
	return(Rslt);
	}

NumUnderlen = 0;
NumStartLess0 = 0;
NumFiltChrom = 0;
NumProcessed = 0;
while((Rslt=m_pCSV->NextLine()) > 0)	// onto next line containing fields
	{
	NumFields = m_pCSV->GetCurFields();
	if(NumFields < 7)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expected 7+ fields in '%s', GetCurFields() returned '%d'",pszRegionFile,NumFields);
		Rslt = eBSFerrFieldCnt;
		break;
		}
	if(!NumProcessed && m_pCSV->IsLikelyHeaderLine())
		continue;
	NumProcessed += 1;

	m_pCSV->GetText(4,&pszChrom);

			// to be included?
	bProcChrom = false;
	for(IncChromIdx = 0; IncChromIdx < NumIncludeChroms; IncChromIdx++)
		{
#ifdef _WIN32
		if(IncludeChromsRE[IncChromIdx]->Match((const STRCHAR *)pszChrom,&mc))
#else
		if(!regexec(&IncludeChromsRE[IncChromIdx],pszChrom,1,&mc,0))
#endif
			{									
			bProcChrom = true;
			break;
			}
		}

			// if not explicitly included then check if to be excluded?
	if(!bProcChrom && !NumIncludeChroms)
		{
		bProcChrom = true;
		for(ExclChromIdx = 0; ExclChromIdx < NumExcludeChroms; ExclChromIdx++)
			{
#ifdef _WIN32	
			if(ExcludeChromsRE[ExclChromIdx]->Match((const STRCHAR *)pszChrom,&mc))
#else
			if(!regexec(&ExcludeChromsRE[ExclChromIdx],pszChrom,1,&mc,0))
#endif
				{
				bProcChrom = false;
				break;
				}
			}
		}

	if(!bProcChrom)
		{
		NumFiltChrom += 1;
		continue;
		}

	m_pCSV->GetInt(7,&Len);
	m_pCSV->GetInt(1,&SrcID);
	m_pCSV->GetText(2,&pszElType);
	m_pCSV->GetText(3,&pszRefSpecies);
	m_pCSV->GetInt(5,&StartLoci);
	m_pCSV->GetInt(6,&EndLoci);

	bOnMinStrand = false;
	if(NumFields >= 8)					// check if strand has been specified
		{
		m_pCSV->GetText(8,&pszStrand);
		if(pszStrand[0] == '-')		// assume anything other than '-' is on the plus strand
			bOnMinStrand = true;
		}
	Instances = 1;
	if(NumFields >= 11)					// check if likely short reads loci with associated reads count
		m_pCSV->GetInt(11,&Instances);

	if((Rslt = AdjLoci(bOnMinStrand,		// true if element is on '-' strand
		&StartLoci,		// starting loci on '+' strand
		&EndLoci,			// ending loci on '+' strand
		OfsLoci,			// offset loci by
		DeltaLen,			// adjust length by
		MinLength,			// must be at least this length
		TruncLength)) < 1)	// truncate to this length
		{
		if(Rslt == -1)
			NumUnderlen += 1;
		else
			if(Rslt == -2)
				NumStartLess0 += 1;
		continue;
		}

	if((Rslt=AddRegion(SrcID,Instances,pszChrom,StartLoci,EndLoci,bOnMinStrand))<0)
		{
		Reset();
		return(Rslt);
		}
	}
delete m_pCSV;
m_pCSV = NULL;
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Accepted %d regions from %d loaded - %d chrom filtered, %d too short, %d adjusted loci before chrom starts",m_NumRegions,NumProcessed,NumFiltChrom,NumUnderlen,NumStartLess0);
if(m_NumRegions > 1)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Sorting regions...");
	qsort(m_pRegions,m_NumRegions,sizeof(tsRegion),SortRegions);
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Finished sorting regions...");
	}
return(m_NumRegions);
}


// LocateBlockMatches
// Locates block subsequence matches w/o any substitutions or microindels
int
LocateBlockMatches(int SubSeqLen)		// subsequences are of this length
{
int BlockIdx;							// index into m_pSeqBlocks[] of current block being processed
tsSeqBlock *pCurBlock;					// cuurent block being processed
unsigned int CurTargEntry;						// current suffix array entry
etSeqBase BlockSeq[cMaxSeqBlockLen];	// current block sequence
etSeqBase RCBlockSeq[cMaxSeqBlockLen];	// reverse complement of current block sequence
int ProbeIdx;							// index into BlockSeq[] and RCBlockSeq[] of cuurent subsequence being processed
int NumMatches;							// sum of matches for all bases in current subsequence
int SeqIdx;

UINT8 *pSeqVal;
etSeqBase *pSeq;

int AccumReadHits;

unsigned int CurHitLoci;
unsigned int HitLoci;
unsigned int HitChrom;
INT64 HitIdx;
char HitStrand;

CurTargEntry = 0;

// iterate over all blocks
pCurBlock = m_pSeqBlocks;
for(BlockIdx = 0; BlockIdx < m_NumBlocks; BlockIdx++,pCurBlock++)
	{
	pSeqVal = pCurBlock->Seq;
	pSeq = BlockSeq;
	for(SeqIdx = 0; SeqIdx < pCurBlock->SeqLen; SeqIdx++,pSeq++,pSeqVal++)
		*pSeq = (*pSeqVal & 0x07);
	memcpy(RCBlockSeq,BlockSeq,pCurBlock->SeqLen);
	CSeqTrans::ReverseComplement(pCurBlock->SeqLen,RCBlockSeq);

	// iterate all subsequences of length SubSeqLen contained within curent block
	for(ProbeIdx = 0; ProbeIdx <= (pCurBlock->SeqLen - SubSeqLen); ProbeIdx++)
		{
		// if any contained N's in subsequence then treat as though there were multiple matches against this subsequence 
		pSeqVal = &pCurBlock->Seq[ProbeIdx];
		pSeq = &BlockSeq[ProbeIdx];
		NumMatches = 0;
		for(SeqIdx = 0; SeqIdx < SubSeqLen; SeqIdx++,pSeq++,pSeqVal++)
			{
			if(*pSeq >= eBaseN)
				{
				*pSeqVal |= 0x010;	// mark
				break;
				}
			NumMatches += (*pSeqVal & 0x060) >> 5;
			}
		if(SeqIdx != SubSeqLen || NumMatches >= (SubSeqLen * 2))	// if any eBaseN's or all subsequence bases matched at least twice then no point in checking for exact matches on this sequence
			continue;

		// initially try finding a match against the '+' strand
		AccumReadHits = 0;	// will contain num hits
		HitIdx = 0;
		HitStrand = '?';
		while((HitIdx = m_pSfxArray->IterateExacts(&BlockSeq[ProbeIdx],			// probe
					 SubSeqLen,						// probe length
					 HitIdx,						// 0 if starting new sequence, otherwise set to return value of previous successful iteration return
					 &CurTargEntry,
					 &CurHitLoci)) > 0)				// if match then where to return loci
			{
			HitStrand = '+';
			HitLoci = CurHitLoci;
			HitChrom = CurTargEntry;
			if(AccumReadHits++ >= 1)				// getting too many hits?
				break;
			}

		// need to also check on the reverse strand if no hits on the positive strand
		if(AccumReadHits <= 1)	// not interested in multiple hits....
			{
			HitIdx = 0;
			while((HitIdx = m_pSfxArray->IterateExacts( &RCBlockSeq[ProbeIdx],			// probe
						 SubSeqLen,						// probe length
						 HitIdx,						// 0 if starting new sequence, otherwise set to return value of previous successful iteration return
						 &CurTargEntry,
						 &CurHitLoci)) > 0)				// if match then where to return loci
				{
				HitStrand = '-';
				HitLoci = CurHitLoci;
				HitChrom = CurTargEntry;
				if(AccumReadHits++ >= 1)					
					break;
				}
			}
		if(!AccumReadHits)
			continue;

		if(AccumReadHits == 1 && pCurBlock->ChromID == HitChrom)
			if(pCurBlock->StartLoci + ProbeIdx == HitLoci)
			continue;

		pSeqVal = &pCurBlock->Seq[ProbeIdx];
		for(SeqIdx = 0; SeqIdx < SubSeqLen; SeqIdx++,pSeqVal++)
			{
			NumMatches = (*pSeqVal & 0x060) >> 5;
			*pSeqVal &= ~0x060;
			NumMatches += AccumReadHits;
			if(NumMatches > 2)
				NumMatches = 2;
			*pSeqVal |= NumMatches << 5;
			}
		}
	}
return(eBSFSuccess);
}

// SortRegions
// sort regions by ChromID-->StartLoci-->EndLoci ascending
int SortRegions( const void *arg1, const void *arg2)
{
tsRegion *pR1 = (tsRegion *)arg1;
tsRegion *pR2 = (tsRegion *)arg2;
if(pR1->ChromID < pR2->ChromID)
	return(-1);
if(pR1->ChromID > pR2->ChromID)
	return(1);
if(pR1->StartLoci < pR2->StartLoci)
	return(-1);
if(pR1->StartLoci > pR2->StartLoci)
	return(1);
return(pR1->EndLoci - pR2->EndLoci);
}