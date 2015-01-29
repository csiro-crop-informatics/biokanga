// genhyperconserved.cpp : Defines the entry point for the console application.
// Locates all hyper and ultraconserved elements according to parameterised filtering critera
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

const unsigned int cProgVer = 303;		// increment with each release


const int cMAtotMaxSeqAlignLen = 0x0fffffff; // total (over all aligned species) max seq length that can be buffered in concatenated seqs

const int cMinCoreLen = 1;				// allow core lengths to be specified down to cMinCoreLen
const int cDfltMinCoreLen= 50;			// if core lengths not specified then default to cDfltMinCoreLen
const int cMaxMinCoreLen = 10000;		// minimum core lengths can be specified upto this length

const int cMinIdentity = 50;			// accept down to 50% identity
const int cMaxIdentity = 100;			// can't do any better than 100% identity!
const int cDfltIdentity = 90;			// use 90% identity if non specified and hypercores are to be processed
const int cRandWalk100Score = 10000;	// random walk score for 100% identity

const int cMaxMismatches = 500;			// total number of mismatches allowed in any hypercore before terminating that hypercore
const int cDfltMaxMismatches = 100;		// default if non specified

const int cMinMismatchHistLen = 10;		// minimal sized window over which mismatch history is maintained
const int cMaxMismatchHistLen = 500;	// maximal sized window over which mismatch history is maintained

const int cMaxIncludeFiles = 10;		// maximun number of include region filter files
const int cMaxExcludeFiles = 10;		// maximun number of exclude region filter files

const int cMaxIncludeChroms = 20;		// max number of include chromosomes regular expressions
const int cMaxExcludeChroms = 20;		// max number of exclude chromosomes regular expressions

const int cMinMatchDistSegments  = 4;	// min match distribution profile segments (applies if eProcModeOutspecies)
const int cDfltMatchDistSegments = 10;	// default match distribution profile segments (applies if eProcModeOutspecies)
const int cMaxMatchDistSegments  = 100;	// max match distribution profile segments (applies if eProcModeOutspecies)


typedef enum eProcMode {
	eProcModeStandard = 0,				// default processing
	eProcModeSummary,					// summary processing
	eProcModeOutspecies					// same as default but with outspecies species additional detail
} etProcMode;

CStopWatch gStopWatch;
CDiagnostics gDiagnostics;				// for writing diagnostics messages to log file
char gszProcName[_MAX_FNAME];			// process name


typedef struct TAG_sLenRangeClass {
	int ID;					// uniquely identifies this range
	int Min;				// minimum length in this range
	int Max;				// maximum length in this range
	const char *pszDescr;			// descriptive text
	}tsLenRangeClass;


typedef struct TAG_sExcludeSpeciesChrom {
	tSpeciesID SpeciesID;		// which species
	tChromID ChromID;			// chromosome is not to be processed
	} tsExcludeSpeciesChrom;


typedef struct TAG_sDistSeg {
	int Matches;			// number of exact matches in this segment
	int Mismatches;			// number of mismatches
	int InDels;				// total number of InDel bases
	int Unaligned;			// number of unaligned bases
	} tsDistSeg;

typedef struct TAG_sProcParams 
	{
	int ProcMode;					// processing mode 0: default, 1: summary stats, 2: outspecies processing
	bool bStatsAvail;
	char *pszSpeciesList;			// comma separated species list starting with reference species
	int NumSpeciesList;				// number of species in species list
	int RefSpeciesIdx;				// current index into pSeqs for the reference species
	char szSpecies[cMaxAlignedSpecies][cMaxDatasetSpeciesChrom];	// species names of interest - other species are sloughed
									// only alignments with species included will be processed
									// first species is the reference species
	int NumCoreSpecies;				// number of core species (in species list priority order) required in an alignment
	int MinAlignSpecies;			// minimum number of species required in an alignment - will be >= NumCoreSpecies
	int MaxNumStatsSpecies;			// max number of species to accumulate stats for
	int NumSpeciesInAlignment;		// actual number of sequences in current alignment
	int MaxAlignIdxSpecies;			// pSeqs[MaxAlignIdxSpecies-1] is last actual alignment sequence
	char szRefChrom[cMaxDatasetSpeciesChrom]; // current reference chromosome
	int RefChromID;					// current reference chromosome identifier
	bool bFiltLoConfidence;			// true if low confidence subseqences to be filtered out
	bool bFilt1stLast;				// true if treat 1st and last subsequences as being low confidence
	int MinIdent;					// treat subsequences of less than this identity as being low confidence
	int MinSubSeqLen;				// subsequences of less than this length are treated as being low confidence
	bool bInDelsAsMismatches;		// treat InDels as if mismatches (one mismatch == one InDel column)
	bool bSloughRefInDels;			// slough columns in which the reference base is an InDel
	int	WindowSize;					// sampling window size
	int NxtOutputOffset;			// when to next output results - i.e end of current window
	int NumDistSegs;				// number of match distribution profile segments
	int SeqOfs;						// offset into sequences
	int SeqLen;						// sequence length
	etSeqBase *pSeqs[cMaxAlignedSpecies];  // ptrs to each sequence
	int MaxSeqAlignLen;				// max length alignment which can be processed (how much mem was alloc'd to pSeq[n])			
	int MinHyperLen;				// minimum hyper core length required
	int MinUltraLen;				// minimum ultra core length required
	int MinCoreLen;					// maximum of either MinHyperLen or MinUltraLen
	int MaxHyperMismatches;			// hyper cores can have upto this number of total mismatches
	int VMismatches;				// number of mismatches in an alignment col to count as a missmatch against MaxMismatches
	int *pCntStepCnts;				// array of stats counters
	int NumCnts;					// number of steps in pCntStepCnts
	int Regions;					// number of regions per step
	CBEDfile *pBiobed;				// if not NULL then opened biobed file for regional characteristics
	int BEDChromID;					// BED chromosome identifier corresponding to RefChromID
	int NumIncludes;				// number of biobed files containing regions to include
	int NumExcludes;				// number of biobed files containing regions to exclude
	CBEDfile *pIncludes[cMaxIncludeFiles];	// if opened biobed files for regions to include - all other regions are to be excluded
	CBEDfile *pExcludes[cMaxExcludeFiles];	// if opened biobed files for regions to exclude 
	int UpDnStreamLen;				// up/dn stream regional length when characterising
	bool bMultipleFeatBits;			// if false then stats only generated if a single feature bit is set - e.g if both exons and introns overlapped then no stat generated
	int hRsltsFile;					// write stats results into this CSV file
	int hCoreCSVRsltsFile;			// write hypercore loci into this CSV file
	int NumIncludeChroms;			// number of chromosomes explicitly defined to be included
	char **ppszIncludeChroms;		// ptr to array of reg expressions defining chroms to include - overides exclude
	int NumExcludeChroms;			// number of chromosomes explicitly defined to be excluded
	char **ppszExcludeChroms;		// ptr to array of reg expressions defining chroms to include
#ifdef _WIN32
	Regexp *IncludeChromsRE[cMaxIncludeChroms];	// compiled regular expressions
	Regexp *ExcludeChromsRE[cMaxExcludeChroms];
#else
	regex_t IncludeChromsRE[cMaxIncludeChroms];	// compiled regular expressions
	regex_t ExcludeChromsRE[cMaxExcludeChroms];
#endif

	bool bAllUniqueSpeciesSeqs;		// true if all species sequences must not overlap any other sequence
	int NumUniqueSpeciesSeqs;		// number of species in UniqueSpeciesSeqs
	char UniqueSpeciesSeqs[cMaxAlignedSpecies][cMaxDatasetSpeciesChrom];	// species for which sequences in alignment blocks must not overlap with sequence in any other block


	int MinIdentity;				// minimum required identity (50-100%) from which MismatchScore and MatchScore were derived  
	int MismatchScore;				// score to reduce random walk score by for mismatches - processing terminates if score drops below 0
	int MatchScore;					// score to increase random walk score by for matches - random walk score limited to cRandWalk100Score;
	} tsProcParams; 

int ParseNumSpecies(char *pszSpeciesList,tsProcParams *pProcParams);
int ParseUniqueSpeciesSeqs(char *pszUniqueSpeciesSeqs,tsProcParams *pProcParams);

bool ProcAlignBlock(int RefChromID,int RefChromOfs,int AlignLen,tsProcParams *pProcParams);
bool ProcAlignBlockSummary(int RefChromID,int RefChromOfs,int AlignLen,tsProcParams *pProcParams);
bool ProcessAlignment(int RefChromID,int ChromOfs,int SeqIdx,int SubSeqLen,tsProcParams *pProcParams);
int ProcessSubSeq(int RefChromID,int ChromOfs,int SeqIdx,int MaxLen,tsProcParams *pProcParams);
bool ChkOutputResults(const char *pszChrom, int ChromOffset, tsProcParams *pProcParams,bool bGenEmptyRows,bool bGenEmptyWindows);
bool OutputResults(const char *pszChrom, int ChromOffset, tsProcParams *pProcParams,bool bGenEmptyRows);
bool ChkOutputSummaryResults(const char *pszChrom, int ChromOffset, tsProcParams *pProcParams,bool bGenEmptyRows,bool bGenEmptyWindows);
bool OutputSummaryResults(const char *pszChrom, int ChromOffset, tsProcParams *pProcParams,bool bGenEmptyRows);
bool OutputHypercore(const char *pszChrom, int ChromStartOffset, int ChromEndOffset, int FeatureBits,
 				int OGUnaligned,int OGMatches,int OGMismatches,int OGInDels,tsDistSeg SegCnts[],tsProcParams *pProcParams);
char *ChkSpeciesChromWellFormed(char *pszSpeciesChroms);
int ReportSummary(int RefChromID,int RefChromOfs,int ProcMode,tsProcParams *pProcParams);
bool IncludeFilter(int RefChromID,int SubRefOfs,int SubRefEndOfs,tsProcParams *pProcParams);
int TrimQuotes(char *pszTxt);
int										// returned blockid to next start loading from
LoadContiguousBlocks(int RefSpeciesID,	// reference species identifier
 			   int  BlockID,			// which block to initially start loading from
			   bool *pbLoaded,			// returned indicator as to if any loaded blocks meet processing requirements
			   int *pRefChromID,		// returned reference chromosome identifier 
			   char *pRefStrand,		// returned reference strand
			   int *pRefAlignLen,		// returned alignment (incl InDels) length
   			   int *pRefChromOfs,		// returned alignment start offset
   			   int *pSpeciesIDs,		// input - species of interest identifier array
			   CMAlignFile *pAlignments,
			   tsProcParams *pProcParams);

int
Process(bool bTargDeps,				// true if process only if any independent src files newer than target
		int ProcMode,					// processing mode 0: default, 1: summary stats only
				 char *pszInputFile,		// bio multialignment (.algn) file to process
					char *pszOutputFile,	// where to write out stats
					char *pszOutputCoreFile, // where to write out the hypercore loci 
					char *pszBiobedFile,	// biobed file containing regional features - exons, introns etc
					int NumIncludeFiles,	// number of include region files
					char **ppszIncludeFiles,// biobed files containing regions to include - default is to exclude none
					int NumExcludeFiles,	// number of exclude region files
					char **ppszExcludeFiles,// biobed file containing regions to exclude - default is to include all
 					char *pszUniqueSpeciesSeqs, // ignore alignment block if these species sequences are not unique
					int WindowSize,			// sampling window size
					int NumCoreSpecies,		// number of core species to be in alignment
					int MinAlignSpecies,	// minimum number of species required in an alignment (includes core species)
					char *pszSpeciesList,	// space or comma separated list of species, priority ordered
					int MinHyperLen,		// minimum hyper core length required
					int MinUltraLen,		// minimum ultra core length required
					int MaxHyperMismatches,	// hyper cores can have upto this number of total mismatches
					int VMismatches,		// number of mismatches in an alignment col to count as a missmatch against MaxMismatches
					int MinIdentity,		// minimum identity required when processing hyperconserved
					int NumDistSegs,		// number of match distribution profile segments
					int RegLen,				// regulatory region length - up/dn stream of 5/3' 
					bool bMultipleFeatBits,	// if true then accept alignments in which multiple feature bits are set
					bool bInDelsAsMismatches, // treat InDels as if mismatches
					bool bSloughRefInDels,	// slough columns in which the reference base is an InDel
					bool bFiltLoConfidence,	// filter out low confidence subsequences
					bool bFilt1stLast,		// treat 1st and last subsequences as being low confidence
					int MinSubSeqLen,		// subsequences of less than this length are treated as being low confidence
					int MinIdent,			// treat subsequences of less than this identity as being low confidence
					int NumIncludeChroms,	// number of chromosomes explicitly defined to be included
					char **ppszIncludeChroms,	// ptr to array of reg expressions defining chroms to include - overides exclude
					int NumExcludeChroms,	// number of chromosomes explicitly defined to be excluded
					char **ppszExcludeChroms);	// ptr to array of reg expressions defining chroms to exclude


int NormaliseInDelColumns(tsProcParams *pProcParams,int AlignLen);

// length range classes
tsLenRangeClass LenRangeClasses[] = {
	{1,0,4,"0-4"},
	{2,5,9,"5-9"},
	{3,10,14,"10-14"},
	{4,15,19,"15-19"},
	{5,20,29,"20-29"},
    {6,30,49,"30-49"},
	{7,50,74,"50-74"},
	{8,75,99,"75-99"},
	{9,100,124,"100-124"},
	{10,125,149,"125-149"},
	{11,150,174,"150-174"},
	{12,175,199,"175-199"},
	{13,200,249,"200-249"},
	{14,250,299,"250-299"},
	{15,300,349,"300-349"},
	{16,350,399,"350-399"},
	{17,400,449,"400-449"},
	{18,450,499,"450-499"},
	{19,500,599,"500-599"},
	{20,600,699,"600-699"},
	{21,700,799,"700-799"},
	{22,800,899,"800-899"},
	{23,900,999,"900-999"},
	{24,1000,1249,"1000-1249"},
	{25,1250,1499,"1250-1499"},
	{26,1500,1749,"1500-1749"},
	{27,1750,1999,"1750-1999"},
	{28,2000,INT_MAX,"2000+"}
};
const int cLenRanges = sizeof(LenRangeClasses)/sizeof(tsLenRangeClass);		  // number of length range classes

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
	return _T("genhyperconserved");
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

int Rslt;
bool bTargDeps;						// true if only process if independent files newer than target

int iProcMode;
int iNumCoreSpecies;
int iMinAlignSpecies;
bool bMultipleFeatBits;
bool bInDelsAsMismatches;
bool bSloughRefInDels;
bool bFiltLoConfidence;

int LenReq;
int Idx;
char szOutputFile[_MAX_PATH];
char szInputFile[_MAX_PATH];
char szOutCoreFile[_MAX_PATH];
char szInputBiobedFile[_MAX_PATH];
int NumIncludeChroms;
char *pszIncludeChroms[cMaxIncludeChroms];
int NumExcludeChroms;
char *pszExcludeChroms[cMaxExcludeChroms];
int NumIncludeFiles;
char *pszIncludeFiles[cMaxExcludeFiles];
int NumExcludeFiles;
char *pszExcludeFiles[cMaxIncludeFiles];
int iWindowSize;
char szSpeciesList[512];
int iNumSpeciesList;
int iRegLen;
int	iMinHyperLen;		// minimum hyper core length required
int iMinUltraLen;		// minimum ultra core length required
int iMaxHyperMismatches;// hyper cores can have upto this number of total mismatches
int	iMinIdentity;		// minimum identity (50-100) required for hypercores
int iDistSegs;			// match distribution profile segments
char szUniqueSpeciesSeqs[512];	// slough blocks in which these species sequences are not unique

// command line args
struct arg_lit  *help    = arg_lit0("hH","help",                "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_int *ScreenLogLevel=arg_int0("S", "ScreenLogLevel",	"<int>","Level of diagnostics written to logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",	"diagnostics log file");

struct arg_int  *RegLen = arg_int0("L","updnstream","<int>",	"length of 5'up or 3'down  stream regulatory region length (default = 2000) 0..1000000");
struct arg_file *InFile = arg_file1("i",NULL,"<file>",			"input from .algn file");
struct arg_file *OutFile = arg_file1("o",NULL,"<file>",			"output to statistics file as CSV");
struct arg_file *OutCoreFile = arg_file0("O",NULL,"<file>",		"output hypercore loci to file as CSV");
struct arg_file *InBedFile = arg_file0("b","bed","<file>",		"characterise regions from biobed file");
struct arg_int  *NumCoreSpecies = arg_int0("m","numcorespecies","<int>", "number of core species in alignment, default is num of species in -s<specieslist>");
struct arg_int  *MinAlignSpecies = arg_int0("M","minalignedspecies","<int>", "min number of species required in alignment (includes core species), default is num of species in -s<specieslist>");
struct arg_str  *SpeciesList = arg_str1("s","species","<string>","species list, ordered by processing priority");
struct arg_lit  *MultipleFeatBits  = arg_lit0("q","multiplefeatbits",	"single featbit (default) or multiple featbits allowed");
struct arg_lit  *InDelsAsMismatches  = arg_lit0("j","indelsasmismatches",	"treat InDels same as mismatches (default is to terminate element if InDel)");
struct arg_lit  *SloughRefInDels  = arg_lit0("k","sloughrefindels",	"slough columns in which the ref base is an InDel (default is to terminate element if InDel)");
struct arg_lit  *FiltLoConfidence = arg_lit0("l","filt",		"filter out low confidence subsequences,( slough subsequence < 15mer and first/last subsequence, subsequence must start/end on identical, identity > 70)");
struct arg_int  *MinHyperLen = arg_int0("N","minhyperlen","<int>",		"minimum (default = 50 or minultralen) required hypercore length (0,10..1000) - will be maximally extended");
struct arg_int  *MinUltraLen = arg_int0("n","minultralen","<int>",		"minimum (default = 50) required ultra length (0,10..1000) - will be maximally extended");
struct arg_int  *MaxHyperMismatches = arg_int0("X","maxmismatches","<int>",	"total number (default = 100) of mismatches allowed in any hypercores (0..500)");
struct arg_int  *MinIdentity = arg_int0("y","minidentity","<int>",	"minimum percentage identity (default = 90) required in hypercore (50-100)");
struct arg_lit  *ChromPer = arg_lit0("c","chromper",			"generate stats for each chromosome -= default is for complete genome");
struct arg_int  *ProcMode = arg_int0("x","procmode","<int>",	"processing mode 0:default, 1:summary, 2:outspecies");

struct arg_str  *UniqueSpeciesSeqs = arg_str0("u","uniquespecieseqs","<string>","Ignore alignment blocks in which these species are not unique (use 'any' for all, or '*' for -s<specieslist>)");

struct arg_int  *DistSegs = arg_int0("d","distsegs","<int>",	"number of match distribution segments (default 10) used when processing outspecies mode");

struct arg_file *ExcludeFile = arg_filen("E","exclude","<file>",0,cMaxExcludeFiles,	"exclude all regions in biobed file from processing ");
struct arg_file *IncludeFile = arg_filen("I","include","<file>",0,cMaxExcludeFiles,	"include all regions (unless specific regions excluded) in biobed file");
struct arg_int  *WindowSize = arg_int0("w","windowsize","<int>","if non-zero then sets fixed size window (in Knt) used to sample along genome");
struct arg_str  *IncludeChroms = arg_strn("z","chromeinclude","<string>",0,cMaxIncludeChroms,"regular expressions defining species.chromosomes to include (overrides exclude) from processing");
struct arg_str  *ExcludeChroms = arg_strn("Z","chromexclude","<string>",0,cMaxExcludeChroms,"regular expressions defining species.chromosomes to exclude from processing");
struct arg_lit  *TargDeps = arg_lit0("D","TargDep",				"Generate target file only if missing or older than any of the independent source files");

struct arg_end *end = arg_end(20);

void *argtable[] = {help,version,FileLogLevel,ScreenLogLevel,LogFile,
					
					InFile,OutFile,OutCoreFile,InBedFile,
					ProcMode,UniqueSpeciesSeqs,
					NumCoreSpecies,MinAlignSpecies,SpeciesList,MultipleFeatBits,
					InDelsAsMismatches,SloughRefInDels,FiltLoConfidence,
					MinHyperLen,MinUltraLen,MaxHyperMismatches,MinIdentity,RegLen,
					ChromPer,DistSegs,ExcludeFile,IncludeFile,WindowSize,IncludeChroms,ExcludeChroms,TargDeps,
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

	iProcMode = ProcMode->count ? ProcMode->ival[0] : eProcModeStandard;
	if(iProcMode < eProcModeStandard || iProcMode > eProcModeOutspecies)
		{
		printf("\nError: Requested processing mode '-x%d' not supported",iProcMode);
		exit(1);
		}

	if(iProcMode == eProcModeSummary)
		{	
		szOutCoreFile[0] = '\0';		// where to write out the hypercore loci 
		iMinHyperLen =0;				// minimum hyper core length required
		iMinUltraLen =0;				// minimum ultra core length required
		iMaxHyperMismatches=0;			// hyper cores can have upto this number of total mismatches
		iMinIdentity=0;					// minimum identity required when processing hyperconserved
		bInDelsAsMismatches=false;		// treat InDels as if mismatches
		bSloughRefInDels =false;		// slough columns in which the reference base is an InDel
		bFiltLoConfidence=false;		// filter out low confidence subsequences
		}
	else	// standard or outspecies processing mode
		{
		if(!MinUltraLen->count && !MinHyperLen->count)	// default will be for ultra processing if neither ultras or hypers not specified
			iMinUltraLen = cDfltMinCoreLen;
		else
			{
			iMinUltraLen = MinUltraLen->count ? MinUltraLen->ival[0] : 0;
			if(iMinUltraLen != 0)
				{
				if(iMinUltraLen < 0)
					{
					printf("\nspecified minimum hyper length '-n%d' < 0, assuming you are not requiring ultra processing",iMinUltraLen);
					iMinUltraLen = 0;
					}
				else
					{
					if(iMinUltraLen < cMinCoreLen)
						{
						printf("\nSpecified minimum ultra length was '-n%d' < %d, assuming you meant '-n%d'",iMinUltraLen,cMinCoreLen,cMinCoreLen);
						iMinUltraLen = cMinCoreLen;
						}
					else
						if(iMinUltraLen > cMaxMinCoreLen)
							{
							printf("\nSpecified minimum ultra length was '-n%d' > %d, assuming you meant '-n%d'",iMinUltraLen,cMaxMinCoreLen,cMaxMinCoreLen);
							iMinUltraLen = cMaxMinCoreLen;
							}
					}
				}
			}
	

		iMinHyperLen = MinHyperLen->count ? MinHyperLen->ival[0] : 0;	// default is for ultra processing only
		if(iMinHyperLen <= 0)
			{
			if(iMinUltraLen == 0)
				{
				printf("\nNothing to do, neither ultra '-n<len>' or hyper '-N<len>' core lengths specified, or lengths are both 0!");
				exit(1);
				}
			if(iMinHyperLen < 0)
				printf("\nminimum hyper length '-n%d' less than 0, assuming you are only processing for ultras",iMinHyperLen);
			iMinHyperLen = 0;
			}

		if(iMinHyperLen > 0)	// Note: will only process minimum required identity if hypercore length was specified > 0
			{
			if(iMinHyperLen < cMinCoreLen)
				{
				printf("\nSpecified minimum hyper length '-N%d' < %d, assuming you meant '-N%d'",iMinHyperLen,cMinCoreLen,cMinCoreLen);
				iMinHyperLen = cMinCoreLen;
				}
			else
				if(iMinHyperLen > cMaxMinCoreLen)
					{
					printf("\nSpecified minimum hyper length '-N%d' > %d, assuming you meant '-N%d'",iMinHyperLen,cMaxMinCoreLen,cMaxMinCoreLen);
					iMinHyperLen = cMaxMinCoreLen;
					}

		if(iMinHyperLen < iMinUltraLen)
			{
			printf("\nSpecified minimum hyper length '-N%d' < minimum ultra length '-n%d', assuming you meant '-N%d'",iMinHyperLen,iMinUltraLen,iMinUltraLen);
			iMinHyperLen = iMinUltraLen;
			}

		iMaxHyperMismatches = MaxHyperMismatches->count ? MaxHyperMismatches->ival[0] : cDfltMaxMismatches;
		if(iMaxHyperMismatches < 0)
			{
			printf("\nMaximun allowed Hypercore mismatches '-X%d' less than 0, assuming you meant to use '-y%d'",iMaxHyperMismatches,cDfltMaxMismatches);
			iMaxHyperMismatches = cDfltMaxMismatches;
			}
		else
			if(iMaxHyperMismatches > cMaxMismatches)
				{
				printf("\nRequested Hypercore mismatches '-X%d' more than than allowed '-X%d', assuming you meant to use '-X%d'",iMaxHyperMismatches,cMaxMismatches,cMaxMismatches);
				iMaxHyperMismatches = cDfltMaxMismatches;
				}

		iMinIdentity = MinIdentity->count ? MinIdentity->ival[0] : cDfltIdentity;
		if(iMinIdentity < cMinIdentity)
			{
			printf("\nMinimum identity specified '-y%d' less than %d%%, assuming you meant to use '-y%d'",iMinIdentity,cMinIdentity,cMinIdentity);
			iMinIdentity = cMinIdentity;
			}
		else
			if(iMinIdentity > 100)
				{
				printf("\nMinimum identity specified '-y%d' more than %d, assuming you meant to use '-y%d'",iMinIdentity,cMaxIdentity,cMaxIdentity);
				iMinIdentity = cMaxIdentity;
				}
			}
		else
			{
			iMaxHyperMismatches = 0;
			iMinIdentity = 100;
			}

		if(OutCoreFile->count)
			strcpy(szOutCoreFile,OutCoreFile->filename[0]);
		else
			szOutCoreFile[0] = '\0';

		bFiltLoConfidence = FiltLoConfidence->count ? true : false;
		bInDelsAsMismatches = InDelsAsMismatches->count ? true : false;
		bSloughRefInDels = SloughRefInDels->count ? true : false;
		}
	

	bMultipleFeatBits = MultipleFeatBits->count ? true : false;

	strcpy(szInputFile,InFile->filename[0]);
	strcpy(szOutputFile,OutFile->filename[0]);

	if(InBedFile->count)
		strcpy(szInputBiobedFile,InBedFile->filename[0]);
	else
		szInputBiobedFile[0] = '\0';

	NumIncludeFiles = IncludeFile->count;
	for(Idx=0;Idx < IncludeFile->count; Idx++)
		{
		LenReq = (int)strlen(IncludeFile->filename[Idx]);
		pszIncludeFiles[Idx] = new char [LenReq+1];
		strcpy(pszIncludeFiles[Idx],IncludeFile->filename[Idx]);
		TrimQuotes(pszIncludeFiles[Idx]);
		}
	NumExcludeFiles = ExcludeFile->count;
	for(Idx=0;Idx < ExcludeFile->count; Idx++)
		{
		LenReq = (int)strlen(ExcludeFile->filename[Idx]);
		pszExcludeFiles[Idx] = new char [LenReq+1];
		strcpy(pszExcludeFiles[Idx],ExcludeFile->filename[Idx]);
		TrimQuotes(pszExcludeFiles[Idx]);
		}

	NumIncludeChroms = IncludeChroms->count;
	for(Idx=0;Idx < IncludeChroms->count; Idx++)
		{
		LenReq = (int)strlen(IncludeChroms->sval[Idx]);
		pszIncludeChroms[Idx] = new char [LenReq+1];
		strcpy(pszIncludeChroms[Idx],IncludeChroms->sval[Idx]);
		TrimQuotes(pszIncludeChroms[Idx]);
		}

	NumExcludeChroms = ExcludeChroms->count;
	for(Idx=0;Idx < ExcludeChroms->count; Idx++)
		{
		LenReq = (int)strlen(ExcludeChroms->sval[Idx]);
		pszExcludeChroms[Idx] = new char [LenReq+1];
		strcpy(pszExcludeChroms[Idx],ExcludeChroms->sval[Idx]);
		TrimQuotes(pszExcludeChroms[Idx]);
		}

	iWindowSize = WindowSize->count ? WindowSize->ival[0] : 0;
	if(iWindowSize < 0)
		{
		printf("\nAllowed sampling result window size '-w%d' less than 0, assuming you meant to use '-w0'",iWindowSize);
		iWindowSize = 0;
		}

	if(ChromPer->count && iWindowSize == 0) // a hack to force per chromsome stats generation
		iWindowSize = 0x1fffffff;

	iRegLen = RegLen->count ? RegLen->ival[0] : cDfltRegLen;
	if(iRegLen < cMinRegLen)
		{
		printf("\nRegulatory region length '-L%d' less than minimum %d, assuming you meant to use '-L%d'",iRegLen,cMinRegLen,cMinRegLen);
		iRegLen = cMinRegLen;
		}
	else
		{
		if(iRegLen > cMaxRegLen)
			{
			printf("\nRegulatory region length '-L%d' more than maximum %d, assuming you meant to use '-L%d'",iRegLen,cMaxRegLen,cMaxRegLen);
			iRegLen = cMaxRegLen;
			}
		}
	
	if(UniqueSpeciesSeqs->count)
		{
		strcpy(szUniqueSpeciesSeqs,UniqueSpeciesSeqs->sval[0]);
		TrimQuotes(szUniqueSpeciesSeqs);
		}
	else
		szUniqueSpeciesSeqs[0] = '\0';

	if(!SpeciesList->count)
		{
		printf("\nError: Species list '-s<specieslist>' is empty\n");
		exit(1);
		}
	strcpy(szSpeciesList,SpeciesList->sval[0]);
	TrimQuotes(szSpeciesList);

	iNumSpeciesList = ParseNumSpecies(szSpeciesList,NULL);
	if(iProcMode == eProcModeOutspecies)
		{
		if(iNumSpeciesList < 3)
			{
			printf("\nError: At least three (2 core + outspecies) species must be specified in Outspecies mode\n");
			exit(1);
			}
		iDistSegs = DistSegs->count ? DistSegs->ival[0] : cDfltMatchDistSegments;
		if(iDistSegs < cMinMatchDistSegments)
			{
			printf("\nWarning: too few distribution segments specified with '-d%d', assume you meant %d segments", iDistSegs,cMinMatchDistSegments);
			iDistSegs = iMinUltraLen;
			}
		else
			if(iDistSegs > min(iMinUltraLen,cMaxMatchDistSegments))
				{
				printf("\nWarning: too many distribution segments specified with '-d%d', assume you meant %d segments", iDistSegs,min(iMinUltraLen,cMaxMatchDistSegments));
				iDistSegs = min(iMinUltraLen,cMaxMatchDistSegments);
				}
		}
	else
		{
		if(iNumSpeciesList < 2)
			{
			printf("\nError: At least two species must be specified\n");
			exit(1);
			}
		iDistSegs = 0;
		}

	if(iProcMode != eProcModeOutspecies)
		{
		iNumCoreSpecies = NumCoreSpecies->count ? NumCoreSpecies->ival[0] : iNumSpeciesList;
		if(iNumCoreSpecies < 2)
			iNumCoreSpecies = 2;
		if(iNumCoreSpecies > iNumSpeciesList)
			{
			printf("NumCoreSpecies %d is more than number of species actually specified %d in species list",iNumCoreSpecies,iNumSpeciesList);
			exit(1);
			}

		iMinAlignSpecies = MinAlignSpecies->count ? MinAlignSpecies->ival[0] : iNumCoreSpecies;
		if(iMinAlignSpecies < iNumCoreSpecies || iMinAlignSpecies > iNumSpeciesList)
			{
			printf("MinAlignSpecies %d must be in range %d..%d",iMinAlignSpecies,iNumCoreSpecies,iNumSpeciesList);
			exit(1);
			}
		}
	else
		{
		iNumCoreSpecies = iNumSpeciesList - 1;
		iMinAlignSpecies = iNumSpeciesList;
		}

	if (TargDeps->count > 0)
		bTargDeps = true;
	else
		bTargDeps = false;

		// now that command parameters have been parsed then initialise diagnostics log system
	if(!gDiagnostics.Open(szLogFile,(etDiagLevel)iScreenLogLevel,(etDiagLevel)iFileLogLevel,true))
		{
		printf("\nError: Unable to start diagnostics subsystem.");
		if(szLogFile[0] != '\0')
			printf(" Most likely cause is that logfile '%s' can't be opened/created",szLogFile);
		exit(1);
		}

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Version: %d.%2.2d Processing parameters:",cProgVer/100,cProgVer%100);
	if(bTargDeps)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Generate target file only if missing or older than any of the independent source files");
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"bio multialignment (.algn) file to process: '%s'",szInputFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"where to write out stats: '%s'",szOutputFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"biobed file containing regional features: '%s'",szInputBiobedFile);
	switch(iProcMode) {
		case eProcModeStandard:
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Processing mode: standard");
			break;

		case eProcModeSummary:
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Processing mode: summary");
			break;

		case eProcModeOutspecies:
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Processing mode: outspecies");
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Match distribution profile segments: %d",iDistSegs);
			break;
		}

	if(szUniqueSpeciesSeqs[0] != '\0')
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Ignore alignment blocks where sequences not unique for: '%s'",szUniqueSpeciesSeqs);

	if(iProcMode != eProcModeSummary)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"where to write out the hypercore loci: '%s'",szOutCoreFile);
	for(Idx = 0; Idx < NumIncludeFiles; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"biobed file containing regions to include: '%s'",pszIncludeFiles[Idx]); 
	for(Idx = 0; Idx < NumExcludeFiles; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"biobed file containing regions to exclude: '%s'",pszExcludeFiles[Idx]);
	if(iWindowSize)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"sampling window size: %d",iWindowSize);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"number of core species to be in alignment: %d",iNumCoreSpecies);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"minimum number of species to be in alignment: %d",iMinAlignSpecies);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"species list: '%s'",	szSpeciesList);
	if(iProcMode != eProcModeSummary)
		{
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"minimum hyper core length required: %d",iMinHyperLen);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"minimum ultra core length required: %d",iMinUltraLen);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"hyper cores can have upto this number of total mismatches: %d",iMaxHyperMismatches);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"minimum identity required when processing hyperconserved: %d",iMinIdentity);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"treat InDels as if mismatches: %s",bInDelsAsMismatches ? "yes" : "no");
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"slough columns in which the reference base is an InDel: %s",bSloughRefInDels ? "yes" : "no");
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"filter out low confidence subsequences: %s",bFiltLoConfidence ? "yes" : "no");
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"treat 1st and last subsequences as being low confidence: %s",bFiltLoConfidence ? "yes" : "no");
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"subsequences of less than this length are treated as being low confidence: %d",bFiltLoConfidence ? 15 : 0);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"treat subsequences of less than this identity as being low confidence: %d",bFiltLoConfidence? 70 : 0);
		}
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"regulatory region length: %d",iRegLen);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"accept alignments in which multiple feature bits are set: %s",bMultipleFeatBits ? "yes" : "no");
	for(Idx = 0; Idx < NumIncludeChroms; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"reg expressions defining chroms to include: '%s'",pszIncludeChroms[Idx]);
	for(Idx = 0; Idx < NumExcludeChroms; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"reg expressions defining chroms to exclude: '%s'",pszExcludeChroms[Idx]); 

	gStopWatch.Start();
#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	Rslt = Process(bTargDeps,				// true if process only if any independent src files newer than target
					iProcMode,		// processing mode 0: default, 1: summary stats only, 2: outspecies processing
					szInputFile,		// bio multialignment (.algn) file to process
					szOutputFile,		// where to write out stats
					szOutCoreFile,		// where to write out the hypercore loci 
					szInputBiobedFile,	// biobed file containing regional features - exons, introns etc
					NumIncludeFiles,	// number of include region files
					pszIncludeFiles,	// biobed files containing regions to include - default is to exclude none
					NumExcludeFiles,	// number of exclude region files
					pszExcludeFiles,	// biobed file containing regions to exclude - default is to include all
 					szUniqueSpeciesSeqs, // ignore alignment blocks if these species sequences are not unique
					iWindowSize,		// sampling window size
					iNumCoreSpecies,	// number of core species to be in alignment
					iMinAlignSpecies,	// minimum number of species to be in alignment
					szSpeciesList,		// space or comma separated list of species, priority ordered
					iMinHyperLen,		// minimum hyper core length required
					iMinUltraLen,		// minimum ultra core length required
					iMaxHyperMismatches,	// hyper cores can have upto this number of total mismatches
					1,					// number of mismatches in an alignment col to count as a missmatch against MaxMismatches
					iMinIdentity,		// minimum identity required when processing hyperconserved
					iDistSegs,			// number of match distribution profile segments
					iRegLen,			// regulatory region length - up/dn stream of 5/3' 
					bMultipleFeatBits,	// if true then accept alignments in which multiple feature bits are set
					bInDelsAsMismatches, // treat InDels as if mismatches
					bSloughRefInDels,		// slough columns in which the reference base is an InDel
					bFiltLoConfidence,		// filter out low confidence subsequences
					bFiltLoConfidence,		// treat 1st and last subsequences as being low confidence
					bFiltLoConfidence ? 15 : 0, // subsequences of less than this length are treated as being low confidence
					bFiltLoConfidence? 70 : 0, // treat subsequences of less than this identity as being low confidence
					NumIncludeChroms,	// number of chromosomes explicitly defined to be included
					pszIncludeChroms,	// ptr to array of reg expressions defining chroms to include - overides exclude
					NumExcludeChroms,	// number of chromosomes explicitly defined to be excluded
					pszExcludeChroms);	// ptr to array of reg expressions defining chroms to be excluded

	gStopWatch.Stop();
	Rslt = Rslt < 0 ? 1 : 0;
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Exit Code: %d Total processing time: %s",Rslt,gStopWatch.Read());
	exit(Rslt);
	}
else
	{
	arg_print_errors(stdout,end,gszProcName);
	arg_print_syntax(stdout,argtable,"\nend of help\n");
	exit(1);
	}
}

// TrimQuotes
// Removes any leading and trailing quotes from pszTxt
int
TrimQuotes(char *pszTxt)
{
char *pDst = pszTxt;
char Chr;
int Len = 0;
while((Chr = *pszTxt++))	
	{
	if((!Len || *pszTxt == '\0') && (Chr == '"' || Chr == '\''))
		continue;
	*pDst++ = Chr;
	Len++;
	}
*pDst = Chr;
return(Len);
}

// GetLengthRangeClass
// Returns ptr to length range class for specified Length, or NULL if can't classify into a range
tsLenRangeClass *GetLengthRangeClass(int Length)
{
int Idx;
tsLenRangeClass *pRange = LenRangeClasses;
for(Idx = 0; Idx < cLenRanges; Idx++,pRange++)
	if(Length >= pRange->Min && Length <= pRange->Max)
		return(pRange);
return(NULL);
}

// GetRangeClass
// Returns ptr to length range class for specified range identifier, or NULL if can't classify into a range
tsLenRangeClass *GetRangeClass(int RangeID)
{
if(RangeID < 1 || RangeID > cLenRanges)
	return(NULL);
return(&LenRangeClasses[RangeID-1]);
}

// ParseNumSpecies
// Initialises pProcParams with parsed species names in space or comma delimited list ptd at by pszSpeciesList
// Returns number of species
int
ParseNumSpecies(char *pszSpeciesList,tsProcParams *pProcParams)
{
// parse out species list
char Chr;
char *pSpecies;
int NumSpecies = 0;
bool InToken = false;
if(pszSpeciesList == NULL || *pszSpeciesList == '\0')
	return(0);

while(Chr = *pszSpeciesList++)
	{
	if(Chr == '"' || Chr == '\'') // change any single or double quotes into spaces
		Chr = ' ';
	if(isspace(Chr) || Chr==',')
		{
		if(!InToken)			// slough whitespace or ',' if not inside a token parse
			continue;
		InToken = false;
		pszSpeciesList[-1] = '\0';
		if(pProcParams != NULL)
			{
			strncpy(pProcParams->szSpecies[NumSpecies],pSpecies,cMaxDatasetSpeciesChrom);
			pProcParams->szSpecies[NumSpecies][cMaxDatasetSpeciesChrom-1] = '\0';
			}
		pszSpeciesList[-1] = Chr;
		NumSpecies++;
		if(NumSpecies >= cMaxAlignedSpecies)
			break;
		continue;
		}
	if(!InToken)			// if not already inside token then start token 
		{
		pSpecies = pszSpeciesList-1;
		InToken = true;
		}
	}
if(InToken)
	{
	pszSpeciesList[-1] = '\0';
	if(pProcParams != NULL)
		{
		strncpy(pProcParams->szSpecies[NumSpecies],pSpecies,cMaxDatasetSpeciesChrom);
		pProcParams->szSpecies[NumSpecies][cMaxDatasetSpeciesChrom-1] = '\0';
		}
	pszSpeciesList[-1] = Chr;
	NumSpecies++;
	}
return(NumSpecies);
}

// ParseUniqueSpeciesSeqs
// Initialises pProcParams with parsed species names in space or comma delimited list ptd at by pszSpeciesList
// Returns number of species
int
ParseUniqueSpeciesSeqs(char *pszUniqueSpeciesSeqs,tsProcParams *pProcParams)
{
// parse out species list for which sequences must be unique
char Chr;
char *pSpecies;
int NumSpecies = 0;
bool InToken = false;

pProcParams->bAllUniqueSpeciesSeqs = false;
pProcParams->NumUniqueSpeciesSeqs = 0;
if(pszUniqueSpeciesSeqs == NULL || *pszUniqueSpeciesSeqs == '\0')
	return(0);

if(!stricmp(pszUniqueSpeciesSeqs,"all"))
	{
	pProcParams->bAllUniqueSpeciesSeqs = true;
	return(cMaxAlignedSpecies);
	}

while(Chr = *pszUniqueSpeciesSeqs++)
	{
	if(Chr == '"' || Chr == '\'') // change any single or double quotes into spaces
		Chr = ' ';
	if(isspace(Chr) || Chr==',')
		{
		if(!InToken)			// slough whitespace or ',' if not inside a token parse
			continue;
		InToken = false;
		pszUniqueSpeciesSeqs[-1] = '\0';
		if(pProcParams != NULL)
			{
			strncpy(pProcParams->UniqueSpeciesSeqs[NumSpecies],pSpecies,cMaxDatasetSpeciesChrom);
			pProcParams->UniqueSpeciesSeqs[NumSpecies][cMaxDatasetSpeciesChrom-1] = '\0';
			}
		pszUniqueSpeciesSeqs[-1] = Chr;
		NumSpecies++;
		if(NumSpecies >= cMaxAlignedSpecies)
			break;
		continue;
		}
	if(!InToken)			// if not already inside token then start token 
		{
		pSpecies = pszUniqueSpeciesSeqs-1;
		InToken = true;
		}
	}
if(InToken)
	{
	pszUniqueSpeciesSeqs[-1] = '\0';
	if(pProcParams != NULL)
		{
		strncpy(pProcParams->UniqueSpeciesSeqs[NumSpecies],pSpecies,cMaxDatasetSpeciesChrom);
		pProcParams->UniqueSpeciesSeqs[NumSpecies][cMaxDatasetSpeciesChrom-1] = '\0';
		}
	pszUniqueSpeciesSeqs[-1] = Chr;
	NumSpecies++;
	}
pProcParams->NumUniqueSpeciesSeqs = NumSpecies;
return(NumSpecies);
}




const int cMaxExcludeHistory = 100;
typedef struct TAG_sExcludeEl {
	struct TAG_sExcludeEl *pNext;
	struct TAG_sExcludeEl *pPrev;
	int SpeciesID;		// identifies species
	int ChromID;		// identifies chromosome
	bool bExclude;		// true if to be excluded, false if not
	} tsExcludeEl;

int gNumExcludeEls = 0;		// current number of elements in gExcludeChroms
tsExcludeEl *gpMRA = NULL;		// pts to most recently accessed or added
tsExcludeEl *gpLRA = NULL;		// pts to least recently accessed
tsExcludeEl gExcludeChroms[cMaxExcludeHistory];

// AddExcludeHistory
// Adds a tsExcludeEl to the cached history
// The newly added element will be the MRA
bool
AddExcludeHistory(int SpeciesID,int ChromID,bool bExclude)
{
tsExcludeEl *pEl;
if(gNumExcludeEls < cMaxExcludeHistory)
	pEl = &gExcludeChroms[gNumExcludeEls++];
else
	{
	pEl = gpLRA;		// reuse the least recently accessed element
	gpLRA = pEl->pPrev;	// 2nd LRA now becomes the LRA
	gpLRA->pNext = NULL;
	}
if(gpMRA != NULL)
	gpMRA->pPrev = pEl;
pEl->pNext = gpMRA;
pEl->pPrev = NULL;
gpMRA = pEl;
if(gpLRA == NULL)
	gpLRA = pEl;
pEl->bExclude = bExclude;
pEl->SpeciesID = SpeciesID;
pEl->ChromID = ChromID;
return(bExclude);
}

// LocateExclude
// Locates - starting from the MRA - a tsExcludeEl which matches on SpeciesID and ChromID
// If matches then this tsExcludeEl is made the MRA
// Returns ptr to matching tsExcludeEl if located or NULL
tsExcludeEl *
LocateExclude(int SpeciesID,int ChromID)
{
tsExcludeEl *pEl;
pEl = gpMRA;
while(pEl != NULL)
	{
	if(pEl->SpeciesID == SpeciesID && pEl->ChromID == ChromID)
		{
		if(gNumExcludeEls ==1 || gpMRA == pEl)	// if only, or already the MRA then no need for any relinking 
			return(pEl);

		if(gpLRA == pEl)						// if was the LRA then the 2nd LRA becomes the LRA
			{
			gpLRA = pEl->pPrev;
			gpLRA->pNext = NULL;
			}
		else									// not the LRA, and not the MRA
			{
			pEl->pPrev->pNext = pEl->pNext;
			pEl->pNext->pPrev = pEl->pPrev;
			}
		gpMRA->pPrev = pEl;
		pEl->pNext = gpMRA;
		pEl->pPrev = NULL;
		gpMRA = pEl;
		return(pEl);
		}
	pEl = pEl->pNext;
	}
return(NULL);
}


// ExcludeThisChrom
// Returns true if SpeciesID.ChromosomeID is to be excluded from processing
// ExcludeThisChrom
bool
ExcludeThisChrom(CMAlignFile *pAlignments,int SpeciesID,int ChromID,tsProcParams *pProcParams)
{
#ifdef _WIN32
RegexpMatch mc;
#else
regmatch_t mc;
#endif

char *pszSpecies;
char *pszChrom;
char szSpeciesChrom[200];
int Idx;
if(!pProcParams->NumExcludeChroms && !pProcParams->NumIncludeChroms)
	return(false);
// check if this species and chromosome are already known to be included/excluded
tsExcludeEl *pEl;
if((pEl = LocateExclude(SpeciesID,ChromID))!=NULL)
	return(pEl->bExclude);
// haven't seen this species or chromosome before - or else they have been discarded from history...
pszSpecies = pAlignments->GetSpeciesName(SpeciesID);
pszChrom = pAlignments->GetChromName(ChromID);
sprintf(szSpeciesChrom,"%s.%s",pszSpecies,pszChrom);
// to be included?
for(Idx = 0; Idx < pProcParams->NumIncludeChroms; Idx++)
	{
#ifdef _WIN32
	if(pProcParams->IncludeChromsRE[Idx]->Match(szSpeciesChrom,&mc))
#else
	if(!regexec(&pProcParams->IncludeChromsRE[Idx],szSpeciesChrom,1,&mc,0))
#endif
		return(AddExcludeHistory(SpeciesID,ChromID,false));
	}

// to be excluded?
for(Idx = 0; Idx < pProcParams->NumExcludeChroms; Idx++)
#ifdef _WIN32	
	if(pProcParams->ExcludeChromsRE[Idx]->Match(szSpeciesChrom,&mc))
#else
	if(!regexec(&pProcParams->ExcludeChromsRE[Idx],szSpeciesChrom,1,&mc,0))
#endif
		return(AddExcludeHistory(SpeciesID,ChromID,true));
// if not explicitly included or excluded then default is to assume include
return(AddExcludeHistory(SpeciesID,ChromID,false));
}

int 
ProcessAlignments(char *pszMAF,			 // source bioseq multialignment file
				  tsProcParams *pProcParams) // processing parameters
{
int RefSpeciesID;
int RefChromID;
int BEDChromID;
int PrevRefChromID;
int PrevDispRefChromID;
char *pszRefChrom;
int RefChromOfs;
char RefStrand;
int SpeciesIDs[cMaxAlignedSpecies];
int Rslt;
int RefAlignLen;
CMAlignFile *pAlignments;
int Idx;
int CurBlockID;
bool bLoaded;

if((pAlignments = new CMAlignFile())==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create new instance of CMAlignFile");
	return(eBSFerrObj);
	}

if((Rslt=pAlignments->Open(pszMAF))!=eBSFSuccess)
	{
	while(pAlignments->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,pAlignments->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open MAF file %s\n",pszMAF);
	return(Rslt);
	}

// ensure all species are represented in multispecies alignment file plus get their species identifiers
for(Idx = 0; Idx < pProcParams->NumSpeciesList; Idx++)
	{
	if((Rslt = SpeciesIDs[Idx] = pAlignments->LocateSpeciesID(pProcParams->szSpecies[Idx]))<1)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Species '%s' not represented in %s",pProcParams->szSpecies[Idx],pszMAF);
		Rslt = pAlignments->GetNumSpecies();
		for(Idx=1;Idx<=Rslt;Idx++)
			gDiagnostics.DiagOut(eDLFatal,gszProcName," Represented species: %s",pAlignments->GetSpeciesName(Idx));
		delete pAlignments;
		return(eBSFerrEntry);
		}
	if(!Idx)	// reference species is always the first species in the species name list
		RefSpeciesID = SpeciesIDs[0];
	}

if(pProcParams->bAllUniqueSpeciesSeqs)
	pAlignments->SetAllConfFilt(true);
else
	{
	pAlignments->SetAllConfFilt(false);
	if(pProcParams->NumUniqueSpeciesSeqs)
		{
		for(Idx = 0; Idx < pProcParams->NumUniqueSpeciesSeqs; Idx++)
			{
			if((Rslt = (int)pAlignments->LocateSpeciesID(pProcParams->UniqueSpeciesSeqs[Idx]))<1)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unique sequence species '%s' not represented in %s",pProcParams->UniqueSpeciesSeqs[Idx],pszMAF);
				Rslt = pAlignments->GetNumSpecies();
				for(Idx=1;Idx<=Rslt;Idx++)
					gDiagnostics.DiagOut(eDLFatal,gszProcName," Represented species: %s",pAlignments->GetSpeciesName(Idx));
				delete pAlignments;
				return(eBSFerrEntry);
				}
			pAlignments->SetConfFilt((tSpeciesID)Rslt,true);
			}
		}
	}

// iterate over reference blocks which are sorted by chrom then offset
CurBlockID = 0;
PrevRefChromID = 0;
PrevDispRefChromID = 0;
BEDChromID = 0;
pProcParams->RefSpeciesIdx = 0;		// reference sequence will always be 1st
pProcParams->NxtOutputOffset = pProcParams->WindowSize;

while(CurBlockID >= 0 && ((CurBlockID =						// returned blockid to next start loading from
	LoadContiguousBlocks(RefSpeciesID,	// reference species identifier
 			   CurBlockID,			// which block to initially start loading from
			   &bLoaded,			// returned indicator as to if any loaded blocks meet processing requirements
			   &RefChromID,			// returned reference chromosome identifier 
			   &RefStrand,			// returned reference strand
			   &RefAlignLen,		// returned alignment (incl InDels) length
   			   &RefChromOfs,		// returned alignment start offset
   			   SpeciesIDs,			// input - species of interest identifier array
			   pAlignments,
			   pProcParams)) > 0 || (CurBlockID == eBSFerrAlignBlk && RefAlignLen > 0)))
	{
	if(RefChromID > 0 && RefChromID != PrevDispRefChromID)
		{
		if(PrevDispRefChromID > 0)
			{
			if(pProcParams->ProcMode == eProcModeSummary)
				ChkOutputSummaryResults(pProcParams->szRefChrom, -1,pProcParams,false,false);
			else
				ChkOutputResults(pProcParams->szRefChrom, -1,pProcParams,false,false);
			}
		PrevDispRefChromID = RefChromID;
		pszRefChrom = pAlignments->GetChromName(RefChromID);

		strcpy(pProcParams->szRefChrom,pszRefChrom);
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing chromosome %s",pszRefChrom);
		pProcParams->NxtOutputOffset = pProcParams->WindowSize;
		if(pProcParams->pBiobed != NULL)
			pProcParams->BEDChromID = pProcParams->pBiobed->LocateChromIDbyName(pProcParams->szRefChrom);
		else
			pProcParams->BEDChromID = 0;
		pProcParams->RefChromID = RefChromID;
		}

	if(!bLoaded)
		continue;

		// not interested if the alignment length would be too short unless extended stats being generated
	if(pProcParams->ProcMode != eProcModeSummary && RefAlignLen < pProcParams->MinCoreLen)
		continue;
	if(pProcParams->bFiltLoConfidence && RefAlignLen < pProcParams->MinSubSeqLen)
		continue;

		// here we need to normalise the alignments so that there will be no case of all InDels in
		// any column which can occur if none of the processed species is not the reference species
	if((RefAlignLen = NormaliseInDelColumns(pProcParams,RefAlignLen))< pProcParams->MinCoreLen)
		continue;
	if(pProcParams->bFiltLoConfidence && RefAlignLen < pProcParams->MinSubSeqLen)
		continue;

	if(RefChromID != PrevRefChromID)
		PrevRefChromID = RefChromID;

	if(pProcParams->ProcMode == eProcModeSummary)
		{
		ChkOutputSummaryResults(pProcParams->szRefChrom, RefChromOfs,pProcParams,false,false);
		if(!ProcAlignBlockSummary(RefChromID,RefChromOfs,RefAlignLen,pProcParams))
			break;
		}
	else
		{
		ChkOutputResults(pProcParams->szRefChrom, RefChromOfs,pProcParams,false,false);
		if(!ProcAlignBlock(RefChromID,RefChromOfs,RefAlignLen,pProcParams))
			break;
		}
	}


delete pAlignments;
return(eBSFSuccess);
}

int										// returned blockid to next start loading from
LoadContiguousBlocks(int RefSpeciesID,	// reference species identifier
 			   int  BlockID,			// which block to initially start loading from
			   bool *pbLoaded,			// returned indicator as to if any loaded blocks meet processing requirements
			   int *pRefChromID,		// returned reference chromosome identifier 
			   char *pRefStrand,		// returned reference strand
			   int *pRefAlignLen,		// returned alignment (incl InDels) length
   			   int *pRefChromOfs,		// returned alignment start offset
   			   int *pSpeciesIDs,		// input - species of interest identifier array
			   CMAlignFile *pAlignments,
			   tsProcParams *pProcParams)
{
int CurBlockID;
int PrvBlockID;
int bFirst;
int CurNumAlignments;
int NumSpecies;
int MaxAlignIdxSpecies;
int Idx;
etSeqBase *pSeq;

int	 CurRefChromID;
char CurRefStrand;
int	 CurRefChromOfs;
int	 CurRefChromEndOfs;
int  CurRefAlignLen;

int	 RefChromID;
char RefStrand;
int	 RefChromOfs;
int	 RefChromEndOfs;
int  RefAlignLen;

int RelSpeciesID;
int RelChromID;

bFirst = true;
CurBlockID = BlockID;
PrvBlockID = CurBlockID;
pProcParams->NumSpeciesInAlignment = 0;
pProcParams->MaxAlignIdxSpecies = 0;
RefAlignLen=   0;
CurRefChromID = 0;
CurRefStrand = '+';
CurRefChromOfs = -1;
CurRefChromEndOfs = -1;
CurRefAlignLen = 0;

while(CurBlockID >= 0 && ((CurBlockID = pAlignments->NxtBlock(CurBlockID)) > 0))
	{
	CurRefChromID  = pAlignments->GetRelChromID(CurBlockID,RefSpeciesID);
	CurRefStrand   = pAlignments->GetStrand(CurBlockID,RefSpeciesID);
	CurRefChromOfs = pAlignments->GetRelChromOfs(CurBlockID,RefSpeciesID);
	CurRefChromEndOfs = pAlignments->GetRelChromEndOfs(CurBlockID,RefSpeciesID);
	CurRefAlignLen = pAlignments->GetAlignLen(CurBlockID,RefSpeciesID);
	
			// terminate with blocks in which there are less species than the minimum number we need!
	if((CurNumAlignments = pAlignments->GetNumSpecies(CurBlockID)) < pProcParams->NumCoreSpecies)
		break;

	// check if this aligned block would overflow alloc'd memory for holding concatenated alignments
	if((CurRefAlignLen + RefAlignLen) > pProcParams->MaxSeqAlignLen)
		{
		CurBlockID = PrvBlockID; // force re-read of block next time
		break;
		}

	if(bFirst) 
		{
		RefChromID =   CurRefChromID;
		RefStrand  =   CurRefStrand;
		RefChromOfs=   CurRefChromOfs;
		RefChromEndOfs=CurRefChromOfs;
		}
	else
		{
		if(CurRefChromID != RefChromID || 
		   CurRefStrand != RefStrand || 
		   (CurRefChromOfs != (RefChromEndOfs + 1)))
			{
			CurBlockID = PrvBlockID; // force re-read of block next time
			break;
			}
		}

		// iterate over species aligned in current block, species are assumed to be in priority order
	for(NumSpecies = Idx = 0; Idx < pProcParams->NumSpeciesList; Idx++) 
		{
		RelSpeciesID = pSpeciesIDs[Idx];
		RelChromID = pAlignments->GetRelChromID(CurBlockID,RelSpeciesID);

		if(RelChromID < 1)			// assume < 1 is because species not in alignment
			{
			if(Idx < pProcParams->NumCoreSpecies)	// first pParams->NumCoreSpecies must always be present
				break;
			memset(&pProcParams->pSeqs[Idx][RefAlignLen],eBaseUndef,CurRefAlignLen);
			continue;		
			}

		if(ExcludeThisChrom(pAlignments,RelSpeciesID,RelChromID,pProcParams))	// should this chromosome be accepted for processing or sloughed?
			{
			if(Idx < pProcParams->NumCoreSpecies)	// first pParams->NumCoreSpecies must always be present
				break;
			memset(&pProcParams->pSeqs[Idx][RefAlignLen],eBaseUndef,CurRefAlignLen);
			continue;	
			}

			// get each sequence
		if((pSeq = pAlignments->GetSeq(CurBlockID,RelSpeciesID))==NULL) 
			{
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"LoadContiguousBlocks: Unexpected missing sequence in alignments");
			break;
			}
		
		memcpy(&pProcParams->pSeqs[Idx][RefAlignLen],pSeq,CurRefAlignLen);
		NumSpecies++;
		MaxAlignIdxSpecies = Idx + 1;
		}
		
	// check if required minimum number of species were in multialignment
	if(NumSpecies < pProcParams->MinAlignSpecies)	
		break;

	if(NumSpecies > pProcParams->NumSpeciesInAlignment)
		pProcParams->NumSpeciesInAlignment = NumSpecies;
	if(MaxAlignIdxSpecies > pProcParams->MaxAlignIdxSpecies)
		pProcParams->MaxAlignIdxSpecies = MaxAlignIdxSpecies;
	RefAlignLen += CurRefAlignLen;
	RefChromEndOfs=CurRefChromEndOfs;
	PrvBlockID = CurBlockID; 
	bFirst = false;
	}

if(bFirst) 
	{
	*pRefChromID  =   CurRefChromID;
	*pRefStrand   =   CurRefStrand;
	*pRefChromOfs =   CurRefChromOfs;
	*pRefAlignLen =   CurRefAlignLen;
	*pbLoaded = false;
	}
else
	{
	*pRefChromID    =RefChromID;
	*pRefStrand     =RefStrand;
	*pRefChromOfs   =RefChromOfs;
	*pRefAlignLen   =RefAlignLen;
	*pbLoaded = true;
	}
return(CurBlockID);
}



// NormaliseInDelColumns
// Because multialignments may be the result of merged alignments resulting in InDels being generated back into the 
// existing reference sequence then if only a subset of species are being processed there could be all InDels in any column
// of the subset sequences. In addition, where there are optional species then there could be eBaseUndefs in InDel columns
// This function will delete all columns which only contain InDels and/or eBaseUndef
// Returns the subset sequence length after eBaseInDel/eBaseUndef columns have been deleted
int
NormaliseInDelColumns(tsProcParams *pProcParams,int AlignLen)
{
etSeqBase *pSeq;
etSeqBase *pSrcSeq;
etSeqBase SeqBase;
int NumVInDels;
int SeqIdx;
int VIdx;
int FirstInDel=-1;
int NormAlignLen=0;
for(SeqIdx = 0; SeqIdx < AlignLen; SeqIdx++)
	{
	NumVInDels = 0;					// assume no InDels in column
	for(VIdx = 0; VIdx < pProcParams->MaxAlignIdxSpecies; VIdx++)
		{
		pSeq = pProcParams->pSeqs[VIdx];
		SeqBase = pSeq[SeqIdx] & ~cRptMskFlg;
		if(SeqBase == eBaseInDel || SeqBase == eBaseUndef)
			NumVInDels++;
		else
			break;				
		}

	if(NumVInDels == pProcParams->MaxAlignIdxSpecies)	// all were InDels or undefined?
		{
		// mark col for deletion
		for(VIdx = 0; VIdx < pProcParams->MaxAlignIdxSpecies; VIdx++)
			{
			pSeq = pProcParams->pSeqs[VIdx];
			pSeq[SeqIdx] = (etSeqBase)0x0ff;
			}
		if(FirstInDel == -1)		// note idx of first InDel
			FirstInDel = SeqIdx;
		}
	else
		NormAlignLen++;				// accept this column
	}
if(NormAlignLen == AlignLen)	// if no columns to delete
	return(AlignLen);

// have at least one column which is all InDels and is to be deleted
for(VIdx = 0; VIdx < pProcParams->MaxAlignIdxSpecies; VIdx++)
	{
	pSeq = pProcParams->pSeqs[VIdx];
	pSeq = &pSeq[FirstInDel];
	pSrcSeq = pSeq;
	for(SeqIdx = FirstInDel; SeqIdx < AlignLen; SeqIdx++,pSrcSeq++)
		{
		if(*pSrcSeq != (etSeqBase)0x0ff)
			*pSeq++ = *pSrcSeq;
		}
	}
return(NormAlignLen);
}



// ProcAlignBlock
// Process an alignment block which may contain aligned subsequences meeting processing requirements
// Any subsequence (bounded by InDels if bInDelsAsMismatches, or end of block) of longer than MinLen is passed on
// to ProcessAlignment() which will process that subsequence looking for hyperconserved cores 
bool 
ProcAlignBlock(int RefChromID,	    // reference chromosome
			   int RefChromOfs,		// offset along reference chromosome 
			   int AlignLen,		// alignment length incl any InDels
			   tsProcParams *pProcParams) // global processing parameters
{
int CurSeqLen;
int CurFiltSeqLen;	
int SubRefOfs;
int SubSeqStartIdx;
int NumIdents;
int SeqIdx;
int NumSubSeqs;
int CurNumSeqs;
int VIdx;
bool bAllIdentical; 
etSeqBase *pSeq;
etSeqBase RefBase;
etSeqBase SeqBase;
int NumVInDels;
bool bRefInDel;


if(pProcParams->bFiltLoConfidence && pProcParams->bFilt1stLast)
	NumSubSeqs = CUtility::GetNumSubseqs(AlignLen,		// alignment length incl InDels
							pProcParams->MinAlignSpecies,
							pProcParams->pSeqs);
else
	NumSubSeqs = 0;
CurNumSeqs=0;
CurSeqLen = 0;
CurFiltSeqLen = 0;
SubRefOfs = RefChromOfs;
SubSeqStartIdx = 0;
//if(RefChromID == 1 && RefChromOfs > 10673120 && RefChromOfs < 10690000)
//	printf("\nAt ChromOfs: %d",RefChromOfs);

for(SeqIdx = 0; SeqIdx < AlignLen; SeqIdx++)
	{
	bAllIdentical = true;			// assume all identical in column
	NumVInDels = 0;					// assume no InDel in column
	bRefInDel = false;				// and specifically that the reference base is not an InDel 
	for(VIdx = 0; VIdx < pProcParams->MinAlignSpecies; VIdx++)
		{
		pSeq = pProcParams->pSeqs[VIdx];
		SeqBase = pSeq[SeqIdx] & ~cRptMskFlg;
		if(VIdx == pProcParams->RefSpeciesIdx)
			RefBase = SeqBase;
		if(SeqBase == eBaseInDel)
			{
			NumVInDels++;
			bAllIdentical = false;
			if(VIdx == pProcParams->RefSpeciesIdx)
				bRefInDel = true;
			}
		else
			if(SeqBase == eBaseN || SeqBase != RefBase)
				bAllIdentical = false;
		}
	// if all bases in column were InDels then alert user as this is a dataset processing error
	// should have been filtered out before this function is called..
	if(NumVInDels == pProcParams->MinAlignSpecies && pProcParams->MinAlignSpecies == pProcParams->NumSpeciesInAlignment)
		gDiagnostics.DiagOut(eDLWarn,gszProcName,"Unexpected, all aligned bases were InDels at reference offset %d",RefChromOfs);

	// if in column any base was eBaseInDel then NumVInDels will be > 0
	// if the reference base was an InDel then bRefInDel will be true
	// if in column all identical and a,c,g or t then bAllIdentical will be true
	// if any base was mismatched or not a,c,g or t then bAllIdentical will be false
	if((pProcParams->ProcMode != eProcModeOutspecies && NumVInDels) && (!pProcParams->bInDelsAsMismatches &&
		(!pProcParams->bSloughRefInDels || (pProcParams->bSloughRefInDels && !bRefInDel))))
		{
		if(CurFiltSeqLen >= pProcParams->MinCoreLen)	// has to be of at least MinLen to be worth further processing
			{
			CurNumSeqs+=1;
			if(CurFiltSeqLen >= pProcParams->MinSubSeqLen)			
				{
				if(pProcParams->bFiltLoConfidence)					// when filtering simply slough any low confidence subsequences
					{
					if(pProcParams->bFilt1stLast && (CurNumSeqs==1 || CurNumSeqs == NumSubSeqs))	// don't process first and last subsequence
						CurFiltSeqLen = 0;
					if(CurFiltSeqLen &&  pProcParams->MinIdent > 0 && ((NumIdents * 100)/CurFiltSeqLen) < pProcParams->MinIdent)
						CurFiltSeqLen = 0;
					}
				if(CurFiltSeqLen >= pProcParams->MinCoreLen)
					ProcessAlignment(RefChromID,SubRefOfs,SubSeqStartIdx,CurFiltSeqLen,pProcParams);
				}
			}
		CurSeqLen = 0;
		CurFiltSeqLen = 0;
		if(RefBase != eBaseInDel)
			RefChromOfs++;
		NumIdents = 0;
		SubSeqStartIdx = SeqIdx + 1;	// mark where next subsequence could start 
		SubRefOfs = RefChromOfs;		// chromosomal offset at which that next subsequence starts
		continue;
		}

	// bInDelsAsMismatches or no InDels in current aligned column - could still be mismatches - ,
	// continue to accept column as part of subsequence
	if(CurSeqLen == 0)		// if first base of putative subsequence...
		{
		if(pProcParams->bFiltLoConfidence && !bAllIdentical) // when filtering, must start on an identical base
			{
			if(RefBase != eBaseInDel)		// advance chromosomal offset  only if not a InDel
				RefChromOfs++;
			SubRefOfs = RefChromOfs;		// chromosomal offset at which that next subsequence starts
			continue;
			}
		NumIdents = 0;
		SubSeqStartIdx = SeqIdx;			// mark where this subsequence has started 
		SubRefOfs = RefChromOfs;			// chromosomal offset at which this subsequence started
		}
	CurSeqLen++;
	if(bAllIdentical)						// mark last identical...
		NumIdents++;
	CurFiltSeqLen = CurSeqLen;
	if(RefBase != eBaseInDel)				// advance chromosomal offset  only if not a InDel
		RefChromOfs++;
	}

// no more bases in this block
if(CurFiltSeqLen >= pProcParams->MinCoreLen)	// has to be at least MinLen to be worth further processing
	{
	if(CurFiltSeqLen >= pProcParams->MinSubSeqLen)			
		{
		if(pProcParams->bFiltLoConfidence) // when filtering simply slough any low confidence subsequences
			{
			if(pProcParams->bFilt1stLast && (!CurNumSeqs || CurNumSeqs == (NumSubSeqs-1)))	// don't process first and last subsequence
				CurFiltSeqLen = 0;
			if(CurFiltSeqLen &&  pProcParams->MinIdent > 0 && ((NumIdents * 100)/CurFiltSeqLen) < pProcParams->MinIdent)
				CurFiltSeqLen = 0;
			}
		if(CurFiltSeqLen >= pProcParams->MinCoreLen)
			ProcessAlignment(RefChromID,SubRefOfs,SubSeqStartIdx,CurFiltSeqLen,pProcParams);
		}
	}

return(true);
}

// ProcAlignBlockSummary
// Process an alignment block which may contain aligned subsequences meeting processing requirements
// Generates stats on the sequences within this alignment block
// Stats: Total number of mismatches
//		  Total number of exact matches
//		  Total numer of ref InDels
//		  Total number of rel InDels
// Above is characterised by region
bool 
ProcAlignBlockSummary(int RefChromID,	    // reference chromosome
			   int RefChromOfs,		// offset along reference chromosome 
			   int AlignLen,		// alignment length incl any InDels
			   tsProcParams *pProcParams) // global processing parameters
{
int CurSeqLen;
int CurFiltSeqLen;	
int SubRefOfs;
int SubSeqStartIdx;
int SeqIdx;
int NumSubSeqs;
int CurNumSeqs;
int VIdx;

etSeqBase *pSeq;
etSeqBase RefBase;
etSeqBase RelBase;
etSeqBase SeqBase;
int NumInDels;
int NumMatching;

if(pProcParams->bFiltLoConfidence && pProcParams->bFilt1stLast)
	NumSubSeqs = CUtility::GetNumSubseqs(AlignLen,		// alignment length incl InDels
							pProcParams->MinAlignSpecies,
							pProcParams->pSeqs);
else
	NumSubSeqs = 0;
CurNumSeqs=0;
CurSeqLen = 0;
CurFiltSeqLen = 0;
SubRefOfs = RefChromOfs;
SubSeqStartIdx = 0;

for(SeqIdx = 0; SeqIdx < AlignLen; SeqIdx++)
	{
	NumInDels = 0;
	NumMatching = 1;	// 1st base always matches itself!
	for(VIdx = 0; VIdx < pProcParams->MinAlignSpecies; VIdx++)
		{
		pSeq = pProcParams->pSeqs[VIdx];
		SeqBase = pSeq[SeqIdx] & ~cRptMskFlg;
		if(VIdx == pProcParams->RefSpeciesIdx)
			RefBase = SeqBase;

		if(SeqBase == eBaseInDel)
			NumInDels++;
		if(!VIdx)
			RelBase = SeqBase;
		else
			if(RelBase == SeqBase)
				NumMatching++;
		}

	// if all bases in column were InDels then alert user as this is a dataset processing error
	// should have been filtered out before this function is called..
	if(NumInDels ==  pProcParams->MinAlignSpecies && pProcParams->MinAlignSpecies == pProcParams->NumSpeciesInAlignment)
		{
		gDiagnostics.DiagOut(eDLWarn,gszProcName,"Unexpected, all aligned bases were InDels at reference offset %d",RefChromOfs);
		return(false);
		}
	if(RefBase == eBaseInDel)
		{
		ReportSummary(RefChromID,RefChromOfs,0,pProcParams);
		continue;	// don't advance chromofs
		}
	if(NumInDels)	// ref wasn't InDel so if any InDels then they must be in a relative species
		ReportSummary(RefChromID,RefChromOfs,1,pProcParams);
	else
		{
		if((RefBase <= eBaseT) && NumMatching == pProcParams->MinAlignSpecies)		
			ReportSummary(RefChromID,RefChromOfs,2,pProcParams); // match
		else						
			ReportSummary(RefChromID,RefChromOfs,3,pProcParams); // at least one mismatch
		}
	RefChromOfs++;
	}
return(true);
}

int
ReportSummary(int RefChromID,int RefChromOfs,int ProcMode,tsProcParams *pProcParams)
{
int FeatureBits;
int RegionIdx;
int BitMsk;
int SpliceSiteOverlaps;
int FeatIdx;
int *pStep = pProcParams->pCntStepCnts;

// ensure that the hypercore is in an included region and not part of an excluded region
if(!IncludeFilter(RefChromID,RefChromOfs,RefChromOfs,pProcParams))
		return(eBSFSuccess);

pStep += ProcMode * pProcParams->Regions;
	
if(pProcParams->pBiobed != NULL)
	{
	if(pProcParams->BEDChromID > 0)
		FeatureBits = pProcParams->pBiobed->GetFeatureBits(pProcParams->BEDChromID,RefChromOfs,RefChromOfs,cRegionFeatBits,pProcParams->UpDnStreamLen);
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
				if(RegionIdx)		// if already have feature
					return(eBSFSuccess);	// although was sequence of interest, more than one feature bit so can't contribute to stats
				RegionIdx = FeatIdx;
				if(pProcParams->bMultipleFeatBits)	// if multiple features allowed then don't check for any additional
					break;
				}
			}
		}

	// need to remap RegionIdx when incr counts so regions are displayed in a 5'->3' logical order
	switch(RegionIdx) {
		case 0: case 2: case 4: case 6:		// IG,5'UTR, Introns and 3'DS
			pStep[RegionIdx] += 1;	
			break;
		case 1:								// CDS
			pStep[3] += 1;
			break;
		case 3:								// 3'UTR
			pStep[5] += 1;
			break;
		case 5:								// 5'US
			pStep[1] += 1;
			break;
		}
			
	if(RegionIdx != 0)		// if not intergenic then check for splice sites
		{
		SpliceSiteOverlaps = pProcParams->pBiobed->GetSpliceSiteBits(pProcParams->BEDChromID,RefChromOfs,RefChromOfs,cMinSpliceOverlap);
		if(SpliceSiteOverlaps & cIntronExonSpliceSite)
			pStep[7]++;
		if(SpliceSiteOverlaps & cExonIntronSpliceSite)
			pStep[8]++;
		}
	pProcParams->bStatsAvail = true;
	ChkOutputSummaryResults(pProcParams->szRefChrom,RefChromOfs,pProcParams,false,false);	// output results as may be appropriate
	}
else
	{
	*pStep += 1;
	ChkOutputSummaryResults(pProcParams->szRefChrom,RefChromOfs,pProcParams,false,false);	// output results as may be appropriate
	}
return(eBSFSuccess);
}



// OpenBedfile
// Attempts to open specified bedfile
// Returns ptr to opened bedfile or NULL
CBEDfile *
OpenBedfile(char *pToOpen)
{
int Rslt;
CBEDfile *pBed;
if(pToOpen != NULL && pToOpen[0] != '\0')
	{
	if((pBed = (CBEDfile *)new CBEDfile())==NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CBEDfile '%s'",pToOpen);
		return(NULL);
		}

	if((Rslt = pBed->Open(pToOpen,eBTAnyBed))!=eBSFSuccess)
		{
		while(pBed->NumErrMsgs())
			gDiagnostics.DiagOut(eDLFatal,gszProcName,pBed->GetErrMsg());
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open biobed file %s",pToOpen);
		delete pBed;
		return(NULL);
		}
	return(pBed);
	}
return(NULL);
}

//CloseBedfiles
//Closes and deletes all created and opened Biobed files
bool
CloseBedfiles(tsProcParams *pProcParams)
{
int Idx;
for(Idx=0; Idx < pProcParams->NumIncludes; Idx++)
	{
	if(pProcParams->pIncludes[Idx] != NULL)
		delete pProcParams->pIncludes[Idx];
	pProcParams->pIncludes[Idx] = NULL;
	}
for(Idx=0; Idx < pProcParams->NumExcludes; Idx++)
	{
	if(pProcParams->pExcludes[Idx] != NULL)
		delete pProcParams->pExcludes[Idx];
	pProcParams->pExcludes[Idx] = NULL;
	}
if(pProcParams->pBiobed != NULL)
	delete pProcParams->pBiobed;
pProcParams->pBiobed = NULL;
return(true);
}

// IncludeFilter
// Returns true if subsequence can be processed or false if subsequence is not to be processed
// To be processed a subsequence -
// A) If no Include BED files specified then subsequence is assumed included unless excluded by following rule B).
//    If at least one Include BED file was specified then if subsequence not in any of the Include files then false is returned
// B) If no Exclude files specified then subsequence is returned as Ok (true) to process. If subsequence not in any of the Exclude files
//    then subsequence is returned as Ok (true) to process, otherwise false is returned
bool
IncludeFilter(int RefChromID,int SubRefOfs,int SubRefEndOfs,tsProcParams *pProcParams)
{
int Idx;
int BEDChromID;
if(pProcParams->NumIncludes)
	{
	for(Idx = 0; Idx < pProcParams->NumIncludes; Idx++)
		{
		if(pProcParams->pIncludes[Idx] == NULL) // should'nt ever be NULL but...
			continue;
		if((BEDChromID = pProcParams->pIncludes[Idx]->LocateChromIDbyName(pProcParams->szRefChrom))<1)
			continue;
		if(pProcParams->pIncludes[Idx]->InAnyFeature(BEDChromID,SubRefOfs,SubRefEndOfs))
			break;
		}
	if(Idx == pProcParams->NumIncludes)
		return(false);
	}

if(pProcParams->NumExcludes)
	{
	for(Idx = 0; Idx < pProcParams->NumExcludes; Idx++)
		{
		if(pProcParams->pExcludes[Idx] == NULL) // should'nt ever be NULL but...
			continue;
		if((BEDChromID = pProcParams->pExcludes[Idx]->LocateChromIDbyName(pProcParams->szRefChrom))<1)
			continue;
		if(pProcParams->pExcludes[Idx]->InAnyFeature(BEDChromID,SubRefOfs,SubRefEndOfs))
			return(false);
		}
	}
return(true);
}

// Process
// Create alignment stats from files in specified source directory
int
Process(bool bTargDeps,				// true if process only if any independent src files newer than target
		int ProcMode,				// processing mode 0: default, 1: summary stats only, 2: outspecies
				 char *pszInputFile,		// bio multialignment (.algn) file to process
					char *pszOutputFile,	// where to write out stats
					char *pszOutputCoreFile, // where to write out the hypercore loci 
					char *pszBiobedFile,	// biobed file containing regional features - exons, introns etc
					int NumIncludeFiles,	// number of include region files
					char **ppszIncludeFiles,	// biobed files containing regions to include - default is to exclude none
					int NumExcludeFiles,	// number of exclude region files
					char **ppszExcludeFiles,	// biobed file containing regions to exclude - default is to include all
 					char *pszUniqueSpeciesSeqs, // ignore alignment block if these species sequences are not unique
					int WindowSize,			// sampling window size
					int NumCoreSpecies,		// number of core species to be in alignment
					int MinAlignSpecies,	// minimum number of species required in an alignment (includes core species)
					char *pszSpeciesList,	// space or comma separated list of species, priority ordered
					int MinHyperLen,		// minimum hyper core length required
					int MinUltraLen,		// minimum ultra core length required
					int MaxHyperMismatches,	// hyper cores can have upto this number of total mismatches
					int VMismatches,		// number of mismatches in an alignment col to count as a missmatch against MaxMismatches
					int MinIdentity,		// minimum identity required when processing hyperconserved
					int NumDistSegs,		// number of match distribution profile segments
					int RegLen,				// regulatory region length - up/dn stream of 5/3' 
					bool bMultipleFeatBits,	// if true then accept alignments in which multiple feature bits are set
					bool bInDelsAsMismatches, // treat InDels as if mismatches
					bool bSloughRefInDels,			// slough columns in which the reference base is an InDel
					bool bFiltLoConfidence,		// filter out low confidence subsequences
					bool bFilt1stLast,			// treat 1st and last subsequences as being low confidence
					int MinSubSeqLen,			// subsequences of less than this length are treated as being low confidence
					int MinIdent,				// treat subsequences of less than this identity as being low confidence
					int NumIncludeChroms,	// number of chromosomes explicitly defined to be included
					char **ppszIncludeChroms,	// ptr to array of reg expressions defining chroms to include - overides exclude
					int NumExcludeChroms,	// number of chromosomes explicitly defined to be excluded
					char **ppszExcludeChroms)	// ptr to array of reg expressions defining chroms to include
{
int Rslt;
int Idx;
int *pCntStepCnts;
int NumCnts;
int Regions;
int UpDnStreamLen;
CBEDfile *pBiobed = NULL;
tsProcParams ProcParams;
int	MismatchScore;
int	MatchScore;
char szCSVSpecies[2048];				// to hold comma separated species list
char szInaccessible[_MAX_PATH];

if(bTargDeps && (Rslt = CUtility::Chk2TargDepend(szInaccessible,_MAX_PATH,pszOutputFile,pszOutputCoreFile,pszBiobedFile,NULL)) <= 0)
	{
	if(Rslt)
		{
		gDiagnostics.DiagOutMsgOnly(eDLFatal,"Unable to access source file '%s'",szInaccessible);
		return(Rslt);
		}
	for(Idx=0;Idx<NumIncludeFiles; Idx++)
		{
		Rslt = CUtility::Chk2TargDepend(szInaccessible,_MAX_PATH,pszOutputFile,pszOutputCoreFile,ppszIncludeFiles[Idx],NULL);
		if(Rslt < 0)
			{
			gDiagnostics.DiagOutMsgOnly(eDLFatal,"Unable to access source include file '%s'",szInaccessible);
			return(Rslt);
			}
		}

	for(Idx=0;Idx<NumExcludeFiles; Idx++)
		{
		Rslt = CUtility::Chk2TargDepend(szInaccessible,_MAX_PATH,pszOutputFile,pszOutputCoreFile,ppszExcludeFiles[Idx],NULL);
		if(Rslt < 0)
			{
			gDiagnostics.DiagOutMsgOnly(eDLFatal,"Unable to access source exclude file '%s'",szInaccessible);
			return(Rslt);
			}
		}


	}

#ifndef _WIN32
int RegErr;				// regular expression parsing error
char szRegErr[128];			// to hold RegErr as textual representation ==> regerror();
#endif

memset(&ProcParams,0,sizeof(tsProcParams));

// parse out species which must have unique alignment block sequences
if(!stricmp(pszUniqueSpeciesSeqs,"*"))
	ParseUniqueSpeciesSeqs(pszSpeciesList,&ProcParams);
else
	ParseUniqueSpeciesSeqs(pszUniqueSpeciesSeqs,&ProcParams);

// parse out species list
ProcParams.NumSpeciesList = ParseNumSpecies(pszSpeciesList,&ProcParams);
szCSVSpecies[0]='\0';
for(Idx = 0; Idx < ProcParams.NumSpeciesList; Idx++)
	{
	if(Idx > 0)
		strcat(szCSVSpecies,",");
	strcat(szCSVSpecies,ProcParams.szSpecies[Idx]);
	}
ProcParams.pszSpeciesList = szCSVSpecies;

for(Idx=0;Idx<NumIncludeFiles; Idx++)
	{
	if((ProcParams.pIncludes[Idx] = OpenBedfile(ppszIncludeFiles[Idx]))==NULL)
		{
		CloseBedfiles(&ProcParams);
		return(eBSFerrObj);
		}
	ProcParams.NumIncludes++;
	}

for(Idx=0;Idx<NumExcludeFiles; Idx++)
	{
	if((ProcParams.pExcludes[Idx] = OpenBedfile(ppszExcludeFiles[Idx]))==NULL)
		{
		CloseBedfiles(&ProcParams);
		return(eBSFerrObj);
		}
	ProcParams.NumExcludes++;
	}


if(pszBiobedFile != NULL && pszBiobedFile[0] != '\0')
	{
	if((ProcParams.pBiobed = OpenBedfile(pszBiobedFile))==NULL)
		{
		CloseBedfiles(&ProcParams);
		return(eBSFerrObj);
		}
	Regions = cNumCharRegs;	
	UpDnStreamLen = RegLen;
	}
else
	{
	Regions = 1;
	UpDnStreamLen = 0;
	ProcParams.pBiobed = NULL;
	}

if(ProcMode == eProcModeSummary)
	NumCnts=Regions * 4;	// allows for RefIndels,RelInDels,Matches and Mismatch counts per region
else
	NumCnts=Regions * cLenRanges;	
if((pCntStepCnts = new int[NumCnts])==NULL)
	{
	CloseBedfiles(&ProcParams);
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"nUnable to allocate memory (%d bytes) for holding alignment statistics",sizeof(int) * NumCnts);
	return(eBSFerrMem);
	}
memset(pCntStepCnts,0,NumCnts * sizeof(int));

ProcParams.bInDelsAsMismatches = bInDelsAsMismatches;
ProcParams.bSloughRefInDels = bSloughRefInDels;			// slough columns in which the reference base is an InDel

ProcParams.bFiltLoConfidence = bFiltLoConfidence;
ProcParams.bFilt1stLast = bFilt1stLast;
ProcParams.MinIdent = MinIdent;
ProcParams.MinSubSeqLen = MinSubSeqLen;
ProcParams.ProcMode = ProcMode;
ProcParams.MinAlignSpecies = MinAlignSpecies;
ProcParams.NumCoreSpecies = NumCoreSpecies;
ProcParams.NumSpeciesInAlignment = 0;
ProcParams.MaxAlignIdxSpecies = 0;
ProcParams.RefSpeciesIdx = 0;
ProcParams.pCntStepCnts = pCntStepCnts;
ProcParams.NumDistSegs = NumDistSegs;
ProcParams.UpDnStreamLen = UpDnStreamLen;
ProcParams.bMultipleFeatBits = bMultipleFeatBits;
ProcParams.Regions = Regions;
ProcParams.NumCnts = NumCnts;
ProcParams.MinUltraLen = MinUltraLen;
ProcParams.MinHyperLen = MinHyperLen;
#ifdef _WIN32
ProcParams.MinCoreLen = max(MinUltraLen,MinHyperLen);
#else
ProcParams.MinCoreLen = MinUltraLen > MinHyperLen ? MinUltraLen : MinHyperLen;
#endif

ProcParams.MaxHyperMismatches = MaxHyperMismatches;
ProcParams.MinIdentity = MinIdentity;		// minimum identity required when processing hyperconserved

if(ProcMode != eProcModeSummary && MinIdentity < 100)
	{
	MismatchScore = (cRandWalk100Score - 1) / (100 - MinIdentity);	   // decrease score by this much on each mismatch
	MatchScore = cRandWalk100Score / MinIdentity;			   // increase score by this much each time there is a match
	}
else
	{
	MismatchScore = 0;		// 100% identity is an ultra!
	MatchScore = 0;
	}

ProcParams.MismatchScore = MismatchScore;
ProcParams.MatchScore = MatchScore;
ProcParams.VMismatches = VMismatches;
ProcParams.NumIncludeChroms = NumIncludeChroms;		// number of chromosomes explicitly defined to be included
ProcParams.ppszIncludeChroms = ppszIncludeChroms;	// ptr to array of reg expressions defining chroms to include - overides exclude
ProcParams.NumExcludeChroms = NumExcludeChroms;		// number of chromosomes explicitly defined to be excluded
ProcParams.ppszExcludeChroms = ppszExcludeChroms;	// ptr to array of reg expressions defining chroms to include

#ifdef _WIN32
try {
	for(Idx=0;Idx < NumIncludeChroms;Idx++)
		{
		ProcParams.IncludeChromsRE[Idx] = new Regexp();
		ProcParams.IncludeChromsRE[Idx]->Parse(ppszIncludeChroms[Idx],false);	// note case insensitive
		}
	}
catch(...)
	{
	delete pCntStepCnts;
	CloseBedfiles(&ProcParams);
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to process include regexpr chrom '%s'",ppszIncludeChroms[Idx]);
	return(eBSFerrMem);
	}
try {
	for(Idx=0;Idx < NumExcludeChroms;Idx++)
		{
		ProcParams.ExcludeChromsRE[Idx] = new Regexp();
		ProcParams.ExcludeChromsRE[Idx]->Parse(ppszExcludeChroms[Idx],false);	// note case insensitive
		}
	}
catch(...)
	{
	delete pCntStepCnts;
	CloseBedfiles(&ProcParams);
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to process exclude regexpr chrom '%s'",ppszExcludeChroms[Idx]);
	return(eBSFerrMem);
	}

#else
for(Idx=0;Idx < NumIncludeChroms;Idx++)
	{

	RegErr=regcomp(&ProcParams.IncludeChromsRE[Idx],ppszIncludeChroms[Idx],REG_EXTENDED | REG_ICASE);	// note case insensitive
	if(RegErr)
		{
		regerror(RegErr,&ProcParams.IncludeChromsRE[Idx],szRegErr,sizeof(szRegErr));
		delete pCntStepCnts;
		CloseBedfiles(&ProcParams);
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to process include chrom '%s' error: %s",ppszIncludeChroms[Idx],szRegErr);
		return(eBSFerrMem);
		}
	}
for(Idx=0;Idx < NumExcludeChroms;Idx++)
	{
	RegErr = regcomp(&ProcParams.ExcludeChromsRE[Idx],ppszExcludeChroms[Idx],REG_EXTENDED | REG_ICASE);	// note case insensitive
	if(RegErr)
		{
		regerror(RegErr,&ProcParams.ExcludeChromsRE[Idx],szRegErr,sizeof(szRegErr));
		delete pCntStepCnts;
		CloseBedfiles(&ProcParams);
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to process exclude chrom '%s' error: %s",ppszExcludeChroms[Idx],szRegErr);
		return(eBSFerrMem);
		}


	}
#endif

// determine max aligned sequence length for any single species which can be handled
ProcParams.MaxSeqAlignLen = cMAtotMaxSeqAlignLen/ProcParams.NumSpeciesList;
for(Idx = 0; Idx < ProcParams.NumSpeciesList; Idx++)
	{
	if((ProcParams.pSeqs[Idx] = new unsigned char [ProcParams.MaxSeqAlignLen])==NULL)
		{
		delete pCntStepCnts;
	    CloseBedfiles(&ProcParams);
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate memory (%d bytes) for holding species sequences",ProcParams.MaxSeqAlignLen);
		return(eBSFerrMem);
		}
	}

#ifdef _WIN32
if((ProcParams.hRsltsFile = open(pszOutputFile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
if((ProcParams.hRsltsFile = open(pszOutputFile,O_RDWR | O_CREAT |O_TRUNC, S_IREAD | S_IWRITE ))==-1)
#endif
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create %s - %s",pszOutputFile,strerror(errno));
	delete pCntStepCnts;
    CloseBedfiles(&ProcParams);
	return(eBSFerrCreateFile);
	}

if(pszOutputCoreFile != NULL && pszOutputCoreFile[0] != '\0')
	{
#ifdef _WIN32
	if((ProcParams.hCoreCSVRsltsFile = open(pszOutputCoreFile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
	if((ProcParams.hCoreCSVRsltsFile = open(pszOutputCoreFile, O_RDWR | O_CREAT | O_TRUNC, S_IREAD | S_IWRITE ))==-1)
#endif
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create %s - %s",pszOutputCoreFile,strerror(errno));
		close(ProcParams.hRsltsFile);
		delete pCntStepCnts;
		CloseBedfiles(&ProcParams);
		return(eBSFerrCreateFile);
		}
	}
else
	ProcParams.hCoreCSVRsltsFile = -1;

Rslt = ProcessAlignments(pszInputFile,&ProcParams);
if(Rslt == eBSFSuccess && !WindowSize)
	{
	if(ProcMode == eProcModeSummary)
		OutputSummaryResults( "Alignment",0,&ProcParams,false);
	else
		OutputResults( "Genome",0,&ProcParams,false);
	}
if(pCntStepCnts != NULL)
	delete(pCntStepCnts);
if(ProcParams.hRsltsFile != -1)
	close(ProcParams.hRsltsFile);
if(ProcParams.hCoreCSVRsltsFile != -1)
	close(ProcParams.hCoreCSVRsltsFile);
CloseBedfiles(&ProcParams);
for(Idx = 0; Idx < ProcParams.NumSpeciesList; Idx++)
	if(ProcParams.pSeqs[Idx] != NULL)
		delete ProcParams.pSeqs[Idx];

#ifdef _WIN32
for(Idx=0;Idx < NumIncludeChroms;Idx++)
	if(ProcParams.IncludeChromsRE[Idx] != NULL)
		delete ProcParams.IncludeChromsRE[Idx];
for(Idx=0;Idx < NumExcludeChroms;Idx++)
	if(ProcParams.ExcludeChromsRE[Idx] != NULL)
		delete ProcParams.ExcludeChromsRE[Idx];
#endif

return(Rslt);
}



// process alignment
// Process an alignment subsequence which may contain mismatches or 'N's but will only contain InDels if bInDelsAsMismatches is true 
// The subsequence length will be of at least MinLen long
bool
ProcessAlignment(int RefChromID,		// reference chromosome
				 int ChromOfs,			// chromosome offset (0..n) at which this subsequence started
				 int SeqIdx,			// index (0..n) into pProcParms->pSeqs[] which alignment subsequence starts	
 				 int SubSeqLen,			// subsequence length
				 tsProcParams *pProcParams)
{
int NxtSeqIdx;
etSeqBase *pRefBase;
pRefBase = pProcParams->pSeqs[pProcParams->RefSpeciesIdx];
pRefBase += SeqIdx;
while(SubSeqLen >= pProcParams->MinCoreLen)
	{
		// from current SubRefOfs maximally extend to right allowing for MaxMismatches
//	if(RefChromID == 1 && (ChromOfs > 10673740 && ChromOfs < 10674000))
//		printf("\nNow at ChromOfs: %d",ChromOfs);
	NxtSeqIdx = ProcessSubSeq(RefChromID,ChromOfs,SeqIdx,SubSeqLen,pProcParams);
	if(NxtSeqIdx == -1)		// -1 flags that there is no more processing required for this alignment
		break;
	SubSeqLen -= (NxtSeqIdx - SeqIdx);	// update length of remaining sequence in this alignment
	if(SubSeqLen < pProcParams->MinCoreLen)
		break;

	// determine chrom offset at which next
	while(SeqIdx < NxtSeqIdx)
		{
		if((*pRefBase++ & ~cRptMskFlg) != eBaseInDel)
			ChromOfs++;
		SeqIdx++;
		}
	}
return(true);
}

// Process subsequence which starts immediately following the end of the previous subsequence of start of an alignment block
// and which ends immediately prior to the next InDel - if bInDelsAsMismatches false -  or end of alignment block
// Returns the SubSeqStartIdx at which to start a subsequent ProcessSubSeq processing call
int	
ProcessSubSeq(int RefChromID,		// reference chromosome identifier
			  int ChromOfs,			// chromosome offset (0..n) at which this subsequence starts
			  int SeqIdx,			// index into pProcParams->pSeq[] at which this subsequence starts
			  int MaxLen,			// maximum length of this subsequence
			  tsProcParams *pProcParams)
{
tsLenRangeClass *pRange;
int VIdx;
int *pStep;
etSeqBase *pVSeq;
etSeqBase VBase;
etSeqBase RefBase;
int RegionIdx;
int FeatureBits;
int BitMsk;
int NxtSeqIdx =-1;
int FeatIdx;
int TotNumMismatches = 0;
int ChromOfsEnd = ChromOfs;
int HypercoreLen = 0;			// current reference species core length (includes any refseq InDels into length)
int VMismatches;
int VIndels;
int SpliceSiteOverlaps;
int IdentityScore = cRandWalk100Score;
int CurUltraCoreLen = 0;		// current core with no mismatches
int MaxUltraCoreLen = 0;		// longest core with no mismatches encountered
int RefHyperCoreLen = 0;		// curent reference species core length (excludes any refseq InDels from length)
int WinIdx = 0;
int NumCoreSpecies;
tsDistSeg OGDistProfile[cMaxMatchDistSegments];

etSeqBase OGBase;
int OGMatches;
int OGMismatches;
int OGInDels;
int OGUnaligned;

OGMatches = 0;
OGMismatches = 0;
OGInDels = 0;
OGUnaligned = 0;

if(pProcParams->ProcMode == eProcModeOutspecies)
	NumCoreSpecies = pProcParams->MinAlignSpecies;
else
	NumCoreSpecies = pProcParams->NumSpeciesInAlignment;

for(HypercoreLen = 0; HypercoreLen < MaxLen ; HypercoreLen++, WinIdx++)
	{
	// determine in the current alignment column the number of mismatches
	VMismatches = 0;
	VIndels = 0;

	for(VIdx = 0; VIdx < NumCoreSpecies; VIdx++)
		{
		pVSeq = pProcParams->pSeqs[VIdx];
		VBase = pVSeq[HypercoreLen+SeqIdx] & ~cRptMskFlg;
		
		if(VIdx == pProcParams->RefSpeciesIdx) // if base for reference sequence
			{
			RefBase = VBase;
			// slough columns in which the ref base is an InDel?
			if(RefBase == eBaseInDel)
				{
				if(pProcParams->bSloughRefInDels)	
					break;
				}
			if(VBase == eBaseInDel)
				VIndels++;
			if(VBase == eBaseN || VBase == eBaseInDel)	  // indeterminate ref seq bases and InDels are 
				VMismatches++;							  // treated as if mismatch
			continue;
			}

		// not reference base so can check for missmatch against relative species base
		if(VBase == eBaseInDel)
				VIndels++;
		if(VBase == eBaseInDel || VBase == eBaseN || // treat InDels and indeterminate bases as if mismatches
		    RefBase != VBase)			   // any bases which don't match the reference base are mismatches
			VMismatches++;
		}

		// processed column of bases
		// should slough columns in which the ref base is an InDel?
	if(RefBase == eBaseInDel && pProcParams->bSloughRefInDels)	
		continue;

		// if processing outspecies mode and core column consisted of InDels then
		// incr InDel count if there is an alignment onto the out species
	if(VIndels == NumCoreSpecies && pProcParams->ProcMode == eProcModeOutspecies)
		{
		if(pProcParams->NumSpeciesInAlignment > NumCoreSpecies)
			OGInDels += 1;
		continue;
		}

	// number of mismatches in column now known - if more than allowed then characterise as a hyperconserved mismatch
	if(VMismatches >= pProcParams->VMismatches)
		{
		CurUltraCoreLen = 0;						// any mismatch terminates current ultracore
		if(pProcParams->MismatchScore)				// if scoring mismatches for hypers, will be 0 for ultras... 
			{
			IdentityScore -= pProcParams->MismatchScore;
			// check if the identity has dropped below minimum required
			if(IdentityScore <= 0)
				break;
			}
		if(NxtSeqIdx == -1)						// if 1st mismatch then note where next search for hypercore should start from
			NxtSeqIdx = HypercoreLen+SeqIdx+1;  // if current hypercore not accepted
				// this many total mismatches allowed?
		if(++TotNumMismatches > pProcParams->MaxHyperMismatches)
			break;
		}
	else	// no mismatch
		{
		CurUltraCoreLen++;							// current ultracore increases as bases identical
		if(CurUltraCoreLen > MaxUltraCoreLen)		// is this the longest ultracore in current subsequence?
			MaxUltraCoreLen = CurUltraCoreLen;
		if(pProcParams->MinUltraLen && CurUltraCoreLen >= pProcParams->MinUltraLen) // when it's an ultra of minimal length
			IdentityScore = cRandWalk100Score;			// then assume backup to 100% identity
		else										// else current ultra not at least minimal length
			{
			IdentityScore += pProcParams->MatchScore;	// increment by random walk match score
			if(IdentityScore > cRandWalk100Score)		// can't do any better than 100%...
				IdentityScore = cRandWalk100Score;
			}
		}

	if(pProcParams->ProcMode == eProcModeOutspecies && pProcParams->NumSpeciesInAlignment > pProcParams->MinAlignSpecies)
		{
		pVSeq = pProcParams->pSeqs[pProcParams->NumSpeciesInAlignment-1];
		OGBase = pVSeq[HypercoreLen+SeqIdx] & ~cRptMskFlg;
		if(OGBase == eBaseUndef)
			OGUnaligned += 1;
		else
			{
			if(RefBase == OGBase && RefBase != eBaseInDel)
				OGMatches+= 1;
			else
				if(RefBase == eBaseInDel || OGBase == eBaseInDel)
					OGInDels += 1;
				else
					OGMismatches += 1;
			}
		}

	if(RefBase != eBaseInDel)
		RefHyperCoreLen++;
	}

// have a sequence which we may be interested in
if(HypercoreLen==MaxLen)// if core was terminated because subsequence ended then no point in subsequent
	NxtSeqIdx = -1;		// checking for more hypercores in same subsequence as these would be internal to this core

// RefHyperCoreLen holds the reference species hypercore length excluding any InDels on the ref species sequence
if(MaxUltraCoreLen < pProcParams->MinUltraLen || 
    RefHyperCoreLen < pProcParams->MinHyperLen)	// if less than required then return where to start next core from
	return(NxtSeqIdx);

ChromOfsEnd = ChromOfs+RefHyperCoreLen-1;	// determine reference chromosomal offset at which sequence ends

// ensure that the hypercore is in an included region and not part of an excluded region
if(!IncludeFilter(RefChromID,ChromOfs,ChromOfsEnd,pProcParams))
		return(NxtSeqIdx);

// classify range
if((pRange = GetLengthRangeClass(RefHyperCoreLen))==NULL)
	{
	gDiagnostics.DiagOut(eDLWarn,gszProcName,"Unexpected length range classification error - length = %d",RefHyperCoreLen);
	return(NxtSeqIdx);
	}

pStep = pProcParams->pCntStepCnts;
pStep += (pRange->ID-1) * pProcParams->Regions;
RegionIdx = 0;
FeatureBits = 0;
SpliceSiteOverlaps = 0;
if(pProcParams->pBiobed != NULL)
	{
	if(pProcParams->BEDChromID > 0)
		FeatureBits = pProcParams->pBiobed->GetFeatureBits(pProcParams->BEDChromID,ChromOfs,ChromOfsEnd,cRegionFeatBits,pProcParams->UpDnStreamLen);
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
				if(RegionIdx)		// if already have feature
					return(NxtSeqIdx);	// although was sequence of interest, more than one feature bit so can't contribute to stats
				RegionIdx = FeatIdx;
				if(pProcParams->bMultipleFeatBits)	// if multiple features allowed then don't check for any additional
					break;
				}
			}
		}

		// need to remap RegionIdx when incr counts so regions are displayed in a 5'->3' logical order
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
	pStep[RegionIdx] += 1;	

	if(RegionIdx != 0)		// if not intergenic then check for splice sites
		{
		SpliceSiteOverlaps = pProcParams->pBiobed->GetSpliceSiteBits(pProcParams->BEDChromID,ChromOfs,ChromOfsEnd,cMinSpliceOverlap);
		if(SpliceSiteOverlaps & cIntronExonSpliceSite)
			pStep[7]++;
		if(SpliceSiteOverlaps & cExonIntronSpliceSite)
			pStep[8]++;
		}
	pProcParams->bStatsAvail = true;
	ChkOutputResults(pProcParams->szRefChrom,ChromOfs,pProcParams,false,false);	// output results as may be appropriate
	}
else
	{
	*pStep += 1;
	ChkOutputResults(pProcParams->szRefChrom,ChromOfs,pProcParams,false,false);	// output results as may be appropriate
	}

memset(OGDistProfile,0,sizeof(OGDistProfile));
if(pProcParams->ProcMode == eProcModeOutspecies)
	{
	if(pProcParams->NumSpeciesInAlignment == pProcParams->MinAlignSpecies)
		OGUnaligned = RefHyperCoreLen;
	else
		{
		int VMatches;
		int CoreOfs;
		
		int CoreLoci = 0;
		int CoreLenLeft = RefHyperCoreLen;
		int NumSegsLeft = pProcParams->NumDistSegs;
		int SegLoci = CoreLenLeft / NumSegsLeft--;
		int OGSegIdx = 0;
		for(CoreOfs = 0; CoreOfs < HypercoreLen ; CoreOfs++)
			{
			// determine in the current alignment column the number of matches
			VMatches = 0;
			for(VIdx = 0; VIdx < NumCoreSpecies; VIdx++)
				{
				pVSeq = pProcParams->pSeqs[VIdx];
				VBase = pVSeq[CoreOfs+SeqIdx] & ~cRptMskFlg;
				if(VIdx == pProcParams->RefSpeciesIdx) // if base for reference sequence
					{
					RefBase = VBase;
					if(VBase == eBaseN || VBase == eBaseInDel)	  // indeterminate ref seq bases and InDels are 
						break;									  // treated as if mismatch
					}
				else
					if(VBase != RefBase)
						break;
				VMatches += 1;
				}

			if(VMatches == NumCoreSpecies)				// if all core bases in column match then check on outspecies
				{
				pVSeq = pProcParams->pSeqs[pProcParams->NumSpeciesInAlignment-1];
				OGBase = pVSeq[CoreOfs+SeqIdx] & ~cRptMskFlg;
				
				if(OGBase == RefBase)					// profile outgroup bases
					OGDistProfile[OGSegIdx].Matches += 1;
				else
					{
					if(OGBase == eBaseUndef)
						OGDistProfile[OGSegIdx].Unaligned += 1;
					else
						{
						if(OGBase == eBaseInDel)
							OGDistProfile[OGSegIdx].InDels += 1;
						else
							OGDistProfile[OGSegIdx].Mismatches += 1;
						}
					}
				CoreLoci += 1;
				CoreLenLeft -= 1;
				if(NumSegsLeft && CoreLoci >= SegLoci)
					{
					OGSegIdx += 1;
					SegLoci += CoreLenLeft / NumSegsLeft--;
					}
				}
			}
		}
	}

if(!OutputHypercore(pProcParams->szRefChrom,ChromOfs,ChromOfsEnd,FeatureBits | SpliceSiteOverlaps,OGUnaligned,OGMatches,OGMismatches,OGInDels,OGDistProfile,pProcParams))
	return(-1);
if(NxtSeqIdx != -1)
	NxtSeqIdx = HypercoreLen+SeqIdx+1;
return(NxtSeqIdx);
}

// ChkOutputResults
bool
ChkOutputResults(const char *pszChrom, int ChromOffset, tsProcParams *pProcParams,bool bGenEmptyRows,bool bGenEmptyWindows)
{
if(!pProcParams->WindowSize)
	return(false);
while(ChromOffset < 0 || ChromOffset > pProcParams->NxtOutputOffset)
	{
	if(bGenEmptyWindows || pProcParams->bStatsAvail)
		{
		OutputResults(pszChrom, pProcParams->NxtOutputOffset, pProcParams,bGenEmptyRows);
		if(pProcParams->bStatsAvail)
			memset(pProcParams->pCntStepCnts,0,pProcParams->NumCnts * pProcParams->Regions * sizeof(int));
		}
	pProcParams->bStatsAvail = false;
	pProcParams->NxtOutputOffset += pProcParams->WindowSize;
	if(ChromOffset < 0)
		break;
	}
return(true);
}



bool	
OutputResults(const char *pszChrom, int ChromOffset, tsProcParams *pProcParams,bool bGenEmptyRows)
{
tsLenRangeClass *pRange;
static bool bOutputHdrFirst = true;
char szLineBuff[2058];
int Idx;
int Steps;
int Instances;
int Len;
int *pStep;
pStep = pProcParams->pCntStepCnts;

if(bOutputHdrFirst)
	{
	bOutputHdrFirst = false;
	if(pProcParams->Regions != 9)
		Len = sprintf(szLineBuff,"\"LenRange\",\"Mismatches\",\"TotInstances\"");
	else
		Len = sprintf(szLineBuff,"\"LenRange\",\"Mismatches\",\"TotInstances\",\"IG\",\"US\",\"5'UTR\",\"CDS\",\"INTRON\",\"3'UTR\",\"DS\",\"5'ExSplice\",\"3'ExSplice\"");
	CUtility::SafeWrite(pProcParams->hRsltsFile,szLineBuff,Len);
	}
pStep = pProcParams->pCntStepCnts;
for(Idx = 0; Idx < cLenRanges; Idx++, pStep += pProcParams->Regions)
	{
	if(pProcParams->Regions > 1)
		{
		for(Instances = Steps = 0; Steps < (pProcParams->Regions - 2); Steps++)
			Instances += pStep[Steps];
		}
	else
		Instances = pStep[0];
	pRange = GetRangeClass(Idx+1);
	Len = sprintf(szLineBuff,"\n\"%s\",%d,%d",
						pRange->pszDescr,pProcParams->MaxHyperMismatches,Instances);
	if(pProcParams->Regions > 1)
		for(Steps = 0; Steps < pProcParams->Regions; Steps++)
			Len += sprintf(&szLineBuff[Len],",%d",pStep[Steps]);
	CUtility::SafeWrite(pProcParams->hRsltsFile,szLineBuff,Len);
	}
return(true);
}


bool	
OutputHypercore(const char *pszChrom, int ChromStartOffset, int ChromEndOffset, 
				int FeatureBits,		// feature bits over lapped
				int OGUnaligned,		// number of unaligned bases in outspecies
				int OGMatches,			// number of matching bases in outspecies
				int OGMismatches,		// number of mismatched bases in outspecies
				int OGInDels,			// number of InDels in outspecies
				tsDistSeg SegCnts[],	// array of segment profile counts	
				tsProcParams *pProcParams)
{
static int CoreID = 0;
int Rslt;
char szLineBuff[4096];
int Len;
if(pProcParams->hCoreCSVRsltsFile != -1)
	{
	Len = sprintf(szLineBuff,"%d,\"%s\",\"%s\",\"%s\",%d,%d,%d,\"%s\",%d",
		++CoreID,pProcParams->MinHyperLen ? "hypercore" : "ultracore",
		pProcParams->szSpecies[pProcParams->RefSpeciesIdx],
		pszChrom,ChromStartOffset,ChromEndOffset,ChromEndOffset-ChromStartOffset+1,
		pProcParams->pszSpeciesList,FeatureBits & (cAnyFeatBits | cOverlaysSpliceSites));
	if(pProcParams->ProcMode == eProcModeOutspecies)
		{
		Len += sprintf(&szLineBuff[Len],",%d,%d,%d,%d",OGUnaligned,OGMatches,OGMismatches,OGInDels);
		for(int Idx = 0; Idx < pProcParams->NumDistSegs; Idx++)
			Len += sprintf(&szLineBuff[Len],",%d,%d,%d,%d",SegCnts[Idx].Matches,SegCnts[Idx].Mismatches,SegCnts[Idx].InDels,SegCnts[Idx].Unaligned);
		}
	Len += sprintf(&szLineBuff[Len],"\n");
	if((Rslt=write(pProcParams->hCoreCSVRsltsFile,szLineBuff,Len))!=Len)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Write to loci file failed - %s",strerror(errno));
		return(false);
		}
	if(CoreID == 1 || !(CoreID % 500))
		_commit(pProcParams->hCoreCSVRsltsFile);
	}
return(true);
}

// ChkOutputSummaryResults
bool
ChkOutputSummaryResults(const char *pszChrom, int ChromOffset, tsProcParams *pProcParams,bool bGenEmptyRows,bool bGenEmptyWindows)
{
if(!pProcParams->WindowSize)
	return(false);
while(ChromOffset < 0 || ChromOffset > pProcParams->NxtOutputOffset)
	{
	if(bGenEmptyWindows || pProcParams->bStatsAvail)
		{
		OutputSummaryResults(pszChrom, pProcParams->NxtOutputOffset, pProcParams,bGenEmptyRows);
		if(pProcParams->bStatsAvail)
			memset(pProcParams->pCntStepCnts,0,pProcParams->NumCnts * pProcParams->Regions * sizeof(int));
		}
	pProcParams->bStatsAvail = false;
	pProcParams->NxtOutputOffset += pProcParams->WindowSize;
	if(ChromOffset < 0)
		break;
	}
return(true);
}



bool	
OutputSummaryResults(const char *pszChrom, int ChromOffset, tsProcParams *pProcParams,bool bGenEmptyRows)
{
static bool bOutputHdrFirst = true;
static char *pszClass[] = {(char *)"RefInDel",(char *)"RelIndel",(char *)"Match",(char *)"Missmatch"};
char szLineBuff[2058];
int Idx;
int Steps;
int Len;
int *pStep;
INT64 Total;
INT64 SumRefBases;
INT64 SumRelBases;

if(bOutputHdrFirst)
	{
	bOutputHdrFirst = false;
	Len = sprintf(szLineBuff,"\"Class\",\"TotInstances\"");
	if(pProcParams->Regions == 9)
		Len += sprintf(&szLineBuff[Len],",\"IG\",\"US\",\"5'UTR\",\"CDS\",\"INTRON\",\"3'UTR\",\"DS\",\"5'ExSplice\",\"3'ExSplice\"");
	CUtility::SafeWrite(pProcParams->hRsltsFile,szLineBuff,Len);
	}

pStep = pProcParams->pCntStepCnts;
for(Idx = 0; Idx < 4; Idx++,pStep += pProcParams->Regions)
	{
	for(Total = Steps = 0; Steps < pProcParams->Regions; Steps++)
		Total += pStep[Steps];
	switch(Idx) {
		case 0:
			SumRelBases = Total;
			break;
		case 1:
			SumRefBases = Total;
			break;
		default:
			SumRelBases += Total;
			SumRefBases += Total;
			break;
		}
#ifdef _WIN32
	Len = sprintf(szLineBuff,"\n\"%s\",%I64d",pszClass[Idx],Total);
#else
	Len = sprintf(szLineBuff,"\n\"%s\",%lld",pszClass[Idx],Total);
#endif
	if(pProcParams->Regions > 1)
		{
		for(Steps = 0; Steps < pProcParams->Regions; Steps++)
			Len += sprintf(&szLineBuff[Len],",%d",pStep[Steps]);
		}
	CUtility::SafeWrite(pProcParams->hRsltsFile,szLineBuff,Len);
	}

pStep = pProcParams->pCntStepCnts;
for(Idx = 0; Idx < 2; Idx++)
	{
	switch(Idx) {	// reference sequence
		case 0:
#ifdef _WIN32
			Len=sprintf(szLineBuff,"\n\"RefBases\",%I64d",SumRefBases);
#else
			Len=sprintf(szLineBuff,"\n\"RefBases\",%lld",SumRefBases);
#endif
			break;
		case 1:		// relative sequence
#ifdef _WIN32
			Len=sprintf(szLineBuff,"\n\"RelBases\",%I64d",SumRelBases);
#else
			Len=sprintf(szLineBuff,"\n\"RelBases\",%lld",SumRelBases);
#endif
			break;
		}

	if(pProcParams->Regions > 1)
		{
		for(Steps = 0; Steps < pProcParams->Regions; Steps++)
			{
			switch(Idx) {
				case 0:	// reference sequence
					SumRefBases = pStep[pProcParams->Regions + Steps] + 
						pStep[(pProcParams->Regions * 2) + Steps] + 
						pStep[(pProcParams->Regions * 3) + Steps];
					break;
				case 1:	// relative sequence
					SumRefBases = pStep[Steps] + 
						pStep[(pProcParams->Regions * 2) + Steps] + 
						pStep[(pProcParams->Regions * 3) + Steps];
					break;
				}
#ifdef _WIN32
			Len += sprintf(&szLineBuff[Len],",%I64d",SumRefBases);
#else
			Len += sprintf(&szLineBuff[Len],",%lld",SumRefBases);
#endif
			}
		}
	CUtility::SafeWrite(pProcParams->hRsltsFile,szLineBuff,Len);
	}
return(true);
}



