// genlAlignLoci2core.cpp : Defines the entry point for the console application.


#include "stdafx.h"
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#if _WIN32

#include "../libbiokanga/commhdrs.h"
#else
#include "../libbiokanga/commhdrs.h"
#endif

const unsigned int cProgVer = 308;			// increment with each release

const int cMAFReqSpecies = 2;				// number of species to be specified
const int cMAtotMaxSeqAlignLen = 0x0ffffff; // total (over all aligned species) max seq length that can be buffered in concatenated seqs

const int cMinCoreLen = 4;				// allow core lengths to be specified down to cMinCoreLen
const int cDfltCoreLen= 10;				// if core lengths not specified then default to cDfltMinCoreLen
const int cMaxCoreLen = 1000000;		// minimum core lengths can be specified upto this length

const int cMaxIncludeFiles = 10;		// maximun number of include region filter files
const int cMaxExcludeFiles = 10;		// maximun number of exclude region filter files

const int cMaxIncludeChroms = 20;		// max number of include chromosomes regular expressions
const int cMaxExcludeChroms = 20;		// max number of exclude chromosomes regular expressions

const int cCoreLociAllocChunk = 10000;	// grow core loci array by this many elements
const int cCoreChromAllocChunk = 1000;	// grow core chrom array by this many elements

const int cMinMatchDistSegments  = 4;	// min match distribution profile segments
const int cDfltMatchDistSegments = 10;	// default match distribution profile segments
const int cMaxMatchDistSegments  = 100;	// max match distribution profile segments

// default blastz score for gaps and non-aligned bases
const int cInDelAffine = -400;		// InDel gap opening score
const int cInDelExtn = -30;			// InDel gap extension score
const int cUnaligned = -1000;		// if outgroup base is unaligned
const int cScoreAA = 91;
const int cScoreAC = -114;
const int cScoreAG = -31;
const int cScoreAT = -123;
const int cScoreAN = -100;

const int cScoreCA = -114;
const int cScoreCC = 100;
const int cScoreCG = -125;
const int cScoreCT = -31;
const int cScoreCN = -100;

const int cScoreGA = -31;
const int cScoreGC = -125;
const int cScoreGG = 100;
const int cScoreGT = -114;
const int cScoreGN = -100;

const int cScoreTA = -123;
const int cScoreTC = -31;
const int cScoreTG = -114;
const int cScoreTT = 91;
const int cScoreTN = -100;

const int cScoreNA = -100;
const int cScoreNC = -100;
const int cScoreNG = -100;
const int cScoreNT = -100;
const int cScoreNN = -100;

// processing mode
typedef enum eProcMode {
	eProcStandard						// standard processing is to output each ref loci + rel alignment distribution
} etProcMode;

// expected input loci file type
typedef enum eLociFileType {
	eLociFileTypeCSV = 0,				// default processing is from CSV file
	eLociFileTypeUCSCBED,				// or from raw UCSC BED file
	eLociFileTypeBIOBED					// or from biobed feature file
} etLociFileType;

typedef struct TAG_sExcludeSpeciesChrom {
	int SpeciesID;			// which species
	int ChromID;			// chromosome is not to be processed
	} tsExcludeSpeciesChrom;

typedef struct TAG_sCoreLoci {
	int RefID;								// reference identifier from loci CSV or BED 
	int ChromID;
	char Strand;
	int StartLoci;
	int EndLoci;
} tsCoreLoci;

typedef struct TAG_sCoreChrom {
	int ChromID;							// unique identifier for this chromosome
	tsCoreLoci *pCoreLoci;					// pts to first core on this chromosome
	UINT16 Hash;							// hash over chromosome name
	char szChrom[cMaxDatasetSpeciesChrom];	// core chromosome
} tsCoreChrom;


typedef struct TAG_sLenRangeClass {
	int ID;					// uniquely identifies this range
	int Min;				// minimum length in this range
	int Max;				// maximum length in this range
	char *pszDescr;			// descriptive text
	}tsLenRangeClass;

typedef struct TAG_sDistSeg {
	int Matches;			// number of exact matches in this segment
	int Mismatches;			// number of mismatches
	int InDels;				// total number of InDel bases
	int Unaligned;			// number of unaligned bases
	} tsDistSeg;

typedef struct TAG_sMatrixScore {
	int Unaligned;			// unaligned base
	int InDelAffine;		// InDel opening score
	int InDelExtn;			// InDel extension score
	int Score[5][5];		// to include scoring for indeterminate nucleotides eBaseN, row == ref, col == outgroup	
} tsMatrixScore;


typedef struct TAG_sProcParams 
	{
	etProcMode ProcMode;			// processing mode
	int MaxNumSpecies;				// maximum number of species required in alignment block (2, ref + rel)
	int MinNumSpecies;				// minimum number of species required in alignment block (1, ref)

	int NumBlockSpecies;			// number of species in currently processed block(s)
	int AlignBlockID;				// current alignment block identifier
	char *pszSpeciesList;			// comma separated species list starting with reference species 
	int RefSpeciesIdx;				// current index into pSeqs for the reference species
	char szSpecies[cMaxAlignedSpecies][cMaxDatasetSpeciesChrom];	// species names of interest - other species are sloughed
									// only alignments with species included will be processed
									// first species is the reference species
	int MaxAlignIdxSpecies;			// pSeqs[MaxAlignIdxSpecies-1] is last actual alignment sequence
	char szRefChrom[cMaxDatasetSpeciesChrom]; // current reference chromosome
	int RefChromID;					// current reference chromosome identifier
	int NxtOutputOffset;			// when to next output results - i.e end of current window
	int SeqOfs;						// offset into sequences
	int SeqLen;						// sequence length
	etSeqBase *pSeqs[cMaxAlignedSpecies];  // ptrs to each alignment sequence

	char szRsltsFile[_MAX_PATH];	// output results file to create/write
	int hRsltsFile;					// opened file handle for output results file

	char szSummaryFile[_MAX_PATH];	// output summary file to create/write
	int hSummaryFile;				// opened file handle for summary file

	int MaxSeqAlignLen;				// max length alignment which can be processed (how much mem was alloc'd to pSeq[n])			
	int NumIncludes;				// number of biobed files containing regions to include
	int NumExcludes;				// number of biobed files containing regions to exclude
	CBEDfile *pIncludes[cMaxIncludeFiles];	// if opened biobed files for regions to include - all other regions are to be excluded
	CBEDfile *pExcludes[cMaxExcludeFiles];	// if opened biobed files for regions to exclude 
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
	etLociFileType LociFileType;	// expected loci file type
	char szLociFile[_MAX_PATH];		// name of file containing core loci
	char szFeaturesFile[_MAX_PATH];	// features charaterisation biobed file
	CBEDfile *pBEDfeatures;			// if not NULL then opened biobed file with regional characterisations
	int BEDChromID;					// BEDfeatures chromosome identifier corresponding to RefChromID
	int UpDnStreamLen;				// up/dn stream regional length when characterising
	bool bMultipleFeatBits;			// if false then stats only generated if a single feature bit is set - e.g if both exons and introns overlapped then no stat generated

	CBEDfile *pBEDloci;				// if not NULL then opened biobed file with core loci
	CCSVFile *pCSVloci;				// if not NULL then opened CSV file with core loci

	tsCoreLoci *pCoreLocs;			// pts to mem allocd to hold array of core loci
	int MaxCoreLocs;				// max allocd core locii
	int NumCoreLocs;				// number of core loci ptd at by pCoreLocs

	int CachCoreChromID;			// identifier of last retrieved core chromosome from pCoreChroms
	tsCoreChrom *pCoreChroms;		// pts to mem allocd to hold array of core chromosome names
	int MaxCoreChroms;				// max allocd core chroms
	int NumCoreChroms;				// number of core chroms ptd at by pCoreLocs

	int	MinCoreLen;					// minimum core length required
	int	MaxCoreLen;					// maximum core length required

	int NumDistSegs;				// number of match distribution profile segments

	int *pCntStepCnts;				// array of stats counters
	int NumCnts;					// number of steps in pCntStepCnts
	int Regions;					// number of regions per step

	CFilterRefIDs *pFilterRefIDs;	// used when filtering by RefID

	tsMatrixScore MatrixScore;		// nucleotide scoring matrix

} tsProcParams; 

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


// length range classes
tsLenRangeClass LenRangeClasses[] = {
	{1,0,4,(char *)"0-4"},
	{2,5,9,(char *)"5-9"},
	{3,10,14,(char *)"10-14"},
	{4,15,19,(char *)"15-19"},
	{5,20,29,(char *)"20-29"},
    {6,30,49,(char *)"30-49"},
	{7,50,74,(char *)"50-74"},
	{8,75,99,(char *)"75-99"},
	{9,100,124,(char *)"100-124"},
	{10,125,149,(char *)"125-149"},
	{11,150,174,(char *)"150-174"},
	{12,175,199,(char *)"175-199"},
	{13,200,249,(char *)"200-249"},
	{14,250,299,(char *)"250-299"},
	{15,300,349,(char *)"300-349"},
	{16,350,399,(char *)"350-399"},
	{17,400,449,(char *)"400-449"},
	{18,450,499,(char *)"450-499"},
	{19,500,599,(char *)"500-599"},
	{20,600,699,(char *)"600-699"},
	{21,700,799,(char *)"700-799"},
	{22,800,899,(char *)"800-899"},
	{23,900,999,(char *)"900-999"},
	{24,1000,1249,(char *)"1000-1249"},
	{25,1250,1499,(char *)"1250-1499"},
	{26,1500,1749,(char *)"1500-1749"},
	{27,1750,1999,(char *)"1750-1999"},
	{28,2000,INT_MAX,(char *)"2000+"}
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
	return _T("genlalignloci2core");
}
// end of str library required code
#endif

int
Process(etProcMode ProcMode,						// processing mode
				etLociFileType LociFileType,		// expected loci file type
				 char *pszInputFile,				// core loci (CSV or Biobed .bsb) file to process
		 		 char *pszFilterRefIDFile,			// exclude any RefIDs in this filter file 
 				 char *pszMAFFile,					// multiple alignment file
					char *pszRsltsFile,				// file to write results into
					char *pszSummaryFile,			// file to write summary into
					char *pszFeaturesFile,			// feature charaterisation biobed file
					int NumDistSegs,				// number of match distribution profile segments
					int RegLen,						// regulatory region length - up/dn stream of 5/3' 
					bool bMultipleFeatBits,			// if true then accept alignments in which multiple feature bits are set
					int NumIncludeFiles,			// number of include region files
					char **ppszIncludeFiles,		// biobed files containing regions to include - default is to exclude none
					int NumExcludeFiles,			// number of exclude region files
					char **ppszExcludeFiles,		// biobed file containing regions to exclude - default is to include all
					char *pszSpeciesList,			// space or comma separated list of species, priority ordered
					int	MinCoreLen,					// minimum core length required
					int	MaxCoreLen,					// maximum core length required
					int NumIncludeChroms,			// number of chromosomes explicitly defined to be included
					char **ppszIncludeChroms,		// ptr to array of reg expressions defining chroms to include - overides exclude
					int NumExcludeChroms,			// number of chromosomes explicitly defined to be excluded
					char **ppszExcludeChroms,		// ptr to array of reg expressions defining chroms to include
					char *pszScoreMatrix);			// scoring matrix file to use


char *ProcMode2Txt(etProcMode ProcMode);
int ParseNumSpecies(char *pszSpeciesList,tsProcParams *pProcParams);
int TrimQuotes(char *pszTxt);
bool IncludeFilter(int RefChromID,int SubRefOfs,int SubRefEndOfs,tsProcParams *pProcParams);
bool CleanupResources(tsProcParams *pProcParams);
CBEDfile *OpenBedfile(char *pToOpen,bool bReportErrs);
bool ExcludeThisChrom(CMAlignFile *pAlignments,int SpeciesID,int ChromID,tsProcParams *pProcParams);
tsExcludeEl *LocateExclude(int SpeciesID,int ChromID);
bool AddExcludeHistory(int SpeciesID,int ChromID,bool bExclude);
int AddLoci(int RefID,char *pszChrom,char Strand,int StartLoci,int Len,tsProcParams *pParams);
int AddChrom(char *pszChrom,tsProcParams *pParams);
int LoadCSVloci(tsProcParams *pParams);
int LoadBIOBEDloci(tsProcParams *pParams);
int LoadUCSCBEDloci(char *pszFile,tsProcParams *pParams);
int ProcessLoci(tsProcParams *pParams);

bool	
OutputLociCore(tsCoreLoci *pCore,		// current core
			   char *pszChrom, int ChromStartOffset, int ChromEndOffset, 
				int FeatureBits,		// feature bits over lapped
				int OGUnaligned,		// number of unaligned bases in outspecies
				int OGMatches,			// number of matching bases in outspecies
				int OGMismatches,		// number of mismatched bases in outspecies
				int OGInDels,			// number of InDels in outspecies
				int Score,				// alignment score
				int RefScore,			// background ref seq score
				tsDistSeg SegCnts[],	// array of segment profile counts	
				tsProcParams *pProcParams);


bool OutputSummary(tsProcParams *pProcParams);
int ProcLociSeq(tsCoreLoci *pCore,		// current core
				int StartSeqIdx,			// index into pProcParams->pSeqs[0][] at which loci mapped alignment starts 
		   int EndSeqIdx,				// index into pProcParams->pSeqs[0][] at which loci mapped alignment ends 
		   char *pszRefChrom,			// reference chromosome
		   int StartLoci,				// start loci on chromosome
		   int EndLoci,					// end loci on chromosome
		   tsProcParams *pProcParams); // global processing parameters


int 
ProcessAlignments(char *pszMAF,			 // source bioseq multialignment file
				  tsProcParams *pProcParams); // processing parameters
int NormaliseInDelColumns(tsProcParams *pProcParams,int AlignLen);
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
bool 
ProcAlignBlock(int RefChromID,	    // reference chromosome
			   int RefChromOfs,		// offset along reference chromosome 
			   int AlignLen,		// alignment length incl any InDels
			   tsProcParams *pProcParams); // global processing parameters

tsCoreLoci *
GetFirstLociOverlaps(char *pszChrom,			// reference species chromosome
   				int RefChromID,			// chrom identifier
				char RefStrand,			// strand
				int RefChromOfs,		// start loci
				int RefChromEndOfs,		// end loci
				tsProcParams *pProcParams); // global processing parameters

tsCoreLoci *GetLociOverlaps(char *pszChrom,			// reference species chromosome
				int RefChromID,			// chrom identifier
				char RefStrand,			// strand
				int RefChromOfs,		// start loci
				int RefChromEndOfs,		// end loci
				tsCoreLoci *pPrevCoreLoci, // previous loci or NULL if none
				tsProcParams *pProcParams); // global processing parameters

UINT16 GenNameHash(char *pszName);

static int CompareCoreEls( const void *arg1, const void *arg2);


CStopWatch gStopWatch;
CDiagnostics gDiagnostics;				// for writing diagnostics messages to log file
char gszProcName[_MAX_FNAME];			// process name

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

int iProcMode;
int iLociFileType;
int Idx;
char szRsltsFile[_MAX_PATH];
char szSummaryFile[_MAX_PATH];
char szInputFile[_MAX_PATH];
char szMAFFile[_MAX_PATH];
char szFeaturesFile[_MAX_PATH];
char szFilterRefIDFile[_MAX_PATH];  // exclude any RefIDs in this filter file

char szScoreMatrix[_MAX_PATH];		// alignment scoring matrix file

int NumIncludeChroms;
char *pszIncludeChroms[cMaxIncludeChroms];
int NumExcludeChroms;
char *pszExcludeChroms[cMaxExcludeChroms];
int NumIncludeFiles;
char *pszIncludeFiles[cMaxIncludeFiles];
int NumExcludeFiles;
char *pszExcludeFiles[cMaxExcludeFiles];
char szSpeciesList[cMaxAlignedSpecies*cMaxDatasetSpeciesChrom];
int iNumSpecies;

int iRegLen;
bool bMultipleFeatBits;
int iDistSegs;			// match distribution profile segments

int	iMinCoreLen;		// minimum core length required
int iMaxCoreLen;		// maximum core length required

int LenFileList;

// command line args
struct arg_lit  *help    = arg_lit0("hH","help",					"print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",				"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",			"<int>","Level of diagnostics written to logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_int *ScreenLogLevel=arg_int0("S", "ScreenLogLevel",		"<int>","Level of diagnostics written to screen 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",			"diagnostics log file");

struct arg_file *InFile = arg_file1("i",NULL,"<file>",				"input loci of interest from .csv or BED file");
struct arg_int  *LociFileType = arg_int0("l","locitype","<int>",	"loci file type 0:CSV loci, 1: raw UCSC BED, 2:Biobed BED (default: CSV)");
struct arg_file *FilterRefIDFile = arg_file0("X",NULL,"<file>",		"filter out any loci (must be CSV file) with RefIDs in this filter file");
struct arg_file *MAFFile = arg_file1("I",NULL,"<file>",				"input from bioseq multiple alignment file");
struct arg_int  *ProcMode = arg_int0("x","procmode","<int>",		"mode 0: default");
struct arg_file *OutFile = arg_file1("o",NULL,"<file>",				"output to results file");
struct arg_file *SummaryFile = arg_file0("O",NULL,"<file>",			"output to summary file");

struct arg_str  *SpeciesList = arg_str1("s","species","<string>",	"species list, must contain ref + rel species with ref the first");

struct arg_int  *MinCoreLen = arg_int0("m","mincorelen","<int>",	"only process elements >= this length (default= 20)");
struct arg_int  *MaxCoreLen = arg_int0("M","maxcorelen","<int>",	"only process elements <= this length (default = 1000000)");

struct arg_file *FeaturesFile = arg_file0("b",NULL,"<file>",			"Characterise features from this biobed file");

struct arg_lit  *MultipleFeatBits  = arg_lit0("q","multiplefeatbits",	"single featbit (default) or multiple featbits allowed");
struct arg_int  *RegLen = arg_int0("L","updnstream","<int>",	"length of 5'up or 3'down  stream regulatory region length (default = 2000) 0..1000000");
struct arg_int  *DistSegs = arg_int0("d","distsegs","<int>",	"number of match distribution segments (default 10) used when processing outspecies mode");

struct arg_file *ExcludeFile = arg_filen("E","exclude","<file>",0,cMaxExcludeFiles,	"exclude all regions in biobed file from processing ");
struct arg_file *IncludeFile = arg_filen("I","include","<file>",0,cMaxIncludeFiles,	"include all regions (unless specific regions excluded) in biobed file");
struct arg_str  *IncludeChroms = arg_strn("z","chromeinclude","<string>",0,cMaxIncludeChroms,"low priority - regular expressions defining species.chromosomes to include for processing");
struct arg_str  *ExcludeChroms = arg_strn("Z","chromexclude","<string>",0,cMaxExcludeChroms,"high priority - regular expressions defining species.chromosomes to exclude from processing");

struct arg_file *ScoreMatrix = arg_file0("scorematrix",NULL,"<file>",	"outgroup score matrix");

struct arg_end *end = arg_end(20);

void *argtable[] = {help,version,FileLogLevel,ScreenLogLevel,LogFile,
					InFile,FilterRefIDFile,MAFFile,OutFile,SummaryFile,MultipleFeatBits,RegLen,DistSegs,
					SpeciesList,MinCoreLen,MaxCoreLen,
					ProcMode,LociFileType,FeaturesFile,ExcludeFile,IncludeFile,IncludeChroms,ExcludeChroms,
					ScoreMatrix,
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
		printf("\n      To invoke this parameter file then precede its name with '@'");
		printf("\n      e.g. %s @myparams.txt\n\n",gszProcName);
		exit(1);
        }

    /* special case: '--version' takes precedence error reporting */
if (version->count > 0)
        {
		printf("\n%s: Version: %d.%2.2d\n",gszProcName,cProgVer/100,cProgVer%100);
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

	iProcMode = ProcMode->count ? ProcMode->ival[0] : eProcStandard;
	if(iProcMode < eProcStandard || iProcMode > eProcStandard)
		{
		printf("\nError: Requested processing mode '-x%d' not supported",iProcMode);
		exit(1);
		}

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


	iMinCoreLen = MinCoreLen->count ? MinCoreLen->ival[0] : cDfltCoreLen;
	if(iMinCoreLen < cMinCoreLen)
		{
		printf("\nSpecified minimum core length '-n%d' < %d, assuming you meant %d",iMinCoreLen,cMinCoreLen,cMinCoreLen);
		iMinCoreLen = cMinCoreLen;
		}
	else
		{
		if(iMinCoreLen > cMaxCoreLen)
			{
			printf("\nSpecified minimum core length '-n%d' > %d, assuming you meant %d",iMinCoreLen,cMaxCoreLen,cMaxCoreLen);
			iMinCoreLen = cMaxCoreLen;
			}
		}


	iMaxCoreLen = MaxCoreLen->count ? MaxCoreLen->ival[0] : cMaxCoreLen;
	if(iMaxCoreLen < iMinCoreLen)
		{
		printf("\nSpecified maximum hyper length '-n%d' < %d, assuming you meant %d",iMaxCoreLen,iMinCoreLen,iMinCoreLen);
		iMaxCoreLen = iMinCoreLen;
		}
	else
		{
		if(iMaxCoreLen > cMaxCoreLen)
			{
			printf("\nSpecified maximum core length was '-n%d' > %d, assuming you meant '-n%d'",iMaxCoreLen,cMaxCoreLen,cMaxCoreLen);
			iMaxCoreLen = cMaxCoreLen;
			}
		}


	bMultipleFeatBits = MultipleFeatBits->count ? true : false;

	iDistSegs = DistSegs->count ? DistSegs->ival[0] : cDfltMatchDistSegments;
	if(iDistSegs < cMinMatchDistSegments)
		{
		printf("\nWarning: too few distribution segments specified with '-d%d', assume you meant %d segments", iDistSegs,cMinMatchDistSegments);
		iDistSegs = cMinMatchDistSegments;
		}
	else
		if(iDistSegs > min(iMinCoreLen,cMaxMatchDistSegments))
			{
			printf("\nWarning: too many distribution segments specified with '-d%d', assume you meant %d segments", iDistSegs,min(iMinCoreLen,cMaxMatchDistSegments));
			iDistSegs = min(iMinCoreLen,cMaxMatchDistSegments);
			}

	iLociFileType = LociFileType->count ? LociFileType->ival[0] : eLociFileTypeCSV;
	if(iLociFileType < eLociFileTypeCSV || iLociFileType > eLociFileTypeBIOBED)
		{
		printf("\nError: Requested loci file type '-l%d' not supported",iLociFileType);
		exit(1);
		}
	else
		if(InFile->count == 0)
			{
			printf("\nError: Input loci file '-i<file>' must be specified");
			exit(1);
			}

	if(FilterRefIDFile->count)
		{
		if(iLociFileType != eLociFileTypeCSV)
			{
			printf("\nError: Filtering by RefID only available if input loci file type is CSV");
			exit(1);
			}
		strncpy(szFilterRefIDFile,FilterRefIDFile->filename[0],sizeof(szFilterRefIDFile));
		szFilterRefIDFile[sizeof(szFilterRefIDFile)-1] = '\0';
		}
	else
		szFilterRefIDFile[0] = '\0';

	if(FeaturesFile->count)
		{
		strncpy(szFeaturesFile,FeaturesFile->filename[0],sizeof(szFeaturesFile));
		szFeaturesFile[sizeof(szFeaturesFile)-1] = '\0';
		}
	else
		szFeaturesFile[0] = '\0';

	if(ScoreMatrix->count)
		{
		strncpy(szScoreMatrix,ScoreMatrix->filename[0],sizeof(szScoreMatrix));
		szScoreMatrix[sizeof(szScoreMatrix)-1] = '\0';
		}
	else
		szScoreMatrix[0] = '\0';

	if(SummaryFile->count)
		{
		if(szFeaturesFile[0] == '\0')
			{
			printf("\nError: Features file '-b<featurebedfile>' must be specified if summary '-O%s' requested",SummaryFile->filename[0]);
			exit(1);
			}
		strncpy(szSummaryFile,SummaryFile->filename[0],sizeof(szSummaryFile));
		szSummaryFile[sizeof(szSummaryFile)-1] = '\0';
		}
	else
		szSummaryFile[0] = '\0';


	strncpy(szInputFile,InFile->filename[0],sizeof(szInputFile));
	szInputFile[sizeof(szInputFile)-1] = '\0';

	strncpy(szRsltsFile,OutFile->filename[0],sizeof(szRsltsFile));
	szRsltsFile[sizeof(szRsltsFile)-1] = '\0';
	strncpy(szMAFFile,MAFFile->filename[0],sizeof(szMAFFile));
	szMAFFile[sizeof(szMAFFile)-1] = '\0';

	NumIncludeFiles = IncludeFile->count;
	for(Idx=0;Idx < IncludeFile->count; Idx++)
		{
		LenFileList = (int)strlen(IncludeFile->filename[Idx]);
		pszIncludeFiles[Idx] = new char [LenFileList+1];
		strcpy(pszIncludeFiles[Idx],IncludeFile->filename[Idx]);
		TrimQuotes(pszIncludeFiles[Idx]);
		}
	NumExcludeFiles = ExcludeFile->count;
	for(Idx=0;Idx < ExcludeFile->count; Idx++)
		{
		LenFileList = (int)strlen(ExcludeFile->filename[Idx]);
		pszExcludeFiles[Idx] = new char [LenFileList+1];
		strcpy(pszExcludeFiles[Idx],ExcludeFile->filename[Idx]);
		TrimQuotes(pszExcludeFiles[Idx]);
		}

	NumIncludeChroms = IncludeChroms->count;
	for(Idx=0;Idx < IncludeChroms->count; Idx++)
		{
		LenFileList = (int)strlen(IncludeChroms->sval[Idx]);
		pszIncludeChroms[Idx] = new char [LenFileList+1];
		strcpy(pszIncludeChroms[Idx],IncludeChroms->sval[Idx]);
		TrimQuotes(pszIncludeChroms[Idx]);
		}

	NumExcludeChroms = ExcludeChroms->count;
	for(Idx=0;Idx < ExcludeChroms->count; Idx++)
		{
		LenFileList = (int)strlen(ExcludeChroms->sval[Idx]);
		pszExcludeChroms[Idx] = new char [LenFileList+1];
		strcpy(pszExcludeChroms[Idx],ExcludeChroms->sval[Idx]);
		TrimQuotes(pszExcludeChroms[Idx]);
		}

	if(!SpeciesList->count)
		{
		printf("\nError: Species list '-s<specieslist>' is empty\n");
		exit(1);
		}
	strcpy(szSpeciesList,SpeciesList->sval[0]);
	TrimQuotes(szSpeciesList);

	iNumSpecies = ParseNumSpecies(szSpeciesList,NULL);
	if(iNumSpecies < cMAFReqSpecies)
		{
		printf("Error: Species list %s",iNumSpecies < 0 ? "is unable to be parsed" : "must contain 2 species");
		exit(1);
		}
	else
		if(iNumSpecies > cMAFReqSpecies)
			{
			printf("Error: Species list '%s' contains more than 2 species",szSpeciesList);
			exit(1);
			}

		// now that command parameters have been parsed then initialise diagnostics log system
	if(!gDiagnostics.Open(szLogFile,(etDiagLevel)iScreenLogLevel,(etDiagLevel)iFileLogLevel,true))
		{
		printf("\nError: Unable to start diagnostics subsystem.");
		if(szLogFile[0] != '\0')
			printf(" Most likely cause is that logfile '%s' can't be opened/created",szLogFile);
		exit(1);
		}

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Version: %d.%2.2d Processing parameters:",cProgVer/100,cProgVer%100);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Processing mode: %s",ProcMode2Txt((etProcMode)iProcMode));
	switch(iLociFileType) {
		case eLociFileTypeCSV:
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Loci CSV input file %s",szInputFile);
			break;

		case eLociFileTypeUCSCBED:
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Loci raw UCSC BED input file %s",szInputFile);
			break;

		case eLociFileTypeBIOBED:
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Loci Biobed input file %s",szInputFile);
			break;
		}

	if(szFilterRefIDFile[0])
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Exclude any RefIDs in this filter file: '%s'",szFilterRefIDFile);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Results file: '%s'",szRsltsFile);
	if(szSummaryFile[0] != '\0')
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Summary file: '%s'",szSummaryFile);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Input (.algn) multialignment file to process: '%s'",szMAFFile);

	for(Idx = 0; Idx < NumIncludeFiles; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"biobed file containing regions to include: '%s'",pszIncludeFiles[Idx]); 
	for(Idx = 0; Idx < NumExcludeFiles; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"biobed file containing regions to exclude: '%s'",pszExcludeFiles[Idx]); 
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Number of species: %d",iNumSpecies);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"species list: '%s'",	szSpeciesList);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"minimum core length required: %d",iMinCoreLen);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"maximum core length required: %d",iMaxCoreLen);
	if(szFeaturesFile[0]!= '\0')
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Characterise features from: '%s'",szFeaturesFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"regulatory region length: %d",iRegLen);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"accept alignments in which multiple feature bits are set: %s",bMultipleFeatBits ? "yes" : "no");
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Match distribution profile segments: %d",iDistSegs);

	for(Idx = 0; Idx < NumIncludeChroms; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"reg expressions defining chroms to include: '%s'",pszIncludeChroms[Idx]);
	for(Idx = 0; Idx < NumExcludeChroms; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"reg expressions defining chroms to exclude: '%s'",pszExcludeChroms[Idx]); 

	if(szScoreMatrix[0] != '\0')
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Nucleotide scoring matrix file : '%s'",szScoreMatrix);

	gStopWatch.Start();
#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	Rslt = Process((etProcMode)iProcMode,				// processing mode
					(etLociFileType)iLociFileType,		// expected input file type
					szInputFile,		// core loci (CSV or Biobed .bsb) file to process
					szFilterRefIDFile, // exclude any RefIDs in this filter file
					szMAFFile,			// multiple alignment file
					szRsltsFile,		// where to write out results
					szSummaryFile,		// where to write out summary
					szFeaturesFile,		// feature charaterisation biobed file
					iDistSegs,			// number of match distribution profile segments
					iRegLen,			// regulatory region length - up/dn stream of 5/3' 
					bMultipleFeatBits,	// if true then accept alignments in which multiple feature bits are set
					NumIncludeFiles,	// number of include region files
					pszIncludeFiles,	// biobed files containing regions to include - default is to exclude none
					NumExcludeFiles,	// number of exclude region files
					pszExcludeFiles,	// biobed file containing regions to exclude - default is to include all
					szSpeciesList,		// space or comma separated list of species, priority ordered
					iMinCoreLen,		// minimum core length required
					iMaxCoreLen,		// maximum core length required
					NumIncludeChroms,	// number of chromosomes explicitly defined to be included
					pszIncludeChroms,	// ptr to array of reg expressions defining chroms to include - overides exclude
					NumExcludeChroms,	// number of chromosomes explicitly defined to be excluded
					pszExcludeChroms,	// ptr to array of reg expressions defining chroms to include
					szScoreMatrix);		// scoring matrix to use
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

char *
ProcMode2Txt(etProcMode ProcMode)
{
switch(ProcMode) {
	case eProcStandard:
		return((char *)"ref loci + rel alignment distribution");
	default:
		break;
	}
return((char *)"Unsupported");
}

char *
TrimWhitespace(char *pTxt)
{
char *pStart;
char Chr;
	// strip leading whitespace
while(Chr = *pTxt++)
	if(!isspace(Chr))
			break;
if(Chr == '\0')					// empty line?
	return(pTxt-1);
pStart = pTxt-1;
while(Chr = *pTxt)			// fast forward to line terminator
	pTxt++;
pTxt-=1;
while(Chr = *pTxt--)
	if(!isspace(Chr))
		break;
pTxt[2] = '\0';
return(pStart);
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

int
LoadScoreMatrix(char *pszScoreMatrix,tsProcParams *pProcParams)
{
int Rslt;
CCSVFile *pCSV;
int NumLines;
int NumFields;

if((pCSV = new CCSVFile())==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"\nUnable to instantiate CCSVFile");
	CleanupResources(pProcParams);
	return(eBSFerrObj);
	}

if((Rslt = pCSV->Open(pszScoreMatrix))!= eBSFSuccess)
	{
	while(pCSV->NumErrMsgs())
			gDiagnostics.DiagOut(eDLFatal,gszProcName,pCSV->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"\nUnable to open loci file '%s'",pszScoreMatrix);
	delete pCSV;
	return(Rslt);
	}

NumLines = 0;
while((Rslt=pCSV->NextLine()) > 0)	// onto next line containing fields
	{
	NumFields = pCSV->GetCurFields();
	if(NumLines == 0)
		{
		if(NumFields < 3)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expected at least 3 fields in '%s', GetCurFields() returned '%d'",pszScoreMatrix,NumFields);
			delete pCSV;
			return(eBSFerrFieldCnt);
			}

		pCSV->GetInt(1,&pProcParams->MatrixScore.Unaligned);
		pCSV->GetInt(2,&pProcParams->MatrixScore.InDelAffine);
		pCSV->GetInt(3,&pProcParams->MatrixScore.InDelExtn);
		}
	else
		{
		if(NumFields < 5)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expected at least 5 fields in '%s', GetCurFields() returned '%d'",pszScoreMatrix,NumFields);
			delete pCSV;
			return(eBSFerrFieldCnt);
			}

		pCSV->GetInt(1,&pProcParams->MatrixScore.Score[NumLines-1][0]);
		pCSV->GetInt(2,&pProcParams->MatrixScore.Score[NumLines-1][1]);
		pCSV->GetInt(3,&pProcParams->MatrixScore.Score[NumLines-1][2]);
		pCSV->GetInt(4,&pProcParams->MatrixScore.Score[NumLines-1][3]);
		pCSV->GetInt(5,&pProcParams->MatrixScore.Score[NumLines-1][4]);
		}
	NumLines += 1;
	}
delete pCSV;
return(eBSFSuccess);
}


int
Process(etProcMode ProcMode,				// processing mode
		etLociFileType LociFileType,		// expected loci file type
				 char *pszInputFile,		// core loci (CSV or Biobed .bsb) file to process
		 		char *pszFilterRefIDFile, // exclude any RefIDs in this filter file
				 char *pszMAFFile,			// multiple alignment file
					char *pszRsltsFile,		// where to write out results
					char *pszSummaryFile,	// where to write out summary
					char *pszFeaturesFile,	// region characterisation biobed file
					int NumDistSegs,		// number of match distribution profile segments
					int RegLen,				// regulatory region length - up/dn stream of 5/3' 
					bool bMultipleFeatBits,	// if true then accept alignments in which multiple feature bits are set
					int NumIncludeFiles,	// number of include region files
					char **ppszIncludeFiles,// biobed files containing regions to include - default is to exclude none
					int NumExcludeFiles,	// number of exclude region files
					char **ppszExcludeFiles,// biobed file containing regions to exclude - default is to include all
					char *pszSpeciesList,	// space or comma separated list of species, priority ordered
					int	MinCoreLen,		// minimum core length required
					int	MaxCoreLen,		// maximum core length required
					int NumIncludeChroms,	// number of chromosomes explicitly defined to be included
					char **ppszIncludeChroms,	// ptr to array of reg expressions defining chroms to include - overides exclude
					int NumExcludeChroms,	// number of chromosomes explicitly defined to be excluded
					char **ppszExcludeChroms,	// ptr to array of reg expressions defining chroms to include
					char *pszScoreMatrix)			// scoring matrix file to use
{
int Rslt;
int Idx;


tsProcParams ProcParams;
char szCSVSpecies[512];				// to hold comma separated species list

#ifndef _WIN32
int RegErr;				// regular expression parsing error
char szRegErr[128];			// to hold RegErr as textual representation ==> regerror();
#endif

memset(&ProcParams,0,sizeof(tsProcParams));

ProcParams.hRsltsFile = -1;
ProcParams.hSummaryFile = -1;
ProcParams.ProcMode = ProcMode;
ProcParams.LociFileType = LociFileType;
ProcParams.MinCoreLen = MinCoreLen;
ProcParams.MaxCoreLen = MaxCoreLen;

if(pszScoreMatrix != NULL && pszScoreMatrix[0] != '\0')
	{
	if((Rslt = LoadScoreMatrix(pszScoreMatrix,&ProcParams))!=eBSFSuccess)
		return(Rslt);
	}
else
	{
	// default scoring is same as blastz scoring matrix
	ProcParams.MatrixScore.Unaligned = cUnaligned;
	ProcParams.MatrixScore.InDelAffine = cInDelAffine;
	ProcParams.MatrixScore.InDelExtn = cInDelExtn;
	ProcParams.MatrixScore.Score[eBaseA][eBaseA] = cScoreAA;
	ProcParams.MatrixScore.Score[eBaseA][eBaseC] = cScoreAC;
	ProcParams.MatrixScore.Score[eBaseA][eBaseG] = cScoreAG;
	ProcParams.MatrixScore.Score[eBaseA][eBaseT] = cScoreAT;
	ProcParams.MatrixScore.Score[eBaseA][eBaseN] = cScoreAN;
												   
	ProcParams.MatrixScore.Score[eBaseC][eBaseA] = cScoreCA;
	ProcParams.MatrixScore.Score[eBaseC][eBaseC] = cScoreCC;
	ProcParams.MatrixScore.Score[eBaseC][eBaseG] = cScoreCG;
	ProcParams.MatrixScore.Score[eBaseC][eBaseT] = cScoreCT;
	ProcParams.MatrixScore.Score[eBaseC][eBaseN] = cScoreCN;
												   
	ProcParams.MatrixScore.Score[eBaseG][eBaseA] = cScoreGA;
	ProcParams.MatrixScore.Score[eBaseG][eBaseC] = cScoreGC;
	ProcParams.MatrixScore.Score[eBaseG][eBaseG] = cScoreGG;
	ProcParams.MatrixScore.Score[eBaseG][eBaseT] = cScoreGT;
	ProcParams.MatrixScore.Score[eBaseG][eBaseN] = cScoreGN;
												   
	ProcParams.MatrixScore.Score[eBaseT][eBaseA] = cScoreTA;
	ProcParams.MatrixScore.Score[eBaseT][eBaseC] = cScoreTC;
	ProcParams.MatrixScore.Score[eBaseT][eBaseG] = cScoreTG;
	ProcParams.MatrixScore.Score[eBaseT][eBaseT] = cScoreTT;
	ProcParams.MatrixScore.Score[eBaseT][eBaseN] = cScoreTN;
												   
	ProcParams.MatrixScore.Score[eBaseN][eBaseA] = cScoreNA;
	ProcParams.MatrixScore.Score[eBaseN][eBaseC] = cScoreNC;
	ProcParams.MatrixScore.Score[eBaseN][eBaseG] = cScoreNG;
	ProcParams.MatrixScore.Score[eBaseN][eBaseT] = cScoreNT;
	ProcParams.MatrixScore.Score[eBaseN][eBaseN] = cScoreNN;
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Matrix Scores - Unaligned: %d, InDelAffine: %d, InDelExtn: %d",
	ProcParams.MatrixScore.Unaligned,ProcParams.MatrixScore.InDelAffine,ProcParams.MatrixScore.InDelExtn);

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Matrix Scores -  A\tC\tG\tT\tN");
for(Idx = 0; Idx < 5; Idx++)
	{	
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"                 %c\t%d\t%d\t%d\t%d\t%d",
			CSeqTrans::MapBase2Ascii(Idx),
				ProcParams.MatrixScore.Score[Idx][0],
				ProcParams.MatrixScore.Score[Idx][1],
				ProcParams.MatrixScore.Score[Idx][2],
				ProcParams.MatrixScore.Score[Idx][3],
				ProcParams.MatrixScore.Score[Idx][4]);
	}

ProcParams.bMultipleFeatBits = bMultipleFeatBits;	// if false then stats only generated if a single feature bit is set - e.g if both exons and introns overlapped then no stat generated
ProcParams.NumDistSegs = NumDistSegs;

strcpy(ProcParams.szLociFile,pszInputFile);
strcpy(ProcParams.szRsltsFile,pszRsltsFile);
if(pszFeaturesFile != NULL && pszFeaturesFile[0] != '\0')
	strcpy(ProcParams.szFeaturesFile,pszFeaturesFile);
if(pszSummaryFile != NULL && pszSummaryFile[0] != '\0')
	strcpy(ProcParams.szSummaryFile,pszSummaryFile);

ProcParams.NumIncludeChroms = NumIncludeChroms;		// number of chromosomes explicitly defined to be included
ProcParams.ppszIncludeChroms = ppszIncludeChroms;	// ptr to array of reg expressions defining chroms to include - overides exclude
ProcParams.NumExcludeChroms = NumExcludeChroms;		// number of chromosomes explicitly defined to be excluded
ProcParams.ppszExcludeChroms = ppszExcludeChroms;	// ptr to array of reg expressions defining chroms to include

if(pszFilterRefIDFile != NULL && pszFilterRefIDFile[0] != '\0')
	{
	ProcParams.pFilterRefIDs = new CFilterRefIDs;
	if((Rslt=ProcParams.pFilterRefIDs->Open(pszFilterRefIDFile)) < 0)
		{
		while(ProcParams.pFilterRefIDs->NumErrMsgs())
			gDiagnostics.DiagOut(eDLFatal,gszProcName,ProcParams.pFilterRefIDs->GetErrMsg());
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open/parse input filter RefID file %s",pszFilterRefIDFile);
		return(Rslt);
		}
	}

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
	CleanupResources(&ProcParams);
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
	CleanupResources(&ProcParams);
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
		CleanupResources(&ProcParams);
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
		CleanupResources(&ProcParams);
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to process exclude chrom '%s' error: %s",ppszExcludeChroms[Idx],szRegErr);
		return(eBSFerrMem);
		}
	}
#endif


// parse out species list
ProcParams.MaxNumSpecies = ParseNumSpecies(pszSpeciesList,&ProcParams); // only ever interested in ref + rel
ProcParams.MinNumSpecies = 1;	// ref species must always be present in block
szCSVSpecies[0]='\0';
for(Idx = 0; Idx < ProcParams.MaxNumSpecies; Idx++)
	{
	if(Idx > 0)
		strcat(szCSVSpecies,",");
	strcat(szCSVSpecies,ProcParams.szSpecies[Idx]);
	}
ProcParams.pszSpeciesList = szCSVSpecies;

ProcParams.Regions = 1;
ProcParams.UpDnStreamLen = 0;
ProcParams.pBEDfeatures = NULL;

if(pszFeaturesFile != NULL && pszFeaturesFile[0] != '\0')
	{
	if((ProcParams.pBEDfeatures = OpenBedfile(pszFeaturesFile,true))==NULL)
		{
		CleanupResources(&ProcParams);
		return(eBSFerrObj);
		}
	if(!ProcParams.pBEDfeatures->ContainsGeneDetail())
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Genomic feature file '%s' does not contain gene detail",pszFeaturesFile);
		CleanupResources(&ProcParams);
		return(eBSFerrFileType);
		}

	if(pszSummaryFile != NULL && pszSummaryFile[0] != '\0')
		{
#ifdef _WIN32
		ProcParams.hSummaryFile = open(pszSummaryFile, O_CREATETRUNC);
#else
		if((ProcParams.hSummaryFile = open64(pszSummaryFile,O_RDWR | O_CREAT,S_IREAD | S_IWRITE))!=-1)
			if(ftruncate(ProcParams.hSummaryFile,0)){};
#endif

		if(ProcParams.hSummaryFile == -1)					// check if file open succeeded
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create/open summary file '%s' - %s",pszSummaryFile,strerror(errno));
			CleanupResources(&ProcParams);
			return(eBSFerrCreateFile);
			}

		ProcParams.Regions = cNumCharRegs;	
		ProcParams.UpDnStreamLen = RegLen;
		if((ProcParams.pCntStepCnts = new int[ProcParams.Regions * cLenRanges])==NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate memory (%d bytes) for holding alignment statistics",sizeof(int) * ProcParams.Regions * cLenRanges);
			CleanupResources(&ProcParams);
			return(eBSFerrMem);
			}
		memset(ProcParams.pCntStepCnts,0,ProcParams.Regions * cLenRanges * sizeof(int));
		}
	}

for(Idx=0;Idx<NumIncludeFiles; Idx++)
	{
	if((ProcParams.pIncludes[Idx] = OpenBedfile(ppszIncludeFiles[Idx],true))==NULL)
		{
		CleanupResources(&ProcParams);
		return(eBSFerrObj);
		}
	ProcParams.NumIncludes++;
	}

for(Idx=0;Idx<NumExcludeFiles; Idx++)
	{
	if((ProcParams.pExcludes[Idx] = OpenBedfile(ppszExcludeFiles[Idx],true))==NULL)
		{
		CleanupResources(&ProcParams);
		return(eBSFerrObj);
		}
	ProcParams.NumExcludes++;
	}


switch(LociFileType) {
	case eLociFileTypeCSV:
		if((ProcParams.pCSVloci = new CCSVFile())==NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"\nUnable to instantiate CCSVFile");
			CleanupResources(&ProcParams);
			return(eBSFerrObj);
			}

		if((Rslt = ProcParams.pCSVloci->Open(pszInputFile))!= eBSFSuccess)
			{
			while(ProcParams.pCSVloci->NumErrMsgs())
				gDiagnostics.DiagOut(eDLFatal,gszProcName,ProcParams.pCSVloci->GetErrMsg());
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"\nUnable to open loci file '%s'",pszInputFile);
			CleanupResources(&ProcParams);
			return(Rslt);
			}

		if((Rslt = LoadCSVloci(&ProcParams)) < 1)
			{
			CleanupResources(&ProcParams);
			return(Rslt);
			}
		delete ProcParams.pCSVloci;
		ProcParams.pCSVloci = NULL;
		break;

	case eLociFileTypeUCSCBED:
		if((Rslt = LoadUCSCBEDloci(pszInputFile,&ProcParams)) < 1)
			{
			CleanupResources(&ProcParams);
			return(Rslt);
			}
		break;

	case eLociFileTypeBIOBED:
		if((ProcParams.pBEDloci = OpenBedfile(pszInputFile,true))==NULL)
			{
			CleanupResources(&ProcParams);
			return(eBSFerrObj);
			}

		if((Rslt = LoadBIOBEDloci(&ProcParams)) < 1)
			{
			CleanupResources(&ProcParams);
			return(Rslt);
			}
		delete ProcParams.pBEDloci;
		ProcParams.pBEDloci = NULL;
		break;
	}

// determine max aligned sequence length for any single species which can be handled
ProcParams.MaxSeqAlignLen = cMAtotMaxSeqAlignLen/ProcParams.MaxNumSpecies;
for(Idx = 0; Idx < ProcParams.MaxNumSpecies; Idx++)
	{
	if((ProcParams.pSeqs[Idx] = new unsigned char [ProcParams.MaxSeqAlignLen])==NULL)
		{
	    CleanupResources(&ProcParams);
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate memory (%d bytes) for holding species sequences",ProcParams.MaxSeqAlignLen);
		return(eBSFerrMem);
		}
	}

#ifdef _WIN32
	ProcParams.hRsltsFile = open(pszRsltsFile, O_CREATETRUNC);
#else
     if((ProcParams.hRsltsFile = open64(pszRsltsFile,O_RDWR | O_CREAT,S_IREAD | S_IWRITE))!=-1)
		 if(ftruncate(ProcParams.hRsltsFile,0)){};
#endif

if(ProcParams.hRsltsFile == -1)					// check if file open succeeded
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create/open output file '%s' - %s",pszRsltsFile,strerror(errno));
	CleanupResources(&ProcParams);
	return(eBSFerrCreateFile);
	}

Rslt = ProcessAlignments(pszMAFFile,	// source bioseq multialignment file
				  &ProcParams); // processing parameters

OutputSummary(&ProcParams);

CleanupResources(&ProcParams);
return(Rslt);
}


int
AddChrom(char *pszChrom,tsProcParams *pParams)
{
tsCoreChrom *pChrom;
int Idx;
UINT16 Hash;
Hash = GenNameHash(pszChrom);

if(pParams->pCoreChroms != NULL && pParams->NumCoreChroms)
	{
	if(pParams->CachCoreChromID > 0)
		{
		pChrom = &pParams->pCoreChroms[pParams->CachCoreChromID-1];
		if((pChrom->Hash == Hash) && !stricmp(pszChrom,pChrom->szChrom))
			return(pChrom->ChromID);
		}

	pChrom = pParams->pCoreChroms;
	for(Idx = 0; Idx < pParams->NumCoreChroms; Idx++,pChrom++)
		{
		if((pChrom->Hash == Hash) && !stricmp(pszChrom,pChrom->szChrom))
			return(pParams->CachCoreChromID = pChrom->ChromID);
		}
	}

if(pParams->pCoreChroms == NULL || pParams->NumCoreChroms == pParams->MaxCoreChroms)
	{
	if((pChrom = new tsCoreChrom[pParams->MaxCoreChroms + cCoreChromAllocChunk])==NULL)
		{
		return(eBSFerrMem);
		}
	if(pParams->pCoreChroms != NULL)
		{
		memcpy(pChrom,pParams->pCoreChroms,sizeof(tsCoreChrom) * pParams->NumCoreChroms);
		pParams->pCoreChroms = pChrom;
		pParams->MaxCoreChroms += cCoreChromAllocChunk;
		}
	else
		{
		pParams->MaxCoreChroms = cCoreChromAllocChunk;
		pParams->NumCoreLocs = 0;
		}
	pParams->pCoreChroms = pChrom;
	}
pChrom = &pParams->pCoreChroms[pParams->NumCoreChroms++];
pChrom->ChromID = pParams->NumCoreChroms;
pChrom->Hash = Hash;
pChrom->pCoreLoci = NULL;
strcpy(pChrom->szChrom,pszChrom);
return(pParams->CachCoreChromID = pChrom->ChromID);
}


int
AddLoci(int RefID,char *pszChrom,char Strand,int StartLoci,int Len,tsProcParams *pParams)
{
tsCoreLoci *pNxtCore;
int ChromID;

if(Len < pParams->MinCoreLen || Len > pParams->MaxCoreLen)
	return(0);

if((ChromID = AddChrom(pszChrom,pParams)) < 1)
	return(ChromID);

if(pParams->pCoreLocs == NULL || (pParams->NumCoreLocs + 1) == pParams->MaxCoreLocs)
	{
	if((pNxtCore = new tsCoreLoci[pParams->MaxCoreLocs + cCoreLociAllocChunk])==NULL)
		{
		return(eBSFerrMem);
		}
	if(pParams->pCoreLocs != NULL)
		{
		if(pParams->NumCoreLocs)
			memcpy(pNxtCore,pParams->pCoreLocs,sizeof(tsCoreLoci) * pParams->NumCoreLocs);
		delete pParams->pCoreLocs;
		}
	else
		{
		pParams->MaxCoreLocs = 0;
		pParams->NumCoreLocs = 0;
		}
	pParams->MaxCoreLocs += cCoreLociAllocChunk;
	pParams->pCoreLocs = pNxtCore;
	}
pNxtCore = &pParams->pCoreLocs[pParams->NumCoreLocs++];
pNxtCore->RefID = RefID;
pNxtCore->ChromID = ChromID;
pNxtCore->Strand = Strand;
pNxtCore->StartLoci = StartLoci;
pNxtCore->EndLoci = StartLoci + Len - 1;
return(pParams->NumCoreLocs);
}

int
LoadCSVloci(tsProcParams *pParams)
{
int NumFields;
int Rslt;
int RefID;
char *pszChrom;
int StartLoci;
int EndLoci;
int Len;
int NumElsRead;		// number of elements before filtering
int NumElsAccepted;	// number of elements accepted after filtering
int NumFiltRefIDs;	// number of elements filtered out because of matching RefID
int NumFiltLen;		// number of elements filtered out because of length

CCSVFile *pCSV = pParams->pCSVloci;

NumElsRead =0;		// number of elements before filtering
NumElsAccepted =0;	// number of elements accepted after filtering
NumFiltRefIDs =0;	// number of elements filtered out because of matching RefID
NumFiltLen=0;		// number of elements filtered out because of length

while((Rslt=pCSV->NextLine()) > 0)	// onto next line containing fields
	{
	NumFields = pCSV->GetCurFields();
	if(NumFields < 7)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expected at least 7 fields in '%s', GetCurFields() returned '%d' at input line '%d'",pParams->szLociFile,NumFields,pCSV->GetLineNumber());
		return(eBSFerrFieldCnt);
		}
	NumElsRead += 1;

	pCSV->GetInt(1,&RefID);
	if(pParams->pFilterRefIDs != NULL && pParams->pFilterRefIDs->Locate(RefID))
		{	
		NumFiltRefIDs += 1;
		continue;
		}

	pCSV->GetInt(7,&Len);
	if(Len < pParams->MinCoreLen || Len > pParams->MaxCoreLen)
		{
		NumFiltLen += 1;
		continue;
		}

	pCSV->GetText(4,&pszChrom);
	pCSV->GetInt(5,&StartLoci);
	pCSV->GetInt(6,&EndLoci);

	if((Rslt = AddLoci(RefID,pszChrom,'+',StartLoci,Len,pParams)) < 0)
		return(Rslt);
	NumElsAccepted += 1;
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"%s: Read: %d, Accepted: %d, FiltRefIDs: %d, FiltLen: %d",
	pParams->szLociFile,NumElsRead,NumElsAccepted,NumFiltRefIDs,NumFiltLen);

return(ProcessLoci(pParams));
}

int
LoadBIOBEDloci(tsProcParams *pParams)
{
int Rslt;
CBEDfile *pBED = pParams->pBEDloci;
int CurFeatureID;
int NumFilteredLoci;
char szChrom[cMaxDatasetSpeciesChrom];
int StartLoci;
int EndLoci;
int FeatureLen;
char Strand;

CurFeatureID = 0;
NumFilteredLoci = 0;
while((CurFeatureID = pBED->GetNextFeatureID(CurFeatureID)) > 0)
	{
	pBED->GetFeature(CurFeatureID,NULL,szChrom,&StartLoci,&EndLoci,0,&Strand);
	FeatureLen = (EndLoci - StartLoci) + 1;
	if(FeatureLen < pParams->MinCoreLen || FeatureLen > pParams->MaxCoreLen)
		continue;

	if((Rslt = AddLoci(NumFilteredLoci+1,szChrom,'+',StartLoci,FeatureLen,pParams)) < 0)
		return(Rslt);

	NumFilteredLoci += 1;
	}

return(ProcessLoci(pParams));
}


int
LoadUCSCBEDloci(char *pszFileName,tsProcParams *pParams)
{
FILE *pBEDStream;
int LineNum;
char szLineBuff[cLineBuffLen];
char szChrom[51];
int StartLoci;
int EndLoci;
char szName[cMaxFeatNameLen+1];
int SupInfoStart;
char szSuppInfo[cLineBuffLen];
int Score;
char Strand;
int Cnt;
int FeatNum;
int FeatureLen;
int Rslt;
char *pTxt;

if(pszFileName == NULL || *pszFileName == '\0')
	return(eBSFerrParams);
if((pBEDStream = fopen(pszFileName,"r"))==NULL)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Unable to open BED format file %s error: %s",pszFileName,strerror(errno));
	return(eBSFerrOpnFile);
	}
LineNum = 0;
FeatNum = 0;

Rslt = eBSFerrParse; // assume the worst..., no features parsed
while(fgets(szLineBuff,sizeof(szLineBuff),pBEDStream)!= NULL)
	{
	pTxt = TrimWhitespace(szLineBuff);
	if(*pTxt=='\0')	// simply slough lines which were just whitespace
		continue;
	Cnt = sscanf(pTxt," %50s %d %d %50s %d %c %n",
			szChrom,&StartLoci,&EndLoci,szName,&Score,&Strand,&SupInfoStart);
	if(Cnt < 3)		// must be rubbish or a comment on this line
		{
		if(!FeatNum)	
			continue;	// no features processed, hope just a comment and keep trying!
		Rslt = eBSFerrParse; // if at least one feature processed then treat rubbish as error
		break;
		}
		
	switch(Cnt) {
		case 3:		// name, score,strand,suppinfo missing 
			sprintf(szName,"feat%d",FeatNum+1);
			Score = 0; Strand = '+'; szSuppInfo[0] = '\0';
			break;
		case 4:		// score,strand,suppinfo missing 
			Score = 0; Strand = '+'; szSuppInfo[0] = '\0';
			break;
		case 5:		// strand,suppinfo missing
			Strand = '+'; szSuppInfo[0] = '\0';
			break;
		default:	// suppinfo may be present but of no interest!
			break;
		}

	FeatureLen = EndLoci - StartLoci;
	if(FeatureLen < pParams->MinCoreLen || FeatureLen > pParams->MaxCoreLen)
		continue;
			
	if(Strand == '.')
		Strand = '+';

	if((Rslt = AddLoci(++FeatNum,szChrom,Strand,StartLoci,FeatureLen,pParams)) < 0)
		{
		fclose(pBEDStream);
		return(Rslt);
		}
	}
fclose(pBEDStream);
return(ProcessLoci(pParams));
}

int
ProcessLoci(tsProcParams *pParams)
{
int CurChromID;
int ChromIdx;
int LociIdx;
tsCoreLoci *pLoci;
tsCoreChrom *pChrom;

if(pParams->pCoreLocs == NULL || pParams->NumCoreLocs < 1)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"No accepted loci to process in '%s'",pParams->szLociFile);
	return(-1);
	}

// sort cores by chrom.start.end
if(pParams->NumCoreLocs > 1)
	qsort(pParams->pCoreLocs,pParams->NumCoreLocs,sizeof(tsCoreLoci),CompareCoreEls);

// add a sentenil element as last element (allocation ensures at least 1 free element!) 
pParams->pCoreLocs[pParams->NumCoreLocs].ChromID = -1;

// process core chroms setting ptrs to first core in each chrom
pLoci = pParams->pCoreLocs;
pChrom = pParams->pCoreChroms;
CurChromID = -1;
for(ChromIdx = LociIdx = 0; LociIdx < pParams->NumCoreLocs; LociIdx++,pLoci++)
	{
	if(pLoci->ChromID != CurChromID)
		{
		pChrom->pCoreLoci = pLoci;
		CurChromID = pLoci->ChromID;
		pChrom += 1;
		}
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Parsed and accepted %d loci for processing from '%s'",LociIdx,pParams->szLociFile);
return(LociIdx);
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

// to be excluded?
for(Idx = 0; Idx < pProcParams->NumExcludeChroms; Idx++)
#ifdef _WIN32	
	if(pProcParams->ExcludeChromsRE[Idx]->Match(szSpeciesChrom,&mc))
#else
	if(!regexec(&pProcParams->ExcludeChromsRE[Idx],szSpeciesChrom,1,&mc,0))
#endif
		return(AddExcludeHistory(SpeciesID,ChromID,true));

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


return(AddExcludeHistory(SpeciesID,ChromID,pProcParams->NumIncludeChroms ? true : false));
}

// OpenBedfile
// Attempts to open specified bedfile
// Returns ptr to opened bedfile or NULL
CBEDfile *
OpenBedfile(char *pToOpen,bool bReportErrs)
{
int Rslt;
CBEDfile *pBed;
if(pToOpen != NULL && pToOpen[0] != '\0')
	{
	if((pBed = (CBEDfile *)new CBEDfile())==NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"\nUnable to instantiate CBEDfile '%s'",pToOpen);
		return(NULL);
		}

	if((Rslt = pBed->Open(pToOpen,eBTAnyBed))!=eBSFSuccess)
		{
		if(bReportErrs)
			{
			while(pBed->NumErrMsgs())
				gDiagnostics.DiagOut(eDLFatal,gszProcName,pBed->GetErrMsg());
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open biobed file %s",pToOpen);
			}
		delete pBed;
		return(NULL);
		}
	return(pBed);
	}
return(NULL);
}

//CleanupResources
//Closes and deletes all opened resources including files
bool
CleanupResources(tsProcParams *pProcParams)
{
int Idx;

if(pProcParams->hRsltsFile != -1)
	{
	close(pProcParams->hRsltsFile);
	pProcParams->hRsltsFile = -1;
	}

if(pProcParams->hSummaryFile != -1)
	{
	close(pProcParams->hSummaryFile);
	pProcParams->hSummaryFile = -1;
	}

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

#ifdef _WIN32
for(Idx=0;Idx < pProcParams->NumIncludeChroms;Idx++)
	if(pProcParams->IncludeChromsRE[Idx] != NULL)
		{
		delete pProcParams->IncludeChromsRE[Idx];
		pProcParams->IncludeChromsRE[Idx] = NULL;
		}
for(Idx=0;Idx < pProcParams->NumExcludeChroms;Idx++)
	if(pProcParams->ExcludeChromsRE[Idx] != NULL)
		{
		delete pProcParams->ExcludeChromsRE[Idx];
		pProcParams->ExcludeChromsRE[Idx] = NULL;
		}
pProcParams->NumIncludeChroms = 0;
pProcParams->NumExcludeChroms = 0;
#endif

if(pProcParams->pFilterRefIDs != NULL)
	delete pProcParams->pFilterRefIDs;
pProcParams->pFilterRefIDs = NULL;

if(pProcParams->pBEDfeatures != NULL)
	delete pProcParams->pBEDfeatures;
pProcParams->pBEDfeatures = NULL;

if(pProcParams->pBEDloci != NULL)
	delete pProcParams->pBEDloci;
pProcParams->pBEDloci = NULL;
if(pProcParams->pCSVloci != NULL)
	delete pProcParams->pCSVloci;
pProcParams->pCSVloci = NULL;
if(pProcParams->pCoreLocs != NULL)
	delete pProcParams->pCoreLocs;
pProcParams->pCoreLocs = NULL;
pProcParams->MaxCoreLocs = 0;
pProcParams->NumCoreLocs = 0;
if(pProcParams->pCoreChroms != NULL)
	delete pProcParams->pCoreChroms;
pProcParams->pCoreChroms = NULL;
pProcParams->MaxCoreChroms = 0;
pProcParams->NumCoreChroms = 0;
pProcParams->CachCoreChromID = 0;

for(Idx = 0; Idx < pProcParams->MaxNumSpecies; Idx++)
	{
	if(pProcParams->pSeqs[Idx] != NULL)
		{
		delete pProcParams->pSeqs[Idx];
		pProcParams->pSeqs[Idx] = NULL;
		}
	}
pProcParams->SeqLen = 0;

if(pProcParams->pCntStepCnts != NULL)
	{
	delete pProcParams->pCntStepCnts;
	pProcParams->pCntStepCnts = NULL;
	}
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


int 
ProcessAlignments(char *pszMAF,			 // source bioseq multialignment file
				  tsProcParams *pProcParams) // processing parameters
{
int RefSpeciesID;
int RelSpeciesID;
int RefChromID;
int PrevRefChromID;
int PrevDispRefChromID;
char *pszRefChrom;
int StartLoci;
int	EndLoci;
int SpeciesIDs[cMaxAlignedSpecies];
int Rslt;
int RefAlignLen;
CMAlignFile *pAlignments;
int Idx;
int ExpAlignLen;
int NormAlignLen;

if((pAlignments = new CMAlignFile())==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create new instance of CMAlignFile");
	return(eBSFerrObj);
	}

if((Rslt=pAlignments->Open(pszMAF))!=eBSFSuccess)
	{
	while(pAlignments->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,pAlignments->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open MAF file '%s'",pszMAF);
	return(Rslt);
	}

// ensure all species are represented in multispecies alignment file plus get their species identifiers
for(Idx = 0; Idx < pProcParams->MaxNumSpecies; Idx++)
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
	}

// iterate over core loci
Rslt = eBSFSuccess;
PrevRefChromID = 0;
PrevDispRefChromID = 0;
RefSpeciesID = SpeciesIDs[0];
RelSpeciesID = SpeciesIDs[1];
pProcParams->MaxAlignIdxSpecies = 2;
tsCoreLoci *pCoreEl = pProcParams->pCoreLocs;
for(Idx = 0; Idx < pProcParams->NumCoreLocs; Idx++,pCoreEl++)
	{
	StartLoci = pCoreEl->StartLoci;
	EndLoci = pCoreEl->EndLoci;
	ExpAlignLen = 1 + EndLoci - StartLoci;
	if(ExpAlignLen < pProcParams->MinCoreLen)
		continue;

	pszRefChrom = pProcParams->pCoreChroms[pCoreEl->ChromID-1].szChrom;
	RefChromID = pAlignments->LocateChromID(RefSpeciesID,pszRefChrom);
	if(RefChromID < 1)		// if core chromosome not in alignment
		continue;
	if(ExcludeThisChrom(pAlignments,RefSpeciesID,RefChromID,pProcParams))
		continue;

	RefAlignLen = pAlignments->LoadAlignment(RefChromID,StartLoci,EndLoci,RelSpeciesID,pProcParams->MaxSeqAlignLen,pProcParams->pSeqs[0],pProcParams->pSeqs[1],eBaseUndef,eBaseUndef);
	if(RefAlignLen < 1 + pCoreEl->EndLoci - pCoreEl->StartLoci)
		gDiagnostics.DiagOut(eDLWarn,gszProcName,"Returned core element length of %d less than expected %d",RefAlignLen,1 + pCoreEl->EndLoci - pCoreEl->StartLoci);

	if(RefChromID != PrevRefChromID)
		{
		PrevRefChromID = RefChromID;
		strcpy(pProcParams->szRefChrom,pszRefChrom);
		pProcParams->RefChromID = RefChromID;
		if(pProcParams->pBEDfeatures != NULL)
			pProcParams->BEDChromID = pProcParams->pBEDfeatures->LocateChromIDbyName(pProcParams->szRefChrom);
		else
			pProcParams->BEDChromID = 0;
		}
	pProcParams->RefChromID = RefChromID;

		// here we need to normalise the alignments so that there will be no case of all InDels in
		// any column which can occur if none of the processed species is not the reference species
	NormAlignLen = NormaliseInDelColumns(pProcParams,RefAlignLen);

	if(PrevDispRefChromID != RefChromID)
		{
		PrevDispRefChromID = RefChromID;
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing chromosome %s",pszRefChrom);
		}

	if((Rslt = ProcLociSeq(pCoreEl,		// current core
			0,							// index into pProcParams->pSeqs[0][] at which loci mapped alignment starts 
		   NormAlignLen-1,				// index into pProcParams->pSeqs[0][] at which loci mapped alignment ends 
		   pszRefChrom,					// reference chromosome
		   StartLoci,					// start loci on chromosome
		   EndLoci,						// end loci on chromosome
		   pProcParams)) < 0)			// global processing parameters
		   break;
	}

delete pAlignments;
return(Rslt);
}


// NormaliseInDelColumns
// Because multialignments may be the result of merged alignments resulting in InDels being generated back into the 
// existing reference sequence then if only a subset of species are being processed there could be all InDels in any column
// of the subset sequences. 
// This function will delete all columns which only contain InDels or in which the reference species has an indel and rel species a undefined base
// Returns the subset sequence length after eBaseInDel columns have been deleted
int
NormaliseInDelColumns(tsProcParams *pProcParams,int AlignLen)
{
etSeqBase *pSrcRefSeq;
etSeqBase *pDstRefSeq;
etSeqBase *pSrcRelSeq;
etSeqBase *pDstRelSeq;
etSeqBase RefBase;
etSeqBase RelBase;
int SeqIdx;
int NormAlignLen=0;

pSrcRefSeq = pDstRefSeq = pProcParams->pSeqs[0];
pSrcRelSeq = pDstRelSeq = pProcParams->pSeqs[1];
for(SeqIdx = 0; SeqIdx < AlignLen; SeqIdx++,pSrcRefSeq++,pSrcRelSeq++)
	{
	RefBase = *pSrcRefSeq & ~cRptMskFlg;
	RelBase = *pSrcRelSeq & ~cRptMskFlg;
	if(RefBase == eBaseInDel && (RelBase == eBaseInDel || RelBase == eBaseUndef))
		continue;
	*pDstRefSeq++ = *pSrcRefSeq;
	*pDstRelSeq++ = *pSrcRelSeq;
	NormAlignLen++;
	}
return(NormAlignLen);
}




int
ProcLociSeq(tsCoreLoci *pCore,		// current core
			int StartSeqIdx,				// index into pProcParams->pSeqs[0][] at which loci mapped alignment starts 
		   int EndSeqIdx,				// index into pProcParams->pSeqs[0][] at which loci mapped alignment ends 
		   char *pszRefChrom,			// reference chromosome
		   int StartLoci,				// start loci on chromosome
		   int EndLoci,					// end loci on chromosome
		   tsProcParams *pProcParams) // global processing parameters
{
tsLenRangeClass *pRange;
int RefSeqLen;
int CoreLen;
int FeatureBits;
int BitMsk;
int SpliceSiteOverlaps;
int FeatIdx;
int RegionIdx;
int *pStep;
int GapLen;
int Score;
int RefScore;
int CoreOfs;
int OGMatches;
int OGMismatches;
int OGInDels;
int OGUnaligned;
etSeqBase OGBase;
etSeqBase RefBase;
etSeqBase *pVSeq;

tsDistSeg OGDistProfile[cMaxMatchDistSegments];


if(pProcParams->hRsltsFile <= 0)
	return(0);

RefSeqLen = 1 + EndLoci - StartLoci;
CoreLen = 1 + EndSeqIdx - StartSeqIdx;

if(pProcParams->pCntStepCnts != NULL)
	{
	if((pRange = GetLengthRangeClass(RefSeqLen))==NULL)
		{
		gDiagnostics.DiagOut(eDLWarn,gszProcName,"Unexpected length range classification error - length = %d",RefSeqLen);
		return(-1);
		}
	pStep = pProcParams->pCntStepCnts;
	pStep += (pRange->ID-1) * pProcParams->Regions;
	}
else
	{
	pStep = NULL;
	pRange = NULL;
	}
RegionIdx = 0;
FeatureBits = 0;
SpliceSiteOverlaps = 0;
if(pProcParams->pBEDfeatures != NULL)
	{
	if(pProcParams->BEDChromID > 0)
		FeatureBits = pProcParams->pBEDfeatures->GetFeatureBits(pProcParams->BEDChromID,StartLoci,EndLoci,cRegionFeatBits,pProcParams->UpDnStreamLen);
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
					return(0);		// although was sequence of interest, more than one feature bit so can't contribute to stats
				RegionIdx = FeatIdx;
				if(pProcParams->bMultipleFeatBits)	// if multiple features allowed then don't check for any additional
					break;
				}
			}
		}

	if(pStep != NULL)
		{
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
		}

	if(FeatureBits != 0)		// if not intergenic then check for splice sites
		SpliceSiteOverlaps = pProcParams->pBEDfeatures->GetSpliceSiteBits(pProcParams->BEDChromID,StartLoci,EndLoci,cMinSpliceOverlap);
	if(pStep != NULL)
		{
		if(SpliceSiteOverlaps & cIntronExonSpliceSite)
			pStep[7]++;
		if(SpliceSiteOverlaps & cExonIntronSpliceSite)
			pStep[8]++;
		}
	}

memset(OGDistProfile,0,sizeof(OGDistProfile));

OGUnaligned = 0;
OGMatches = 0;
OGMismatches = 0;
OGInDels = 0;
	
Score = 0;
GapLen = 0;
int CoreLoci = 0;
int CoreLenLeft = RefSeqLen;
int NumSegsLeft = pProcParams->NumDistSegs;
int SegLoci = RefSeqLen / NumSegsLeft--;
int OGSegIdx = 0;
for(CoreOfs = 0; CoreOfs < CoreLen ; CoreOfs++)
	{
	pVSeq = pProcParams->pSeqs[0];
	RefBase = pVSeq[CoreOfs+StartSeqIdx] & ~cRptMskFlg;
	pVSeq = pProcParams->pSeqs[1];
	OGBase = pVSeq[CoreOfs+StartSeqIdx] & ~cRptMskFlg;

	if(RefBase == eBaseInDel || OGBase == eBaseInDel)	// handle InDels first
		{
		if(!GapLen)			// if opening a gap
			Score += pProcParams->MatrixScore.InDelAffine;		
		else
			Score += pProcParams->MatrixScore.InDelExtn; // existing gap so extension score
		GapLen += 1;
		if(RefBase > eBaseN)
			continue;
		}
	else
		GapLen = 0;
	
	if(RefBase == eBaseUndef)
		{
		OGUnaligned += 1;
		OGDistProfile[OGSegIdx].Unaligned += 1;
		Score += pProcParams->MatrixScore.Unaligned;
		}
	else
		{
		if(RefBase == eBaseN)			// if ref base unknown then treat as as mismatch				
			{
			OGMismatches += 1;
			OGDistProfile[OGSegIdx].Mismatches += 1;
			}
		else
			{
			if(OGBase == RefBase)		// profile outgroup bases
				{
				OGMatches += 1;
				OGDistProfile[OGSegIdx].Matches += 1;
				}
			else
				{
				if(OGBase == eBaseUndef)
					{
					OGUnaligned += 1;
					OGDistProfile[OGSegIdx].Unaligned += 1;
					Score += pProcParams->MatrixScore.Unaligned;
					}
				else
					{
					if(OGBase == eBaseInDel)
						{
						OGInDels += 1;
						OGDistProfile[OGSegIdx].InDels += 1;
						}
					else
						{
						OGMismatches += 1;
						OGDistProfile[OGSegIdx].Mismatches += 1;
						}
					}
				}
			}
		}

	if(OGBase <= eBaseN)
		Score += pProcParams->MatrixScore.Score[RefBase][OGBase];
	
	CoreLoci += 1;
	CoreLenLeft -= 1;
	if(NumSegsLeft && CoreLoci >= SegLoci)
		{
		OGSegIdx += 1;
		SegLoci += CoreLenLeft / NumSegsLeft--;
		}
	}
	

// calc background score as sum of refbase scores
pVSeq = pProcParams->pSeqs[0];
RefScore = 0;
for(CoreOfs = 0; CoreOfs < CoreLen ; CoreOfs++)
	{
	RefBase = pVSeq[CoreOfs+StartSeqIdx] & ~cRptMskFlg;
	if(RefBase <= eBaseN)
		RefScore += pProcParams->MatrixScore.Score[RefBase][RefBase];
	// currently the ref score does not account for any InDels
	}
OutputLociCore(pCore,pProcParams->szRefChrom,StartLoci,EndLoci,FeatureBits | SpliceSiteOverlaps,OGUnaligned,OGMatches,OGMismatches,OGInDels,Score,RefScore,OGDistProfile,pProcParams);
return(0);
}


// GenNameHash
// Generates a 16bit hash on specified lowercased name
UINT16 
GenNameHash(char *pszName)
{
unsigned long hash = 5381;
char Chr;
while (Chr = *pszName++)
	hash = ((hash << 5) + hash) + tolower(Chr);
return ((UINT16)hash);
}


tsCoreChrom *
LocateCoreChrom(char *pszChrom,
				tsProcParams *pProcParams) // global processing parameters
{
tsCoreChrom *pChrom;
UINT16 Hash;
int Idx;
if(pProcParams->pCoreChroms == NULL || pProcParams->NumCoreChroms < 1)
	return(NULL);
Hash = GenNameHash(pszChrom);

// quick check against cached chrom name
if(pProcParams->CachCoreChromID > 0)
	{
	pChrom = &pProcParams->pCoreChroms[pProcParams->CachCoreChromID-1];
	if(Hash == pChrom->Hash && !stricmp(pChrom->szChrom,pszChrom))
		return(pChrom);
	}
pChrom = pProcParams->pCoreChroms;
for(Idx = 0; Idx < pProcParams->NumCoreChroms; Idx++,pChrom++)
	{
	if(Hash == pChrom->Hash && !stricmp(pChrom->szChrom,pszChrom))
		{
		pProcParams->CachCoreChromID = pChrom->ChromID;
		return(pChrom);
		}
	}
return(NULL);
}


tsCoreLoci AcceptAllLoci;

tsCoreLoci *
GetLociOverlaps(char *pszChrom,			// reference species chromosome
				int RefChromID,			// chrom identifier
				char RefStrand,			// strand
				int RefChromOfs,		// start loci
				int RefChromEndOfs,		// end loci
				tsCoreLoci *pPrevCoreLoci, //
				tsProcParams *pProcParams) // global processing parameters
{
tsCoreLoci *pLoci;

if(pPrevCoreLoci == NULL)
	return(GetFirstLociOverlaps(pszChrom,RefChromID,RefStrand,RefChromOfs,RefChromEndOfs,pProcParams));

pLoci = pPrevCoreLoci + 1;
while(pLoci->ChromID > 0 && pLoci->ChromID == pPrevCoreLoci->ChromID)
	{
	if(pLoci->StartLoci > RefChromEndOfs)
		break;

	if(pLoci->Strand == RefStrand && pLoci->StartLoci <= RefChromEndOfs && pLoci->EndLoci >= RefChromOfs)
		return(pLoci);
	pLoci += 1;
	}
return(NULL);
}

bool	
OutputLociCore(tsCoreLoci *pCore,		// current core
			   char *pszChrom, int ChromStartOffset, int ChromEndOffset, 
				int FeatureBits,		// feature bits over lapped
				int OGUnaligned,		// number of unaligned bases in outspecies
				int OGMatches,			// number of matching bases in outspecies
				int OGMismatches,		// number of mismatched bases in outspecies
				int OGInDels,			// number of InDels in outspecies
				int Score,				// alignment score
				int RefScore,			// background ref seq score
				tsDistSeg SegCnts[],	// array of segment profile counts	
				tsProcParams *pProcParams)
{
char szLineBuff[4096];
int Len;
if(pProcParams->hRsltsFile != -1)
	{
	Len = sprintf(szLineBuff,"%d,\"LociCore\",\"%s\",\"%s\",%d,%d,%d,\"%s\",%d",
		pCore->RefID,
		pProcParams->szSpecies[pProcParams->RefSpeciesIdx],
		pszChrom,ChromStartOffset,ChromEndOffset,1 + ChromEndOffset-ChromStartOffset,
		pProcParams->pszSpeciesList,FeatureBits & (cAnyFeatBits | cOverlaysSpliceSites));

	Len += sprintf(&szLineBuff[Len],",%d,%d,%d,%d,%d,%d",OGUnaligned,OGMatches,OGMismatches,OGInDels,Score,RefScore);
	for(int Idx = 0; Idx < pProcParams->NumDistSegs; Idx++)
		Len += sprintf(&szLineBuff[Len],",%d,%d,%d,%d",SegCnts[Idx].Matches,SegCnts[Idx].Mismatches,SegCnts[Idx].InDels,SegCnts[Idx].Unaligned);
	Len += sprintf(&szLineBuff[Len],"\n");
	CUtility::SafeWrite(pProcParams->hRsltsFile,szLineBuff,Len);
	}
return(true);
}

bool	
OutputSummary(tsProcParams *pProcParams)
{
tsLenRangeClass *pRange;
static bool bOutputHdrFirst = true;
char szLineBuff[cMaxReadLen+1];
int Idx;
int Steps;
int Instances;
int Len;
int *pStep;
if(pProcParams->hSummaryFile != -1 && pProcParams->pCntStepCnts != NULL)
	{
	pStep = pProcParams->pCntStepCnts;

	if(bOutputHdrFirst)
		{
		bOutputHdrFirst = false;
		Len = sprintf(szLineBuff,"\"LenRange\",\"Mismatches\",\"TotInstances\",\"IG\",\"US\",\"5'UTR\",\"CDS\",\"INTRON\",\"3'UTR\",\"DS\",\"5'ExSplice\",\"3'ExSplice\"");
		CUtility::SafeWrite(pProcParams->hSummaryFile,szLineBuff,Len);
		}
	pStep = pProcParams->pCntStepCnts;
	for(Idx = 0; Idx < cLenRanges; Idx++, pStep += pProcParams->Regions)
		{
		for(Instances = Steps = 0; Steps < (pProcParams->Regions - 2); Steps++)
			Instances += pStep[Steps];
		pRange = GetRangeClass(Idx+1);
		Len = sprintf(szLineBuff,"\n\"%s\",0,%d",
							pRange->pszDescr,Instances);
		for(Steps = 0; Steps < pProcParams->Regions; Steps++)
				Len += sprintf(&szLineBuff[Len],",%d",pStep[Steps]);
		CUtility::SafeWrite(pProcParams->hSummaryFile,szLineBuff,Len);
		}
	}
return(true);
}



tsCoreLoci *
GetFirstLociOverlaps(char *pszChrom,	// reference species chromosome
   				int RefChromID,			// chrom identifier
				char RefStrand,			// strand
				int RefChromOfs,		// start loci
				int RefChromEndOfs,		// end loci
				tsProcParams *pProcParams) // global processing parameters
{
tsCoreChrom *pChrom;
tsCoreLoci *pLoci;

if((pChrom = LocateCoreChrom(pszChrom,pProcParams))==NULL)
   return(NULL);
if((pLoci = pChrom->pCoreLoci)==NULL)
	return(NULL);
while(pLoci->ChromID > 0 && pLoci->ChromID == pChrom->ChromID)
	{
	if(pLoci->StartLoci > RefChromEndOfs)
		break;
	if(pLoci->Strand == RefStrand && pLoci->StartLoci <= RefChromEndOfs && pLoci->EndLoci >= RefChromOfs)
		return(pLoci);
	pLoci += 1;
	}
return(NULL);
}



// CompareCoreEls
// Used to sort core loci elements 
int 
CompareCoreEls( const void *arg1, const void *arg2)
{
tsCoreLoci *pEl1 = (tsCoreLoci *)arg1;
tsCoreLoci *pEl2 = (tsCoreLoci *)arg2;
if(pEl1->ChromID == pEl2->ChromID)
	{
	if(pEl1->StartLoci < pEl2->StartLoci)
		return(-1);
	if(pEl1->StartLoci > pEl2->StartLoci)
		return(1);
	if(pEl1->EndLoci < pEl2->EndLoci) 
		return(-1);									
	if(pEl1->EndLoci > pEl2->EndLoci)
		return(1);
	return(0);
	}
return(pEl1->ChromID < pEl2->ChromID ? -1 : 1);
}


