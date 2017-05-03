// genalignstats.cpp : Defines the entry point for the console application.
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

const unsigned int cProgVer = 302;		// increment with each release



const int cMAFmaxSpecies = 50;			// maximum number of species that can be handled
const int cMAtotMaxSeqAlignLen = 0x03ffffff; // total (over all aligned species) max seq length which can be buffered

const int cMaxIncludeFiles = 10;		// maximun number of include region filter files
const int cMaxExcludeFiles = 10;		// maximun number of exclude region filter files

const int cMaxIncludeChroms = 20;		// max number of include chromosomes regular expressions
const int cMaxExcludeChroms = 20;		// max number of exclude chromosomes regular expressions

typedef enum eProcMode {
	eProcModeStandard = 0,				// default processing - report aligned columns (can contain InDels) vs aligned and identical 
	eProcModeNoInDels,					// report aligned columns (InDels not allowed) vs aligned and identical
	eProcModePairwise					// report pairwise substitutions (relative to reference species)
} etProcMode;

typedef struct TAG_sExcludeSpeciesChrom {
	int SpeciesID;						// which species
	int ChromID;						// chromosome is not to be processed
	} tsExcludeSpeciesChrom;

typedef struct TAG_sBaseStats {
	int MonoCnts[4][4];					// counts for ref base to rel base substitutions as [RefBase][RelBase]
} tsBaseStats;

typedef struct TAG_sRegionStats {
	tsBaseStats Regions[9];
} tsRegionStats;


typedef struct TAG_sColStats {
	int Aligned[7];					// counts for number of aligned columns (a,c,g,t,N,-) where at least 1 species is <= eBaseN
	int Identical[7];				// counts for number of aligned columns in which all species are identical to refspecies
} tsColStats;

typedef struct TAG_sColRegionStats {
	tsColStats Regions[9];
} tsColRegionStats;

typedef struct TAG_sProcParams 
	{
	int ProcMode;							// processing mode 0: default
	bool bStatsAvail;

	int NumSpecies;						// number of species in szSpecies[]
	char *pszSpeciesList;				// comma separated species list starting with the reference species
	int RefSpeciesIdx;					// current index into pSeqs for the reference species
	char szSpecies[cMAFmaxSpecies][cMaxDatasetSpeciesChrom];	// species names of interest - other species are sloughed
										// only alignments with species included will be processed
										// first species is the reference species
	int NumSpeciesAligned;				// actual number of sequences in current alignment
	int MinNumSpeciesAligned;			// minimum number of species (from species list) required in alignment block
	char szRefChrom[cMaxDatasetSpeciesChrom];
	int RefChromID;						// current reference chromosome identifier
	int SeqOfs;							// offset into sequences
	int SeqLen;							// sequence length
	etSeqBase *pSeqs[cMAFmaxSpecies];	// ptrs to each sequence
	int MaxSeqAlignLen;					// max length alignment which can be processed (how much mem was alloc'd to pSeq[n])			
	int Regions;						// number of regions per step

	tsColRegionStats ColStats;			// stats for number of aligned columns (to ref species) and cols which are identical over all species

	tsRegionStats RegionStats[cMAFmaxSpecies];	// ref to rel species counts for each region

	CBEDfile *pBiobed;					// if not NULL then opened biobed file for regional characteristics
	int BEDChromID;						// BED chromosome identifier corresponding to RefChromID
	int NumIncludes;					// number of biobed files containing regions to include
	int NumExcludes;					// number of biobed files containing regions to exclude
	CBEDfile *pIncludes[cMaxIncludeFiles];	// if opened biobed files for regions to include - all other regions are to be excluded
	CBEDfile *pExcludes[cMaxExcludeFiles];	// if opened biobed files for regions to exclude 
	int UpDnStreamLen;					// up/dn stream regional length when characterising
	bool bMultipleFeatBits;				// if false then stats only generated if a single feature bit is set - e.g if both exons and introns overlapped then no stat generated
	int hRsltsFile;						// write stats results into this CSV file
	int NumIncludeChroms;				// number of chromosomes explicitly defined to be included
	char **ppszIncludeChroms;			// ptr to array of reg expressions defining chroms to include - overides exclude
	int NumExcludeChroms;				// number of chromosomes explicitly defined to be excluded
	char **ppszExcludeChroms;			// ptr to array of reg expressions defining chroms to include
#ifdef _WIN32
	Regexp *IncludeChromsRE[cMaxIncludeChroms];	// compiled regular expressions
	Regexp *ExcludeChromsRE[cMaxExcludeChroms];
#else
	regex_t IncludeChromsRE[cMaxIncludeChroms];	// compiled regular expressions
	regex_t ExcludeChromsRE[cMaxExcludeChroms];
#endif
	} tsProcParams; 

int
Process(int ProcMode,				// processing mode 0
				 char *pszInputFile,		// bio multialignment (.algn) file to process
					char *pszOutputFile,	// where to write out stats
					char *pszBiobedFile,	// biobed file containing regional features - exons, introns etc
					int NumIncludeFiles,	// number of include region files
					char **ppszIncludeFiles,	// biobed files containing regions to include - default is to exclude none
					int NumExcludeFiles,	// number of exclude region files
					char **ppszExcludeFiles,	// biobed file containing regions to exclude - default is to include all
					int MinNumAlignSpecies,	// minimum number of species (from species list) which must be in any aligned block before that block is processed
 					char *pszSpeciesList,	// space or comma separated list of species
					int RegLen,				// regulatory region length - up/dn stream of 5/3' 
					bool bMultipleFeatBits,	// if true then accept alignments in which multiple feature bits are set
					int NumIncludeChroms,	// number of chromosomes explicitly defined to be included
					char **ppszIncludeChroms,	// ptr to array of reg expressions defining chroms to include - overides exclude
					int NumExcludeChroms,	// number of chromosomes explicitly defined to be excluded
					char **ppszExcludeChroms);	// ptr to array of reg expressions defining chroms to include

int ParseNumSpecies(char *pszSpeciesList,tsProcParams *pProcParams);
int TrimQuotes(char *pszToTrim);
bool ProcAlignBlockSummary(int RefChromID,int RefChromOfs,int AlignLen,tsProcParams *pProcParams);
bool ProcessAlignment(int RefChromID,int ChromOfs,int SeqIdx,int SubSeqLen,tsProcParams *pProcParams);
bool IncludeFilter(int SubRefOfs,int SubRefEndOfs,tsProcParams *pProcParams);

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

int ProcAlignBlock(int RefChromID,int RefChromOfs,int RefAlignLen,tsProcParams *pProcParams);
int ReportResults(tsProcParams *pProcParams);
char *MapRegion2Txt(int RegionIdx);

CStopWatch gStopWatch;
CDiagnostics gDiagnostics;				// for writing diagnostics messages to log file
char gszProcName[_MAX_FNAME];			// process name

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
	return _T("genalignstats");
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
int Rslt = 0;

int iProcMode;
bool bMultipleFeatBits;
int Idx;
int LenReq;
char szOutputFile[_MAX_PATH];
char szInputFile[_MAX_PATH];
char szInputBiobedFile[_MAX_PATH];

int NumIncludeChroms;
char *pszIncludeChroms[cMaxIncludeChroms];
int NumExcludeChroms;
char *pszExcludeChroms[cMaxExcludeChroms];
int NumIncludeFiles;
char *pszIncludeFiles[cMaxExcludeFiles];
int NumExcludeFiles;
char *pszExcludeFiles[cMaxIncludeFiles];
int iMinNumSpecies;
char szSpeciesList[512];
int iNumSpecies;
int iRegLen;

// command line args
struct arg_lit  *help    = arg_lit0("hH","help",                "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_int *ScreenLogLevel=arg_int0("S", "ScreenLogLevel",	"<int>","Level of diagnostics written to logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"diagnostics log file");

struct arg_int  *ProcMode = arg_int0("m","procmode","<int>",	"processing mode 0: aligned vs identical columns (InDels allowed), 1: aligned vs identical columns (no InDels), 2: pairwise (relative to ref species) substitutions");
struct arg_int  *RegLen = arg_int0("L","updnstream","<int>",	"length of 5'up or 3'down  stream regulatory region length (default = 2000) 0..1000000");
struct arg_file *InFile = arg_file1("i",NULL,"<file>",			"input from .algn files");
struct arg_file *OutFile = arg_file1("o",NULL,"<file>",			"output to statistics file as CSV");
struct arg_file *InBedFile = arg_file0("b","bed","<file>",		"characterise regions from biobed file");
struct arg_int  *MinNumSpecies = arg_int0("M","minspecies","<int>",	"minimum number of species (from species list) in any alignment block (default 2)");
struct arg_str  *SpeciesList = arg_str1("s","species","<string>","comma or space delimited species list, first species must be the reference species");
struct arg_lit  *MultipleFeatBits  = arg_lit0("q","multiplefeatbits",	"single featbit (default) or multiple featbits allowed");
struct arg_file *ExcludeFile = arg_filen("E","exclude","<file>",0,cMaxExcludeFiles,	"exclude all regions in biobed file from processing ");
struct arg_file *IncludeFile = arg_filen("I","include","<file>",0,cMaxExcludeFiles,	"include all regions (unless specific regions excluded) in biobed file");
struct arg_str  *IncludeChroms = arg_strn("z","chromeinclude","<string>",0,cMaxIncludeChroms,"low priority - regular expressions defining species.chromosomes to include for processing");
struct arg_str  *ExcludeChroms = arg_strn("Z","chromexclude","<string>",0,cMaxExcludeChroms,"high priority - regular expressions defining species.chromosomes to exclude from processing");

struct arg_end *end = arg_end(20);

void *argtable[] = {help,version,FileLogLevel,ScreenLogLevel,LogFile,
					ProcMode,RegLen,
					InFile,OutFile,InBedFile,
					MinNumSpecies,SpeciesList,MultipleFeatBits,
					ExcludeFile,IncludeFile,IncludeChroms,ExcludeChroms,
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

	iProcMode = ProcMode->count ? ProcMode->ival[0] : 0;
	if(iProcMode < eProcModeStandard || iProcMode > eProcModePairwise)
		{
		printf("\nError: Requested processing mode '-m%d' not supported",iProcMode);
		exit(1);
		}


	bMultipleFeatBits = MultipleFeatBits->count ? true : false;

	strncpy(szInputFile,InFile->filename[0],_MAX_PATH);
	szInputFile[_MAX_PATH-1] = '\0';
	strncpy(szOutputFile,OutFile->filename[0],_MAX_PATH);
	szOutputFile[_MAX_PATH-1] = '\0';

	if(InBedFile->count)
		{
		strncpy(szInputBiobedFile,InBedFile->filename[0],_MAX_PATH);
		szInputBiobedFile[_MAX_PATH-1] = '\0';
		}
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
	
	if(!SpeciesList->count)
		{
		printf("\nError: Species list '-s<specieslist>' is empty\n");
		exit(1);
		}

	strcpy(szSpeciesList,SpeciesList->sval[0]);
	TrimQuotes(szSpeciesList);
	iNumSpecies = ParseNumSpecies(szSpeciesList,NULL);
	if(iNumSpecies < 2)
		{
		printf("\nAt least 2 species must be specified in species list '-s%s'",szSpeciesList);
		exit(1);
		}

	iMinNumSpecies = MinNumSpecies->count ? MinNumSpecies->ival[0] : 2;
	if(iMinNumSpecies < 2 || iMinNumSpecies > iNumSpecies)
		{
		printf("\nMinimum number of species in alignment block '-M%d' must be in range 2..%d",iMinNumSpecies,iNumSpecies);
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
	switch(iProcMode) {
		case eProcModeStandard:
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Processing mode: Aligned columns vs identical, can contain InDels");
			break;
		case eProcModeNoInDels:
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Processing mode: Aligned columns vs identical, no InDels allowed");
			break;
		case eProcModePairwise:
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Processing mode: pairwise substitutions (relative to reference species)");
			break;
		}
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"bio multialignment (.algn) file to process: '%s'",szInputFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"where to write out stats: '%s'",szOutputFile);
	if(szInputBiobedFile[0]!='\0')
		{
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"biobed file containing regional features: '%s'",szInputBiobedFile);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"accept alignments in which multiple feature bits are set: %s",bMultipleFeatBits ? "yes" : "no");
		}

	for(Idx = 0; Idx < NumIncludeFiles; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"biobed file containing regions to include: '%s'",pszIncludeFiles[Idx]); 
	for(Idx = 0; Idx < NumExcludeFiles; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"biobed file containing regions to exclude: '%s'",pszExcludeFiles[Idx]);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"minimum number of species in alignment block: %d",iMinNumSpecies);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"species list: '%s'",	szSpeciesList);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"regulatory region length: %d",iRegLen);

	for(Idx = 0; Idx < NumIncludeChroms; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"reg expressions defining chroms to include: '%s'",pszIncludeChroms[Idx]);
	for(Idx = 0; Idx < NumExcludeChroms; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"reg expressions defining chroms to exclude: '%s'",pszExcludeChroms[Idx]); 

	gStopWatch.Start();
#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif

	Rslt = Process(iProcMode,			// processing mode 0
				 szInputFile,		// bio multialignment (.algn) file to process
				 szOutputFile,		// where to write out stats
				 szInputBiobedFile,	// biobed file containing regional features - exons, introns etc
				 NumIncludeFiles,	// number of include region files
				 pszIncludeFiles,	// biobed files containing regions to include - default is to exclude none
				 NumExcludeFiles,	// number of exclude region files
				 pszExcludeFiles,	// biobed file containing regions to exclude - default is to include all
				 iMinNumSpecies,	// minimum number of species (from species list) required in alignment block
 				 szSpeciesList,		// space or comma separated list of species, priority ordered
				 iRegLen,			// regulatory region length - up/dn stream of 5/3' 
				 bMultipleFeatBits,	// if true then accept alignments in which multiple feature bits are set
				 NumIncludeChroms,	// number of chromosomes explicitly defined to be included
				 pszIncludeChroms,	// ptr to array of reg expressions defining chroms to include - overides exclude
				 NumExcludeChroms,	// number of chromosomes explicitly defined to be excluded
				 pszExcludeChroms);	// ptr to array of reg expressions defining chroms to include

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
		if(NumSpecies >= cMAFmaxSpecies)
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

return(AddExcludeHistory(SpeciesID,ChromID,pProcParams->NumIncludeChroms > 0 ? true : false));
}

int 
ProcessAlignments(char *pszMAF,			 // source multialignment file
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
int SpeciesIDs[cMAFmaxSpecies];
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
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open MAF file '%s'",pszMAF);
	return(Rslt);
	}

// ensure all requested species are represented in multispecies alignment file plus get their species identifiers
for(Idx = 0; Idx < pProcParams->NumSpecies; Idx++)
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


// iterate over reference blocks which are sorted by chrom then offset
CurBlockID = 0;
PrevRefChromID = 0;
PrevDispRefChromID = 0;
BEDChromID = 0;
pProcParams->RefSpeciesIdx = 0;		// reference sequence will always be 1st
memset(&pProcParams->RegionStats,0,sizeof(pProcParams->RegionStats));

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
		PrevDispRefChromID = RefChromID;
		pszRefChrom = pAlignments->GetChromName(RefChromID);

		strcpy(pProcParams->szRefChrom,pszRefChrom);
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing chromosome %s",pszRefChrom);
		if(pProcParams->pBiobed != NULL)
			pProcParams->BEDChromID = pProcParams->pBiobed->LocateChromIDbyName(pProcParams->szRefChrom);
		else
			pProcParams->BEDChromID = 0;
		pProcParams->RefChromID = RefChromID;
		}

	if(!bLoaded)
		continue;

	if(RefChromID != PrevRefChromID)
		PrevRefChromID = RefChromID;

	if(ProcAlignBlock(RefChromID,RefChromOfs,RefAlignLen,pProcParams) < eBSFSuccess)
		break;
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
int NumAlignedSpecies;
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

	// should this reference chromosome be accepted for processing or sloughed?
	if(ExcludeThisChrom(pAlignments,RefSpeciesID,CurRefChromID,pProcParams))	
		break;

			// terminate with blocks in which there is no alignment to at least a minimum number of species
	if((CurNumAlignments = pAlignments->GetNumSpecies(CurBlockID)) < pProcParams->MinNumSpeciesAligned)
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

		// iterate over all species
		// if species not in current block then set it's sequence to unaligned
	for(NumAlignedSpecies = Idx = 0; Idx < pProcParams->NumSpecies; Idx++) 
		{
		RelSpeciesID = pSpeciesIDs[Idx];
		RelChromID = pAlignments->GetRelChromID(CurBlockID,RelSpeciesID);
		if(RelChromID < 1)			// assume < 1 because species not in alignment
			{
			memset(&pProcParams->pSeqs[Idx][RefAlignLen],eBaseUndef,CurRefAlignLen);
			continue;		
			}

		if(ExcludeThisChrom(pAlignments,RelSpeciesID,RelChromID,pProcParams))	// should this chromosome be accepted for processing or sloughed?
			{
			memset(&pProcParams->pSeqs[Idx][RefAlignLen],eBaseUndef,CurRefAlignLen);
			continue;		
			}

			// get sequence
		if((pSeq = pAlignments->GetSeq(CurBlockID,RelSpeciesID))==NULL) // some species are not present in all blocks
			{
			memset(&pProcParams->pSeqs[Idx][RefAlignLen],eBaseUndef,CurRefAlignLen);
			continue;		
			}

		memcpy(&pProcParams->pSeqs[Idx][RefAlignLen],pSeq,CurRefAlignLen);
		NumAlignedSpecies++;
		}
		
	// check if required minimum number of species required were in multialignment
	if(NumAlignedSpecies < pProcParams->MinNumSpeciesAligned)	
		break;

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
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"\nUnable to instantiate CBEDfile '%s'",pToOpen);
		return(NULL);
		}

	if((Rslt = pBed->Open(pToOpen,eBTAnyBed))!=eBSFSuccess)
		{
		while(pBed->NumErrMsgs())
			gDiagnostics.DiagOut(eDLFatal,gszProcName,pBed->GetErrMsg());
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open biobed file '%s'",pToOpen);
		delete pBed;
		return(NULL);
		}
	return(pBed);
	}
return(NULL);
}

//Cleanup
//Closes and deletes all resources 
bool
Cleanup(tsProcParams *pProcParams)
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

if(pProcParams->hRsltsFile != -1)
	{
	close(pProcParams->hRsltsFile);
	pProcParams->hRsltsFile = -1;
	}

for(Idx = 0; Idx < pProcParams->NumSpecies; Idx++)
	{
	if(pProcParams->pSeqs[Idx] != NULL)
		delete pProcParams->pSeqs[Idx];
	pProcParams->pSeqs[Idx] = NULL;
	}
pProcParams->NumSpecies = 0;

#ifdef _WIN32
for(Idx=0;Idx < pProcParams->NumIncludeChroms;Idx++)
	{
	if(pProcParams->IncludeChromsRE[Idx] != NULL)
		delete pProcParams->IncludeChromsRE[Idx];
	pProcParams->IncludeChromsRE[Idx] = NULL;
	}
pProcParams->NumIncludeChroms = 0;

for(Idx=0;Idx < pProcParams->NumExcludeChroms;Idx++)
	{
	if(pProcParams->ExcludeChromsRE[Idx] != NULL)
		delete pProcParams->ExcludeChromsRE[Idx];
	pProcParams->ExcludeChromsRE[Idx] = NULL;
	}
pProcParams->NumExcludeChroms = 0;
#endif

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
IncludeFilter(int SubRefOfs,int SubRefEndOfs,tsProcParams *pProcParams)
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
// Generate stats from alignment file
int
Process(int ProcMode,							// processing mode 0
				 char *pszInputFile,		// bio multialignment (.algn) file to process
					char *pszOutputFile,	// where to write out stats
					char *pszBiobedFile,	// biobed file containing regional features - exons, introns etc
					int NumIncludeFiles,	// number of include region files
					char **ppszIncludeFiles,	// biobed files containing regions to include - default is to exclude none
					int NumExcludeFiles,	// number of exclude region files
					char **ppszExcludeFiles,	// biobed file containing regions to exclude - default is to include all
					int MinNumAlignSpecies,	// minimum number of species which must be in any aligned block before that block is processed
					char *pszSpeciesList,	// space or comma separated list of species
					int RegLen,				// regulatory region length - up/dn stream of 5/3' 
					bool bMultipleFeatBits,	// if true then accept alignments in which multiple feature bits are set
					int NumIncludeChroms,	// number of chromosomes explicitly defined to be included
					char **ppszIncludeChroms,	// ptr to array of reg expressions defining chroms to include - overides exclude
					int NumExcludeChroms,	// number of chromosomes explicitly defined to be excluded
					char **ppszExcludeChroms)	// ptr to array of reg expressions defining chroms to include
{
int Rslt;
int NumSpecies;
int Idx;
int Regions;
int UpDnStreamLen;
CBEDfile *pBiobed = NULL;
tsProcParams ProcParams;
//int	MismatchScore;
//int	MatchScore;
char szCSVSpecies[512];				// to hold comma separated species list

#ifndef _WIN32
int RegErr;				// regular expression parsing error
char szRegErr[128];			// to hold RegErr as textual representation ==> regerror();
#endif

memset(&ProcParams,0,sizeof(tsProcParams));

// parse out species list
NumSpecies = ParseNumSpecies(pszSpeciesList,&ProcParams);
szCSVSpecies[0]='\0';
for(Idx = 0; Idx < NumSpecies; Idx++)
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
		Cleanup(&ProcParams);
		return(eBSFerrObj);
		}
	ProcParams.NumIncludes++;
	}

for(Idx=0;Idx<NumExcludeFiles; Idx++)
	{
	if((ProcParams.pExcludes[Idx] = OpenBedfile(ppszExcludeFiles[Idx]))==NULL)
		{
		Cleanup(&ProcParams);
		return(eBSFerrObj);
		}
	ProcParams.NumExcludes++;
	}


if(pszBiobedFile != NULL && pszBiobedFile[0] != '\0')
	{
	if((ProcParams.pBiobed = OpenBedfile(pszBiobedFile))==NULL)
		{
		Cleanup(&ProcParams);
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

ProcParams.ProcMode = ProcMode;
ProcParams.NumSpecies = NumSpecies;
ProcParams.MinNumSpeciesAligned = MinNumAlignSpecies;			
ProcParams.RefSpeciesIdx = 0;

ProcParams.UpDnStreamLen = UpDnStreamLen;
ProcParams.bMultipleFeatBits = bMultipleFeatBits;
ProcParams.Regions = Regions;

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
	Cleanup(&ProcParams);
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
	Cleanup(&ProcParams);
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
		Cleanup(&ProcParams);
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
		Cleanup(&ProcParams);
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to process exclude chrom '%s' error: %s",ppszExcludeChroms[Idx],szRegErr);
		return(eBSFerrMem);
		}


	}
#endif

// determine max aligned sequence length for any single species which can be handled
ProcParams.MaxSeqAlignLen = cMAtotMaxSeqAlignLen/NumSpecies;
for(Idx = 0; Idx < NumSpecies; Idx++)
	{
	if((ProcParams.pSeqs[Idx] = new unsigned char [ProcParams.MaxSeqAlignLen])==NULL)
		{
	    Cleanup(&ProcParams);
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
    Cleanup(&ProcParams);
	return(eBSFerrCreateFile);
	}

Rslt = ProcessAlignments(pszInputFile,&ProcParams);
if(Rslt == eBSFSuccess)
	ReportResults(&ProcParams);

if(ProcParams.hRsltsFile != -1)
	{
	close(ProcParams.hRsltsFile);
	ProcParams.hRsltsFile = -1;
	}

Cleanup(&ProcParams);
return(Rslt);
}

int
ReportResults(tsProcParams *pProcParams)
{
int LineID;
int SpeciesIdx;
int RegionIdx;
int RefBase;
int RelBase;
int TotSpeciesBases;
int LineLen;
tsBaseStats *pStats;
tsColStats *pColStats;
tsRegionStats *pRegionStats;
char szLineBuff[0x03fff];

LineID = 0;
LineLen = 0;

if(pProcParams->ProcMode > eProcModePairwise)
	{
	for(SpeciesIdx = 1; SpeciesIdx < pProcParams->NumSpecies; SpeciesIdx++)
		{
		pRegionStats = &pProcParams->RegionStats[SpeciesIdx]; 
		pStats = &pRegionStats->Regions[0];

		for(RegionIdx = 0; RegionIdx < pProcParams->Regions; RegionIdx++,pStats++)
			{
			for(RefBase = 0; RefBase < 4; RefBase++)
				{
				for(TotSpeciesBases = RelBase = 0; RelBase < 4; RelBase++)
					TotSpeciesBases += pStats->MonoCnts[RefBase][RelBase];

				for(RelBase = 0; RelBase < 4; RelBase++)
					{
					LineLen += sprintf(&szLineBuff[LineLen],"%d,\"%s\",\"%s\",\"%s\",\"%c->%c\",%d,%1.8f\n",
							++LineID,pProcParams->szSpecies[0],pProcParams->szSpecies[SpeciesIdx],
							pProcParams->Regions == 1 ? "All" : MapRegion2Txt(RegionIdx),
							CSeqTrans::MapBase2Ascii(RefBase), CSeqTrans::MapBase2Ascii(RelBase), pStats->MonoCnts[RefBase][RelBase],
							TotSpeciesBases == 0 ? 0.0 : (1.0 * pStats->MonoCnts[RefBase][RelBase])/TotSpeciesBases);
					if(LineLen > ((sizeof(szLineBuff) * 9) / 10))
						{
						CUtility::SafeWrite(pProcParams->hRsltsFile,szLineBuff,LineLen);
						LineLen = 0;
						}
					}
				}
			}
		}
	}
else
	{
	// report aligned vs identical columns
	pColStats = &pProcParams->ColStats.Regions[0];
	for(RegionIdx = 0; RegionIdx < pProcParams->Regions; RegionIdx++,pColStats++)
		{
		for(RefBase = 0; RefBase <= eBaseInDel; RefBase++)
			{
			if(RefBase == eBaseUndef)
				continue;
			LineLen += sprintf(&szLineBuff[LineLen],"%d,\"%s\",\"%c\",%d,%d\n",
					++LineID,
					pProcParams->Regions == 1 ? "All" : MapRegion2Txt(RegionIdx),
					CSeqTrans::MapBase2Ascii(RefBase),pColStats->Aligned[RefBase],pColStats->Identical[RefBase]);
			if(LineLen > ((sizeof(szLineBuff) * 9) / 10))
				{
				CUtility::SafeWrite(pProcParams->hRsltsFile,szLineBuff,LineLen);
				LineLen = 0;
				}
			}
		}
	}

if(LineLen)
	CUtility::SafeWrite(pProcParams->hRsltsFile,szLineBuff,LineLen);
return(eBSFSuccess);
}

// MapLoci2RegionIdx
// Maps chromosome loci to a region index in 5'->3' logical order e.g. Intergenic = 0,...3' downstream
// Defaults to intergenic if unable to map loci to bedfile features
// Returns -1 if loci overlaps more than 1 feature and multiple feature overlaps not allowed 
const int cRegionFeatIdxMsk = 0x07;		// mask to retain intergenic...3'ds region index, removing splice site flags
// splice site flags are cIntronExonSpliceSite and cExonIntronSpliceSite

int
MapLoci2RegionIdx(int Loci,tsProcParams *pProcParams)
{
int RegionIdx;
int FeatureBits;
int SpliceSiteOverlaps;
int BitMsk;
int FeatIdx;

RegionIdx = 0;
FeatureBits = 0;
SpliceSiteOverlaps = 0;
if(pProcParams->pBiobed != NULL)
	{
	if(pProcParams->BEDChromID > 0)
		FeatureBits = pProcParams->pBiobed->GetFeatureBits(pProcParams->BEDChromID,Loci,Loci,cRegionFeatBits,pProcParams->UpDnStreamLen);
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
					return(-1);	// although was base of interest, more than one feature bit so can't contribute to stats
				RegionIdx = FeatIdx;
				if(pProcParams->bMultipleFeatBits)	// if multiple features allowed then don't check for any additional
					break;
				}
			}
		}

		// need to remap RegionIdx to 5'->3' logical order
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
	
	if(RegionIdx != 0)		// if not intergenic then check for splice sites
		{
		SpliceSiteOverlaps = pProcParams->pBiobed->GetSpliceSiteBits(pProcParams->BEDChromID,Loci,Loci,cMinSpliceOverlap);
		if(SpliceSiteOverlaps & cIntronExonSpliceSite)
			RegionIdx |= cIntronExonSpliceSite;
		if(SpliceSiteOverlaps & cExonIntronSpliceSite)
			RegionIdx |= cExonIntronSpliceSite;
		}
	}
return(RegionIdx);
}

char *
MapRegion2Txt(int RegionIdx)
{
switch(RegionIdx) {
	case 0:
		return((char *)"Intergenic");
	case 1:
		return((char *)"5'Upstream");
	case 2:
		return((char *)"5'UTR");
	case 3:
		return((char *)"CDS");
	case 4:
		return((char *)"Intron");
	case 5:
		return((char *)"3'UTR");
	case 6:
		return((char *)"3'Dnstream");
	case 7:
		return((char *)"5'Splice");
	case 8:
		return((char *)"3'Splice");
	default:
		break;
	}
return((char *)"Unrecognised region");
}

int
ProcAlignBlock(int RefChromID,int RefChromOfs,int RefAlignLen,tsProcParams *pProcParams)
{
int Idx;
int SpeciesIdx;
int CurRefChromOfs;
int RegionIdx;
tsBaseStats *pStats;
tsColStats *pColStats;
tsRegionStats *pRegionStats;
etSeqBase *pRefBase;
etSeqBase *pRelBase;
etSeqBase RefBase;
etSeqBase RelBase;
int AlignedBases;
int AlignedExactBases;
int NumInDelBases;

if(pProcParams->ProcMode > eProcModePairwise)
	{
	for(SpeciesIdx = 1; SpeciesIdx < pProcParams->NumSpecies; SpeciesIdx++)
		{
		pRegionStats = &pProcParams->RegionStats[SpeciesIdx];
		pRefBase = pProcParams->pSeqs[0];
		pRelBase = pProcParams->pSeqs[SpeciesIdx];
		CurRefChromOfs = RefChromOfs;
		for(Idx = 0; Idx < RefAlignLen; Idx++,pRefBase++,pRelBase++)
			{
			RefBase = (*pRefBase & ~cRptMskFlg);
			if(RefBase <= eBaseT && ((RelBase = (*pRelBase & ~cRptMskFlg)) <= eBaseT) &&
				IncludeFilter(CurRefChromOfs,CurRefChromOfs,pProcParams))
				{
				RegionIdx = MapLoci2RegionIdx(CurRefChromOfs,pProcParams);
				if(RegionIdx >= 0)
					{
					pStats = &pRegionStats->Regions[RegionIdx & cRegionFeatIdxMsk];
					pStats->MonoCnts[RefBase][RelBase] += 1;
					if(RegionIdx & cIntronExonSpliceSite)
						{
						pStats = &pRegionStats->Regions[7];
						pStats->MonoCnts[RefBase][RelBase] += 1;
						}
					if(RegionIdx & cExonIntronSpliceSite)
						{
						pStats = &pRegionStats->Regions[8];
						pStats->MonoCnts[RefBase][RelBase] += 1;
						}
					}
				}
			if(RefBase <= eBaseN)
				CurRefChromOfs += 1;
			}
		}
	}
else
	{
	CurRefChromOfs = RefChromOfs;
	pRefBase = pProcParams->pSeqs[0];
	for(Idx = 0; Idx < RefAlignLen; Idx++,pRefBase++)
		{
		RefBase = *pRefBase & ~cRptMskFlg;
		if(RefBase <= eBaseN)
			CurRefChromOfs += 1;

		if((RefBase > eBaseN) && (pProcParams->ProcMode == eProcModeNoInDels || RefBase != eBaseInDel))
			continue;

		if(!IncludeFilter(CurRefChromOfs-1,CurRefChromOfs-1,pProcParams))
			continue;

		RegionIdx = MapLoci2RegionIdx(CurRefChromOfs-1,pProcParams);
		if(RegionIdx < 0)
			continue;
		
		AlignedBases = 1;	
		if(RefBase <= eBaseT)
			AlignedExactBases = 1;
		else
			AlignedExactBases = 0;
		if(RefBase == eBaseInDel)
			NumInDelBases = 1;
		else
			NumInDelBases = 0;

		for(SpeciesIdx = 1; SpeciesIdx < pProcParams->NumSpecies; SpeciesIdx++)
			{
			pRelBase = pProcParams->pSeqs[SpeciesIdx];
			pRelBase += Idx;
			RelBase = (*pRelBase & ~cRptMskFlg);
			if(RelBase <= eBaseN || (pProcParams->ProcMode == eProcModeStandard && RelBase == eBaseInDel))
				{
				AlignedBases += 1;
				if(RefBase <= eBaseT && RelBase == RefBase)				
					AlignedExactBases += 1;
				if(RelBase == eBaseInDel)
					NumInDelBases += 1;
				}
			}

		if(AlignedBases >= pProcParams->MinNumSpeciesAligned && AlignedBases > NumInDelBases) 
			{
			pColStats = &pProcParams->ColStats.Regions[RegionIdx & cRegionFeatIdxMsk];
			pColStats->Aligned[RefBase] += 1;
			if(AlignedExactBases == AlignedBases)
				pColStats->Identical[RefBase] += 1;

			if(RegionIdx & cIntronExonSpliceSite)
				{
				pColStats = &pProcParams->ColStats.Regions[7];
				pColStats->Aligned[RefBase] += 1;
				if(AlignedExactBases == AlignedBases)
					pColStats->Identical[RefBase] += 1;
				}
			if(RegionIdx & cExonIntronSpliceSite)
				{
				pColStats = &pProcParams->ColStats.Regions[8];
				pColStats->Aligned[RefBase] += 1;
				if(AlignedExactBases == AlignedBases)
					pColStats->Identical[RefBase] += 1;
				}
			}
		}
	}

return(eBSFSuccess);
}

