// genultras.cpp : Defines the entry point for the console application.
// Locates all ultraconserved elements according to parameterised filtering critera
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

const unsigned int cProgVer = 305;		// increment with each release


const int cMAtotMaxSeqAlignLen = 0x0fffffff; // total (over all aligned species) max seq length that can be buffered in concatenated seqs

const int cMinUltraLen = 4;				// allow core lengths to be specified down to cMinCoreLen
const int cDfltUltraLen= 50;			// if core lengths not specified then default to cDfltMinCoreLen
const int cMaxUltraLen = 10000;		// minimum core lengths can be specified upto this length

const int cMaxIncludeFiles = 10;		// maximun number of include region filter files
const int cMaxExcludeFiles = 10;		// maximun number of exclude region filter files

const int cMaxIncludeChroms = 20;		// max number of include chromosomes regular expressions
const int cMaxExcludeChroms = 20;		// max number of exclude chromosomes regular expressions

typedef enum eProcMode {
	eProcModeStandard = 0				// default processing
} etProcMode;

CStopWatch gStopWatch;
CDiagnostics gDiagnostics;				// for writing diagnostics messages to log file
char gszProcName[_MAX_FNAME];			// process name


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
	int RefSpeciesIdx;				// current index into pSeqs for the reference species
	char szSpecies[cMaxAlignedSpecies][cMaxDatasetSpeciesChrom];	// species names of interest - other species are sloughed
									// only alignments with species included will be processed
									// first species is the reference species
	int MaxNumSpecies;				// max number of species to be processed (number of species in species list)
	int NumCoreSpecies;				// number of core species required in an alignment (1st NumCoreSpecies in species list must be present in alignment)
	int MinNumSpecies;				// minimum number of species including core species to be in alignment (will always be at least NumCoreSpecies)
	int NumSpecies;					// actual number of sequences in current alignment
	int MaxAlignIdxSpecies;			// pSeqs[MaxAlignIdxSpecies-1] is last actual alignment sequence
	char szRefChrom[cMaxDatasetSpeciesChrom]; // current reference chromosome
	int RefChromID;					// current reference chromosome identifier
	int	WindowSize;					// sampling window size
	int NxtOutputOffset;			// when to next output results - i.e end of current window
	int SeqOfs;						// offset into sequences
	int SeqLen;						// sequence length
	etSeqBase *pSeqs[cMaxAlignedSpecies];  // ptrs to each sequence
	int MaxSeqAlignLen;				// max length alignment which can be processed (how much mem was alloc'd to pSeq[n])			
	int MinUltraLen;				// minimum ultra core length required
	CBEDfile *pBiobed;				// if not NULL then opened biobed file for regional characteristics
	int BEDChromID;					// BED chromosome identifier corresponding to RefChromID
	int NumIncludes;				// number of biobed files containing regions to include
	int NumExcludes;				// number of biobed files containing regions to exclude
	CBEDfile *pIncludes[cMaxIncludeFiles];	// if opened biobed files for regions to include - all other regions are to be excluded
	CBEDfile *pExcludes[cMaxExcludeFiles];	// if opened biobed files for regions to exclude 
	int UpDnStreamLen;				// up/dn stream regional length when characterising
	int hRsltsFile;					// write ultra loci into this CSV file
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

	} tsProcParams; 

int ParseNumSpecies(char *pszSpeciesList,tsProcParams *pProcParams);
int ParseUniqueSpeciesSeqs(char *pszUniqueSpeciesSeqs,tsProcParams *pProcParams);

bool ProcAlignBlock(int RefChromID,int RefChromOfs,int AlignLen,tsProcParams *pProcParams);
bool ProcessAlignment(int RefChromID,int ChromOfs,int SeqIdx,int SubSeqLen,tsProcParams *pProcParams);
int ProcessSubSeq(int RefChromID,int ChromOfs,int SeqIdx,int MaxLen,tsProcParams *pProcParams);
bool OutputUltraCore(char *pszChrom, int ChromStartOffset, int ChromEndOffset, int FeatureBits,tsProcParams *pProcParams);
char *ChkSpeciesChromWellFormed(char *pszSpeciesChroms);
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
Process(int ProcMode,				// processing mode 0: default
				 char *pszInputFile,		// bio multialignment (.algn) file to process
					char *pszOutputFile,	// where to write out loci
					char *pszBiobedFile,	// biobed file containing regional features - exons, introns etc
					int NumIncludeFiles,	// number of include region files
					char **ppszIncludeFiles,// biobed files containing regions to include - default is to exclude none
					int NumExcludeFiles,	// number of exclude region files
					char **ppszExcludeFiles,// biobed file containing regions to exclude - default is to include all
					int NumCoreSpecies,		// number of core species required to be in alignment
					int MinNumSpecies,		// minimum number of species including NumCoreSpecies to be in alignment
					char *pszSpeciesList,	// space or comma separated list of species, priority ordered
					int MinUltraLen,		// minimum ultra core length required
					int RegLen,				// regulatory region length - up/dn stream of 5/3' 
					int NumIncludeChroms,	// number of chromosomes explicitly defined to be included
					char **ppszIncludeChroms,	// ptr to array of reg expressions defining chroms to include - overides exclude
					int NumExcludeChroms,	// number of chromosomes explicitly defined to be excluded
					char **ppszExcludeChroms);	// ptr to array of reg expressions defining chroms to include


int NormaliseInDelColumns(tsProcParams *pProcParams,int AlignLen);


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
	return _T("genultras");
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

int iProcMode;
int iMinNumSpecies;
int iNumCoreSpecies;

int LenReq;
int Idx;
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
char szSpeciesList[512];
int iNumSpecies;
int iRegLen;
int iMinUltraLen;				// minimum ultra core length required

// command line args
struct arg_lit  *help    = arg_lit0("hH","help",                "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_int *ScreenLogLevel=arg_int0("S", "ScreenLogLevel",	"<int>","Level of diagnostics written to logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",	"diagnostics log file");

struct arg_int  *ProcMode = arg_int0("x","procmode","<int>",	"processing mode 0:default");
struct arg_int  *RegLen = arg_int0("L","updnstream","<int>",	"length of 5'up or 3'down  stream regulatory region length (default = 2000) 0..1000000");
struct arg_file *InFile = arg_file1("i",NULL,"<file>",			"input from .algn file");
struct arg_file *OutFile = arg_file0("o",NULL,"<file>",		    "output hypercore loci to file as CSV");
struct arg_file *InBedFile = arg_file0("b","bed","<file>",		"characterise regions from biobed file");
struct arg_int  *NumCoreSpecies = arg_int0("m","corespecies","<int>", "number of core species required in alignment, default is num of species in -s<specieslist>");
struct arg_int  *MinNumSpecies = arg_int0("M","minspecies","<int>", "minimum number ( >= NumCoreSpecies) of species required in alignment, default is number species specified - min number of core species");
struct arg_str  *SpeciesList = arg_str0("s","species","<string>","species list, ordered by processing priority (if not specified then all species in alignments");
struct arg_int  *MinUltraLen = arg_int0("n","minultralen","<int>","minimum (default = 50) required ultra length (10..1000) - will be maximally extended");


struct arg_file *ExcludeFile = arg_filen("E","exclude","<file>",0,cMaxExcludeFiles,	"exclude all regions in biobed file from processing ");
struct arg_file *IncludeFile = arg_filen("I","include","<file>",0,cMaxExcludeFiles,	"include all regions (unless specific regions excluded) in biobed file");
struct arg_str  *IncludeChroms = arg_strn("z","chromeinclude","<string>",0,cMaxIncludeChroms,"low priority - regular expressions defining species.chromosomes to include for processing");
struct arg_str  *ExcludeChroms = arg_strn("Z","chromexclude","<string>",0,cMaxExcludeChroms,"high priority - regular expressions defining species.chromosomes to exclude from processing");

struct arg_end *end = arg_end(20);

void *argtable[] = {help,version,FileLogLevel,ScreenLogLevel,LogFile,
					InFile,OutFile,InBedFile,
					ProcMode,NumCoreSpecies,MinNumSpecies,SpeciesList,
					MinUltraLen,RegLen,
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
	if(iProcMode < eProcModeStandard || iProcMode > eProcModeStandard)
		{
		printf("\nError: Requested processing mode '-x%d' not supported",iProcMode);
		exit(1);
		}

	iMinUltraLen = MinUltraLen->count ? MinUltraLen->ival[0] : cDfltUltraLen;
	if(iMinUltraLen < cMinUltraLen || iMinUltraLen > cMaxUltraLen)
		{
		printf("\nSpecified minimum ultra length was '-n%d', must be in range %d..%d",iMinUltraLen,cMinUltraLen,cMaxUltraLen);
		exit(1);;
		}

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
		printf("\nWarning: Species list empty, defaulting to all species in alignment file\n");
		szSpeciesList[0] = '\0';
		iNumSpecies = 0;
		}
	else
		{
		strcpy(szSpeciesList,SpeciesList->sval[0]);
		TrimQuotes(szSpeciesList);

		iNumSpecies = ParseNumSpecies(szSpeciesList,NULL);
	
		if(iNumSpecies < 2)
			{
			printf("\nError: At least two species must be specified\n");
			exit(1);
			}
		}
	
	iNumCoreSpecies = NumCoreSpecies->count ? NumCoreSpecies->ival[0] : iNumSpecies;
	if((iNumSpecies > 1 && iNumCoreSpecies < 1) || (iNumSpecies > 1 && iNumCoreSpecies > iNumSpecies))
		{
		printf("NumCoreSpecies '-m%d' must be in range 1..%d (number of species in species list)",iNumCoreSpecies,iNumSpecies);
		exit(1);
		}

	iMinNumSpecies = MinNumSpecies->count ? MinNumSpecies->ival[0] : iNumSpecies;
	if((iNumSpecies > 1 && iMinNumSpecies < iNumCoreSpecies) || (iNumSpecies > 1 && iMinNumSpecies > iNumSpecies))
		{
		printf("Error: MinNumSpecies '-M%d' must be in range %d..%d",iMinNumSpecies,iNumCoreSpecies,iNumSpecies);
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
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"bio multialignment (.algn) file to process: '%s'",szInputFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"where to write out stats: '%s'",szOutputFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"biobed file containing regional features: '%s'",szInputBiobedFile);
	switch(iProcMode) {
		case eProcModeStandard:
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Processing mode: standard");
			break;

		}

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"where to write out the ultracore loci: '%s'",szOutputFile);
	for(Idx = 0; Idx < NumIncludeFiles; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"biobed file containing regions to include: '%s'",pszIncludeFiles[Idx]); 
	for(Idx = 0; Idx < NumExcludeFiles; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"biobed file containing regions to exclude: '%s'",pszExcludeFiles[Idx]);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"number of core species to be in alignment: %d",iNumCoreSpecies);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Number of species including core species to be in alignment: %d",iMinNumSpecies);
	if(iNumSpecies > 1)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"species (%d) list: '%s'",iNumSpecies,szSpeciesList);
	else
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"species (%d) list: '%s'",iNumSpecies,"all species in alignments");
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"minimum ultra core length required: %d",iMinUltraLen);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"regulatory region length: %d",iRegLen);
	for(Idx = 0; Idx < NumIncludeChroms; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"reg expressions defining chroms to include: '%s'",pszIncludeChroms[Idx]);
	for(Idx = 0; Idx < NumExcludeChroms; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"reg expressions defining chroms to exclude: '%s'",pszExcludeChroms[Idx]); 

	gStopWatch.Start();
#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	Rslt = Process(iProcMode,		// processing mode 0: default
					szInputFile,		// bio multialignment (.algn) file to process
					szOutputFile,		// where to write out loci
					szInputBiobedFile,	// biobed file containing regional features - exons, introns etc
					NumIncludeFiles,	// number of include region files
					pszIncludeFiles,	// biobed files containing regions to include - default is to exclude none
					NumExcludeFiles,	// number of exclude region files
					pszExcludeFiles,	// biobed file containing regions to exclude - default is to include all
					iNumCoreSpecies,	// number of core species required to be in alignment
					iMinNumSpecies,		// minimum number of species including NumCoreSpecies to be in alignment
					szSpeciesList,		// space or comma separated list of species, priority ordered
					iMinUltraLen,		// minimum ultra core length required
					iRegLen,			// regulatory region length - up/dn stream of 5/3' 
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
	if(!Idx)	// reference species is always the first species in the species name list
		RefSpeciesID = SpeciesIDs[0];
	}


pAlignments->SetAllConfFilt(false);

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

		// not interested if the alignment length would be too short
	if(RefAlignLen < pProcParams->MinUltraLen)
		continue;

		// here we need to normalise the alignments so that there will be no case of all InDels in
		// any column which can occur if none of the processed species is not the reference species
	if((RefAlignLen = NormaliseInDelColumns(pProcParams,RefAlignLen))< pProcParams->MinUltraLen)
		continue;

	if(RefChromID != PrevRefChromID)
		PrevRefChromID = RefChromID;

	if(!ProcAlignBlock(RefChromID,RefChromOfs,RefAlignLen,pProcParams))
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
pProcParams->NumSpecies = 0;
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
	if((CurNumAlignments = pAlignments->GetNumSpecies(CurBlockID)) < pProcParams->MinNumSpecies)
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
	for(NumSpecies = Idx = 0; Idx < pProcParams->MaxNumSpecies; Idx++) 
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
		
	// check if required minimum number of species required were in multialignment
	if(NumSpecies < pProcParams->MinNumSpecies)	
		break;

	if(NumSpecies > pProcParams->NumSpecies)
		pProcParams->NumSpecies = NumSpecies;
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


NumSubSeqs = 0;
CurNumSeqs=0;
CurSeqLen = 0;
CurFiltSeqLen = 0;
SubRefOfs = RefChromOfs;
SubSeqStartIdx = 0;

for(SeqIdx = 0; SeqIdx < AlignLen; SeqIdx++)
	{
	bAllIdentical = true;			// assume all identical in column
	NumVInDels = 0;					// assume no InDel in column
	bRefInDel = false;				// and specifically that the reference base is not an InDel 
	for(VIdx = 0; VIdx < pProcParams->MinNumSpecies; VIdx++)
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

	// if in column any base was eBaseInDel then NumVInDels will be > 0
	// if the reference base was an InDel then bRefInDel will be true
	// if in column all identical and a,c,g or t then bAllIdentical will be true
	// if any base was mismatched or not a,c,g or t then bAllIdentical will be false
	if(NumVInDels)
		{
		if(CurFiltSeqLen >= pProcParams->MinUltraLen)	// has to be of at least MinLen to be worth further processing
			ProcessAlignment(RefChromID,SubRefOfs,SubSeqStartIdx,CurFiltSeqLen,pProcParams);
		CurSeqLen = 0;
		CurFiltSeqLen = 0;
		if(RefBase != eBaseInDel)
			RefChromOfs++;
		NumIdents = 0;
		SubSeqStartIdx = SeqIdx + 1;	// mark where next subsequence could start 
		SubRefOfs = RefChromOfs;		// chromosomal offset at which that next subsequence starts
		continue;
		}

	// no InDels in current aligned column - could still be mismatches - ,
	// continue to accept column as part of subsequence
	if(CurSeqLen == 0)		// if first base of putative subsequence...
		{
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
if(CurFiltSeqLen >= pProcParams->MinUltraLen)	// has to be at least MinUltraLen to be worth further processing
	ProcessAlignment(RefChromID,SubRefOfs,SubSeqStartIdx,CurFiltSeqLen,pProcParams);

return(true);
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
Process(int ProcMode,				// processing mode 0: default
				 char *pszInputFile,		// bio multialignment (.algn) file to process
					char *pszOutputFile,	// where to write out loci
					char *pszBiobedFile,	// biobed file containing regional features - exons, introns etc
					int NumIncludeFiles,	// number of include region files
					char **ppszIncludeFiles,// biobed files containing regions to include - default is to exclude none
					int NumExcludeFiles,	// number of exclude region files
					char **ppszExcludeFiles,// biobed file containing regions to exclude - default is to include all
					int NumCoreSpecies,		// number of core of species to be in alignment
					int MinNumSpecies,		// minimum number including core species to be in alignment
					char *pszSpeciesList,	// space or comma separated list of species, priority ordered
					int MinUltraLen,		// minimum ultra core length required
					int RegLen,				// regulatory region length - up/dn stream of 5/3' 
					int NumIncludeChroms,	// number of chromosomes explicitly defined to be included
					char **ppszIncludeChroms,	// ptr to array of reg expressions defining chroms to include - overides exclude
					int NumExcludeChroms,	// number of chromosomes explicitly defined to be excluded
					char **ppszExcludeChroms)	// ptr to array of reg expressions defining chroms to include
{
int Rslt;
int Idx;

int Regions;
int UpDnStreamLen;
CBEDfile *pBiobed = NULL;
tsProcParams ProcParams;
int NumSpecies;
char szCSVSpecies[2048];				// to hold comma separated species list

#ifndef _WIN32
int RegErr;				// regular expression parsing error
char szRegErr[128];			// to hold RegErr as textual representation ==> regerror();
#endif

memset(&ProcParams,0,sizeof(tsProcParams));

if(pszSpeciesList[0] == '\0')
	{
	CMAlignFile *pAlignments;
	int NumAlignedSpecies;
	char *pszSpecies;
	if((pAlignments = new CMAlignFile())==NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create new instance of CMAlignFile");
		return(eBSFerrObj);
		}
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Opening alignment file to get list of aligned species...");
	if((Rslt=pAlignments->Open(pszInputFile))!=eBSFSuccess)
		{
		while(pAlignments->NumErrMsgs())
			gDiagnostics.DiagOut(eDLFatal,gszProcName,pAlignments->GetErrMsg());
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open MAF file %s\n",pszInputFile);
		return(Rslt);
		}

	// extract species from multispecies alignment file
	szCSVSpecies[0]='\0';
	NumAlignedSpecies = pAlignments->GetNumSpecies();
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Alignment file contains %d aligned species...",NumAlignedSpecies);
	for(Idx = 0; Idx < NumAlignedSpecies; Idx++)
		{
		pszSpecies = pAlignments->GetSpeciesName(Idx+1);
		if(Idx > 0)
			strcat(szCSVSpecies,",");
		strcat(szCSVSpecies,pszSpecies);
		}
	ProcParams.pszSpeciesList = szCSVSpecies;
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Aligned species: '%s'",szCSVSpecies);
	pAlignments->Reset();
	delete pAlignments;
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Closed alignment file");

	NumSpecies = NumAlignedSpecies;

	if(NumCoreSpecies == 0 || NumCoreSpecies > NumSpecies)
		NumCoreSpecies = NumSpecies;
	if(MinNumSpecies == 0 || MinNumSpecies > NumSpecies)
		MinNumSpecies = NumSpecies;
	NumSpecies = ParseNumSpecies(szCSVSpecies,&ProcParams);
	}
else
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

ProcParams.ProcMode = ProcMode;
ProcParams.MinNumSpecies = MinNumSpecies;
ProcParams.NumCoreSpecies = NumCoreSpecies;
ProcParams.NumSpecies = NumSpecies;
ProcParams.MaxNumSpecies = NumSpecies;
ProcParams.MaxAlignIdxSpecies = 0;
ProcParams.RefSpeciesIdx = 0;
ProcParams.UpDnStreamLen = UpDnStreamLen;
ProcParams.MinUltraLen = MinUltraLen;
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
		CloseBedfiles(&ProcParams);
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
	    CloseBedfiles(&ProcParams);
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate memory (%d bytes) for holding species sequences",ProcParams.MaxSeqAlignLen);
		return(eBSFerrMem);
		}
	}

if(pszOutputFile != NULL && pszOutputFile[0] != '\0')
	{
#ifdef _WIN32
	if((ProcParams.hRsltsFile = open(pszOutputFile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
	if((ProcParams.hRsltsFile = open(pszOutputFile, O_RDWR | O_CREAT | O_TRUNC, S_IREAD | S_IWRITE ))==-1)
#endif
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create %s - %s",pszOutputFile,strerror(errno));
		close(ProcParams.hRsltsFile);
		CloseBedfiles(&ProcParams);
		return(eBSFerrCreateFile);
		}
	}
else
	ProcParams.hRsltsFile = -1;

Rslt = ProcessAlignments(pszInputFile,&ProcParams);

if(ProcParams.hRsltsFile != -1)
	close(ProcParams.hRsltsFile);

CloseBedfiles(&ProcParams);
for(Idx = 0; Idx < NumSpecies; Idx++)
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
// Process an alignment subsequence which may contain mismatches or 'N's but will not contain InDels 
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
while(SubSeqLen >= pProcParams->MinUltraLen)
	{
		// from current SubRefOfs maximally extend to right allowing for MaxMismatches
	NxtSeqIdx = ProcessSubSeq(RefChromID,ChromOfs,SeqIdx,SubSeqLen,pProcParams);
	if(NxtSeqIdx == -1)		// -1 flags that there is no more processing required for this alignment
		break;
	SubSeqLen -= (NxtSeqIdx - SeqIdx);	// update length of remaining sequence in this alignment
	if(SubSeqLen < pProcParams->MinUltraLen)
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
// and which ends immediately prior to the next InDel  or end of alignment block
// Returns the SubSeqStartIdx at which to start a subsequent ProcessSubSeq processing call
int	
ProcessSubSeq(int RefChromID,		// reference chromosome identifier
			  int ChromOfs,			// chromosome offset (0..n) at which this subsequence starts
			  int SeqIdx,			// index into pProcParams->pSeq[] at which this subsequence starts
			  int MaxLen,			// maximum length of this subsequence
			  tsProcParams *pProcParams)
{
int VIdx;

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
int CurUltraCoreLen = 0;		// current core with no mismatches
int MaxUltraCoreLen = 0;		// longest core with no mismatches encountered
int RefHyperCoreLen = 0;		// curent reference species core length (excludes any refseq InDels from length)
int WinIdx = 0;
int NumSpecies;

NumSpecies = pProcParams->NumSpecies;

for(HypercoreLen = 0; HypercoreLen < MaxLen ; HypercoreLen++, WinIdx++)
	{
	// determine in the current alignment column the number of mismatches
	VMismatches = 0;
	VIndels = 0;

	for(VIdx = 0; VIdx < NumSpecies; VIdx++)
		{
		pVSeq = pProcParams->pSeqs[VIdx];
		VBase = pVSeq[HypercoreLen+SeqIdx] & ~cRptMskFlg;
		
		if(VIdx == pProcParams->RefSpeciesIdx) // if base for reference sequence
			{
			RefBase = VBase;
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
	// number of mismatches in column now known - if more than allowed then characterise as a hyperconserved mismatch
	if(VMismatches >= 1)
		{
		CurUltraCoreLen = 0;						// any mismatch terminates current ultracore
		if(NxtSeqIdx == -1)						// if 1st mismatch then note where next search for hypercore should start from
			NxtSeqIdx = HypercoreLen+SeqIdx+1;  // if current hypercore not accepted
				// this many total mismatches allowed?
		break;
		}
	else	// no mismatch
		{
		CurUltraCoreLen++;							// current ultracore increases as bases identical
		if(CurUltraCoreLen > MaxUltraCoreLen)		// is this the longest ultracore in current subsequence?
			MaxUltraCoreLen = CurUltraCoreLen;
		}

	if(RefBase != eBaseInDel)
		RefHyperCoreLen++;
	}

// have a sequence which we may be interested in
if(HypercoreLen==MaxLen)// if core was terminated because subsequence ended then no point in subsequent
	NxtSeqIdx = -1;		// checking for more hypercores in same subsequence as these would be internal to this core

// RefHyperCoreLen holds the reference species hypercore length excluding any InDels on the ref species sequence
if(MaxUltraCoreLen < pProcParams->MinUltraLen)	// if less than required then return where to start next core from
	return(NxtSeqIdx);

ChromOfsEnd = ChromOfs+RefHyperCoreLen-1;	// determine reference chromosomal offset at which sequence ends

// ensure that the hypercore is in an included region and not part of an excluded region
if(!IncludeFilter(RefChromID,ChromOfs,ChromOfsEnd,pProcParams))
		return(NxtSeqIdx);

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
	
	if(RegionIdx != 0)		// if not intergenic then check for splice sites
		{
		SpliceSiteOverlaps = pProcParams->pBiobed->GetSpliceSiteBits(pProcParams->BEDChromID,ChromOfs,ChromOfsEnd,cMinSpliceOverlap);
		}
	}

OutputUltraCore(pProcParams->szRefChrom,ChromOfs,ChromOfsEnd,FeatureBits | SpliceSiteOverlaps,pProcParams);
if(NxtSeqIdx != -1)
	NxtSeqIdx = HypercoreLen+SeqIdx+1;
return(NxtSeqIdx);
}



bool	
OutputUltraCore(char *pszChrom, int ChromStartOffset, int ChromEndOffset, 
				int FeatureBits,		// feature bits over lapped
				tsProcParams *pProcParams)
{
static int CoreID = 0;
char szLineBuff[4096];
int Len;
if(pProcParams->hRsltsFile != -1)
	{
	Len = sprintf(szLineBuff,"%d,\"%s\",\"%s\",\"%s\",%d,%d,%d,\"%s\",%d",
		++CoreID,"ultracore",
		pProcParams->szSpecies[pProcParams->RefSpeciesIdx],
		pszChrom,ChromStartOffset,ChromEndOffset,ChromEndOffset-ChromStartOffset+1,
		pProcParams->pszSpeciesList,FeatureBits & (cAnyFeatBits | cOverlaysSpliceSites));
	Len += sprintf(&szLineBuff[Len],"\n");
	CUtility::SafeWrite(pProcParams->hRsltsFile,szLineBuff,Len);
	}
return(true);
}


