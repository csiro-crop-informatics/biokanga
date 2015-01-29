// quickcount.cpp : Defines the entry point for the console application
// Counts the NMer distributions from mononucleotide upto user specified NMer limit <= cMaxNMerLen
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

const unsigned int cProgVer = 106;		// increment program version with each release

const int cMaxNMerLen = 15;				// maximum length N-Mer (don't change unless code dependencies checked and plenty of memory!)
const int cDfltNMerLen= 5;				// default N-Mer length

const int cMaxIncludeFiles = 10;		// maximum number of include region filter files
const int cMaxExcludeFiles = 10;		// maximum number of exclude region filter files

const int cMaxIncludeChroms = 20;		// max number of include chromosomes regular expressions
const int cMaxExcludeChroms = 20;		// max number of exclude chromosomes regular expressions

const int cChromSeqLen = 0x03fffff;		// process sequences in this size chunks

const int cDensityMB = 1000000;			// normalise instance density to number of instances per million nucleotides

// processing mode
typedef enum eProcMode {
	eProcModeNMerDistAllSeqs,			// process N-Mer frequency distribution over all sequences
	eProcModeNMerDistPerSeq,			// process per sequence N-Mer frequency distribution
	eProcModeNMerDistNorm				// process N-Mer frequence distribution over all sequences normalising for sequence lengths
} etProcMode;


const int cMaxExcludeHistory = 100;
typedef struct TAG_sExcludeEl {
	struct TAG_sExcludeEl *pNext;
	struct TAG_sExcludeEl *pPrev;
	int ChromID;		// identifies chromosome
	bool bExclude;		// true if to be excluded, false if not
	} tsExcludeEl;

#pragma pack(1)
typedef struct TAG_sNormDensity {
	UINT32 Instances;			// absolute number of instances for current chromosome (later used to hold sequence identifier)
	double SumDensities;		// densities for each chromosome summed (later used to hold means)
	double SumDensitiesSquared;	// densities for each chromosome squared and then summed (later used to hold stddev)
	double COV;					// holds coefficient of variation
} tsNormDensity;
#pragma pack()

typedef struct TAG_sProcParams 
	{
	etProcMode ProcMode;			// processing mode
	int MaxNMerLen;					// max N-Mer to process
	int MinNMerLen;					// min N-Mer to process
	char szCurChrom[cMaxDatasetSpeciesChrom];
	int CurChromID;					// current chromosome identifier
	etSeqBase *pSeq;				// ptr to sequence being processed
	int CurSeqLen;					// actual current sequence length
	int AllocdSeqLen;				// how many bytes of memory were alloc'd to pSeq
	int *pCntStepCnts;				// array of stats counters organised in [Step][Region] order

	int NumSeqs;					// number of sequences or chromsomes processed
	int NumNormCnts;				// number of normalised stats counters
	tsNormDensity *pNormCnts;			// array of normalised stats counters

	int NumCntSteps;				// number of steps in pCntStepCnts
	int NumRegions;					// number of regions per step
	int CurRegion;					// NMer is clasified into this region
	int CntStepOfs[cMaxNMerLen];	// ofs into pCntStepCnts for each NMer length
	
	int BEDChromID;					// BED chromosome identifier corresponding to CurChromID
	int NumIncludes;				// number of biobed files containing regions to include
	int NumExcludes;				// number of biobed files containing regions to exclude
	CBEDfile *pIncludes[cMaxIncludeFiles];	// if opened biobed files for regions to include - all other regions are to be excluded
	CBEDfile *pExcludes[cMaxExcludeFiles];	// if opened biobed files for regions to exclude 
	int UpDnStreamLen;				// up/dn stream regional length when characterising
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
	char szInFile[_MAX_PATH];			// assembly or bioseq .fa
	CBioSeqFile *pBioSeqFile;			// bioseq class instance containing assembly or fasta sequences
	CFasta *pFasta;						// Fasta instance if fasta file as input
	char szOutFile[_MAX_PATH];			// where to write out stats
	int hRsltsFile;						// write results into this CSV file
	char szBiobedFile[_MAX_PATH];		// biobed file containing regional features - exons, introns etc
	CBEDfile *pBiobed;					// if not NULL then opened biobed file for regional characteristics
	} tsProcParams; 

int TrimQuotes(char *pszTxt);
int Process(etProcMode ProcMode,		// processing mode
			int MaxNMerLen,				// max length N-Mer to process
			int MinNMerLen,				// min N-Mer to process
			char *pszInFile,			// input assembly or bioseq .fa
			char *pszOutFile,			// where to write out stats
			char *pszBiobedFile,		// biobed file containing regional features - exons, introns etc
			int	NumIncludeFiles,		// number of include region files
			char **ppszIncludeFiles,	// biobed files containing regions to include - default is to exclude none
			int	NumExcludeFiles,		// number of exclude region files
			char **ppszExcludeFiles,	// biobed file containing regions to exclude - default is to include all
			int		RegLen,				// regulatory region length - up/dn stream of 5/3' 
			int		NumIncludeChroms,	// number of chromosomes explicitly defined to be included
			char **ppszIncludeChroms,	// ptr to array of reg expressions defining chroms to include - overides exclude
			int		NumExcludeChroms,	// number of chromosomes explicitly defined to be excluded
			char **ppszExcludeChroms);	// ptr to array of reg expressions defining chroms to include
	
int TrimQuotes(char *pszTxt);
CBEDfile *OpenBedfile(char *pToOpen);
bool CleanupResources(tsProcParams *pProcParams);
int OutputResults(tsProcParams *pProcParams);
int OutputDensityResults(tsProcParams *pProcParams);
int GenBioseqFreqCounts(char *pszInFile,tsProcParams *pProcParams);
int ProcessSequence(tsProcParams *pProcParams);
bool ExcludeThisChrom(tsProcParams *pProcParams);
tsExcludeEl *LocateExclude(int ChromID);
bool AddExcludeHistory(int ChromID,bool bExclude);
char *StepIdx2Seq(int SeqLen,int SeqIdx);
int GenSeqIdx(int SeqLen,etSeqBase *pSeq);
bool IncludeFilter(int SubRefOfs,int SubRefEndOfs,tsProcParams *pProcParams);
int ProcessSequenceDensity(tsProcParams *pProcParams);
static int SortNMerSumDensities( const void *arg1, const void *arg2);
static int SortNMerSumDensitiesSquared( const void *arg1, const void *arg2);
static int SortNMerCOV( const void *arg1, const void *arg2);


CStopWatch gStopWatch;
CDiagnostics gDiagnostics;				// for writing diagnostics messages to log file
char gszProcName[_MAX_FNAME];			// process name

int gNumExcludeEls = 0;					// current number of elements in gExcludeChroms
tsExcludeEl *gpMRA = NULL;				// pts to most recently accessed or added
tsExcludeEl *gpLRA = NULL;				// pts to least recently accessed
tsExcludeEl gExcludeChroms[cMaxExcludeHistory];


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
	return _T("quickcount");
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
int Idx;
int LenReq;
int Rslt;
int iProcMode;
int iMaxNMerLen;
int iMinNMerLen;
char szInFile[_MAX_PATH];
char szOutFile[_MAX_PATH];
char szInBedFile[_MAX_PATH];
int iRegLen;
int NumIncludeChroms;
char *pszIncludeChroms[cMaxIncludeChroms];
int NumExcludeChroms;
char *pszExcludeChroms[cMaxExcludeChroms];
int NumIncludeFiles;
char *pszIncludeFiles[cMaxExcludeFiles];
int NumExcludeFiles;
char *pszExcludeFiles[cMaxIncludeFiles];

// command line args
struct arg_lit  *help    = arg_lit0("hH","help",                "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_int *ScreenLogLevel=arg_int0("S", "ScreenLogLevel",	"<int>","Level of diagnostics written to logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",	"diagnostics log file");

struct arg_int  *ProcMode = arg_int0("m","mode","<int>",		"processing mode: 0 = N-Mer distribution, 1 = N-Mer distribution per sequence, 2 = N-Mer normalised extreme distributions");
struct arg_file *InFile = arg_file1("i","infile","<file>",		"assembly file to process (either fasta or bioseq)");
struct arg_int *MinNMerLen=arg_int0("l", "minnmerlen","<int>",  "process N-Mers from this length up (default is 1)");
struct arg_int *MaxNMerLen=arg_int0("L", "maxnmerlen","<int>",  "process N-Mers upto and including this length (default is 5)");
struct arg_file *OutFile = arg_file1("o",NULL,"<file>",			"output results to CSV file");
struct arg_file *InBedFile = arg_file0("b","bed","<file>",		"characterise regions from biobed file");
struct arg_int  *RegLen = arg_int0("L","updnstream","<int>",	"length of 5'up or 3'down  stream regulatory region length (default = 2000) 0..1000000");
struct arg_file *ExcludeFile = arg_filen("E","exclude","<file>",0,cMaxExcludeFiles,	"exclude all regions in biobed file from processing ");
struct arg_file *IncludeFile = arg_filen("I","include","<file>",0,cMaxExcludeFiles,	"include all regions (unless specific regions excluded) in biobed file");
struct arg_str  *IncludeChroms = arg_strn("z","chromeinclude","<string>",0,cMaxIncludeChroms,"low priority - regular expressions defining species.chromosomes to include for processing");
struct arg_str  *ExcludeChroms = arg_strn("Z","chromexclude","<string>",0,cMaxExcludeChroms,"high priority - regular expressions defining species.chromosomes to exclude from processing");
struct arg_end *end = arg_end(20);

void *argtable[] = {help,version,FileLogLevel,ScreenLogLevel,LogFile,
						ProcMode,MinNMerLen,MaxNMerLen,InFile,OutFile,InBedFile,RegLen,ExcludeFile,IncludeFile,
						IncludeChroms,ExcludeChroms,
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
		printf("\n%s Version: %d.%2.2d\n",gszProcName,cProgVer/100,cProgVer%100);
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

	iProcMode = ProcMode->count ? ProcMode->ival[0] : eProcModeNMerDistAllSeqs;
	if(iProcMode < eProcModeNMerDistAllSeqs || iProcMode > eProcModeNMerDistNorm)
		{
		printf("Error: Unsupported processing mode '-m%d'\n",iProcMode);
		exit(1);
		}

	iMinNMerLen = MinNMerLen->count ? MinNMerLen->ival[0] : 1;
	if(iMinNMerLen < 1 || iMinNMerLen > cMaxNMerLen)
		{
		printf("Error: min N-Mer length requested '-l%d' outside of range %d..%d\n",iMinNMerLen,1,cMaxNMerLen);
		exit(1);
		}
	if(iMinNMerLen > cDfltNMerLen)
		iMaxNMerLen = iMinNMerLen;
	else
		iMaxNMerLen = cDfltNMerLen;



	if(iProcMode == eProcModeNMerDistNorm)
		{
		printf("\nIn normalised distribution mode '-m2' the maximum NMer length will be set to be same as minimim NMer length");
		iMaxNMerLen = iMinNMerLen;
		}
	else
		{
		iMaxNMerLen = MaxNMerLen->count ? MaxNMerLen->ival[0] : iMaxNMerLen;
		if(iMaxNMerLen < iMinNMerLen || iMaxNMerLen > cMaxNMerLen)
			{
			printf("Error: max N-Mer length requested '-n%d' outside of range %d..%d\n",iMaxNMerLen,iMinNMerLen,cMaxNMerLen);
			exit(1);
			}
		}

	if(InBedFile->count)
		strcpy(szInBedFile,InBedFile->filename[0]);
	else
		szInBedFile[0] = '\0';

	
	strncpy(szInFile,InFile->filename[0],sizeof(szInFile));
	szInFile[sizeof(szInFile)-1] = '\0';
	strncpy(szOutFile,OutFile->filename[0],sizeof(szOutFile));
	szOutFile[sizeof(szOutFile)-1] = '\0';

	iRegLen = RegLen->count ? RegLen->ival[0] : cDfltRegLen;
	if(iRegLen < cMinRegLen)
		{
			printf("Warning: Regulatory region length '-L%d' less than minimum %d, assuming you meant to use '-L%d'\n",iRegLen,cMinRegLen,cMinRegLen);
		iRegLen = cMinRegLen;
		}
	else
		{
		if(iRegLen > cMaxRegLen)
			{
				printf("Warning: Regulatory region length '-L%d' more than maximum %d, assuming you meant to use '-L%d'\n",iRegLen,cMaxRegLen,cMaxRegLen);
			iRegLen = cMaxRegLen;
			}
		}

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
		case eProcModeNMerDistAllSeqs:			
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Processing Mode: N-Mer frequency distribution over all sequences");
			break;
		case eProcModeNMerDistPerSeq:			
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Processing Mode: N-Mer frequency distribution for each sequences");
			break;
		case eProcModeNMerDistNorm:			
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"process N-Mer density distribution over all sequences or chromosomes");
			break;
		}

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Minimum N-Mer length: %d",iMinNMerLen);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Maximum N-Mer length: %d",iMaxNMerLen);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Input bioseq assembly/sequence file: '%s'",szInFile);
	if(szInBedFile[0] != '\0')
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Input biobed region file: '%s'",szInBedFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Results output to file: '%s'",szOutFile);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"regulatory region length: %d",iRegLen);
	for(Idx = 0; Idx < NumIncludeChroms; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"reg expressions defining chroms to include: '%s'",pszIncludeChroms[Idx]);
	for(Idx = 0; Idx < NumExcludeChroms; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"reg expressions defining chroms to exclude: '%s'",pszExcludeChroms[Idx]); 

	gStopWatch.Start();
#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	Rslt = Process((etProcMode)iProcMode,	// processing mode 0: default
					iMinNMerLen,		// min length N-Mer to process
					iMaxNMerLen,		// max length N-Mer to process
					szInFile,			// assembly or bioseq .fa
					szOutFile,			// where to write out stats
					szInBedFile,		// biobed file containing regional features - exons, introns etc
					NumIncludeFiles,	// number of include region files
					pszIncludeFiles,	// biobed files containing regions to include - default is to exclude none
					NumExcludeFiles,	// number of exclude region files
					pszExcludeFiles,	// biobed file containing regions to exclude - default is to include all
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
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open biobed file %s",pToOpen);
		delete pBed;
		return(NULL);
		}
	return(pBed);
	}
return(NULL);
}

//CleanupResources
//Closes and deletes all created and opened Biobed files
bool
CleanupResources(tsProcParams *pProcParams)
{
int Idx;
if(pProcParams->pBioSeqFile != NULL)
	{
	delete pProcParams->pBioSeqFile;
	pProcParams->pBioSeqFile = NULL;
	}

if(pProcParams->pFasta != NULL)
	{
	delete pProcParams->pFasta;
	pProcParams->pFasta = NULL;
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
if(pProcParams->pBiobed != NULL)
	delete pProcParams->pBiobed;
pProcParams->pBiobed = NULL;

if(pProcParams->pCntStepCnts != NULL)
	{
	delete pProcParams->pCntStepCnts;
	pProcParams->pCntStepCnts = NULL;
	}

if(pProcParams->pNormCnts != NULL)
	{
	delete pProcParams->pNormCnts;
	pProcParams->pNormCnts = NULL;
	}

if(pProcParams->pSeq != NULL)
	{
	delete pProcParams->pSeq;
	pProcParams->pSeq = NULL;
	}
pProcParams->CurSeqLen = 0;
pProcParams->AllocdSeqLen = 0;

if(pProcParams->hRsltsFile != -1)
	{
	close(pProcParams->hRsltsFile);
	pProcParams->hRsltsFile = -1;
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
		if((BEDChromID = pProcParams->pIncludes[Idx]->LocateChromIDbyName(pProcParams->szCurChrom))<1)
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
		if((BEDChromID = pProcParams->pExcludes[Idx]->LocateChromIDbyName(pProcParams->szCurChrom))<1)
			continue;
		if(pProcParams->pExcludes[Idx]->InAnyFeature(BEDChromID,SubRefOfs,SubRefEndOfs))
			return(false);
		}
	}
return(true);
}

int
GenSeqIdx(int SeqLen,etSeqBase *pSeq)
{
int Idx;
int SeqIdx;
int Base;
for(Idx=SeqIdx=0; Idx < SeqLen; Idx++,pSeq++)
	{
	Base = *pSeq & ~cRptMskFlg;
	if(Base > eBaseT)
		return(-1);
	SeqIdx <<= 2;
	SeqIdx |= Base;
	}
return(SeqIdx);
}

char *
StepIdx2Seq(int SeqLen,int SeqIdx)
{
static char szSeqBuff[256];
char *pChr;
int Base;
int Idx;

szSeqBuff[SeqLen] = '\0';
pChr = &szSeqBuff[SeqLen-1];
for(Idx=0; Idx < SeqLen; Idx++,pChr--)
	{
	Base = SeqIdx & 0x03;
	switch(Base) {
		case 0: *pChr = 'a'; break;
		case 1: *pChr = 'c'; break;
		case 2: *pChr = 'g'; break;
		case 3: *pChr = 't'; break;
		}
	SeqIdx >>= 2;
	}
return(szSeqBuff);
}




// AddExcludeHistory
// Adds a tsExcludeEl to the cached history
// The newly added element will be the MRA
bool
AddExcludeHistory(int ChromID,bool bExclude)
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
pEl->ChromID = ChromID;
return(bExclude);
}

// LocateExclude
// Locates - starting from the MRA - a tsExcludeEl which matches on SpeciesID and ChromID
// If matches then this tsExcludeEl is made the MRA
// Returns ptr to matching tsExcludeEl if located or NULL
tsExcludeEl *
LocateExclude(int ChromID)
{
tsExcludeEl *pEl;
pEl = gpMRA;
while(pEl != NULL)
	{
	if(pEl->ChromID == ChromID)
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
ExcludeThisChrom(tsProcParams *pProcParams)
{
#ifdef _WIN32
RegexpMatch mc;
#else
regmatch_t mc;
#endif

int Idx;
if(!pProcParams->NumExcludeChroms && !pProcParams->NumIncludeChroms)
	return(false);
// check if this chromosome are already known to be included/excluded
tsExcludeEl *pEl;
if((pEl = LocateExclude(pProcParams->CurChromID))!=NULL)
	return(pEl->bExclude);
// haven't seen this chromosome before - or else they have been discarded from history...

// to be excluded?
for(Idx = 0; Idx < pProcParams->NumExcludeChroms; Idx++)
#ifdef _WIN32	
	if(pProcParams->ExcludeChromsRE[Idx]->Match(pProcParams->szCurChrom,&mc))
#else
	if(!regexec(&pProcParams->ExcludeChromsRE[Idx],pProcParams->szCurChrom,1,&mc,0))
#endif
		return(AddExcludeHistory(pProcParams->CurChromID,true));

// to be included?
for(Idx = 0; Idx < pProcParams->NumIncludeChroms; Idx++)
	{
#ifdef _WIN32
	if(pProcParams->IncludeChromsRE[Idx]->Match(pProcParams->szCurChrom,&mc))
#else
	if(!regexec(&pProcParams->IncludeChromsRE[Idx],pProcParams->szCurChrom,1,&mc,0))
#endif
		return(AddExcludeHistory(pProcParams->CurChromID,false));
	}

return(AddExcludeHistory(pProcParams->CurChromID,pProcParams->NumIncludeChroms > 0 ? true : false));
}


int
ProcessSequence(tsProcParams *pProcParams)
{
int *pCnt;
int Idx;
int FeatureBits;
int Region;
int SeqIdx;
int CurSeqIdxLen;
bool bChkInclude;
bool bChkRegion;

for(Idx = 0; Idx < pProcParams->CurSeqLen; Idx++)
	{
	bChkInclude = true;
	bChkRegion = true;
	for(CurSeqIdxLen = pProcParams->MaxNMerLen; CurSeqIdxLen >= pProcParams->MinNMerLen; CurSeqIdxLen--)
		{
		if((Idx + CurSeqIdxLen) > pProcParams->CurSeqLen)
			continue;
		SeqIdx = GenSeqIdx(CurSeqIdxLen,&pProcParams->pSeq[Idx]);
		if(SeqIdx < 0)		// < 0 generally indicates that sequences contains 'N' or some other none a,c,g,t base
			continue;

		// ensure that the hypercore is in an included region and not part of an excluded region
		if(bChkInclude && !IncludeFilter(Idx,Idx+CurSeqIdxLen-1,pProcParams))
			continue;
		bChkInclude = false;	// if longest N-Mer is included then shorter will always be!

		// which region does this N-Mer overlap?
		if(bChkRegion)
			{
			if(pProcParams->BEDChromID > 0)
				FeatureBits = pProcParams->pBiobed->GetFeatureBits(pProcParams->BEDChromID,Idx,Idx+CurSeqIdxLen-1,cRegionFeatBits,pProcParams->UpDnStreamLen);
			else
				FeatureBits = 0;
			Region = pProcParams->pBiobed->MapFeatureBits2Idx(FeatureBits);
			if(Region == 0 || Region == 4)	// if longer N-Mer contained within IG or Intron then shorter will also be
				bChkRegion = false;
			}

		pProcParams->CurRegion = Region;
		pCnt = pProcParams->pCntStepCnts;
		pCnt += pProcParams->CntStepOfs[CurSeqIdxLen-1];
		pCnt[(SeqIdx * pProcParams->NumRegions) + Region] += 1;
		}
	}
return(eBSFSuccess);
}

int
ProcessSequenceDensity(tsProcParams *pProcParams)
{
int Idx;
int SeqIdx;
for(Idx = 0; Idx <= (pProcParams->CurSeqLen - pProcParams->MaxNMerLen); Idx++)
	{
	SeqIdx = GenSeqIdx(pProcParams->MaxNMerLen,&pProcParams->pSeq[Idx]);
	if(SeqIdx < 0)		// < 0 generally indicates that sequences contains 'N' or some other none a,c,g,t base
		continue;

	// ensure that the hypercore is in an included region and not part of an excluded region
	if(!IncludeFilter(Idx,Idx+pProcParams->MaxNMerLen-1,pProcParams))
		continue;
	pProcParams->pNormCnts[SeqIdx].Instances += 1;
	}
return(eBSFSuccess);
}


int
GenBioseqFreqCountsBioseq(tsProcParams *pProcParams)
{
int Rslt;
int CurEntryID;
int CurSeqLen;

CurEntryID = 0;
while((CurEntryID = pProcParams->pBioSeqFile->Next(CurEntryID))>0)
	{
	pProcParams->CurChromID = CurEntryID;
	pProcParams->pBioSeqFile->GetName(CurEntryID,sizeof(pProcParams->szCurChrom),pProcParams->szCurChrom);
	if(ExcludeThisChrom(pProcParams))	// should this chromosome be accepted for processing or sloughed?
		{
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing %s...Skipped",pProcParams->szCurChrom);
		continue;	
		}
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing %s...",pProcParams->szCurChrom);
	CurSeqLen = pProcParams->pBioSeqFile->GetDataLen(CurEntryID);
	if(CurSeqLen > pProcParams->AllocdSeqLen)
		{
		delete pProcParams->pSeq;
		if((pProcParams->pSeq = new unsigned char [CurSeqLen+0x07fff])==NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to alloc memory (%d requested) for holding raw sequence data",CurSeqLen+0x07fff);
			Rslt = eBSFerrMem;
			break;
			}
		pProcParams->AllocdSeqLen = CurSeqLen+0x07fff;
		}

	if((Rslt=pProcParams->pBioSeqFile->GetData(CurEntryID,eSeqBaseType,0,pProcParams->pSeq,CurSeqLen))<eBSFSuccess)
		{
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Unable to retrieve sequence for %s...Skipped",pProcParams->szCurChrom);
		continue;
		}
	if(pProcParams->pBiobed != NULL)
		{
		pProcParams->BEDChromID = pProcParams->pBiobed->LocateChromIDbyName(pProcParams->szCurChrom);
		if(pProcParams->BEDChromID < 1)
			{
			gDiagnostics.DiagOut(eDLWarn,gszProcName,"\nUnable to locate chromosome '%s' in biobed file, unable to generate regional stats...\n",pProcParams->szCurChrom);
			pProcParams->BEDChromID = 0;
			}
		}
	else
		pProcParams->BEDChromID = -1;
	
	if(pProcParams->ProcMode == eProcModeNMerDistPerSeq) 
		memset(pProcParams->pCntStepCnts,0,pProcParams->NumCntSteps * sizeof(int));
	else
		if(pProcParams->ProcMode == eProcModeNMerDistNorm)
			{
			tsNormDensity *pNormCnt = pProcParams->pNormCnts;
			for(int NormIdx = 0; NormIdx < pProcParams->NumNormCnts; NormIdx++,pNormCnt++)
				pNormCnt->Instances = 0;
			}

	pProcParams->CurSeqLen = CurSeqLen;
	pProcParams->NumSeqs += 1;

	if(pProcParams->ProcMode == eProcModeNMerDistNorm)
		{
		double Density;
		if((Rslt=ProcessSequenceDensity(pProcParams)) < 0)
			break;

		tsNormDensity *pNormCnt = pProcParams->pNormCnts;
		for(int NormIdx = 0; NormIdx < pProcParams->NumNormCnts; NormIdx++,pNormCnt++)
			{
			Density = ((double)cDensityMB * pNormCnt->Instances) / CurSeqLen;
			pNormCnt->SumDensities += Density;
			pNormCnt->SumDensitiesSquared += Density * Density;
			pNormCnt->Instances = 0;
			}
		}
	else
		if((Rslt=ProcessSequence(pProcParams)) < 0)
			break;

	if(pProcParams->ProcMode == eProcModeNMerDistPerSeq && 
		(Rslt = OutputResults(pProcParams)) < 0)
		break;
	}

return(Rslt);
}


int
GenBioseqFreqCountsFasta(tsProcParams *pProcParams)
{
int Rslt;
int CurEntryID;
int NumSeqs;
int DisplayEveryN;
int MaxSeqLen;
int CurSeqLen;
bool bDescriptor;
char szDescription[cBSFDescriptionSize];
CFasta *pFasta = pProcParams->pFasta;

// ensure allocd for largest sequence in fasta
pFasta->Reset();
bDescriptor = false;
MaxSeqLen = 0;
NumSeqs = 0;
while((Rslt = CurSeqLen = pFasta->ReadSequence()) > eBSFSuccess)
	{
	if(CurSeqLen == eBSFFastaDescr)		// just read a descriptor line
		continue;
	if(CurSeqLen > MaxSeqLen)
		MaxSeqLen = CurSeqLen;
	NumSeqs += 1;
	}
if(pProcParams->pSeq == NULL || pProcParams->AllocdSeqLen < MaxSeqLen)
	{
	if(pProcParams->pSeq != NULL)
		delete pProcParams->pSeq;
	pProcParams->AllocdSeqLen = 0;
	if((pProcParams->pSeq = new unsigned char [MaxSeqLen])==NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate memory (%d bytes) for seq processing ",MaxSeqLen);
		return(eBSFerrMem);
		}
	pProcParams->AllocdSeqLen = MaxSeqLen;
	}

pFasta->Reset();

if(NumSeqs <= 30)
	DisplayEveryN = 1;
else
	DisplayEveryN = (2 * NumSeqs) / 30;

CurEntryID = 0;
bDescriptor = false;
while((Rslt = CurSeqLen = pFasta->ReadSequence(pProcParams->pSeq,MaxSeqLen)) > eBSFSuccess)
	{
	if(Rslt == eBSFFastaDescr)		// just read a descriptor line
		{
		pFasta->ReadDescriptor(szDescription,cBSFDescriptionSize);
		bDescriptor = true;
		continue;
		}
	else									// just read sequence
		if(!bDescriptor)					// if there was no descriptor then dummy up one...
			{
			sprintf(szDescription,"Probe%d",CurEntryID+1);
			bDescriptor = true;
			}
	strncpy(pProcParams->szCurChrom,szDescription,cMaxDatasetSpeciesChrom);
	pProcParams->szCurChrom[cMaxDatasetSpeciesChrom-1] = '\0';
	pProcParams->CurChromID = ++CurEntryID;

	if(CurEntryID == 1 || !(CurEntryID % DisplayEveryN))
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing %s...",pProcParams->szCurChrom);
	pProcParams->BEDChromID = -1;

	if(pProcParams->ProcMode == eProcModeNMerDistPerSeq) 
		memset(pProcParams->pCntStepCnts,0,pProcParams->NumCntSteps * sizeof(int));
	else
		if(pProcParams->ProcMode == eProcModeNMerDistNorm)
			{
			tsNormDensity *pNormCnt = pProcParams->pNormCnts;
			for(int NormIdx = 0; NormIdx < pProcParams->NumNormCnts; NormIdx++,pNormCnt++)
				pNormCnt->Instances = 0;
			}

	pProcParams->CurSeqLen = CurSeqLen;
	pProcParams->NumSeqs += 1;
	if((Rslt=ProcessSequence(pProcParams)) < 0)
		break;

	if(pProcParams->ProcMode == eProcModeNMerDistNorm)
		{
		double Density;
		tsNormDensity *pNormCnt = pProcParams->pNormCnts;
		for(int NormIdx = 0; NormIdx < pProcParams->NumNormCnts; NormIdx++,pNormCnt++)
			{
			Density = ((double)cDensityMB * pNormCnt->Instances) / CurSeqLen;
			pNormCnt->SumDensities += Density;
			pNormCnt->SumDensitiesSquared += Density * Density;
			pNormCnt->Instances = 0;
			}
		}

	if(pProcParams->ProcMode == eProcModeNMerDistPerSeq && 
		(Rslt = OutputResults(pProcParams)) < 0)
		break;
	}

return(Rslt);
}

int
GenBioseqFreqCounts(char *pszInFile,					// either bioseq or multifasta
					tsProcParams *pProcParams)
{
int Rslt;
char szTitleLine[512];
int Len;

if((pProcParams->pSeq = new unsigned char [cChromSeqLen])==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate memory (%d bytes) for seq processing ",cChromSeqLen);
	return(eBSFerrMem);
	}
pProcParams->AllocdSeqLen = cChromSeqLen;
if(pProcParams->pCntStepCnts != NULL)
	memset(pProcParams->pCntStepCnts,0,pProcParams->NumCntSteps * sizeof(int));

// Initially try opening as fasta, if that fails then open as bioseq
pProcParams->pFasta = new CFasta();
if((Rslt=pProcParams->pFasta->Open(pszInFile)) == eBSFSuccess)
	Rslt = GenBioseqFreqCountsFasta(pProcParams);
else
	{
	if(Rslt != eBSFerrNotFasta)		
		{
		while(pProcParams->pFasta->NumErrMsgs())
			gDiagnostics.DiagOut(eDLFatal,gszProcName,pProcParams->pFasta->GetErrMsg());
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open file '%s' for processing",pszInFile);
		}
	else
		{
		pProcParams->pFasta->ClearErrs();
		pProcParams->pBioSeqFile = new CBioSeqFile();
		if((Rslt=pProcParams->pBioSeqFile->Open(pszInFile,cBSFTypeAny,false))==eBSFSuccess)
			Rslt = GenBioseqFreqCountsBioseq(pProcParams);
		if(Rslt < eBSFSuccess)
			{
			while(pProcParams->pBioSeqFile->NumErrMsgs())
				gDiagnostics.DiagOut(eDLFatal,gszProcName,pProcParams->pBioSeqFile->GetErrMsg());
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open file '%s' for processing",pszInFile);
			}
		}
	}

if(Rslt >= eBSFSuccess)
	{
	if(pProcParams->ProcMode == eProcModeNMerDistAllSeqs)
		{
		Len = sprintf(szTitleLine,"\"Genome\",\"N-Mer\",\"Ref\",\"Oligo\",\"Total\",\"IG\",\"5'US\",\"5'UTR\",\"CDS\",\"INTRON\",\"3'UTR\",\"3'DS\",\"5'ExSplice\",\"3'ExSplice\"\n");
		CUtility::SafeWrite(pProcParams->hRsltsFile,szTitleLine,Len);
		Rslt = OutputResults(pProcParams);
		}
	else
		if(pProcParams->ProcMode == eProcModeNMerDistNorm)
			OutputDensityResults(pProcParams);
	}

if(pProcParams->pFasta != NULL)
	{
	delete pProcParams->pFasta;
	pProcParams->pFasta = NULL;
	}
if(pProcParams->pBioSeqFile != NULL)
	{
	delete pProcParams->pBioSeqFile;
	pProcParams->pBioSeqFile = NULL;
	}
if(pProcParams->pSeq != NULL)
	{
	delete pProcParams->pSeq;
	pProcParams->pSeq = NULL;
	}
pProcParams->AllocdSeqLen = 0;
return(Rslt);
}



int Process(etProcMode ProcMode,		// processing mode 0: N-Mer frequency distribution
			int MinNMerLen,				// min length N-Mer to process
			int MaxNMerLen,				// max length N-Mer to process
			char *pszInFile,			// assembly or bioseq .fa
			char *pszOutFile,			// where to write out stats
			char *pszBiobedFile,		// biobed file containing regional features - exons, introns etc
			int	NumIncludeFiles,		// number of include region files
			char **ppszIncludeFiles,	// biobed files containing regions to include - default is to exclude none
			int	NumExcludeFiles,		// number of exclude region files
			char **ppszExcludeFiles,	// biobed file containing regions to exclude - default is to include all
			int		RegLen,				// regulatory region length - up/dn stream of 5/3' 
			int		NumIncludeChroms,	// number of chromosomes explicitly defined to be included
			char **ppszIncludeChroms,	// ptr to array of reg expressions defining chroms to include - overides exclude
			int		NumExcludeChroms,	// number of chromosomes explicitly defined to be excluded
			char **ppszExcludeChroms)	// ptr to array of reg expressions defining chroms to include
{
int Rslt;
int Idx;
int NumCnts;
int RegionCnts;
int NumRegions;
int UpDnStreamLen;
tsProcParams ProcParams;

#ifndef _WIN32
int RegErr;				// regular expression parsing error
char szRegErr[128];		// to hold RegErr as textual representation ==> regerror();
#endif

memset(&ProcParams,0,sizeof(tsProcParams));
if(pszInFile != NULL && pszInFile[0] != '\0')
	strcpy(ProcParams.szInFile,pszInFile);

if(pszBiobedFile != NULL && pszBiobedFile[0] != '\0')
	strcpy(ProcParams.szBiobedFile,pszBiobedFile);
if(pszOutFile != NULL && pszOutFile[0] != '\0')
	strcpy(ProcParams.szOutFile,pszOutFile);

for(Idx=0;Idx<NumIncludeFiles; Idx++)
	{
	if((ProcParams.pIncludes[Idx] = OpenBedfile(ppszIncludeFiles[Idx]))==NULL)
		{
		CleanupResources(&ProcParams);
		return(eBSFerrObj);
		}
	ProcParams.NumIncludes++;
	}

for(Idx=0;Idx<NumExcludeFiles; Idx++)
	{
	if((ProcParams.pExcludes[Idx] = OpenBedfile(ppszExcludeFiles[Idx]))==NULL)
		{
		CleanupResources(&ProcParams);
		return(eBSFerrObj);
		}
	ProcParams.NumExcludes++;
	}


if(pszBiobedFile != NULL && pszBiobedFile[0] != '\0')
	{
	if((ProcParams.pBiobed = OpenBedfile(pszBiobedFile))==NULL)
		{
		CleanupResources(&ProcParams);
		return(eBSFerrObj);
		}
	NumRegions = cNumCharRegs;	// intergenic,5'upstream,5'utr,cds,introns,3'utr,3'dnstream,3'spliceoverlap,5'spliceoverlap
	UpDnStreamLen = RegLen;
	}
else
	{
	NumRegions = 1;
	UpDnStreamLen = 0;
	ProcParams.pBiobed = NULL;
	}
NumCnts = 0;
RegionCnts = 0;


	
if(ProcMode != eProcModeNMerDistNorm)
	{
	for(Idx = 0; Idx < MaxNMerLen; Idx++)
		{
		ProcParams.CntStepOfs[Idx] = NumCnts * NumRegions;
		if(!Idx)
			RegionCnts = 4;
		else
			RegionCnts *= 4;
		NumCnts += RegionCnts;
		}
	NumCnts *= NumRegions;
	if((ProcParams.pCntStepCnts = new int[NumCnts])==NULL)
		{
		CleanupResources(&ProcParams);
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"nUnable to allocate memory (%d bytes) for holding composition statistics",sizeof(int) * NumCnts);
		return(eBSFerrMem);
		}
	memset(ProcParams.pCntStepCnts,0,NumCnts * sizeof(int));
	}
else
	{
	NumCnts = 4 << (MaxNMerLen - 1) * 2;
	if((ProcParams.pNormCnts = new tsNormDensity[NumCnts])==NULL)
		{
		CleanupResources(&ProcParams);
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"nUnable to allocate memory (%d bytes) for holding composition statistics",sizeof(tsNormDensity) * NumCnts);
		return(eBSFerrMem);
		}
	memset(ProcParams.pNormCnts,0,NumCnts * sizeof(tsNormDensity));
	ProcParams.NumNormCnts = NumCnts;
	}

ProcParams.ProcMode = ProcMode;
ProcParams.MinNMerLen = MinNMerLen;
ProcParams.MaxNMerLen = MaxNMerLen;
ProcParams.UpDnStreamLen = UpDnStreamLen;
ProcParams.NumRegions = NumRegions;
ProcParams.CurRegion = 0;
ProcParams.NumCntSteps = NumCnts;
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

#ifdef _WIN32
if((ProcParams.hRsltsFile = open(pszOutFile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
if((ProcParams.hRsltsFile = open(pszOutFile,O_RDWR | O_CREAT |O_TRUNC, S_IREAD | S_IWRITE ))==-1)
#endif
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create %s - %s",pszOutFile,strerror(errno));
    CleanupResources(&ProcParams);
	return(eBSFerrCreateFile);
	}

Rslt = GenBioseqFreqCounts(pszInFile,&ProcParams);

CleanupResources(&ProcParams);
return(Rslt);
}


int
OutputResults(tsProcParams *pProcParams)
{
char szLineBuff[2048];
int *pCnt;
int *pTmp;
int Tot;
int Len;
int NumSteps;
int StepIdx;
int Idx;
int Region;
NumSteps = 0;
for(Idx = 0; Idx < pProcParams->MaxNMerLen; Idx++)
	{
	if(Idx == 0)
		NumSteps = 1;
	NumSteps *= 4;

	if(Idx < (pProcParams->MinNMerLen -1))
		continue;

	pCnt = &pProcParams->pCntStepCnts[pProcParams->CntStepOfs[Idx]];
	for(StepIdx = 0; StepIdx < NumSteps; StepIdx++)
		{
		Tot = 0;		
		pTmp = pCnt;
		for(Region=0; Region < pProcParams->NumRegions; Region++,pCnt++)
			Tot += *pCnt;
		if(Tot == 0)
			continue;
		pCnt = pTmp;
		Len = sprintf(szLineBuff,"\"%s\",%d,%d,\"%s\",%d",
			pProcParams->ProcMode == eProcModeNMerDistAllSeqs ? pProcParams->szInFile : pProcParams->szCurChrom,
			Idx+1,StepIdx+1,StepIdx2Seq(Idx+1,StepIdx),Tot);
		if(pProcParams->NumRegions > 1)
			for(Region=0; Region < pProcParams->NumRegions; Region++,pCnt++)
				Len += sprintf(&szLineBuff[Len],",%d",*pCnt);
		else
			pCnt++;
		Len += sprintf(&szLineBuff[Len],"\n");
		CUtility::SafeWrite(pProcParams->hRsltsFile,szLineBuff,Len);
		}
	}
return(eBSFSuccess);
}

int
OutputDensityResults(tsProcParams *pProcParams)
{
tsNormDensity *pCurDensity;
int MaxLim;
int NumReported;
int Idx;
int Len;
char szLineBuff[16000];

double Deviation;
tsNormDensity *pNormCnt = pProcParams->pNormCnts;
for(int NormIdx = 0; NormIdx < pProcParams->NumNormCnts; NormIdx++,pNormCnt++)
	{
	if(pNormCnt->SumDensitiesSquared > 0.0)
		{
		Deviation = ((pNormCnt->SumDensitiesSquared - ((pNormCnt->SumDensities * pNormCnt->SumDensities) / pProcParams->NumSeqs)))/(pProcParams->NumSeqs-1);
		pNormCnt->SumDensities /= pProcParams->NumSeqs; // gives mean density 
		pNormCnt->SumDensitiesSquared = sqrt(Deviation); // gives stddev
		pNormCnt->COV = (100.0 * pNormCnt->SumDensitiesSquared) / pNormCnt->SumDensities; // gives COV
		}
	pNormCnt->Instances = NormIdx;					 // identifies which NMer
	}
	
// interested in COVs as these show the density variations
qsort(pProcParams->pNormCnts,pProcParams->NumNormCnts,sizeof(tsNormDensity),SortNMerCOV);


// report on most constant (lowest COV) 20000 first 
MaxLim = min(20000,pProcParams->NumNormCnts);
pCurDensity = pProcParams->pNormCnts;
Len = 0;
NumReported = 0;
for(Idx = 0; Idx < pProcParams->NumNormCnts; Idx++,pCurDensity++)
	{
	if(pCurDensity->SumDensitiesSquared > 0.0)
		{
		Len += sprintf(&szLineBuff[Len],"%d,\"minCOV\",%1.8f,%1.8f,%1.8f,\"%s\"\n",
			Idx+1,pCurDensity->COV,pCurDensity->SumDensitiesSquared,pCurDensity->SumDensities,StepIdx2Seq(pProcParams->MaxNMerLen,pCurDensity->Instances));
		if(Len > (sizeof(szLineBuff) * 9) / 10)
			{
			CUtility::SafeWrite(pProcParams->hRsltsFile,szLineBuff,Len);
			Len = 0;
			}
		if(++NumReported > MaxLim)
			break;
		}
	}

// report on most variable (highest COV) next
pCurDensity = &pProcParams->pNormCnts[pProcParams->NumNormCnts-MaxLim];
for(Idx = 0; Idx < MaxLim; Idx++,pCurDensity++)
	{
	if(pCurDensity->SumDensitiesSquared > 0.0)
		{
		Len += sprintf(&szLineBuff[Len],"%d,\"maxstddev\",%1.8f,%1.8f,%1.8f,\"%s\"\n",
			pProcParams->NumNormCnts-MaxLim+Idx+1,pCurDensity->COV,pCurDensity->SumDensitiesSquared,pCurDensity->SumDensities,StepIdx2Seq(pProcParams->MaxNMerLen,pCurDensity->Instances));
		if(Len > (sizeof(szLineBuff) * 9) / 10)
			{
			CUtility::SafeWrite(pProcParams->hRsltsFile,szLineBuff,Len);
			Len = 0;
			}
		}
	}

if(Len)
	CUtility::SafeWrite(pProcParams->hRsltsFile,szLineBuff,Len);
	
return(eBSFSuccess);
}


// SortNMerSumDensities
// Used to sort NMers by SumDensities (mean density over all chromosomes) 
static int 
SortNMerSumDensities( const void *arg1, const void *arg2)
{
tsNormDensity *pEl1 = (tsNormDensity *)arg1;
tsNormDensity *pEl2 = (tsNormDensity *)arg2;

if(pEl1->SumDensities < pEl2->SumDensities)
	return(-1);
if(pEl1->SumDensities > pEl2->SumDensities)
	return(1);
return(0);
}

// SortNMerSumDensitiesSquared
// Used to sort NMers by SumDensitiesSquared ( stddev)
static int 
SortNMerSumDensitiesSquared( const void *arg1, const void *arg2)
{
tsNormDensity *pEl1 = (tsNormDensity *)arg1;
tsNormDensity *pEl2 = (tsNormDensity *)arg2;

if(pEl1->SumDensitiesSquared < pEl2->SumDensitiesSquared)
	return(-1);
if(pEl1->SumDensitiesSquared > pEl2->SumDensitiesSquared)
	return(1);
return(0);
}

// SortNMerCOV
// Used to sort NMers by COV
static int 
SortNMerCOV( const void *arg1, const void *arg2)
{
tsNormDensity *pEl1 = (tsNormDensity *)arg1;
tsNormDensity *pEl2 = (tsNormDensity *)arg2;

if(pEl1->COV < pEl2->COV)
	return(-1);
if(pEl1->COV > pEl2->COV)
	return(1);
return(0);
}