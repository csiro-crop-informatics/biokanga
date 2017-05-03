/*
 * CSIRO Open Source Software License Agreement (variation of the BSD / MIT License)
 * Copyright (c) 2017, Commonwealth Scientific and Industrial Research Organisation (CSIRO) ABN 41 687 119 230.
 * See LICENCE for the complete licence information (https://github.com/csiro-crop-informatics/biokanga/LICENSE)
 * Contact: Alex Whan <alex.whan@csiro.au>
 */
#include "stdafx.h"

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#if _WIN32
#include <process.h>
#include "../libbiokanga/commhdrs.h"
#else
#include <sys/mman.h>
#include <pthread.h>
#include "../libbiokanga/commhdrs.h"
#endif

#include "biokanga.h"

#include "../libbiokanga/bgzf.h"
#include "SQLitePSL.h"
#include "Blitz.h"

int
Process(etBLZPMode PMode,				// processing mode
		char *pszExprName,				// experiment name
		char *pszExprDescr,				// experiment description
		char *pszParams,				// string containing blitz parameters
		etBLZSensitivity Sensitivity,	// sensitivity 0 - standard, 1 - high, 2 - very high, 3 - low sensitivity
		eALStrand AlignStrand,			// align on to watson, crick or both strands of target
		int MismatchScore,				// decrease score by this for each mismatch bp
		int ExactMatchScore,			// increase score by this for each exactly matching bp
		int GapOpenScore,				// decrease score by this for each gap open
		int  CoreLen,					// use this core length as the exactly matching seed length to be 5' and 3' extended
		int  CoreDelta,					// offset cores by this many bp
		int MaxExtnScoreThres,			// terminate overlap extension if curent extension score more than this; if mismatch then extension score += 2, if match and score > 0 then score -= 1 
		int MaxOccKMerDepth,			// maximum depth to explore over-occurring core K-mers
		int  MinPathScore,				// only report alignment paths on any target sequence if the path score is >= this minimum score
		int QueryLenAlignedPct,			// only report alignment paths if the percentage of total aligned bases to the query sequence length is at least this percentage (1..100)
		int  MaxPathsToReport,			// report at most this many alignment paths for any query
		etBLZRsltsFomat RsltsFormat,	// output results format
		char *pszInputFile,				// name of input file containing query sequences
		char *pszSfxFile,				// target as suffix array
		char *pszOutFile,				// where to write alignments
		int NumThreads);				// number of worker threads to use


#ifdef _WIN32
int Blitz(int argc, char* argv[])
{
// determine my process name
_splitpath(argv[0],NULL,NULL,gszProcName,NULL);
#else
int
Blitz(int argc, char** argv)
{
// determine my process name
CUtility::splitpath((char *)argv[0],NULL,gszProcName);
#endif

int iFileLogLevel;			// level of file diagnostics
int iScreenLogLevel;		// level of file diagnostics
char szLogFile[_MAX_PATH];	// write diagnostics to this file
int Rslt = 0;   			// function result code >= 0 represents success, < 0 on failure

int PMode;					// processing mode
int FMode;					// format output mode

int NumberOfProcessors;		// number of installed CPUs
int NumThreads;				// number of threads (0 defaults to number of CPUs)
int Sensitivity;			// sensitivity 0 - standard, 1 - high, 2 - very high, 3 - low sensitivity (default is 0)

int MismatchScore;			// decrease score by this for each mismatch bp
int ExactMatchScore;		// increase score by this for each exactly matching bp
int GapOpenScore;			// decrease score by this for each gap open

int CoreLen;				// use this core length as the exactly matching seed length to be 5' and 3' extended whilst no more than m_MaxSubRate
int CoreDelta;				// offset cores by this many bp
int MaxOccKMerDepth;		// maximum depth to explore over-occurring core K-mers
int MaxExtnScoreThres;		// terminate overlap extension if curent extension score more than this; if mismatch then extension score += 2, if match and score > 0 then score -= 1 

int  MinPathScore;			// only report alignment paths on any target sequence if the path score is >= this minimum score
int QueryLenAlignedPct;		// only report alignment paths if the percentage of total aligned bases to the query sequence length is at least this percentage (1..100)
int  MaxPathsToReport;		// report at most this many alignment paths for any query
int AlignStrand;			// align on to watson, crick or both strands of target
char szRsltsFile[_MAX_PATH];			// results to this file
char szTargFile[_MAX_PATH];				// align against this target suffix array genome file

char szInputFile[_MAX_PATH];		// input file containing sequences to be aligned

char szSQLiteDatabase[_MAX_PATH];	// results summaries to this SQLite file
char szExperimentName[cMaxDatasetSpeciesChrom+1];	// experiment name
char szExperimentDescr[1000];		// describes experiment
char szBlitzParams[2000];			// to hold Blitz parameters
//
struct arg_lit  *help    = arg_lit0("h","help",                 "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to screen and logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"diagnostics log file");

struct arg_int *pmode = arg_int0("m","mode","<int>",		    "alignment processing mode: 0 - standard");
struct arg_int *format = arg_int0("M","format","<int>",		    "output format: 0 - PSL, 1 - PSLX, 2 - MAF, 3 - BED, 4 - SQLite (default 0 - PSL)");
struct arg_file *inputfile = arg_file1("i","in","<file>",		"input sequences to align from this file");

struct arg_int  *alignstrand = arg_int0("Q","alignstrand","<int>", "align to this strand: 0 either, 1 Watson '+', 2 Crick '-' (default is to align to either strand)");
struct arg_file *sfxfile = arg_file1("I","sfx","<file>",		"align against this suffix array (kangax generated) file");
struct arg_file *outfile = arg_file1("o","out","<file>",		"output alignments to this file");

struct arg_int *maxocckmerdepth = arg_int0("k","maxocckmerdepth","<int>",	"maximum depth to explore over-occurring seed K-mers (default is 0 for auto, range 100 to 20000)");
struct arg_int *sensitivity = arg_int0("s","sensitivity","<int>",	"sensitivity 0 - standard, 1 - high, 2 - very high, 3 - low sensitivity (default is 0)");
struct arg_int *maxextnscorethres = arg_int0("e","extnscorethres","<int>",	"extension score threshold, core overlap extensions with extension score above this threshold are terminated;\n\t\t\t\textension score += 2 if mismatch, extension score -= 1 if match and score > 0 (default is 12)");

struct arg_int *mismatchscore = arg_int0("j","mismatchscore","<int>",	"penalise score for bp mismatches (default is 2, range 1..50)");
struct arg_int *exactmatchscore = arg_int0("J","exactmatchscore","<int>",	"score exact bp matching (default is 1, range 1..50)");
struct arg_int *gapopenscore = arg_int0("g","gapopenscore","<int>",	"penalise score for gap openings (default is 5, range 1..50)");

struct arg_int *coredelta = arg_int0("c","coredelta","<int>",	"core (seed) delta (default is 0 for auto, range 1..50)");
struct arg_int *corelen = arg_int0("C","corelen","<int>",		"core (seed) length (default is 0 for auto, range 5..50)");

struct arg_int *minpathscore = arg_int0("p","minpathscore","<int>",		"minimum alignment path score (default is 0 for auto, range 50..50000)");
struct arg_int *querylendpct = arg_int0("a","querylendpct","<int>",		"minimum required percentage of query sequence aligned (default is 25, range 1 to 100)");

struct arg_int *maxpathstoreport = arg_int0("P","maxpathstoreport","<int>",	"report at most this many highest scored alignment paths for each query (default is 10)");


struct arg_int *threads = arg_int0("T","threads","<int>",		"number of processing threads 0..128 (defaults to 0 which sets threads to number of CPU cores)");

struct arg_file *summrslts = arg_file0("q","sumrslts","<file>",		"Output results summary to this SQLite3 database file");
struct arg_str *experimentname = arg_str0("w","experimentname","<str>",		"experiment name SQLite3 database file");
struct arg_str *experimentdescr = arg_str0("W","experimentdescr","<str>",	"experiment description SQLite3 database file");

struct arg_end *end = arg_end(200);

void *argtable[] = {help,version,FileLogLevel,LogFile,
					summrslts,experimentname,experimentdescr,
					pmode,sensitivity,alignstrand,mismatchscore,exactmatchscore,gapopenscore,coredelta,corelen,maxocckmerdepth,maxextnscorethres,minpathscore,querylendpct,maxpathstoreport,format,inputfile,sfxfile,outfile,threads,
					end};

char **pAllArgs;
int argerrors;
argerrors = CUtility::arg_parsefromfile(argc,(char **)argv,&pAllArgs);
if(argerrors >= 0)
	argerrors = arg_parse(argerrors,pAllArgs,argtable);

/* special case: '--help' takes precedence over error reporting */
if (help->count > 0)
        {
		printf("\n%s %s %s, Version %s\nOptions ---\n", gszProcName,gpszSubProcess->pszName,gpszSubProcess->pszFullDescr,cpszProgVer);
        arg_print_syntax(stdout,argtable,"\n");
        arg_print_glossary(stdout,argtable,"  %-25s %s\n");
		printf("\nNote: Parameters can be entered into a parameter file, one parameter per line.");
		printf("\n      To invoke this parameter file then precede its name with '@'");
		printf("\n      e.g. %s %s @myparams.txt\n",gszProcName,gpszSubProcess->pszName);
		printf("\nPlease report any issues regarding usage of %s at https://github.com/csiro-crop-informatics/biokanga/issues\n\n",gszProcName);
		return(1);
        }

    /* special case: '--version' takes precedence error reporting */
if (version->count > 0)
        {
		printf("\n%s %s Version %s\n",gszProcName,gpszSubProcess->pszName,cpszProgVer);
		return(1);
        }

if (!argerrors)
	{
	if(FileLogLevel->count && !LogFile->count)
		{
		printf("\nError: FileLogLevel '-f%d' specified but no logfile '-F<logfile>\n'",FileLogLevel->ival[0]);
		exit(1);
		}

	iScreenLogLevel = iFileLogLevel = FileLogLevel->count ? FileLogLevel->ival[0] : eDLInfo;
	if(iFileLogLevel < eDLNone || iFileLogLevel > eDLDebug)
		{
		printf("\nError: FileLogLevel '-l%d' specified outside of range %d..%d\n",iFileLogLevel,eDLNone,eDLDebug);
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

	// now that log parameters have been parsed then initialise diagnostics log system
	if(!gDiagnostics.Open(szLogFile,(etDiagLevel)iScreenLogLevel,(etDiagLevel)iFileLogLevel,true))
		{
		printf("\nError: Unable to start diagnostics subsystem\n");
		if(szLogFile[0] != '\0')
			printf(" Most likely cause is that logfile '%s' can't be opened/created\n",szLogFile);
		exit(1);
		}

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Subprocess %s Version %s starting",gpszSubProcess->pszName,cpszProgVer);
	gExperimentID = 0;
	gProcessID = 0;
	gProcessingID = 0;
	szSQLiteDatabase[0] = '\0';
	szExperimentName[0] = '\0';
	szExperimentDescr[0] = '\0';


	if(experimentname->count)
		{
		strncpy(szExperimentName,experimentname->sval[0],sizeof(szExperimentName));
		szExperimentName[sizeof(szExperimentName)-1] = '\0';
		CUtility::TrimQuotedWhitespcExtd(szExperimentName);
		CUtility::ReduceWhitespace(szExperimentName);
		}

	if(strlen(szExperimentName) < 1)
		strcpy(szExperimentName,"N/A");

	if(experimentdescr->count)
		{
		strncpy(szExperimentDescr,experimentdescr->sval[0],sizeof(szExperimentDescr)-1);
		szExperimentDescr[sizeof(szExperimentDescr)-1] = '\0';
		CUtility::TrimQuotedWhitespcExtd(szExperimentDescr);
		CUtility::ReduceWhitespace(szExperimentDescr);
		}
	if(strlen(szExperimentDescr) < 1)
		strcpy(szExperimentDescr,"N/A");

	gExperimentID = 0;
	gProcessID = 0;
	gProcessingID = 0;
	szSQLiteDatabase[0] = '\0';

	if(summrslts->count)
		{
		strncpy(szSQLiteDatabase,summrslts->filename[0],sizeof(szSQLiteDatabase)-1);
		szSQLiteDatabase[sizeof(szSQLiteDatabase)-1] = '\0';
		CUtility::TrimQuotedWhitespcExtd(szSQLiteDatabase);
		if(strlen(szSQLiteDatabase) < 1)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: After removal of whitespace, no SQLite database specified with '-q<filespec>' option");
			return(1);
			}

		gExperimentID = gSQLiteSummaries.StartExperiment(szSQLiteDatabase,false,true,szExperimentName,szExperimentName,szExperimentDescr);
		if(gExperimentID < 1)
			return(1);
		gProcessID = gSQLiteSummaries.AddProcess((char *)gpszSubProcess->pszName,(char *)gpszSubProcess->pszName,(char *)gpszSubProcess->pszFullDescr);
		if(gProcessID < 1)
			return(1);
		gProcessingID = gSQLiteSummaries.StartProcessing(gExperimentID,gProcessID,(char *)cpszProgVer);
		if(gProcessingID < 1)
			return(1);
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Initialised SQLite database '%s' for results summary collection",szSQLiteDatabase);
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"SQLite database experiment identifier for '%s' is %d",szExperimentName,gExperimentID);
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"SQLite database process identifier for '%s' is %d",(char *)gpszSubProcess->pszName,gProcessID);
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"SQLite database processing instance identifier is %d",gProcessingID);
		}
	else
		szSQLiteDatabase[0] = '\0';

	// ensure all filenames are initialised in case not user specified
	szRsltsFile[0] = '\0';
	szTargFile[0] = '\0';
	szInputFile[0] = '\0';

	PMode = pmode->count ? pmode->ival[0] : (int)eBLZPMdefault;
	if(PMode < eBLZPMdefault || PMode >= eBLZPMplaceholder)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Processing mode '-m%d' specified outside of range %d..%d\n",PMode,eBLZPMdefault,(int)eBLZPMplaceholder-1);
		exit(1);
		}

	Sensitivity = sensitivity->count ? sensitivity->ival[0] : (int)eBLZSdefault;
	if(Sensitivity < eBLZSdefault || Sensitivity >= eBLZSplaceholder)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Sensitivity mode '-s%d' specified outside of range %d..%d\n",Sensitivity,eBLZSdefault,(int)eBLZSplaceholder-1);
		exit(1);
		}

	QueryLenAlignedPct = querylendpct->count ? querylendpct->ival[0] : cDfltMinQueryLenAlignedPct;
	if(QueryLenAlignedPct < 1 || QueryLenAlignedPct > 100)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Percentage of query sequence aligned '-a%d' specified outside of range 1..100\n",QueryLenAlignedPct);
		exit(1);
		}

	CoreLen = 0;
	CoreDelta = 0;
	MinPathScore = 0;

	AlignStrand = (eALStrand)(alignstrand->count ? alignstrand->ival[0] : eALSboth);
	if(AlignStrand < eALSboth || AlignStrand >= eALSnone)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Aligned to strand '-Q%d' specified outside of range %d..%d\n",AlignStrand,eALSboth,(int)eALSnone-1);
		exit(1);
		}


	FMode = format->count ? format->ival[0] : (int)eBLZRsltsPSL;
	if(FMode < eBLZRsltsPSL || FMode >= eBLZRsltsplaceholder)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Output format mode '-M%d' specified outside of range %d..%d\n",FMode,eBLZRsltsPSL,(int)eBLZRsltsplaceholder-1);
		exit(1);
		}

	CoreLen = corelen->count ?  corelen->ival[0] : CoreLen;
	if(CoreLen != 0 && (CoreLen < cMinCoreLen) || CoreLen > cMaxCoreLen)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: core or seed length '-s%d' specified outside of range %d..%d\n",CoreLen,cMinCoreLen,cMaxCoreLen);
		exit(1);
		}

	CoreDelta = coredelta->count ?  coredelta->ival[0] : CoreDelta;
	if((CoreDelta != 0 && CoreDelta < cMinCoreDelta) || CoreDelta > cMaxCoreDelta)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: core delta '-s%d' specified outside of range %d..%d\n",CoreDelta,cMinCoreDelta,cMaxCoreDelta);
		exit(1);
		}

	MismatchScore = mismatchscore->count ?  mismatchscore->ival[0] : cDfltMismatchScore;
	if(MismatchScore < cMinMismatchScore || MismatchScore > cMaxMismatchScore)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: mismatch penalty '-s%d' specified outside of range %d..%d\n",MismatchScore,cMinMismatchScore,cMaxMismatchScore);
		exit(1);
		}
	ExactMatchScore = exactmatchscore->count ?  exactmatchscore->ival[0] : cDfltExactMatchScore;
	if(ExactMatchScore < cMinExactMatchScore || ExactMatchScore > cMaxExactMatchScore)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: exact match score '-s%d' specified outside of range %d..%d\n",ExactMatchScore,cMinExactMatchScore,cMaxExactMatchScore);
		exit(1);
		}
	GapOpenScore = gapopenscore->count ?  gapopenscore->ival[0] : cDfltGapOpenScore;
	if(GapOpenScore < cMinGapOpenScore || GapOpenScore > cMaxGapOpenScore)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: gap open penalty '-s%d' specified outside of range %d..%d\n",GapOpenScore,cMinGapOpenScore,cMaxGapOpenScore);
		exit(1);
		}

	MaxOccKMerDepth = maxocckmerdepth->count ?  maxocckmerdepth->ival[0] : 0;
	if(MaxOccKMerDepth != 0 && (MaxOccKMerDepth < cMinOccKMerDepth) || MaxOccKMerDepth > cMaxOccKMerDepth)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: maximum depth to explore over-occurring seed K-mers '-k%d' specified outside of range %d..%d\n",MaxOccKMerDepth,cMinOccKMerDepth,cMaxOccKMerDepth);
		exit(1);
		}

	MaxExtnScoreThres = maxextnscorethres->count ?  maxextnscorethres->ival[0] : -1;	// -1 used when autodeterming the threshold
	if(MaxExtnScoreThres != -1 && (MaxExtnScoreThres < 0 || MaxExtnScoreThres > cMaxMaxExtnScoreThres))
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: maximum extension score threshold '-e%d' specified outside of range %d..%d\n",MaxExtnScoreThres,0,cMaxMaxExtnScoreThres);
		exit(1);
		}

	MinPathScore = minpathscore->count ?  minpathscore->ival[0] : MinPathScore;
	if(MinPathScore != 0 && (MinPathScore < cMinPathScore) || MinPathScore > cMaxPathScore)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: minimum path score '-p%d' specified outside of range %d..%d\n",MinPathScore,cMinPathScore,cMaxPathScore);
		exit(1);
		}

	MaxPathsToReport = maxpathstoreport->count ?  maxpathstoreport->ival[0] : cDfltMaxPathsToReport;
	if(MaxPathsToReport < 1)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: maximum number of highest scoring paths per query '-P%d' must be at least 1\n",MaxPathsToReport);
		exit(1);
		}

// show user current resource limits
#ifndef _WIN32
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Resources: %s",CUtility::ReportResourceLimits());
#endif

#ifdef _WIN32
	SYSTEM_INFO SystemInfo;
	GetSystemInfo(&SystemInfo);
	NumberOfProcessors = SystemInfo.dwNumberOfProcessors;
#else
	NumberOfProcessors = sysconf(_SC_NPROCESSORS_CONF);
#endif
	int MaxAllowedThreads = min(cMaxWorkerThreads,NumberOfProcessors);	// limit to be at most cMaxWorkerThreads
	if((NumThreads = threads->count ? threads->ival[0] : MaxAllowedThreads)==0)
		NumThreads = MaxAllowedThreads;
	if(NumThreads < 0 || NumThreads > MaxAllowedThreads)
		{
		gDiagnostics.DiagOut(eDLWarn,gszProcName,"Warning: Number of threads '-T%d' specified was outside of range %d..%d",NumThreads,1,MaxAllowedThreads);
		gDiagnostics.DiagOut(eDLWarn,gszProcName,"Warning: Defaulting number of threads to %d",MaxAllowedThreads);
		NumThreads = MaxAllowedThreads;
		}

	strcpy(szTargFile,sfxfile->filename[0]);
	strcpy(szInputFile,inputfile->filename[0]);
	strcpy(szRsltsFile,outfile->filename[0]);

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing parameters:");

	const char *pszDescr;
	switch(PMode) {
		case eBLZPMdefault:
		default:
			pszDescr = "Standard alignment processing";
			break;
		}
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Alignment processing is : '%s'",pszDescr);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Experiment name : '%s'",szExperimentName);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Experiment description : '%s'",szExperimentDescr);

	switch(Sensitivity) {
		case eBLZSdefault:
			pszDescr = "Standard alignment sensitivity";
			break;
		case eBLZSMoreSens:
			pszDescr = "High alignment sensitivity";
			break;
		case eBLZSUltraSens:
			pszDescr = "Very high alignment sensitivity - caution: very slow";
			break;
		default:
			pszDescr = "Less sensitive alignment (quicker)";
			break;
		}
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Sensitivity is : '%s'",pszDescr);
	switch(AlignStrand) {
		case eALSboth:
			pszDescr = "Watson '+' and Crick '-' strands";
			break;
		case eALSWatson:
			pszDescr = "Watson '+' strand only";
			break;
		case eALSCrick:
			pszDescr = "Crick '-' strand only";
			break;
		}

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Mismatch score penalty : %d",MismatchScore);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Exact match score : %d",ExactMatchScore);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Gap open score penalty : %d",GapOpenScore);

	if(MaxExtnScoreThres == -1)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Core extension score threshold : Auto");
	else
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Core extension score threshold : %d",MaxExtnScoreThres);
	if(CoreLen == 0)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Core length : Auto");
	else
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Core length : %d",CoreLen);
	if(CoreDelta == 0)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Core delta : Auto");
	else
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Core delta : %d",CoreDelta);
	if(MaxOccKMerDepth == 0)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Maximum depth to explore over-occurring seed K-mers : Auto");
	else
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Maximum depth to explore over-occurring seed K-mers : %d",MaxOccKMerDepth);
	if(MinPathScore == 0)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Minimum path score : Auto");
	else
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Minimum path score : %d",MinPathScore);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Minimum percentage of query sequence aligned : %d",QueryLenAlignedPct);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"maximum number of highest scoring paths per query : %d",MaxPathsToReport);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"alignments are to : %s",pszDescr);

	switch(FMode) {
		case eBLZRsltsPSL:
			pszDescr = "PSL";
			break;
		case eBLZRsltsPSLX:		
			pszDescr = "PSLX";
			break;
		case eBLZRsltsMAF:
			pszDescr = "MAF";
			break;
		case eBLZRsltsBED:
			pszDescr = "BED";
			break;
		case eBLZRsltsSQLite:
			pszDescr = "SQLite";
			break;
		}

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"output format is : '%s'",pszDescr);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"input query sequences file: '%s'",szInputFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"input target sequence(s) suffix array file: '%s'",szTargFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"output results file: '%s'",szRsltsFile);


	if(szExperimentName[0] != '\0')
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"This processing reference: %s",szExperimentName);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"number of threads : %d",NumThreads);

	if(gExperimentID > 0)
		{
		int ParamID;
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szLogFile),"log",szLogFile);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(PMode),"mode",&PMode);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(Sensitivity),"sensitivity",&Sensitivity);

		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(MismatchScore),"mismatchscore",&MismatchScore);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(ExactMatchScore),"exactmatchscore",&ExactMatchScore);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(GapOpenScore),"gapopenscore",&GapOpenScore);

		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(FMode),"format",&FMode);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(AlignStrand),"alignstrand",&AlignStrand);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(CoreLen),"corelen",&CoreLen);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(CoreDelta),"coredelta",&CoreDelta);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(MaxExtnScoreThres),"extnscorethres",&MaxExtnScoreThres);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(MaxOccKMerDepth),"maxocckmerdepth",&MaxOccKMerDepth);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(MinPathScore),"minpathscore",&MinPathScore);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(QueryLenAlignedPct),"querylendpct",&QueryLenAlignedPct);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(MaxPathsToReport),"maxpathstoreport",&MaxPathsToReport);

		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szTargFile),"sfx",szTargFile);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szRsltsFile),"out",szRsltsFile);
		
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(NumThreads),"threads",&NumThreads);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(NumberOfProcessors),"cpus",&NumberOfProcessors);

		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szSQLiteDatabase),"sumrslts",szSQLiteDatabase);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szExperimentName),"experimentname",szExperimentName);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szExperimentDescr),"experimentdescr",szExperimentDescr);
		}

	sprintf(szBlitzParams,"mode: %d sensitivity: %d mismatchscore: %d exactmatchscore: %d gapopenscore: %d alignstrand: %d corelen: %d coredelta: %d extnscorethres: %d maxocckmerdepth: %d minpathscore: %d querylendpct: %d maxpathstoreport: %d",
							PMode, Sensitivity, MismatchScore, ExactMatchScore, GapOpenScore,AlignStrand,CoreLen,CoreDelta,MaxExtnScoreThres,MaxOccKMerDepth,MinPathScore,QueryLenAlignedPct,MaxPathsToReport);

#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	gStopWatch.Start();
	Rslt = Process((etBLZPMode)PMode,szExperimentName,szExperimentDescr,szBlitzParams,(etBLZSensitivity)Sensitivity,(eALStrand)AlignStrand,MismatchScore,ExactMatchScore,GapOpenScore,CoreLen, CoreDelta,MaxExtnScoreThres,MaxOccKMerDepth,MinPathScore,QueryLenAlignedPct,MaxPathsToReport,(etBLZRsltsFomat)FMode,szInputFile,szTargFile,szRsltsFile,NumThreads);
	Rslt = Rslt >=0 ? 0 : 1;
	if(gExperimentID > 0)
		{
		if(gProcessingID)
			gSQLiteSummaries.EndProcessing(gProcessingID,Rslt);
		gSQLiteSummaries.EndExperiment(gExperimentID);
		}
	gStopWatch.Stop();

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Exit code: %d Total processing time: %s",Rslt,gStopWatch.Read());
	exit(Rslt);
	}
else
	{
    printf("\n%s %s %s, Version %s\n", gszProcName,gpszSubProcess->pszName,gpszSubProcess->pszFullDescr,cpszProgVer);
	arg_print_errors(stdout,end,gszProcName);
	arg_print_syntax(stdout,argtable,"\nUse '-h' to view option and parameter usage\n");
	exit(1);
	}
return 0;
}


int
Process(etBLZPMode PMode,				// processing mode
		char *pszExprName,				// experiment name
		char *pszExprDescr,				// experiment description
		char *pszParams,				// string containing blitz parameters
		etBLZSensitivity Sensitivity,	// sensitivity 0 - standard, 1 - high, 2 - very high, 3 - low sensitivity
		eALStrand AlignStrand,			// align on to watson, crick or both strands of target
		int MismatchScore,				// decrease score by this for each mismatch bp
		int ExactMatchScore,			// increase score by this for each exactly matching bp
		int GapOpenScore,				// decrease score by this for each gap open
		int  CoreLen,					// use this core length as the exactly matching seed length to be 5' and 3' extended whilst no more than m_MaxSubRate
		int  CoreDelta,					// offset cores by this many bp
		int MaxExtnScoreThres,			// terminate overlap extension if curent extension score more than this; if mismatch then extension score += 2, if match and score > 0 then score -= 1 
		int MaxOccKMerDepth,			// maximum depth to explore over-occurring core K-mers
		int  MinPathScore,				// only report alignment paths on any target sequence if the path score is >= this minimum score
		int QueryLenAlignedPct,				// only report alignment paths if the percentage of total aligned bases to the query sequence length is at least this percentage (1..100)
		int  MaxPathsToReport,			// report at most this many alignment paths for any query
		etBLZRsltsFomat RsltsFormat,	// output results format
		char *pszInputFile,				// name of input file containting query sequences
		char *pszSfxFile,				// target as suffix array
		char *pszOutFile,				// where to write alignments
		int NumThreads)					// number of worker threads to use
{
int Rslt;
CBlitz *pBlitzer;

if((pBlitzer = new CBlitz)==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: Unable to instantiate CAligner");
	return(eBSFerrObj);
	}
Rslt = pBlitzer->Process(PMode,pszExprName,pszExprDescr,pszParams,Sensitivity,AlignStrand,MismatchScore,ExactMatchScore,GapOpenScore,CoreLen,CoreDelta,MaxExtnScoreThres,MaxOccKMerDepth,MinPathScore,QueryLenAlignedPct,MaxPathsToReport,RsltsFormat,pszInputFile,pszSfxFile,pszOutFile,NumThreads);
delete pBlitzer;
return(Rslt);
}

CBlitz::CBlitz()
{
Init();
}

CBlitz::~CBlitz()
{
Reset(false);
}


void
CBlitz::Init(void)
{
m_hInFile = -1;
m_hOutFile = -1;

m_TotSeqIDs = 0;
m_NumQuerySeqs = 0;
m_NxtQuerySeqIdx = 0;
m_AllocdQuerySeqs = 0;
m_pQuerySeqs = NULL;
m_pSQLitePSL=NULL;
m_ProcMode = eBLZPMdefault;
m_Sensitivity = eBLZSdefault;				
m_AlignStrand = eALSboth;	
m_ExtnScoreThres = cDfltMaxExtnScoreThres;
m_CoreLen = cDfltCoreLen;
m_CoreDelta = (cDfltCoreLen+1)/2;
m_MaxOccKMerDepth = cDfltSensCoreIters;
m_QueryLenAlignedPct = cDfltMinQueryLenAlignedPct;
m_MinPathScore = cDfltPathScore;
m_MaxPathsToReport = cDfltMaxPathsToReport;
m_AlignPathID = 0;
m_RsltsFormat = eBLZRsltsPSL;	
m_pszInputFile = NULL;	
m_pszSfxFile = NULL;		
m_pszOutFile = NULL;		
m_pSfxArray = NULL;
m_pszLineBuff = NULL;
m_szLineBuffIdx = 0;
m_ReportedPaths = 0;
m_QueriesPaths = 0;
m_NumQueriesProc = 0;
m_NumThreads = 0;		
m_bMutexesCreated = false;
m_TermBackgoundThreads = 0;
}

void
CBlitz::Reset(bool bSync)			// if bSync true then fsync before closing output file handles
{
m_TermBackgoundThreads = 0x01;	// need to require any background threads to self-terminate
if(m_hInFile != -1)
	{
	close(m_hInFile);
	m_hInFile = -1;
	}
if(m_hOutFile != -1)
	{
	if(bSync)
#ifdef _WIN32
		_commit(m_hOutFile);
#else
		fsync(m_hOutFile);
#endif
	close(m_hOutFile);
	m_hOutFile = -1;
	}

if(m_pSQLitePSL != NULL)
	{
	delete m_pSQLitePSL;
	m_pSQLitePSL = NULL;
	}

if(m_pszLineBuff != NULL)
	{
	delete m_pszLineBuff;
	m_pszLineBuff = NULL;
	}

if(m_pSfxArray != NULL)
	{
	delete m_pSfxArray;
	m_pSfxArray = NULL;
	}

if(m_pQuerySeqs != NULL)
	{
	tsQuerySeq *pQuerySeq;
	if(m_NumQuerySeqs)
		{
		do {
			pQuerySeq = &m_pQuerySeqs[m_NxtQuerySeqIdx++];
			if(pQuerySeq->pQuerySeq != NULL)
				delete pQuerySeq->pQuerySeq;
			if(m_NxtQuerySeqIdx == m_AllocdQuerySeqs)
				m_NxtQuerySeqIdx = 0;
			}
		while(m_NumQuerySeqs--);
		}
	delete m_pQuerySeqs;
	m_pQuerySeqs = NULL;
	}

DeleteMutexes();

Init();
m_TermBackgoundThreads = 0x0;	// can startup any background thread processing
}



int
CBlitz::CreateMutexes(void)
{
if(m_bMutexesCreated)
	return(eBSFSuccess);

#ifdef _WIN32
InitializeSRWLock(&m_hRwLock);
#else
if(pthread_rwlock_init (&m_hRwLock,NULL)!=0)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: unable to create rwlock");
	return(eBSFerrInternal);
	}
#endif

#ifdef _WIN32
if((m_hMtxIterReads = CreateMutex(NULL,false,NULL))==NULL)
	{
#else
if(pthread_mutex_init (&m_hMtxIterReads,NULL)!=0)
	{
	pthread_rwlock_destroy(&m_hRwLock);
#endif
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: unable to create mutex");
	return(eBSFerrInternal);
	}

#ifdef _WIN32
if((m_hMtxMHReads = CreateMutex(NULL,false,NULL))==NULL)
	{
#else
if(pthread_mutex_init (&m_hMtxMHReads,NULL)!=0)
	{
#endif
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: unable to create mutex");
#ifdef _WIN32
	CloseHandle(m_hMtxIterReads);
#else
	pthread_rwlock_destroy(&m_hRwLock);
	pthread_mutex_destroy(&m_hMtxIterReads);
#endif
	return(eBSFerrInternal);
	}

m_bMutexesCreated = true;
return(eBSFSuccess);
}

void
CBlitz::DeleteMutexes(void)
{
if(!m_bMutexesCreated)
	return;
#ifdef _WIN32
CloseHandle(m_hMtxIterReads);
CloseHandle(m_hMtxMHReads);
#else
pthread_mutex_destroy(&m_hMtxIterReads);
pthread_mutex_destroy(&m_hMtxMHReads);
pthread_rwlock_destroy(&m_hRwLock);
#endif
m_bMutexesCreated = false;
}

void
CBlitz::AcquireSerialise(void)
{
#ifdef _WIN32
WaitForSingleObject(m_hMtxIterReads,INFINITE);
#else
pthread_mutex_lock(&m_hMtxIterReads);
#endif
}

void
CBlitz::ReleaseSerialise(void)
{
#ifdef _WIN32
ReleaseMutex(m_hMtxIterReads);
#else
pthread_mutex_unlock(&m_hMtxIterReads);
#endif
}

void
CBlitz::AcquireSerialiseMH(void)
{
#ifdef _WIN32
WaitForSingleObject(m_hMtxMHReads,INFINITE);
#else
pthread_mutex_lock(&m_hMtxMHReads);
#endif
}

void
CBlitz::ReleaseSerialiseMH(void)
{
#ifdef _WIN32
ReleaseMutex(m_hMtxMHReads);
#else
pthread_mutex_unlock(&m_hMtxMHReads);
#endif
}

void
CBlitz::AcquireLock(bool bExclusive)
{
#ifdef _WIN32
if(bExclusive)
	AcquireSRWLockExclusive(&m_hRwLock);
else
	AcquireSRWLockShared(&m_hRwLock);
#else
if(bExclusive)
	pthread_rwlock_wrlock(&m_hRwLock);
else
	pthread_rwlock_rdlock(&m_hRwLock);
#endif
}

void
CBlitz::ReleaseLock(bool bExclusive)
{

#ifdef _WIN32
if(bExclusive)
	ReleaseSRWLockExclusive(&m_hRwLock);
else
	ReleaseSRWLockShared(&m_hRwLock);
#else
pthread_rwlock_unlock(&m_hRwLock);
#endif
}

int
CBlitz::Process(etBLZPMode PMode,		// processing mode
		char *pszExprName,				// experiment name
		char *pszExprDescr,				// experiment description
		char *pszParams,				// string containing blitz parameters
		etBLZSensitivity Sensitivity,	// sensitivity 0 - standard, 1 - high, 2 - very high, 3 - low sensitivity
		eALStrand AlignStrand,			// align on to watson, crick or both strands of target
		int MismatchScore,				// decrease score by this for each mismatch bp
		int ExactMatchScore,			// increase score by this for each exactly matching bp
		int GapOpenScore,				// decrease score by this for each gap open
		int  CoreLen,					// use this core length (0 if determined from total target sequence length) as the exactly matching seed length to be 5' and 3' extended whilst no more than m_MaxSubRate
		int  CoreDelta,					// offset cores by this many bp
		int MaxExtnScoreThres,			// terminate overlap extension if curent extension score more than this; if mismatch then extension score += 2, if match and score > 0 then score -= 1 
		int MaxOccKMerDepth,			// maximum depth to explore over-occurring core K-mers
		int  MinPathScore,				// only report alignment paths on any target sequence if the path score is >= this minimum score
		int QueryLenAlignedPct,			// only report alignment paths if the percentage of total aligned bases to the query sequence length is at least this percentage (1..100)
		int  MaxPathsToReport,			// report at most this many alignment paths for any query
		etBLZRsltsFomat RsltsFormat,	// output results format
		char *pszInputFile,				// name of input file containing query sequences
		char *pszSfxFile,				// target as suffix array
		char *pszOutFile,				// where to write alignments
		int NumThreads)					// number of worker threads to use
{
int Rslt;
Init();

m_ProcMode = PMode;
m_Sensitivity = Sensitivity;				
m_AlignStrand = AlignStrand;	
m_MismatchScore = MismatchScore;
m_ExactMatchScore = ExactMatchScore;
m_GapOpenScore = GapOpenScore;
m_CoreLen = CoreLen;
m_CoreDelta = CoreDelta;
m_MaxOccKMerDepth = MaxOccKMerDepth;		
m_MinPathScore = MinPathScore;	
m_ExtnScoreThres = MaxExtnScoreThres;
m_QueryLenAlignedPct = QueryLenAlignedPct;
m_MaxPathsToReport = MaxPathsToReport;
m_RsltsFormat = RsltsFormat;	
m_pszInputFile = pszInputFile;	
m_pszSfxFile = pszSfxFile;		
m_pszOutFile = pszOutFile;		
m_NumThreads = NumThreads;	

if((m_pszLineBuff = new char [cAlignRprtBufferSize]) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Failed to allocate memory for alignment report buffering");
	Reset(false);
	return(eBSFerrMem);
	}

if((m_pQuerySeqs = new tsQuerySeq [cMaxReadAheadQuerySeqs]) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Failed to allocate memory for query sequences");
	Reset(false);
	return(eBSFerrMem);
	}
m_AllocdQuerySeqs = cMaxReadAheadQuerySeqs;
m_NumQuerySeqs = 0;
m_NxtQuerySeqIdx = 0;
m_TotSeqIDs = 0;

if(CreateMutexes()!=eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Failed to create thread synchronisation mutexes");
	Reset(false);
	return(cBSFSyncObjErr);
	}

m_mtqsort.SetMaxThreads(NumThreads);	
// open bioseq file containing suffix array for targeted assembly to align reads against
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Loading suffix array file '%s'", pszSfxFile);
if((m_pSfxArray = new CSfxArrayV3()) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CSfxArrayV2");
	Reset(false);
	return(eBSFerrObj);
	}
if((Rslt=m_pSfxArray->Open(pszSfxFile,false,false,false))!=eBSFSuccess)
	{
	while(m_pSfxArray->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,m_pSfxArray->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open input bioseq suffix array file '%s'",pszSfxFile);
	Reset(false);
	return(Rslt);
	}

// report to user some sfx array metadata as conformation the targeted assembly is correct
strcpy(m_szTargSpecies,m_pSfxArray->GetDatasetName());
tsSfxHeaderV3 SfxHeader;
m_pSfxArray->GetSfxHeader(&SfxHeader);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Genome Assembly Name: '%s' Descr: '%s' Title: '%s' Version: %d",
					 m_szTargSpecies,SfxHeader.szDescription,SfxHeader.szTitle,SfxHeader.Version);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Assembly has blocks: %d, max block size: %llu",SfxHeader.NumSfxBlocks,SfxHeader.SfxBlockSize);

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Loading genome assembly suffix array...");
if((Rslt=m_pSfxArray->SetTargBlock(1))<0)
	{
	while(m_pSfxArray->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,m_pSfxArray->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: unable to load genome assembly suffix array");
	return(Rslt);
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Genome assembly suffix array loaded");
// from the total sequence length then determine the core length to use
// the autodetermined core length is that length such that on average there would be expected ~ 128 to 256 copies of the core sequence in a random sequence of same length as targeted genome
UINT64 TotSeqsLen = m_pSfxArray->GetTotSeqsLen();
int AutoCoreLen = 10;
while(TotSeqsLen >>= 2)
	AutoCoreLen++;

if(CoreLen == 0)
	{
	CoreLen = AutoCoreLen;

	switch(Sensitivity) {
		case eBLZSdefault:			// default sensitivity
			break;
		case eBLZSMoreSens:			// more sensitive - slower
			CoreLen -= 2;
			break;
		case eBLZSUltraSens:		// ultra sensitive - much slower
			CoreLen -= 4;
			break;
		case eBLZSLessSens:			// less sensitive - quicker
		default:
			CoreLen += 4;
		}
	if(CoreLen > cMaxCoreLen)
		CoreLen = cMaxCoreLen;
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Using core length : %d",CoreLen);
	m_CoreLen = CoreLen;
	}
m_MinExtdCoreLen = max(30,m_CoreLen * 3);

if(CoreDelta == 0)
	{
	switch(Sensitivity) {
		case eBLZSdefault:			// default sensitivity
			CoreDelta = (CoreLen+1) / 2;
			break;
		case eBLZSMoreSens:			// more sensitive - slower
			CoreDelta = (CoreLen+2) / 3;
			break;
		case eBLZSUltraSens:		// ultra sensitive - much slower
			CoreDelta = (CoreLen+3) / 4;
			break;
		case eBLZSLessSens:			// less sensitive - quicker
		default:
			CoreDelta = CoreLen;
		}
	if(CoreDelta == 0)
		CoreDelta = 1;
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Using core delta : %d",CoreDelta);
	m_CoreDelta = CoreDelta;
	}

if(MinPathScore == 0)
	{
	MinPathScore = max(10 * CoreLen,cMinPathScore);
	if(MinPathScore > cMaxPathScore)
		MinPathScore = cMaxPathScore;	
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Using minimum path score : %d",MinPathScore);
	m_MinPathScore = MinPathScore;
	}

m_MaxIter = MaxOccKMerDepth;

// restrict the max core iterations and substitutions thresholding according to the requested sensitivity
switch(Sensitivity) {
	case eBLZSdefault:			// default sensitivity
		if(!MaxOccKMerDepth)
			m_MaxIter = cDfltSensCoreIters;
		if(MaxExtnScoreThres == -1)
			m_ExtnScoreThres = cDfltMaxExtnScoreThres;
		break;
	case eBLZSMoreSens:			// more sensitive - slower
		if(!MaxOccKMerDepth)
			m_MaxIter = cMoreSensCoreIters;
		if(MaxExtnScoreThres == -1)
			m_ExtnScoreThres = cDfltMaxExtnScoreThres + 1;
		break;
	case eBLZSUltraSens:			// ultra sensitive - much slower
		if(!MaxOccKMerDepth)
			m_MaxIter = cUltraSensCoreIters;
		if(MaxExtnScoreThres == -1)
			m_ExtnScoreThres = cDfltMaxExtnScoreThres + 2;
		break;
	case eBLZSLessSens:			// less sensitive - quicker
	default:
		if(!MaxOccKMerDepth)
			m_MaxIter = cMinSensCoreIters;
		if(MaxExtnScoreThres == -1)
			m_ExtnScoreThres = cDfltMaxExtnScoreThres - 1;
		break;
	}

if(MaxExtnScoreThres == -1)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Using overlap extension threshold score : %d",m_ExtnScoreThres);
	MaxExtnScoreThres = m_ExtnScoreThres;
	}

if(MaxOccKMerDepth == 0)
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Using maximum depth to explore over-occurring seed K-mers : %d",m_MaxIter);

m_pSfxArray->SetMaxIter(m_MaxIter);

if(CoreLen <= cMaxKmerLen)
	{
	if((Rslt = m_pSfxArray->InitOverOccKMers(CoreLen,m_MaxIter))!=eBSFSuccess)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Failed to initialise for over occurring K-mers");
		Reset(false);
		return(Rslt);
		}
	}


if(RsltsFormat == eBLZRsltsSQLite)
	{
	if((m_pSQLitePSL = new CSQLitePSL)==NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate instance of CSQLitePSL");
		return(eBSFerrObj);
		}

	if(m_pSQLitePSL->CreateDatabase(pszOutFile,false)==NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create SQLite database '%s'",pszOutFile);
		delete m_pSQLitePSL;
		m_pSQLitePSL = NULL;
		return(eBSFerrObj);
		}

	if((m_ExprID = m_pSQLitePSL->CreateExperiment(pszExprName,pszOutFile,pszInputFile,pszSfxFile,pszExprDescr,pszParams,0)) < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to initialise SQLite database '%s' with experiment details",pszOutFile);
		delete m_pSQLitePSL;
		m_pSQLitePSL = NULL;
		return(m_ExprID);
		}

	if((Rslt=m_pSQLitePSL->BeginPopulatingTables())!=eBSFSuccess)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to begin SQLite database '%s' table populating",pszOutFile);
		delete m_pSQLitePSL;
		m_pSQLitePSL = NULL;
		return(Rslt);
		}

	// add summary instances for all the target sequences
	char szSeqIdent[100];
	UINT32 SeqLen;
	int NumEntryIDs;
	int CurEntryID;
	NumEntryIDs = m_pSfxArray->GetNumEntries();
	for(CurEntryID = 1; CurEntryID <= NumEntryIDs; CurEntryID+=1)
		{
		m_pSfxArray->GetIdentName(CurEntryID,sizeof(szSeqIdent)-1,szSeqIdent);
		SeqLen = m_pSfxArray->GetSeqLen(CurEntryID);
		m_pSQLitePSL->AddAlignSummary(m_ExprID,NULL,0,szSeqIdent,SeqLen);		
		}
	}


// reads are loaded asynchronously to the alignment processing
if((Rslt=InitLoadQuerySeqs()) < eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Failed to load reads");
	Reset(false);
	return(Rslt);
	}

if(RsltsFormat != eBLZRsltsSQLite)
	{
#ifdef _WIN32
	m_hOutFile = open(pszOutFile,( O_WRONLY | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC),(_S_IREAD | _S_IWRITE));
#else
if((m_hOutFile = open(pszOutFile,O_WRONLY | O_CREAT,S_IREAD | S_IWRITE))!=-1)
	if(ftruncate(m_hOutFile,0)!=0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate %s - %s",pszOutFile,strerror(errno));
		Reset(false);
		return(eBSFerrCreateFile);
		}
#endif
	if(m_hOutFile < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate output file '%s'",pszOutFile);
		Reset(false);
		return(eBSFerrCreateFile);
		}

	// write out format specific headers
	switch(RsltsFormat) {
		case eBLZRsltsPSL:
			m_szLineBuffIdx = sprintf(m_pszLineBuff,"psLayout version 3\nGenerated by %s %s, Version %s\n",gszProcName,gpszSubProcess->pszName,cpszProgVer);
			m_szLineBuffIdx += sprintf(&m_pszLineBuff[m_szLineBuffIdx],"match	mis- 	rep. 	N's	Q gap	Q gap	T gap	T gap	strand	Q        	Q   	Q    	Q  	T        	T   	T    	T  	block	blockSizes 	qStarts	 tStarts\n");
			m_szLineBuffIdx += sprintf(&m_pszLineBuff[m_szLineBuffIdx],"     	match	match	   	count	bases	count	bases	      	name     	size	start	end	name     	size	start	end	count\n");
			m_szLineBuffIdx += sprintf(&m_pszLineBuff[m_szLineBuffIdx],"---------------------------------------------------------------------------------------------------------------------------------------------------------------\n");
			break;
		case eBLZRsltsPSLX:
			m_szLineBuffIdx = sprintf(m_pszLineBuff,"psLayout version 3\nGenerated by %s %s, Version %s\n",gszProcName,gpszSubProcess->pszName,cpszProgVer);
			m_szLineBuffIdx += sprintf(&m_pszLineBuff[m_szLineBuffIdx],"match	mis- 	rep. 	N's	Q gap	Q gap	T gap	T gap	strand	Q        	Q   	Q    	Q  	T        	T   	T    	T  	block	blockSizes 	qStarts	 tStarts\n");
			m_szLineBuffIdx += sprintf(&m_pszLineBuff[m_szLineBuffIdx],"     	match	match	   	count	bases	count	bases	      	name     	size	start	end	name     	size	start	end	count\n");
			m_szLineBuffIdx += sprintf(&m_pszLineBuff[m_szLineBuffIdx],"---------------------------------------------------------------------------------------------------------------------------------------------------------------\n");
			break;
		case eBLZRsltsMAF:
			m_szLineBuffIdx = sprintf(m_pszLineBuff,"##maf version=1 scoring=blitz");
			break;
		case eBLZRsltsBED:
			m_szLineBuffIdx = sprintf(m_pszLineBuff,"track type=bed name=\"Blitz\" description=\"Biokanga Blitz alignments\"\n");
			break;
		}
	CUtility::SafeWrite(m_hOutFile,m_pszLineBuff,m_szLineBuffIdx);
	}

m_szLineBuffIdx = 0;

InitQuerySeqThreads(NumThreads,cNumAllocdAlignNodes);	

// pickup the query sequence loader thread, if the alignment processing threads all finished then the loader thread should also have finished
#ifdef _WIN32
if(m_hThreadLoadQuerySeqs != NULL)
	{
	while(WAIT_TIMEOUT == WaitForSingleObject(m_hThreadLoadQuerySeqs, 5000))
		{
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Progress: waiting on query sequence loader thread to complete");
		}
	CloseHandle(m_hThreadLoadQuerySeqs);
	m_hThreadLoadQuerySeqs = NULL;
	}
#else
if(m_ThreadLoadQuerySeqsID != 0)
	{
	struct timespec ts;
	int JoinRlt;
	clock_gettime(CLOCK_REALTIME, &ts);
	ts.tv_sec += 5;
	while((JoinRlt = pthread_timedjoin_np(m_ThreadLoadQuerySeqsID, NULL, &ts)) != 0)
		{
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Progress: waiting on query sequence loader thread to complete");
		ts.tv_sec += 60;
		}
	}
#endif

// Checking here that the reads were all loaded w/o any major dramas!
if(m_LoadQuerySeqsRslt < 0)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Progress: problem loading");
	Reset(false);
	return(m_LoadQuerySeqsRslt);
	}

if(m_hOutFile != -1)
	{
	if(m_szLineBuffIdx > 0)
		CUtility::SafeWrite(m_hOutFile,m_pszLineBuff,m_szLineBuffIdx);
#ifdef _WIN32
	_commit(m_hOutFile);
#else
	fsync(m_hOutFile);
#endif
	close(m_hOutFile);
	m_hOutFile = -1;
	}

if(RsltsFormat == eBLZRsltsSQLite)
	{
	if((Rslt=m_pSQLitePSL->EndPopulatingTables())!=eBSFSuccess)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to complete SQLite database '%s' table populating",pszOutFile);
		delete m_pSQLitePSL;
		m_pSQLitePSL = NULL;
		return(Rslt);
		}
	delete m_pSQLitePSL;
	m_pSQLitePSL = NULL;
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"SQLite database ready for use");
	}

Reset(false);
return(Rslt);
}




#ifdef _WIN32
unsigned __stdcall AlignQuerySeqsThread(void * pThreadPars)
#else
void *AlignQuerySeqsThread(void * pThreadPars)
#endif
{
int Rslt;
tsThreadQuerySeqsPars *pPars = (tsThreadQuerySeqsPars *)pThreadPars;			// makes it easier not having to deal with casts!
CBlitz *pBlitzer = (CBlitz *)pPars->pThis;

Rslt = pBlitzer->ProcAlignQuerySeqs(pPars);
pPars->Rslt = Rslt;
#ifdef _WIN32
_endthreadex(0);
return(eBSFSuccess);
#else
pthread_exit(NULL);
#endif
}

int
CBlitz::InitQuerySeqThreads(int NumThreads,			// use this many threads
							int AlignNodes)			// each thread is allocatd this many subsequence alignment nodes
{
tsThreadQuerySeqsPars *pThreads;
tsThreadQuerySeqsPars *pThread;
m_QueriesPaths = 0;
m_ReportedPaths = 0;

if ((pThreads = new tsThreadQuerySeqsPars[NumThreads]) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "ProcessReadsetDist: Memory allocation for thread context failed");
	return(eBSFerrMem);
	}
memset(pThreads,0,sizeof(tsThreadQuerySeqsPars) * NumThreads);
int ThreadIdx;
pThread = pThreads;
for (ThreadIdx = 1; ThreadIdx <= NumThreads; ThreadIdx++, pThread++)
	{
	pThread->ThreadIdx = ThreadIdx;
	pThread->pThis = this;
	pThread->pAllocdAlignNodes = new tsQueryAlignNodes [cNumAllocdAlignNodes];
	pThread->ppFirst2Rpts = new tsQueryAlignNodes * [cNumAllocdAlignNodes];
	pThread->NumAllocdAlignNodes = cNumAllocdAlignNodes;

#ifdef _WIN32
	pThread->threadHandle = (HANDLE)_beginthreadex(NULL, 0x0fffff, AlignQuerySeqsThread, pThread, 0, &pThread->threadID);
#else
	pThread->threadRslt = pthread_create(&pThread->threadID, NULL, AlignQuerySeqsThread, pThread);
#endif
	}

// allow threads a few seconds to startup
#ifdef _WIN32
Sleep(2000);
#else
sleep(2);
#endif
UINT32 ReportedPaths;
UINT32 PrevReportedPaths = 0;
UINT32 QueriesPaths;
UINT32 PrevQueriesPaths = 0;
UINT32 NumQueriesProc;
UINT32 PrevNumQueriesProc = 0;


gDiagnostics.DiagOut(eDLInfo,gszProcName,"Progress: Generated 0 alignment paths for 0 query sequences from 0 processed");
pThread = pThreads;
for (ThreadIdx = 0; ThreadIdx < NumThreads; ThreadIdx++, pThread++)
	{
#ifdef _WIN32
	while (WAIT_TIMEOUT == WaitForSingleObject(pThread->threadHandle, 60000))
		{
		AcquireSerialise();
		ReportedPaths = m_ReportedPaths;
		QueriesPaths = m_QueriesPaths;
		NumQueriesProc = m_NumQueriesProc;
		if(ReportedPaths > PrevReportedPaths || QueriesPaths > PrevQueriesPaths || NumQueriesProc > PrevNumQueriesProc)
			{
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"Progress: Generated %u alignment paths for %u query sequences from %d processed",ReportedPaths,QueriesPaths,NumQueriesProc);
			PrevReportedPaths = ReportedPaths;
			PrevQueriesPaths = QueriesPaths;
			PrevNumQueriesProc = NumQueriesProc;
			}
		ReleaseSerialise();
		};
	CloseHandle(pThread->threadHandle);
#else
	struct timespec ts;
	int JoinRlt;
	clock_gettime(CLOCK_REALTIME, &ts);
	ts.tv_sec += 60;
	while ((JoinRlt = pthread_timedjoin_np(pThread->threadID, NULL, &ts)) != 0)
		{
		AcquireSerialise();
		ReportedPaths = m_ReportedPaths;
		QueriesPaths = m_QueriesPaths;
		NumQueriesProc = m_NumQueriesProc;
		if(ReportedPaths > PrevReportedPaths || QueriesPaths > PrevQueriesPaths || NumQueriesProc > PrevNumQueriesProc)
			{
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"Progress: Generated %u alignment paths for %u query sequences from %d processed",ReportedPaths,QueriesPaths,NumQueriesProc);
			PrevReportedPaths = ReportedPaths;
			PrevQueriesPaths = QueriesPaths;
			PrevNumQueriesProc = NumQueriesProc;
			}
		ReleaseSerialise();
		ts.tv_sec += 60;
		}
#endif
	if(pThread->pAllocdAlignNodes != NULL)
		{
		delete pThread->pAllocdAlignNodes;
		pThread->pAllocdAlignNodes = NULL; 
		}
	if(pThread->ppFirst2Rpts != NULL)
		{
		delete pThread->ppFirst2Rpts;
		pThread->ppFirst2Rpts = NULL; 
		}
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Completed reporting %u alignment paths for %u query sequences from %d processed",m_ReportedPaths,m_QueriesPaths,m_NumQueriesProc);

if(pThreads != NULL)
	delete pThreads;
return(0);
}


int
CBlitz::ProcAlignQuerySeqs(tsThreadQuerySeqsPars *pPars) 
{
int NumQueriesProc;
int PrevNumQueriesProc;
int NumMatches;
int QuerySeqLen;
int SeqID;
UINT8 *pQuerySeq;
char szQuerySeqIdent[cMaxQuerySeqIdentLen+1];

NumQueriesProc = 0;
PrevNumQueriesProc = 0;
while((pQuerySeq = DequeueQuerySeq(120,sizeof(szQuerySeqIdent),&SeqID,szQuerySeqIdent,&QuerySeqLen))!=NULL)
	{
	NumQueriesProc += 1;
	NumMatches = m_pSfxArray->LocateQuerySeqs(SeqID,pQuerySeq,QuerySeqLen,m_ExtnScoreThres,m_CoreLen,m_CoreDelta,m_AlignStrand,m_MinExtdCoreLen,pPars->NumAllocdAlignNodes,pPars->pAllocdAlignNodes,m_MaxIter);
	if(NumMatches)
		{
		if(NumMatches > 1)	// sorting by TargSeqID.QueryID.FlgStrand.TargStartOfs.QueryStartOfs
			qsort(pPars->pAllocdAlignNodes,NumMatches,sizeof(tsQueryAlignNodes),SortQueryAlignNodes);

		Report(m_MinPathScore,m_MaxPathsToReport,szQuerySeqIdent,QuerySeqLen,pQuerySeq,NumMatches,pPars->pAllocdAlignNodes,pPars->ppFirst2Rpts);

		AcquireSerialise();
		m_QueriesPaths += 1;
		m_NumQueriesProc += NumQueriesProc - PrevNumQueriesProc;
		ReleaseSerialise();
		PrevNumQueriesProc = NumQueriesProc;
		}
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Progress: Thread %d completed, processed %d query sequences",pPars->ThreadIdx,NumQueriesProc);
if(NumQueriesProc > PrevNumQueriesProc)
	{
	AcquireSerialise();
	m_NumQueriesProc += NumQueriesProc - PrevNumQueriesProc;
	ReleaseSerialise();
	}
return(0);
}

// locates and identifies the highest scoring smith-waterman path
// returns the highest score for all paths explored starting at the current node
// NOTE: recursive!!!
int											// returned best score for paths starting at pAlignNodes[ExploreNodeIdx]
CBlitz::HighScoreSW(UINT32 QueryLen,			// query length
			UINT32 TargSeqLen,					// targeted sequence length
 			bool bStrand,						// scoring for series on this strand - false if sense, true if antisense
			UINT32 ExploreNodeIdx,				// node to be explored for maximally scored path
			UINT32 NumNodes,					// total number of alignment nodes 
			tsQueryAlignNodes *pAlignNodes)	// alignment nodes
{
int BestHighScore;
UINT32 NodeIdx;
int CurNodeScore;
int GapScore;
int GapLen;
int TargGapLen;
int QueryGapLen;
int PutHighScore;
tsQueryAlignNodes *pCurNode;
tsQueryAlignNodes *pExploreNode;
pCurNode = &pAlignNodes[ExploreNodeIdx - 1];

if(pCurNode->FlgScored)
	return(pCurNode->HiScore);

// score for exactly matching bp
CurNodeScore = ((pCurNode->AlignLen - pCurNode->NumMismatches) * m_ExactMatchScore);
// penalise score for mismatches
if(pCurNode->NumMismatches)
	{
	CurNodeScore -= (pCurNode->NumMismatches * m_MismatchScore);
	if(CurNodeScore < 0)		// much easier in subsequent processing to not worry about negative scores!
		CurNodeScore = 0;
	}

pCurNode->HiScore = CurNodeScore;
pExploreNode = pAlignNodes;
BestHighScore = 0;
for(NodeIdx = 1; NodeIdx <= NumNodes; NodeIdx++,pExploreNode++)
	{
	if(ExploreNodeIdx == NodeIdx)
		continue;
	if(pExploreNode->Flg2Rpt || pExploreNode->FlgStrand != (bStrand ? 1 : 0))		// skip if already path to be reported, or if not requested strand
		continue;

	if(pExploreNode->TargStartOfs < (pCurNode->TargStartOfs + pCurNode->AlignLen - cMaxOverlapFloat))	// allowing for possible overlaps on the target sequence
		continue;
	if(pExploreNode->TargStartOfs > (pCurNode->TargStartOfs + pCurNode->AlignLen + cGapMaxLength))	// if gap too large then not on same path
		continue;
	if(pExploreNode->QueryStartOfs < (pCurNode->QueryStartOfs + pCurNode->AlignLen - cMaxOverlapFloat))   // allowing for possible overlaps on on the query sequence
		continue;
	if(pExploreNode->QueryStartOfs > (pCurNode->QueryStartOfs + cGapMaxLength))   // if gap too large then not on same path
		continue;
	QueryGapLen = abs((int)(pExploreNode->QueryStartOfs - (pCurNode->QueryStartOfs + pCurNode->AlignLen)));
	TargGapLen = abs((int)(pExploreNode->TargStartOfs - (pCurNode->TargStartOfs + pCurNode->AlignLen)));
	GapLen = (int)sqrt(((double)QueryGapLen * QueryGapLen) + ((double)TargGapLen * TargGapLen));
	GapScore = 1 + ((GapLen / 10) * cGapExtendCost);
	if(GapScore > cGapExtendCostLimit)
		GapScore = cGapExtendCostLimit;
	GapScore += m_GapOpenScore;

	if(pExploreNode->FlgScored)
		PutHighScore = pExploreNode->HiScore;
	else
		PutHighScore = HighScoreSW(QueryLen,TargSeqLen,bStrand,NodeIdx,NumNodes,pAlignNodes);
	PutHighScore += CurNodeScore;
	PutHighScore -= GapScore;
	if(PutHighScore < 0)
		PutHighScore = 0;		// much easier in subsequent processing to not worry about negative scores!
	if(PutHighScore > BestHighScore)
		{
		BestHighScore = PutHighScore;
		pCurNode->HiScore = BestHighScore; 
		pCurNode->HiScorePathNextIdx = NodeIdx;
		}
	}
pCurNode->FlgScored = 1;
return(pCurNode->HiScore);
}

// expectation is that nodes will have been sorted in TargSeqID.QueryID.FlgStrand.QueryStartOfs.TargStartOfs ascending order
// essentially is dynamic programming (a.k smith-waterman) using nodes instead of the sequences as the nodes already contain
// matching + mismatches along the diagonals
int		// returns count of paths meeting scoring threshold
CBlitz::IdentifyHighScorePaths(UINT32 QueryLen,			// query length
			UINT32 TargSeqLen,					// targeted sequence length
			bool bStrand,						// reporting for paths on this strand - false if sense, true if antisense
			UINT32 NumNodes,					// reporting on this number of nodes starting from StartNodeIdx
			UINT32 StartNodeIdx,				// report for nodes starting at this node index (1..NumNodes) which is expected to be the first alignment node of a new target sequence
			tsQueryAlignNodes *pAlignNodes,		// alignment nodes
			int MinPathScore,					// only interested in paths having at least this score
			int  MaxPathsToReport)				// report at most this many alignment paths for any query
{
int PutBestHighScore;
int BestHighScore;
int BestHighScoreNodeIdx;
UINT32 CurNodeIdx;
int NumPutPaths2Rpt;
tsQueryAlignNodes *pCurNode;
tsQueryAlignNodes *pAlignSeqNodes;

pAlignSeqNodes =  &pAlignNodes[StartNodeIdx-1];
NumPutPaths2Rpt = 0;
do {
	pCurNode = pAlignSeqNodes;
	for(CurNodeIdx = 1; CurNodeIdx <= NumNodes; CurNodeIdx++,pCurNode++)
		{
		if(pCurNode->FlgStrand != (bStrand ? 1 : 0) || pCurNode->Flg2Rpt == 1) // path scoring is strand specific, and retain scores and flags if this node already in a putative path to be reported
			continue;
		pCurNode->FlgFirst2tRpt = 0;
		pCurNode->Flg2Rpt = 0;
		pCurNode->FlgScored = 0;
		pCurNode->HiScore = 0;
		pCurNode->HiScorePathNextIdx = 0;
		}
	BestHighScore = MinPathScore - 1;
	BestHighScoreNodeIdx = 0;
	pCurNode = pAlignSeqNodes;
	for(CurNodeIdx = 1; CurNodeIdx <= NumNodes; CurNodeIdx++,pCurNode++)
		{
		if(pCurNode->Flg2Rpt || pCurNode->FlgStrand != (bStrand ? 1 : 0))	// once a node marked for reporting in a given path then can't be reported in any other path 
			continue;
		// get best score for any path originating from pCurNode, and retain if higher than any other originating nodes highest path score
		if((PutBestHighScore=HighScoreSW(QueryLen,TargSeqLen,bStrand,CurNodeIdx,NumNodes,pAlignSeqNodes)) > BestHighScore)
			{
			BestHighScore = PutBestHighScore;		// best thus far
			BestHighScoreNodeIdx = CurNodeIdx;		// record which node is the originating node
			}
		} 

	// if best path score meets the minimum required then mark 1st node as being the first and all nodes on path to be putatively reported
	if(BestHighScore >= MinPathScore)
		{
		UINT32 CurPathNodeIdx = BestHighScoreNodeIdx;
		UINT32 PathAlignedLen = 0;

		do {
			pCurNode = &pAlignSeqNodes[CurPathNodeIdx-1];
			PathAlignedLen += pCurNode->AlignLen;
			CurPathNodeIdx = pCurNode->HiScorePathNextIdx;
			}
		while(CurPathNodeIdx != 0);
		if(((PathAlignedLen * 100) / QueryLen) >= (UINT32)m_QueryLenAlignedPct)
			{
			CurPathNodeIdx = BestHighScoreNodeIdx;
			do {
				pCurNode = &pAlignSeqNodes[CurPathNodeIdx-1];
				if(CurPathNodeIdx == BestHighScoreNodeIdx)
					pCurNode->FlgFirst2tRpt=1;			// 1st node on path
				else
					pCurNode->FlgFirst2tRpt = 0;
				pCurNode->Flg2Rpt = 1;					// marking all nodes on path as being reportable
				CurPathNodeIdx = pCurNode->HiScorePathNextIdx;
				if(CurPathNodeIdx)
					pCurNode->HiScorePathNextIdx = CurPathNodeIdx + StartNodeIdx - 1;
				}
			while(CurPathNodeIdx != 0);
			NumPutPaths2Rpt += 1;
			}
		else
			BestHighScore = 0;
		}
	}
while(BestHighScore >= MinPathScore && NumPutPaths2Rpt < MaxPathsToReport);

return(NumPutPaths2Rpt);
}


int				// number of blocks processed
CBlitz::BlocksAlignStats(UINT32 *pMatches,			// returned number of bases that match that aren't repeats
					UINT32 *pmisMatches,			// returned number of bases that don't match
					UINT32 *prepMatches,			// returned number of bases that match but are part of repeats
					UINT32 *pnCount,				// returned number of 'N' bases
				char  Strand,						// query sequence strand, '+' or '-')
				UINT8 *pQuerySeq,					// the query sequence
				UINT32 qSize,						// Query sequence size
				UINT32 TargSeqID,					// CSfxArray sequence identifier
				UINT32 tSize,						// Target sequence size 
				UINT32 TargPathStartOfs,			// at this starting offset
				UINT32 TargPathEndOfs,				// ending at this offset (inclusive)
				UINT32 NumPathNodes,				// number of alignment nodes in alignment path
				int SortedPathIdx,
				tsQueryAlignNodes *pAlignNodes,		// alignment nodes
				tsQueryAlignNodes **ppFirst2Rpts) // allocated to hold ptrs to alignment nodes which are marked as being FlgFirst2tRpt
{
UINT32 Matches;
UINT32 misMatches;
UINT32 nCount;
UINT32 TotMatches;
UINT32 TotMisMatches;
UINT32 TotNCount;
UINT32 Idx;
int NumBlocks;
UINT32 MaxBlockSize;
tsQueryAlignNodes *pCurNode;
UINT8 *pTSeq;
UINT8 *pQSeq;
etSeqBase *pTargBase;
etSeqBase *pQueryBase;

if(pMatches != NULL)
	*pMatches = 0;
if(pmisMatches != NULL)
	*pmisMatches = 0;
if(prepMatches != NULL)
	*prepMatches = 0;
if(pnCount != NULL)
	*pnCount = 0;

// determine maximal sized block to allocate for
MaxBlockSize = 0;
pCurNode = ppFirst2Rpts[SortedPathIdx];
do {
	if(pCurNode->AlignLen > MaxBlockSize)
		MaxBlockSize = pCurNode->AlignLen;
	if(pCurNode->HiScorePathNextIdx > 0)
		pCurNode = &pAlignNodes[pCurNode->HiScorePathNextIdx-1];
	else
		pCurNode = NULL;
	}
while(pCurNode != NULL);

if((pTSeq = new UINT8 [MaxBlockSize+10])==NULL)
	return(eBSFerrMem);
if((pQSeq = new UINT8 [MaxBlockSize+10])==NULL)
	{
	delete pTSeq;
	return(eBSFerrMem);
	}

// iterate each block accruing mismatches and number of indeterminates
NumBlocks = 0;
TotMatches = 0;
TotMisMatches = 0;
TotNCount =  0;
pCurNode = ppFirst2Rpts[SortedPathIdx];
do {
	Matches = 0;
	misMatches = 0;
	nCount = 0;
	NumBlocks += 1;
	m_pSfxArray->GetSeq(TargSeqID,pCurNode->TargStartOfs,pTSeq,pCurNode->AlignLen);
	if(Strand == '-')
		{
		memcpy(pQSeq,&pQuerySeq[qSize - (pCurNode->QueryStartOfs + pCurNode->AlignLen)],pCurNode->AlignLen);
		CSeqTrans::ReverseComplement(pCurNode->AlignLen,pQSeq);
		}
	else
		memcpy(pQSeq,&pQuerySeq[pCurNode->QueryStartOfs],pCurNode->AlignLen);

	pTargBase = pTSeq;
	pQueryBase = pQSeq;

	for(Idx = 0; Idx < pCurNode->AlignLen; Idx++,pTargBase++,pQueryBase++)
		{
		if((*pQueryBase & 0x07) > eBaseT ||  (*pTargBase & 0x07) > eBaseT)
			{
			nCount += 1;
			misMatches += 1;
			continue;
			}
		if(*pQueryBase == *pTargBase)
			Matches += 1;
		else
			misMatches += 1;
		}

	TotMatches += Matches;
	TotMisMatches += misMatches;
	TotNCount += nCount;
	pCurNode->HiScore = pCurNode->AlignLen - misMatches;

	if(pCurNode->HiScorePathNextIdx > 0)
		pCurNode = &pAlignNodes[pCurNode->HiScorePathNextIdx-1];
	else
		pCurNode = NULL;
	}
while(pCurNode != NULL);

delete pTSeq;
delete pQSeq;
if(pMatches != NULL)
	*pMatches = TotMatches;
if(pmisMatches != NULL)
	*pmisMatches = TotMisMatches;
if(pnCount != NULL)
	*pnCount = TotNCount;
return(NumBlocks);
}


// reporting alignment as SQLite PSL format
int 
CBlitz::ReportAsSQLitePSL(UINT32 Matches,				// Number of bases that match that aren't repeats
					UINT32 misMatches,			// Number of bases that don't match
					UINT32 repMatches,			// Number of bases that match but are part of repeats
					UINT32 nCount,				// Number of 'N' bases
					UINT32	qNumInsert,			// Number of inserts in query
					UINT32 qBaseInsert,			// Number of bases inserted in query
					UINT32 tNumInsert,			// Number of inserts in target
					UINT32 tBaseInsert,			// Number of bases inserted in target
					char  Strand,				// query sequence strand, '+' or '-'
					char *pszQuerySeqIdent,     // this query sequence
					UINT32 qSize,				// Query sequence size
					UINT32 qStart,				// Alignment start position in query
					UINT32 qEnd,				// Alignment end position in query
					char *pszTargName,			// aligning to this target
					UINT32 tSize,				// Target sequence size 
					UINT32 TargPathStartOfs,	// at this starting offset
					UINT32 TargPathEndOfs,		// ending at this offset (inclusive)
					UINT32 NumPathNodes,		// number of alignment nodes in alignment path
					int SortedPathIdx,
					tsQueryAlignNodes *pAlignNodes,		// alignment nodes
					tsQueryAlignNodes **ppFirst2Rpts) // allocated to hold ptrs to alignment nodes which are marked as being FlgFirst2tRpt
{
int Score;				// Alignment score (using Blat pslScore() function)
int Identity;           // Alignment identity (using Blat 100.0 - pslCalcMilliBad(psl, TRUE) * 0.1)
int StrandQStart;
int StrandQEnd;
char szStrand[2];
UINT32 *pBlockLens;
UINT32 *pQueryBlockStarts;
UINT32 *pTargBlockStarts;
UINT32 BlockLens[5000];
UINT32 QueryBlockStarts[5000];
UINT32 TargBlockStarts[5000];

tsQueryAlignNodes *pCurNode;

StrandQStart = Strand == '+' ? qStart : qSize - (qEnd + 1),
StrandQEnd = Strand == '+' ? qEnd+1 : qSize - qStart,

szStrand[0] = Strand;
szStrand[1] = '\0';
// block sizes and starts
pCurNode = ppFirst2Rpts[SortedPathIdx];
pBlockLens = BlockLens;
pQueryBlockStarts = QueryBlockStarts;
pTargBlockStarts = TargBlockStarts;
do {
	*pBlockLens++ = pCurNode->AlignLen;
	*pQueryBlockStarts++ = pCurNode->QueryStartOfs;
	*pTargBlockStarts++ = pCurNode->TargStartOfs;
	if(pCurNode->HiScorePathNextIdx > 0)
		pCurNode = &pAlignNodes[pCurNode->HiScorePathNextIdx-1];
	else
		pCurNode = NULL;
	}
while(pCurNode != NULL);

Score = m_pSQLitePSL->pslScore(Matches,misMatches,repMatches,qNumInsert,tNumInsert,szStrand,tSize,TargPathStartOfs,TargPathEndOfs+1,NumPathNodes,(int *)BlockLens,(int *)QueryBlockStarts,(int *)TargBlockStarts);
double pslIdent = (double)m_pSQLitePSL->pslCalcMilliBad(Matches,misMatches,repMatches,qNumInsert,tNumInsert,qSize,StrandQStart,StrandQEnd,szStrand,tSize,TargPathStartOfs,TargPathEndOfs,NumPathNodes,(int *)BlockLens,(int *)QueryBlockStarts,(int *)TargBlockStarts,true); 
Identity = (int)(100.0 - pslIdent * 0.1);
AcquireLock(true);

m_pSQLitePSL->AddAlignment(m_ExprID,Score,Identity,Matches,misMatches,repMatches,nCount,qNumInsert,qBaseInsert,tNumInsert,tBaseInsert,szStrand,pszQuerySeqIdent,qSize,StrandQStart,StrandQEnd,pszTargName,tSize,TargPathStartOfs,TargPathEndOfs+1,NumPathNodes,(int *)BlockLens,(int *)QueryBlockStarts,(int *)TargBlockStarts);

m_ReportedPaths += 1;				
ReleaseLock(true);
return(NumPathNodes);
}

// reporting alignment as PSL format
int 
CBlitz::ReportAsPSL(UINT32 Matches,				// Number of bases that match that aren't repeats
					UINT32 misMatches,			// Number of bases that don't match
					UINT32 repMatches,			// Number of bases that match but are part of repeats
					UINT32 nCount,				// Number of 'N' bases
					UINT32	qNumInsert,			// Number of inserts in query
					UINT32 qBaseInsert,			// Number of bases inserted in query
					UINT32 tNumInsert,			// Number of inserts in target
					UINT32 tBaseInsert,			// Number of bases inserted in target
					char  Strand,				// query sequence strand, '+' or '-'
					char *pszQuerySeqIdent,     // this query sequence
					UINT32 qSize,				// Query sequence size
					UINT32 qStart,				// Alignment start position in query
					UINT32 qEnd,				// Alignment end position in query
					char *pszTargName,			// aligning to this target
					UINT32 tSize,				// Target sequence size 
					UINT32 TargPathStartOfs,	// at this starting offset
					UINT32 TargPathEndOfs,		// ending at this offset (inclusive)
					UINT32 NumPathNodes,		// number of alignment nodes in alignment path
					int SortedPathIdx,
					tsQueryAlignNodes *pAlignNodes,		// alignment nodes
					tsQueryAlignNodes **ppFirst2Rpts) // allocated to hold ptrs to alignment nodes which are marked as being FlgFirst2tRpt
{
tsQueryAlignNodes *pCurNode;

AcquireSerialise();
if(m_szLineBuffIdx > (cAlignRprtBufferSize * 9) / 10)
	{
	CUtility::SafeWrite(m_hOutFile,m_pszLineBuff,m_szLineBuffIdx);
	m_szLineBuffIdx = 0;
	}

m_szLineBuffIdx += sprintf(&m_pszLineBuff[m_szLineBuffIdx],"%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%c\t%s\t%u\t%u\t%u\t%s\t%u\t%u\t%u\t%u\t",
					Matches,misMatches,repMatches,nCount,qNumInsert,qBaseInsert,tNumInsert,tBaseInsert,Strand,
					pszQuerySeqIdent,qSize,
					Strand == '+' ? qStart : qSize - (qEnd + 1),
					Strand == '+' ? qEnd+1 : qSize - qStart,
					pszTargName,tSize,TargPathStartOfs,TargPathEndOfs+1,NumPathNodes);

// block sizes
pCurNode = ppFirst2Rpts[SortedPathIdx];
do {
	m_szLineBuffIdx += sprintf(&m_pszLineBuff[m_szLineBuffIdx],"%u,",pCurNode->AlignLen);
	if(pCurNode->HiScorePathNextIdx > 0)
		pCurNode = &pAlignNodes[pCurNode->HiScorePathNextIdx-1];
	else
		pCurNode = NULL;
	}
while(pCurNode != NULL);
m_pszLineBuff[m_szLineBuffIdx++] = '\t';

// query starts
pCurNode = ppFirst2Rpts[SortedPathIdx];
do {
	m_szLineBuffIdx += sprintf(&m_pszLineBuff[m_szLineBuffIdx],"%u,", pCurNode->QueryStartOfs);
	if(pCurNode->HiScorePathNextIdx > 0)
		pCurNode = &pAlignNodes[pCurNode->HiScorePathNextIdx-1];
	else
		pCurNode = NULL;
	}
while(pCurNode != NULL);
m_pszLineBuff[m_szLineBuffIdx++] = '\t';

// target starts
pCurNode = ppFirst2Rpts[SortedPathIdx];
do {
	m_szLineBuffIdx += sprintf(&m_pszLineBuff[m_szLineBuffIdx],"%u,",pCurNode->TargStartOfs);
	if(pCurNode->HiScorePathNextIdx > 0)
		pCurNode = &pAlignNodes[pCurNode->HiScorePathNextIdx-1];
	else
		pCurNode = NULL;
	}
while(pCurNode != NULL);
m_szLineBuffIdx += sprintf(&m_pszLineBuff[m_szLineBuffIdx],"\n");	
if(m_szLineBuffIdx > (cAlignRprtBufferSize * 9) / 10)
	{
	CUtility::SafeWrite(m_hOutFile,m_pszLineBuff,m_szLineBuffIdx);
	m_szLineBuffIdx = 0;
	}
m_ReportedPaths += 1;				
ReleaseSerialise();
return(NumPathNodes);
}

// reporting alignment as PSLX format
int 
CBlitz::ReportAsPSLX(UINT32 Matches,				// Number of bases that match that aren't repeats
					UINT32 misMatches,			// Number of bases that don't match
					UINT32 repMatches,			// Number of bases that match but are part of repeats
					UINT32 nCount,				// Number of 'N' bases
					UINT32	qNumInsert,			// Number of inserts in query
					UINT32 qBaseInsert,			// Number of bases inserted in query
					UINT32 tNumInsert,			// Number of inserts in target
					UINT32 tBaseInsert,			// Number of bases inserted in target
					char  Strand,				// query sequence strand, '+' or '-'
					char *pszQuerySeqIdent,     // this query sequence
					UINT32 qSize,				// Query sequence size
					UINT8 *pQuerySeq,			// the query sequence
					UINT32 qStart,				// Alignment start position in query
					UINT32 qEnd,				// Alignment end position in query
					UINT32 TargSeqID,			// CSfxArray sequence identifier
					char *pszTargName,			// aligning to this target
					UINT32 tSize,				// Target sequence size 
					UINT32 TargPathStartOfs,	// at this starting offset
					UINT32 TargPathEndOfs,		// ending at this offset (inclusive)
					UINT32 NumPathNodes,		// number of alignment nodes in alignment path
					int SortedPathIdx,
					tsQueryAlignNodes *pAlignNodes,		// alignment nodes
					tsQueryAlignNodes **ppFirst2Rpts) // allocated to hold ptrs to alignment nodes which are marked as being FlgFirst2tRpt
{
UINT32 MaxBlockSize;
etSeqBase *pCurSeq;
tsQueryAlignNodes *pCurNode;

pCurSeq = NULL;
AcquireSerialise();
if(m_szLineBuffIdx > (cAlignRprtBufferSize * 9) / 10)
	{
	CUtility::SafeWrite(m_hOutFile,m_pszLineBuff,m_szLineBuffIdx);
	m_szLineBuffIdx = 0;
	}

m_szLineBuffIdx += sprintf(&m_pszLineBuff[m_szLineBuffIdx],"%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%c\t%s\t%u\t%u\t%u\t%s\t%u\t%u\t%u\t%u\t",
					Matches,misMatches,repMatches,nCount,qNumInsert,qBaseInsert,tNumInsert,tBaseInsert,Strand,
					pszQuerySeqIdent,qSize,
					Strand == '+' ? qStart : qSize - (qEnd + 1),
					Strand == '+' ? qEnd+1 : qSize - qStart,
					pszTargName,tSize,TargPathStartOfs,TargPathEndOfs+1,NumPathNodes);

// block sizes
MaxBlockSize = 0;
pCurNode = ppFirst2Rpts[SortedPathIdx];
do {
	m_szLineBuffIdx += sprintf(&m_pszLineBuff[m_szLineBuffIdx],"%u,",pCurNode->AlignLen);
	if(pCurNode->AlignLen > MaxBlockSize)
		MaxBlockSize = pCurNode->AlignLen;
	if(pCurNode->HiScorePathNextIdx > 0)
		pCurNode = &pAlignNodes[pCurNode->HiScorePathNextIdx-1];
	else
		pCurNode = NULL;
	}
while(pCurNode != NULL);
m_pszLineBuff[m_szLineBuffIdx++] = '\t';

// query starts
pCurNode = ppFirst2Rpts[SortedPathIdx];
do {
	m_szLineBuffIdx += sprintf(&m_pszLineBuff[m_szLineBuffIdx],"%u,", pCurNode->QueryStartOfs);
	if(pCurNode->HiScorePathNextIdx > 0)
		pCurNode = &pAlignNodes[pCurNode->HiScorePathNextIdx-1];
	else
		pCurNode = NULL;
	}
while(pCurNode != NULL);
m_pszLineBuff[m_szLineBuffIdx++] = '\t';

// target starts
pCurNode = ppFirst2Rpts[SortedPathIdx];
do {
	m_szLineBuffIdx += sprintf(&m_pszLineBuff[m_szLineBuffIdx],"%u,",pCurNode->TargStartOfs);
	if(pCurNode->HiScorePathNextIdx > 0)
		pCurNode = &pAlignNodes[pCurNode->HiScorePathNextIdx-1];
	else
		pCurNode = NULL;
	}
while(pCurNode != NULL);

// allocate to hold the maximally sized sequence block
if((pCurSeq = new UINT8 [10 + TargPathEndOfs - TargPathStartOfs])==NULL)	// inplace translation to ascii so allow a few extra for terminagting '\0'	
	{
	ReleaseSerialise();
	return(eBSFerrMem);
	}

// get query sequences for each block and report these
// query starts
pCurNode = ppFirst2Rpts[SortedPathIdx];
do {
	if((m_szLineBuffIdx + pCurNode->AlignLen) > (cAlignRprtBufferSize * 9) / 10)
		{
		CUtility::SafeWrite(m_hOutFile,m_pszLineBuff,m_szLineBuffIdx);
		m_szLineBuffIdx = 0;
		}
	if(Strand == '-')
		{
		memcpy(pCurSeq,&pQuerySeq[qSize - (pCurNode->QueryStartOfs + pCurNode->AlignLen)],pCurNode->AlignLen);
		CSeqTrans::ReverseComplement(pCurNode->AlignLen,pCurSeq);
		}
	else
		memcpy(pCurSeq,&pQuerySeq[pCurNode->QueryStartOfs],pCurNode->AlignLen);
	m_szLineBuffIdx += sprintf(&m_pszLineBuff[m_szLineBuffIdx],"%s,",CSeqTrans::MapSeq2Ascii(pCurSeq,pCurNode->AlignLen,(char *)pCurSeq));
	if(pCurNode->HiScorePathNextIdx > 0)
		pCurNode = &pAlignNodes[pCurNode->HiScorePathNextIdx-1];
	else
		pCurNode = NULL;
	}
while(pCurNode != NULL);
m_pszLineBuff[m_szLineBuffIdx++] = '\t';


// get target sequences for each block and report these
pCurNode = ppFirst2Rpts[SortedPathIdx];
do {
	if((m_szLineBuffIdx + pCurNode->AlignLen) > (cAlignRprtBufferSize * 9) / 10)
		{
		CUtility::SafeWrite(m_hOutFile,m_pszLineBuff,m_szLineBuffIdx);
		m_szLineBuffIdx = 0;
		}
	m_pSfxArray->GetSeq(TargSeqID,pCurNode->TargStartOfs,pCurSeq,pCurNode->AlignLen);
	m_szLineBuffIdx += sprintf(&m_pszLineBuff[m_szLineBuffIdx],"%s,",CSeqTrans::MapSeq2Ascii(pCurSeq,pCurNode->AlignLen,(char *)pCurSeq));
	if(pCurNode->HiScorePathNextIdx > 0)
		pCurNode = &pAlignNodes[pCurNode->HiScorePathNextIdx-1];
	else
		pCurNode = NULL;
	}
while(pCurNode != NULL);

m_szLineBuffIdx += sprintf(&m_pszLineBuff[m_szLineBuffIdx],"\n");	
if(m_szLineBuffIdx > (cAlignRprtBufferSize * 9) / 10)
	{
	CUtility::SafeWrite(m_hOutFile,m_pszLineBuff,m_szLineBuffIdx);
	m_szLineBuffIdx = 0;
	}
m_ReportedPaths += 1;				
ReleaseSerialise();
if(pCurSeq!=NULL)
	delete pCurSeq;
return(NumPathNodes);
}

// reporting alignment as MAF format
int 
CBlitz::ReportAsMAF(int PathScore,				// score for this path
					UINT32 Matches,				// Number of bases that match that aren't repeats
					UINT32 misMatches,			// Number of bases that don't match
					UINT32 repMatches,			// Number of bases that match but are part of repeats
					UINT32 nCount,				// Number of 'N' bases
					UINT32	qNumInsert,			// Number of inserts in query
					UINT32 qBaseInsert,			// Number of bases inserted in query
					UINT32 tNumInsert,			// Number of inserts in target
					UINT32 tBaseInsert,			// Number of bases inserted in target
					char  Strand,				// query sequence strand, '+' or '-'
					char *pszQuerySeqIdent,     // this query sequence
					UINT32 qSize,				// Query sequence size
					UINT8 *pQuerySeq,			// the query sequence
					UINT32 qStart,				// Alignment start position in query
					UINT32 qEnd,				// Alignment end position in query
					UINT32 TargSeqID,			// CSfxArray sequence identifier
					char *pszTargName,			// aligning to this target
					UINT32 tSize,				// Target sequence size 
					UINT32 TargPathStartOfs,	// at this starting offset
					UINT32 TargPathEndOfs,		// ending at this offset (inclusive)
					UINT32 NumPathNodes,		// number of alignment nodes in alignment path
					int SortedPathIdx,
					tsQueryAlignNodes *pAlignNodes,		// alignment nodes
					tsQueryAlignNodes **ppFirst2Rpts) // allocated to hold ptrs to alignment nodes which are marked as being FlgFirst2tRpt
{
UINT32 MaxBlockSize;
tsQueryAlignNodes *pCurNode;
UINT8 *pCurSeq;

// // allocate buffering for maximal sized alignment block
MaxBlockSize = 0;
pCurNode = ppFirst2Rpts[SortedPathIdx];
do {
	if(pCurNode->AlignLen > MaxBlockSize)
		MaxBlockSize = pCurNode->AlignLen;
	if(pCurNode->HiScorePathNextIdx > 0)
		pCurNode = &pAlignNodes[pCurNode->HiScorePathNextIdx-1];
	else
		pCurNode = NULL;
	}
while(pCurNode != NULL);

if((pCurSeq = new UINT8 [10 + MaxBlockSize])==NULL)	// inplace translation to ascii so allow a few extra for terminagting '\0'	
	return(eBSFerrMem);

AcquireSerialise();
if(m_szLineBuffIdx > (cAlignRprtBufferSize * 9) / 10)
	{
	CUtility::SafeWrite(m_hOutFile,m_pszLineBuff,m_szLineBuffIdx);
	m_szLineBuffIdx = 0;
	}

pCurNode = ppFirst2Rpts[SortedPathIdx];
do {
	if((m_szLineBuffIdx + pCurNode->AlignLen) > (cAlignRprtBufferSize * 9) / 10)
		{
		CUtility::SafeWrite(m_hOutFile,m_pszLineBuff,m_szLineBuffIdx);
		m_szLineBuffIdx = 0;
		}	
	m_szLineBuffIdx += sprintf(&m_pszLineBuff[m_szLineBuffIdx],"\na\tscore=%d",pCurNode->HiScore);
	m_szLineBuffIdx += sprintf(&m_pszLineBuff[m_szLineBuffIdx],"\ns\t%s\t%u\t%u\t%c\t%u\t",pszTargName,pCurNode->TargStartOfs,pCurNode->AlignLen,'+',tSize);
	m_pSfxArray->GetSeq(TargSeqID,pCurNode->TargStartOfs,pCurSeq,pCurNode->AlignLen);
	m_szLineBuffIdx += sprintf(&m_pszLineBuff[m_szLineBuffIdx],"%s",CSeqTrans::MapSeq2Ascii(pCurSeq,pCurNode->AlignLen,(char *)pCurSeq));

	m_szLineBuffIdx += sprintf(&m_pszLineBuff[m_szLineBuffIdx],"\ns\t%s\t%u\t%u\t%c\t%u\t",pszQuerySeqIdent,pCurNode->TargStartOfs,pCurNode->AlignLen,'+',qSize);
	if(Strand == '-')
		{
		memcpy(pCurSeq,&pQuerySeq[qSize - (pCurNode->QueryStartOfs + pCurNode->AlignLen)],pCurNode->AlignLen);
		CSeqTrans::ReverseComplement(pCurNode->AlignLen,pCurSeq);
		}
	else
		memcpy(pCurSeq,&pQuerySeq[pCurNode->QueryStartOfs],pCurNode->AlignLen);

	m_szLineBuffIdx += sprintf(&m_pszLineBuff[m_szLineBuffIdx],"%s\n",CSeqTrans::MapSeq2Ascii(pCurSeq,pCurNode->AlignLen,(char *)pCurSeq));
	if(pCurNode->HiScorePathNextIdx > 0)
		pCurNode = &pAlignNodes[pCurNode->HiScorePathNextIdx-1];
	else
		pCurNode = NULL;
	}
while(pCurNode != NULL);
m_ReportedPaths += 1;
ReleaseSerialise();
if(pCurSeq != NULL)
	delete pCurSeq;

return(NumPathNodes);
}

// reporting alignment as BED format
int 
CBlitz::ReportAsBED(char *pszQuerySeqIdent,     // this query sequence
					char  Strand,				// query sequence strand, '+' or '-'
					UINT32 AlignScore,			// alignment has this score
					char *pszTargName,			// aligning to this target
					UINT32 NumPathNodes,		// number of alignment nodes in alignment path
					UINT32 TargPathStartOfs,	// at this starting offset
					UINT32 TargPathEndOfs,		// ending at this offset (inclusive)
					int SortedPathIdx,
					tsQueryAlignNodes *pAlignNodes,		// alignment nodes
					tsQueryAlignNodes **ppFirst2Rpts) // allocated to hold ptrs to alignment nodes which are marked as being FlgFirst2tRpt
{
tsQueryAlignNodes *pCurNode;
AlignScore = (UINT32)(10 * sqrt(AlignScore));
if(AlignScore > 1000)
	AlignScore = 1000;
AcquireSerialise();
if(m_szLineBuffIdx > (cAlignRprtBufferSize * 9) / 10)
	{
	CUtility::SafeWrite(m_hOutFile,m_pszLineBuff,m_szLineBuffIdx);
	m_szLineBuffIdx = 0;
	}

m_szLineBuffIdx += sprintf(&m_pszLineBuff[m_szLineBuffIdx],"%s\t%d\t%d\t%s\t%d\t%c\t%d\t%d\t0\t%d\t",
					pszTargName,TargPathStartOfs,TargPathEndOfs+1,pszQuerySeqIdent,AlignScore,Strand,TargPathStartOfs,TargPathEndOfs+1,NumPathNodes);
pCurNode = ppFirst2Rpts[SortedPathIdx];
do {
	m_szLineBuffIdx += sprintf(&m_pszLineBuff[m_szLineBuffIdx],"%u,",pCurNode->AlignLen);
	if(pCurNode->HiScorePathNextIdx > 0)
		pCurNode = &pAlignNodes[pCurNode->HiScorePathNextIdx-1];
	else
		pCurNode = NULL;
	}
while(pCurNode != NULL);
m_pszLineBuff[m_szLineBuffIdx++] = '\t';
pCurNode = ppFirst2Rpts[SortedPathIdx];
do {
	m_szLineBuffIdx += sprintf(&m_pszLineBuff[m_szLineBuffIdx],"%u",pCurNode->TargStartOfs - TargPathStartOfs);
	if(pCurNode->HiScorePathNextIdx > 0)
		{
		m_pszLineBuff[m_szLineBuffIdx++] = ',';
		pCurNode = &pAlignNodes[pCurNode->HiScorePathNextIdx-1];
		}
	else
		pCurNode = NULL;
	}
while(pCurNode != NULL);
m_szLineBuffIdx += sprintf(&m_pszLineBuff[m_szLineBuffIdx],"\n");
if(m_szLineBuffIdx > (cAlignRprtBufferSize * 9) / 10)
	{
	CUtility::SafeWrite(m_hOutFile,m_pszLineBuff,m_szLineBuffIdx);
	m_szLineBuffIdx = 0;
	}
m_ReportedPaths += 1;
ReleaseSerialise();
return(NumPathNodes);
}


int											// returns number of reported paths
CBlitz::Report(int MinPathScore,			// only report paths having at least this minimum score
				int  MaxPathsToReport,		// report at most this many alignment paths for any query
				char *pszQuerySeqIdent,		// query sequence 
				UINT32 QueryLen,			// query length
				UINT8 *pQuerySeq,			// the query sequence
				UINT32 NumNodes,			// number of alignment nodes
				tsQueryAlignNodes *pAlignNodes,		// alignment nodes
				tsQueryAlignNodes **ppFirst2Rpts)	// allocated to hold ptrs to alignment nodes which are marked as being FlgFirst2tRpt
{
tsQueryAlignNodes *pCurNode;
UINT32 MaxAlignLen;
UINT32 CurTargSeqID;
UINT32 CurTargMatchLen;
UINT32 MaxTargMatchLen;
UINT32 CurTargMatchNodes;
UINT32 MaxTargMatchNodes;
UINT32 MaxNodeAlignLen;
UINT32 MaxNodeAlignLenMMs;

UINT32 NodeIdx;
UINT32 CurPathNodeIdx;
int SortedPathIdx;
UINT32 QueryPathStartOfs;
UINT32 TargPathStartOfs;
UINT32 QueryPathEndOfs;
UINT32 TargPathEndOfs;
UINT32 TotalMMs;
UINT32 NumPathNodes;
UINT32	qNumInsert;			// Number of inserts in query
UINT32 qBaseInsert;			// Number of bases inserted in query
UINT32 tNumInsert;			// Number of inserts in target
UINT32 tBaseInsert;			// Number of bases inserted in target
int TargGap;
int QueryGap;

UINT32 TotalAlignLen;
bool bStrand;
int PathHiScore;
UINT32 TargSeqLen;
char szTargName[100];

pCurNode = pAlignNodes;
MaxAlignLen = 0;
CurTargSeqID = 0;
CurTargMatchLen = 0;
MaxTargMatchLen = 0;
CurTargMatchNodes = 0;
MaxTargMatchNodes = 0;
MaxNodeAlignLen = 0;
MaxNodeAlignLenMMs = 0;

// nodes will have been sorted in TargSeqID.QueryID.FlgStrand.QueryStartOfs.TargStartOfs ascending order
UINT32 StartTargNodeIdx;
UINT32 EndTargNodeIdx;
UINT32 NumPutPaths;

NumPutPaths = 0;
CurTargSeqID = pCurNode->TargSeqID;
CurTargMatchNodes = 0;
for(StartTargNodeIdx=EndTargNodeIdx=1; EndTargNodeIdx <= NumNodes; EndTargNodeIdx++, pCurNode++)
	{
	if(pCurNode->TargSeqID != CurTargSeqID || EndTargNodeIdx == NumNodes)				// onto a new target sequence or last node?
		{
		if(EndTargNodeIdx == NumNodes)
			CurTargMatchNodes += 1;
		TargSeqLen = m_pSfxArray->GetSeqLen(CurTargSeqID);
		if(m_AlignStrand != eALSCrick)
 			NumPutPaths += IdentifyHighScorePaths(QueryLen,TargSeqLen,false,CurTargMatchNodes,StartTargNodeIdx,pAlignNodes,MinPathScore,MaxPathsToReport);	// sense/sense paths
		if(m_AlignStrand != eALSWatson)
			NumPutPaths += IdentifyHighScorePaths(QueryLen,TargSeqLen,true,CurTargMatchNodes,StartTargNodeIdx,pAlignNodes,MinPathScore,MaxPathsToReport);     // antisense/sense paths
		CurTargSeqID = pCurNode->TargSeqID;
		CurTargMatchNodes = 0;
		StartTargNodeIdx = EndTargNodeIdx;
		}
	CurTargMatchNodes += 1;
	}
if(NumPutPaths == 0)
	return(0);

UINT32 NumHeadNodes = 0;
pCurNode = pAlignNodes;
for(NodeIdx=1; NodeIdx <= NumNodes; NodeIdx++, pCurNode++)
	{
	if(pCurNode->FlgFirst2tRpt)
		ppFirst2Rpts[NumHeadNodes++] = pCurNode;
	}

if(NumHeadNodes > 1)			// sort by highest scoring path descending
	qsort(ppFirst2Rpts,NumHeadNodes,sizeof(tsQueryAlignNodes *),SortHighScoreDescend);

for(SortedPathIdx=0; SortedPathIdx < min((int)NumHeadNodes,MaxPathsToReport); SortedPathIdx++)
	{
	pCurNode = ppFirst2Rpts[SortedPathIdx];
	bStrand = pCurNode->FlgStrand ? true : false;
	PathHiScore = pCurNode->HiScore;
	m_pSfxArray->GetIdentName(pCurNode->TargSeqID,sizeof(szTargName),szTargName);
	TargSeqLen = m_pSfxArray->GetSeqLen(pCurNode->TargSeqID);
	QueryPathStartOfs = pCurNode->QueryStartOfs;
	TargPathStartOfs = pCurNode->TargStartOfs;
	TotalMMs = pCurNode->NumMismatches;
	TotalAlignLen = pCurNode->AlignLen;
	QueryPathEndOfs = pCurNode->QueryStartOfs + pCurNode->AlignLen - 1;
	TargPathEndOfs = pCurNode->TargStartOfs + pCurNode->AlignLen - 1;
	NumPathNodes = 1;
	qNumInsert = 0;
	qBaseInsert = 0;
	tNumInsert = 0;
	tBaseInsert = 0;

	tsQueryAlignNodes *pTT = pCurNode;
	UINT32 QueryPathEndOfsTT = QueryPathEndOfs;
	UINT32 TargPathEndOfsTT = TargPathEndOfs;
	int DeltaGap;
	while((CurPathNodeIdx = pTT->HiScorePathNextIdx) != 0)
		{
		pTT = &pAlignNodes[CurPathNodeIdx-1];
		QueryGap = pTT->QueryStartOfs - QueryPathEndOfsTT - 1;
		TargGap = pTT->TargStartOfs - TargPathEndOfsTT - 1;	
		if(QueryGap < 0 || TargGap < 0)
			{
			if(QueryGap >= 0)	// if query gap is >= 0 then targ gap must be < 0
				DeltaGap = abs(TargGap);
			else
				{
				if(TargGap >= 0)	// if targ gap >= 0 then query gap must be < 0
					DeltaGap = abs(QueryGap);
				else    // else both query and targ gap must have been < 0
					{
					DeltaGap = abs(min(QueryGap,TargGap));
					}
				}
			pTT->QueryStartOfs += DeltaGap;
			pTT->TargStartOfs += DeltaGap;
			pTT->AlignLen -= DeltaGap;
			}
		
		QueryPathEndOfsTT = pTT->QueryStartOfs + pTT->AlignLen - 1;
		TargPathEndOfsTT = pTT->TargStartOfs + pTT->AlignLen - 1;
		}

	while((CurPathNodeIdx = pCurNode->HiScorePathNextIdx) != 0)
		{
		pCurNode = &pAlignNodes[CurPathNodeIdx-1];
		QueryGap = pCurNode->QueryStartOfs - QueryPathEndOfs - 1;
		TargGap = pCurNode->TargStartOfs - TargPathEndOfs - 1;	
		TotalMMs += pCurNode->NumMismatches;
		TotalAlignLen += pCurNode->AlignLen;
		QueryPathEndOfs = pCurNode->QueryStartOfs + pCurNode->AlignLen - 1;
		TargPathEndOfs = pCurNode->TargStartOfs + pCurNode->AlignLen - 1;
		if(QueryGap > 0)
			{
			qNumInsert += 1;
			qBaseInsert += QueryGap;
			}
		if(TargGap > 0)
			{
			tNumInsert += 1;
			tBaseInsert += TargGap;		
			}
		NumPathNodes += 1;
		}

	// recalc total matches, mismatches, indeterminates as the block lengths and starting offsets may have been modified
	UINT32 Matches,MisMatches,NumbNs;
	BlocksAlignStats(&Matches,&MisMatches,NULL,&NumbNs,bStrand == true ? '-' : '+',pQuerySeq,QueryLen,pCurNode->TargSeqID,TargSeqLen,TargPathStartOfs,TargPathEndOfs,NumPathNodes,SortedPathIdx,pAlignNodes,ppFirst2Rpts); 

	switch(m_RsltsFormat) {
		case eBLZRsltsSQLite:
				ReportAsSQLitePSL(Matches,MisMatches,0,NumbNs,						// TotalAlignLen -TotalMMs,TotalMMs,0,0,
									qNumInsert,								// qNumInsert, Number of inserts in query
									qBaseInsert,							// qBaseInsert, Number of bases inserted in query
									tNumInsert,								// tNumInsert, Number of inserts in target
									tBaseInsert,							// tBaseInsert, Number of bases inserted in target
									bStrand == true ? '-' : '+',
									pszQuerySeqIdent,QueryLen,QueryPathStartOfs,QueryPathEndOfs,
									szTargName,TargSeqLen,TargPathStartOfs,TargPathEndOfs,
									NumPathNodes,SortedPathIdx,pAlignNodes,ppFirst2Rpts);
			break;



		case eBLZRsltsPSL:
			ReportAsPSL(Matches,MisMatches,0,NumbNs,						// TotalAlignLen -TotalMMs,TotalMMs,0,0,
									qNumInsert,								// qNumInsert, Number of inserts in query
									qBaseInsert,							// qBaseInsert, Number of bases inserted in query
									tNumInsert,								// tNumInsert, Number of inserts in target
									tBaseInsert,							// tBaseInsert, Number of bases inserted in target
									bStrand == true ? '-' : '+',
									pszQuerySeqIdent,QueryLen,QueryPathStartOfs,QueryPathEndOfs,
									szTargName,TargSeqLen,TargPathStartOfs,TargPathEndOfs,
									NumPathNodes,SortedPathIdx,pAlignNodes,ppFirst2Rpts);
			break;

		case eBLZRsltsPSLX:
			ReportAsPSLX(Matches,MisMatches,0,NumbNs,						// TotalAlignLen -TotalMMs,TotalMMs,0,0,
									qNumInsert,								// qNumInsert, Number of inserts in query
									qBaseInsert,							// qBaseInsert, Number of bases inserted in query
									tNumInsert,								// tNumInsert, Number of inserts in target
									tBaseInsert,							// tBaseInsert, Number of bases inserted in target
									bStrand == true ? '-' : '+',
									pszQuerySeqIdent,QueryLen,pQuerySeq,QueryPathStartOfs,QueryPathEndOfs,
									pCurNode->TargSeqID,
									szTargName,TargSeqLen,TargPathStartOfs,TargPathEndOfs,
									NumPathNodes,SortedPathIdx,pAlignNodes,ppFirst2Rpts);
			break;

		case eBLZRsltsMAF:
			ReportAsMAF(PathHiScore,Matches,MisMatches,0,NumbNs,						// TotalAlignLen -TotalMMs,TotalMMs,0,0,
									qNumInsert,								// qNumInsert, Number of inserts in query
									qBaseInsert,							// qBaseInsert, Number of bases inserted in query
									tNumInsert,								// tNumInsert, Number of inserts in target
									tBaseInsert,							// tBaseInsert, Number of bases inserted in target
									bStrand == true ? '-' : '+',
									pszQuerySeqIdent,QueryLen,pQuerySeq,QueryPathStartOfs,QueryPathEndOfs,
									pCurNode->TargSeqID,
									szTargName,TargSeqLen,TargPathStartOfs,TargPathEndOfs,
									NumPathNodes,SortedPathIdx,pAlignNodes,ppFirst2Rpts);
			break;


		case eBLZRsltsBED:
			ReportAsBED(pszQuerySeqIdent,bStrand == true ? '-' : '+',PathHiScore,szTargName,NumPathNodes,TargPathStartOfs,TargPathEndOfs,SortedPathIdx,pAlignNodes,ppFirst2Rpts);
			break;
		}
	}
return(SortedPathIdx);
}

#ifdef _WIN32
unsigned __stdcall LoadQuerySeqsFileThread(void * pThreadPars)
#else
void *LoadQuerySeqsFileThread(void * pThreadPars)
#endif
{
int Rslt;
tsLoadQuerySeqsThreadPars *pPars = (tsLoadQuerySeqsThreadPars *)pThreadPars;			// makes it easier not having to deal with casts!
CBlitz *pBlitzer = (CBlitz *)pPars->pThis;

Rslt = pBlitzer->ProcLoadQuerySeqsFile(pPars);
pPars->Rslt = Rslt;
#ifdef _WIN32
_endthreadex(0);
return(eBSFSuccess);
#else
pthread_exit(NULL);
#endif
}

int
CBlitz::InitLoadQuerySeqs(void)
{
tsLoadQuerySeqsThreadPars ThreadPars;

// initiate loading the reads
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Loading query sequences from file...");
m_ThreadLoadQuerySeqsRslt = -1;

ThreadPars.pRslt = &m_ThreadLoadQuerySeqsRslt;
ThreadPars.pThis = this;
ThreadPars.Rslt = 0;

#ifdef _WIN32
m_hThreadLoadQuerySeqs = ThreadPars.threadHandle = (HANDLE)_beginthreadex(NULL,0x0fffff,LoadQuerySeqsFileThread,&ThreadPars,0,&m_ThreadLoadQuerySeqsID);
#else
int ThreadRslt = ThreadPars.threadRslt = pthread_create (&m_ThreadLoadQuerySeqsID , NULL , LoadQuerySeqsFileThread , &ThreadPars );
#endif

// wait a few seconds, if major problems with loading reads then should show very quickly
#ifdef _WIN32
if(WAIT_TIMEOUT != WaitForSingleObject(m_hThreadLoadQuerySeqs, 3000))
	{
	CloseHandle(m_hThreadLoadQuerySeqs);
	m_hThreadLoadQuerySeqs = NULL;
	return(m_ThreadLoadQuerySeqsRslt);
	}
#else
struct timespec ts;
int JoinRlt;
clock_gettime(CLOCK_REALTIME, &ts);
ts.tv_sec += 3;
if((JoinRlt = pthread_timedjoin_np(m_ThreadLoadQuerySeqsID, NULL, &ts)) == 0)
	{
	m_ThreadLoadQuerySeqsID = 0;
	return(m_ThreadLoadQuerySeqsRslt);
	}
#endif
return(eBSFSuccess);
}

int									// returned enqueued query identifier
CBlitz::EnqueueQuerySeq(char *pszQueryIdent,    // query identifier parsed from fasta descriptor
			int QuerySeqLen,		// query sequence length
			UINT8 *pQuerySeq)       // query sequence
{
int Idx;
int SeqID;
UINT8 *pSeq;
tsQuerySeq *psQuery;
if(m_TermBackgoundThreads != 0)	// need to immediately self-terminate?
	return(0);
if((pSeq = new UINT8 [QuerySeqLen]) == NULL)
	return(eBSFerrMem);
memcpy(pSeq,pQuerySeq,QuerySeqLen);

// any room left in query sequence queue?
while(1) {
	AcquireLock(true);
	if((m_NumQuerySeqs + 1) < m_AllocdQuerySeqs)
		break;
	ReleaseLock(true);
	if(m_TermBackgoundThreads != 0)	// need to immediately self-terminate?
		{
		delete pSeq;
		return(0);
		}
#ifdef _WIN32
	Sleep(1000);
#else
	sleep(1);
#endif
	}
Idx = (m_NxtQuerySeqIdx + m_NumQuerySeqs) % m_AllocdQuerySeqs;
psQuery = &m_pQuerySeqs[Idx];
SeqID = ++m_TotSeqIDs;
psQuery->SeqID = SeqID;
psQuery->pQuerySeq = pSeq;
psQuery->QuerySeqLen = QuerySeqLen; 
strncpy(psQuery->szQueryIdent,pszQueryIdent,cMaxQuerySeqIdentLen);
psQuery->szQueryIdent[cMaxQuerySeqIdentLen] = '\0';
m_NumQuerySeqs += 1;
ReleaseLock(true);
return(SeqID);
}

UINT8 *										// returned dequeued sequence, caller is responsible for deleting memory allocated to hold the returned sequence (delete pRetSeq;)
CBlitz::DequeueQuerySeq(int WaitSecs,		// if no sequences available to be dequeued then wait at most this many seconds for a sequence to become available
			int MaxLenQueryIdent,			// maximum length query identifier
			int *pSeqID,					// returned sequence identifier
			char *pszQueryIdent,			// where to return query identifier
			int *pQuerySeqLen)				// where to return query sequence length
{
bool bAllQuerySeqsLoaded;
UINT8 *pSeq;
tsQuerySeq *psQuery;

// any sequences available to be dequeued?
if(WaitSecs < 0)
	WaitSecs = 0;
while(1) {
	if(m_TermBackgoundThreads != 0)	// need to immediately self-terminate?
		return(NULL);
	AcquireLock(true);
	if(m_NumQuerySeqs)
		break;
	bAllQuerySeqsLoaded = m_bAllQuerySeqsLoaded;
	ReleaseLock(true);
	if(bAllQuerySeqsLoaded || WaitSecs-- <= 0 || m_TermBackgoundThreads != 0)	// need to immediately self-terminate?
		return(NULL);
#ifdef _WIN32
	Sleep(1000);
#else
	sleep(1);
#endif
	}
psQuery = &m_pQuerySeqs[m_NxtQuerySeqIdx++];
if(m_NxtQuerySeqIdx == m_AllocdQuerySeqs)
	m_NxtQuerySeqIdx = 0;
*pSeqID = psQuery->SeqID;
pSeq = psQuery->pQuerySeq;
psQuery->pQuerySeq = NULL;
*pQuerySeqLen = psQuery->QuerySeqLen; 
strncpy(pszQueryIdent,psQuery->szQueryIdent,MaxLenQueryIdent);
pszQueryIdent[MaxLenQueryIdent-1] = '\0';
m_NumQuerySeqs -= 1;

if(m_RsltsFormat == eBLZRsltsSQLite)
	m_pSQLitePSL->AddAlignSummary(m_ExprID,pszQueryIdent,*pQuerySeqLen,NULL,0);
ReleaseLock(true);
return(pSeq);
}

int
CBlitz::ProcLoadQuerySeqsFile(tsLoadQuerySeqsThreadPars *pPars)
{
CFasta Fasta;
unsigned char *pSeqBuff;
unsigned char *pMskBase;
UINT32 MskIdx;
size_t BuffOfs;
size_t AllocdBuffSize;
size_t AvailBuffSize;
char szName[cBSFSourceSize];
char szDescription[cBSFDescriptionSize];
UINT32 SeqLen;
int Descrlen;
bool bFirstEntry;
bool bEntryCreated;
bool bTruncSeq;
int Rslt;
int SeqID;
int *pRslt = pPars->pRslt;
AcquireLock(true);
m_bAllQuerySeqsLoaded = false;
m_LoadQuerySeqsRslt = eBSFSuccess;		// presumed success, changed if any processing errors
ReleaseLock(true);

if((Rslt=Fasta.Open(m_pszInputFile,true))!=eBSFSuccess)
	{
	if(Rslt != eBSFerrNotFasta)
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcLoadQuerySeqsFile: Unable to open '%s' [%s] %s",m_pszInputFile,Fasta.ErrText((teBSFrsltCodes)Rslt),Fasta.GetErrMsg());
	AcquireLock(true);
	m_bAllQuerySeqsLoaded = true;
	m_LoadQuerySeqsRslt = (teBSFrsltCodes)Rslt;		
	ReleaseLock(true);
	return(Rslt);
	}

// note malloc is used as can then simply realloc to expand as may later be required
AllocdBuffSize = (size_t)cAllocQuerySeqLen;
if((pSeqBuff = (unsigned char *)malloc(AllocdBuffSize)) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcLoadQuerySeqsFile:- Unable to allocate memory (%u bytes) for sequence buffer",(UINT32)cAllocQuerySeqLen);
	Fasta.Close();
	*pRslt = eBSFerrMem;
	AcquireLock(true);
	m_bAllQuerySeqsLoaded = true;
	m_LoadQuerySeqsRslt = (teBSFrsltCodes)Rslt;		
	ReleaseLock(true);
	return(eBSFerrMem);
	}
AvailBuffSize = cAllocQuerySeqLen;

bFirstEntry = true;
bEntryCreated = false;
bTruncSeq = false;
SeqID = 0;
BuffOfs = 0;
while((Rslt = SeqLen = Fasta.ReadSequence(&pSeqBuff[BuffOfs],(int)min(AvailBuffSize,(size_t)cMaxQuerySeqLen),true,false)) > eBSFSuccess)
	{
	if(m_TermBackgoundThreads != 0)	// requested to immediately self-terminate?
		{
		Rslt = eBSErrSession;
		break;
		}

	if(SeqLen == eBSFFastaDescr)		// just read a descriptor line
		{
		SeqID++;
		if(bEntryCreated)				// add any previous entry
			{
			if((Rslt=EnqueueQuerySeq(szName,(int)BuffOfs,pSeqBuff)) <= eBSFSuccess)
				break;
			}
		Descrlen = Fasta.ReadDescriptor(szDescription,cBSFDescriptionSize);
		// An assumption - will one day bite real hard - is that the
		// fasta descriptor line starts with some form of unique identifier.
		// Use this identifier as the entry name.
		if(sscanf(szDescription," %s[ ,]",szName)!=1)
			sprintf(szName,"%s.%d",m_pszInputFile,++SeqID);

		bFirstEntry = false;
		bEntryCreated = true;
		bTruncSeq = false;
		BuffOfs = 0;
		continue;
		}
	else
		if(bFirstEntry)	// if there was no descriptor then dummy up one...
			{
			SeqID++;
			sprintf(szName,"%s.%d",m_pszInputFile,SeqID);
			strcpy(szDescription,"No Description provided");
			bFirstEntry = false;
			bEntryCreated = true;
			bTruncSeq = false;
			}
	if(bTruncSeq)
		continue;

	// remove any repeat masking flags
	pMskBase = &pSeqBuff[BuffOfs];
	int SeqNs = 0;
	for(MskIdx = 0; MskIdx < SeqLen; MskIdx++,pMskBase++)
		*pMskBase &= ~cRptMskFlg;

	BuffOfs += SeqLen;

	if(BuffOfs > cMaxQuerySeqLen)	// truncate at cMaxQuerySeqLen
		{
		gDiagnostics.DiagOut(eDLWarn,gszProcName,"ProcLoadQuerySeqsFile:- Truncating overlength query sequence '%s' to %d",szName,cMaxQuerySeqLen);
		BuffOfs = cMaxQuerySeqLen;
		AvailBuffSize = AllocdBuffSize - BuffOfs;
		bTruncSeq = true;
		continue;
		}
	AvailBuffSize -= SeqLen;

	if(AvailBuffSize < (size_t)(cAllocQuerySeqLen / 2))
		{
		size_t NewSize = (size_t)cAllocQuerySeqLen + AllocdBuffSize;
		unsigned char *pTmp;
		if((pTmp = (unsigned char *)realloc(pSeqBuff,NewSize))==NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcLoadQuerySeqsFile:- Unable to reallocate memory (%u bytes) for sequence buffer",(UINT32)NewSize);
			Rslt = eBSFerrMem;
			break;
			}
		pSeqBuff = pTmp;
		AllocdBuffSize = NewSize;
		AvailBuffSize = AllocdBuffSize - BuffOfs;
		}
	}
if(Rslt < eBSFSuccess && Rslt != eBSErrSession)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcLoadQuerySeqsFile: Parsing errors");
	while(Fasta.NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,Fasta.GetErrMsg());
	}

if(Rslt >= eBSFSuccess && bEntryCreated && BuffOfs > 0)			// last entry
	Rslt=EnqueueQuerySeq(szName,(int)BuffOfs,pSeqBuff);
if(Rslt > eBSFSuccess)
	Rslt = eBSFSuccess;
if(pSeqBuff != NULL)
	free(pSeqBuff);
*pRslt = Rslt;
AcquireSerialise();
m_bAllQuerySeqsLoaded = true;
m_LoadQuerySeqsRslt = (teBSFrsltCodes)Rslt;
ReleaseSerialise();
return(Rslt);
}

// SortQueryAlignNodes
// Sort alignment nodes by TargSeqID.QueryID.FlgStrand.QueryStartOfs.TargStartOfs ascending
int
CBlitz::SortQueryAlignNodes(const void *arg1, const void *arg2)
{
tsQueryAlignNodes *pEl1 = (tsQueryAlignNodes *)arg1;
tsQueryAlignNodes *pEl2 = (tsQueryAlignNodes *)arg2;

if(pEl1->TargSeqID > pEl2->TargSeqID)
	return(1);
if(pEl1->TargSeqID < pEl2->TargSeqID)
	return(-1);
if(pEl1->QueryID > pEl2->QueryID)
	return(1);
if(pEl1->QueryID < pEl2->QueryID)
	return(-1);
if(pEl1->FlgStrand > pEl2->FlgStrand)
	return(1);
if(pEl1->FlgStrand < pEl2->FlgStrand)
	return(-1);
if(pEl1->QueryStartOfs > pEl2->QueryStartOfs)
	return(1);
if(pEl1->QueryStartOfs < pEl2->QueryStartOfs)
	return(-1);
if(pEl1->TargStartOfs > pEl2->TargStartOfs)
	return(1);
if(pEl1->TargStartOfs < pEl2->TargStartOfs)
	return(-1);
return(0);
}

// SortHighScoreDescend
// Sort alignment nodes which are the first in path by score descending
int
CBlitz::SortHighScoreDescend(const void *arg1, const void *arg2)
{
tsQueryAlignNodes *pEl1 = *(tsQueryAlignNodes **)arg1;
tsQueryAlignNodes *pEl2 = *(tsQueryAlignNodes **)arg2;
if(pEl1->HiScore < pEl2->HiScore)
	return(1);
if(pEl1->HiScore > pEl2->HiScore)
	return(-1);
if(pEl1->TargSeqID > pEl2->TargSeqID)
	return(1);
if(pEl1->TargSeqID < pEl2->TargSeqID)
	return(-1);
return(0);
}




