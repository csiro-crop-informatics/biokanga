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

#include "AlignsBootstrap.h"

int
Process(ePMBSAlign PMode,					// bootstrap processing mode
				int RandSeed,				// if > 0 then random generator seed, otherwise current time used as the seed
				bool bSenseOnly,			// true if to align sense only, false if to align both sense and antisense
				int MaxSubs,				// allowing at most this many subs as percentage of query length before accepting alignment
				int NumBootstraps,			// number of bootstrap iterations, excludes initial original query sequences aligned onto initial target sequences
				bool bWORreplacement,		// sample without replacement, default is to sample with replacement
				bool bNoOverlaps,			// sample without overlap, default is to sample allowing overlaps
				char *pszQuerySeqsFile,		// fasta file containing initial query sequences from which to derive query length distributions
				char *pszTargSeqsFile,		// fasta file containing initial target sequences from which to derive target length distributions
				char *pszQueryAssembFile,	// file containing fasta assembly to be bootstrap sampled for query sequences with same length distributions as sequences in pszQuerySeqsFile 
				char *pszTargAssembFile,	// file containing fasta assembly to be bootstrap sampled for target sequences with same length distributions as sequences in pszTargSeqsFile
				char *pszQRsltsFile,		// summary number of query hits onto at least one target bootstrap results to this file 
				char *pszTRsltsFile,		// summary number of targets hit by at least one query bootstrap results to this file  
			    int NumThreads);			// number of worker threads to use



#ifdef _WIN32
int AlignsBootstrap(int argc, char* argv[])
{
// determine my process name
_splitpath(argv[0],NULL,NULL,gszProcName,NULL);
#else
int
AlignsBootstrap(int argc, char** argv)
{
// determine my process name
CUtility::splitpath((char *)argv[0],NULL,gszProcName);
#endif

int iFileLogLevel;			// level of file diagnostics
int iScreenLogLevel;		// level of file diagnostics
char szLogFile[_MAX_PATH];	// write diagnostics to this file
int Rslt = 0;   			// function result code >= 0 represents success, < 0 on failure

ePMBSAlign PMode;			// processing mode
int RandSeed;				// if > 0 then random generator seed, otherwise current time used as the seed

bool bWORreplacement;	// sample without replacement, default is to sample with replacement
bool bNoOverlaps;		// sample without overlap, default is to sample allowing overlaps

int NumberOfProcessors;		// number of installed CPUs
int NumThreads;				// number of threads (0 defaults to number of CPUs)

bool bSenseOnly;			// true if to align sense only, false if to align both sense and antisense
int MaxSubs;				// allowing at most this many subs as percentage of query length before accepting alignment

int NumBootstraps;			// number of bootstrap iterations, excludes initial original query sequences aligned onto initial target sequences
char szQuerySeqsFile[_MAX_PATH];		// fasta file containing initial query sequences from which to derive query length distributions
char szTargSeqsFile[_MAX_PATH];		// fasta file containing initial target sequences from which to derive target length distributions
char szQueryAssembFile[_MAX_PATH];	// file containing fasta assembly to be bootstrap sampled for query sequences with same length distributions as sequences in pszQuerySeqsFile 
char szTargAssembFile[_MAX_PATH];	// file containing fasta assembly to be bootstrap sampled for target sequences with same length distributions as sequences in pszTargSeqsFile
char szQRsltsFile[_MAX_PATH];		// summary number of query hits onto at least one target bootstrap results to this file 
char szTRsltsFile[_MAX_PATH];		// summary number of targets hit by at least one query bootstrap results to this file   

char szSQLiteDatabase[_MAX_PATH];	// results summaries to this SQLite file
char szExperimentName[cMaxDatasetSpeciesChrom+1];			// experiment name
char szExperimentDescr[1000];		// describes experiment

//
struct arg_lit  *help    = arg_lit0("h","help",                 "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to screen and logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"diagnostics log file");

struct arg_int *pmode = arg_int0("m","mode","<int>",		    "bootstrap processing mode: 0 - standard, 1 - report bootstrap sequences used in each iteration");

struct arg_int *randseed = arg_int0("r","randseed","<int>",		    "if > 0 then random generator seed, otherwise current time used as the seed (default 0)");


struct arg_int *numbootstraps = arg_int0("b","numbootstraps","<int>",	"number of bootstrap iterations - excluding initial original query seqs aligned to original target seqs (default 1000, range 1..10000)");

struct arg_int *maxsubs = arg_int0("s","maxsubs","<int>",	"allowing at most this many subs as percentage of query length before accepting alignment (default 0, range 0..50");

struct arg_file *queryseqsfile = arg_file1("p","queryseqsfile","<file>",		"fasta file containing initial query sequences from which to derive query length distributions");
struct arg_file *queryassembfile = arg_file1("P","queryassembfile","<file>",	"file containing fasta assembly to be bootstrap sampled for query sequences with same length distributions as sequences in queryseqsfile");

struct arg_file *targseqsfile = arg_file1("i","targseqsfile","<file>",		"fasta file containing initial target sequences from which to derive target length distributions");
struct arg_file *targassembfile = arg_file1("I","targassembfile","<file>",	"file containing fasta assembly to be bootstrap sampled for target sequences with same length distributions as sequences in targassembfile");

struct arg_lit  *senseonly    = arg_lit0("a","senseonly",                 "align sense only, default is to align both sense and antisense");

struct arg_lit  *woreplacement    = arg_lit0("x","woreplacement",    "sample without replacement, default is to sample with replacement");
struct arg_lit  *nooverlaps    = arg_lit0("X","nooverlaps",         "sample without overlap, default is to sample allowing overlaps");


struct arg_file *qrsltsfile = arg_file1("o","qrsltsfile","<file>",		"summary number of query hits onto at least one target bootstrap results to this file");
struct arg_file *trsltsfile = arg_file1("O","trsltsfile","<file>",		"summary number of targets hit by at least one query bootstrap results to this file");

struct arg_int *threads = arg_int0("T","threads","<int>",		"number of processing threads 0..128 (defaults to 0 which sets threads to number of CPU cores)");

struct arg_file *summrslts = arg_file0("q","sumrslts","<file>",		"Output results summary to this SQLite3 database file");
struct arg_str *experimentname = arg_str0("w","experimentname","<str>",		"experiment name SQLite3 database file");
struct arg_str *experimentdescr = arg_str0("W","experimentdescr","<str>",	"experiment description SQLite3 database file");

struct arg_end *end = arg_end(200);

void *argtable[] = {help,version,FileLogLevel,LogFile,
					summrslts,experimentname,experimentdescr,
					pmode,randseed,senseonly,maxsubs,numbootstraps,woreplacement,nooverlaps,queryseqsfile,queryassembfile,targseqsfile,targassembfile,qrsltsfile,trsltsfile,threads,
					end};

char **pAllArgs;
int argerrors;
argerrors = CUtility::arg_parsefromfile(argc,(char **)argv,&pAllArgs);
if(argerrors >= 0)
	argerrors = arg_parse(argerrors,pAllArgs,argtable);

/* special case: '--help' takes precedence over error reporting */
/* special case: '--help' takes precedence over error reporting */
if (help->count > 0)
        {
		printf("\n%s %s %s, Version %s\nOptions ---\n", gszProcName,gpszSubProcess->pszName,gpszSubProcess->pszFullDescr,cpszProgVer);
        arg_print_syntax(stdout,argtable,"\n");
        arg_print_glossary(stdout,argtable,"  %-25s %s\n");
		printf("\nNote: Parameters can be entered into a parameter file, one parameter per line.");
		printf("\n      To invoke this parameter file then precede it's name with '@'");
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

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Subprocess %s Version %s starting",gszProcName,cpszProgVer);
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
	else
		szExperimentName[0] = '\0';

	gExperimentID = 0;
	gProcessID = 0;
	gProcessingID = 0;
	szSQLiteDatabase[0] = '\0';
	szExperimentDescr[0] = '\0';

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

		if(strlen(szExperimentName) < 1)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: After removal of whitespace, no SQLite experiment name specified with '-w<str>' option");
			return(1);
			}
		if(experimentdescr->count)
			{
			strncpy(szExperimentDescr,experimentdescr->sval[0],sizeof(szExperimentDescr)-1);
			szExperimentDescr[sizeof(szExperimentDescr)-1] = '\0';
			CUtility::TrimQuotedWhitespcExtd(szExperimentDescr);
			}
		if(strlen(szExperimentDescr) < 1)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: After removal of whitespace, no SQLite experiment description specified with '-W<str>' option");
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
		{
		szSQLiteDatabase[0] = '\0';
		szExperimentDescr[0] = '\0';
		}

	PMode = (ePMBSAlign)(pmode->count ? pmode->ival[0] : (int)ePMBSAdefault);
	if(PMode < ePMBSAdefault || PMode >= ePMBSAPlaceholder)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Processing mode '-m%d' specified outside of range %d..%d\n",PMode,ePMBSAdefault,(int)ePMBSAPlaceholder-1);
		exit(1);
		}


	bWORreplacement=woreplacement->count ? true : false;
	bNoOverlaps=nooverlaps->count ? true : false;

	RandSeed = randseed->count ? randseed->ival[0] : 0;
	if(RandSeed <= 0)
		{
#ifdef _WIN32
		INT64 Now;		
		QueryPerformanceCounter((LARGE_INTEGER *)&Now);
		RandSeed = (int)(Now & 0x07fffffff);
#else
		struct timeval TimeNow;
		gettimeofday(&TimeNow,NULL);
		RandSeed = (int)(TimeNow.tv_usec & 0x07fffffff);
#endif
		}

	bSenseOnly = senseonly->count ? true : false;

	NumBootstraps = numbootstraps->count ? numbootstraps->ival[0] : cDfltBootstraps;
	if(NumBootstraps < cMinBootstraps || NumBootstraps > cMaxBootstraps)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Number of query bootstrap iterations '-b%d' specified outside of range %d..%d\n",NumBootstraps,cMinBootstraps,cMaxBootstraps);
		exit(1);
		}

	MaxSubs = maxsubs->count ? maxsubs->ival[0] : 0;
	if(MaxSubs < 0 || MaxSubs > 50)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Maximum allowed subs in accepted alignments '-s%d' specified outside of range %d..%d\n",MaxSubs,0,50);
		exit(1);
		}

	strcpy(szQuerySeqsFile,queryseqsfile->filename[0]);
	strcpy(szTargSeqsFile,targseqsfile->filename[0]);
	strcpy(szQueryAssembFile,queryassembfile->filename[0]);
	strcpy(szTargAssembFile,targassembfile->filename[0]);
	strcpy(szQRsltsFile,qrsltsfile->filename[0]);
	strcpy(szTRsltsFile,trsltsfile->filename[0]);

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

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing parameters:");

	const char *pszDescr;
	switch(PMode) {
		case ePMSAreportseqs:
			pszDescr = "Standard bootstrap processing with reporting of sequences used in bootstraps excluding original sequences";
			break;

		case ePMBSAdefault:
		default:
			pszDescr = "Standard bootstrap processing";
			break;
		}
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Processing is : '%s'",pszDescr);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Sampling %s replacement", bWORreplacement ? "without" : "with");
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Sampling %s allowing overlaps", bNoOverlaps ? "without" : "with");

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Random number generator using seed: %d",RandSeed);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Aligning: '%s'",bSenseOnly ? "Sense only" : "Sense and antisense");

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Number of bootstrap iterations: %d",NumBootstraps);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Max allowed subs in accepted alignments as percentage of query length: %d",MaxSubs);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"input query sequences file: '%s'",szQuerySeqsFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"input target sequences file: '%s'",szTargSeqsFile);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"input query assembly file: '%s'",szQueryAssembFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"input target assembly file: '%s'",szTargAssembFile);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"output query results file: '%s'",szQRsltsFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"output target results file: '%s'",szTRsltsFile);

	if(szExperimentName[0] != '\0')
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"This processing reference: %s",szExperimentName);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"number of threads : %d",NumThreads);

	if(gExperimentID > 0)
		{
		int ParamID;
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szLogFile),"log",szLogFile);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(PMode),"mode",&PMode);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(NumBootstraps),"numbootstraps",&NumBootstraps);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(RandSeed),"randseed",&RandSeed);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(MaxSubs),"maxsubs",&MaxSubs);

		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(bWORreplacement),"woreplacement",&bWORreplacement);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(bNoOverlaps),"nooverlaps",&bNoOverlaps);

		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szQuerySeqsFile),"queryseqsfile",szQuerySeqsFile);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szTargSeqsFile),"targseqsfile",szTargSeqsFile);

		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szQueryAssembFile),"queryassembfile",szQueryAssembFile);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szTargAssembFile),"targassembfile",szTargAssembFile);

		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szQRsltsFile),"qrsltsfile",szQRsltsFile);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szTRsltsFile),"trsltsfile",szTRsltsFile);
		
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(NumThreads),"threads",&NumThreads);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(NumberOfProcessors),"cpus",&NumberOfProcessors);

		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szSQLiteDatabase),"sumrslts",szSQLiteDatabase);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szExperimentName),"experimentname",szExperimentName);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szExperimentDescr),"experimentdescr",szExperimentDescr);
		}


#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	gStopWatch.Start();
	Rslt = Process((ePMBSAlign)PMode,RandSeed,bSenseOnly,MaxSubs,NumBootstraps,bWORreplacement,bNoOverlaps,szQuerySeqsFile,szTargSeqsFile,szQueryAssembFile,szTargAssembFile,szQRsltsFile,szTRsltsFile,NumThreads);
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
Process(ePMBSAlign PMode,					// bootstrap processing mode
				int RandSeed,				// if > 0 then random generator seed, otherwise current time used as the seed
				bool bSenseOnly,			// true if to align sense only, false if to align both sense and antisense
				int MaxSubs,				// allowing at most this many subs as percentage of query length before accepting alignment
				int NumBootstraps,			// number of bootstrap iterations, excludes initial original query sequences aligned onto initial target sequences
				bool bWORreplacement,		// sample without replacement, default is to sample with replacement
				bool bNoOverlaps,			// sample without overlap, default is to sample allowing overlaps
				char *pszQuerySeqsFile,		// file containing initial query sequences from which to derive query length distributions - can be fasta, bed, or csv 
				char *pszTargSeqsFile,		// file containing initial target sequences from which to derive target length distributions - can be fasta, bed, or csv
				char *pszQueryAssembFile,	// file containing assembly to be bootstrap sampled for query sequences with same length distributions as sequences in pszQuerySeqsFile 
				char *pszTargAssembFile,	// file containing assembly to be bootstrap sampled for target sequences with same length distributions as sequences in pszTargSeqsFile
				char *pszQRsltsFile,		// summary number of query hits onto target bootstrap results to this file 
				char *pszTRsltsFile,		// summary number of targets hit by at least one query bootstrap results to this file  
			    int NumThreads)				// number of worker threads to use
{
int Rslt;
CAlignsBootstrap *pBootstrapper;
if((pBootstrapper = new CAlignsBootstrap) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: Unable to instantiate instance of CAlignsBootstrap");
	return(eBSFerrObj);
	}
Rslt = pBootstrapper->Process(PMode,RandSeed,bSenseOnly,MaxSubs,NumBootstraps,bWORreplacement,bNoOverlaps,pszQuerySeqsFile,pszTargSeqsFile,pszQueryAssembFile,pszTargAssembFile,pszQRsltsFile,pszTRsltsFile,NumThreads);
delete pBootstrapper;
return(Rslt);
}



CAlignsBootstrap::CAlignsBootstrap()
{
memset(m_Seqs,0,sizeof(m_Seqs));
memset(m_WorkerInstances,0,sizeof(m_WorkerInstances));
m_hCSVQRslts = -1;
m_hCSVTRslts = -1;
m_pQueryHits = NULL;
m_pszQRsltsBuff = NULL;
m_pszTRsltsBuff = NULL;
m_pRandomMersenne = NULL;
}


CAlignsBootstrap::~CAlignsBootstrap()
{
Reset();
}

void
CAlignsBootstrap::Reset(void)
{
int SrcIdx;
tsSeqAllocs *pSeq;

if(m_hCSVQRslts != -1)
	{
#ifdef _WIN32
	_commit(m_hCSVQRslts);
#else
	fsync(m_hCSVQRslts);
#endif
	close(m_hCSVQRslts);
	m_hCSVQRslts = -1;
	}
if(m_hCSVTRslts != -1)
	{
#ifdef _WIN32
	_commit(m_hCSVTRslts);
#else
	fsync(m_hCSVTRslts);
#endif
	close(m_hCSVTRslts);
	m_hCSVTRslts = -1;
	}
if(m_pszQRsltsBuff != NULL)
	{
	free(m_pszQRsltsBuff);
	m_pszQRsltsBuff = NULL;
	}
if(m_pszTRsltsBuff != NULL)
	{
	free(m_pszTRsltsBuff);
	m_pszTRsltsBuff = NULL;
	}

pSeq = m_Seqs;
for(SrcIdx = 0; SrcIdx < ePMBSSrcPlaceholder; SrcIdx++, pSeq++)
	{
	if(pSeq->pSeqs != NULL)
		{
#ifdef _WIN32
		free(pSeq->pSeqs);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
		if(pSeq->pSeqs != MAP_FAILED)
			munmap(pSeq->pSeqs,pSeq->AllocdSeqsSize);
#endif
		pSeq->pSeqs = NULL;
		}

	if(pSeq->pSeqBlocks != NULL)
		{
#ifdef _WIN32
		free(pSeq->pSeqBlocks);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
		if(pSeq->pSeqBlocks != MAP_FAILED)
			munmap(pSeq->pSeqBlocks,pSeq->SeqBlocksAllocSize);
#endif
		pSeq->pSeqBlocks = NULL;
		}

	if(pSeq->pSeqDescrs != NULL)
		{
#ifdef _WIN32
		free(pSeq->pSeqDescrs);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
		if(pSeq->pSeqDescrs != MAP_FAILED)
			munmap(pSeq->pSeqDescrs,pSeq->AllocdSeqDescrsSize);
#endif
		pSeq->pSeqDescrs = NULL;
		}
	}
memset(m_Seqs,0,sizeof(m_Seqs));

if(m_pQueryHits != NULL)
	{
	free(m_pQueryHits);
	m_pQueryHits = NULL;
	}

if(m_pRandomMersenne != NULL)
	{
	delete m_pRandomMersenne;
	m_pRandomMersenne = NULL;
	}
m_PMode = ePMBSAdefault;
m_RandSeed = 0;
m_bSenseOnly = false;
m_bWORreplacement = false;
m_bNoOverlaps = false;
m_MaxSubs = 1;
m_NumBootstraps = cDfltBootstraps;
m_NumThreads = cMaxWorkerThreads;
m_NumWorkerInsts = 0;
m_CompletedWorkerInsts = 0;
m_AlignReqID = 0;
m_TermAllThreads = 0;
m_CurBootstrap = 0;
memset(m_WorkerInstances,0,sizeof(m_WorkerInstances));
}

int
CAlignsBootstrap::Init(int RandSeed)		// if > 0 then random generator seed, otherwise current time used as the seed
{
size_t memreq;
int SrcIdx;
tsSeqAllocs *pSeq;

Reset();
pSeq = m_Seqs;
for(SrcIdx = 0; SrcIdx < ePMBSSrcPlaceholder; SrcIdx++,pSeq++)
	{
	pSeq->SeqSrc = (ePMBSSeqSrc)SrcIdx;

	memreq = ((size_t)cNumSeqDescrs * (sizeof(tsSeqDescr) + cAvgSeqDescrLen));
#ifdef _WIN32
	pSeq->pSeqDescrs = (tsSeqDescr *) malloc(memreq);	// initial and perhaps the only allocation
	if(pSeq->pSeqDescrs == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Init: Memory allocation of %lld bytes for sequence descriptors failed - %s",(INT64)memreq,strerror(errno));
		return(eBSFerrMem);
		}
#else
	// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
	pSeq->pSeqDescrs = (tsSeqDescr *)mmap(NULL,memreq, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
	if(pSeq->pSeqDescrs == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Init: Memory allocation of %lld bytes for sequence descriptors failed - %s",(INT64)memreq,strerror(errno));
		pSeq->pSeqDescrs = NULL;
		return(eBSFerrMem);
		}
#endif
	pSeq->AllocdSeqDescrsSize = memreq;
	pSeq->UsedSeqDescrsSize = 0;
	pSeq->NumSeqDescrs = 0;

	memreq = ((size_t)cDfltSeqBlocks * sizeof(tsSeqBlock));
#ifdef _WIN32
	pSeq->pSeqBlocks = (tsSeqBlock *) malloc(memreq);	// initial and perhaps the only allocation
	if(pSeq->pSeqBlocks == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Init: Memory allocation of %lld bytes for sequence blocks failed - %s",(INT64)memreq,strerror(errno));
		Reset();
		return(eBSFerrMem);
		}
#else
	// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
	pSeq->pSeqBlocks = (tsSeqBlock *)mmap(NULL,memreq, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
	if(pSeq->pSeqBlocks == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Init: Memory allocation of %lld bytes for sequence blocks failed - %s",(INT64)memreq,strerror(errno));
		pSeq->pSeqBlocks = NULL;
		Reset();
		return(eBSFerrMem);
		}
#endif
	pSeq->SeqBlocksAllocSize = memreq;
	pSeq->UsedSeqBlocksSize = 0;
	pSeq->UsedSeqBlocks = 0;


	memreq = ((size_t)cDfltSeqAlloc * sizeof(etSeqBase));
#ifdef _WIN32
	pSeq->pSeqs = (etSeqBase *) malloc(memreq);	// initial and perhaps the only allocation
	if(pSeq->pSeqs == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Init: Memory allocation of %lld bytes for sequence bases failed - %s",(INT64)memreq,strerror(errno));
		Reset();
		return(eBSFerrMem);
		}
#else
	// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
	pSeq->pSeqs = (etSeqBase *)mmap(NULL,memreq, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
	if(pSeq->pSeqs == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Init: Memory allocation of %lld bytes for sequence bases failed - %s",(INT64)memreq,strerror(errno));
		pSeq->pSeqs = NULL;
		Reset();
		return(eBSFerrMem);
		}
#endif
	pSeq->AllocdSeqsSize = memreq;
	pSeq->UsedSeqsSize = 0;
	pSeq->NumSeqs = 0;
	}

if(RandSeed <= 0)
	{
#ifdef _WIN32
	INT64 Now;		
	QueryPerformanceCounter((LARGE_INTEGER *)&Now);
	RandSeed = (int)(Now & 0x07fffffff);
#else
	struct timeval TimeNow;
	gettimeofday(&TimeNow,NULL);
	RandSeed = (int)(TimeNow.tv_usec & 0x07fffffff);
#endif
	}

m_pRandomMersenne = new CRandomMersenne(RandSeed);
return(eBSFSuccess);
}

int
CAlignsBootstrap::AddSeq(ePMBSSeqSrc SeqSrc,     // descriptor source - 0: query seqs, 1: target sequences, 2: query assembly, 3: target assembly
			 char *pszDescr,		// '\0' terminated sequence descriptor
			 UINT32 SeqLen,			// sequence is this length
			 UINT8 *pSeqBuff)		// sequence
{
size_t memreq;
int DescrLen;
tsSeqDescr *pSeqDescrs;
tsSeqAllocs *pSeqAllocs;
tsSeqBlock *pSeqBlocks;
etSeqBase *pSeqs;
pSeqAllocs = &m_Seqs[SeqSrc];

DescrLen = min((int)strlen(pszDescr),cTruncSeqDescrLen);
pszDescr[DescrLen] = '\0';

if(pSeqAllocs->AllocdSeqDescrsSize <= (pSeqAllocs->UsedSeqDescrsSize + sizeof(tsSeqDescr) + DescrLen))
	{
	memreq = ((size_t)cNumSeqDescrs * (sizeof(tsSeqDescr) + cAvgSeqDescrLen));
	memreq += pSeqAllocs->AllocdSeqDescrsSize; 
#ifdef _WIN32
	pSeqDescrs = (tsSeqDescr *) realloc(pSeqAllocs->pSeqDescrs,memreq);
#else
	pSeqDescrs = (tsSeqDescr *)mremap(pSeqAllocs->pSeqDescrs,pSeqAllocs->AllocdSeqDescrsSize,memreq,MREMAP_MAYMOVE);
	if(pSeqDescrs == MAP_FAILED)
		pSeqDescrs = NULL;
#endif
	if(pSeqDescrs == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddSeq: Memory reallocation to %lld bytes for sequence descriptors failed - %s",(INT64)memreq,strerror(errno));
		return(eBSFerrMem);
		}
	pSeqAllocs->AllocdSeqDescrsSize = memreq;
	pSeqAllocs->pSeqDescrs = pSeqDescrs;
	}

if(pSeqAllocs->SeqBlocksAllocSize <= pSeqAllocs->UsedSeqBlocksSize + sizeof(tsSeqBlock))
	{
	memreq = (size_t)cDfltSeqBlocks * sizeof(tsSeqBlock);
	memreq += pSeqAllocs->SeqBlocksAllocSize; 
#ifdef _WIN32
	pSeqBlocks = (tsSeqBlock *) realloc(pSeqAllocs->pSeqBlocks,memreq);
#else
	pSeqBlocks = (tsSeqBlock *)mremap(pSeqAllocs->pSeqBlocks,pSeqAllocs->SeqBlocksAllocSize,memreq,MREMAP_MAYMOVE);
	if(pSeqBlocks == MAP_FAILED)
		pSeqBlocks = NULL;
#endif
	if(pSeqBlocks == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddSeq: Memory reallocation to %lld bytes for sequence blocks failed - %s",(INT64)memreq,strerror(errno));
		return(eBSFerrMem);
		}
	pSeqAllocs->SeqBlocksAllocSize = memreq;
	pSeqAllocs->pSeqBlocks = pSeqBlocks;
	}

if(pSeqAllocs->AllocdSeqsSize <= pSeqAllocs->UsedSeqsSize + (sizeof(etSeqBase) * (SeqLen + 1)))
	{
	memreq = (size_t)cDfltSeqAlloc * sizeof(etSeqBase);
	memreq += pSeqAllocs->AllocdSeqsSize + SeqLen; 
#ifdef _WIN32
	pSeqs = (etSeqBase *) realloc(pSeqAllocs->pSeqs,memreq);
#else
	pSeqs = (etSeqBase *)mremap(pSeqAllocs->pSeqs,pSeqAllocs->AllocdSeqsSize,memreq,MREMAP_MAYMOVE);
	if(pSeqs == MAP_FAILED)
		pSeqs = NULL;
#endif
	if(pSeqs == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddSeq: Memory reallocation to %lld bytes for sequence bases failed - %s",(INT64)memreq,strerror(errno));
		return(eBSFerrMem);
		}
	pSeqAllocs->AllocdSeqsSize = memreq;
	pSeqAllocs->pSeqs = pSeqs;
	}

pSeqs = &pSeqAllocs->pSeqs[pSeqAllocs->UsedSeqsSize];
pSeqDescrs = (tsSeqDescr *)((UINT8 *)pSeqAllocs->pSeqDescrs + pSeqAllocs->UsedSeqDescrsSize);
pSeqBlocks = (tsSeqBlock *)((UINT8 *)pSeqAllocs->pSeqBlocks + pSeqAllocs->UsedSeqBlocksSize);

memcpy(pSeqs,pSeqBuff,SeqLen);
pSeqs[SeqLen] = eBaseEOS;

pSeqAllocs->NumSeqs += 1;
pSeqAllocs->UsedSeqBlocks += 1;

strcpy((char *)pSeqDescrs->Descr,pszDescr);
pSeqDescrs->DescrLen = DescrLen;
pSeqDescrs->SeqSrc = SeqSrc;
pSeqDescrs->SeqBlockID = pSeqAllocs->UsedSeqBlocks;

pSeqBlocks->SeqLen = SeqLen;
pSeqBlocks->SmplSeqOfs = pSeqAllocs->UsedSeqsSize;
pSeqBlocks->PopSeqOfs = -1;
pSeqBlocks->SeqSrc = SeqSrc;
pSeqBlocks->SeqBlockID = pSeqAllocs->UsedSeqBlocks;
pSeqBlocks->DescrOfs = pSeqAllocs->UsedSeqDescrsSize;

pSeqAllocs->UsedSeqsSize += SeqLen + 1;
pSeqAllocs->UsedSeqDescrsSize += sizeof(tsSeqDescr) + DescrLen;
pSeqAllocs->UsedSeqBlocksSize += sizeof(tsSeqBlock);

return(eBSFSuccess);
}

int
CAlignsBootstrap::LoadFastaSeqs(ePMBSSeqSrc SeqSrc,     // descriptor source - 0: query seqs, 1: target sequences, 2: query assembly, 3: target assembly
			int MinSeqLen,								// only accepting sequences which are at least this length
			char *pszFastaFile)							// fasta file to load from
{
int Rslt;
UINT32 SeqIdx;
UINT8 *pSeqBase;
size_t BuffOfs;
size_t AllocdBuffSize;
size_t AvailBuffSize;
UINT8 *pSeqBuff;
char szName[cBSFSourceSize];
char szDescription[cBSFDescriptionSize];
UINT32 SeqLen;
int MaxSeqLen;
int Descrlen;
int SeqID;
UINT32 NumSeqsUnderlength;
UINT32 NumSeqsOverlength;
bool bFirstEntry;
bool bEntryCreated;
CFasta Fasta;

switch(SeqSrc) {
	case ePMBSSQuerySeqs:	// sequence loaded from query sequence
		MaxSeqLen = cMaxQuerySeqLen;
		break;

	case ePMBSSTargSeqs:	// sequence loaded from target sequence
		MaxSeqLen = cMaxTargSeqLen;
		break;

	case ePMBSSQueryAssemb:  // sequence loaded from query assembly
	case ePMBSSTargAssemb:  // sequence loaded from target assembly
		MaxSeqLen = cMaxAssembSeqLen;
		break;
	}

if((Rslt=Fasta.Open(pszFastaFile,true,(UINT32)cMaxAssembSeqLen * 2))!=eBSFSuccess)
	{
	if(Rslt != eBSFerrNotFasta)
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadFastaSeqs: Unable to open '%s' [%s] %s",pszFastaFile,Fasta.ErrText((teBSFrsltCodes)Rslt),Fasta.GetErrMsg());
	return(Rslt);
	}

AllocdBuffSize = (size_t)cMaxAssembSeqLen;
// note malloc is used as can then simply realloc to expand as may later be required
if((pSeqBuff = (unsigned char *)malloc(AllocdBuffSize)) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadFastaSeqs:- Unable to allocate memory (%u bytes) for sequence buffer",(UINT32)AllocdBuffSize);
	Fasta.Close();
	return(eBSFerrMem);
	}
AvailBuffSize = AllocdBuffSize;

bFirstEntry = true;
bEntryCreated = false;
SeqID = 0;
BuffOfs = 0;
NumSeqsUnderlength = 0;
NumSeqsOverlength = 0;
while((Rslt = SeqLen = Fasta.ReadSequence(&pSeqBuff[BuffOfs],(int)min(AvailBuffSize,(size_t)cMaxAssembSeqLen),true,false)) > eBSFSuccess)
	{
	if(SeqLen == eBSFFastaDescr)		// just read a descriptor line
		{
		SeqID++;
		if(bEntryCreated)				// add any previous entry
			{
			if(BuffOfs < (size_t)MinSeqLen)
				NumSeqsUnderlength += 1;
			else
				if(BuffOfs > (size_t)MaxSeqLen)
					NumSeqsOverlength += 1;
				else
					if((Rslt=AddSeq(SeqSrc,szName,(UINT32)BuffOfs,pSeqBuff)) < eBSFSuccess)
						{
						gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadFastaSeqs - AddSeq error %d",Rslt);
						break;
						}
			Rslt = eBSFSuccess;
			}
		Descrlen = Fasta.ReadDescriptor(szDescription,cBSFDescriptionSize);
		// An assumption - will one day bite real hard - is that the
		// fasta descriptor line starts with some form of unique identifier.
		// Use this identifier as the entry name.
		if(sscanf(szDescription," %s[ ,]",szName)!=1)
			sprintf(szName,"%s.%d",pszFastaFile,++SeqID);

		bFirstEntry = false;
		bEntryCreated = true;
		BuffOfs = 0;
		continue;
		}
	else
		if(bFirstEntry)	// if there was no descriptor then dummy up one...
			{
			SeqID++;
			sprintf(szName,"%s.%d",pszFastaFile,SeqID);
			strcpy(szDescription,"No Description provided");
			bFirstEntry = false;
			bEntryCreated = true;
			}

	// not interested in any repeat masking
	pSeqBase = pSeqBuff;
	for(SeqIdx = 0; SeqIdx < SeqLen; SeqIdx++)
		*pSeqBase++ &= ~cRptMskFlg;

	BuffOfs += SeqLen;
	AvailBuffSize = AllocdBuffSize - BuffOfs;
	}

if(Rslt >= eBSFSuccess && bEntryCreated && BuffOfs > 0)			// close entry
	{
	if(BuffOfs < (size_t)MinSeqLen)
		{
		NumSeqsUnderlength += 1;
		Rslt = eBSFSuccess;
		}
	else
		if(BuffOfs > (size_t)MaxSeqLen)
			{
			NumSeqsOverlength += 1;
			Rslt = eBSFSuccess;
			}
		else
			{
			if((Rslt=AddSeq(SeqSrc,szName,(UINT32)BuffOfs,pSeqBuff)) < eBSFSuccess)
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadFastaSeqs - AddSeq error %d",Rslt);
			else
				Rslt = eBSFSuccess;
			}
	}
if(pSeqBuff != NULL)
	free(pSeqBuff);
if(NumSeqsUnderlength > 0 || NumSeqsOverlength > 0)
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"LoadFastaSeqs - removed %u with length under %dbp, %u with length over %dbp  ",NumSeqsUnderlength,MinSeqLen,NumSeqsOverlength,MaxSeqLen);
return(Rslt);
}

INT64		// returned random number will be at most 60bits (2^60)
CAlignsBootstrap::GenRand60(INT64 Limit)	// generate random number between 0 and Limit inclusive where Limit is <= 2^60
{
int RandomLo;
int RandomHi;

if(Limit <= 0x7fffffff)
	return((INT64)m_pRandomMersenne->IRandom(0,(int)Limit));

RandomLo = m_pRandomMersenne->IRandom(0,(int)(Limit & 0x3fffffff));
RandomHi = m_pRandomMersenne->IRandom(0,(int)((Limit >> 30) & 0x3fffffff));

return((INT64)RandomHi << 30 | (INT64)RandomLo);
}


int
CAlignsBootstrap::GenBootstrap(UINT32 MaxBootstrapAttempts, // allow at most this many attempts at bootstrapping a set of samples before returning error
					UINT32 MaxSampleAttempts,				// allow at most this many attempts at randomly locating a sample before restarting the bootstrap
					bool bTargs,						// false: generate bootstrap sampling from query assembly sequences, true: bootstrap sampling from target assembling sequences
					bool bWithoutReplacement,			// false: sampling with replacement, true: sampling without replacement (currently not implemented)
					bool bNonOverlapping)				// false: samples may be overlapping, true: samples must be non-overlapping (currently not implemented)
{
UINT8 Base;
UINT32 CurSampleAttempt;
UINT32 CurBootstrapAttempt;
UINT32 SampleIdx;
UINT32 Len;
INT64 PopSeqOfs;
UINT8 *pSeq;
tsSeqBlock *pSeqBlock;
tsSeqAllocs *pPopulation;
tsSeqAllocs *pSample;

if(MaxBootstrapAttempts < 0)
	MaxBootstrapAttempts = 1;
if(MaxSampleAttempts < 0)
	MaxSampleAttempts = 1;

if(bTargs)
	{
	pPopulation = &m_Seqs[ePMBSSTargAssemb];
	pSample = &m_Seqs[ePMBSSTargSeqs];
	}
else
	{
	pPopulation = &m_Seqs[ePMBSSQueryAssemb];
	pSample = &m_Seqs[ePMBSSQuerySeqs];
	}

CurBootstrapAttempt = 0;
do {
	CurBootstrapAttempt += 1;
	// needing to clear any sequence constraint flags from a previous sampling bootstrapping iteration?
	if(bWithoutReplacement || bNonOverlapping)
		{
		pSeq = pPopulation->pSeqs;
		for(Len = 0; Len < pPopulation->UsedSeqsSize-1; Len++, pSeq++)
			*pSeq &= 0x0f;
		}

	pSeqBlock = pSample->pSeqBlocks;
	for(SampleIdx = 0; SampleIdx < pSample->UsedSeqBlocks; SampleIdx++,pSeqBlock++)
		{
		for(CurSampleAttempt = 0; CurSampleAttempt < MaxSampleAttempts; CurSampleAttempt++)
			{
			PopSeqOfs = GenRand60(pPopulation->UsedSeqsSize - 1);
			pSeq = &pPopulation->pSeqs[PopSeqOfs];
			for(Len = 0; Len < pSeqBlock->SeqLen; Len++, pSeq++)
				{
				Base = *pSeq;
				if((Base & 0x0f) > eBaseT || (bWithoutReplacement && Base & 0x10) || (bNonOverlapping && Base & 0x30))
					break;
				}
			if(Len == pSeqBlock->SeqLen)
				{
				pSeq = &pPopulation->pSeqs[PopSeqOfs];
				if(bWithoutReplacement)
					*pSeq |= 0x10;
				if(bNonOverlapping)
					{
					for(Len = 0; Len < pSeqBlock->SeqLen; Len++, pSeq++)
						*pSeq |= 0x20;
					}
				break;
				}
			}
		if(CurSampleAttempt == MaxSampleAttempts)
			break;
		// have an accepted sample which is same length as original
		pSeqBlock->PopSeqOfs = PopSeqOfs;
		}
	if(SampleIdx == pSample->UsedSeqBlocks)	// if able to bootstrap all samples then success!!!
		{
		if(bWithoutReplacement || bNonOverlapping)
			{
			pSeq = pPopulation->pSeqs;
			for(Len = 0; Len < pPopulation->UsedSeqsSize-1; Len++, pSeq++)  // ensure subsequent alignments won't fail because sequence constraint flags were not reset!
				*pSeq &= 0x0f;
			}
		return(SampleIdx);
		}
	}
while(CurBootstrapAttempt < MaxBootstrapAttempts);
return(-1); // failure to bootstrap
}


int			// index, -1 if no matches, at which Query matched onto target with at most MaxSubs
CAlignsBootstrap::AlignQueriesToTargs(bool bSenseOnly,	// true if to align sense only, default is to align both sense and antisense
						UINT32 StartQuerySeqIdx,		// alignments starting with this query sequence
						UINT32 EndQuerySeqIdx,			// through to this query sequence inclusive
						int MaxSubs)				// accepting at most this percentage of bases of query length to be mismatches
{
UINT8 RevCplQBases[cMaxQuerySeqLen+1];
tsSeqAllocs *pQSeqs;
tsSeqAllocs *pQAssemb;	
tsSeqBlock *pQBlock;
tsSeqAllocs *pTSeqs;
tsSeqAllocs *pTAssemb;	
tsSeqBlock *pTBlock;
UINT8 *pQBases;
UINT8 *pQWBase;
int QLen;
UINT8 *pTBases;
UINT8 *pTWBases;
UINT8 *pTWBase;
int TLen;
int BaseIdx;
int WinIdx;
int MMCnt;
int MaxMMCnt;
UINT32 CurQueryIdx;
UINT32 CurTargIdx;
int MaxNumPasses;
int CurPass;
int NumQueryHits;
tsQueryHit *pQueryHits;

pQAssemb = &m_Seqs[ePMBSSQueryAssemb];
pQSeqs = &m_Seqs[ePMBSSQuerySeqs];
pTAssemb = &m_Seqs[ePMBSSTargAssemb];
pTSeqs = &m_Seqs[ePMBSSTargSeqs];

NumQueryHits = 0;
MaxNumPasses = bSenseOnly ? 1 : 2; // either 1 (bSenseOnly true) or 2 passes (bSenseOnly false); first pass will be with query sense, if two passes then second pass will be with query antisense

for(CurTargIdx = 0; CurTargIdx < pTSeqs->NumSeqs; CurTargIdx++)
	{
	pTBlock = &pTSeqs->pSeqBlocks[CurTargIdx];
	TLen = pTBlock->SeqLen;
	if(m_bUseTargBS == false)
		pTBases = &pTSeqs->pSeqs[pTBlock->SmplSeqOfs];
	else
		pTBases = &pTAssemb->pSeqs[pTBlock->PopSeqOfs];
	CurPass = 1;
	do {
		pQueryHits = &m_pQueryHits[StartQuerySeqIdx];
		for(CurQueryIdx = StartQuerySeqIdx; CurQueryIdx <= EndQuerySeqIdx; CurQueryIdx++,pQueryHits++)
			{	
			if(pQueryHits->flgHit == 1 && pQueryHits->MMCnt == 0)
				continue;
			pQBlock = &pQSeqs->pSeqBlocks[CurQueryIdx];
			QLen = pQBlock->SeqLen;
			if(TLen < QLen)
				continue;

			if(m_bUseQueryBS == false)
				pQBases = &pQSeqs->pSeqs[pQBlock->SmplSeqOfs];
			else
				pQBases = &pQAssemb->pSeqs[pQBlock->PopSeqOfs];

			if(CurPass == 2)  // a second pass only if first pass for sense completed and bSenseOnly was false
				{
				memcpy(RevCplQBases,pQBases,QLen);
				CSeqTrans::ReverseComplement(QLen,RevCplQBases);
				pQBases = RevCplQBases;
				}

			MaxMMCnt = (QLen * MaxSubs) / 100;
			if(pQueryHits->flgHit && pQueryHits->MMCnt < MaxMMCnt)
				MaxMMCnt = pQueryHits->MMCnt;	
			pTWBases = pTBases;
			for(WinIdx = 0; WinIdx < (TLen - QLen); WinIdx++,pTWBases++)
				{
				pTWBase = pTWBases;
				pQWBase = pQBases;
				MMCnt = 0;
				for(BaseIdx = 0; BaseIdx < QLen; BaseIdx++,pQWBase++,pTWBase++)
					{
					if(*pQWBase != *pTWBase)
						{
						MMCnt += 1;
						if((pQueryHits->flgHit && MMCnt > pQueryHits->MMCnt) || MMCnt > MaxMMCnt)
							break;
						}
					}
				if(BaseIdx == QLen && (pQueryHits->flgHit == 0 || MMCnt < pQueryHits->MMCnt))
					{
					if(pQueryHits->flgHit == 0)
						{
						NumQueryHits += 1;
						pQueryHits->flgHit = 1;
						}
					pQueryHits->MMCnt = MMCnt;
					pQueryHits->TargOfs = WinIdx;
					pQueryHits->TargIdx = CurTargIdx;
					pQueryHits->flgAntisense = CurPass == 2 ? 1 : 0;
					}
				if(pQueryHits->flgHit == 1 && pQueryHits->MMCnt == 0)
					break;
				}
			}
		}
	while(++CurPass <= MaxNumPasses);
	}
return(NumQueryHits);
}

int
CAlignsBootstrap::AlignBootstrap(int NumRepeats)  // number of times current set of counts are to be reported
{
int RepeatIdx;
UINT32 NumQueryHits;
UINT32 NumTargHits;
UINT32 QueryHitIdx;
UINT32 TargHitIdx;
tsQueryHit *pQueryHit;
tsSeqBlock *pTargHit;
tsSeqAllocs *pTargSeqs;
memset(m_pQueryHits,0,sizeof(tsQueryHit) * m_Seqs[ePMBSSQuerySeqs].NumSeqs);
StartAlignments();
while(!WaitAlignments(60))
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Process: Aligning ...");
		
pTargSeqs = &m_Seqs[ePMBSSTargSeqs];
pTargHit = pTargSeqs->pSeqBlocks;
for(TargHitIdx = 0; TargHitIdx < pTargSeqs->NumSeqs; TargHitIdx++, pTargHit++)
	pTargHit->NumQueryHits = 0;	

NumQueryHits = 0;
pQueryHit = m_pQueryHits;
for(QueryHitIdx = 0; QueryHitIdx < m_Seqs[ePMBSSQuerySeqs].NumSeqs; QueryHitIdx++, pQueryHit++)
	{
	if(pQueryHit->flgHit)
		{
		NumQueryHits += 1;
		pTargHit = &pTargSeqs->pSeqBlocks[pQueryHit->TargIdx];
		pTargHit->NumQueryHits += 1;
		}
	}

NumTargHits = 0;
pTargHit = pTargSeqs->pSeqBlocks;
for(TargHitIdx = 0; TargHitIdx < m_Seqs[ePMBSSTargSeqs].NumSeqs; TargHitIdx++, pTargHit++)
	if(pTargHit->NumQueryHits > 0)
		NumTargHits += 1;

for(RepeatIdx = 0; RepeatIdx < NumRepeats; RepeatIdx++)
	{
	if(m_CurQRsltsOfs + 100 > m_AllocRsltsBuff)
		{
		CUtility::SafeWrite(m_hCSVQRslts,m_pszQRsltsBuff,m_CurQRsltsOfs);
		m_CurQRsltsOfs = 0;
		}
	m_CurQRsltsOfs += sprintf(&m_pszQRsltsBuff[m_CurQRsltsOfs],",%d",NumQueryHits);
	}


for(RepeatIdx = 0; RepeatIdx < NumRepeats; RepeatIdx++)
	{
	if(m_CurTRsltsOfs + 100 > m_AllocRsltsBuff)
		{
		CUtility::SafeWrite(m_hCSVTRslts,m_pszTRsltsBuff,m_CurTRsltsOfs);
		m_CurTRsltsOfs = 0;
		}
	m_CurTRsltsOfs += sprintf(&m_pszTRsltsBuff[m_CurTRsltsOfs],",%d",NumTargHits);
	}

return(0);
}

int
CAlignsBootstrap::ReportBootstrapSeqs(bool bTargSeqs,		// true if target sequences to be reported 
					int Iteration,		// which bootstrap iteration (1..n)
					char *pszSeqsFile)  // write bootstraps into this file, will have bootstrap iteration specific suffix appended
{
int hFile;
int RsltsFileLen;
char szSeqsFile[_MAX_PATH];
tsSeqAllocs *pSeqs;
tsSeqAllocs *pAssemb;	
UINT32 SeqLen;
UINT8 *pBases;
tsSeqBlock *pBlock;

UINT8 *pSeqBuff;
char szLineBuff[100];
int CurSeqBuffIdx;
UINT32 CurIdx;

pSeqBuff = new UINT8 [cMaxTargSeqLen+1];

hFile = -1;
CurSeqBuffIdx = 0;
RsltsFileLen = (int)strlen(pszSeqsFile);
strcpy(szSeqsFile,pszSeqsFile);
if(stricmp(&szSeqsFile[RsltsFileLen-4],".csv")==0)
	RsltsFileLen -= 4;
sprintf(&szSeqsFile[RsltsFileLen],".%cbs.%d.fasta",bTargSeqs ? 't' : 'q',Iteration);

// all sequences loaded with no errors, worth now creating the summary result files
#ifdef _WIN32
if((hFile = open(szSeqsFile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
if((hFile = open(szSeqsFile, O_RDWR | O_CREAT | O_TRUNC, S_IREAD | S_IWRITE ))==-1)
#endif
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create or truncate sampled query sequences file %s error: %s",szSeqsFile,strerror(errno));
	delete pSeqBuff;
	Reset();
	return(eBSFerrCreateFile);
	}

if(bTargSeqs)
	{
	pAssemb = &m_Seqs[ePMBSSTargAssemb];
	pSeqs = &m_Seqs[ePMBSSTargSeqs];
	}
else
	{
	pAssemb = &m_Seqs[ePMBSSQueryAssemb];
	pSeqs = &m_Seqs[ePMBSSQuerySeqs];
	}

for(CurIdx = 0; CurIdx < pSeqs->NumSeqs; CurIdx++)
	{
	pBlock = &pSeqs->pSeqBlocks[CurIdx];
	SeqLen = pBlock->SeqLen;
	pBases = &pAssemb->pSeqs[pBlock->PopSeqOfs];
	CurSeqBuffIdx+=sprintf((char *)&pSeqBuff[CurSeqBuffIdx],">bsseq%d\n",CurIdx+1);
	while(SeqLen > 0)
		{
		CSeqTrans::MapSeq2Ascii(pBases,min(SeqLen,70),szLineBuff);
		pBases += 70;
		CurSeqBuffIdx += sprintf((char *)&pSeqBuff[CurSeqBuffIdx],"%s\n",szLineBuff);
		SeqLen -= min(SeqLen,70);
		if((CurSeqBuffIdx + 1000) > cMaxTargSeqLen)
			{
			CUtility::SafeWrite(hFile,pSeqBuff,CurSeqBuffIdx);
			CurSeqBuffIdx = 0;
			}
		}
	}

if(hFile != -1 && CurSeqBuffIdx > 0)
	{
	CUtility::SafeWrite(hFile,pSeqBuff,CurSeqBuffIdx);
	CurSeqBuffIdx = 0;
	}
if(hFile != -1)
	{
#ifdef _WIN32
	_commit(hFile);
#else
	fsync(hFile);
#endif
	close(hFile);
	hFile = -1;
	}
if(pSeqBuff != NULL)
	delete pSeqBuff;
return(0);
}


int 
CAlignsBootstrap::Process(ePMBSAlign PMode,	// bootstrap processing mode
				int RandSeed,				// if > 0 then random generator seed, otherwise current time used as the seed
				bool bSenseOnly,			// true if to align only sense, default is to align both sense and antisense
				int MaxSubs,				// allowing at most this many subs as percentage of query length before accepting alignment
				int NumBootstraps,			// number of bootstrap iterations, excludes initial original query sequences aligned onto initial target sequences
				bool bWORreplacement,		// sample without replacement, default is to sample with replacement
				bool bNoOverlaps,			// sample without overlap, default is to sample allowing overlaps
				char *pszQuerySeqsFile,		// fasta file containing initial query sequences from which to derive query length distributions
				char *pszTargSeqsFile,		// fasta file containing initial target sequences from which to derive target length distributions
				char *pszQueryAssembFile,	// file containing fasta assembly to be bootstrap sampled for query sequences with same length distributions as sequences in pszQuerySeqsFile 
				char *pszTargAssembFile,	// file containing fasta assembly to be bootstrap sampled for target sequences with same length distributions as sequences in pszTargSeqsFile
				char *pszQRsltsFile,		// summary number of query hits onto at least one target bootstrap results to this file 
				char *pszTRsltsFile,		// summary number of targets hit by at least one query bootstrap results to this file
				int NumThreads)				// number of worker threads to use
{
int Rslt;

int CurIterQueries;

size_t memreq;

if((Rslt = Init(RandSeed))< eBSFSuccess)
	return(Rslt);

m_PMode = PMode;
m_RandSeed = RandSeed;
m_bSenseOnly = bSenseOnly;
m_MaxSubs = MaxSubs;
m_NumBootstraps = NumBootstraps;
m_NumThreads = NumThreads;
m_bWORreplacement = bWORreplacement;
m_bNoOverlaps = bNoOverlaps;

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Process: Loading query sequences from '%s'",pszQuerySeqsFile);
if((Rslt = LoadFastaSeqs(ePMBSSQuerySeqs,10,pszQuerySeqsFile)) < eBSFSuccess)
	{
	Reset();
	return(Rslt);
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Process: Loading target sequences from '%s'",pszTargSeqsFile);
if((Rslt = LoadFastaSeqs(ePMBSSTargSeqs,10,pszTargSeqsFile)) < eBSFSuccess)
	{
	Reset();
	return(Rslt);
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Process: Loading query assembly from '%s'",pszQueryAssembFile);
if((Rslt = LoadFastaSeqs(ePMBSSQueryAssemb,10,pszQueryAssembFile)) < eBSFSuccess)
	{
	Reset();
	return(Rslt);
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Process: Loading target assembly from '%s'",pszTargAssembFile);
if((Rslt = LoadFastaSeqs(ePMBSSTargAssemb,10,pszTargAssembFile)) < eBSFSuccess)
	{
	Reset();
	return(Rslt);
	}

memreq = sizeof(tsQueryHit) * m_Seqs[ePMBSSQuerySeqs].NumSeqs;
if((m_pQueryHits = (tsQueryHit *)malloc(memreq)) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: Unable to allocate memory for %u query hits",m_Seqs[ePMBSSQuerySeqs].NumSeqs);
	Reset();
	return(eBSFerrMem);
	}

m_AllocRsltsBuff = 50000;
m_CurQRsltsOfs = 0;
m_CurTRsltsOfs = 0;

memreq = (size_t)m_AllocRsltsBuff;		
if((m_pszQRsltsBuff = (char *)malloc(memreq)) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: Unable to allocate memory for %u query hit summary results",m_AllocRsltsBuff);
	Reset();
	return(eBSFerrMem);
	}

if((m_pszTRsltsBuff = (char *)malloc(memreq)) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: Unable to allocate memory for %u target hit summary results",m_AllocRsltsBuff);
	Reset();
	return(eBSFerrMem);
	}

// all sequences loaded with no errors, worth now creating the summary result files
#ifdef _WIN32
if((m_hCSVQRslts = open(pszQRsltsFile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
if((m_hCSVQRslts = open(pszQRsltsFile, O_RDWR | O_CREAT | O_TRUNC, S_IREAD | S_IWRITE ))==-1)
#endif
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create or truncate results file %s error: %s",pszQRsltsFile,strerror(errno));
	Reset();
	return(eBSFerrCreateFile);
	}
#ifdef _WIN32
if((m_hCSVTRslts = open(pszTRsltsFile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
if((m_hCSVTRslts = open(pszTRsltsFile, O_RDWR | O_CREAT | O_TRUNC, S_IREAD | S_IWRITE ))==-1)
#endif
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create or truncate results file %s error: %s",pszTRsltsFile,strerror(errno));
	Reset();
	return(eBSFerrCreateFile);
	}

m_CurQRsltsOfs = sprintf(m_pszQRsltsBuff,"\"QSets\"");
for(CurIterQueries = 1; CurIterQueries <= NumBootstraps; CurIterQueries++)
	{
	m_CurQRsltsOfs += sprintf(&m_pszQRsltsBuff[m_CurQRsltsOfs],",\"BS:%d\"",CurIterQueries);
	if(m_CurQRsltsOfs + 100 > m_AllocRsltsBuff)
		{
		CUtility::SafeWrite(m_hCSVQRslts,m_pszQRsltsBuff,m_CurQRsltsOfs);
		m_CurQRsltsOfs = 0;
		}
	}

if(m_CurQRsltsOfs > 0)
	{
	CUtility::SafeWrite(m_hCSVQRslts,m_pszQRsltsBuff,m_CurQRsltsOfs);
	m_CurQRsltsOfs = 0;
	}

m_CurTRsltsOfs = sprintf(m_pszTRsltsBuff,"\"TSets\"");
for(CurIterQueries = 1; CurIterQueries <= NumBootstraps; CurIterQueries++)
	{
	m_CurTRsltsOfs += sprintf(&m_pszTRsltsBuff[m_CurTRsltsOfs],",\"BS:%d\"",CurIterQueries);
	if(m_CurTRsltsOfs + 100 > m_AllocRsltsBuff)
		{
		CUtility::SafeWrite(m_hCSVTRslts,m_pszTRsltsBuff,m_CurTRsltsOfs);
		m_CurTRsltsOfs = 0;
		}
	}

if(m_CurTRsltsOfs > 0)
	{
	CUtility::SafeWrite(m_hCSVTRslts,m_pszTRsltsBuff,m_CurTRsltsOfs);
	m_CurTRsltsOfs = 0;
	}

#ifdef _DEBUG
#ifdef _WIN32
_ASSERTE( _CrtCheckMemory());
#endif
#endif

// create pool of worker threads
if(m_NumThreads > (int)m_Seqs[ePMBSSQuerySeqs].NumSeqs)
	m_NumThreads = (int)m_Seqs[ePMBSSQuerySeqs].NumSeqs;
StartWorkerThreads(m_NumThreads,m_Seqs[ePMBSSQuerySeqs].NumSeqs);

// firstly process Class1 - this is a pseudo bootstrap which is doing alignment of original query sequences onto original target sequences but reporting as though bootstraps
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Process: Aligning original query sequences onto original target sequences ...");
m_CurQRsltsOfs += sprintf(&m_pszQRsltsBuff[m_CurQRsltsOfs],"\n\"QalignT\"");
m_CurTRsltsOfs += sprintf(&m_pszTRsltsBuff[m_CurTRsltsOfs],"\n\"QalignT\"");
m_bUseTargBS = false;
m_bUseQueryBS = false;
AlignBootstrap(m_NumBootstraps);

// next process Class2 - query sequences are bootstrapped and aligned to original targets
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Process: Aligning bootstrapped query sequences onto original target sequences ...");
m_CurQRsltsOfs += sprintf(&m_pszQRsltsBuff[m_CurQRsltsOfs],"\n\"QBSalignT\"");
m_CurTRsltsOfs += sprintf(&m_pszTRsltsBuff[m_CurTRsltsOfs],"\n\"QBSalignT\"");
m_bUseTargBS = false;
m_bUseQueryBS = true;
for(m_CurBootstrap = 0; m_CurBootstrap < m_NumBootstraps; m_CurBootstrap++)
	{
	if((Rslt = GenBootstrap(cDfltBootstrappingAttempts,cDfltSamplingAttempts,false,m_bWORreplacement,m_bNoOverlaps)) < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: Bootstrapping query sequences failed");
		Reset();
		return(Rslt);
		}
	if(m_PMode == ePMSAreportseqs)
		{
		// write bootstrap query sequences to file
		ReportBootstrapSeqs(false,m_CurBootstrap+1,pszQRsltsFile); 
		}
	AlignBootstrap(1);
	}

// next process Class3 - original query sequences are aligned to bootstrapped target sequences
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Process: Aligning original query sequences onto bootstrapped target sequences ...");
m_CurQRsltsOfs += sprintf(&m_pszQRsltsBuff[m_CurQRsltsOfs],"\n\"QalignTBS\"");
m_CurTRsltsOfs += sprintf(&m_pszTRsltsBuff[m_CurTRsltsOfs],"\n\"QalignTBS\"");
m_bUseTargBS = true;
m_bUseQueryBS = false;
for(m_CurBootstrap = 0; m_CurBootstrap < m_NumBootstraps; m_CurBootstrap++)
	{
	if((Rslt = GenBootstrap(cDfltBootstrappingAttempts,cDfltSamplingAttempts,true,m_bWORreplacement,m_bNoOverlaps)) < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: Bootstrapping target sequences failed");
		Reset();
		return(Rslt);
		}
	if(m_PMode == ePMSAreportseqs)
		{
		// write bootstrap target sequences to file
		ReportBootstrapSeqs(true,m_CurBootstrap+1,pszTRsltsFile);
		}
	AlignBootstrap(1);
	}

// finally process Class 4 - bootstrapped query sequences are aligned to bootstrapped target sequences
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Process: Aligning bootstrapped query sequences onto bootstrapped target sequences ...");
m_CurQRsltsOfs += sprintf(&m_pszQRsltsBuff[m_CurQRsltsOfs],"\n\"QBSalignTBS\"");
m_CurTRsltsOfs += sprintf(&m_pszTRsltsBuff[m_CurTRsltsOfs],"\n\"QBSalignTBS\"");
m_bUseTargBS = true;
m_bUseQueryBS = true;
for(m_CurBootstrap = 0; m_CurBootstrap < m_NumBootstraps; m_CurBootstrap++)
	{
	if((Rslt = GenBootstrap(cDfltBootstrappingAttempts,cDfltSamplingAttempts,false,m_bWORreplacement,m_bNoOverlaps)) < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: Bootstrapping query sequences failed");
		Reset();
		return(Rslt);
		}

	if((Rslt = GenBootstrap(cDfltBootstrappingAttempts,cDfltSamplingAttempts,true,m_bWORreplacement,m_bNoOverlaps)) < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: Bootstrapping target sequences failed");
		Reset();
		return(Rslt);
		}

	AlignBootstrap(1);
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Process: Bootstrapping and alignments completed");

#ifdef _DEBUG
#ifdef _WIN32
_ASSERTE( _CrtCheckMemory());
#endif
#endif

if(m_hCSVQRslts != -1 && m_CurQRsltsOfs > 0)
	{
	CUtility::SafeWrite(m_hCSVQRslts,m_pszQRsltsBuff,m_CurQRsltsOfs);
	m_CurQRsltsOfs = 0;
	}
if(m_hCSVQRslts != -1)
	{
#ifdef _WIN32
	_commit(m_hCSVQRslts);
#else
	fsync(m_hCSVQRslts);
#endif
	close(m_hCSVQRslts);
	m_hCSVQRslts = -1;
	}
if(m_hCSVTRslts != -1 && m_CurTRsltsOfs > 0)
	{
	CUtility::SafeWrite(m_hCSVTRslts,m_pszTRsltsBuff,m_CurTRsltsOfs);
	m_CurTRsltsOfs = 0;
	}
if(m_hCSVTRslts != -1)
	{
#ifdef _WIN32
	_commit(m_hCSVTRslts);
#else
	fsync(m_hCSVTRslts);
#endif
	close(m_hCSVTRslts);
	m_hCSVTRslts = -1;
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Process: Bootstrapping completed");
TerminateWorkerThreads();
#ifdef _DEBUG
#ifdef _WIN32
_ASSERTE( _CrtCheckMemory());
#endif
#endif
Reset();
return(Rslt);
} 


#ifdef _WIN32
unsigned __stdcall WorkerInstance(void * pThreadPars)
#else
void *WorkerInstance(void * pThreadPars)
#endif
{
int Rslt;
tsWorkerInstance *pPars = (tsWorkerInstance *)pThreadPars;			// makes it easier not having to deal with casts!
CAlignsBootstrap *pWorkerInstance = (CAlignsBootstrap *)pPars->pThis;

Rslt = pWorkerInstance->ProcWorkerThread(pPars);
pPars->Rslt = Rslt;
#ifdef _WIN32
_endthreadex(0);
return(eBSFSuccess);
#else
pthread_exit(NULL);
#endif
}

// initialise and start pool of worker threads
int
CAlignsBootstrap::StartWorkerThreads(UINT32 NumThreads,		// there are this many threads in pool
									UINT32 NumQuerySeqs)	// which will be processing a total of this many query sequences
{
UINT32 MaxWait;
UINT32 StartQuerySeqIdx;
UINT32 ThreadIdx;
UINT32 StartedInstances;
tsWorkerInstance *pThreadPar;
m_TermAllThreads = 0;

if(NumThreads > NumQuerySeqs)
	NumThreads = NumQuerySeqs;

StartQuerySeqIdx = 0;
pThreadPar = m_WorkerInstances;
for(ThreadIdx = 0; ThreadIdx < NumThreads; ThreadIdx++,pThreadPar++)
	{
	memset(pThreadPar,0,sizeof(tsWorkerInstance));
#ifdef _WIN32
	pThreadPar->threadHandle = NULL;
#else
	pThreadPar->threadID = 0;
#endif
	pThreadPar->AlignReqID = 0;
	pThreadPar->StartQuerySeqIdx = StartQuerySeqIdx;		// alignments starting with this query sequence
	pThreadPar->EndQuerySeqIdx = (UINT32)((((float)NumQuerySeqs * (ThreadIdx+1))/(float)NumThreads)-1.0);			// through to this query sequence inclusive
	StartQuerySeqIdx = pThreadPar->EndQuerySeqIdx+1;
	}


m_NumWorkerInsts = 0;
m_CompletedWorkerInsts = 0;
pThreadPar = m_WorkerInstances;
for (ThreadIdx = 1; ThreadIdx <= NumThreads; ThreadIdx++, pThreadPar++)
	{
	pThreadPar->ThreadIdx = ThreadIdx;
	pThreadPar->pThis = this;
#ifdef _WIN32
	pThreadPar->threadHandle = (HANDLE)_beginthreadex(NULL, 0x0fffff, WorkerInstance, pThreadPar, 0, &pThreadPar->threadID);
#else
	pThreadPar->threadRslt = pthread_create(&pThreadPar->threadID, NULL, WorkerInstance, pThreadPar);
#endif
	}

// allow threads a few seconds to startup
// check if all threads did actually startup; if they did so then m_NumWorkerInsts will have been incremented to NumInstances
MaxWait = 60;		// allowing at most 60 secs for threads to startup
do {
#ifdef WIN32
	Sleep(1000);
	StartedInstances = InterlockedCompareExchange(&m_NumWorkerInsts,NumThreads,NumThreads);
#else
	sleep(1);
	StartedInstances = __sync_val_compare_and_swap (&m_NumWorkerInsts,NumThreads,NumThreads);
#endif
	MaxWait -= 1;
	}
while(StartedInstances != NumThreads && MaxWait > 0);
if(StartedInstances != NumThreads)
	{
	TerminateWorkerThreads();
	StartedInstances = 0;
	}
return(StartedInstances);
}


bool // true if any worker threads in pool to start alignments, false if no worker threads 
CAlignsBootstrap::StartAlignments(void) // signal worker pool of threads that there is a new bootstrap set to be aligned
{
UINT32 StartedInstances;

#ifdef WIN32
StartedInstances = InterlockedCompareExchange(&m_NumWorkerInsts,0,0);
#else
StartedInstances = __sync_val_compare_and_swap (&m_NumWorkerInsts,0,0);
#endif
if(StartedInstances == 0)
	return(false);

m_CompletedWorkerInsts = 0;

#ifdef WIN32
	InterlockedIncrement(&m_AlignReqID);
#else
	__sync_fetch_and_add(&m_AlignReqID,1);
#endif
return(true);
}


bool
CAlignsBootstrap::WaitAlignments(int WaitSecs)	// allow at most this many seconds for pool of worker threads to complete aligning current bootstrap set
{
INT64 WaitMS;
WaitMS = (INT64)WaitSecs * 1000;

UINT32 CompletedInstances;
do {
	CUtility::SleepMillisecs(20);
#ifdef WIN32
	CompletedInstances = InterlockedCompareExchange(&m_CompletedWorkerInsts,m_NumWorkerInsts,m_NumWorkerInsts);
#else
	CompletedInstances = __sync_val_compare_and_swap (&m_CompletedWorkerInsts,m_NumWorkerInsts,m_NumWorkerInsts);
#endif
	WaitMS -= 20;
	}
while(CompletedInstances != m_NumWorkerInsts && WaitMS > 0);

return(CompletedInstances != m_NumWorkerInsts ? false : true);
}

// stop all threads in worker pool
int
CAlignsBootstrap::TerminateWorkerThreads(int WaitSecs)				// alow at most this many seconds before force terminating worker pool threads
{
int NumForceTerminated;
UINT32 Idx;
UINT32 StartedInstances; 
tsWorkerInstance *pThreadPar;
time_t Then;
time_t Now;

#ifdef WIN32
StartedInstances = InterlockedCompareExchange(&m_NumWorkerInsts,0,0);
#else
StartedInstances = __sync_val_compare_and_swap (&m_NumWorkerInsts,0,0);
#endif
if(StartedInstances == 0)
	{
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "TerminateWorkerThreads: No worker threads to terminate");
	return(0);
	}

// request all worker threads to self terminate
#ifdef WIN32
InterlockedCompareExchange(&m_TermAllThreads,1,0);
#else
__sync_val_compare_and_swap (&m_TermAllThreads,0,1);
#endif
gDiagnostics.DiagOut(eDLInfo, gszProcName, "TerminateWorkerThreads: Requesting %u worker threads to terminate",StartedInstances);
Then = time(NULL) + WaitSecs;
NumForceTerminated = 0;
pThreadPar = m_WorkerInstances;
for(Idx = 0; Idx < StartedInstances; Idx++, pThreadPar += 1)
	{
	Now = time(NULL);
	if(Now >= Then)
		Now = 1;
	else
		Now = Then - Now;

#ifdef WIN32
	if(pThreadPar->threadHandle != NULL)
		{
		if(WAIT_TIMEOUT == WaitForSingleObject(pThreadPar->threadHandle, (UINT32)Now * 1000))
			{
			NumForceTerminated += 1;
			TerminateThread(pThreadPar->threadHandle,0);
			}
		pThreadPar->threadHandle = NULL;
		}
#else
	if(pThreadPar->threadID != 0)
		{
		struct timespec ts;
		int JoinRlt;
		void *pExitRslt;
		clock_gettime(CLOCK_REALTIME, &ts);
		ts.tv_sec += Now;
		if ((JoinRlt = pthread_timedjoin_np(pThreadPar->threadID, &pExitRslt, &ts)) != 0)
			{
			gDiagnostics.DiagOut(eDLInfo, gszProcName, "TerminateWorkerThreads: Force terminating thread %u, pthread_timedjoin_np() returned %d",pThreadPar->ThreadIdx,JoinRlt);
			NumForceTerminated += 1;
			pthread_cancel(pThreadPar->threadID);	
			pthread_join(pThreadPar->threadID, NULL);
			}
		pThreadPar->threadID = 0;
		}
#endif
	}

pThreadPar = m_WorkerInstances;
for(Idx = 0; Idx < StartedInstances; Idx++, pThreadPar += 1)
	memset(pThreadPar,0,sizeof(tsWorkerInstance));

m_TermAllThreads = 0;	
return(NumForceTerminated);
}

int
CAlignsBootstrap::ProcWorkerThread(tsWorkerInstance *pThreadPar)	// worker thread parameters
{
UINT32 AlignReqID;

// this thread has started, one more worker thread
#ifdef WIN32
InterlockedIncrement(&m_NumWorkerInsts);
#else
__sync_fetch_and_add(&m_NumWorkerInsts,1);
#endif

while(1) {
	// check if requested to terminate
#ifdef WIN32
	if(InterlockedCompareExchange(&m_TermAllThreads,1,1)==1)
#else
	if(__sync_val_compare_and_swap (&m_TermAllThreads,1,1)==1)
#endif
		break;

	// any work to do? m_AlignReqID will have been incremented by controlling thread when new set of alignments required
#ifdef WIN32
	AlignReqID = InterlockedCompareExchange(&m_AlignReqID,0,0);
#else
	AlignReqID = __sync_val_compare_and_swap (&m_AlignReqID,0,0);
#endif
	if(AlignReqID == pThreadPar->AlignReqID)
		{
		CUtility::SleepMillisecs(20);
		continue;
		}

	pThreadPar->AlignReqID = AlignReqID;
	AlignQueriesToTargs(m_bSenseOnly,			// true if to align sense only, default is to align both sense and antisense
						pThreadPar->StartQuerySeqIdx,		// alignments starting with this query sequence
						pThreadPar->EndQuerySeqIdx,			// through to this query sequence inclusive
						m_MaxSubs);				// accepting at most this percentage of bases of query length to be mismatches

#ifdef WIN32
	InterlockedIncrement(&m_CompletedWorkerInsts);
#else
	__sync_fetch_and_add(&m_CompletedWorkerInsts,1);
#endif
	}

// terminating, one less worker thread
#ifdef WIN32
InterlockedDecrement(&m_NumWorkerInsts);
#else
__sync_fetch_and_sub(&m_NumWorkerInsts,1);
#endif
return(0);
}

