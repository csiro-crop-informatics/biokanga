// SSRdiscovery.cpp : Defines the entry point for the console application.
//
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

#include "SQLiteSummaries.h"
#include "SSRIdentify.h"

const char *cpszProgVer = "0.6.0";		// increment with each release

const int cMaxInFileSpecs = 20;			// allow at most this number of input files

int
Process(int PMode,					// processing mode
		teRptSSRsFromat RptSSRsFormat,	// report SSRs in this file format 
		int MinRepElLen,			// identify repeating elements of this minimum length
		int MaxRepElLen,			// ranging upto this maximum length
		int MinTandemRpts,			// minimum number of tandem repeats
		int MaxTandemRpts,			// maximum number of repeats
		int NumInputFiles,			// number of input filespecs
		char *pszInputFiles[],		// input multifasta files
		char *pszKMerFreqFile,		// optional, output element KMer freq to this file
		char *pszOutFile);			// output SSRs to this file

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
	return _T("SSRdiscovery");
}
// end of str library required code
#endif

CStopWatch gStopWatch;
CDiagnostics gDiagnostics;				// for writing diagnostics messages to log file
char gszProcName[_MAX_FNAME];			// process name

CSQLiteSummaries gSQLiteSummaries;		// for writing processing result summaries to SQLite database
int	gExperimentID = 0;					// SQLite experiment identifier
int gProcessID = 0;						// SQLite process identifier
int	gProcessingID = 0;					// SQLite processing identifier

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

int PMode;					// processing mode

int MinRepElLen;			// identify repeating elements of this minimum length
int MaxRepElLen;			// ranging upto this maximum length
int MinTandemRpts;			// minimum number of tandem repeats
int MaxTandemRpts;			// maximum number of repeats

int NumInputFiles;			// number of intput filespecs
char *pszInputFiles[cMaxInFileSpecs];	// names of input files (wildcards allowed) containing sequences to be processed for SSRs

char szKMerFreqFile[_MAX_PATH];	// optional, output element KMer freq to this file
char szOutFile[_MAX_PATH];		// write identified SSRs to this file

char szSQLiteDatabase[_MAX_PATH];	// results summaries to this SQLite file
char szExperimentName[cMaxDatasetSpeciesChrom+1];			// experiment name
char szExperimentDescr[1000];		// describes experiment

// command line args
struct arg_lit  *help    = arg_lit0("h","help",                 "Print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"Print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to screen and logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"Diagnostics log file");

struct arg_int *mode = arg_int0("m","mode","<int>",								"Processing mode: 0 default");

struct arg_int *minrepellen = arg_int0("k","minrepellen","<int>",				"identify repeating K-mer elements of this minimum length (2 default)");
struct arg_int *maxrepellen = arg_int0("K","maxrepellen","<int>",				"identify repeating K-mer elements of this maximum length (5 default)");
struct arg_int *mintandemrpts = arg_int0("r","mintandemrpts","<int>",			"minimum number of tandem element repeats (5 default)");
struct arg_int *maxtandemrpts = arg_int0("R","maxtandemrpts","<int>",			"maximum number of tandem element repeats (10 default)");

struct arg_file *inputfiles = arg_filen("i","in","<file>",0,cMaxInFileSpecs,	"Input file(s) containing sequences to process for SSRs");
struct arg_file *kmerfreq = arg_file0("O","outkmerfreq","<file>",				"Output K-mer element freq to this file");
struct arg_file *outfile = arg_file1("o","out","<file>",						"Output SSRs to this file");

struct arg_file *summrslts = arg_file0("q","sumrslts","<file>",					"Output results summary to this SQLite3 database file");
struct arg_str *experimentname = arg_str0("w","experimentname","<str>",			"experiment name SQLite3 database file");
struct arg_str *experimentdescr = arg_str0("W","experimentdescr","<str>",		"experiment description SQLite3 database file");

struct arg_end *end = arg_end(40);

void *argtable[] = {help,version,FileLogLevel,LogFile,
					summrslts,experimentname,experimentdescr,
	                mode,minrepellen, maxrepellen, mintandemrpts, maxtandemrpts, inputfiles, kmerfreq, outfile,
					end};

char **pAllArgs;
int argerrors;
argerrors = CUtility::arg_parsefromfile(argc,(char **)argv,&pAllArgs);
if(argerrors >= 0)
	argerrors = arg_parse(argerrors,pAllArgs,argtable);

/* special case: '--help' takes precedence over error reporting */
if (help->count > 0)
        {
		printf("\n%s PE1/PE2 Corelation, Version %s\nOptions ---\n", gszProcName,cpszProgVer);
        arg_print_syntax(stdout,argtable,"\n");
        arg_print_glossary(stdout,argtable,"  %-25s %s\n");
		printf("\nNote: Parameters can be entered into a parameter file, one parameter per line.");
		printf("\n      To invoke this parameter file then precede it's name with '@'");
		printf("\n      e.g. %s @myparams.txt\n",gszProcName);
		printf("\nPlease report any issues regarding usage of %s to stuart.stephen@csiro.au\n\n",gszProcName);
		return(1);
        }

    /* special case: '--version' takes precedence error reporting */
if (version->count > 0)
        {
		printf("\n%s Version %s\n",gszProcName,cpszProgVer);
		return(1);
        }

if (!argerrors)
	{
	if(FileLogLevel->count && !LogFile->count)
		{
		printf("\nError: FileLogLevel '-f%d' specified but no logfile '-F<logfile>\n'",FileLogLevel->ival[0]);
		return(1);
		}

	iScreenLogLevel = iFileLogLevel = FileLogLevel->count ? FileLogLevel->ival[0] : eDLInfo;
	if(iFileLogLevel < eDLNone || iFileLogLevel > eDLDebug)
		{
		printf("\nError: FileLogLevel '-l%d' specified outside of range %d..%d\n",iFileLogLevel,eDLNone,eDLDebug);
		return(1);
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
		return(1);
		}
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Process %s Version %s starting",gszProcName,cpszProgVer);

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
		gProcessID = gSQLiteSummaries.AddProcess((char *)gszProcName,(char *)gszProcName,(char *)szExperimentDescr);
		if(gProcessID < 1)
			return(1);
		gProcessingID = gSQLiteSummaries.StartProcessing(gExperimentID,gProcessID,(char *)cpszProgVer);
		if(gProcessingID < 1)
			return(1);
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Initialised SQLite database '%s' for SSR results summary",szSQLiteDatabase);
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"SQLite database experiment identifier for '%s' is %d",szExperimentName,gExperimentID);
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"SQLite database process identifier for '%s' is %d",(char *)gszProcName,gProcessID);
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"SQLite database processing instance identifier is %d",gProcessingID);
		}
	else
		{
		szSQLiteDatabase[0] = '\0';
		szExperimentDescr[0] = '\0';
		}

	PMode = mode->count ? mode->ival[0] : 0;
	if(PMode < 0 || PMode > 1)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Processing mode '-m%d' must be in range 0..0",PMode);
		return(1);
		}


	MinRepElLen = minrepellen->count ? minrepellen->ival[0] : cDfltMinRepElLen;
	if(MinRepElLen < cMinRepElLen || MinRepElLen > cMaxRepElLen)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Repeating K-mer element minimum length '-k%d' must be in range %d..%d",MinRepElLen,cMinRepElLen,cMaxRepElLen);
		return(1);
		}

	MaxRepElLen = maxrepellen->count ? maxrepellen->ival[0] : max(cDfltMaxRepElLen,MinRepElLen);
	if(MaxRepElLen < MinRepElLen || MaxRepElLen > cMaxRepElLen)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Repeating K-mer element maximum length '-K%d' must be in range %d..%d",MaxRepElLen,MinRepElLen,cMaxRepElLen);
		return(1);
		}

	MinTandemRpts = mintandemrpts->count ? mintandemrpts->ival[0] : cDfltMinTandemRpts;
	if(MinTandemRpts < cMinTandemRpts || MinTandemRpts > cMaxTandemRpts)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Repeating K-mer element tandem min repeats '-r%d' must be in range %d..%d",MinTandemRpts,cMinTandemRpts,cMaxTandemRpts);
		return(1);
		}

	MaxTandemRpts = maxtandemrpts->count ? maxtandemrpts->ival[0] : max(cDfltMaxTandemRpts,MinTandemRpts);
	if(MaxTandemRpts < MinTandemRpts || MaxTandemRpts > cMaxTandemRpts)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Repeating K-mer element tandem max repeats '-R%d' must be in range %d..%d",MaxTandemRpts,MinTandemRpts,cMaxTandemRpts);
		return(1);
		}

	int Idx;
	for(NumInputFiles=Idx=0;NumInputFiles < cMaxInFileSpecs && Idx < inputfiles->count; Idx++)
		{
		pszInputFiles[Idx] = NULL;
		if(pszInputFiles[NumInputFiles] == NULL)
			pszInputFiles[NumInputFiles] = new char [_MAX_PATH];
		strncpy(pszInputFiles[NumInputFiles],inputfiles->filename[Idx],_MAX_PATH);
		pszInputFiles[NumInputFiles][_MAX_PATH-1] = '\0';
		CUtility::TrimQuotedWhitespcExtd(pszInputFiles[NumInputFiles]);
		if(pszInputFiles[NumInputFiles][0] != '\0')
			NumInputFiles++;
		}

	if(!NumInputFiles)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: After removal of whitespace, no input file(s) specified with '-i<filespec>' option)\n");
		exit(1);
		}

	if(kmerfreq->count)
		{
		if(!(MinRepElLen == MaxRepElLen && MinRepElLen <= 10))
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: SSR element frequency counting not supported for range of element lengths or if element length > 10\n");
			exit(1);
			}
		strcpy(szKMerFreqFile,kmerfreq->filename[0]);
		CUtility::TrimQuotedWhitespcExtd(szKMerFreqFile);
		}
	else
		szKMerFreqFile[0] = '\0';

	strcpy(szOutFile,outfile->filename[0]);
	CUtility::TrimQuotedWhitespcExtd(szOutFile);


	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing parameters:");
	const char *pszDescr;
	switch(PMode) {
		case 0:
			pszDescr = "Identify SSRs";
			break;
		}

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Processing mode is : '%s'",pszDescr);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"identify repeating K-mer elements of this minimum length: %d",MinRepElLen);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"identify repeating K-mer elements of this maximum length: %d",MaxRepElLen);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"minimum number of tandem element repeats: %d",MinTandemRpts);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"maximum number of tandem element repeats: %d",MaxTandemRpts);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Experiment name : '%s'",szExperimentName);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Experiment description : '%s'",szExperimentDescr);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"SQLite database file: '%s'",szOutFile);

	for(Idx = 0; Idx < NumInputFiles; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Input file (%d): '%s'",Idx+1,pszInputFiles[Idx]);
	if(szKMerFreqFile[0] != '\0')
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"SSR element K-Mer freq file: '%s'",szKMerFreqFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"SSRs to file: '%s'",szOutFile);

	if(gExperimentID > 0)
		{
		int ParamID;
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szLogFile),"log",szLogFile);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,sizeof(PMode),"mode",&PMode);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,sizeof(MinRepElLen),"minrepellen",&MinRepElLen);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,sizeof(MaxRepElLen),"maxrepellen",&MaxRepElLen);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,sizeof(MinTandemRpts),"mintandemrpts",&MinTandemRpts);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,sizeof(MaxTandemRpts),"maxtandemrpts",&MaxTandemRpts);

	    ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,sizeof(NumInputFiles),"NumInputFiles",&NumInputFiles);
		for(Idx=0; Idx < NumInputFiles; Idx++)
			ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(pszInputFiles[Idx]),"in",pszInputFiles[Idx]);

		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szOutFile),"out",szOutFile);
		if(szKMerFreqFile[0] != '\0')
			ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szKMerFreqFile),"outkmerfreq",szKMerFreqFile);

		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szSQLiteDatabase),"sumrslts",szSQLiteDatabase);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szExperimentName),"experimentname",szExperimentName);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szExperimentDescr),"experimentdescr",szExperimentDescr);
		}


#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	gStopWatch.Start();

	Rslt = Process(PMode,eRFCsv,MinRepElLen,MaxRepElLen,MinTandemRpts,MaxTandemRpts,NumInputFiles,pszInputFiles,szKMerFreqFile,szOutFile);
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
    printf("\n%s, Version %s\n", gszProcName,cpszProgVer);
	arg_print_errors(stdout,end,gszProcName);
	arg_print_syntax(stdout,argtable,"\nUse '-h' to view option and parameter usage\n");
	exit(1);
	}
return 0;
}

int
Process(int PMode,					// processing mode
		teRptSSRsFromat RptSSRsFormat,	// report SSRs in this file format 
		int MinRepElLen,			// identify repeating elements of this minimum length
		int MaxRepElLen,			// ranging upto this maximum length
		int MinTandemRpts,			// minimum number of tandem repeats
		int MaxTandemRpts,			// maximum number of repeats
		int NumInputFiles,			// number of input filespecs
		char *pszInputFiles[],		// input multifasta files
		char *pszKMerFreqFile,		// optional, output element KMer freq to this file
		char *pszOutFile)			// output SSRs to this file
{
int Rslt;
CSSRIdentify CSSRIdentify;
Rslt = CSSRIdentify.Process(PMode,RptSSRsFormat,MinRepElLen,MaxRepElLen,MinTandemRpts,MaxTandemRpts,NumInputFiles,pszInputFiles,pszKMerFreqFile,pszOutFile);
return(Rslt);
}




