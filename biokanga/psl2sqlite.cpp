// psl2sqlite.cpp : Defines the entry point for the console application.
// Process BLAT psl output format into SQLite


#include "stdafx.h"
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#if _WIN32

#include "../libbiokanga/commhdrs.h"
#else
#include "../libbiokanga/commhdrs.h"
#endif

#include "biokanga.h"
#include "SQLitePSL.h"

int 
Process(int PMode,					// processing mode, 0 to delete any existing then create new SQLite, 1 to append to existing SQLite
		int MinIdentity,			// minimum required identity
		int MinScore,				// minimum required score
		int MinMatches,				// minimum required base matches
		char *pszDatabase,			// SQLite database file
		char *pszExprName,			// experiment name
		char *pszPSLFile,			// alignments were parsed from this BLAT generated PSL file
		char *pszQueryFile,			// Blat'd query sequences in this file
		char *pszTargetFile,		// against targeted sequences in this file
		char *pszExprDescr,			// describes experiment
		char *pszBlatParams);		// Blat parameters used

int TrimQuotes(char *pszTxt);

#ifdef _WIN32
int PSL2SQLite(int argc, char* argv[])
{
// determine my process name
_splitpath(argv[0],NULL,NULL,gszProcName,NULL);
#else
int
PSL2SQLite(int argc, char** argv)
{
// determine my process name
CUtility::splitpath((char *)argv[0],NULL,gszProcName);
#endif
int iScreenLogLevel;		// level of screen diagnostics
int iFileLogLevel;			// level of file diagnostics
char szLogFile[_MAX_PATH];	// write diagnostics to this file

int Rslt;
int TrimLen;

int PMode;                              // processing mode
int MinIdentity;						// minimum required identity
int MinScore;							// minimum required score
int MinMatches;							// minimum required base matches

char szExprName[cMaxIdntNameLen];		// name of this experiment
char szExprDescr[cMaxIdntDescrLen];		// describes experiment
char szBlatParams[cMaxIdntDescrLen];	// Blat parameters used
char szQueryFile[_MAX_PATH];	// Blat query file name
char szTargetFile[_MAX_PATH];   // Blat target file name

char szPSLinFile[_MAX_PATH];	// process from these psl files
char szOutFile[_MAX_PATH];		// output into this SQLite file

// command line args
struct arg_lit  *help    = arg_lit0("h","help",                "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to screen and logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"diagnostics log file");

struct arg_int *pmode = arg_int0("m","mode","<int>",              "processing mode - 0 create new, 1 append (default 0)");
struct arg_int *minidentity = arg_int0("j","minidentity","<int>", "minimum required identity (default 80, range 25..100)");
struct arg_int *minscore = arg_int0("s","minscore","<int>",       "minimum required score (default 200, min 25)");
struct arg_int *minmatches = arg_int0("r","minmatches","<int>",   "minimum required base matches (default 100, min 25)");


struct arg_file *infile = arg_file1("i","infile","<file>",		"load alignments from this Blat generated psl file");
struct arg_file *outfile = arg_file1("o","outfile","<file>",	"store alignments into this SQLite file");

struct arg_file *queryfile = arg_file1("q","queryfile","<file>",	"Blat input query sequences file");
struct arg_file *targetfile = arg_file1("t","targetfile","<file>",	"Blat input target sequences file");

struct arg_str *exprname = arg_str1("e","experiment","<str>",	"name of experiment");
struct arg_str *exprdescr = arg_str0("E","description","<str>",	"describes experiment");
struct arg_str *blatparams = arg_str0("b","parameters","<str>",	"Blat parameters used");

struct arg_end *end = arg_end(40);

void *argtable[] = {help,version,FileLogLevel,LogFile,
					pmode,minidentity,minscore,minmatches,infile,outfile,queryfile,targetfile,exprname,exprdescr,blatparams,
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
		printf("\n      To invoke this parameter file then precede it's name with '@'");
		printf("\n      e.g. %s %s @myparams.txt\n",gszProcName,gpszSubProcess->pszName);
		printf("\nPlease report any issues regarding usage of %s to stuart.stephen@csiro.au\n\n",gszProcName);
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
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Subprocess %s Version %s starting",gpszSubProcess->pszName,cpszProgVer);

	PMode = pmode->count ? pmode->ival[0] : 0;
	if(PMode < 0 || PMode > 1)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Processing mode '-m%d' must be in range 0..1",PMode);
		exit(1);
		}

	MinIdentity = minidentity->count ? minidentity->ival[0] : 80;
	if(MinIdentity < 25 || MinIdentity > 100)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Minimum required identity '-j%d' must be in range 25..100",MinIdentity);
		exit(1);
		}
	MinScore = minscore->count ? minscore->ival[0] : 200;
	if(MinScore < 25)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Minimum required score '-s%d' must be at least 25",MinScore);
		exit(1);
		}
	MinMatches = minmatches->count ? minmatches->ival[0] : 100;
	if(MinMatches < 25)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Minimum required base matches '-r%d' must be at least 25",MinMatches);
		exit(1);
		}


	strcpy(szPSLinFile,infile->filename[0]);
	strcpy(szOutFile,outfile->filename[0]);

	strcpy(szQueryFile,queryfile->filename[0]);
	strcpy(szTargetFile,targetfile->filename[0]);

	strcpy(szExprName,exprname->sval[0]);
	TrimLen = TrimQuotes(szExprName);
	if(TrimLen < 1)
		strcpy(szExprName,"Not specified");

	if(exprdescr->count)
		strcpy(szExprDescr,exprdescr->sval[0]);
	else
		strcpy(szExprDescr,"N/A");
	TrimLen = TrimQuotes(szExprDescr);
	if(TrimLen < 1)
		strcpy(szExprDescr,"N/A");

	if(blatparams->count)
		strcpy(szBlatParams,blatparams->sval[0]);
	else
		strcpy(szBlatParams,"N/A");
	TrimLen = TrimQuotes(szBlatParams);
	if(TrimLen < 1)
		strcpy(szBlatParams,"N/A");


	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing parameters:");
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing mode: %d",PMode);
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Minimum identity: %d",MinIdentity);
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Minimum score: %d",MinScore);
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Minimum base matches: %d",MinMatches);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Input psl file: '%s'",szPSLinFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Output SQLite file:   '%s'",szOutFile);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Experiment name:   '%s'",szExprName);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Blat query sequences file:   '%s'",szQueryFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Blat target sequences file:   '%s'",szTargetFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Experiment descrption:   '%s'",szExprDescr);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Blat parameters used:   '%s'",szBlatParams);
	
	gStopWatch.Start();
#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	Rslt = Process(PMode,MinIdentity,MinScore,MinMatches,szOutFile,szExprName,szPSLinFile,szQueryFile,szTargetFile,szExprDescr,szBlatParams);
	gStopWatch.Stop();
	Rslt = Rslt < 0 ? 1 : 0;
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Exit Code: %d Total processing time: %s",Rslt,gStopWatch.Read());
	exit(Rslt);
	}
else
	{
	printf("\n%s PSL to SQLite, Version %s\n", gszProcName,cpszProgVer);
	arg_print_errors(stdout,end,gszProcName);
	arg_print_syntax(stdout,argtable,"\nUse '-h' to view option and parameter usage\n");
	return(1);	}
return 0;
}

int 
Process(int PMode,					// processing mode, 0 to delete any existing then create new SQLite, 1 to append to existing SQLite
		int MinIdentity,			// minimum required identity
		int MinScore,				// minimum required score
		int MinMatches,				// minimum required base matches
		char *pszDatabase,			// SQLite database file
		char *pszExprName,			// experiment name
		char *pszPSLFile,			// alignments were parsed from this BLAT generated PSL file
		char *pszQueryFile,			// Blat'd query sequences in this file
		char *pszTargetFile,		// against targeted sequences in this file
		char *pszExprDescr,			// describes experiment
		char *pszBlatParams)		// Blat parameters used
{
int Rslt;
CSQLitePSL *pSQLitePSL;

if((pSQLitePSL = new CSQLitePSL)==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate instance of CSQLitePSL");
	return(eBSFerrObj);
	}
Rslt = pSQLitePSL->ProcessPSL2SQLite(PMode,MinIdentity,MinScore,MinMatches,pszDatabase,pszExprName,pszPSLFile,pszQueryFile,pszTargetFile,pszExprDescr,pszBlatParams,0);
delete pSQLitePSL;
return(Rslt);
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


