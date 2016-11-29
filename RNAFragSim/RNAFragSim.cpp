// RNAFragSim.cpp : Defines the entry point for the console application.
// RNA sequencing fragmentation simulator
// Given a list of transcripts and their occurances attempts to fragment these assuming each
// base has equal probabilty of being a fragmentation site taking into account 5' and 3' steric occlusion zones

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

#include "./Transcriptome.h"

const char *cpszProgVer = "0.0.5";		// increment with each release

int Process(int TargFragLen,			// Fragment until <= mean fragment length
			int Steric5Len,				// length of 5' steric occlusion region
			int Steric3Len,				// length of 3' steric occlusion region
			int MinFragLen,				// minimum fragment length
			int MaxFragLen,				// maximum fragment length
			char *pszTransFile,			// input transcribed region (CSV) file
			char *pszFragFile,			// output fragment counts (CSV) file
			char *pszDistFile);			// output fragment length distributions (CSV) file

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

int TargFragLen;	    	// fragment transcripts until fragments are <= this mean length
int Steric5Len;				// length of 5' steric occlusion region
int Steric3Len;				// length of 3' steric occlusion region
int MinFragLen;				// report minimum fragment length
int MaxFragLen;				// report maximum fragment length

char szTransFile[_MAX_PATH]; // input transcribed region (CSV) file
char szFragFile[_MAX_PATH];	 // output fragment counts (CSV) file
char szDistFile[_MAX_PATH];	 // output length distribution (CSV) file

// command line args
struct arg_lit  *help    = arg_lit0("hH","help",                "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"diagnostics log file");

struct arg_file *transfile = arg_file1("i","intrans","<file>",	"input transcribed region (CSV) file");
struct arg_file *fragfile = arg_file1("o","outfrags","<file>",	"output transcribed region fragment counts (CSV) file");
struct arg_file *distfile = arg_file1("O","outdist","<file>",	"output fragment length distribution (CSV) file");

struct arg_int  *targfraglen = arg_int0("n","targfraglen","<int>", "fragment until this mean fragment length (default 300, range 100..1000)");
struct arg_int  *steric5len = arg_int0("s","steric5len","<int>",   "length of 5' steric occlusion zone (default = 30, range 5..100");
struct arg_int  *steric3len = arg_int0("S","steric3len","<int>",   "length of 3' steric occlusion zone (default = 30, range 5..100");
struct arg_int  *minfraglen = arg_int0("l","minfraglen","<int>",   "report minimum fragment length (default 125, range 100..1000)");
struct arg_int  *maxfraglen = arg_int0("L","maxfraglen","<int>",   "report maximum fragment length (default 175, range (minfraglen+10)..1100)");

struct arg_end *end = arg_end(20);

void *argtable[] = {help,version,FileLogLevel,LogFile,
					targfraglen, steric5len, steric3len, minfraglen, maxfraglen, transfile, fragfile,distfile,
					end};

char **pAllArgs;
int argerrors;
argerrors = CUtility::arg_parsefromfile(argc,(char **)argv,&pAllArgs);
if(argerrors >= 0)
	argerrors = arg_parse(argerrors,pAllArgs,argtable);

    /* special case: '--help' takes precedence over error reporting */
if (help->count > 0)
        {
		printf("\n%s simulate RNA transcript fragmentation process, Version %s\nOptions ---\n", gszProcName,cpszProgVer);
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
		printf("\n%s Version %s\n",gszProcName,cpszProgVer);
		exit(1);
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

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Version: %s",cpszProgVer);


	TargFragLen = targfraglen->count ? targfraglen->ival[0] : cDfltTargFragLen;
	if(TargFragLen < 100 || TargFragLen > 1000)
		{
		printf("\nError: Targeted mean fragmentation length '-n%d' must be in range %d..%d",TargFragLen,100,1000);
		exit(1);
		}

	Steric5Len = steric5len->count ? steric5len->ival[0] : cDfltStericOcc5Len;
	if(Steric5Len < 5 || Steric5Len > 100)
		{
		printf("\nError: 5' steric occlusion zone length '-s%d' must be in range %d..%d",Steric5Len,5,100);
		exit(1);
		}

	Steric3Len = steric3len->count ? steric3len->ival[0] : cDfltStericOcc3Len;
	if(Steric3Len < 5 || Steric3Len > 100)
		{
		printf("\nError: 3' steric occlusion region length '-S%d' must be in range %d..%d",Steric3Len,5,100);
		exit(1);
		}

	MinFragLen = minfraglen->count ? minfraglen->ival[0] : cDfltMinFragLen;
	if(MinFragLen < 100 || MinFragLen > 1000)
		{
		printf("\nError: Minimum reported fragment length '-l%d' must be in range %d..%d",MinFragLen,100,1000);
		exit(1);
		}

	MaxFragLen = maxfraglen->count ? maxfraglen->ival[0] : cDfltMaxFragLen;
	if(MaxFragLen < (MinFragLen+10) || MaxFragLen > 1100)
		{
		printf("\nError: Maximum reported fragment length '-L%d' must be in range %d..%d",MaxFragLen,(MinFragLen+10),1100);
		exit(1);
		}


	strncpy(szTransFile,transfile->filename[0],_MAX_PATH);
	szTransFile[_MAX_PATH-1] = '\0';
	

	strncpy(szFragFile,fragfile->filename[0],_MAX_PATH);
	szFragFile[_MAX_PATH-1] = '\0';

	strncpy(szDistFile,distfile->filename[0],_MAX_PATH);
	szDistFile[_MAX_PATH-1] = '\0';

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing parameters:");

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Input transcribed regions from file : '%s'",szTransFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Output transcribed regions fragment counts to file : '%s'",szFragFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Output transcribed regions fragment length distributions to file : '%s'",szDistFile);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Fragment until <= mean fragment length: %d",TargFragLen);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"5' steric occlusion region length: %d",Steric5Len);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"3' steric occlusion region length: %d",Steric3Len);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Report on fragments >= length: %d",MinFragLen);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Report on fragments <= length: %d",MaxFragLen);


#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	gStopWatch.Start();
	Rslt = Process(TargFragLen,Steric5Len,Steric3Len,MinFragLen,MaxFragLen,szTransFile,szFragFile,szDistFile);

	gStopWatch.Stop();
	Rslt = Rslt >=0 ? 0 : 1;
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Exit code: %d Total processing time: %s",Rslt,gStopWatch.Read());
	exit(Rslt);
	}
else
	{
	printf("\n%s simulate RNA transcript fragmentation process, Version %s\n",gszProcName,cpszProgVer);
	arg_print_errors(stdout,end,gszProcName);
	arg_print_syntax(stdout,argtable,"\nUse '-h' to view option and parameter usage\n");
	exit(1);
	}
}


int Process(int TargFragLen,			// Fragment until <= mean fragment length
			int Steric5Len,				// length of 5' steric occlusion region
			int Steric3Len,				// length of 3' steric occlusion region
			int MinFragLen,				// minimum fragment length
			int MaxFragLen,				// maximum fragment length
			char *pszTransFile,			// input transcribed region (CSV) file
			char *pszFragFile,			// output fragment counts (CSV) file
			char *pszDistFile)			// output fragment length distributions (CSV) file
{
int Rslt;
CTranscriptome Transcriptome;

Rslt = Transcriptome.Fragmentate(pszTransFile,pszFragFile,pszDistFile,TargFragLen,10,MinFragLen,MaxFragLen,Steric5Len,Steric3Len);

return(Rslt);
}