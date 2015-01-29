// genreadsde.cpp : Defines the entry point for the console application.
// Generate reads differential expression - or k-mer differential expression where k is the read length
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

const char *cpszProgVer = "0.0.2";		// increment with each release


teBSFrsltCodes Process(etPRRMode PMode,		// processing mode
   		UINT32 NumReadsLimit,				// limit processing to this many reads
		int Trim5,							// trim this many bases off leading sequence 5' end
		int Trim3,							// trim this many bases off trailing sequence 3' end
		int	MinSampleCnts,					// minimum sample counts
		int	MinTotalCnts,					// minimum total samples counts
		int NumInputFileSpecs,				// number of input file specs 
		char *pszInfileSpecs[],				// names of inputs file (wildcards allowed unless in dump mode) containing raw reads
		char *pszOutFile);					// output into this file only if not NULL or not '\0' 


CStopWatch gStopWatch;
CDiagnostics gDiagnostics;				// for writing diagnostics messages to log file
char gszProcName[_MAX_FNAME];			// process name
bool gbActivity;						// used to determine if activity messages vi printf's can be used - output activity if eDLInfo or less

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
	return _T("kangarde");
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
int iScreenLogLevel;		// level of file diagnostics
int iFileLogLevel;			// level of file diagnostics
char szLogFile[_MAX_PATH];	// write diagnostics to this file
int Rslt = 0;   			// function result code >= 0 represents success, < 0 on failure
int Idx;					// general iteration indexer

etPRRMode PMode;				// processing mode

int NumInputFiles;			// number of input files
char *pszInfileSpecs[cRRMaxInFileSpecs];  // input (wildcards allowed if single ended) raw sequencer read files
char szOutfile[_MAX_PATH];  // output rds file or csv if user requested stats

int NumReadsLimit;			// only process this many reads, useful whilst testing workflow
int Trim5;					// trim this many bases off leading sequence 5' end
int Trim3;					// trim this many bases off trailing sequence 3' end
int MinSampleCnts;			// only report if at least one sample has this minimum number of counts and 
int MinTotalCnts;			// only report if total number of counts is at least this

// command line args
struct arg_lit  *help    = arg_lit0("hH","help",                "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to screen and logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"diagnostics log file");
struct arg_int *pmode = arg_int0("m","mode","<int>",		    "processing mode: 0 - single end reads processing (default = 0)");
struct arg_int *trim5 = arg_int0("t","trim5","<int>",		    "trim this many bases off leading sequence 5' end (default = 0)");
struct arg_int *trim3 = arg_int0("T","trim3","<int>",		    "trim this many bases off trailing sequence 3' end (default = 0)");
struct arg_int *minsamplecnts = arg_int0("s","minsamplecnts","<int>","any sample must have at least this number of cnts (default = 5)");
struct arg_int *mintotalcnts = arg_int0("S","mintotalcnts","<int>",	 "total of all samples must have at least this number of cnts (default = 10)");

struct arg_file *inreadfiles = arg_filen("i","in","<file>",0,cRRMaxInFileSpecs,"input from these raw sequencer read files, wildcards allowed if single ended");
struct arg_file *outfile = arg_file1("o","out","<file>",		"output read frequencies to this file");
struct arg_int *readslimit = arg_int0("n","numreadslimit","<int>","limit number of reads in each input file to this many - 0 (default) if no limit");

struct arg_end *end = arg_end(20);

void *argtable[] = {help,version,FileLogLevel,LogFile,
					pmode,trim5,trim3,minsamplecnts,mintotalcnts,inreadfiles,outfile,readslimit,
					end};

char **pAllArgs;
int argerrors;
argerrors = CUtility::arg_parsefromfile(argc,(char **)argv,&pAllArgs);
if(argerrors >= 0)
	argerrors = arg_parse(argerrors,pAllArgs,argtable);

/* special case: '--help' takes precedence over error reporting */
if (help->count > 0)
        {
		printf("\n%s the K-mer Adaptive Next Generation Aligner Reads preprocessor, Version %s\nOptions ---\n", gszProcName,cpszProgVer);
        arg_print_syntax(stdout,argtable,"\n");
        arg_print_glossary(stdout,argtable,"  %-25s %s\n");
		printf("\nNote: Parameters can be entered into a parameter file, one parameter per line.");
		printf("\n      To invoke this parameter file then precede it's name with '@'");
		printf("\n      e.g. %s @myparams.txt\n",gszProcName);
		printf("\nPlease report any issues regarding usage of %s to stuart.stephen@csiro.au\n\n",gszProcName);
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
	gbActivity = true;
	if(FileLogLevel->count && !LogFile->count)
		{
		printf("\nError: FileLogLevel '-f%d' specified but no logfile '-F<logfile>'",FileLogLevel->ival[0]);
		exit(1);
		}

	iScreenLogLevel = iFileLogLevel = FileLogLevel->count ? FileLogLevel->ival[0] : eDLInfo;
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

			// now that log parameters have been parsed then initialise diagnostics log system
	if(!gDiagnostics.Open(szLogFile,(etDiagLevel)iScreenLogLevel,(etDiagLevel)iFileLogLevel,true))
		{
		printf("\nError: Unable to start diagnostics subsystem\n");
		if(szLogFile[0] != '\0')
			printf(" Most likely cause is that logfile '%s' can't be opened/created\n",szLogFile);
		exit(1);
		}
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Version: %s",cpszProgVer);


	PMode = (etPRRMode)(pmode->count ? pmode->ival[0] : ePMRRNewSingle);
	if(PMode < 0 || PMode >= ePMRRplaceholder)
		{
		printf("\nError: Processing mode '-m%d' specified outside of range %d..%d",PMode,0,(int)ePMRRplaceholder-1);
		exit(1);
		}

	
	Trim5 = trim5->count ? trim5->ival[0] : 0;
	if(Trim5 < 0 || Trim5 > 99)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: 5' end trim '-t%d' specified outside of range 0..99",Trim5);
		exit(1);
		}
	Trim3 = trim3->count ? trim3->ival[0] : 0;
	if(Trim3 < 0 || Trim3 > 99)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: 3' end trim '-t%d' specified outside of range 0..99",Trim3);
		exit(1);
		}

	MinSampleCnts = minsamplecnts->count ? minsamplecnts->ival[0] : 5;
	if(MinSampleCnts < 1 || MinSampleCnts > 100000)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: minimum sample counts '-s%d' specified outside of range 1..100000",MinSampleCnts);
		exit(1);
		}

	MinTotalCnts = mintotalcnts->count ? mintotalcnts->ival[0] : 10;
	if(MinTotalCnts < 1 || MinTotalCnts > 1000000)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: minimum total samples counts '-S%d' specified outside of range %d..1000000",MinSampleCnts,MinTotalCnts);
		exit(1);
		}


	NumReadsLimit = readslimit->count ? readslimit->ival[0] : 0;
	NumReadsLimit &= 0x7fffffff;	// simply to ensure limit can't be later taken as a negative value...

	if(!inreadfiles->count)
		{
		printf("\nError: No input files specified with with '-i<filespec>' option)");
		exit(1);
		}

	for(NumInputFiles=Idx=0;NumInputFiles < cRRMaxInFileSpecs && Idx < inreadfiles->count; Idx++)
		{
		pszInfileSpecs[Idx] = NULL;
		if(pszInfileSpecs[NumInputFiles] == NULL)
			pszInfileSpecs[NumInputFiles] = new char [_MAX_PATH];
		strncpy(pszInfileSpecs[NumInputFiles],inreadfiles->filename[Idx],_MAX_PATH);
		pszInfileSpecs[NumInputFiles][_MAX_PATH-1] = '\0';
		CUtility::TrimQuotedWhitespcExtd(pszInfileSpecs[NumInputFiles]);
		if(pszInfileSpecs[NumInputFiles][0] != '\0')
			NumInputFiles++;
		}

	if(!NumInputFiles)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: After removal of whitespace, no input file%s specified with '-i<filespec>' option)", PMode >= ePMRRStats ? "" : "(s)");
		exit(1);
		}

	strncpy(szOutfile,outfile->filename[0],_MAX_PATH);
	szOutfile[_MAX_PATH-1] = '\0';
	CUtility::TrimQuotedWhitespcExtd(szOutfile);
	if(szOutfile[0] == '\0')
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: No output file ('-o<file> option') specified after removal of leading/trailing quotes and whitespace");
		exit(1);
		}

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing parameters:");

	const char *pszProcMode;
	switch(PMode) {
		case ePMRRNewSingle:
			pszProcMode = "Process raw reads for differential expression";
			break;

		};

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Processing mode: '%s'",pszProcMode);


	if(NumReadsLimit > 0)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"limiting processing to the first %d reads",NumReadsLimit);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Trim %d bases from 5' end and %d bases from 3' end",Trim5,Trim3);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Minimum sample counts: %d",MinSampleCnts);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Minimum total samples counts: %d",MinTotalCnts);

	for(Idx=0; Idx < NumInputFiles; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"input raw reads files (%d): '%s'",Idx+1,pszInfileSpecs[Idx]);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"output file to create: '%s'", szOutfile);
#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	gStopWatch.Start();
	Rslt = Process(PMode,			// processing mode
					NumReadsLimit,				// limit processing to this many reads
					Trim5,						// trim this many bases off leading sequence 5' end
					Trim3,						// trim this many bases off trailing sequence 3' end
					MinSampleCnts,				// minimum sample counts
					MinTotalCnts,				// minimum total samples counts
					NumInputFiles,				// number of input file specs
					pszInfileSpecs,				// names of inputs file (wildcards allowed unless in dump mode) containing raw reads
					szOutfile);					// output into this file

	gStopWatch.Stop();
	Rslt = Rslt >=0 ? 0 : 1;
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Exit code: %d Total processing time: %s",Rslt,gStopWatch.Read());
	exit(Rslt);
	}
else
	{
	printf("\n%s the K-mer Adaptive Next Generation Aligner Reads preprocessor, Version %s\n", gszProcName,cpszProgVer);
	arg_print_errors(stdout,end,gszProcName);
	arg_print_syntax(stdout,argtable,"\nUse '-h' to view option and parameter usage\n");
	exit(1);
	}
return 0;
}

teBSFrsltCodes
Process(etPRRMode PMode,						// processing mode
		UINT32 NumReadsLimit,				// limit processing to this many reads
		int Trim5,							// trim this many bases off leading sequence 5' end
		int Trim3,							// trim this many bases off trailing sequence 3' end
		int	MinSampleCnts,					// minimum sample counts
		int	MinTotalCnts,					// minimum total samples counts
		int NumInputFileSpecs,				// number of input file specs 
		char *pszInfileSpecs[],				// names of inputs file (wildcards allowed unless in dump mode) containing raw reads
		char *pszOutFile)					// output into this file
{
teBSFrsltCodes Rslt;
CProcRawReads RawReads;
RawReads.ReportActivity(true);

switch(PMode) {
	case ePMRRNewSingle:
		Rslt = RawReads.LoadAndProcessReadsDE(PMode,NumReadsLimit,Trim5,Trim3,MinSampleCnts,MinTotalCnts,NumInputFileSpecs,pszInfileSpecs,pszOutFile);
		break;

	}
return(Rslt);
}
