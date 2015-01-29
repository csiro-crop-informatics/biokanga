// denovo.cpp : Defines the entry point for the console application.
// Currently this code is targeting pair end reads in which the 3' sequence is antisense to the 5' sequence and
// there is an overlap. An output multifasta file is generated in which the 3' sequences will be revcpl'd and overlaid
// on to their 5' respective prime reads and output as merged sequences. 
//
// 1.0.0	 inital release with fastq as an output option
// 1.0.1     modified allowed overlay substitutions to allow for higher sequencer error rates at the 3' end of reads
// 1.0.3     fix for potential buffer overflow causing trunctation of output quality scores
// 1.0.4     changed orphan output processing

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

#include "./MergeReadPairs.h"

const char *cpszProgVer = "1.0.4";			// inital release with fastq as an output option

const int cDfltOverlapPercSubs = 5;			// default allowed substitutions as a percentage of the overlap
const int cDfltRequiredOverlap = 10;		// default paired ends must overlap by at least this many bases
const int cMaxMinRequiredOverlap  = 200;	// max paired ends must overlap by at least this many bases


const int cMaxInFileSpecs = 10;				// allow at most 10 files for input per end

const int cMaxWorkerThreads = 64;			// allow up to this many worker threads



int
Process(etPMode PMode,				// processing mode
	    etOFormat OFormat,			// output file format
		int MinOverlap,				// overlaps must be of at least this length
		int MaxSubPerc,				// and percentage of substitutions in overlap must be no more than this
		int NumInPE5Files,			// number of input single ended or 5' end if paired end reads file
		char **pszInPE5Files,		// input single ended or 5' end if paired end reads files
		int NumInPE3Files,			// number of input input 3' end if paired end reads files
		char **pszInPE3Files,		// input 3' end if paired end reads files
		char *pszOutCtgsFile);		// output file

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
	return _T("denovo");
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

int Idx;

etPMode PMode;				// processing mode
etOFormat OFormat;			// output file format

int MinOverlap;				// overlaps must be of at least this length
int MaxSubPerc;				// and percentage of substitutions in overlap must be no more than this

int NumInPE5Files;								// number of input single ended or 5' end if paired end reads file
char *pszInPE5Files[cMaxInFileSpecs];			// input single ended or 5' end if paired end reads files
int NumInPE3Files;								// number of input input 3' end if paired end reads files
char *pszInPE3Files[cMaxInFileSpecs];			// input 3' end if paired end reads files
char szOutCtgsFile[_MAX_PATH];					// output assembled contigs file

// command line args
struct arg_lit  *help    = arg_lit0("hH","help",                "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to screen and logfile fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"diagnostics log file");

struct arg_int *pmode = arg_int0("m","mode","<int>",		    "merge processing : 0 - output overlaped reads only, 1 - output overlape and orphan reads combined, 2 - output overlaped and orphan reads separately (default is 0)");
struct arg_int *oformat = arg_int0("M","oformat","<int>",		"output format as : 0 - auto, 1 - fasta, 2 - fastq (default 0, uses input file format)");

struct arg_int *minoverlap = arg_int0("l","minoverlap","<int>",	"paired end 3' reads must overlap onto 5' reads by at least this number of base (minimum 1, default 10)");
struct arg_int *naxsubperc = arg_int0("s","maxsubperc","<int>",	"allow at most this percentage of substitutions in the paired end overlaps (maximum 10, default 5)");

struct arg_file *inpe5files = arg_filen("i","inpe5","<file>",1,cMaxInFileSpecs, "input P1 5' end raw read files (wildcards not allowed, fasta or fastq)");
struct arg_file *inpe3files = arg_filen("I","inpe3","<file>",1,cMaxInFileSpecs, "input P2 3' end raw read files (wildcards not allowed, fasta or fastq)");
struct arg_file *outctgsfile = arg_file1("o","outctgsfile","<file>", "output merged pair sequences to this file");

struct arg_end *end = arg_end(20);

void *argtable[] = {help,version,FileLogLevel,LogFile,
					pmode,oformat,minoverlap,naxsubperc,
					inpe5files,inpe3files,outctgsfile,
					end};

char **pAllArgs;
int argerrors;
argerrors = CUtility::arg_parsefromfile(argc,(char **)argv,&pAllArgs);
if(argerrors >= 0)
	argerrors = arg_parse(argerrors,pAllArgs,argtable);

/* special case: '--help' takes precedence over error reporting */
if (help->count > 0)
        {
		printf("\n%s Overlap paired end merger, Version %s\nOptions ---\n", gszProcName,cpszProgVer);
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

	PMode = (etPMode)(pmode->count ? pmode->ival[0] : ePMdefault);
	if(PMode < ePMdefault || PMode >= ePMplaceholder)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Processing mode '-m%d' specified outside of range %d..%d",PMode,0,(int)ePMplaceholder-1);
		exit(1);
		}

	OFormat = (etOFormat)(oformat->count ? oformat->ival[0] : eOFauto);
	if(OFormat < eOFauto || OFormat >= eOFplaceholder)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Output format '-M%d' specified outside of range %d..%d",OFormat,0,(int)eOFplaceholder-1);
		exit(1);
		}

	MinOverlap = minoverlap->count ? minoverlap->ival[0] : cDfltRequiredOverlap;
	if(MinOverlap < cMinOverlapLen || MinOverlap >= cMaxMinRequiredOverlap)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Minimum required paired end overlap '-S%d' specified outside of range %d..%d",MinOverlap,cMinOverlapLen,cMaxMinRequiredOverlap);
		exit(1);
		}

	MaxSubPerc = naxsubperc->count ? naxsubperc->ival[0] : cDfltOverlapPercSubs;
	if(MaxSubPerc < 0 || MaxSubPerc > cMaxOverlapPercSubs)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Maximum accepted overlap substitution percentage '-s%d' specified outside of range %d..%d",MaxSubPerc,0,cMaxOverlapPercSubs);
		exit(1);
		}

	for(NumInPE5Files=Idx=0;NumInPE5Files < cMaxInFileSpecs && Idx < inpe5files->count; Idx++)
		{
		pszInPE5Files[Idx] = NULL;
		if(pszInPE5Files[NumInPE5Files] == NULL)
			pszInPE5Files[NumInPE5Files] = new char [_MAX_PATH];
		strncpy(pszInPE5Files[NumInPE5Files],inpe5files->filename[Idx],_MAX_PATH);
		pszInPE5Files[NumInPE5Files][_MAX_PATH-1] = '\0';
		CUtility::TrimQuotedWhitespcExtd(pszInPE5Files[NumInPE5Files]);
		if(pszInPE5Files[NumInPE5Files][0] != '\0')
			NumInPE5Files++;
		}

	if(!NumInPE5Files)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: After removal of whitespace, no input P1 5' end file(s) specified with '-i<filespec>' option)\n");
		exit(1);
		}

	NumInPE3Files = 0;
	for(Idx=0;NumInPE3Files < cMaxInFileSpecs && Idx < inpe3files->count; Idx++)
		{
		pszInPE3Files[Idx] = NULL;
		if(pszInPE3Files[NumInPE3Files] == NULL)
			pszInPE3Files[NumInPE3Files] = new char [_MAX_PATH];
		strncpy(pszInPE3Files[NumInPE3Files],inpe3files->filename[Idx],_MAX_PATH);
		pszInPE3Files[NumInPE3Files][_MAX_PATH-1] = '\0';
		CUtility::TrimQuotedWhitespcExtd(pszInPE3Files[NumInPE3Files]);
		if(pszInPE3Files[NumInPE3Files][0] != '\0')
			NumInPE3Files++;
		}

	if(!NumInPE3Files)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: After removal of whitespace, no input P2 3' end file(s) specified with '-I<filespec>' option)\n");
		exit(1);
		}

	if(NumInPE5Files != NumInPE3Files)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Number (%d) of P1 5' input files must equal number (%d) of P2 3' input files\n",NumInPE5Files,NumInPE3Files);
		exit(1);
		}

	strcpy(szOutCtgsFile,outctgsfile->filename[0]);					// output file

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing parameters:");

	const char *pszProcMode;
	switch(PMode) {
		case ePMdefault:
			pszProcMode = "Overlaped read sequences file output";
			break;

		case ePMcombined:
			pszProcMode = "Overlaped and orphan read sequences combined file output";
			break;

		case ePMseparate:
			pszProcMode = "Overlaped and orphan read sequences separate output files";
			break;

		}
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Processing mode: '%s'",pszProcMode);

	switch(OFormat) {
		case eOFauto:
			pszProcMode = "Same format as input files";
			break;

		case eOFfasta:
			pszProcMode = "Basespase as fasta";
			break;

		case eOFfastq:
			pszProcMode = "Basespace with quality scores as fastq";
			break;

		}
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Output format: '%s'",pszProcMode);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Minimum overlap length: %d",MinOverlap);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Maximum percentage of substitutions in overlap: %d%%",MaxSubPerc);
	

	for(Idx=0; Idx < NumInPE5Files; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"input P1 5' reads files (%d): '%s'",Idx+1,pszInPE5Files[Idx]);
	for(Idx=0; Idx < NumInPE3Files; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"input P2 3' reads files (%d): '%s'",Idx+1,pszInPE3Files[Idx]);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"output assembled contigs file to create: '%s'", szOutCtgsFile);

#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	gStopWatch.Start();
	Rslt = Process(PMode,OFormat,MinOverlap,MaxSubPerc,NumInPE5Files,pszInPE5Files,NumInPE3Files,pszInPE3Files,szOutCtgsFile);
	gStopWatch.Stop();
	Rslt = Rslt >=0 ? 0 : 1;
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Exit code: %d Total processing time: %s",Rslt,gStopWatch.Read());
	exit(Rslt);
	}
else
	{
	printf("\n%s Overlap paired end merger, Version %s\n",gszProcName,cpszProgVer);
	arg_print_errors(stdout,end,gszProcName);
	arg_print_syntax(stdout,argtable,"\nUse '-h' to view option and parameter usage\n");
	exit(1);
	}
return 0;
}


int
MergeOverlaps(etPMode PMode,				// processing mode 
			  etOFormat OFormat,			// output file format
			  bool bAppendOut,				// if true then append if output files exist, otherwise trunctate existing output files 
			  int MergedSeqID,				// starting merged sequence identifier
			  int MinOverlap,				// reads must overlap by at least this many bases
              int MaxSubs,					// and overlap can have at most this many substitutions
			  char *pszIn5ReadsFile,		// 5' reads are in this file
			  char *pszIn3ReadsFile,		// 3' reads are in this file
			  char *pszMergeOutFile)		// write merged overlaping reads to this file
{
int Rslt;
CMergeReadPairs *pMergeReads = NULL;
if((pMergeReads = new CMergeReadPairs()) == NULL)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Unable to instantiate CMergeReads");
	return(-1);
	}
Rslt = pMergeReads->MergeOverlaps(PMode,OFormat,bAppendOut,MergedSeqID,MinOverlap,MaxSubs,pszIn5ReadsFile,pszIn3ReadsFile,pszMergeOutFile);

return(Rslt);
}

int
Process(etPMode PMode,				// processing mode
	    etOFormat OFormat,			// output file format
		int MinOverlap,				// overlaps must be of at least this length
		int MaxSubPerc,				// and percentage of substitutions in overlap must be no more than this
		int NumInPE5Files,			// number of input single ended or 5' end if paired end reads file
		char **pszInPE5Files,		// input single ended or 5' end if paired end reads files
		int NumInPE3Files,			// number of input input 3' end if paired end reads files
		char **pszInPE3Files,		// input 3' end if paired end reads files
		char *pszOutCtgsFile)		// output file
{
int Rslt;
int Idx;
int GenMergeSeqID;
GenMergeSeqID = 0;
for(Idx = 0; Idx < NumInPE5Files; Idx++)
	{
	if((Rslt = MergeOverlaps(PMode,OFormat,GenMergeSeqID == 0 ? false : true,GenMergeSeqID,MinOverlap,MaxSubPerc,pszInPE5Files[Idx],pszInPE3Files[Idx],pszOutCtgsFile)) < 0)
		{
		break;
		}
	GenMergeSeqID += Rslt;
	}

return(Rslt);
}