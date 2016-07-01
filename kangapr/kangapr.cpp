// kangapr.cpp : Defines the entry point for the console application.
// Kanga NG preprocessing NG reads
// Targeted towards error correction and filtering of raw reads
// Required functionality
// Trim fixed number of bases from 5' and/or 3' read ends
// filter out under and/or overlength reads
// flter out reads with more than acceptable number of indeterminate bases
// Trim to a maximum length
// Trim adaptor bases
// Trim low quality 5' and 3' bases
// Error correct reads
// Validate (if paired end) the correct read ordering and that there are none orphan single ended
//
// Release 1.0.1   changed _MAX_PATH to 260 to maintain compatibility with windows

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

#include "kangapr.h"
#include "SampleReads.h"


const char *cpszProgVer = "1.0.1";		// increment with each release

int	TrimSeqChrs(char *pszTxt);	// trims quote marks, space/tabs and validates sequence as only containing acgtu
int Process(etPMode PMode,				// process filtering mode
			bool bSOLiD,				// true if expecting SOLiD colorspace reads
		    int SampleOfs,				// start at this read instance (1 = first read)
			int SampleNth,				// sample every Nth from SampleOfs
			int MaxSamples,				// sample at most this number of reads
			int NumInputP1Files,		// number of input P1, or combined P1 and P2 files
			char **ppszInP1Files,		// input P1 files
			int NumInputP2Files,		// number of input P2 files
			char **ppszInP2Files,		// input P2 files
			char *pszOutP1File,			// output 5' sampled reads to this file
			char *pszOutP2File);			// output 3' sampled reads to this file

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
	return _T("kangapr");
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

int iFileLogLevel;			// level of file diagnostics
int iScreenLogLevel;		// level of file diagnostics
char szLogFile[_MAX_PATH];	// write diagnostics to this file
int Rslt = 0;   			// function result code >= 0 represents success, < 0 on failure
int Idx;

etPMode PMode;				// process filtering mode
bool bSOLiD;				// true if expecting SOLiD colorspace reads

int SampleOfs;				// start at this read instance (1 = first read)

int SampleNth;				// sample every Nth from SampleOfs

int MaxSamples;				// sample at most this number of reads

int NumInputP1Files;			// number of input P1 files
char *pszInP1Files[cRRMaxInFileSpecs];  // input P1 files

int NumInputP2Files;					// number of input P2 files
char *pszInP2Files[cRRMaxInFileSpecs];  // input P2 files

char szOutP1File[_MAX_PATH];			// output 5' sampled reads to this file
char szOutP2File[_MAX_PATH];			// output 3' sampled reads to this file

// command line args
struct arg_lit  *help    = arg_lit0("h","help",                 "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to screen and logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"diagnostics log file");

struct arg_int *pmode = arg_int0("m","mode","<int>",		    "Processing mode: 0 SE reads in P1 input files, 1 PE reads in P1 and P2 input files (default SE reads)");
struct arg_lit  *solid = arg_lit0("C","colorspace",             "Input reads are colorspace (SOLiD)");


struct arg_int *sampleof = arg_int0("s","sampleof","<int>",		"Start sampling from this read (default 1st)");
struct arg_int *samplenth = arg_int0("S","samplenth","<int>",	"Sample every Nth read (1 default)");


struct arg_int *maxsamples = arg_int0("M","Naxsamples","<int>",	"Max number of samples (0 default, sample until last read)");

struct arg_file *inp1files = arg_filen("i","in","<file>",1,cRRMaxInFileSpecs,"Load P1 fasta/fastq sequences from file(s) NO wildcards allowed");
struct arg_file *inp2files = arg_file0("I","in","<file>",		"Load P2 fasta/fastq sequences from file(s), NO wildcards allowed");

struct arg_file *outp1file = arg_file1("o","out","<file>",		"output P1 sampled reads to this file");
struct arg_file *outp2file = arg_file0("O","out","<file>",		"output P2 sampled reads to this file");




struct arg_end *end = arg_end(20);

void *argtable[] = {help,version,FileLogLevel,LogFile,
					pmode,solid,sampleof,samplenth,maxsamples,inp1files,inp2files,outp1file,outp2file,
					end};

char **pAllArgs;
int argerrors;
argerrors = CUtility::arg_parsefromfile(argc,(char **)argv,&pAllArgs);
if(argerrors >= 0)
	argerrors = arg_parse(argerrors,pAllArgs,argtable);

/* special case: '--help' takes precedence over error reporting */
if (help->count > 0)
        {
		printf("\n%s Kanga NGS Preprocess Reads, Version %s\nOptions ---\n", gszProcName,cpszProgVer);
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

	PMode = (etPMode)(pmode->count ? pmode->ival[0] : ePMDefault);
	if(PMode < ePMDefault || PMode >= ePMplaceholder)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Processing mode '-m%d' specified outside of range %d..%d\n",PMode,ePMDefault,(int)ePMplaceholder-1);
		exit(1);
		}

	bSOLiD = solid->count ? true : false;

	SampleOfs = sampleof->count ? sampleof->ival[0] : 1;
	if(SampleOfs < 1 || (UINT32)SampleOfs > 0x07fffffff)
		{




		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Sample reads offset '-s%d' must be in range 1..%d",SampleOfs,0x07fffffff);
		exit(1);
		}


	SampleNth = samplenth->count ? samplenth->ival[0] : 1;
	if(SampleNth < 1 || (UINT32)SampleNth > 0x07fffffff)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Sample Nth '-S%d' must be in range 1..%d",SampleOfs,0x07fffffff);
		exit(1);
		}


	MaxSamples = maxsamples->count ? maxsamples->ival[0] : 0;
	if(MaxSamples < 0 || (UINT32)MaxSamples > 0x07fffffff)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Max samples '-M%d' must be in range 0..%d",MaxSamples,0x07fffffff);

		exit(1);
		}

	if(!inp1files->count)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: No input P1 fasta/fastq file(s) specified with with '-i<filespec>' option)");
		exit(1);
		}

	for(NumInputP1Files=Idx=0;NumInputP1Files < cRRMaxInFileSpecs && Idx < inp1files->count; Idx++)
		{
		pszInP1Files[Idx] = NULL;
		if(pszInP1Files[NumInputP1Files] == NULL)
			pszInP1Files[NumInputP1Files] = new char [_MAX_PATH];
		strncpy(pszInP1Files[NumInputP1Files],inp1files->filename[Idx],_MAX_PATH);
		pszInP1Files[NumInputP1Files][_MAX_PATH-1] = '\0';
		CUtility::TrimQuotedWhitespcExtd(pszInP1Files[NumInputP1Files]);
		if(pszInP1Files[NumInputP1Files][0] != '\0')
			NumInputP1Files++;
		}

	if(!NumInputP1Files)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: After removal of whitespace, no input P1 fasta/fastq file(s) specified with '-i<filespec>' option");
		exit(1);
		}

	strncpy(szOutP1File,outp1file->filename[0],_MAX_PATH);
	szOutP1File[_MAX_PATH-1] = '\0';
	NumInputP2Files = 0;
	pszInP2Files[0] = NULL;

	if(PMode == ePMPairedEnd)
		{
		if(!inp2files->count)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: P2 paired end read files expected to be specified with '-I<readsfile'");
			exit(1);
			}

		for(NumInputP2Files=Idx=0;NumInputP2Files < cRRMaxInFileSpecs && Idx < inp2files->count; Idx++)
			{
			pszInP2Files[Idx] = NULL;
			if(pszInP2Files[NumInputP2Files] == NULL)
				pszInP2Files[NumInputP2Files] = new char [_MAX_PATH];
			strncpy(pszInP2Files[NumInputP2Files],inp2files->filename[Idx],_MAX_PATH);
			pszInP2Files[NumInputP2Files][_MAX_PATH-1] = '\0';
			CUtility::TrimQuotedWhitespcExtd(pszInP2Files[NumInputP2Files]);
			if(pszInP2Files[NumInputP2Files][0] != '\0')
				NumInputP2Files++;
			}

		if(!NumInputP2Files)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: After removal of whitespace, no input P2 fasta/fastq file(s) specified with '-I<filespec>' option");
			exit(1);
			}

		if(NumInputP1Files != NumInputP2Files)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Number of P1 files (%d) specified with '-i<P1filespec>' not same as number of P2 files (%d) specified with '-I<P2filespec>' option",NumInputP1Files,NumInputP2Files);
			exit(1);
			}


		if(!outp2file->count)
		{

			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Paired end sampling requested but no PE2 output file specified");
		exit(1);
			}
		strncpy(szOutP2File,outp2file->filename[0],_MAX_PATH);
		szOutP2File[_MAX_PATH-1] = '\0';
		}
	else
		szOutP2File[0] = '\0';

// show user current resource limits
#ifndef _WIN32
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Resources: %s",CUtility::ReportResourceLimits());
#endif

	// report back to user the parameter settings
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing parameters:");
	const char *pszDescr;

	switch(PMode) {
		case ePMPairedEnd:
			pszDescr = "paired end reads expected in separate P1 and P2 file(s)";
			break;

		default:
			PMode = ePMDefault;
			pszDescr = "single ended reads expected in P1 input file(s)";
			break;
		}

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Processing mode is : '%s'",pszDescr);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Reads are expected to be in : '%s'",bSOLiD ? "Colorspace" : "Basespace");



	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Sample input sequences starting from : %d",SampleOfs);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Sample every Nth sequence : %d",SampleNth);
	if(MaxSamples == 0)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Sample at most this many sequences : No limit");
	else
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Sample at most this many sequences : %d",MaxSamples);


	for(Idx=0; Idx < NumInputP1Files; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Input P1 sequences file (%d) : '%s'",Idx+1,pszInP1Files[Idx]);

	if(PMode == ePMPairedEnd)
		{
		for(Idx=0; Idx < NumInputP2Files; Idx++)
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Input P2 sequences file (%d) : '%s'",Idx+1,pszInP2Files[Idx]);
		}

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Output sampled P1 reads file: '%s'",szOutP1File);
	if(PMode == ePMPairedEnd)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Output sampled P2 reads file: '%s'",szOutP2File);


#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	gStopWatch.Start();

	Rslt = Process(PMode,					// process filtering mode
					bSOLiD,					// true if expecting SOLiD colorspace reads
					SampleOfs,				// start at this read instance (1 = first read)
					SampleNth,				// sample every Nth from SampleOfs
					MaxSamples,				// sample at most this number of reads
					NumInputP1Files,		// number of input P1, or combined P1 and P2 files
					pszInP1Files,			// input P1 files
					NumInputP2Files,		// number of input P2 files
					pszInP2Files,			// input P2 files
					szOutP1File,			// output 5' sampled reads to this file
					szOutP2File);			// output 3' sampled reads to this file

	gStopWatch.Stop();
	Rslt = Rslt >=0 ? 0 : 1;
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Exit code: %d Total processing time: %s",Rslt,gStopWatch.Read());
	exit(Rslt);
	}
else
	{
	printf("\n%s Kanga NGS Preprocess Reads, Version %s\n",gszProcName,cpszProgVer);
	arg_print_errors(stdout,end,gszProcName);
	arg_print_syntax(stdout,argtable,"\nUse '-h' to view option and parameter usage\n");
	exit(1);
	}
return 0;

}

// TrimSeqChrs
// Removes any quotes and whitespace (space and tabs only) from pszTxt
// Also checks for legal base chars 'acgt'
// Returns -1 if illegal char, 0 if empty sequence, otherwise the sequence length
int
TrimSeqChrs(char *pszTxt)
{
char *pDst = pszTxt;
char Chr;
int Len = 0;
while((Chr = *pszTxt++))
	{
	if(Chr == '\0' || (Chr == '\t'  || Chr == ' '  || Chr == '"' || Chr == '\''))
		continue;
	switch(Chr) {
		case 'a': case 'A':
			Chr = 'a';
			break;
		case 'c': case 'C':
			Chr = 'c';
			break;
		case 'g': case 'G':
			Chr = 'g';
			break;
		case 't': case 'T': case 'u': case 'U':
			Chr = 't';
			break;
		default:
			return(-1);
		}
	*pDst++ = Chr;
	Len++;
	}
*pDst = Chr;
return(Len);
}


int Process(etPMode PMode,				// process filtering mode
			bool bSOLiD,				// true if expecting SOLiD colorspace reads
		    int SampleOfs,				// start at this read instance (1 = first read)
			int SampleIncr,				// sample every Nth from SampleOfs
			int SampleMax,				// sample at most this number of reads
			int NumInputP1Files,		// number of input P1, or combined P1 and P2 files
			char **ppszInP1Files,		// input P1 files
			int NumInputP2Files,		// number of input P2 files
			char **ppszInP2Files,		// input P2 files
			char *pszOutP1File,			// output 5' sampled reads to this file
			char *pszOutP2File)			// output 3' sampled reads to this file
{
int Rslt;
int NumSampled;
int FileIdx;
char *pszInP1File;
char *pszInP2File;
CSampleReads SampleReads;

NumSampled = 0;
for(FileIdx = 0; FileIdx < NumInputP1Files; FileIdx++)
	{
	pszInP1File = ppszInP1Files[FileIdx];
	if(PMode == ePMPairedEnd)
		pszInP2File = ppszInP2Files[FileIdx];
	else
		pszInP2File = NULL;
	if((Rslt = SampleReads.GenSamples(PMode,FileIdx == 0 ? false : true,SampleOfs,SampleIncr,SampleMax,pszInP1File,pszInP2File,pszOutP1File,pszOutP2File)) != eBSFSuccess)
		break;
//	NumSampled = SampleReads.NumSamples();
	if(SampleMax > 0)
		{
		SampleMax -= NumSampled;
		if(SampleMax == 0)
			break;
		}
//	SampleOfs = SampleReads.NumBeforeNxtSample();
	}
return(0);
}

