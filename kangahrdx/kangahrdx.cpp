// kangacrdx.cpp : Defines the entry point for the console application.
//
// Release 0.0.9	Restructured build directories to reflect application name
// Release 0.0.13   Reduced default homozygotic rate down to 3
// Release 0.0.14   New option allowing user to specify the contig descriptor line prefix string (defaults as 'HrdxCtg%d')
// Release 0.0.15   Can handle Can handle tsRawReadV6 kangar preprocessed reads
// Release 1.11.0   public binary release
// Release 1.12.0   public binary release
// Release 1.12.2   increased limit on number of bases compared when sorting sequences
// Release 1.12.3   changed _MAX_PATH to 260 to maintain compatibility with windows
// Release 2.0.0    removed limits on assembly size

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

#include "./HomozyReduce.h"

const char *cpszProgVer = "2.0.1";		// increment with each release

int
Process(etPMode PMode,					// processing mode
		bool bStrand,					// strand specific homozygous region reduction - homozygous regions between any two contigs must be in same orientation
		int MaxHomozySubs,				// characterise as homozygotic if substitution rate between regions <= this rate per 100bp
		int MinHomozyLen,				// homozygotic regions to be at least this length
		int MinHetrozyLen,				// island (marked as homozygotic either flank) hetrozygotic regions must be at least this length otherwise treat as homozygotic
		int MaxNs,						// filter out input contigs having higher than this number of indeterminate bases per 100bp (default is 1, range 0..10)
		int Trim5,						// trim this number of 5' bases from input contigs sequences (default is 0, range 0..20)
		int Trim3,						// trim this number of 3' bases from input contigs sequences (default is 0, range 0..20)
		int MinSeqLen,		            // filter out contig sequences (after any trimming) which are less than this length (default is 50bp, range 30..10000)
		int MinCtgLen,					// filter out homozygotic region reduced contigs of less than this length
		char *pszCtgDescr,				// contig descriptor prefix
		int NumThreads,					// number of worker threads to use
    	int NumInputFiles,				// number of input file specs
		char *pszInfileSpecs[],			// names of input files
		char *pszOutFile);				// homozygotic region reduced contigs written to this file

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
	return _T("kangacrdx");
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
int Rslt = 0;   			// function result code >= 0 represents success, < 0 on failure
int Idx;

etPMode PMode;				// processing mode

int NumberOfProcessors;		// number of installed CPUs
int NumThreads;				// number of threads (0 defaults to number of CPUs)
bool bStrand;				// strand specific homozygous region reduction - homozygous regions between any two contigs must be in same orientation

int MaxNs;					// filter out input sequences having higher than this number of indeterminate bases per 100bp (default is 1, range 0..10)
int Trim5;					// trim this number of 5' bases from input sequences (default is 0, range 0..20)
int Trim3;					// trim this number of 3' bases from input sequences (default is 0, range 0..20)
int MinSeqLen;              // filter out input sequences (after any trimming) which are less than this length (default is 50bp, range 30..10000)
int MinCtgLen;				// filter out homozygotic region reduced contigs of less than this length
int MaxHomozySubs;			// characterise as homozygotic if substitution rate between regions <= this rate per 100bp
int MinHomozyLen;			// homozygotic regions to be at least this length
int MinHetrozyLen;			// island (marked as homozygotic either flank) hetrozygotic regions must be at least this length otherwise treat as homozygotic
char szCtgDescr[80];		// contig descriptor prefix

char szRsltsFile[_MAX_PATH];	//  homozygotic region reduced contigs written to this file

int NumInputFiles;			// number of input files
char *pszInfileSpecs[cRRMaxInFileSpecs];  // input read files

// command line args
struct arg_lit  *help    = arg_lit0("h","help",                 "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to screen and logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"diagnostics log file");

struct arg_int *pmode = arg_int0("m","mode","<int>",		    "Shared homozygous region identification processing sensitivity: 0 - standard sensitivity, 1 - more sensitive (slower), 2 - ultra sensitive (slowest), 3 - less sensitive (quicker)");

struct arg_file *infiles = arg_filen("i","in","<file>",1,cRRMaxInFileSpecs,"Load contig sequences from file(s), must be fasta or fastq, wildcards allowed");

struct arg_file *outfile = arg_file1("o","out","<file>",		"Output homozygous region reduced contigs to this file");

struct arg_int *maxns = arg_int0("n","indeterminates","<int>",  "filter out input sequences having higher than this percentage of indeterminate bases (default is 1, range 0..5)");
struct arg_int *trim5 = arg_int0("x","trim5","<int>",			"trim this number of 5' bases from input sequences (default is 0, range 0..20)");
struct arg_int *trim3 = arg_int0("X","trim3","<int>",			"trim this number of 3' bases from input sequences (default is 0, range 0..20)");
struct arg_int *minseqlen= arg_int0("l","minlen","<int>",       "filter out input sequences (after any trimming) which are less than this length (default is 50bp, range 30..10000)");
struct arg_int *minctglen= arg_int0("L","minctglen","<int>",    "filter out homozygous region reduced contigs which which are less than this length (default is 100bp, range 30..10000)");


struct arg_int *threads = arg_int0("T","threads","<int>",		"number of processing threads 0..128 (defaults to 0 which sets threads to number of CPU cores)");
struct arg_lit  *strand   = arg_lit0("S","strand",              "strand specific homozygous region reduction - identified homozygous regions between any two contigs must be same orientation");

struct arg_int  *maxhomozysubs = arg_int0("z","maxhomozysubs","<int>","characterise region as homozygous if differs by at most this base rate per 100 from any other region (default is 3%%, range 0..7%%");
struct arg_int  *minhomozylen = arg_int0("Z","minhomozylen","<int>","characterise region as homozygous if region at least this length (default is 75, range 30..10000)");
struct arg_int  *minhetrozylen = arg_int0("k","minhetrozylen","<int>","treat hetrozygotic region as homozygous if less than this length and when flanked by homozygous region (default is 30, range 10..10000)");

struct arg_str *ctgdescr = arg_str0("c","ctgdescr","<string>",	"contig identifer descriptor prefix (default is 'KHrdxCtg')");

struct arg_end *end = arg_end(20);

void *argtable[] = {help,version,FileLogLevel,LogFile,
					pmode,strand,maxns,trim5,trim3,minseqlen,maxhomozysubs,minhomozylen,minctglen,minhetrozylen,ctgdescr,infiles,outfile,
					threads,
					end};

char **pAllArgs;
int argerrors;
argerrors = CUtility::arg_parsefromfile(argc,(char **)argv,&pAllArgs);
if(argerrors >= 0)
	argerrors = arg_parse(argerrors,pAllArgs,argtable);

/* special case: '--help' takes precedence over error reporting */
if (help->count > 0)
        {
		printf("\n%s the K-mer Adaptive Next Generation Assembler Homozygous Region Reducer, Version %s\nOptions ---\n", gszProcName,cpszProgVer);
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
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Processing mode '-m%d' specified outside of range %d..%d\n",PMode,ePMdefault,(int)ePMplaceholder-1);
		exit(1);
		}

	bStrand = strand->count ? true : false;

	strncpy(szRsltsFile,outfile->filename[0],_MAX_PATH);
	szRsltsFile[_MAX_PATH-1] = '\0';

	if(!infiles->count)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: No input file(s) specified with with '-i<filespec>' option)");
		exit(1);
		}

	for(NumInputFiles=Idx=0;NumInputFiles < cRRMaxInFileSpecs && Idx < infiles->count; Idx++)
		{
		pszInfileSpecs[Idx] = NULL;
		if(pszInfileSpecs[NumInputFiles] == NULL)
			pszInfileSpecs[NumInputFiles] = new char [_MAX_PATH];
		strncpy(pszInfileSpecs[NumInputFiles],infiles->filename[Idx],_MAX_PATH);
		pszInfileSpecs[NumInputFiles][_MAX_PATH-1] = '\0';
		CUtility::TrimQuotedWhitespcExtd(pszInfileSpecs[NumInputFiles]);
		if(pszInfileSpecs[NumInputFiles][0] != '\0')
			NumInputFiles++;
		}

	if(!NumInputFiles)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: After removal of whitespace, no input file(s) specified with '-i<filespec>' option");
		exit(1);
		}


	MinCtgLen = minctglen->count ? minctglen->ival[0] : 100;
	if(MinCtgLen < 0 || MinCtgLen > 10000)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: minimum contig length '-L%d' must be in range 30..10000",MinCtgLen);
		exit(1);
		}

	MaxHomozySubs = maxhomozysubs->count ? maxhomozysubs->ival[0] : 3;
	if(MaxHomozySubs < 0 || MaxHomozySubs > 7)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: max homozygotic sub rate '-z%d' must be in range 0..7",MaxHomozySubs);
		exit(1);
		}

	MinHomozyLen = minhomozylen->count ? minhomozylen->ival[0] : 75;
	if(MinHomozyLen < 30 || MinHomozyLen > 10000)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: minimum homozygotic region length '-Z%d' must be in range 30..10000",MinHomozyLen);
		exit(1);
		}

	MinHetrozyLen = minhetrozylen->count ? minhetrozylen->ival[0] : 30;
	if(MinHetrozyLen < 30 || MinHetrozyLen > 10000)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: minimum hetrozygotic region length '-Z%d' must be in range 30..10000",MinHetrozyLen);
		exit(1);
		}

	Trim5 = trim5->count ? trim5->ival[0] : 0;
	if(Trim5 < 0 || Trim5 > 20)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: 5' trim '-x%d' must be in range 0..20 bases",Trim5);
		exit(1);
		}

	Trim3 = trim3->count ? trim3->ival[0] : 0;
	if(Trim3 < 0 || Trim3 > 20)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: 3' trim '-x%d' must be in range 0..20 bases",Trim3);
		exit(1);
		}

	MaxNs = maxns->count ? maxns->ival[0] : 1;
	if(MaxNs < 0 || MaxNs > 5)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: max percentage of indeterminate bases '-n%d' must be in range 0..5",MaxNs);
		exit(1);
		}

	MinSeqLen = minseqlen->count ? minseqlen->ival[0] : 50;
	if(MinSeqLen < 30 || MinSeqLen > 10000)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Minimum sequence length '-l%d' must be in range 30..10000 bases",MinSeqLen);
		exit(1);
		}

	if(ctgdescr->count)
		{
		strncpy(szCtgDescr,ctgdescr->sval[0],sizeof(szCtgDescr));
		szCtgDescr[sizeof(szCtgDescr)-1] = '\0';
		CUtility::TrimQuotedWhitespcExtd(szCtgDescr);
		CUtility::ReduceWhitespace(szCtgDescr);
		}
	else
		strcpy(szCtgDescr,"KHrdxCtg");

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

// show user current resource limits
#ifndef _WIN32
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Resources: %s",CUtility::ReportResourceLimits());
#endif

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing parameters:");
	const char *pszDescr;

	switch(PMode) {
		case ePMdefault:
			pszDescr = "Standard shared homozygous region identification sensitivity";
			break;
		case ePMMoreSens:
			pszDescr = "More sensitive shared homozygous region identification (slower)";
			break;
		case ePMUltraSens:
			pszDescr = "Ultra sensitive shared homozygous region identification (very slow)";
			break;
		case ePMLessSens:
		default:
			pszDescr = "Less sensitive shared homozygous region identification (quicker)";
			break;
		}

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Processing mode is : '%s'",pszDescr);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Strand specific homozygous region reduction: '%s'",bStrand ? "Yes" : "No");

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Accept input sequences if percentage of indeterminate bases (Ns) no more than : %d",MaxNs);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Trim 5' input sequences by : %d",Trim5);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Trim 3' input sequences by : %d",Trim3);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Accept input sequences (after any trim) if at least this length : %d",MinSeqLen);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Characterise as homozygotic if substitution rate between regions <= this rate per 100bp: %d%%",MaxHomozySubs);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Characterise as homozygotic if identified region at least this length: %d",MinHomozyLen);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Characterise as homozygotic if hetrozygotic region flanked by homozygotic region(s) less than this length: %d",MinHetrozyLen);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Prefix assembled contig descriptors with : '%s'",szCtgDescr);

	for(Idx=0; Idx < NumInputFiles; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Input sequences file (%d) : '%s'",Idx+1,pszInfileSpecs[Idx]);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Filter out homozygous region reduced contigs if length less than: %d",MinCtgLen);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Output homozygous region reduced contigs to file: '%s'",szRsltsFile);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"number of threads : %d",NumThreads);

#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	gStopWatch.Start();
	Rslt = Process(PMode,bStrand,MaxHomozySubs,MinHomozyLen,MinHetrozyLen,MaxNs,Trim5,Trim3,MinSeqLen,MinCtgLen,szCtgDescr,NumThreads,NumInputFiles,pszInfileSpecs,szRsltsFile);
	gStopWatch.Stop();
	Rslt = Rslt >=0 ? 0 : 1;
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Exit code: %d Total processing time: %s",Rslt,gStopWatch.Read());
	exit(Rslt);
	}
else
	{
	printf("\n%s the K-mer Adaptive Next Generation Assembler Homozygous Region Reducer, Version %s\n",gszProcName,cpszProgVer);
	arg_print_errors(stdout,end,gszProcName);
	arg_print_syntax(stdout,argtable,"\nUse '-h' to view option and parameter usage\n");
	exit(1);
	}
return 0;
}

int
Process(etPMode PMode,					// processing mode
		bool bStrand,					// strand specific homozygous region reduction - homozygous regions between any two contigs must be in same orientation
		int MaxHomozySubs,				// characterise as homozygotic if substitution rate between regions <= this rate per 100bp
		int MinHomozyLen,				// homozygotic regions to be at least this length
		int MinHetrozyLen,				// island (marked as homozygotic either flank) hetrozygotic regions must be at least this length otherwise treat as homozygotic
		int MaxNs,						// filter out input contigs having higher than this number of indeterminate bases per 100bp (default is 1, range 0..10)
		int Trim5,						// trim this number of 5' bases from input contigs sequences (default is 0, range 0..20)
		int Trim3,						// trim this number of 3' bases from input contigs sequences (default is 0, range 0..20)
		int MinSeqLen,		            // filter out contig sequences (after any trimming) which are less than this length (default is 50bp, range 30..10000)
		int MinCtgLen,					// filter out homozygotic region reduced contigs of less than this length
		char *pszCtgDescr,				// contig descriptor prefix
		int NumThreads,					// number of worker threads to use
    	int NumInputFiles,				// number of input file specs
		char *pszInfileSpecs[],			// names of input files
		char *pszOutFile)				// homozygotic region reduced contigs written to this file
{
int Rslt;
int Idx;
int NumInputFilesProcessed;
int TotNumContigsAccepted;
char *pszInFile;
CHomozyReduce *pHomozyReduce;
pHomozyReduce = new CHomozyReduce();
pHomozyReduce->SetNumThreads(NumThreads);
pHomozyReduce->Reset(false);
pHomozyReduce->SetCtgDescr(pszCtgDescr);

#ifdef CHECKXFORMID
UINT64 XFormID;
UINT32 Idy;
UINT32 FormID;
for(Idy = 0; Idy < 0x80000000; Idy++)
	{
	XFormID = pHomozyReduce->IDtoXForm(Idy);
	XFormID &= 0x0ffffffffff; // 5 bytes...
	FormID = pHomozyReduce->XFormToID(XFormID);
	if(FormID != Idy)
		printf("\n A problem...");
	}
#endif


NumInputFilesProcessed = 0;
TotNumContigsAccepted = 0;
CSimpleGlob glob(SG_GLOB_FULLSORT);

// now load the raw contig sequences applying any trimming and length filtering
for(Idx = 0; Idx < NumInputFiles; Idx++)
	{
	glob.Init();
	if(glob.Add(pszInfileSpecs[Idx]) < SG_SUCCESS)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to glob '%s",pszInfileSpecs[Idx]);
		return(eBSFerrOpnFile);	// treat as though unable to open file
		}

	if(glob.FileCount() <= 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to locate any source contig sequences file matching '%s",pszInfileSpecs[Idx]);
		continue;
		}
	Rslt = eBSFSuccess;
	for (int FileID = 0; Rslt >= eBSFSuccess &&  FileID < glob.FileCount(); ++FileID)
		{
		pszInFile = glob.File(FileID);
		NumInputFilesProcessed += 1;
		if((Rslt = pHomozyReduce->LoadContigs(MaxNs,Trim5,Trim3,MinSeqLen,NumInputFilesProcessed,pszInFile)) < eBSFSuccess)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: Load failed for source contig sequences file '%s'\n",pszInFile);
			pHomozyReduce->Reset(false);
			delete pHomozyReduce;
			return((teBSFrsltCodes)Rslt);
			}
		TotNumContigsAccepted += (int)Rslt;
		}
	}

if(NumInputFilesProcessed == 0)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: No contig sequence files were loaded - nothing to do\n");
	pHomozyReduce->Reset(true);
	delete pHomozyReduce;
	return(eBSFerrNoEntries);
	}

if(TotNumContigsAccepted < 1)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: No contig sequences were loaded - nothing to homozygosity reduce\n");
	pHomozyReduce->Reset(true);
	delete pHomozyReduce;
	return(eBSFerrNoEntries);
	}

if((Rslt = pHomozyReduce->ReduceHomozygosity(PMode,pszOutFile,bStrand,MaxHomozySubs,MinHomozyLen,MinHetrozyLen,MinCtgLen)) < eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to generate assembly");
	pHomozyReduce->Reset(false);
	}
else
	pHomozyReduce->Reset(true);

delete pHomozyReduce;
return(Rslt);
}



