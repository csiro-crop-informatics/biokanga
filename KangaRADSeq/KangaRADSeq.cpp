// KangaRADSeq.cpp : Defines the entry point for the console application.
// Processes either single or paired end RADseq sequenced reads

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

#include "StackSeqs.h"

const char *cpszProgVer = "0.0.2";		// increment with each release

CStopWatch gStopWatch;
CDiagnostics gDiagnostics;				// for writing diagnostics messages to log file
char gszProcName[_MAX_FNAME];			// process name



int Process(etPMode PMode,				// processing sensitivity mode
	int MaxNs,							// filter out input sequences having higher than this number of indeterminate bases per 100bp (default is 1, range 0..10)
	int Trim5,							// trim this number of 5' bases from input sequences (default is 0, range 0..20)
	int Trim3,							// trim this number of 3' bases from input sequences (default is 0, range 0..20)
	int MinSeqLen,						// filter out input sequences (after any trimming) which are less than this length (default is 50bp, range 30..10000)

	int P1StackEnd,						// P1 stack maximum end float (default is 1, range 0..10)");
	int P1StackDepth,					// P1 stack minimum depth (default is 10, range 5..1000)");
	int P1StackSubRate,					// P1 stack maximum substitution rate (default is 1%, range 0..10%)");

	int P2MinOvrl,						// P2 minimum read overlap (default is 30% of read length, range 10..100)");
	int P2MaxOvrlSubRate,				// P2 maximum read overlap substitution rate (default is 5%%, range 0..10%%)");
	char *pszCtgDescr,					// generated contig descriptor prefix 
	int NumThreads,						// number of worker threads to use
	int NumInputP1Files,				// number of input P1 file specs
	char *pszInP1Files[],				// names of input files
	int NumInputP2Files,				// number of input P2 files
	char *pszInP2Files[],				// input P2 files
	char *pszOutCtgsFile,				// assembled contigs written to this file
	char *pszOutVCFfile);				// output Variant Call Format (VCF 4.1) to this file

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
	return _T("kangaRADSeq");
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

etPMode PMode;				// processing sensitivity mode
int NumberOfProcessors;		// number of installed CPUs
int NumThreads;				// number of threads (0 defaults to number of CPUs)

int MaxNs;					// filter out input sequences having higher than this number of indeterminate bases per 100bp (default is 1, range 0..10)
int Trim5;					// trim this number of 5' bases from input sequences (default is 0, range 0..20)
int Trim3;					// trim this number of 3' bases from input sequences (default is 0, range 0..20)
int MinSeqLen;              // filter out input sequences (after any trimming) which are less than this length (default is 50bp, range 30..10000)

int P1StackEnd;				// P1 stack maximum end float (default is 1, range 0..10)");
int P1StackDepth;			// P1 stack minimum depth (default is 10, range 5..1000)");
int P1StackSubRate;			// P1 stack maximum substitution rate (default is 1%, range 0..10%)");

int P2MinOvrl;				// P2 minimum read overlap (default is 30% of read length, range 10..100)");
int P2MaxOvrlSubRate;		// P2 maximum read overlap substitution rate (default is 5%%, range 0..10%%)");


char szCtgDescr[80];			// generated contig descriptor prefix 
char szOutCtgsFile[_MAX_PATH];	// assembled contigs written to this file
char szOutVCFfile[_MAX_PATH];	// output Variant Call Format (VCF 4.1) to this file

int NumInputP1Files;			// number of input P1 files
char *pszInP1Files[cRRMaxInFileSpecs];  // input P1 files

int NumInputP2Files;			// number of input P2 files
char *pszInP2Files[cRRMaxInFileSpecs];  // input P2 files

// command line args
struct arg_lit  *help    = arg_lit0("h","help",                 "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to screen and logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"diagnostics log file");

struct arg_int *pmode = arg_int0("m","mode","<int>",		    "RAD processing sensitivity: 0 - standard sensitivity, 1 - more sensitive (slower), 2 - ultra sensitive (slowest), 3 - less sensitive (quicker)");

struct arg_int *p1stackend   = arg_int0("z","p1stackend","<int>", "P1 stack maximum end float (default is 5% of mean P1 read length, range 0..20bp)");
struct arg_int *p1stackdepth = arg_int0("Z","p1stackdepth","<int>", "P1 stack minimum depth (default is 10, range 5..1000)");
struct arg_int *p1stacksubrate = arg_int0("s","p1stacksubrate","<int>", "P1 stack maximum substitution rate (default is 1%%, range 0..10%%)");

struct arg_int *p2minovrl = arg_int0("y","p2minovrl","<int>",    "P2 minimum read overlap (default is 50% of mean P2 read length, range 10..100bp)");
struct arg_int *p2maxovrlsubrate = arg_int0("S","p2maxovrlsubrate","<int>",    "P2 maximum read overlap substitution rate (default is 1%%, range 0..10%%)");


struct arg_file *inp1files = arg_filen("i","in","<file>",1,cRRMaxInFileSpecs,"Load P1 fasta/fastq sequences from file(s) wildcards allowed");

struct arg_file *inp2files = arg_filen("I","in","<file>",1,cRRMaxInFileSpecs,"Load P2 fasta/fastq sequences from file(s), wildcards allowed");

struct arg_file *outctgsfile = arg_file1("o","out","<file>",	"Output assembled contigs to this file");
struct arg_file *outvcffile = arg_file0("O","variants","<file>",	"Output Variant Call Format (VCF 4.1) to this file");

struct arg_int *maxns = arg_int0("n","indeterminates","<int>",  "filter out input sequences having higher than this percentage of indeterminate bases (default is 1, range 0..5)");
struct arg_int *trim5 = arg_int0("x","trim5","<int>",			"trim this number of 5' bases from input sequences (default is 0, range 0..20)");
struct arg_int *trim3 = arg_int0("X","trim3","<int>",			"trim this number of 3' bases from input sequences (default is 0, range 0..20)");
struct arg_int *minseqlen= arg_int0("l","minlen","<int>",       "filter out input sequences (after any trimming) which are less than this length (default is 50bp, range 30..1000)");
struct arg_int *minctglen= arg_int0("L","minctglen","<int>",    "filter out assembled contigs which which are less than this length (default is 100bp, range 30..1000)");

struct arg_int *threads = arg_int0("T","threads","<int>",		"number of processing threads 0..128 (defaults to 0 which sets threads to number of CPU cores)");

struct arg_str *ctgdescr = arg_str0("c","ctgdescr","<string>",	"contig identifer descriptor prefix (default is 'KRADCtg')");


struct arg_end *end = arg_end(20);

void *argtable[] = {help,version,FileLogLevel,LogFile,
					pmode,maxns,trim5,trim3,minseqlen,p1stackend,p1stackdepth,p1stacksubrate,p2minovrl,p2maxovrlsubrate,minctglen,inp1files,inp2files,outctgsfile,outvcffile,ctgdescr,
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
		printf("\n%s the K-mer Adaptive Next Generation Assembler, Version %s\nOptions ---\n", gszProcName,cpszProgVer);
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


	if(ctgdescr->count)
		{
		strncpy(szCtgDescr,ctgdescr->sval[0],sizeof(szCtgDescr));
		szCtgDescr[sizeof(szCtgDescr)-1] = '\0';
		CUtility::TrimQuotedWhitespcExtd(szCtgDescr);
		CUtility::ReduceWhitespace(szCtgDescr);
		}
	else
		strcpy(szCtgDescr,"KDNACtg");

	strncpy(szOutCtgsFile,outctgsfile->filename[0],_MAX_PATH);
	szOutCtgsFile[_MAX_PATH-1] = '\0';

	if(outvcffile->count)
		{
		strncpy(szOutVCFfile,outvcffile->filename[0],_MAX_PATH);
		szOutCtgsFile[_MAX_PATH-1] = '\0';
		}
	else
		szOutCtgsFile[0] = '\0';

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


	if(inp2files->count)
		{
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
		}
	else
		{
		NumInputP2Files = 0;
		pszInP2Files[0] = NULL;
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
	if(MinSeqLen < 30 || MinSeqLen > 1000)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Minimum sequence length '-l%d' must be in range 30..1000 bases",MinSeqLen);
		exit(1);
		}

	P1StackEnd = p1stackend->count ? p1stackend->ival[0] : cDfltP1StackEnd;
	if(P1StackEnd != cDfltP1StackEnd && (P1StackEnd < cMinP1StackEnd || P1StackEnd > cMaxP1StackEnd))
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: P1 stack maximum end float '-z%d' must be in range %d..%d bases",P1StackEnd,cMinP1StackEnd,cMaxP1StackEnd);
		exit(1);
		}

	P1StackDepth = p1stackdepth->count ? p1stackdepth->ival[0] : cDfltP1StackDepth;
	if(P1StackDepth < cMinP1StackDepth || P1StackDepth > cMaxP1StackDepth)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: P1 stack minimum depth '-z%d' must be in range %d..%d bases",P1StackDepth,cMinP1StackDepth,cMaxP1StackDepth);
		exit(1);
		}

	P1StackSubRate = p1stacksubrate->count ? p1stacksubrate->ival[0] : cDfltP1StackSubRate;
	if(P1StackSubRate < cMinP1StackSubRate || P1StackSubRate > cMaxP1StackSubRate)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: P1 stack maximum substitution rate '-z%d' must be in range %d..%d%%",P1StackSubRate,cMinP1StackSubRate,cMaxP1StackSubRate);
		exit(1);
		}


	if(NumInputP2Files > 0)
		{
		P2MinOvrl = p2minovrl->count ? p2minovrl->ival[0] : cDfltP2MinOvrl;
		if(P2MinOvrl != -1 && (P2MinOvrl < cMinP2MinOvrl || P2MinOvrl > cMaxP2MinOvrl))
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: P2 minimum read overlap '-z%d' must be in range %d..%d bases",P2MinOvrl,cMinP2MinOvrl,cMaxP2MinOvrl);
			exit(1);
			}

		P2MaxOvrlSubRate = p2maxovrlsubrate->count ? p2maxovrlsubrate->ival[0] : cDfltP2MaxOvrlSubRate;
		if(P2MaxOvrlSubRate < cMinP2MaxOvrlSubRate || P2MaxOvrlSubRate > cMaxP2MaxOvrlSubRate)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: P2 maximum read overlap substitution rate '-z%d' must be in range %d..%d%%",P2MaxOvrlSubRate,cMinP2MaxOvrlSubRate,cMaxP2MaxOvrlSubRate);
			exit(1);
			}
		}
	else
		{
		P2MinOvrl = 0;
		P2MaxOvrlSubRate = 0;
		}

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
		case ePMDefault:
			pszDescr = "Standard sensitivity";
			break;
		case ePMMoreSens:
			pszDescr = "More sensitive (slower)";
			break;
		case ePMUltraSens:
			pszDescr = "Ultra sensitive (very slow)";
			break;
		case ePMLessSens:
		default:
			pszDescr = "Less sensitive (quicker)";
			break;
		}

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Processing mode is : '%s'",pszDescr);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Accept input sequences if percentage of indeterminate bases (Ns) no more than : %d",MaxNs);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Trim 5' input sequences by : %d",Trim5);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Trim 3' input sequences by : %d",Trim3);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Accept input sequences (after any trim) if at least this length : %d",MinSeqLen);
	
	if(P1StackEnd == -1)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"P1 stacks built with maximum end float: 5%% of mean P1 read lengths");
	else
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"P1 stacks built with maximum end float: %d bases",P1StackEnd);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"P1 stacks built with stacked read mismatch rate at most : %d%%",P1StackSubRate);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"P1 stacks to contain at least this number of stacked reads : %d",P1StackDepth);

	if(NumInputP2Files)
		{
		if(P2MinOvrl == -1)
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"P2 contigs built with minimum read overlaps: 30%% of mean P2 read lengths");
		else
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"P2 contigs built with minimum read overlaps: %d bases",P2MinOvrl);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"P2 contigs built with overlaying read mismatch rate at most : %d%%",P2MaxOvrlSubRate);
		}

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Prefix assembled contig descriptors with : '%s'",szCtgDescr);

	for(Idx=0; Idx < NumInputP1Files; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Input P1 sequences file (%d) : '%s'",Idx+1,pszInP1Files[Idx]);

	for(Idx=0; Idx < NumInputP2Files; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Input P2 sequences file (%d) : '%s'",Idx+1,pszInP2Files[Idx]);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Output assembled contigs to file: '%s'",szOutCtgsFile);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"number of threads : %d",NumThreads);

#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	gStopWatch.Start();
    Rslt = Process(PMode,				// processing sensitivity mode
			MaxNs,						// filter out input sequences having higher than this number of indeterminate bases per 100bp (default is 1, range 0..10)
			Trim5,						// trim this number of 5' bases from input sequences (default is 0, range 0..20)
			Trim3,						// trim this number of 3' bases from input sequences (default is 0, range 0..20)
			MinSeqLen,					// filter out input sequences (after any trimming) which are less than this length (default is 50bp, range 30..10000)
			P1StackEnd,					// P1 stack maximum end float (default is 1, range 0..10)");
			P1StackDepth,				// P1 stack minimum depth (default is 10, range 5..1000)");
			P1StackSubRate,				// P1 stack maximum substitution rate (default is 1%, range 0..10%)");
			P2MinOvrl,					// P2 minimum read overlap (default is 30% of read length, range 10..100)");
			P2MaxOvrlSubRate,			// P2 maximum read overlap substitution rate (default is 5%%, range 0..10%%)");
			szCtgDescr,					// generated contig descriptor prefix 
			NumThreads,					// number of worker threads to use
			NumInputP1Files,			// number of input P1 file specs
			pszInP1Files,				// names of input files
			NumInputP2Files,			// number of input P2 files
			pszInP2Files,				// input P2 files
			szOutCtgsFile,				// assembled contigs written to this file
			szOutVCFfile);				// output Variant Call Format (VCF 4.1) to this file

	gStopWatch.Stop();
	Rslt = Rslt >=0 ? 0 : 1;
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Exit code: %d Total processing time: %s",Rslt,gStopWatch.Read());
	exit(Rslt);
	}
else
	{
	printf("\n%s the K-mer Adaptive Next Generation Assembler, Version %s\n",gszProcName,cpszProgVer);
	arg_print_errors(stdout,end,gszProcName);
	arg_print_syntax(stdout,argtable,"\nUse '-h' to view option and parameter usage\n");
	exit(1);
	}
return 0;
}

// basic assumption is that the paired ends are in the same order in their respective paired end files
// so when reads are loaded from file P1 and discarded for whatever reason, then the corresponding read can be discarded from P2
int Process(etPMode PMode,				// processing sensitivity mode
	int MaxNs,							// filter out input sequences having higher than this number of indeterminate bases per 100bp (default is 1, range 0..10)
	int Trim5,							// trim this number of 5' bases from input sequences (default is 0, range 0..20)
	int Trim3,							// trim this number of 3' bases from input sequences (default is 0, range 0..20)
	int MinSeqLen,						// filter out input sequences (after any trimming) which are less than this length (default is 50bp, range 30..10000)

	int P1StackEnd,						// P1 stack maximum end float (default is 1, range 0..10)");
	int P1StackDepth,					// P1 stack minimum depth (default is 10, range 5..1000)");
	int P1StackSubRate,					// P1 stack maximum substitution rate (default is 1%, range 0..10%)");

	int P2MinOvrl,						// P2 minimum read overlap (default is 30% of read length, range 10..100)");
	int P2MaxOvrlSubRate,				// P2 maximum read overlap substitution rate (default is 5%%, range 0..10%%)");
	char *pszCtgDescr,					// generated contig descriptor prefix 
	int NumThreads,						// number of worker threads to use
	int NumInputP1Files,				// number of input P1 file specs
	char *pszInP1Files[],				// names of input files
	int NumInputP2Files,				// number of input P2 files
	char *pszInP2Files[],				// input P2 files
	char *pszOutCtgsFile,				// assembled contigs written to this file
	char *pszOutVCFfile)				// output Variant Call Format (VCF 4.1) to this file
{
int Rslt;
int Idx;
int NumInputFilesProcessed;
UINT32 TotNumReadsAccepted;
int NumP1Files;
int NumP2Files;

char *pszP1InFile;
char *pszP2InFile;
CStackSeqs *pStackSeqs;
pStackSeqs = new CStackSeqs();
pStackSeqs->Init();
pStackSeqs->SetNumThreads(NumThreads);
pStackSeqs->Reset(false);
pStackSeqs->SetCtgDescr(pszCtgDescr);

NumInputFilesProcessed = 0;
TotNumReadsAccepted = 0;
CSimpleGlob P1Glob(SG_GLOB_FULLSORT);
CSimpleGlob P2Glob(SG_GLOB_FULLSORT);

// now load the raw read sequences applying any trimming and colorspace processing
P1Glob.Init();
for(Idx = 0; Idx < NumInputP1Files; Idx++)
	{
	if(P1Glob.Add(pszInP1Files[Idx]) < SG_SUCCESS)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to glob P1 '%s",pszInP1Files[Idx]);
		return(eBSFerrOpnFile);	// treat as though unable to open file
		}
	}

if((NumP1Files = P1Glob.FileCount()) <= 0)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to locate any source reads file matching '%s",pszInP1Files[Idx]);
	return(eBSFerrOpnFile);	// treat as though unable to open file
	}

if(NumInputP2Files)
	{
	P2Glob.Init();
	for(Idx = 0; Idx < NumInputP2Files; Idx++)
		{
		if(P2Glob.Add(pszInP2Files[Idx]) < SG_SUCCESS)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to glob P2 '%s",pszInP2Files[Idx]);
			return(eBSFerrOpnFile);	// treat as though unable to open file
			}
		}

	if((NumP2Files = P2Glob.FileCount()) <= 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to locate any source reads file matching '%s",pszInP2Files[Idx]);
		return(eBSFerrOpnFile);	// treat as though unable to open file
		}

	if(NumP1Files != NumP2Files)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expected same number of P1 (%d) and P2 (%d) files",NumP1Files,NumP2Files);
		return(eBSFerrOpnFile);	// treat as though unable to open file
		}
	pStackSeqs->SetPairedEndProc(true);
	}
else
	pStackSeqs->SetPairedEndProc(false);


Rslt = eBSFSuccess;
for (int FileID = 0; Rslt >= eBSFSuccess &&  FileID < NumP1Files; ++FileID)
	{
	pszP1InFile = P1Glob.File(FileID);
	if(NumP2Files)
		pszP2InFile = P2Glob.File(FileID);
	else
		pszP2InFile = NULL;
	NumInputFilesProcessed += 1;
	if((Rslt = pStackSeqs->LoadRawReads(MaxNs,Trim5,Trim3,MinSeqLen,NumInputFilesProcessed,pszP1InFile,pszP2InFile)) < eBSFSuccess)
		{
		if(NumP2Files)
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: Load failed for input sequences file '%s'\n",pszP1InFile);
		else
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: Load failed for input sequences files P1 '%s' and P2 '%s'\n",pszP1InFile,pszP2InFile);
		pStackSeqs->Reset(false);
		delete pStackSeqs;
		return((teBSFrsltCodes)Rslt);
		}
	TotNumReadsAccepted += (int)Rslt;
	}

if(NumInputFilesProcessed == 0)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: No sequence files were loaded - nothing to assemble\n");
	pStackSeqs->Reset(true);
	delete pStackSeqs;
	return(eBSFerrNoEntries);
	}

if(TotNumReadsAccepted < 1)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: No sequence reads were accepted - nothing to assemble\n");
	pStackSeqs->Reset(true);
	delete pStackSeqs;
	return(eBSFerrNoEntries);
	}

if((Rslt = pStackSeqs->Process(PMode,						// processing sensitivity mode
							P1StackEnd,						// P1 stack maximum end float (default is 1, range 0..10)");
							P1StackDepth,					// P1 stack minimum depth (default is 10, range 5..1000)");
							P1StackSubRate,					// P1 stack maximum substitution rate (default is 1%, range 0..10%)");
							P2MinOvrl,						// P2 minimum read overlap (default is 30% of read length, range 10..100)");
							P2MaxOvrlSubRate,				// P2 maximum read overlap substitution rate (default is 5%%, range 0..10%%)");
							pszCtgDescr,					// generated contig descriptor prefix 
							pszOutCtgsFile,					// assembled contigs written to this file
							pszOutVCFfile,					// output Variant Call Format (VCF 4.1) to this file
							NumThreads)) < eBSFSuccess)		// max number of worker threads to use
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to generate assembly");
	pStackSeqs->Reset(false);
	}
else
	pStackSeqs->Reset(true);

delete pStackSeqs;
return(Rslt);
}