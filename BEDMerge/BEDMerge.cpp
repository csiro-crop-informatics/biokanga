// BEDMerge.cpp : Defines the entry point for the console application.
//


#include "stdafx.h"
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#if _WIN32
#include "../libbiokanga/commhdrs.h"
#else
#include "../libbiokanga/commhdrs.h"
#endif

#include "./ChromFeatures.h"

const char *cpszProgVer = "1.1.2";		// increment with each release

const int cMinMergeLen  = 5;				// after merging, features must be of at least this length
const int cDfltMergeLen = 20;				// default merge feature length
const int cMaxMergeLen = 1024000;			// max merge feature length



const int cRsltsLineLen = 4096;			// buffer size for results loci lines

// processing modes
typedef enum eProcMode {
	ePMDefault,						// default processing is for strand independent union of features
	ePMStrandDep,					// strand dependent union of features
	ePMplaceholder					// used to set the enumeration range
} etProcMode;

int
Process(etProcMode ProcMode,	// processing mode
		char Strand,			// process for this strand only
		int MinLen,				// generated merged features must be at least this length
		int JoinLen,			// join features if separated by at most this length
		etBEDRegion Region,		// process for these regions only
		int NumIncludeChroms,		 // number of chromosome regular expressions to include
		char **ppszIncludeChroms,	 // array of include chromosome regular expressions
		int NumExcludeChroms,		 // number of chromosome expressions to exclude
		char **ppszExcludeChroms,	 // array of exclude chromosome regular expressions
		int NumInFiles,			// number of input files containing features to merge
		char *pszInFiles[],		// BED files containing elements to be combined
		char *pszOutFile);		// write combined elements into this BED file

int TrimQuotes(char *pszTxt);
char *ProcMode2Txt(etProcMode ProcMode);
char *Region2Txt(etBEDRegion Region);

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
	return _T("BEDMerge");
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

int Rslt;
int Idx;
int iProcMode;
char Strand;							// process features on this strand
int MinLen;
int JoinLen;
int LenFileList;
int NumIncludeChroms;
char *pszIncludeChroms[cMaxIncludeChroms];
int NumExcludeChroms;
char *pszExcludeChroms[cMaxExcludeChroms];
etBEDRegion Region;				// process for this functional region only
int NumInFiles;							// number of control input files
char *pszInFileSpecs[cMaxNumBedFiles];  // input (wildcards allowed) BED files

char szOutFile[_MAX_PATH];	// write merged to this file

// command line args
struct arg_lit  *help    = arg_lit0("h","help",                 "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to screen and logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"diagnostics log file");
struct arg_int  *ProcMode = arg_int0("m","mode","<int>",		"merge processing mode: 0 - Strand independent union, 1 - Strand dependent union");
struct arg_int  *strand = arg_int0("s","strand","<int>",		"filter for this strand: 0 - any, 1 - Watson '+', 2 - Crick '-' (default is any)");
struct arg_int *region = arg_int0("r","genomicregion","<int>",	"Retain annotated regions 0:ALL,1:Intergenic,2:Exons,3:Introns,4:CDS,5:UTRs,6:5'UTR,7:3'UTR (default = ALL)");
struct arg_int  *minlen = arg_int0("l","minlen","<int>",		"generated merged features must be of at least this length (default is 20)");
struct arg_int  *joinlen = arg_int0("j","joinlen","<int>",		"merge features which are separated by at most this length, 0 if no merging (default is 1)");
struct arg_str  *excludechroms = arg_strn("Z","chromexclude","<string>",0,cMaxExcludeChroms,"regular expressions defining chromosomes to always exclude");
struct arg_str  *includechroms = arg_strn("z","chromeinclude","<string>",0,cMaxIncludeChroms,"regular expressions defining chromosomes to explicitly include if not already excluded");
struct arg_file *infiles = arg_filen("i","srcfiles","<file>",0,  cMaxNumBedFiles,"merge features contained in these BED files (wildcards alllowed)");
struct arg_file *OutFile = arg_file1("o",NULL,"<file>",			"output merged features to this BED file");
struct arg_end *end = arg_end(20);

void *argtable[] = {help,version,FileLogLevel,LogFile,
					ProcMode,strand,minlen,joinlen,region,excludechroms,includechroms,infiles,OutFile,
					end};

char **pAllArgs;
int argerrors;
argerrors = CUtility::arg_parsefromfile(argc,(char **)argv,&pAllArgs);
if(argerrors >= 0)
	argerrors = arg_parse(argerrors,pAllArgs,argtable);

    /* special case: '--help' takes precedence over error reporting */
if (help->count > 0)
        {
		printf("\n%s BED Merge Blocks, Version %s\nOptions ---\n", gszProcName,cpszProgVer);
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

	iProcMode = ProcMode->count ? ProcMode->ival[0] : ePMDefault;
	if(iProcMode < ePMDefault || iProcMode >= ePMplaceholder)
		{
		printf("Error: Processing mode '-p%d' is not in range 0..%d",iProcMode,ePMplaceholder-1);
		exit(1);
		}

	MinLen = minlen->count ? minlen->ival[0] : cDfltMergeLen;
	if(MinLen < cMinMergeLen || MinLen > cMaxMergeLen)
		{
		printf("Error: minimum length '-l%d' is not in range %d..%d",MinLen,cMinMergeLen,cMaxMergeLen);
		exit(1);
		}


	JoinLen = joinlen->count ? joinlen->ival[0] : 1;
	if(JoinLen < 0 || JoinLen > cMaxMergeLen)
		{
		printf("Error: Join length '-j%d' is not in range 0..%d",JoinLen,cMaxMergeLen);
		exit(1);
		}


	Strand = strand->count ? strand->ival[0] : 0;
	if(Strand < 0 || Strand > 2)
		{
		printf("\nError: Strand '-s%d' specified outside of range 0..2",Strand);
		exit(1);
		}

	switch(Strand) {
		case 1: Strand = (int)'+'; break;
		case 2: Strand = (int)'-'; break;
		case 0: Strand = (int)'*'; break;
		}

	Region = (etBEDRegion)(region->count ? region->ival[0] : eMEGRAny);	// default as being any region
	if(Region < eMEGRAny || Region > eMEG3UTR)
		{
		printf("\nSpecified region '-g%d' outside of range 0..%d",Region,eMEG3UTR);
		exit(1);
		}

	NumIncludeChroms = includechroms->count;
	for(Idx=0;Idx < includechroms->count; Idx++)
		{
		LenFileList = (int)strlen(includechroms->sval[Idx]);
		pszIncludeChroms[Idx] = new char [LenFileList+1];
		strcpy(pszIncludeChroms[Idx],includechroms->sval[Idx]);
		TrimQuotes(pszIncludeChroms[Idx]);
		}

	NumExcludeChroms = excludechroms->count;
	for(Idx=0;Idx < excludechroms->count; Idx++)
		{
		LenFileList = (int)strlen(excludechroms->sval[Idx]);
		pszExcludeChroms[Idx] = new char [LenFileList+1];
		strcpy(pszExcludeChroms[Idx],excludechroms->sval[Idx]);
		TrimQuotes(pszExcludeChroms[Idx]);
		}

	if(!infiles->count)
		{
		printf("\nError: No input file(s) specified with with '-i<filespec>' option)");
		exit(1);
		}

	for(NumInFiles=Idx=0;NumInFiles < cMaxNumBedFiles && Idx < infiles->count; Idx++)
		{
		pszInFileSpecs[Idx] = NULL;
		if(pszInFileSpecs[NumInFiles] == NULL)
			pszInFileSpecs[NumInFiles] = new char [_MAX_PATH];
		strncpy(pszInFileSpecs[NumInFiles],infiles->filename[Idx],_MAX_PATH);
		pszInFileSpecs[NumInFiles][_MAX_PATH-1] = '\0';
		CUtility::TrimQuotedWhitespcExtd(pszInFileSpecs[NumInFiles]);
		if(pszInFileSpecs[NumInFiles][0] != '\0')
			NumInFiles++;
		}

	if(!NumInFiles)
		{
		printf("\nError: After removal of whitespace, no input file(s) specified with '-i<filespec>' option)");
		exit(1);
		}

	strncpy(szOutFile,OutFile->filename[0],_MAX_PATH);
	szOutFile[_MAX_PATH-1] = '\0';


			// now that command parameters have been parsed then initialise diagnostics log system
	if(!gDiagnostics.Open(szLogFile,(etDiagLevel)iScreenLogLevel,(etDiagLevel)iFileLogLevel,true))
		{
		printf("\nError: Unable to start diagnostics subsystem.");
		if(szLogFile[0] != '\0')
			printf(" Most likely cause is that logfile '%s' can't be opened/created",szLogFile);
		exit(1);
		}

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Version: %s Processing parameters:",cpszProgVer);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Processing Mode: %d (%s)",iProcMode,ProcMode2Txt((etProcMode)iProcMode));
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"process for this strand only: '%c'",(char)Strand);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"generated merged features will be at least this length: %d",MinLen);
	if(JoinLen == 0)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"keep separate features, no merging");
	else
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"merge features separated by upto: %d",JoinLen);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Retain Region: %s",Region2Txt((etBEDRegion)Region));
	for(Idx = 0; Idx < NumIncludeChroms; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"reg expressions defining chroms to include: '%s'",pszIncludeChroms[Idx]);
	for(Idx = 0; Idx < NumExcludeChroms; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"reg expressions defining chroms to exclude: '%s'",pszExcludeChroms[Idx]); 
	for(Idx=0;Idx < NumInFiles; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"merge features from this input BED file (%d): '%s'",Idx+1,pszInFileSpecs[Idx]);		
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Output merged features into BED file: '%s'",szOutFile);

	gStopWatch.Start();
#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif	
	Rslt = Process((etProcMode)iProcMode,Strand,MinLen,JoinLen,Region,
			NumIncludeChroms,	// number of chromosome regular expressions to include
			pszIncludeChroms,	// array of include chromosome regular expressions
			NumExcludeChroms,	// number of chromosome expressions to exclude
			pszExcludeChroms,	// array of exclude chromosome regular expressions
		NumInFiles,pszInFileSpecs,szOutFile);
	gStopWatch.Stop();
	Rslt = Rslt >=0 ? 0 : 1;
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Exit code: %d Total processing time: %s",Rslt,gStopWatch.Read());
	exit(Rslt);
	}
else
	{
	printf("\n%s BED Merge Blocks, Version %s\n",gszProcName,cpszProgVer);
	arg_print_errors(stdout,end,gszProcName);
	arg_print_syntax(stdout,argtable,"\nUse '-h' to view option and parameter usage\n");
	exit(1);
	}
}

char *
Region2Txt(etBEDRegion Region)
{
switch(Region) {
	case eMEGRAny:		// process any region
		return((char *)"All");

	case eMEGRIntergenic:	// only process intergenic
		return((char *)"Intergenic");

	case eMEGRExons:	// only process exons
		return((char *)"EXONS");

	case eMEGRIntrons:	// only process introns
		return((char *)"INTRONS");

	case eMEGRCDS:		// only process CDSs
		return((char *)"CDS");

	case eMEGUTR:		// only process UTRs
		return((char *)"UTR");

	case eMEG5UTR:		// only process 5'UTRs
		return((char *)"5'UTR");

	case eMEG3UTR:		// only process 3'UTRs
		return((char *)"3'UTR");

	default:
		break;
	}
return((char *)"Unsupported");
}

char *
ProcMode2Txt(etProcMode ProcMode)
{
switch(ProcMode) {
	case ePMDefault:
		return((char *)"Strand independent union of overlapping features");
	case ePMStrandDep:
		return((char *)"Strand dependent union of overlapping features");
	default:
		break;
	}
return((char *)"Unrecognised");
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

CChromFeatures *m_pChromFeatures;

void
Init(void)
{
m_pChromFeatures = NULL;
}

void
Reset(void)
{
if(m_pChromFeatures != NULL)
	{
	delete m_pChromFeatures;
	m_pChromFeatures = NULL;
	}
}

int
Process(etProcMode ProcMode,	// processing mode
		char Strand,			// process for this strand only
		int MinLen,				// generated merged features must be at least this length
		int JoinLen,			// join features if separated by at most this length
		etBEDRegion Region,		// process for these regions only
		int NumIncludeChroms,		 // number of chromosome regular expressions to include
		char **ppszIncludeChroms,	 // array of include chromosome regular expressions
		int NumExcludeChroms,		 // number of chromosome expressions to exclude
		char **ppszExcludeChroms,	 // array of exclude chromosome regular expressions
		int NumInFiles,			// number of input files containing features to merge
		char *pszInFileSpecs[],	// BED files containing elements to be combined
		char *pszOutFile)		// write combined elements into this BED file
{
teBSFrsltCodes Rslt;
char *pszInfile;
int NumBEDsLoaded = 0;
int Idx;			// general processing iteration index

Init();

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Output merge file created/truncated: '%s'",pszOutFile);

m_pChromFeatures = new CChromFeatures;

if(NumIncludeChroms > 0 || NumExcludeChroms > 0)
	{
	if((Rslt=(teBSFrsltCodes)m_pChromFeatures->SetChromFilters(NumIncludeChroms,ppszIncludeChroms,NumExcludeChroms,ppszExcludeChroms)) != eBSFSuccess)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to initialise chrom include/exclude REs");
		Reset();
		return(Rslt);
		}
	}

// process the control files
CSimpleGlob glob(SG_GLOB_FULLSORT);
for(Idx = 0; Idx < NumInFiles; Idx++)
	{
	glob.Init();
	if(glob.Add((const char *)pszInFileSpecs[Idx]) < SG_SUCCESS)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to glob '%s",pszInFileSpecs[Idx]);
		Reset();
		return(eBSFerrOpnFile);	// treat as though unable to open file
		}

	if(glob.FileCount() <= 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to locate any source file matching '%s",pszInFileSpecs[Idx]);
		continue;
		}

	Rslt = eBSFSuccess;
	for (int FileID = 0; Rslt >= eBSFSuccess &&  FileID < glob.FileCount(); ++FileID)
		{
		pszInfile = glob.File(FileID);
		NumBEDsLoaded += 1;
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing input BED file %d '%s'",NumBEDsLoaded,pszInfile);
		
		Rslt = m_pChromFeatures->LoadBED(Region,Strand,pszInfile);
		if(Rslt != eBSFSuccess)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: Load failed for input BED file '%s'",pszInfile);
			Reset();
			return(Rslt);
			}
		}
	}
if(!NumBEDsLoaded)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: Unable to load any BED input files");
	Reset();
	return(eBSFerrOpnFile);
	}

Rslt = (teBSFrsltCodes)m_pChromFeatures->MergeFeatures(pszOutFile,ProcMode == ePMStrandDep ? true : false,JoinLen,MinLen);
if(Rslt >= 0)
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Process: %d merged features written to BED file '%s'",(int)Rslt,pszOutFile);
else
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: Merge failed");
Reset();
return(Rslt);
}


