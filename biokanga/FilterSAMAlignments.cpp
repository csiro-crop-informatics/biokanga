/*
 * CSIRO Open Source Software License Agreement (GPLv3)
 * Copyright (c) 2017, Commonwealth Scientific and Industrial Research Organisation (CSIRO) ABN 41 687 119 230.
 * See LICENSE for the complete license information (https://github.com/csiro-crop-informatics/biokanga/LICENSE)
 * Contact: Alex Whan <alex.whan@csiro.au>
 */

// FilterSAMAlignments.cpp : Defines the entry point for the console application.

#include "stdafx.h"
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "../libbiokanga/commhdrs.h"

#include "biokanga.h"

#include "FilterSAMAlignments.h"


int Process(int NumIncludeChroms,		// number of retained chromosomes regular expressions
	char **ppszIncludeChroms,	// array of include chromosome regular expressions
	int NumExcludeChroms,		// number of chromosome expressions to explicitly exclude
	char **ppszExcludeChroms,	// array of exclude chromosome regular expressions
	char *pszInFile,			// input file containing alignments to be filtered
	char *pszOutFile);			// write filtered alignments to this output file)

int TrimREQuotes(char *pszTxt);

#ifdef _WIN32
int FilterSAMAlignments(int argc, char* argv[])
{
// determine my process name
_splitpath(argv[0],NULL,NULL,gszProcName,NULL);
#else
int
FilterSAMAlignments(int argc, char** argv)
{
// determine my process name
CUtility::splitpath((char *)argv[0],NULL,gszProcName);
#endif
int iScreenLogLevel;		// level of screen diagnostics
int iFileLogLevel;			// level of file diagnostics
char szLogFile[_MAX_PATH];	// write diagnostics to this file

int Rslt;
int Idx;
int ReLen;

int PMode;				// processing mode
int NumIncludeChroms;
char *pszIncludeChroms[cMaxIncludeChroms];
int NumExcludeChroms;
char *pszExcludeChroms[cMaxExcludeChroms];

char szInFile[_MAX_PATH];	// input alignments from this file
char szOutFile[_MAX_PATH];	// write filtered alignments to this file

char szSQLiteDatabase[_MAX_PATH];	// results summaries to this SQLite file
char szExperimentName[cMaxDatasetSpeciesChrom+1];			// experiment name
char szExperimentDescr[1000];		// describes experiment

// command line args
struct arg_lit  *help    = arg_lit0("hH","help",                "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"diagnostics log file");

struct arg_int *pmode = arg_int0("m","mode","<int>",		    "processing mode: 0 - default");
struct arg_str  *excludechroms = arg_strn("Z","chromexclude","<string>",0,cMaxExcludeChroms,"regular expressions defining chromosomes to always exclude");
struct arg_str  *includechroms = arg_strn("z","chromeinclude","<string>",0,cMaxIncludeChroms,"regular expressions defining chromosomes to explicitly include if not already excluded");
struct arg_file *infile = arg_file1("i","in","<file>",			"input alignments file to be filtered (SAM/BAM) file");
struct arg_file *outfile = arg_file1("o","output","<file>",		"write accepted alignments to this (SAM/BAM) file");
struct arg_file *summrslts = arg_file0("q","sumrslts","<file>",		"Output results summary to this SQLite3 database file");
struct arg_str *experimentname = arg_str0("w","experimentname","<str>",		"experiment name SQLite3 database file");
struct arg_str *experimentdescr = arg_str0("W","experimentdescr","<str>",	"experiment description SQLite3 database file");

struct arg_end *end = arg_end(200);

void *argtable[] = {help,version,FileLogLevel,LogFile,
					summrslts,experimentname,experimentdescr,
					pmode,excludechroms,includechroms,infile,outfile,
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
		printf("\n      To invoke this parameter file then precede its name with '@'");
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

	if(!gDiagnostics.Open(szLogFile,(etDiagLevel)iScreenLogLevel,(etDiagLevel)iFileLogLevel,true))
		{
		printf("\nError: Unable to start diagnostics subsystem.");
		if(szLogFile[0] != '\0')
			printf(" Most likely cause is that logfile '%s' can't be opened/created",szLogFile);
		exit(1);
		}
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Subprocess %s Version %s starting",gpszSubProcess->pszName,cpszProgVer);
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
	PMode = pmode->count ? pmode->ival[0] : 0;
	if(PMode < 0 || PMode > 0)
		{
		printf("\nError: Processing mode '-m%d' specified outside of range %d..%d",PMode,0,(int)0);
		exit(1);
		}

	NumIncludeChroms = includechroms->count;
	for(Idx=0;Idx < includechroms->count; Idx++)
		{
		ReLen = (int)strlen(includechroms->sval[Idx]);
		pszIncludeChroms[Idx] = new char [ReLen+1];
		strcpy(pszIncludeChroms[Idx],includechroms->sval[Idx]);
		TrimREQuotes(pszIncludeChroms[Idx]);
		}

	NumExcludeChroms = excludechroms->count;
	for(Idx=0;Idx < excludechroms->count; Idx++)
		{
		ReLen = (int)strlen(excludechroms->sval[Idx]);
		pszExcludeChroms[Idx] = new char [ReLen+1];
		strcpy(pszExcludeChroms[Idx],excludechroms->sval[Idx]);
		TrimREQuotes(pszExcludeChroms[Idx]);
		}



	strncpy(szInFile,infile->filename[0],_MAX_PATH);
	szInFile[_MAX_PATH-1] = '\0';

	strncpy(szOutFile,outfile->filename[0],_MAX_PATH);
	szOutFile[_MAX_PATH-1] = '\0';

// show user current resource limits
#ifndef _WIN32
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Resources: %s",CUtility::ReportResourceLimits());
#endif

	const char *pszDescr;
	switch(PMode) {
		case 0:
			pszDescr = "Filtering file alignment chromosomes";
			break;
		}

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Processing mode is : '%s'",pszDescr);

	for(Idx = 0; Idx < NumIncludeChroms; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"reg expressions defining chroms to include: '%s'",pszIncludeChroms[Idx]);
	for(Idx = 0; Idx < NumExcludeChroms; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"reg expressions defining chroms to exclude: '%s'",pszExcludeChroms[Idx]);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Input alignment file: '%s'",szInFile);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Output accepted alignments to file: '%s'",szOutFile);

	if(szExperimentName[0] != '\0')
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"This processing reference: %s",szExperimentName);

#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	gStopWatch.Start();
	Rslt = Process(NumIncludeChroms,pszIncludeChroms,NumExcludeChroms,pszExcludeChroms,szInFile,szOutFile);
	Rslt = Rslt >=0 ? 0 : 1;
	if(gExperimentID > 0)
		{
		if(gProcessingID)
			gSQLiteSummaries.EndProcessing(gExperimentID, gProcessingID,Rslt);
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
}

// TrimREQuotes
// Removes any leading and trailing quotes from pszTxt
int
TrimREQuotes(char *pszTxt)
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

int Process(int NumIncludeChroms,		// number of retained chromosomes regular expressions
	char **ppszIncludeChroms,			// array of include chromosome regular expressions
	int NumExcludeChroms,				// number of chromosome expressions to explicitly exclude
	char **ppszExcludeChroms,			// array of exclude chromosome regular expressions
	char *pszInFile,					// input file containing alignments to be filtered
	char *pszOutFile)					// write filtered alignments to this output file
{
CFilterSAMAlignments FilterSAMAlignments;
return(FilterSAMAlignments.FilterSAMbyChrom(NumIncludeChroms,ppszIncludeChroms,NumExcludeChroms,ppszExcludeChroms,pszInFile,pszOutFile));
}

CFilterSAMAlignments::CFilterSAMAlignments()
{
m_pInBAMfile = NULL;
m_pOutBAMfile = NULL;
}


CFilterSAMAlignments::~CFilterSAMAlignments()
{
}


void
CFilterSAMAlignments::Reset(void)
{
if(m_pInBAMfile != NULL)
	{
	m_pInBAMfile->Close();
	delete m_pInBAMfile;
	m_pInBAMfile = NULL;
	}
if(m_pOutBAMfile != NULL)
	{
	m_pOutBAMfile->Close();
	delete m_pOutBAMfile;
	m_pOutBAMfile = NULL;
	}

m_szFiltChrom[0] = 0;
m_bFiltChrom = false;

m_NumIncludeChroms = 0;			
m_ppszIncludeChroms = NULL;		
m_NumExcludeChroms = 0;			
m_ppszExcludeChroms = NULL;		
}


void
CFilterSAMAlignments::Init(void)
{
Reset();
}


char *
CFilterSAMAlignments::TrimWhitespace(char *pTxt)
{
char *pStart;
char Chr;
	// strip leading whitespace
while(Chr = *pTxt++)
	if(!isspace(Chr))
			break;
if(Chr == '\0')					// empty line?
	return(pTxt-1);
pStart = pTxt-1;
while(Chr = *pTxt)			// fast forward to line terminator
	pTxt++;
pTxt-=1;
while(Chr = *pTxt--)
	if(!isspace(Chr))
		break;
pTxt[2] = '\0';
return(pStart);
}

bool								// true if file to be generated compressed with gzopen/gzwrite/gzclose
CFilterSAMAlignments::FileReqWriteCompr(char *pszFile) // If last 3 chars of file name is ".gz" then this file is assumed to require compression
{
int Len;
if(pszFile == NULL || pszFile[0] == '\0')
	return(false);
if((Len = (int)strlen(pszFile)) < 4)
	return(false);
return(stricmp(".gz",&pszFile[Len-3]) == 0 ? true : false);
}



int			// eBSFSuccess or error code
CFilterSAMAlignments::SetChromFilters(int NumIncludeChroms,		 // number of chromosome regular expressions to include
						char **ppszIncludeChroms,	 // array of include chromosome regular expressions
						int NumExcludeChroms,		 // number of chromosome expressions to exclude
						char **ppszExcludeChroms)	 // array of exclude chromosome regular expressions
{
int Idx;

#ifndef _WIN32
int RegErr;				// regular expression parsing error
char szRegErr[128];			// to hold RegErr as textual representation ==> regerror();
#endif

m_NumIncludeChroms = NumIncludeChroms;		// number of chromosomes explicitly defined to be included
m_ppszIncludeChroms = NumIncludeChroms == 0 ? NULL : ppszIncludeChroms;	// ptr to array of reg expressions defining chroms to include - overides exclude
m_NumExcludeChroms = NumExcludeChroms;		// number of chromosomes explicitly defined to be excluded
m_ppszExcludeChroms =  NumExcludeChroms == 0 ? NULL : ppszExcludeChroms;	// ptr to array of reg expressions defining chroms to include

m_szFiltChrom[0] = 0;
m_bFiltChrom = false;

if(NumIncludeChroms == 0 && NumExcludeChroms == 0)
	return(eBSFSuccess);

#ifdef _WIN32
try {
	for(Idx=0;Idx < NumIncludeChroms;Idx++)
		{
		m_IncludeChromsRE[Idx] = new Regexp();
		m_IncludeChromsRE[Idx]->Parse(ppszIncludeChroms[Idx],false);	// note case insensitive
		}
	}
catch(...)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to process include regexpr chrom '%s'",ppszIncludeChroms[Idx]);
	Reset();
	return(eBSFerrMem);
	}
try {
	for(Idx=0;Idx < NumExcludeChroms;Idx++)
		{
		m_ExcludeChromsRE[Idx] = new Regexp();
		m_ExcludeChromsRE[Idx]->Parse(ppszExcludeChroms[Idx],false);	// note case insensitive
		}
	}
catch(...)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to process exclude regexpr chrom '%s'",ppszExcludeChroms[Idx]);
	Reset();
	return(eBSFerrMem);
	}

#else
for(Idx=0;Idx < NumIncludeChroms;Idx++)
	{

	RegErr=regcomp(&m_IncludeChromsRE[Idx],ppszIncludeChroms[Idx],REG_EXTENDED | REG_ICASE);	// note case insensitive
	if(RegErr)
		{
		regerror(RegErr,&m_IncludeChromsRE[Idx],szRegErr,sizeof(szRegErr));
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to process include chrom '%s' error: %s",ppszIncludeChroms[Idx],szRegErr);
		Reset();
		return(eBSFerrMem);
		}
	}
for(Idx=0;Idx < NumExcludeChroms;Idx++)
	{
	RegErr = regcomp(&m_ExcludeChromsRE[Idx],ppszExcludeChroms[Idx],REG_EXTENDED | REG_ICASE);	// note case insensitive
	if(RegErr)
		{
		regerror(RegErr,&m_ExcludeChromsRE[Idx],szRegErr,sizeof(szRegErr));
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to process exclude chrom '%s' error: %s",ppszExcludeChroms[Idx],szRegErr);
		Reset();
		return(eBSFerrMem);
		}
	}
#endif
return(eBSFSuccess);
}

// ExcludeThisChrom
// Returns true if pszChrom is to be excluded from processing
bool
CFilterSAMAlignments::ExcludeThisChrom(char *pszChrom)
{
#ifdef _WIN32
RegexpMatch mc;
#else
regmatch_t mc;
#endif

int Idx;
if(!m_NumExcludeChroms && !m_NumIncludeChroms)
	return(false);

if(m_szFiltChrom[0] != 0)
	{
	if(!stricmp(m_szFiltChrom,pszChrom))
		return(m_bFiltChrom);
	}
strcpy(m_szFiltChrom,pszChrom);
m_bFiltChrom = false;

// explicitly to be excluded?
for(Idx = 0; Idx < m_NumExcludeChroms; Idx++)
#ifdef _WIN32	
	if(m_ExcludeChromsRE[Idx]->Match(pszChrom,&mc))
#else
	if(!regexec(&m_ExcludeChromsRE[Idx],pszChrom,1,&mc,0))
#endif
		return(m_bFiltChrom = true);

// explicitly to be included?
for(Idx = 0; Idx < m_NumIncludeChroms; Idx++)
	{
#ifdef _WIN32
	if(m_IncludeChromsRE[Idx]->Match(pszChrom,&mc))
#else
	if(!regexec(&m_IncludeChromsRE[Idx],pszChrom,1,&mc,0))
#endif
		return(m_bFiltChrom = false);
	}


// if chromosomes were defined as to explicitly include then this chrom is to be filtered out
m_bFiltChrom = m_NumIncludeChroms > 0 ? true : false;
return(m_bFiltChrom);
}

int								// number of alignments which were retained and written to output file after filtering was applied
CFilterSAMAlignments::FilterSAMbyChrom(int NumIncludeChroms,		// number of retained chromosomes regular expressions
	char **ppszIncludeChroms,	// array of include chromosome regular expressions
	int NumExcludeChroms,		// number of chromosome expressions to explicitly exclude
	char **ppszExcludeChroms,	// array of exclude chromosome regular expressions
	char *pszInFile,			// input file containing alignments to be filtered
	char *pszOutFile)			// write filtered alignments to this output file
{
teBSFrsltCodes Rslt;

Reset();

if((Rslt=(teBSFrsltCodes)SetChromFilters(NumIncludeChroms,ppszIncludeChroms,NumExcludeChroms,ppszExcludeChroms)) != eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to initialise chrom include/exclude REs");
	Reset();
	return(Rslt);
	}

int NumParsedElLines;
int NumAcceptedEls;
int NumUnmappedEls;
int NumMappedChroms;
int LineLen;
char szLine[cMaxReadLen  * 3];				// buffer input lines
bool bFirstAlignment;
tsBAMalign ProvBAMalignment;
tsBAMalign AcceptedBAMalignment;

char szGenome[100];
char szContig[100];

int ContigLen;

int Tmp = 0;

char *pTxt;

	// if SAM output format then could be BGZF compressed BAM; use the extension to determine which...
	// if extension is '.bam' then BGZF compressed BAM, any other extension is for SAM
int Len;
int SAMFormat = 0;
Len = (int)strlen(pszOutFile);
if(Len > 5)
	{
	if(!stricmp(".bam",&pszOutFile[Len-4]))
		SAMFormat = 1;
	}

// open SAM for reading
if(pszInFile == NULL || *pszInFile == '\0')
	{
	Reset();
	return(eBSFerrParams);
	}

if((m_pInBAMfile = new CSAMfile) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"RemapSAMLocii: Unable to instantiate class CSAMfile");
	Reset();
	return(eBSFerrInternal);
	}

if((Rslt = (teBSFrsltCodes)m_pInBAMfile->Open(pszInFile)) != eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"RemapSAMLocii: Unable to open SAM/BAM format file %s",pszInFile);
	delete m_pInBAMfile;
	m_pInBAMfile = NULL;
	Reset();
	return((teBSFrsltCodes)Rslt);
	}

if((m_pOutBAMfile = new CSAMfile) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"RemapSAMLocii: Unable to instantiate class CSAMfile");
	delete m_pInBAMfile;
	m_pInBAMfile = NULL;
	Reset();
	return(eBSFerrInternal);
	}

eSAMFileType FileType;
bool bgzOutFile = FileReqWriteCompr(pszOutFile);

switch(SAMFormat) {
	case 0:			// output SAM
		if(bgzOutFile)
			FileType = eSFTSAMgz;
		else
			FileType = eSFTSAM;
		break;
	case 1:				// output as BAM compressed with bgzf
		FileType = eSFTBAM_BAI;
		break;
	}

if((Rslt = (teBSFrsltCodes)m_pOutBAMfile->Create(FileType,pszOutFile,6,(char *)cpszProgVer)) < eBSFSuccess) // defaulting to compression level 6 if compressed BAM 
	{
	delete m_pInBAMfile;
	m_pInBAMfile = NULL;
	delete m_pOutBAMfile;
	m_pInBAMfile = NULL;
	Reset();
	return(Rslt);
	}

NumParsedElLines = 0;
NumAcceptedEls = 0;
NumUnmappedEls = 0;
bFirstAlignment = true;
Rslt = eBSFSuccess;

NumMappedChroms = 0;
time_t Then = time(NULL);
time_t Now;
while(Rslt >= eBSFSuccess && (LineLen = m_pInBAMfile->GetNxtSAMline(szLine)) > 0)
	{
	NumParsedElLines += 1;
	if(!(NumParsedElLines % 100000) || NumParsedElLines == 1)
		{
		Now = time(NULL);
		if((Now - Then) >= 60)
			{
			if(bFirstAlignment)
				gDiagnostics.DiagOut(eDLInfo,gszProcName,"Parsed %d element lines, accepted %d Chroms for filtering",NumParsedElLines,NumMappedChroms);
			else
				gDiagnostics.DiagOut(eDLInfo,gszProcName,"Parsed %d element lines, unmapped %d, accepted %d",NumParsedElLines,NumUnmappedEls,NumAcceptedEls);
			Then += 60;
			}
		}

	szLine[sizeof(szLine)-1] = '\0';
	pTxt = TrimWhitespace(szLine);
	if(*pTxt=='\0')			// simply slough lines which are just whitespace
		continue;

	if(*pTxt== '@')			// if a reference sequence dictionary entry then check if an entry in the mapping loci BED
		{
		if(pTxt[1] != 'S' && pTxt[2] != 'Q')
			continue;
		if(3 != sscanf(&pTxt[3]," AS:%s SN:%s LN:%d",szGenome,szContig,&ContigLen))
			continue;

		if(ExcludeThisChrom(szContig))
			continue;

		m_pOutBAMfile->AddRefSeq(szGenome,szContig,ContigLen);
		NumMappedChroms += 1;
		continue;
		}	
	if(bFirstAlignment)
		{
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Parsed %d element lines, accepted %d Chroms for filtering",NumParsedElLines,NumMappedChroms);
		m_pOutBAMfile->StartAlignments();	
		bFirstAlignment = false;
		}
	if(!NumMappedChroms)
		break;

	// primary interest is in the reference chromname, startloci, length
	if((Rslt = (teBSFrsltCodes)m_pOutBAMfile->ParseSAM2BAMalign(pTxt,&ProvBAMalignment,NULL)) < eBSFSuccess)
		{
		if(Rslt == eBSFerrFeature)
			{
			NumUnmappedEls += 1;
			Rslt = eBSFSuccess;
			continue;
			}
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"ParseSAM2BAMalign() returned error: %d after parsed %d element lines, unmapped %d, accepted %d",Rslt, NumParsedElLines,NumUnmappedEls,NumAcceptedEls);
		break;
		}

		// check if read has been mapped, if not then slough ...
	if(ProvBAMalignment.refID == -1 || (ProvBAMalignment.flag_nc >> 16) & 0x04 || ProvBAMalignment.cigar[0] == '*')	// set if unmapped or Cigar is unknown
		{
		NumUnmappedEls += 1;
	    continue;
		}

	if(NumAcceptedEls > 0)
		{
		if((Rslt = (teBSFrsltCodes)m_pOutBAMfile->AddAlignment(&AcceptedBAMalignment,false)) < eBSFSuccess)
			{
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"AddAlignment() returned error: %d after parsed %d element lines, unmapped %d, accepted %d",Rslt, NumParsedElLines,NumUnmappedEls,NumAcceptedEls);
			break;
			}
		}
	AcceptedBAMalignment = ProvBAMalignment;
	NumAcceptedEls += 1;
	}

if(Rslt >= eBSFSuccess && NumAcceptedEls > 0)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Parsed %d element lines (final), unmapped %d, accepted %d",NumParsedElLines,NumUnmappedEls,NumAcceptedEls);
	if((Rslt = (teBSFrsltCodes)m_pOutBAMfile->AddAlignment(&AcceptedBAMalignment,true)) < eBSFSuccess)
		{
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"AddAlignment() returned error %d after parsed %d element lines (final), unmapped %d, accepted %d",Rslt,NumParsedElLines,NumUnmappedEls,NumAcceptedEls);
		}
	}
if(Rslt >= eBSFSuccess)
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Completed parsing %d element lines, unmapped %d, accepted %d",NumParsedElLines,NumUnmappedEls,NumAcceptedEls);
if(m_pInBAMfile != NULL)
	{
	m_pInBAMfile->Close();
	delete m_pInBAMfile;
	m_pInBAMfile = NULL;
	}
if(m_pOutBAMfile != NULL)
	{
	m_pOutBAMfile->Close();
	delete m_pOutBAMfile;
	m_pOutBAMfile = NULL;
	}

Reset();
return(Rslt >= 0 ? NumAcceptedEls : Rslt);
}

