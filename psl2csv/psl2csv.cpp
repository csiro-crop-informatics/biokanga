// psl2csv.cpp : Defines the entry point for the console application.
// Process BLAT psl output format into csv


#include "stdafx.h"
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#if _WIN32

#include "../libbiokanga/commhdrs.h"
#else
#include "../libbiokanga/commhdrs.h"
#endif

const unsigned int cProgVer = 209;		// increment with each release

const int cMaxIncludeChroms = 20;		// max number of include chromosomes regular expressions
const int cMaxExcludeChroms = 20;		// max number of exclude chromosomes regular expressions

const int cMaxInBuffSize  = 10000000;	// read in chunks of this size from source psl file
const int cMaxOutBuffSize = 10000000;	// and write in chunks of this size to output csv file

const int cMaxLenPSLline = 32000;		// max length PSL line expected - just a guess!
const int cMaxNumPSLblocks = 1000;		// max number of blocks in any PSL alignment - again just a guess!

typedef struct TAG_sPSLline {
	int Matches;			// number of matches which aren't repeats
	int misMatches;			// number of bases which do not match
	int repMatches;			// number of bases which match but are also repeats
	int nCount;				// number of N bases
	int qNumInsert;			// number of InDel seqs in query
	int qBaseInsert;		// number of bases total in all InDels in query
	int tNumInsert;			// number of InDel seqs in target
	int tBaseInsert;		// number of bases total in all InDels in target
	char strand[3];			// '+' or '-' for query strand, for translated alignments second '+'/'-' is for genomic strand
	char qName[128];		// query sequence name
	int qNameLen;			// query name length 
	int  qSize;				// query sequence size
	int  qStart;			// alignment start psn in query
	int  qEnd;				// alignment end psn in query
	char tName[128];		// target sequence name
	int  tSize;				// target sequence size
	int  tStart;			// alignment start psn in target
	int  tEnd;				// alignment end psn in target
	int  blockCount;		// number of blocks in the alignment
	int  blockSizes[cMaxNumPSLblocks];		// sizes of each block
	int  qStarts[cMaxNumPSLblocks];			// starting psn of each block in query
	int  tStarts[cMaxNumPSLblocks];			// starting psn of each block in target
	} tsPSLline;


const int cMaxExcludeHistory = 100;
typedef struct TAG_sExcludeEl {
	struct TAG_sExcludeEl *pNext;
	struct TAG_sExcludeEl *pPrev;
	char szChrom[cMaxDatasetSpeciesChrom];	// chromosome is not to be processed
	bool bExclude;							// true if to be excluded, false if not
	} tsExcludeEl;

typedef struct TAG_sProcParams 
	{
	int hPSLinFile;				// opened file handle for psl files
	int hOutFile;				// opened file handle for csv file being written
	unsigned char *pInBuffer;	// mem allocd to buffer chars being read from psl
	int NumInBuffer;			// num of chars currently in pInBuffer
	int InBuffIdx;				// index of next char to read from pInBuffer[]
	int PushedBack;				// last pushed back char (only 1 level of pushback supported!)

	char szOutFile[_MAX_PATH];	// output csv file
	char szPSLinFile[_MAX_PATH];   // current PSL file

	char *pszLineBuff;			// allocd PSL line buffer
	int CurLineLen;				// current line length in pszLineBuff
	tsPSLline PSLline;			// parsed PSL alignment

	int NumBlatHits;			// number of blat hits parsed
	int NumFiltHitChroms;			// number of blat hits filtered out because they hit on chroms to exclude
	int NumBlatHitsAccepted;	// number of blat hits accepted (NumBlatHits - NumFiltHitChroms) 

	int NumIncludeChroms;			// number of chromosomes explicitly defined to be included
	char **ppszIncludeChroms;		// ptr to array of reg expressions defining chroms to include - overides exclude
	int NumExcludeChroms;			// number of chromosomes explicitly defined to be excluded
	char **ppszExcludeChroms;		// ptr to array of reg expressions defining chroms to include
#ifdef _WIN32
	Regexp *IncludeChromsRE[cMaxIncludeChroms];	// compiled regular expressions
	Regexp *ExcludeChromsRE[cMaxExcludeChroms];
#else
	regex_t IncludeChromsRE[cMaxIncludeChroms];	// compiled regular expressions
	regex_t ExcludeChromsRE[cMaxExcludeChroms];
#endif

	int gNumExcludeEls;		// current number of elements in gExcludeChroms
	tsExcludeEl *gpMRA;		// pts to most recently accessed or added
	tsExcludeEl *gpLRA;		// pts to least recently accessed
	tsExcludeEl gExcludeChroms[cMaxExcludeHistory];
	} tsProcParams;


int 
Process(int NumIncludeChroms,		// number of chromosome regular expressions to include
		char **ppszIncludeChroms,	// array of include chromosome regular expressions
		int NumExcludeChroms,		// number of chromosome expressions to exclude
		char **ppszExcludeChroms,	// array of exclude chromosome regular expressions
		char *pszPSLinFile,
		char *pszOutFile);
int ProcessPSL(tsProcParams *pParams);
int ProcessPSLline(tsProcParams *pParams);
int	GetNext(tsProcParams *pParams);
int OutputCSV(tsProcParams *pParams);
int TrimQuotes(char *pszTxt);
bool ExcludeThisChrom(char *pszChrom,tsProcParams *pProcParams);

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
	return _T("psl2csv");
}
// end of str library required code
#endif


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
int Idx;
int LenFileList;

char szPSLinFile[_MAX_PATH];	// process from these psl files
char szOutFile[_MAX_PATH];		// output into this csv file
int NumIncludeChroms;
char *pszIncludeChroms[cMaxIncludeChroms];
int NumExcludeChroms;
char *pszExcludeChroms[cMaxExcludeChroms];

// command line args
struct arg_lit  *help    = arg_lit0("hH","help",                "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_int *ScreenLogLevel=arg_int0("S", "ScreenLogLevel",	"<int>","Level of diagnostics written to logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"diagnostics log file");

struct arg_file *InFile = arg_file1("i","infile","<file>",		"process from these psl files");
struct arg_file *OutFile = arg_file1("o","outfile","<file>",	"process into this csv file");

struct arg_str  *ExcludeChroms = arg_strn("Z","chromexclude","<string>",0,cMaxExcludeChroms,"high priority: regular expressions defining target hit chromosomes to exclude");
struct arg_str  *IncludeChroms = arg_strn("z","chromeinclude","<string>",0,cMaxIncludeChroms,"low priority: regular expressions defining target hit chromosomes to include");

struct arg_end *end = arg_end(20);

void *argtable[] = {help,version,FileLogLevel,ScreenLogLevel,LogFile,
					InFile,OutFile,ExcludeChroms,IncludeChroms,
					end};

char **pAllArgs;
int argerrors;
argerrors = CUtility::arg_parsefromfile(argc,(char **)argv,&pAllArgs);
if(argerrors >= 0)
	argerrors = arg_parse(argerrors,pAllArgs,argtable);

    /* special case: '--help' takes precedence over error reporting */
if (help->count > 0)
        {
		printf("Usage: %s ", gszProcName);
        arg_print_syntax(stdout,argtable,"\n");
        arg_print_glossary(stdout,argtable,"  %-25s %s\n");
		exit(1);
        }

    /* special case: '--version' takes precedence error reporting */
if (version->count > 0)
        {
			printf("\n%s Version: %d.%2.2d\n",gszProcName,cProgVer/100,cProgVer%100);
		exit(1);
        }
if (!argerrors)
	{
	iScreenLogLevel = ScreenLogLevel->count ? ScreenLogLevel->ival[0] : eDLInfo;
	if(iScreenLogLevel < eDLNone || iScreenLogLevel > eDLDebug)
		{
		printf("\nError: ScreenLogLevel '-S%d' specified outside of range %d..%d",iScreenLogLevel,eDLNone,eDLDebug);
		exit(1);
		}
	if(FileLogLevel->count && !LogFile->count)
		{
		printf("\nError: FileLogLevel '-f%d' specified but no logfile '-F<logfile>'",FileLogLevel->ival[0]);
		exit(1);
		}

	iFileLogLevel = FileLogLevel->count ? FileLogLevel->ival[0] : eDLInfo;
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

	NumIncludeChroms = IncludeChroms->count;
	for(Idx=0;Idx < IncludeChroms->count; Idx++)
		{
		LenFileList = (int)strlen(IncludeChroms->sval[Idx]);
		pszIncludeChroms[Idx] = new char [LenFileList+1];
		strcpy(pszIncludeChroms[Idx],IncludeChroms->sval[Idx]);
		TrimQuotes(pszIncludeChroms[Idx]);
		}

	NumExcludeChroms = ExcludeChroms->count;
	for(Idx=0;Idx < ExcludeChroms->count; Idx++)
		{
		LenFileList = (int)strlen(ExcludeChroms->sval[Idx]);
		pszExcludeChroms[Idx] = new char [LenFileList+1];
		strcpy(pszExcludeChroms[Idx],ExcludeChroms->sval[Idx]);
		TrimQuotes(pszExcludeChroms[Idx]);
		}

	strcpy(szPSLinFile,InFile->filename[0]);
	strcpy(szOutFile,OutFile->filename[0]);

		// now that command parameters have been parsed then initialise diagnostics log system
	if(!gDiagnostics.Open(szLogFile,(etDiagLevel)iScreenLogLevel,(etDiagLevel)iFileLogLevel,true))
		{
		printf("\nError: Unable to start diagnostics subsystem.");
		if(szLogFile[0] != '\0')
			printf(" Most likely cause is that logfile '%s' can't be opened/created",szLogFile);
		exit(1);
		}

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Version: %d.%2.2d Processing parameters:",cProgVer/100,cProgVer%100);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"input psl files: '%s'",szPSLinFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"output csv file:   '%s'",szOutFile);
	for(Idx = 0; Idx < NumIncludeChroms; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"reg expressions defining chroms to include: '%s'",pszIncludeChroms[Idx]);
	for(Idx = 0; Idx < NumExcludeChroms; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"reg expressions defining chroms to exclude: '%s'",pszExcludeChroms[Idx]); 

	gStopWatch.Start();
#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	Rslt = Process(NumIncludeChroms,pszIncludeChroms,NumExcludeChroms,pszExcludeChroms,szPSLinFile,szOutFile);
	gStopWatch.Stop();
	Rslt = Rslt < 0 ? 1 : 0;
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Exit Code: %d Total processing time: %s",Rslt,gStopWatch.Read());
	exit(Rslt);
	}
else
	{
	arg_print_errors(stdout,end,gszProcName);
	arg_print_syntax(stdout,argtable,"\nend of help\n");
	exit(1);
	}
return 0;
}


int
ProcessThisFile(char *pszInPSL,void *pHandlerParams)
{
int Rslt;
tsProcParams *pParams = (tsProcParams *)pHandlerParams;
strcpy(pParams->szPSLinFile,pszInPSL);

#ifdef _WIN32
if((pParams->hPSLinFile = open(pszInPSL,_O_RDWR | _O_BINARY | _O_SEQUENTIAL))==-1)
#else
if((pParams->hPSLinFile = open(pszInPSL,O_RDWR))==-1)
#endif
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: Unable to open input file for processing - '%s' - %s", 
				pszInPSL,strerror(errno));
	return(eBSFerrOpnFile);
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"ProcessThisFile: Processing input file '%s'",pszInPSL); 

pParams->NumInBuffer = 0;
pParams->InBuffIdx = 0;

Rslt = ProcessPSL(pParams);

if(Rslt >= 0)
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"File '%s' Accepted: %d, Parsed: %d, Filtered because hit to chrom: %d",
					 pszInPSL,pParams->NumBlatHitsAccepted,pParams->NumBlatHits,pParams->NumFiltHitChroms); 
pParams->NumBlatHitsAccepted = 0;
pParams->NumBlatHits = 0;
pParams->NumFiltHitChroms = 0;

close(pParams->hPSLinFile);
pParams->hPSLinFile = -1;

if(Rslt < 0 && pParams->hOutFile != -1)
	{
	close(pParams->hOutFile);
	pParams->hOutFile = -1;
	}

return(Rslt);
}



int
ProcessPSL(tsProcParams *pParams)
{
int Rslt;
bool bInWhiteSpace = false;		

pParams->CurLineLen = 0;
Rslt = 0;
while(Rslt >= 0 && (Rslt = GetNext(pParams)) > 0) {
	switch((char)Rslt) {
		case '\r':			// silently slough CR - must have been generated on windows/msdos machine
			continue;

		case '\n':				// accept linefeeds - both Linux and windows are happy with lines terminated by NL
			bInWhiteSpace = false;
			if(pParams->CurLineLen)
				Rslt = ProcessPSLline(pParams);
			pParams->CurLineLen = 0;
			continue;
		
		default:
			if(isspace(Rslt))
				{
				if(bInWhiteSpace)	// slough multiple whitespace
					continue;
				bInWhiteSpace = true;
				}
			else
				bInWhiteSpace = false;
			pParams->pszLineBuff[pParams->CurLineLen++] = (char)Rslt;
			continue;
			}
	}
if(!Rslt && pParams->CurLineLen)
	Rslt = ProcessPSLline(pParams);
return(Rslt);
}

//
// ProcessPSLline
// Parses PSL line
// 
int
ProcessPSLline(tsProcParams *pParams)
{
int Cnt;
int Psn;
int BlockIdx;
char *pChr;
tsPSLline *pPSL;
pParams->pszLineBuff[pParams->CurLineLen] = '\0';	
pPSL = &pParams->PSLline;
pChr = pParams->pszLineBuff;
while(isspace(*pChr))
	pChr++;
if(!isdigit(*pChr))	// assume header line if not a digit	
	return(0);

Cnt = sscanf(pChr,"%d %d %d %d %d %d %d %d %n",
				&pPSL->Matches,&pPSL->misMatches,&pPSL->repMatches,&pPSL->nCount,&pPSL->qNumInsert,
				&pPSL->qBaseInsert,&pPSL->tNumInsert,&pPSL->tBaseInsert,&Psn);
if(Cnt != 8)
	return(-1);

pChr += Psn;
Cnt = sscanf(pChr,"%s %s %d %d %d %s %d %d %d %d %n",
				pPSL->strand,pPSL->qName,&pPSL->qSize,&pPSL->qStart,&pPSL->qEnd,
				pPSL->tName,&pPSL->tSize,&pPSL->tStart,&pPSL->tEnd,&pPSL->blockCount,&Psn);
if(Cnt != 10)
	return(-1);


pPSL->qNameLen = (int)strlen(pPSL->qName);

pChr += Psn;
for(BlockIdx = 0; BlockIdx < pPSL->blockCount; BlockIdx++)
	{
	Cnt = sscanf(pChr,"%d , %n",&pPSL->blockSizes[BlockIdx],&Psn);
	if(Cnt != 1)
		return(-1);
	pChr += Psn;
	}
for(BlockIdx = 0; BlockIdx < pPSL->blockCount; BlockIdx++)
	{
	Cnt = sscanf(pChr,"%d , %n",&pPSL->qStarts[BlockIdx],&Psn);
	if(Cnt != 1)
		return(-1);
	pChr += Psn;
	}
for(BlockIdx = 0; BlockIdx < pPSL->blockCount; BlockIdx++)
	{
	if(BlockIdx == pPSL->blockCount-1)
		Cnt = sscanf(pChr,"%d",&pPSL->tStarts[BlockIdx]);
	else
		Cnt = sscanf(pChr,"%d , %n",&pPSL->tStarts[BlockIdx],&Psn);
	if(Cnt != 1)
		return(-1);
	pChr += Psn;
	}

pParams->NumBlatHits += 1;

// this where to filter on chromsome (tName)
if(ExcludeThisChrom(pPSL->tName,pParams))
	{
	pParams->NumFiltHitChroms += 1;
	return(0);
	}
//
pParams->NumBlatHitsAccepted += 1;
return(OutputCSV(pParams));
}

int		// 0: EOF -1: error >0 chr
GetNext(tsProcParams *pParams)
{
int Chr;
if((Chr = pParams->PushedBack) > 0)
	{
	pParams->PushedBack = 0;
	return(Chr);
	}
if(pParams->InBuffIdx == -1 || pParams->InBuffIdx >= pParams->NumInBuffer)
	{
	pParams->NumInBuffer = read(pParams->hPSLinFile,pParams->pInBuffer,cMaxInBuffSize);
	if(pParams->NumInBuffer <= 0)
		{
		pParams->InBuffIdx = -1;
		return(pParams->NumInBuffer);
		}
	pParams->InBuffIdx = 0;
	}
return(pParams->pInBuffer[pParams->InBuffIdx++]);
}

void
Cleanup(tsProcParams *pProcParams)
{
if(pProcParams->pInBuffer != NULL)
	{
	delete pProcParams->pInBuffer;
	pProcParams->pInBuffer = NULL;
	}
if(pProcParams->pszLineBuff != NULL)
	{
	delete pProcParams->pszLineBuff;
	pProcParams->pszLineBuff = NULL;
	}
	
if(pProcParams->hOutFile != -1)			// close output file if still opened
	{
	close(pProcParams->hOutFile);
	pProcParams->hOutFile = -1;
	}
if(pProcParams->hPSLinFile != -1)
	{
	close(pProcParams->hPSLinFile);
	pProcParams->hPSLinFile = -1;
	}
}


int 
Process(int NumIncludeChroms,		// number of chromosome regular expressions to include
		char **ppszIncludeChroms,	// array of include chromosome regular expressions
		int NumExcludeChroms,		// number of chromosome expressions to exclude
		char **ppszExcludeChroms,	// array of exclude chromosome regular expressions
		char *pszPSLinFile,
		char *pszOutFile)
{
int Rslt;
int Idx;
tsProcParams ProcParams;		// initialised to hold processing parameters

#ifndef _WIN32
int RegErr;				// regular expression parsing error
char szRegErr[128];			// to hold RegErr as textual representation ==> regerror();
#endif

memset(&ProcParams,0,sizeof(ProcParams));

ProcParams.hPSLinFile = -1;
ProcParams.hOutFile = -1;
ProcParams.InBuffIdx = 0;
ProcParams.NumInBuffer = 0;
ProcParams.PushedBack = 0;
strcpy(ProcParams.szOutFile,pszOutFile);

if((ProcParams.pInBuffer = new unsigned char [cMaxInBuffSize])==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: Unable to allocate memory (%d bytes) for input buffering", 
				cMaxInBuffSize);
	return(eBSFerrMem);
	}


ProcParams.InBuffIdx = 0;
ProcParams.NumInBuffer = 0;
ProcParams.PushedBack = 0;

if((ProcParams.pszLineBuff = new char [cMaxLenPSLline])==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: Unable to allocate memory (%d bytes) for PSL line buffer", 
				cMaxLenPSLline);
	
	Cleanup(&ProcParams);
	return(eBSFerrMem);
	}

ProcParams.NumIncludeChroms = NumIncludeChroms;		// number of chromosomes explicitly defined to be included
ProcParams.ppszIncludeChroms = ppszIncludeChroms;	// ptr to array of reg expressions defining chroms to include - overides exclude
ProcParams.NumExcludeChroms = NumExcludeChroms;		// number of chromosomes explicitly defined to be excluded
ProcParams.ppszExcludeChroms = ppszExcludeChroms;	// ptr to array of reg expressions defining chroms to include

#ifdef _WIN32
try {
	for(Idx=0;Idx < NumIncludeChroms;Idx++)
		{
		ProcParams.IncludeChromsRE[Idx] = new Regexp();
		ProcParams.IncludeChromsRE[Idx]->Parse(ppszIncludeChroms[Idx],false);	// note case insensitive
		}
	}
catch(...)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to process include regexpr chrom '%s'",ppszIncludeChroms[Idx]);
	Cleanup(&ProcParams);
	return(eBSFerrMem);
	}
try {
	for(Idx=0;Idx < NumExcludeChroms;Idx++)
		{
		ProcParams.ExcludeChromsRE[Idx] = new Regexp();
		ProcParams.ExcludeChromsRE[Idx]->Parse(ppszExcludeChroms[Idx],false);	// note case insensitive
		}
	}
catch(...)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to process exclude regexpr chrom '%s'",ppszExcludeChroms[Idx]);
	Cleanup(&ProcParams);
	return(eBSFerrMem);
	}

#else
for(Idx=0;Idx < NumIncludeChroms;Idx++)
	{

	RegErr=regcomp(&ProcParams.IncludeChromsRE[Idx],ppszIncludeChroms[Idx],REG_EXTENDED | REG_ICASE);	// note case insensitive
	if(RegErr)
		{
		regerror(RegErr,&ProcParams.IncludeChromsRE[Idx],szRegErr,sizeof(szRegErr));
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to process include chrom '%s' error: %s",ppszIncludeChroms[Idx],szRegErr);
		Cleanup(&ProcParams);
		return(eBSFerrMem);
		}
	}
for(Idx=0;Idx < NumExcludeChroms;Idx++)
	{
	RegErr = regcomp(&ProcParams.ExcludeChromsRE[Idx],ppszExcludeChroms[Idx],REG_EXTENDED | REG_ICASE);	// note case insensitive
	if(RegErr)
		{
		regerror(RegErr,&ProcParams.ExcludeChromsRE[Idx],szRegErr,sizeof(szRegErr));
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to process exclude chrom '%s' error: %s",ppszExcludeChroms[Idx],szRegErr);
		Cleanup(&ProcParams);
		return(eBSFerrMem);
		}
	}
#endif


#ifdef _WIN32
if((ProcParams.hOutFile = open(pszOutFile,( _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC),(_S_IREAD | _S_IWRITE)))==-1)
#else
if((ProcParams.hOutFile = open(pszOutFile,O_RDWR | O_CREAT |O_TRUNC, S_IREAD | S_IWRITE))==-1)
#endif
{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: Unable to open output file for processing - '%s' - %s", 
			pszOutFile,strerror(errno));
	Cleanup(&ProcParams);
	return(eBSFerrOpnFile);
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"ProcessThisFile: Generating output file '%s'",pszOutFile); 

CSimpleGlob glob(SG_GLOB_FULLSORT);
if (glob.Add(pszPSLinFile) >= SG_SUCCESS)
	{
	Rslt = eBSFSuccess;
	for (int n = 0; Rslt >= eBSFSuccess &&  n < glob.FileCount(); ++n)
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Will process in this order: %d '%s", n+1,glob.File(n));

    for (int n = 0; Rslt >= eBSFSuccess &&  n < glob.FileCount(); ++n)
		Rslt = ProcessThisFile(glob.File(n),&ProcParams);
	}
else
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to glob '%s",pszPSLinFile);
	Rslt = eBSFerrOpnFile;	// treat as though unable to open file
    }

Cleanup(&ProcParams);
return(Rslt);
}

int
OutputCSV(tsProcParams *pParams)
{
static bool bFirst=true;
char szLineBuff[4096];
int BuffIdx;
int BlockIdx;
tsPSLline *pPSL = &pParams->PSLline;

if(bFirst)
	{
	bFirst=false;
	BuffIdx=sprintf(szLineBuff,"\"QName\",\"QLen\",\"QAlignLen\",\"QStart\",\"QEnd\",\"TName\",\"TStrand\",\"TLen\",\"TAlignLen\",\"TStart\",\"TEnd\",\"Matches\",\"Mismatches\",\"MatchNonRepeats\",\"MatchRepeats\",\"NBases\",\"QNumInDels\",\"QInDelsBases\",\"TNumInDels\",\"TInDelsBases\",\"BlockCnt\",\"BlockLens\",\"QBlockStarts\",\"TBlockStarts\"\n");
	CUtility::SafeWrite(pParams->hOutFile,szLineBuff,BuffIdx);
	}

BuffIdx = sprintf(szLineBuff,"\"%s\",%d,%d,%d,%d,\"%s\",\"%s\",%d,%d,%d,%d,",
							pPSL->qName,pPSL->qSize,(pPSL->qEnd - pPSL->qStart),pPSL->qStart,pPSL->qEnd,
							pPSL->tName,pPSL->strand,pPSL->tSize,(pPSL->tEnd - pPSL->tStart),pPSL->tStart,pPSL->tEnd);

BuffIdx += sprintf(&szLineBuff[BuffIdx],"%d,%d,%d,%d,%d,%d,%d,%d,%d,%d",
						pPSL->Matches+pPSL->repMatches,
						pPSL->misMatches,pPSL->Matches,pPSL->repMatches,pPSL->nCount,pPSL->qNumInsert,
						pPSL->qBaseInsert,pPSL->tNumInsert,pPSL->tBaseInsert,pPSL->blockCount);

CUtility::SafeWrite(pParams->hOutFile,szLineBuff,BuffIdx);
BuffIdx=sprintf(szLineBuff,",\"");
for(BlockIdx = 0; BlockIdx < pPSL->blockCount; BlockIdx++)
	BuffIdx += sprintf(&szLineBuff[BuffIdx],"%d,",pPSL->blockSizes[BlockIdx]);
BuffIdx -= 1;		// removes last ','
BuffIdx += sprintf(&szLineBuff[BuffIdx],"\"");
CUtility::SafeWrite(pParams->hOutFile,szLineBuff,BuffIdx);

BuffIdx=sprintf(szLineBuff,",\"");
for(BlockIdx = 0; BlockIdx < pPSL->blockCount; BlockIdx++)
	BuffIdx += sprintf(&szLineBuff[BuffIdx],"%d,",pPSL->qStarts[BlockIdx]);
BuffIdx -= 1;		// removes last ','
BuffIdx += sprintf(&szLineBuff[BuffIdx],"\"");
CUtility::SafeWrite(pParams->hOutFile,szLineBuff,BuffIdx);

BuffIdx=sprintf(szLineBuff,",\"");
for(BlockIdx = 0; BlockIdx < pPSL->blockCount; BlockIdx++)
	BuffIdx += sprintf(&szLineBuff[BuffIdx],"%d,",pPSL->tStarts[BlockIdx]);
BuffIdx -= 1;		// removes last ','
BuffIdx += sprintf(&szLineBuff[BuffIdx],"\"\n");
CUtility::SafeWrite(pParams->hOutFile,szLineBuff,BuffIdx);
return(0);
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

// AddExcludeHistory
// Adds a tsExcludeEl to the cached history
// The newly added element will be the MRA
bool
AddExcludeHistory(char *pszChrom,bool bExclude,tsProcParams *pParams)
{
tsExcludeEl *pEl;
if(pParams->gNumExcludeEls < cMaxExcludeHistory)
	pEl = &pParams->gExcludeChroms[pParams->gNumExcludeEls++];
else
	{
	pEl = pParams->gpLRA;		// reuse the least recently accessed element
	pParams->gpLRA = pEl->pPrev;	// 2nd LRA now becomes the LRA
	pParams->gpLRA->pNext = NULL;
	}
if(pParams->gpMRA != NULL)
	pParams->gpMRA->pPrev = pEl;
pEl->pNext = pParams->gpMRA;
pEl->pPrev = NULL;
pParams->gpMRA = pEl;
if(pParams->gpLRA == NULL)
	pParams->gpLRA = pEl;
pEl->bExclude = bExclude;
strcpy(pEl->szChrom,pszChrom);
return(bExclude);
}

// LocateExclude
// Locates - starting from the MRA - a tsExcludeEl which matches on szChrom
// If matches then this tsExcludeEl is made the MRA
// Returns ptr to matching tsExcludeEl if located or NULL
tsExcludeEl *
LocateExclude(char *pszChrom,tsProcParams *pParams)
{
tsExcludeEl *pEl;
pEl = pParams->gpMRA;
while(pEl != NULL)
	{
	if(!stricmp(pEl->szChrom,pszChrom))
		{
		if(pParams->gNumExcludeEls ==1 || pParams->gpMRA == pEl)	// if only, or already the MRA then no need for any relinking 
			return(pEl);

		if(pParams->gpLRA == pEl)						// if was the LRA then the 2nd LRA becomes the LRA
			{
			pParams->gpLRA = pEl->pPrev;
			pParams->gpLRA->pNext = NULL;
			}
		else									// not the LRA, and not the MRA
			{
			pEl->pPrev->pNext = pEl->pNext;
			pEl->pNext->pPrev = pEl->pPrev;
			}
		pParams->gpMRA->pPrev = pEl;
		pEl->pNext = pParams->gpMRA;
		pEl->pPrev = NULL;
		pParams->gpMRA = pEl;
		return(pEl);
		}
	pEl = pEl->pNext;
	}
return(NULL);
}


// ExcludeThisChrom
// Returns true if szChromosome is to be excluded from processing
// ExcludeThisChrom
bool
ExcludeThisChrom(char *pszChrom,tsProcParams *pProcParams)
{
#ifdef _WIN32
RegexpMatch mc;
#else
regmatch_t mc;
#endif

int Idx;
if(!pProcParams->NumExcludeChroms && !pProcParams->NumIncludeChroms)
	return(false);

// check if this species and chromosome are already known to be included/excluded
tsExcludeEl *pEl;
if((pEl = LocateExclude(pszChrom,pProcParams))!=NULL)
	return(pEl->bExclude);
// haven't seen this chromosome before - or else they have been discarded from history...
// to be excluded?
for(Idx = 0; Idx < pProcParams->NumExcludeChroms; Idx++)
#ifdef _WIN32	
	if(pProcParams->ExcludeChromsRE[Idx]->Match(pszChrom,&mc))
#else
	if(!regexec(&pProcParams->ExcludeChromsRE[Idx],pszChrom,1,&mc,0))
#endif
		return(AddExcludeHistory(pszChrom,true,pProcParams));

// to be included?
for(Idx = 0; Idx < pProcParams->NumIncludeChroms; Idx++)
	{
#ifdef _WIN32
	if(pProcParams->IncludeChromsRE[Idx]->Match(pszChrom,&mc))
#else
	if(!regexec(&pProcParams->IncludeChromsRE[Idx],pszChrom,1,&mc,0))
#endif
		return(AddExcludeHistory(pszChrom,false,pProcParams));
	}

return(AddExcludeHistory(pszChrom,pProcParams->NumIncludeChroms > 0 ? true : false,pProcParams));
}



