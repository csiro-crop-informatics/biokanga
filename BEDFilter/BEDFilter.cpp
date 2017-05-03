// BEDFilter.cpp : Defines the entry point for the console application.
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

const char *cpszProgVer = "0.0.1";		// increment with each release

const int cMinFeatLen  = 1;				// features must be of at least this length
const int cDfltFeatLen = 1;				// default feature length
const int cMaxFeatLen = 1000000000;		// max feature length

const int cMaxExcludeChroms = 20;		// allow upto this many regexpr for specifying chroms to exclude
const int cMaxIncludeChroms = 20;		// allow upto this many regexpr for specifying chroms to include

const size_t cAllocBuffSize = 1000000; // buffer size to allocate for holding formatted output ready to file write

// processing modes
typedef enum eProcMode {
	ePMDefault,						// default processing
	ePMplaceholder					// used to set the enumeration range
} etProcMode;

int
Process(etProcMode ProcMode,	// processing mode
		char Strand,			// process for this strand only
		int MinLen,				// retained features must be at least this length
		int MaxLen,				// retained features must be no longer than this length
		int NumIncludeChroms,	// number of chromosome regular expressions to include
		char **ppszIncludeChroms,	// array of include chromosome regular expressions
		int NumExcludeChroms,		// number of chromosome expressions to exclude
		char **ppszExcludeChroms,	// array of exclude chromosome regular expressions
		char *pszInFile,		// BED file containing elements to be filtered
		char *pszOutFile);		// write combined elements into this BED file

int TrimQuotes(char *pszTxt);
char *ProcMode2Txt(etProcMode ProcMode);

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
	return _T("BEDFilter");
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
int MaxLen;
int LenFileList;
int NumIncludeChroms;
char *pszIncludeChroms[cMaxIncludeChroms];
int NumExcludeChroms;
char *pszExcludeChroms[cMaxExcludeChroms];
char szInFile[_MAX_PATH];  // input BED file

char szOutFile[_MAX_PATH];	// write merged to this file

// command line args
struct arg_lit  *help    = arg_lit0("h","help",                 "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to screen and logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"diagnostics log file");
struct arg_int  *ProcMode = arg_int0("m","mode","<int>",		"processing mode: 0 - default");
struct arg_int  *strand = arg_int0("s","strand","<int>",		"filter for this strand: 0 - any, 1 - Watson '+', 2 - Crick '-' (default is any)");
struct arg_int  *minlen = arg_int0("l","minlen","<int>",		"retained features must be of at least this length (default is 1)");
struct arg_int  *maxlen = arg_int0("L","maxlen","<int>",		"retained features must be no longer than this length (default is 20)");

struct arg_str  *excludechroms = arg_strn("Z","chromexclude","<string>",0,cMaxExcludeChroms,"regular expressions defining chromosomes to always exclude");
struct arg_str  *includechroms = arg_strn("z","chromeinclude","<string>",0,cMaxIncludeChroms,"regular expressions defining chromosomes to explicitly include if not already excluded");
struct arg_file *infile = arg_file1("i","infile","<file>",		"filter features contained in this BED file");
struct arg_file *outfile = arg_file1("o","outfile","<file>",	"output retained features to this BED file");
struct arg_end *end = arg_end(20);

void *argtable[] = {help,version,FileLogLevel,LogFile,
					ProcMode,strand,minlen,maxlen,excludechroms,includechroms,infile,outfile,
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
		printf("\n      To invoke this parameter file then precede its name with '@'");
		printf("\n      e.g. %s @myparams.txt\n",gszProcName);
		printf("\nPlease report any issues regarding usage of %s at https://github.com/csiro-crop-informatics/biokanga/issues\n\n",gszProcName);
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

	MinLen = minlen->count ? minlen->ival[0] : cDfltFeatLen;
	if(MinLen < cMinFeatLen || MinLen > cMaxFeatLen)
		{
		printf("Error: minimum length '-l%d' is not in range %d..%d",MinLen,cMinFeatLen,cMaxFeatLen);
		exit(1);
		}

	MaxLen = maxlen->count ? maxlen->ival[0] : cMaxFeatLen;
	if(MaxLen < MinLen || MaxLen > cMaxFeatLen)
		{
		printf("Error: maximum length '-L%d' is not in range %d..%d",MaxLen,MinLen,cMaxFeatLen);
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


	strncpy(szInFile,infile->filename[0],_MAX_PATH);
	szInFile[_MAX_PATH-1] = '\0';
	strncpy(szOutFile,outfile->filename[0],_MAX_PATH);
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
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"retained features must be at least this length: %d",MinLen);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"retained features must be no longer than this length: %d",MaxLen);
	for(Idx = 0; Idx < NumIncludeChroms; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"reg expressions defining chroms to include: '%s'",pszIncludeChroms[Idx]);
	for(Idx = 0; Idx < NumExcludeChroms; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"reg expressions defining chroms to exclude: '%s'",pszExcludeChroms[Idx]);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Input features to filter from BED file: '%s'",szInFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Output merged features into BED file: '%s'",szOutFile);

	gStopWatch.Start();
#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	Rslt = Process((etProcMode)iProcMode,Strand,MinLen,MaxLen,
			NumIncludeChroms,	// number of chromosome regular expressions to include
			pszIncludeChroms,	// array of include chromosome regular expressions
			NumExcludeChroms,	// number of chromosome expressions to exclude
			pszExcludeChroms,	// array of exclude chromosome regular expressions
		    szInFile,szOutFile);
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
ProcMode2Txt(etProcMode ProcMode)
{
switch(ProcMode) {
	case ePMDefault:
		return((char *)"Default processing");
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

CBEDfile *m_pBedFile; // BED file from which features are to be filtered
int m_hOutFile;		  // file handle for output file contained retained features
int m_BuffIdx;		  // index into pszOutBuff at which to next write
int m_AllocBuffSize;  // pszOutBuff allocated to hold at most this many chars
char *m_pszOutBuff;	  // output buffer

int m_NumIncludeChroms;			// number of chromosomes explicitly defined to be included
char **m_ppszIncludeChroms;		// ptr to array of reg expressions defining chroms to include
int m_NumExcludeChroms;			// number of chromosomes explicitly defined to be excluded
char **m_ppszExcludeChroms;		// ptr to array of reg expressions defining chroms to include
#ifdef _WIN32
Regexp *m_IncludeChromsRE[cMaxIncludeChroms];	// compiled regular expressions
Regexp *m_ExcludeChromsRE[cMaxExcludeChroms];
#else
regex_t m_IncludeChromsRE[cMaxIncludeChroms];	// compiled regular expressions
regex_t m_ExcludeChromsRE[cMaxExcludeChroms];
#endif

char m_szFiltChrom[_MAX_PATH];	// used to cache last chrom processed	
bool m_bFiltChrom;				// and it's filtered status



void
Reset(void)
{
if(m_pBedFile != NULL)
	{
	delete m_pBedFile;
	m_pBedFile = NULL;
	}
if(m_hOutFile != -1)
	{
	close(m_hOutFile);
	m_hOutFile = -1;
	}
if(m_pszOutBuff != NULL)
	{
	delete m_pszOutBuff;
	m_pszOutBuff = NULL;
	}
m_BuffIdx = 0;
m_AllocBuffSize = 0;

m_szFiltChrom[0] = 0;
m_bFiltChrom = false;

m_NumIncludeChroms = 0;			
m_ppszIncludeChroms = NULL;		
m_NumExcludeChroms = 0;			
m_ppszExcludeChroms = NULL;		

}


void
Init(void)
{
m_pBedFile = NULL;
m_pszOutBuff = NULL;
m_hOutFile = -1;
m_ppszIncludeChroms = NULL;
m_ppszExcludeChroms = NULL;
Reset();
}

int
SetChromFilters(int NumIncludeChroms,		 // number of chromosome regular expressions to include
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
m_ppszIncludeChroms = ppszIncludeChroms;	// ptr to array of reg expressions defining chroms to include - overides exclude
m_NumExcludeChroms = NumExcludeChroms;		// number of chromosomes explicitly defined to be excluded
m_ppszExcludeChroms = ppszExcludeChroms;	// ptr to array of reg expressions defining chroms to include

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
ExcludeThisChrom(char *pszChrom)
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



int
Process(etProcMode ProcMode,	// processing mode
		char Strand,			// process for this strand only
		int MinLen,				// retained features must be at least this length
		int MaxLen,				// retained features must be no longer than this length
		int NumIncludeChroms,	// number of chromosome regular expressions to include
		char **ppszIncludeChroms,	 // array of include chromosome regular expressions
		int NumExcludeChroms,		 // number of chromosome expressions to exclude
		char **ppszExcludeChroms,	 // array of exclude chromosome regular expressions
		char *pszInFile,		// BED file containing elements to be filtered
		char *pszOutFile)		// write combined elements into this BED file

{
teBSFrsltCodes Rslt;
int NumBEDsLoaded = 0;
int StartLoci;
int EndLoci;
int FeatLen;
int Score;
char szChrom[128];
char szPrevChrom[128];
char szFeatName[128];
char FeatStrand;
int CurFeatureID;
int NumEls;
int RetainedEls;
int FiltChromCnt;
int FiltStrandCnt;
int FiltMinLenCnt;
int FiltMaxLenCnt;

Init();

if(NumIncludeChroms > 0 || NumExcludeChroms > 0)
	{
	if((Rslt=(teBSFrsltCodes)SetChromFilters(NumIncludeChroms,ppszIncludeChroms,NumExcludeChroms,ppszExcludeChroms)) != eBSFSuccess)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to initialise chrom include/exclude REs");
		Reset();
		return(Rslt);
		}
	}


if((m_pszOutBuff = new char [cAllocBuffSize]) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate %d bytes memory for write buffering",cAllocBuffSize);
	Reset();
	return(eBSFerrMem);
	}
m_AllocBuffSize = cAllocBuffSize;

if((m_pBedFile = new CBEDfile)==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CBEDfile");
	Reset();
	return(eBSFerrObj);
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Starting to load from %s",pszInFile);
if((Rslt=m_pBedFile->Open(pszInFile,eBTAnyBed)) !=eBSFSuccess)
	{
	while(m_pBedFile->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,m_pBedFile->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open '%s' for processing",pszInFile);
	Reset();
	return(eBSFerrOpnFile);
	}

#ifdef _WIN32
if((m_hOutFile = open(pszOutFile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
if((m_hOutFile = open(pszOutFile,O_RDWR | O_CREAT |O_TRUNC, S_IREAD | S_IWRITE))==-1)
#endif
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create or truncate %s - %s",pszOutFile,strerror(errno));
	Reset();
	return(eBSFerrCreateFile);
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Load completed, iterating features and filtering...");

// now iterate over the features, filtering as may be appropriate
szPrevChrom[0] = '\0';
CurFeatureID = 0;
Rslt = eBSFSuccess;

NumEls = 0;
RetainedEls = 0;
FiltChromCnt = 0;
FiltStrandCnt = 0;
FiltMinLenCnt = 0;
FiltMaxLenCnt = 0;

while(Rslt == eBSFSuccess && (CurFeatureID = m_pBedFile->GetNextFeatureID(CurFeatureID)) > 0)
	{
	NumEls += 1;
	m_pBedFile->GetFeature(CurFeatureID,	// feature instance identifier
				szFeatName,				// where to return feature name
				szChrom,				// where to return chromosome name
				&StartLoci,				// where to return feature start on chromosome (0..n) 
				&EndLoci,				// where to return feature end on chromosome
 				&Score,					// where to return score
 				&FeatStrand);				// where to return strand

	// exclude this chrom?

	if(ExcludeThisChrom(szChrom))
		{
		FiltChromCnt += 1;
		continue;
		}
	if(Strand != '*')
		{
		if(FeatStrand != Strand)
			{
			FiltStrandCnt += 1;
			continue;
			}
		}

	FeatLen = 1 + EndLoci - StartLoci;
	if(MinLen > FeatLen)
		{
		FiltMinLenCnt += 1;
		continue;
		}
	if(MaxLen < FeatLen)
		{
		FiltMaxLenCnt += 1;
		continue;
		}		

	if((m_BuffIdx + 0x07fff) > m_AllocBuffSize)
		{
		CUtility::SafeWrite(m_hOutFile,m_pszOutBuff,m_BuffIdx);
		m_BuffIdx = 0;
		}
	Rslt = (teBSFrsltCodes)m_pBedFile->GetBEDFormatedFeature(CurFeatureID,m_AllocBuffSize - m_BuffIdx,&m_pszOutBuff[m_BuffIdx]);
	if(Rslt < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"GetBEDFormatedFeature failed");
		Reset();
		return(Rslt);
		}
	m_BuffIdx += (int)Rslt;
	m_BuffIdx += sprintf(&m_pszOutBuff[m_BuffIdx],"\n");
	RetainedEls += 1;
	Rslt = eBSFSuccess;
	}
if(m_BuffIdx > 0)
	CUtility::SafeWrite(m_hOutFile,m_pszOutBuff,m_BuffIdx);
close(m_hOutFile);
m_hOutFile = -1;

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Retained %d features from %d processed",RetainedEls,NumEls);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Filtered out %d by chrom, %d by strand, %d by minlen and %d by maxlen",FiltChromCnt,FiltStrandCnt,FiltMinLenCnt,FiltMaxLenCnt);
Reset();
return((teBSFrsltCodes)Rslt);
}


