// gennucstats.cpp : Defines the entry point for the console application.
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

#include "IncExclChroms.h"

const unsigned int cProgVer = 101;			// increment with each release
const int cDfltDyadOfs = 73;				// default offset to add to the loci start for the dyad loci
const int cMaxDyadOfs  = 500;				// maximum offset to add to the loci start for the dyad loci
const int cDfltWindDyad = 5;				// accept sample dyads as matching if within this many bases of background reference dyad
const int cMaxWindDyad = 50;				// maximum dyad flanking window

// processing modes
typedef enum TAG_ePMode {		
	ePMNDyadRegionDist,						// dyad regional distribution
	ePMNDyadOverlap,						// dyad overlap probabilities distribution
	ePMplaceholder							// used to set the enumeration range
	} etPMode;

int
Process(etPMode PMode,				// processing mode
		int	BkgDyadOfs,				// offset background loci by this many bases to derive dyad loci
		int	SmplDyadOfs,			// offset sample loci by this many bases to derive dyad loci
		int	WindDyad,				// accept sample dyad as matching background loci if within +/- this many bases
		int	NumIncludeChroms,		// number of chromosome regular expressions to include
		char **ppszIncludeChroms,	// array of include chromosome regular expressions
		int	NumExcludeChroms,		// number of chromosome expressions to exclude
		char **ppszExcludeChroms,	// array of exclude chromosome regular expressions
		char *pszDyadFile,			// CSV file containing background dyad loci
		char *pszSampleFile,		// CSV file containing sample dyad loci
		char *pszRsltsFile); 		// output results to this file



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
	return _T("gennucstats");
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
int LenChromList;


etPMode PMode;				// processing mode

int BkgDyadOfs;				// offset background loci by this many bases to derive dyad loci
int SmplDyadOfs;			// offset sample loci by this many bases to derive dyad loci
int WindDyad;				// accept sample dyad as matching background loci if within +/- this many bases

char szDyadFile[_MAX_PATH];	// CSV file containing background dyad loci
char szSampleFile[_MAX_PATH];	// CSV file containing sample dyad loci
char szRsltsFile[_MAX_PATH];	// output results to this file

int NumIncludeChroms;
char *pszIncludeChroms[cMaxIncludeChroms];
int NumExcludeChroms;
char *pszExcludeChroms[cMaxExcludeChroms];


// command line args
struct arg_lit  *help    = arg_lit0("hH","help",                "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_int *ScreenLogLevel=arg_int0("S", "ScreenLogLevel",	"<int>","Level of diagnostics written to screen 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"diagnostics log file");

struct arg_int *pmode = arg_int0("m","mode","<int>",		    "processing mode: 0 - background dyad regional distribution, 1 - Sample dyads overlap onto background (default = 0)");

struct arg_int *bkgdyadofs = arg_int0("b","bkgdyadofs","<int>",		"offset background loci by this many bases to derive dyad loci (default = 73)");
struct arg_int *smpldyadofs = arg_int0("s","smpldyadofs","<int>",	"offset sample loci by this many bases to derive dyad loci (default = 73)");
struct arg_int *winddyad = arg_int0("w","winddyad","<int>",		"accept sample dyad as matching background loci if within +/- this many bases (default = 5)");

struct arg_file *dyadfile = arg_file1("i","infile","<file>",	"CSV file containing background dyad loci");
struct arg_file *samplefile = arg_file0("I","sample","<file>",	"CSV file containing sample dyad loci");
struct arg_file *rsltsfile = arg_file1("o","outfile","<file>",	"output results to this file");

struct arg_str  *ExcludeChroms = arg_strn("Z","chromexclude","<string>",0,cMaxExcludeChroms,"high priority - regular expressions defining species.chromosomes to exclude");
struct arg_str  *IncludeChroms = arg_strn("z","chromeinclude","<string>",0,cMaxIncludeChroms,"low priority - regular expressions defining species.chromosomes to include");


struct arg_end *end = arg_end(20);

void *argtable[] = {help,version,FileLogLevel,ScreenLogLevel,LogFile,
					pmode,bkgdyadofs,smpldyadofs,winddyad,dyadfile,samplefile,rsltsfile,ExcludeChroms,IncludeChroms,
					end};

char **pAllArgs;
int argerrors;
argerrors = CUtility::arg_parsefromfile(argc,(char **)argv,&pAllArgs);
if(argerrors >= 0)
	argerrors = arg_parse(argerrors,pAllArgs,argtable);

/* special case: '--help' takes precedence over error reporting */
if (help->count > 0)
        {
		printf("\n%s ", gszProcName);
        arg_print_syntax(stdout,argtable,"\n");
        arg_print_glossary(stdout,argtable,"  %-25s %s\n");
		printf("\nNote: Parameters can be entered into a parameter file, one parameter per line.");
		printf("\n      To invoke this parameter file then precede its name with '@'");
		printf("\n      e.g. %s @myparams.txt\n\n",gszProcName);
		exit(1);
        }

    /* special case: '--version' takes precedence error reporting */
if (version->count > 0)
        {
		printf("\n%s Version %d.%2.2d",gszProcName,cProgVer/100,cProgVer%100);
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

	PMode = (etPMode)(pmode->count ? pmode->ival[0] : ePMNDyadRegionDist);
	if(PMode < 0 || PMode >= ePMplaceholder)
		{
		printf("\nError: Processing mode '-m%d' specified outside of range %d..%d",PMode,0,(int)ePMplaceholder-1);
		exit(1);
		}

	BkgDyadOfs = bkgdyadofs->count ? bkgdyadofs->ival[0] : cDfltDyadOfs;
	if(abs(BkgDyadOfs) > cMaxDyadOfs)
		{
		printf("\nError: Dyad offset specified as '-b%d' must be in range +/-0..%d",BkgDyadOfs,cMaxDyadOfs);
		exit(1);
		}

	if(PMode == ePMNDyadOverlap)
		{
		if(!samplefile->count)
			{
			printf("\nError: Expected CSV file containing sample dyad loci to be specified with '-I<file>");
			exit(1);
			}
		strncpy(szSampleFile,samplefile->filename[0],_MAX_PATH);
		szSampleFile[_MAX_PATH-1] = '\0';
		CUtility::TrimQuotedWhitespcExtd(szSampleFile);
		if(szSampleFile[0] == '\0')
			{
			printf("\nError: No sample dyad file ('-I<file> option') specified after removal of leading/trailing quotes and whitespace");
			exit(1);
			}
		SmplDyadOfs = smpldyadofs->count ? smpldyadofs->ival[0] : cDfltDyadOfs;
		if(abs(SmplDyadOfs) > cMaxDyadOfs)
			{
			printf("\nError: Dyad offset specified as '-s%d' must be in range +/-0..%d",SmplDyadOfs,cMaxDyadOfs);
			exit(1);
			}

		WindDyad = winddyad->count ? winddyad->ival[0] : cDfltWindDyad;
		if(WindDyad < 0 || WindDyad > cMaxWindDyad)
			{
			printf("\nError: Dyad matching window flanks specified as '-w%d' must be in range 0..%d",WindDyad,cMaxWindDyad);
			exit(1);
			}
		}
	else
		{
		szSampleFile[0] = '\0';
		SmplDyadOfs = 0;
		WindDyad = 0;
		}

	strncpy(szDyadFile,dyadfile->filename[0],_MAX_PATH);
	szDyadFile[_MAX_PATH-1] = '\0';
	CUtility::TrimQuotedWhitespcExtd(szDyadFile);
	if(szDyadFile[0] == '\0')
		{
		printf("\nError: No dyad file ('-i<file> option') specified after removal of leading/trailing quotes and whitespace");
		exit(1);
		}

	strncpy(szRsltsFile,rsltsfile->filename[0],_MAX_PATH);
	szRsltsFile[_MAX_PATH-1] = '\0';
	CUtility::TrimQuotedWhitespcExtd(szRsltsFile);
	if(szRsltsFile[0] == '\0')
		{
		printf("\nError: No results file ('-o<file> option') specified after removal of leading/trailing quotes and whitespace");
		exit(1);
		}

	NumIncludeChroms = IncludeChroms->count;
	for(Idx=0;Idx < IncludeChroms->count; Idx++)
		{
		LenChromList = (int)strlen(IncludeChroms->sval[Idx]);
		pszIncludeChroms[Idx] = new char [LenChromList+1];
		strcpy(pszIncludeChroms[Idx],IncludeChroms->sval[Idx]);
		CUtility::TrimQuotes(pszIncludeChroms[Idx]);
		}

	NumExcludeChroms = ExcludeChroms->count;
	for(Idx=0;Idx < ExcludeChroms->count; Idx++)
		{
		LenChromList = (int)strlen(ExcludeChroms->sval[Idx]);
		pszExcludeChroms[Idx] = new char [LenChromList+1];
		strcpy(pszExcludeChroms[Idx],ExcludeChroms->sval[Idx]);
		CUtility::TrimQuotes(pszExcludeChroms[Idx]);
		}

			// now that command parameters have been parsed then initialise diagnostics log system
	if(!gDiagnostics.Open(szLogFile,(etDiagLevel)iScreenLogLevel,(etDiagLevel)iFileLogLevel,true))
		{
		printf("\nError: Unable to start diagnostics subsystem.");
		if(szLogFile[0] != '\0')
			printf(" Most likely cause is that logfile '%s' can't be opened/created",szLogFile);
		exit(1);
		}

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Version: %d.%2.2d Processing parameters:",cProgVer/100,cProgVer%100);
	const char *pszProcMode;
	switch(PMode) {
		case ePMNDyadRegionDist:
			pszProcMode = "background dyad regional distribution";
			break;

		case ePMNDyadOverlap:
			pszProcMode = "sample dyads overlap onto background dyads";
			break;

		};
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Processing mode: '%s'",pszProcMode);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"offset background loci by: %d",BkgDyadOfs);
	if(PMode == ePMNDyadOverlap)
		{
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"offset sample loci by: %d",SmplDyadOfs);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"window flanking: %d",WindDyad);
		}

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"input background loci file: '%s'",szDyadFile);

	if(PMode == ePMNDyadOverlap)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"input sample loci file: '%s'",szSampleFile);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"results file: '%s'",szRsltsFile);

	for(Idx = 0; Idx < NumIncludeChroms; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"reg expressions defining chroms to include: '%s'",pszIncludeChroms[Idx]);
	for(Idx = 0; Idx < NumExcludeChroms; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"reg expressions defining chroms to exclude: '%s'",pszExcludeChroms[Idx]); 

#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	gStopWatch.Start();
	Rslt = Process(PMode,					// processing mode
					BkgDyadOfs,				// offset background loci by this many bases to derive dyad loci
					SmplDyadOfs,			// offset sample loci by this many bases to derive dyad loci
					WindDyad,				// accept sample dyad as matching background loci if within +/- this many bases
					NumIncludeChroms,
					pszIncludeChroms,
					NumExcludeChroms,
					pszExcludeChroms,
					szDyadFile,				// CSV file containing background dyad loci
					szSampleFile,			// CSV file containing sample dyad loci
					szRsltsFile);			// output results to this file
	gStopWatch.Stop();
	Rslt = Rslt >=0 ? 0 : 1;
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Exit code: %d Total processing time: %s",Rslt,gStopWatch.Read());
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


const int cDyadNumAlloc = 5000000;	// allocate for dyads in this number increments


typedef struct TAG_sDyadLoci {
	int DyadID;				// uniquely identifies this dyad
	int ChromID;			// dyad is on this chrom
	int Loci;				// and at this loci
	int Map;				// this dyad is mapped to this region
	int MatchID;			// this dyad is matched to this other dyad
	int Rank;				// used to determine the median distance between dyads
} tsDyadLoci;

const int cNumDistBins = 2000;	// bin dyad distances into this many bins
typedef struct TAG_sDyadLocii {
	tsDyadLoci *pDyads;
	int AllocdDyads;
	int NumDyads;
	int NumOverlaps;			// number of overlap dyads
	int SortMode;			// 0 if unsorted, 1 if sorted by ascending chrom.loci, 2 if sorted by ascending rank
	int LociOfs;			// loci offset used to derive dyad locii
	double Mean;			// mean distance between dyads
	double Median;			// median distance between dyads
	double StdDev;			// stddev of distances
	int DistCnts[cNumDistBins];		// distance distribution counts
	} tsDyadLocii;

int DyadicOverlaps(tsDyadLocii *pADyadics,tsDyadLocii *pBDyadics,int WindDyad);
int DyadicDistances(tsDyadLocii *pDyadics);

static int SortDyadsLoci(const void *arg1, const void *arg2);
static int SortDyadsRank(const void *arg1, const void *arg2);

void Init(void);
void Reset(void);
int DumpDistribution(tsDyadLocii *pDyadics,char *pszDyadFile);

int m_hOutFile;					// output results file handle
CIncExclChroms m_IncExclChroms;		// include/exclude chromosome processing classs instance
CCSVFile *m_pCSV;			// file by which background/sample loci is being loaded
tsDyadLocii m_RefDyads;		// background, or reference, dyads
tsDyadLocii m_SampleDyads;	// sample dyads

void
Init(void)
{
m_pCSV = NULL;
m_hOutFile = -1;
memset(&m_RefDyads,0,sizeof(m_RefDyads));
memset(&m_SampleDyads,0,sizeof(m_SampleDyads));
Reset();
}

void
Reset(void)
{
if(m_hOutFile != -1)
	{
	close(m_hOutFile);
	m_hOutFile = -1;
	}
if(m_pCSV != NULL)
	{
	delete m_pCSV;
	m_pCSV = NULL;
	}
if(m_RefDyads.pDyads != NULL)
	{
	free(m_RefDyads.pDyads);
	m_RefDyads.pDyads = NULL;
	}	
if(m_SampleDyads.pDyads != NULL)
	{
	free(m_SampleDyads.pDyads);
	m_SampleDyads.pDyads = NULL;
	}
m_IncExclChroms.Reset();
memset(&m_RefDyads,0,sizeof(m_RefDyads));
memset(&m_SampleDyads,0,sizeof(m_SampleDyads));
}


bool			
AdjLoci(bool bOnMinStrand,		// true if element is on '-' strand
		int *pStartLoci,		// starting loci on '+' strand
		int *pEndLoci,			// ending loci on '+' strand
		int OfsLoci)			// offset loci by
{
int Ofs;
if(*pStartLoci < 0 || *pStartLoci > *pEndLoci)
	return(false);
if(!OfsLoci)
	return(true);

if(bOnMinStrand)
	Ofs = OfsLoci * -1;
else
	Ofs = OfsLoci;
*pStartLoci += Ofs;
*pEndLoci += Ofs;

if((bOnMinStrand && *pEndLoci < 0) || (!bOnMinStrand && *pStartLoci < 0))
	return(false);
return(true);
}

// LoadDyads
// Loads putative dyads into memory applying LociOfs to the start loci to derive the dyad loci
int
LoadDyads(int LociOfs,tsDyadLocii *pDyadLocii,char *pszLociFile)
{
int Rslt;
tsDyadLoci *pDyadLoci;
int NumElsProcessed;
int NumFields;

int SrcID;
char *pszElType;
char *pszRefSpecies;
int StartLoci;
int EndLoci;
int Len;
char *pszChrom;
char *pszStrand;
bool bOnMinStrand;
int ChromID;

int NumExclChroms;
int NumExclOfs;

if(m_pCSV != NULL)
	delete m_pCSV;
if(pDyadLocii->pDyads == NULL)
	pDyadLocii->AllocdDyads = 0;
pDyadLocii->NumDyads = 0;
pDyadLocii->LociOfs = LociOfs;

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Loading dyads from '%s'...",pszLociFile);

if((m_pCSV = new CCSVFile) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CCSVfile");
	return(eBSFerrObj);
	}
if((Rslt=m_pCSV->Open(pszLociFile))!=eBSFSuccess)
	{
	while(m_pCSV->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,m_pCSV->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open file: %s",pszLociFile);
	Reset();
	return(Rslt);
	}

if(pDyadLocii->pDyads == NULL)
	{
	if((pDyadLocii->pDyads = (tsDyadLoci *)malloc(sizeof(tsDyadLoci) * cDyadNumAlloc))==NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to allocate memory (%d bytes) for dyad loci",sizeof(tsDyadLoci) * cDyadNumAlloc);
		Reset();
		return(eBSFerrMem);
		}
	pDyadLocii->AllocdDyads = cDyadNumAlloc;
	}

NumElsProcessed = 0;
NumExclChroms = 0;
NumExclOfs = 0;
while((Rslt=m_pCSV->NextLine()) > 0)	// onto next line containing fields
	{
	if(!NumElsProcessed && m_pCSV->IsLikelyHeaderLine())
		continue;

	NumElsProcessed++;
	NumFields = m_pCSV->GetCurFields();
	if(NumFields < 7)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expected 7+ fields in '%s', GetCurFields() returned '%d'",pszLociFile,NumFields);
		Reset();
		return(eBSFerrFieldCnt);
		}

	m_pCSV->GetText(4,&pszChrom);
	if((Rslt = m_IncExclChroms.IncludeThisChrom(pszChrom))<=0)
		{
		if(Rslt < 0)
			break;
		NumExclChroms += 1;
		continue;
		}
	ChromID = Rslt;
	m_pCSV->GetInt(1,&SrcID);
	m_pCSV->GetText(2,&pszElType);
	m_pCSV->GetText(3,&pszRefSpecies);

	m_pCSV->GetInt(5,&StartLoci);
	m_pCSV->GetInt(6,&EndLoci);
	m_pCSV->GetInt(7,&Len);
	if(StartLoci > EndLoci)
		{
		int tmp = StartLoci;
		StartLoci = EndLoci;
		EndLoci = tmp;
		}
	bOnMinStrand = false;
	if(NumFields >= 8)					// check if strand has been specified
		{
		m_pCSV->GetText(8,&pszStrand);
		if(pszStrand[0] == '-')		// assume anything other than '-' is on the plus strand
			bOnMinStrand = true;
		}
	if(!AdjLoci(bOnMinStrand,&StartLoci,&EndLoci,LociOfs))
		{
		NumExclOfs += 1;
		continue;
		}

	if(pDyadLocii->NumDyads == pDyadLocii->AllocdDyads)	// time to realloc?
		{
		int ReqMem = sizeof(tsDyadLoci) * (cDyadNumAlloc + pDyadLocii->AllocdDyads);
		if((pDyadLoci = (tsDyadLoci *)realloc(pDyadLocii->pDyads,ReqMem))==NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to allocate memory (%d bytes) for dyad loci",ReqMem);
			Reset();
			return(eBSFerrMem);
			}
		pDyadLocii->pDyads = pDyadLoci;
		pDyadLocii->AllocdDyads += cDyadNumAlloc;
		}
	pDyadLoci = &pDyadLocii->pDyads[pDyadLocii->NumDyads++];
	pDyadLoci->ChromID = ChromID;
	pDyadLoci->DyadID = SrcID;
	pDyadLoci->Loci = bOnMinStrand ? EndLoci : StartLoci;
	pDyadLoci->Map = 0;
	}
if(Rslt >= 0)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processed %d elements and accepted %d, filtered out %d bcause of chrom and %d because offset dyad would be before chrom start",
					 NumElsProcessed,pDyadLocii->NumDyads,NumExclChroms,NumExclOfs);
	if(pDyadLocii->NumDyads > 1)
		{
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Sorting %d dyads...",pDyadLocii->NumDyads);
		qsort(pDyadLocii->pDyads,pDyadLocii->NumDyads,sizeof(tsDyadLoci),SortDyadsLoci);
		pDyadLocii->SortMode = 1;
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Sort completed");
		}
	}
return(Rslt >= 0 ? pDyadLocii->NumDyads : Rslt);
}

int
Process(etPMode PMode,				// processing mode
		int	BkgDyadOfs,				// offset background loci by this many bases to derive dyad loci
		int	SmplDyadOfs,			// offset sample loci by this many bases to derive dyad loci
		int	WindDyad,				// accept sample dyad as matching background loci if within +/- this many bases
		int	NumIncludeChroms,		// number of chromosome regular expressions to include
		char **ppszIncludeChroms,	// array of include chromosome regular expressions
		int	NumExcludeChroms,		// number of chromosome expressions to exclude
		char **ppszExcludeChroms,	// array of exclude chromosome regular expressions
		char *pszDyadFile,			// CSV file containing background dyad loci
		char *pszSampleFile,		// CSV file containing sample dyad loci
		char *pszRsltsFile) 		// output results to this file
{
int Rslt;
Init();

#ifdef _WIN32
if((m_hOutFile = open(pszRsltsFile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
if((m_hOutFile = open(pszRsltsFile, O_RDWR | O_CREAT | O_TRUNC, S_IREAD | S_IWRITE ))==-1)
#endif
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create results output file: %s - %s",pszRsltsFile,strerror(errno));
	Reset();
	return(eBSFerrCreateFile);
	}

if((Rslt=m_IncExclChroms.InitChromExclusion(NumIncludeChroms,ppszIncludeChroms,NumExcludeChroms,ppszExcludeChroms))!=eBSFSuccess)
	{
	Reset();
	return(Rslt);
	}
if((Rslt = LoadDyads(BkgDyadOfs,&m_RefDyads,pszDyadFile)) <= 0)
	{
	Reset();
	return(Rslt);
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Generating basic stats on background/reference dyads from '%s'",pszDyadFile);
DyadicDistances(&m_RefDyads);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Generating basic stats on background/reference dyads completed");
if(PMode == ePMNDyadOverlap)
	{
	if((Rslt = LoadDyads(SmplDyadOfs,&m_SampleDyads,pszSampleFile)) <= 0)
		{
		Reset();
		return(Rslt);
		}
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Generating basic stats on sample dyads from '%s'",pszSampleFile);
	DyadicDistances(&m_SampleDyads);
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Generating basic stats on sample dyads completed");
	
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Generating basic stats on dyads overlapping from '%s' onto '%s'",pszSampleFile,pszDyadFile);
	DyadicOverlaps(&m_RefDyads,&m_SampleDyads,WindDyad);
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Generating basic stats on dyads overlapping completed, %d from sample (%d) overlap reference (%d) within +/-%dnt",
															m_RefDyads.NumOverlaps,m_SampleDyads.NumDyads,m_RefDyads.NumDyads,WindDyad);
	}
DumpDistribution(&m_RefDyads,pszDyadFile);
if(PMode == ePMNDyadOverlap)
	DumpDistribution(&m_SampleDyads,pszSampleFile);
Reset();
return(Rslt);
}

int
DumpDistribution(tsDyadLocii *pDyadics,char *pszDyadFile)
{
char szBuff[10000];
int BuffOfs;
int SumCnt;
int Idx;
BuffOfs = sprintf(szBuff,"\"LociFile\",\"LociOfs\",\"NumDyads\",\"MeanDist\",\"MedianDist\",\"StdDevDist\",\"NumOverlaps\"");
BuffOfs += sprintf(&szBuff[BuffOfs],"\n\"%s\",%d,%d,%1.3f,%1.3f,%1.3f,%d",
				   pszDyadFile,pDyadics->LociOfs,pDyadics->NumDyads,pDyadics->Mean,pDyadics->Median,pDyadics->StdDev,pDyadics->NumOverlaps);
BuffOfs += sprintf(&szBuff[BuffOfs],"\n\n\"InterDyadDist\",\"Count\"");
SumCnt = 0;
for(Idx = 1; Idx < cNumDistBins; Idx++)
	{
	if(Idx <= 999)
		{
		BuffOfs += sprintf(&szBuff[BuffOfs],"\n%d,%d",Idx,pDyadics->DistCnts[Idx]);
		continue;
		}
	SumCnt += pDyadics->DistCnts[Idx];

	if(Idx == (cNumDistBins-1))
		{
		BuffOfs += sprintf(&szBuff[BuffOfs],"\n1000+,%d",SumCnt);
		SumCnt = 0;
		}
	}
BuffOfs += sprintf(&szBuff[BuffOfs],"\n\n");
CUtility::SafeWrite(m_hOutFile,szBuff,BuffOfs);
BuffOfs = 0;
return(0);
}	

int
DyadicOverlaps(tsDyadLocii *pADyadics,tsDyadLocii *pBDyadics,int WindDyad)
{
int AIdx;
int BIdx;
tsDyadLoci *pADyad;
tsDyadLoci *pBDyad;

if(pADyadics->SortMode != 1)
	{
	qsort(pADyadics->pDyads,pADyadics->NumDyads,sizeof(tsDyadLoci),SortDyadsLoci);
	pADyadics->SortMode = 1;
	}
if(pBDyadics->SortMode != 1)
	{
	qsort(pBDyadics->pDyads,pBDyadics->NumDyads,sizeof(tsDyadLoci),SortDyadsLoci);
	pBDyadics->SortMode = 1;
	}

pADyadics->NumOverlaps = 0;
pADyad = pADyadics->pDyads;
for(AIdx = 0; AIdx < pADyadics->NumDyads; AIdx++,pADyad++)
	pADyad->MatchID = 0;
pBDyadics->NumOverlaps = 0;
pBDyad = pBDyadics->pDyads;
for(BIdx = 0; BIdx < pBDyadics->NumDyads; BIdx++,pBDyad++)
	pBDyad->MatchID = 0;

pADyad = pADyadics->pDyads;
pBDyad = pBDyadics->pDyads;
AIdx = BIdx = 0;
do	{
	if(pADyad->ChromID > pBDyad->ChromID)
		{
		pBDyad++;
		BIdx++;
		continue;
		}
	else
		if(pADyad->ChromID < pBDyad->ChromID)
			{
			pADyad++;
			AIdx++;
			continue;
			}
	// same chromsome..
	if(pBDyad->Loci > (pADyad->Loci + WindDyad))
		{
		pADyad++;
		AIdx++;
		continue;
		}
	if(pBDyad->Loci < pADyad->Loci - WindDyad)
		{
		pBDyad++;
		BIdx++;
		continue;
		}

	// match
	pADyadics->NumOverlaps += 1;
	pADyad->MatchID = pBDyad->DyadID;
	pBDyad->MatchID = pADyad->DyadID;
	pADyad++;
	pBDyad++;
	AIdx++;
	BIdx++;
	}
while(AIdx < pADyadics->NumDyads && BIdx < pBDyadics->NumDyads);

return(pADyadics->NumOverlaps);
}

// DyadicDistances
// generate stats on distances between dyads
int
DyadicDistances(tsDyadLocii *pDyadics)
{
int Idx;
int NumProcessed;
int CurChromID;
int DistIdx;
double Variance;
tsDyadLoci *pDyad;

memset(pDyadics->DistCnts,0,sizeof(pDyadics->DistCnts));
pDyadics->Mean = 0.0f;
pDyadics->Median = 0.0f;
pDyadics->StdDev = 0.0f;

if(pDyadics->SortMode != 1)
	{
	qsort(pDyadics->pDyads,pDyadics->NumDyads,sizeof(tsDyadLoci),SortDyadsLoci);
	pDyadics->SortMode = 1;
	}

NumProcessed = 0;
CurChromID = -1;
pDyad = pDyadics->pDyads;
for(Idx = 0; Idx < pDyadics->NumDyads; Idx++, pDyad++)
	{
	if(pDyad->ChromID != CurChromID)		// new chrom?
		{
		pDyad->Rank = 0;
		CurChromID = pDyad->ChromID;
		continue;
		}
	DistIdx = pDyad->Loci - pDyad[-1].Loci;
	pDyadics->Mean += DistIdx;
	pDyad->Rank = DistIdx;
	if(DistIdx >= cNumDistBins)
		DistIdx = cNumDistBins-1;
	pDyadics->DistCnts[DistIdx] += 1;
	NumProcessed += 1;
	}
pDyadics->Mean /= NumProcessed;
CurChromID = -1;
pDyad = pDyadics->pDyads;
Variance = 0.0f;
for(Idx = 0; Idx < pDyadics->NumDyads; Idx++, pDyad++)
	{
	if(pDyad->ChromID != CurChromID)		// new chrom?
		{
		CurChromID = pDyad->ChromID;
		continue;
		}
	DistIdx = pDyad->Loci - pDyad[-1].Loci;
	Variance += pow((double)DistIdx - pDyadics->Mean,2.0);
	}
pDyadics->StdDev = sqrt(Variance);

qsort(pDyadics->pDyads,pDyadics->NumDyads,sizeof(tsDyadLoci),SortDyadsRank);
pDyadics->SortMode = 2;
pDyad = &pDyadics->pDyads[(pDyadics->NumDyads+1)/2];
pDyadics->Median = pDyad->Rank;
return(0);
}

// SortDyadsLoci
// Sort by ascending chrom, loci
static int
SortDyadsLoci(const void *arg1, const void *arg2)
{
tsDyadLoci *pEl1 = (tsDyadLoci *)arg1;
tsDyadLoci *pEl2 = (tsDyadLoci *)arg2;

if(pEl1->ChromID < pEl2->ChromID )
		return(-1);
if(pEl1->ChromID > pEl2->ChromID )
	return(1);
if(pEl1->Loci < pEl2->Loci )
		return(-1);
if(pEl1->Loci > pEl2->Loci )
	return(1);
return(0);
}

// SortDyadsRank
// Sort by ascending Rank
static int
SortDyadsRank(const void *arg1, const void *arg2)
{
tsDyadLoci *pEl1 = (tsDyadLoci *)arg1;
tsDyadLoci *pEl2 = (tsDyadLoci *)arg2;

if(pEl1->Rank < pEl2->Rank)
		return(-1);
if(pEl1->Rank > pEl2->Rank)
	return(1);
return(0);
}