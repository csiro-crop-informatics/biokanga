// ufilter.cpp : Defines the entry point for the console application.
// use to filter csv loci files
// filter by chromosome, regions, and comformational twist/groove characteristics
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

const unsigned int cProgVer = 104;		// increment with each release

const unsigned int cMaxProcSeqLen = 2000; // can only handle sequences of upto this length
const double cMinBkgndGroove = 11.00f;  // min background minor groove width
const double cDfltBkgndGroove = 11.20f; // default background minor groove width
const double cMaxBkgndGroove = 11.40f;  // max background minor groove width

const double cMinDyadratio = 1.000f;	// min dyad grooves must be at least this ratio to background
const double cDfltDyadratio = 1.030f;	// default dyad grooves must be at least this ratio to background
const double cMaxDyadratio = 1.050f;	// max dyad grooves must be at least this ratio to background

const double cMinDyad2ratio = 1.000f;	// min immediately flanking grooves must be at least this ratio to background
const double cDfltDyad2ratio = 1.020f;	// default immediately flanking grooves must be at least this ratio to background
const double cMaxDyad2ratio = 1.040f;	// max immediately flanking grooves must be at least this ratio to background

const double cMinDyad3ratio = 1.000f;	// min remainder of flanking grooves must be at least this ration to background
const double cDfltDyad3ratio = 1.015f;	// default remainder of flanking grooves must be at least this ration to background
const double cMaxDyad3ratio = 1.030f;	// max remainder of flanking grooves must be at least this ration to background

const int cAllocOutBuff = 1024000;		// results output buffer allocation size

// processing modes
typedef enum TAG_ePMode {		
	ePMDefault,							// default is to filter on chromosome/region/strand only
	ePMFiltNucConf,						// filter on chromosome/region/strand and also on putative nucleosome conformation
	ePMplaceholder						// used to set the enumeration range
	} etPMode;

#pragma pack(1)
typedef struct TAG_sElement {
	int ElementID;					// uniquely identifies this element
	int SrcID;						// original element source identifier
	int TypeID;					    // original element type
	int SpeciesID;					// original species 
	int NumInstances;				// number of instances of this element (useful to hold number of reads)
	int ChromID;					// element is on this chromosome
	int StartLoci;					// loci at which element starts
	int EndLoci;					// loci at which element ends
	char Strand;					// is on '+' or '-' strand
	int Score;						// score for higest scoring dyad
	int DyadIdx;					// offset from element start of highest scoring dyad
	int NumDyads;					// number of dyads located
} tsElement;
#pragma pack()



int
Process(etPMode PMode,						// processing mode
		char Strand,						// only retain elements on this strand
		int OfsLoci,						// offset element start loci by this many nt
		int DeltaLen,						// change element length by this many nt
		int MinLength,						// slough elements of less than this length
		int TruncLength,					// truncate elements to be a maximum of this length
		double BkgndGroove,					// background minor groove 
		double DyadratioThres,				// dyad grooves must be at least this ratio to background
		double Dyad2ratioThres,				// immediately flanking grooves must be at least this ratio to background
		double Dyad3ratioThres,				// remainder of flanking grooves must be at least this ration to background
		char *pszInGenomefile,				// bioseq genome file
		char *pszInConfFile,				// file containing conformation characteristics
		char *pszLociFile,					// read element loci ffrom this file
		char *pszRsltsFile,					// write filtered and translated element loci to this file
		int	NumIncludeChroms,				// number of chromosome regular expressions to include
		char **ppszIncludeChroms,			// array of include chromosome regular expressions
		int	NumExcludeChroms,				// number of chromosome expressions to exclude
		char **ppszExcludeChroms);			// array of exclude chromosome regular expressions

void Init(void);
void Reset(void);
int AddSpecies(char *pszSpecies);
int AddElType(char *pszElType);
int SortElements(const void *arg1, const void *arg2);

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
	return _T("ufilter");
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

char Strand;				// if '-' or '+' then only process that specific strand
int MinLength;				// elements must be of at least this length
int OfsLoci;				// offset loci by this many bases
int DeltaLen;				// change element lengths by this many bases
int TruncLength;			// and truncate to be no longer than this length

double BkgndGroove;			// background minor groove (if not specified then defaults to 11.12)
double Dyadratio;
double Dyad2ratio;
double Dyad3ratio;

char szLocifile[_MAX_PATH];		// input csv loci file
char szOutfile[_MAX_PATH];      // output loci file
char szGenomefile[_MAX_PATH];	// bioseq genome file
char szConfFile[_MAX_PATH];		// structural conformation characteristics file

int LenChromList;
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

struct arg_int *pmode = arg_int0("m","mode","<int>",		    "processing mode: 0 - filter on chrom/region/strand only, 1 - additionally filter on putative nucleosome conformation (default = 0)");

struct arg_str *strand=arg_str0("s", "strand","<str>",          "process for this strand '+' or '-' only (default is to process both)");

struct arg_int  *minlength = arg_int0("l","minlength","<int>",	"minimum element length (default 30)");
struct arg_int  *trunclength = arg_int0("T","truncatelength","<int>","truncate elements to be no  longer than this length (default 300)");
struct arg_int  *ofsloci   = arg_int0("u","offset","<int>",	    "offset element loci by this many bases, -2048 to +2048 (default 0)");
struct arg_int  *deltalen  = arg_int0("U","deltalen","<int>",	"delta element lengths by this many bases, -2048 to +2048 (default 0)");

struct arg_dbl *bkgndgroove=arg_dbl0("b", "bkgndgroove",		"<double>","background minor groove (default 11.1200)");
struct arg_dbl *dyadratio=arg_dbl0("d", "dyadratio",			"<double>","dyad minor grooves must be at least this ratio to background (default 1.030");
struct arg_dbl *dyad2ratio=arg_dbl0("D", "dyad2ratio",			"<double>","immediately flanking minor grooves must be at least this ratio to background (default 1.020");
struct arg_dbl *dyad3ratio=arg_dbl0("e", "dyad3ratio",			"<double>","remainder flanking minor grooves must be at least this ratio to background (default 1.015");

struct arg_file *inlocifile = arg_file1("i","inloci","<file>",	"sequence loci CSV file");
struct arg_file *outfile = arg_file1("o","out","<file>",		"output filtered loci to this file");

struct arg_file *inGenomefile = arg_file0("g","in","<file>",	"input from this bioseq genome assembly file");
struct arg_file *inConfFile = arg_file0("I","conf","<file>",	"in mode 1, input conformation characteristics from this file");

struct arg_str  *ExcludeChroms = arg_strn("Z","chromexclude","<string>",0,cMaxExcludeChroms,"high priority - regular expressions defining chromosomes to exclude");
struct arg_str  *IncludeChroms = arg_strn("z","chromeinclude","<string>",0,cMaxIncludeChroms,"low priority - regular expressions defining chromosomes to include");
struct arg_end *end = arg_end(20);

void *argtable[] = {help,version,FileLogLevel,ScreenLogLevel,LogFile,
					pmode,strand,minlength,trunclength,ofsloci,deltalen,bkgndgroove,dyadratio,dyad2ratio,dyad3ratio,ExcludeChroms,IncludeChroms,
					inlocifile,outfile,inGenomefile,inConfFile,
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
		printf("\n      To invoke this parameter file then precede it's name with '@'");
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

	PMode = (etPMode)(pmode->count ? pmode->ival[0] : ePMDefault);
	if(PMode < ePMDefault || PMode >= ePMplaceholder)
		{
		printf("\nError: Processing mode '-m%d' specified outside of range %d..%d",PMode,ePMDefault,(int)ePMplaceholder-1);
		exit(1);
		}

	Strand = strand->count ? strand->sval[0][0] : '*';
	if(Strand != '+' && Strand != '-' && Strand != '*')
		{
		printf("\nError: Strand '-s%c' must be specified as '-', '+', or '*'",Strand);
		exit(1);
		}


	OfsLoci = ofsloci->count ? ofsloci->ival[0] : 0;
	if(abs(OfsLoci) > 2048)
		{
		printf("Error: loci offset '-u%d' must be in the range -2048 to +2048",OfsLoci);
		exit(1);
		}

	DeltaLen = deltalen->count ? deltalen->ival[0] : 0;
	if(abs(DeltaLen) > 2048)
		{
		printf("Error: delta length '-U%d' must be in the range -2048 to +2048",DeltaLen);
		exit(1);
		}

	MinLength = minlength->count ? minlength->ival[0] : 30;
	if(MinLength < 10)
		{
		printf("Error: minimum length '-l%d' must be in range 10..%d",MinLength,cMaxProcSeqLen);
		exit(1);
		}

	TruncLength = trunclength->count ? trunclength->ival[0] : 300;
	if(TruncLength < MinLength || TruncLength > cMaxProcSeqLen)
		{
		printf("Error: truncation length '-T%d' must be in range %d..%d",TruncLength,MinLength,cMaxProcSeqLen);
		exit(1);
		}

	if(PMode == ePMFiltNucConf)
		{
		BkgndGroove = bkgndgroove->count ? bkgndgroove->dval[0] : cDfltBkgndGroove;
		if(BkgndGroove != 0.0f && (BkgndGroove < cMinBkgndGroove || BkgndGroove > cMaxBkgndGroove))
			{
			printf("\nError: Bacgound minor groove '-b%1.4f' must be either 0.0 or in range %1.4f to %1.4f",BkgndGroove,cMinBkgndGroove,cMaxBkgndGroove);
			exit(1);
			}

		Dyadratio = dyadratio->count ? dyadratio->dval[0] : cDfltDyadratio;
		if(Dyadratio < cMinDyadratio || Dyadratio > cMaxDyadratio)
			{
			printf("\nError: Centrall dyad threshold ratio '-d%1.4f' must be in range %1.4f to %1.4f",Dyadratio,cMinDyadratio,cMaxDyadratio);
			exit(1);
			}


		Dyad2ratio = dyad2ratio->count ? dyad2ratio->dval[0] : cDfltDyad2ratio;
		if(Dyad2ratio < cMinDyad2ratio || Dyadratio > cMaxDyad2ratio)
			{
			printf("\nError: Centrall dyad threshold ratio '-d%1.4f' must be in range %1.4f to %1.4f",Dyad2ratio,cMinDyad2ratio,cMaxDyad2ratio);
			exit(1);
			}

		Dyad3ratio = dyad3ratio->count ? dyad3ratio->dval[0] : cDfltDyad3ratio;
		if(Dyad3ratio < cMinDyad2ratio || Dyadratio > cMaxDyad2ratio)
			{
			printf("\nError: Centrall dyad threshold ratio '-d%1.4f' must be in range %1.4f to %1.4f",Dyad3ratio,cMinDyad3ratio,cMaxDyad3ratio);
			exit(1);
			}

		if(inConfFile->count)
			strcpy(szConfFile,inConfFile->filename[0]);
		else
			{
			printf("\nError: In mode 1 expected structural conformation file to be specified with '-I<file>'");
			exit(1);
			}
		if(inGenomefile->count)
			strcpy(szGenomefile,inGenomefile->filename[0]);
		else
			{
			printf("\nError: In mode 1 expected bioseq genome assembly file to be specified with '-g<file>'");
			exit(1);
			}
		}
	else
		{
		Dyadratio = cDfltDyadratio;
		Dyad2ratio = cDfltDyad2ratio;
		Dyad3ratio = cDfltDyad3ratio;
		BkgndGroove = cDfltBkgndGroove;
		szGenomefile[0] = '\0';
		szConfFile[0] = '\0';
		}

	strcpy(szLocifile,inlocifile->filename[0]);
	strcpy(szOutfile,outfile->filename[0]);

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
		case ePMDefault:
			pszProcMode = "filter on chromosome/region/strand only";
			break;
		case ePMFiltNucConf:
			pszProcMode = "filter on chromosome/region/strand as well as on putative nucleosome conformation";
			break;
		};
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Processing mode: '%s'",pszProcMode);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Input element loci from csv file: '%s'",szLocifile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Write filtered and translated loci to this csv file: '%s'",szOutfile);
	if(PMode == ePMFiltNucConf)
		{
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Use genome assembly sequences from this bioseq file: '%s'",szGenomefile);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Use structural conformations from this file: '%s'",szConfFile);
		}

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Strand: '%c'",Strand);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Filter out elements with length less than: %d",MinLength);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Offset element loci by: %d",OfsLoci);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Change element lengths by: %d",DeltaLen);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Truncate elements to be no longer than: %d",TruncLength);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Filter out elements with length less than: %d",MinLength);

	if(PMode == ePMFiltNucConf)
		{
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"background minor groove: %1.4f",BkgndGroove);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"dyad minor groove threshold: %1.4f",Dyadratio);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"dyad immediately flanking minor groove threshold: %1.4f",Dyad2ratio);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"dyad remainder flanking minor groove threshold: %1.4f",Dyad3ratio);
		}

	for(Idx = 0; Idx < NumIncludeChroms; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"reg expressions defining chroms to include: '%s'",pszIncludeChroms[Idx]);
	for(Idx = 0; Idx < NumExcludeChroms; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"reg expressions defining chroms to exclude: '%s'",pszExcludeChroms[Idx]); 

#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	gStopWatch.Start();
	Rslt = Process(PMode,			// processing mode
		Strand,						// only retain elements on this strand
		OfsLoci,					// offset element start loci by this many nt
		DeltaLen,					// change element length by this many nt
		MinLength,					// slough elements of less than this length
		TruncLength,				// truncate elements to be a maximum of this length
		BkgndGroove,				// background minor groove 
		Dyadratio,					// dyad grooves must be at least this ratio to background
		Dyad2ratio,					// immediately flanking grooves must be at least this ratio to background
		Dyad3ratio,					// remainder of flanking grooves must be at least this ration to background
		szGenomefile,				// bioseq genome file
		szConfFile,					// file containing conformation characteristics
		szLocifile,					// read element loci ffrom this file
		szOutfile,					// write filtered and translated element loci to this file
		NumIncludeChroms,			// number of chromosome regular expressions to include
		pszIncludeChroms,			// array of include chromosome regular expressions
		NumExcludeChroms,			// number of chromosome expressions to exclude
		pszExcludeChroms);			// array of exclude chromosome regular expressions

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

CBioSeqFile *m_pBioSeqFile;			// holds instantiated genome assembly sequences
CTwister *m_pTwister;				// contains hunter group dsDNA conformational characteristic values
double m_BkgndGroove;				// dyads are relative to this background groove
double m_DyadratioThres;			// dyad grooves must be at least this ratio to background
double m_Dyad2ratioThres;			// immediately flanking grooves must be at least this ratio to background
double m_Dyad3ratioThres;			// remainder of flanking grooves must be at least this ration to background


etPMode m_PMode;					// processing mode

int m_hOutFile;						// file handle for output results
char *m_pszOutBuff;					// used for buffering output results
int m_AllocdOutBuff;				// how many chars have been allocated to m_pszOutBuff
int m_UsedOutBuff;					// how many chars of results currently in m_pszOutBuff

CCSVFile *m_pCSV;					// used whilst loading elements from CSV file
const int cElementAlloc = 1000000;	// allocate for this many tsElement
tsElement *m_pElements;				// pts to loaded elements of interest
int m_AllocdElements;				// memory allocated for this many elements
int m_NumElements;					// this many elements have been loaded

int m_ChromSeqLen;					// currently being processed assembly chromosome sequence length
char m_szCurChrom[cMaxDatasetSpeciesChrom+1]; // chromosome currently being processed

CIncExclChroms m_IncExclChroms;		// include/exclude chromosome processing class instance

int m_NumElTypes;					// current number of element types parsed
char m_szElType[cMaxElTypes][cMaxDatasetSpeciesChrom];
int m_NumSpecies;					// current numer of element species parsed
char m_szSpecies[cMaxSpecies][cMaxDatasetSpeciesChrom];


void
Init(void)
{
m_pCSV = NULL;
m_pBioSeqFile = NULL;
m_pTwister = NULL;
m_pszOutBuff = NULL;
m_pElements = NULL;
m_hOutFile = -1;
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
if(m_pszOutBuff != NULL)
	{
	delete m_pszOutBuff;
	m_pszOutBuff = NULL;
	}
if(m_pBioSeqFile != NULL)
	{
	delete m_pBioSeqFile;
	m_pBioSeqFile = NULL;
	}
if(m_pTwister != NULL)
	{
	delete m_pTwister;
	m_pTwister = NULL;
	}
if(m_pElements != NULL)
	{
	free(m_pElements);				// NOTE - use free() as memory malloc'd/realloc'd
	m_pElements = NULL;
	}
m_IncExclChroms.Reset();
m_AllocdOutBuff = 0;				// how many chars have been allocated to m_pszOutBuff
m_UsedOutBuff = 0;					// how many chars of results currently in m_pszOutBuff
m_AllocdElements = 0;				// memory allocated for this many regions
m_NumElements = 0;					// this many regions have been loaded
m_NumSpecies = 0;
m_NumElTypes = 0;
}


int
AddElement(int SrcID,				// region source identifier
		   int SpeciesID,			// source species
		   int TypeID,				// source type
		  int Instances,			// number of instances to associate with this region
		  int ChromID,				// is on this chrom
		  int StartLoci,			// starting at this loci
		  int EndLoci,				// and ending at this loci
		  bool bOnMinStrand)		// true if on minus strand
{
tsElement *pElement;

if(m_pElements == NULL || m_NumElements == m_AllocdElements)
	{
	if(m_pElements == NULL)
		{
		if((pElement = (tsElement *)malloc(sizeof(tsElement) * cElementAlloc))==NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to malloc memory (%d bytes requested) for region filtering",sizeof(tsElement) * cElementAlloc);
			Reset();
			return(eBSFerrMem);
			}
		m_AllocdElements = 0;
		m_NumElements = 0;
		}
	else
		{
		if((pElement = (tsElement *)realloc(m_pElements,sizeof(tsElement) * (m_AllocdElements+cElementAlloc)))==NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to realloc memory (%d bytes requested) for region filtering",sizeof(tsElement) * (m_AllocdElements+cElementAlloc));
			Reset();
			return(eBSFerrMem);
			}
		}
	m_pElements = pElement;
	m_AllocdElements += cElementAlloc;
	}
pElement = &m_pElements[m_NumElements++];
memset(pElement,0,sizeof(tsElement));
pElement->ElementID = m_NumElements;
pElement->SpeciesID = SpeciesID;
pElement->TypeID = TypeID;
pElement->SrcID = SrcID;
pElement->NumInstances = Instances;
pElement->ChromID = ChromID;
pElement->StartLoci = StartLoci;
pElement->EndLoci = EndLoci;
pElement->Strand = bOnMinStrand ? '-' : '+';
return(m_NumElements);
}


int			// actual length, or if under MinLen then -1, or if start would be < 0 or end < start then returns -2
AdjLoci(bool bOnMinStrand,		// true if element is on '-' strand
		int *pStartLoci,		// starting loci on '+' strand
		int *pEndLoci,			// ending loci on '+' strand
		int OfsLoci,			// offset loci by
		int DeltaLen,			// adjust length by
		int MinLen,				// must be at least this length
		int TruncLen)			// truncate to this length
{
int Ofs;
int DLen;

if(*pStartLoci < 0 || *pStartLoci > *pEndLoci)
	return(-2);

if(OfsLoci == 0 && DeltaLen == 0 && TruncLen == 0)
	return(1 + *pEndLoci -  *pStartLoci);

if(OfsLoci != 0 || DeltaLen != 0)
	{
	if(bOnMinStrand)
		{
		Ofs = OfsLoci * -1;
		DLen = DeltaLen * -1;
		}
	else
		{
		Ofs = OfsLoci;
		DLen = DeltaLen;
		}

	*pStartLoci += Ofs;
	*pEndLoci += Ofs;
	if(bOnMinStrand)
		*pStartLoci += DLen;
	else
		*pEndLoci += DLen;
	}

if(*pStartLoci > *pEndLoci || (!bOnMinStrand && *pStartLoci < 0))
	return(-2);
DLen = 1 + *pEndLoci - *pStartLoci;
if(DLen < MinLen)
	return(-1);
if(DLen > TruncLen)
	{
	if(bOnMinStrand)
		{
		*pStartLoci = 1 + *pEndLoci - TruncLen;
		if(*pStartLoci < 0)
			return(-2);
		}
	else
		*pEndLoci = *pStartLoci + TruncLen - 1;
	DLen= TruncLen;
	}
return(DLen);
}

int
AddSpecies(char *pszSpecies)
{
int Idx;
if(m_NumSpecies > 0)
	{
	for(Idx = 0; Idx <  m_NumSpecies; Idx++)
		if(!stricmp(m_szSpecies[Idx],pszSpecies))
			return(Idx+1);
	}
if(m_NumSpecies == cMaxSpecies)
	return(-1);
strcpy(m_szSpecies[m_NumSpecies++],pszSpecies);
return(m_NumSpecies);
}

int
AddElType(char *pszElType)
{
int Idx;
if(m_NumElTypes > 0)
	{
	for(Idx = 0; Idx <  m_NumElTypes; Idx++)
		if(!stricmp(m_szElType[Idx],pszElType))
			return(Idx+1);
	}
if(m_NumElTypes == cMaxElTypes)
	return(-1);
strcpy(m_szElType[m_NumElTypes++],pszElType);
return(m_NumSpecies);
}

int
LoadElements(char *pszElementFile,		// csv file containing region loci
		int OfsLoci,					// offset region start loci by this many nt
		int DeltaLen,					// change region length by this many nt
		int TruncLength,				// truncate regions to be a maximum of this length
		int MinLength)					// regions to be of this minimum length	
{
int Rslt;
int NumFields;
int NumProcessed;

int Len;
int SrcID;
char *pszElType;
char *pszRefSpecies;
int ChromID;
char *pszChrom;
char szPrevChrom[cMaxDatasetSpeciesChrom+1];
int ChromSeqLen;
char *pszStrand;
int StartLoci;
int EndLoci;
int Instances;
bool bOnMinStrand;
int NumUnderlen;
int NumStartLess0;
int NumEndPastChrom;
int NumFiltChrom;

int SpeciesID;
int TypeID;

if(pszElementFile == NULL || pszElementFile[0] == '\0')
	return(eBSFSuccess);

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Loading element loci from CSV file '%s'",pszElementFile);

m_pCSV = new CCSVFile;
if(m_pCSV == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CCSVfile");
	return(eBSFerrObj);
	}

if((Rslt=m_pCSV->Open(pszElementFile))!=eBSFSuccess)
	{
	while(m_pCSV->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,m_pCSV->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open file: %s",pszElementFile);
	Reset();
	return(Rslt);
	}

NumUnderlen = 0;
NumStartLess0 = 0;
NumEndPastChrom = 0;
NumFiltChrom = 0;
NumProcessed = 0;
ChromID = 0;
szPrevChrom[0] = '\0';
ChromSeqLen = 0;
while((Rslt=m_pCSV->NextLine()) > 0)	// onto next line containing fields
	{
	NumFields = m_pCSV->GetCurFields();
	if(NumFields < 7)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expected 7+ fields in '%s', GetCurFields() returned '%d'",pszElementFile,NumFields);
		Rslt = eBSFerrFieldCnt;
		break;
		}
	if(!NumProcessed && m_pCSV->IsLikelyHeaderLine())
		continue;
	NumProcessed += 1;

	m_pCSV->GetText(4,&pszChrom);
	if((Rslt = m_IncExclChroms.IncludeThisChrom(pszChrom))<=0)
		{
		if(Rslt < 0)
			break;
		NumFiltChrom += 1;
		continue;
		}

	if(m_pBioSeqFile != NULL)
		{
		if(ChromID == 0 || stricmp(pszChrom,szPrevChrom))
			{
			if((ChromID = m_pBioSeqFile->LocateEntryIDbyName(pszChrom)) < 1)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to locate element chromosome '%s'",pszChrom);
				return(eBSFerrChrom);
				}
			ChromSeqLen = m_pBioSeqFile->GetDataLen(ChromID);
			strcpy(szPrevChrom,pszChrom);
			}
		}
	else
		{
		ChromID = Rslt;
		ChromSeqLen = -1;
		}
	m_pCSV->GetInt(7,&Len);
	m_pCSV->GetInt(1,&SrcID);
	m_pCSV->GetText(2,&pszElType);
	m_pCSV->GetText(3,&pszRefSpecies);
	m_pCSV->GetInt(5,&StartLoci);
	m_pCSV->GetInt(6,&EndLoci);

	bOnMinStrand = false;
	if(NumFields >= 8)					// check if strand has been specified
		{
		m_pCSV->GetText(8,&pszStrand);
		if(pszStrand[0] == '-')		// assume anything other than '-' is on the plus strand
			bOnMinStrand = true;
		}
	Instances = 1;
	if(NumFields >= 11)					// check if likely short reads loci with associated reads count
		m_pCSV->GetInt(11,&Instances);

	if((Rslt = AdjLoci(bOnMinStrand,		// true if element is on '-' strand
		&StartLoci,		// starting loci on '+' strand
		&EndLoci,			// ending loci on '+' strand
		OfsLoci,			// offset loci by
		DeltaLen,			// adjust length by
		MinLength,			// must be at least this length
		TruncLength)) < 1)	// truncate to this length
		{
		if(Rslt == -1)
			NumUnderlen += 1;
		else
			if(Rslt == -2)
				NumStartLess0 += 1;
		continue;
		}

	if(ChromSeqLen > 0 && EndLoci >= ChromSeqLen)
		{
		NumEndPastChrom += 1;
		continue;
		}
	if((SpeciesID = AddSpecies(pszRefSpecies))<1)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to add species '%s', can process at most %d",pszRefSpecies,cMaxSpecies);
		Reset();
		return(eBSFerrMaxEntries);
		}

	if((TypeID = AddElType(pszElType))<1)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to add element type '%s', can process at most %d",pszElType,cMaxElTypes);
		Reset();
		return(eBSFerrMaxEntries);
		}

	if((Rslt=AddElement(SrcID,SpeciesID,TypeID,Instances,ChromID,StartLoci,EndLoci,bOnMinStrand))<0)
		{
		Reset();
		return(Rslt);
		}
	}
delete m_pCSV;
m_pCSV = NULL;
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Accepted %d elements from %d loaded - %d chrom filtered, %d too short, %d adjusted loci before chrom starts, %d end loci past end of chrom",m_NumElements,NumProcessed,NumFiltChrom,NumUnderlen,NumStartLess0,NumEndPastChrom);
if(m_NumElements > 1)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Sorting elements...");
	qsort(m_pElements,m_NumElements,sizeof(tsElement),SortElements);
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Finished sorting elements...");
	}
return(m_NumElements);
}

int 
LocateDyads(int SeqLen,				// sequence is this length
			etSeqBase *pSeq,		// interested in dyads for this sequence
			int *pScore,			// score for higest scoring dyad
			int *pDyadIdx,			// index into pSeq[] of highest scoring dyad
			int *pNumDyads)			// number of dyads located
{
int Rslt;
int *pTwist;

double DyadRatio;
double Dyad2Ratio;
double Dyad3Ratio;

int NumSteps;
int *pValues;
int *pConfGroove;
int *pConfTwist;
int ChkGroove[13];		// to hold dyad (ChkGroove[6]) and +/- 6 at approx decimer (depends on twist) offsets 
int DecIdx;				// index into ChkGroove, incr/decr every 360 degree twist
int GrooveCnt;			// to hold number of groove values contributing to current ChkGroove[DecIdx] so average can be calculated
int AccumTwist;			// to hold accumulated twist relative to dyad
int ChkTwist;			// AccumTwist % 360 used to determine if minor groove back on same plane as at dyad

int ConfValues[cMaxProcSeqLen];
int TwistValues[cMaxProcSeqLen];	

int SeqIdx;
int DyadFirstOfs;
int	DyadLastOfs;
int NumDyads;
int DyadIdx;
int Score;
int PrevBestScore;

if(SeqLen > cMaxProcSeqLen || pSeq == NULL)
	return(eBSFerrParams);

if(pScore != NULL)
	*pScore = 0;
if(pDyadIdx != NULL)
	*pDyadIdx = -1;
if(pNumDyads != NULL)
	*pNumDyads = 0;

if((Rslt = m_pTwister->GetSequenceConformation(eSSminorgroove,	// process for this conformational parameter
					  0,								// initial starting offset (0..n) in pSeq
					  0,				                // number of steps (0 for all) to process starting at pSeq[iStartPsn]|pSeq[iStartPsn+1]
  					  SeqLen,				// number of nucleotides
					  pSeq,				// sequence to be processed
					  ConfValues))!=eBSFSuccess)		// where to return conformational values
					{
					gDiagnostics.DiagOut(eDLFatal,"LoadSequences","ProcessSequence failed");
					Reset();
					return(Rslt);
					}

if((Rslt = m_pTwister->GetSequenceConformation(eSStwist,	// process for this conformational parameter
						  0,									// initial starting offset (0..n) in pSeq
						  0,					                // number of steps (0 for all) to process starting at pSeq[iStartPsn]|pSeq[iStartPsn+1]
  						  SeqLen,								// number of nucleotides
						  pSeq,									// sequence to be processed
						  TwistValues))!=eBSFSuccess)			// where to return conformational values
						{
						gDiagnostics.DiagOut(eDLFatal,"LoadSequences","ProcessSequence failed");
						Reset();
						return(Rslt);
						}

NumSteps = SeqLen - 1;
pValues = ConfValues;
pTwist = TwistValues;
DyadIdx = -1;
Score = 0;
PrevBestScore = 0;
NumDyads = 0;
DyadFirstOfs = 68;				// allows for twist as 6 360 degrees of twist can require 68 nt
DyadLastOfs = NumSteps - 68;
if(DyadLastOfs < DyadFirstOfs)
	return(0);

memset(ChkGroove,0,sizeof(ChkGroove));
for(SeqIdx = DyadFirstOfs; SeqIdx <= DyadLastOfs; SeqIdx++)
	{
	DecIdx = 6;
	pConfGroove = &pValues[SeqIdx];
	ChkGroove[DecIdx++] = *pConfGroove++;
	DyadRatio = (double)ChkGroove[6]/(m_BkgndGroove * 10000.0f);
	if(DyadRatio < m_DyadratioThres)
		continue;
	pConfTwist =  &pTwist[SeqIdx+1];
	AccumTwist = *pConfTwist;
	ChkGroove[DecIdx] = 0;
	GrooveCnt = 0;
	
	// iterate over bases to right of putative dyad and every rotation of the dsDNA get the minor groove
	int Bases = 1;
	while(DecIdx <= 12)
		{
		Bases += 1;
		pConfTwist += 1;
		pConfGroove += 1;
		AccumTwist += *pConfTwist;
		ChkTwist = AccumTwist % 3600000;
		if(ChkTwist >= 3300000 || ChkTwist <= 300000)
			{
			ChkGroove[DecIdx] += *pConfGroove;
			GrooveCnt += 1;
			}
		else
			{
			if(GrooveCnt > 0)
				{
				ChkGroove[DecIdx] /= GrooveCnt;
				GrooveCnt = 0;
				if(DecIdx++ < 12)
					ChkGroove[DecIdx]= 0;
				}
			}
		}
	// now iterate over bases to left of putative dyad and every rotation of the dsDNA get the minor groove
	DecIdx = 5;
	pConfGroove =  &pValues[SeqIdx-1];
	pConfTwist =  &pTwist[SeqIdx-1];
	AccumTwist = *pConfTwist;
	ChkGroove[DecIdx] = 0;
	GrooveCnt = 0;
	Bases = 1;
	while(DecIdx >= 0)
		{
		Bases += 1;
		pConfTwist -= 1;
		pConfGroove -= 1;
		AccumTwist += *pConfTwist;
		ChkTwist = AccumTwist % 3600000;
		if(ChkTwist >= 3300000 || ChkTwist <= 300000)
			{
			ChkGroove[DecIdx] += *pConfGroove;
			GrooveCnt += 1;
			}
		else
			{
			if(GrooveCnt > 0)
				{
				ChkGroove[DecIdx] /= GrooveCnt;
				GrooveCnt = 0;
				if(DecIdx-- > 0)
					ChkGroove[DecIdx] = 0;
				}
			}
		}

	Dyad2Ratio = (double)(ChkGroove[5] + ChkGroove[7])/(2*m_BkgndGroove*10000.0f);
	Dyad3Ratio = (double)(ChkGroove[0] + ChkGroove[1] + ChkGroove[2] + ChkGroove[3] + ChkGroove[4] +
				ChkGroove[8] + ChkGroove[9] + ChkGroove[10] + ChkGroove[11] + ChkGroove[12])/(10*m_BkgndGroove*10000.0f);

	
	if(Dyad2Ratio < m_Dyad2ratioThres || Dyad3Ratio < m_Dyad3ratioThres)
		continue;

		// how to really score these dyads? One day real scores will be associated...
	int Score = (int)(1000 * ((DyadRatio - 1.0f) + ((Dyad2Ratio - 1.0f) * 0.85) + ((Dyad3Ratio - 1.0f) * 0.75)));

	if(Score > PrevBestScore)
		{
		if(pScore != NULL)
			*pScore = Score;
		if(pDyadIdx != NULL)
			*pDyadIdx = SeqIdx;
		PrevBestScore = Score;
		}
	NumDyads += 1;
	if(pNumDyads != NULL)
		*pNumDyads = NumDyads;
	}
return(NumDyads);
}


int
Process(etPMode PMode,						// processing mode
		char Strand,						// only retain elements on this strand
		int OfsLoci,						// offset element start loci by this many nt
		int DeltaLen,						// change element length by this many nt
		int MinLength,						// slough elements of less than this length
		int TruncLength,					// truncate elements to be a maximum of this length
		double BkgndGroove,					// background minor groove 
		double DyadratioThres,				// dyad grooves must be at least this ratio to background
		double Dyad2ratioThres,				// immediately flanking grooves must be at least this ratio to background
		double Dyad3ratioThres,				// remainder of flanking grooves must be at least this ration to background
		char *pszInGenomefile,				// bioseq genome file
		char *pszInConfFile,				// file containing conformation characteristics
		char *pszLociFile,					// read element loci ffrom this file
		char *pszRsltsFile,					// write filtered and translated element loci to this file
		int	NumIncludeChroms,				// number of chromosome regular expressions to include
		char **ppszIncludeChroms,			// array of include chromosome regular expressions
		int	NumExcludeChroms,				// number of chromosome expressions to exclude
		char **ppszExcludeChroms)			// array of exclude chromosome regular expressions
{
int Rslt;
int Score;					// score for higest scoring dyad
int DyadIdx;				// index into pSeq[] of highest scoring dyad
int NumDyads;				// number of dyads located
tsElement *pElement;
int ElIdx;
int TotNumDyads;
int NumFilterEnds;
int CurChromID;
int	ChromTotDyads;
int	ChromSeqs;


etSeqBase Sequence[cMaxProcSeqLen];		// element sequence currently being processed for nucleosome dyads
int SeqLen;					// length of element sequence currently being processed

Init();

if((m_pszOutBuff = new char [cAllocOutBuff])==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to alloc memory (%d bytes requested) for bufering output results",cAllocOutBuff);
	Reset();
	return(eBSFerrMem);
	}
m_AllocdOutBuff = cAllocOutBuff;

#ifdef _WIN32
if((m_hOutFile = open(pszRsltsFile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
if((m_hOutFile = open(pszRsltsFile, O_RDWR | O_CREAT | O_TRUNC, S_IREAD | S_IWRITE ))==-1)
#endif
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create output file: %s - %s",pszRsltsFile,strerror(errno));
	Reset();
	return(eBSFerrCreateFile);
	}

if((Rslt=m_IncExclChroms.InitChromExclusion(NumIncludeChroms,ppszIncludeChroms,NumExcludeChroms,ppszExcludeChroms))!=eBSFSuccess)
	{
	Reset();
	return(Rslt);
	}

// if specified, then load the bioseq genome assembly file first as that is needed for the chromosome identifier mapping
if(pszInGenomefile[0] != '\0')
	{
	if((m_pBioSeqFile = new CBioSeqFile) == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CBioSeqFile object");
		Reset();
		return(eBSFerrObj);
		}

	if((Rslt = m_pBioSeqFile->Open(pszInGenomefile))!=eBSFSuccess)
		{
		while(m_pBioSeqFile->NumErrMsgs())
			gDiagnostics.DiagOut(eDLFatal,gszProcName,m_pBioSeqFile->GetErrMsg());
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open genome assembly sequence file '%s'",pszInGenomefile);
		Reset();
		return(Rslt);
		}
	}

// can now load the input elements
if((Rslt = LoadElements(pszLociFile,// csv file containing element loci
		OfsLoci,					// offset element start loci by this many nt
		DeltaLen,					// change element length by this many nt
		TruncLength,				// truncate element to be a maximum of this length
		MinLength)) <= 0)			// element to be of this minimum length	
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to load elements loci file '%s'",pszLociFile);
		Reset();
		return(Rslt);
		}

if(PMode == ePMFiltNucConf)
	{
	// can now apply any conformational characteristics filtering
	if((m_pTwister = new CTwister)==NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,"Process","Unable to create CTwister object");
		Reset();
		return(eBSFerrObj);
		}
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Loading structural characteristic values from %s...",pszInConfFile);
	if((Rslt = m_pTwister->LoadStructParams(pszInConfFile))  < eBSFSuccess)
		{
		while(m_pTwister->NumErrMsgs())
			gDiagnostics.DiagOut(eDLFatal,gszProcName,m_pTwister->GetErrMsg());
		gDiagnostics.DiagOut(eDLFatal,"Process","LoadStructParams(%s) failed",pszInConfFile);
		Reset();
		return(Rslt);
		}

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Now filtering for dyads...",pszInConfFile);

	m_BkgndGroove=BkgndGroove;			// dyads are relative to this background groove
	m_DyadratioThres=DyadratioThres;	// dyad grooves must be at least this ratio to background
	m_Dyad2ratioThres=Dyad2ratioThres;	// immediately flanking grooves must be at least this ratio to background
	m_Dyad3ratioThres=Dyad3ratioThres;	// remainder of flanking grooves must be at least this ration to background

	// iterate over each element, get it's sequence and check if it contains a putative dyad..
	NumFilterEnds = 0;
	TotNumDyads = 0;
	ChromTotDyads = 0;
	ChromSeqs = 0;
	CurChromID = 0;
	pElement = m_pElements;
	for(ElIdx = 0; ElIdx < m_NumElements; ElIdx++,pElement++)
		{
		if(pElement->ChromID != CurChromID)
			{
			if(CurChromID > 0)
				gDiagnostics.DiagOut(eDLInfo,gszProcName,"On chrom %s from %d sequences there are %d putative dyads",m_szCurChrom,ChromSeqs,ChromTotDyads);
			m_pBioSeqFile->GetName(pElement->ChromID,sizeof(m_szCurChrom),m_szCurChrom);
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing %s...",m_szCurChrom);
			m_ChromSeqLen = m_pBioSeqFile->GetDataLen(pElement->ChromID);
			CurChromID = pElement->ChromID;
			ChromTotDyads = 0;
			ChromSeqs = 0;
			}
		ChromSeqs += 1;
		if(pElement->EndLoci >= m_ChromSeqLen)
			{
			if(NumFilterEnds++ < 10)
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Sequence starting at loci %s %d and finishing %d extends past chromosome length of %d - sloughing'",m_szCurChrom,pElement->StartLoci,SeqLen);
			continue;
			}
		SeqLen = 1 + pElement->EndLoci - pElement->StartLoci;
		if(SeqLen > cMaxProcSeqLen)	// should have been picked up earlier - but better safe than core dump :-)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Sequence of length %d is too long",SeqLen);
			Reset();
			return(eBSFerrOfs);
			}
		if((Rslt=m_pBioSeqFile->GetData(pElement->ChromID,eSeqBaseType,pElement->StartLoci,Sequence,SeqLen)) != SeqLen)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Loading sequence of length %d failed from chrom: %s file: '%s'",SeqLen,m_szCurChrom,pszInGenomefile);
			Reset();
			return(Rslt);
			}

	// remove any repeat masking and randomly substitute bases for eBaseN's - not expecting too many of these say's he hopefully!
	int SeqIdx;
	etSeqBase *pSeq = Sequence;
	for(SeqIdx = 0; SeqIdx < SeqLen; SeqIdx++,pSeq++)
		if((*pSeq &= ~cRptMskFlg) > eBaseT)
			*pSeq = rand() % 4;

	Rslt = LocateDyads(SeqLen,		// sequence is of this length
				Sequence,				// interested in dyads for this sequence
				&Score,					// score for higest scoring dyad
				&DyadIdx,				// index into pSeq[] of highest scoring dyad
				&NumDyads);				// number of dyads located
	if(Rslt > 0)
			{
			TotNumDyads += 1;
			ChromTotDyads += 1;
			pElement->Score = Score;
			pElement->DyadIdx = DyadIdx;
			pElement->NumDyads = NumDyads;
			}
		}
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"On chrom %s from %d sequences there are %d putative dyads",m_szCurChrom,ChromSeqs,ChromTotDyads);
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Over all chroms, %d filtered sequences, there are %d putative dyads",m_NumElements,TotNumDyads);
	}

// can now write out those filtered and transformed element loci which are still of interest
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Writing filtered %s to file '%s'",pszRsltsFile,PMode == ePMFiltNucConf ? "loci with putative dyads " : "loci");
pElement = m_pElements;
CurChromID = 0;
for(ElIdx = 0; ElIdx < m_NumElements; ElIdx++,pElement++)
	{
	if(pElement->ChromID != CurChromID)
		{
		if(m_pBioSeqFile != NULL)
			{
			m_pBioSeqFile->GetName(pElement->ChromID,sizeof(m_szCurChrom),m_szCurChrom);
			m_ChromSeqLen = m_pBioSeqFile->GetDataLen(pElement->ChromID);
			CurChromID = pElement->ChromID;
			}
		else
			strcpy(m_szCurChrom,m_IncExclChroms.LocateChrom(pElement->ChromID));
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Writing out filtered elements %s...",m_szCurChrom);
		}
	if(PMode != ePMFiltNucConf || pElement->Score > 0)
		{
		m_UsedOutBuff += sprintf(&m_pszOutBuff[m_UsedOutBuff],"%d,\"%s\",\"%s\",\"%s\",%d,%d,%d,\"%c\",%d,%d,%d,%d\n",
			pElement->SrcID,m_szElType[pElement->TypeID-1],m_szSpecies[pElement->SpeciesID-1],m_szCurChrom,
			pElement->StartLoci,pElement->EndLoci,1 + pElement->EndLoci-pElement->StartLoci,pElement->Strand,
			PMode,pElement->NumDyads,pElement->NumInstances,pElement->Score);

		if((m_UsedOutBuff+1000) > m_AllocdOutBuff)
			{
			if((Rslt=write(m_hOutFile,m_pszOutBuff,m_UsedOutBuff))!=m_UsedOutBuff)
				{
				gDiagnostics.DiagOut(eDLFatal,"Unable to write file %s - error %s",pszRsltsFile,strerror(errno));
				Reset();
				return(eBSFerrWrite);
				}
			m_UsedOutBuff = 0;
			}
		}
	}
if(m_UsedOutBuff > 0)
	if((Rslt=write(m_hOutFile,m_pszOutBuff,m_UsedOutBuff))!=m_UsedOutBuff)
		{
		gDiagnostics.DiagOut(eDLFatal,"Unable to write file %s - error %s",pszRsltsFile,strerror(errno));
		Reset();
		return(eBSFerrWrite);
		}
Reset();
return(0);
}

// SortElements
// sort elements by ChromID-->StartLoci-->EndLoci ascending
int SortElements( const void *arg1, const void *arg2)
{
tsElement *pR1 = (tsElement *)arg1;
tsElement *pR2 = (tsElement *)arg2;
if(pR1->ChromID < pR2->ChromID)
	return(-1);
if(pR1->ChromID > pR2->ChromID)
	return(1);
if(pR1->StartLoci < pR2->StartLoci)
	return(-1);
if(pR1->StartLoci > pR2->StartLoci)
	return(1);
return(pR1->EndLoci - pR2->EndLoci);
}