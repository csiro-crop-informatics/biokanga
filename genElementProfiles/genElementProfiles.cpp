// genElementProfiles.cpp : Defines the entry point for the console application
// generates element frequency profiles ( could be short read alignments ) around transcriptional start sites
// or other user specified annotated regions
// 1.1.2 added option to bin by read starts
// 1.2.0 extended binning option to also allow for unique loci read starts
// 1.2.0 started to add processing for SAM alignment files
// 1.2.2 completed adding processing for SAM alignment files
// 1.2.3 input files can be wildcarded 

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

const char *cpszProgVer = "1.2.4";			// increment with each release

const int cChromSeqReAlloc = 5000000;	// realloc chrom sequence size
const UINT32 cMaxAccumCnt= 0x07ffffff;  // truncate accumulated counts to be at most this value
const UINT32 cRetainCnt  = 0x08000000;  // accumulated counts >= this value are counts marked as being to be retained

const int cMinFlankLen   = 0;			// minimum user specified flank length
const int cDfltFlankLen  = 2000;		// default flank length
const int cMaxFlankLen   = 100000;		// maximum user specified flank length

const int cMinNumBins    = 10;			// minimum user specified number of bins in each flank
const int cDfltNumBins   = 100;			// default number of bins
const int cMaxNumBins    = 2000;		// maximum user specified number of bins

const int cMinKDist      = 0;			// minimum user specified distance of overlapping other features
const int cDfltKDist     = 1000;		// default distance of overlapping other features
const int cMaxKDist      = 100000;		// maximum user specified distance of overlapping other features

const int cMaxInFileSpecs = 20;			// allow for upto this many input wildcarded files

// processing modes
typedef enum TAG_ePMode {
	ePMauto,					// autodetermination of file type
	ePMCsv,						// CSV format
	ePMBed,						// UCSC BED format
	ePMSam,						// SAM format
	ePMplaceholder				// used to set the enumeration range
	} etPMode;

// output format modes
typedef enum TAG_eFMode {
	eFMdefault,					// default is for CSV only
	eFMplaceholder				// used to set the enumeration range
	} etFMode;

// profile relative to
typedef enum ePFLFeature {
	ePFLIsoforms,				// profile gene isoforms
	ePFLTSS,					// profile TSS flanking regions
	ePFLTES,					// profile TES flanking regions
	ePFLplaceholder				// used to set the enumeration range
} etPFLFeature;

// profile overlap policy
typedef enum ePOPpolicy {
	ePOPdefault,				// accept all reads 
	ePOPkdist,					// do not accept reads within -K distance of overlapping other features
	ePOPplaceholder				// used to set the enumeration range
} etPOPpolicy;

typedef enum eReadProfile {
	eRPdensity = 0,				// profile for: 0 read density
	eRPstartcnts,				// profile read start counts 
	eRPstartuniques,			// unique read starts
	eRPplaceholder
	} etReadProfile;

int
Process(etPMode PMode,					// processing mode
		etFMode FMode,					// output format - CSV or BED
		etReadProfile ReadProfile,		// profile for: 0 read density, 1: read starts, 2:  unique read starts (default is for read density)
		int PkWinLen,					// window length over which to calculate the peak counts
		char Strand,					// ROI are on this strand
		etPFLFeature Feature,			// center profiles around this functional feature
		int InterGenicLen,				// intergenic flank length
		int IntraGenicLen,				// intragenic flank length
		int NumBins,					// profile into this many bins intragenic and intergenic flank counts
		etPOPpolicy ProfPol,			// profile overlap policy
		int KDist,						// overlap distance
		char *pszExcludeFile,			// reads exclude features
		char *pszFeatureFile,			// file containing annotation regions of interest to user
		int NumInputFiles,				// number of input files
		char **ppszInFileSpecs,			// input files
		char *pszRsltsFile);			// output CSV profilefile

char *Feature2Txt(etPFLFeature Feature);

CStopWatch gStopWatch;
CDiagnostics gDiagnostics;				// for writing diagnostics messages to log file
char gszProcName[_MAX_FNAME];			// process name

int m_hRsltsFile;						// handle for opened results file

CCSVFile *m_pCSVFile;							// used if processing input CSV file for coverage
CBEDfile *m_pBEDFile;							// used if processing input BED files for coverage
CBEDfile *m_pExclBEDFile;						// used if exclusion BED file was specified by user

const int cMaxChromCov = 1000;					// can handle at most this many chromosomes
const int cAllocCovCnts = 0x07fffff;			// allocate for chrom coverage cnts in this sized increments

typedef struct TAG_sChromCnts {
	char szChrom[cMaxDatasetSpeciesChrom+1];	// coverage is on this chromosome
	int AllocCovCnts;							// allocated to hold cnts for this sized chromosome
	int StartOfs;								// pCovCnts[offset] of first coverage cnt
	int EndOfs;									// pCovCnts[offset] of last coverage cnt
	UINT32 *pCovCnts;							// coverage counts
} tsChromCnts;

tsChromCnts m_ChromCnts[cMaxChromCov];
int m_NumChromsCov = 0;

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
	return _T("genElementProfiles");
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
int Idx;

int Rslt = 0;   			// function result code >= 0 represents success, < 0 on failure

etPMode PMode;				// processing mode
etFMode FMode;				// format output mode

int Strand;						// ROI are on this strand
etPFLFeature Feature;			// center profiles around this functional feature
int InterGenicLen;				// intergenic flank length
int IntraGenicLen;				// intragenic flank length
int NumBins;					// profile into this many bins
int PkWinLen;					// peak window length
etPOPpolicy ProfPol;			// profile overlap policy
int KDist;						// do not accept reads within K distance of overlapping other features
etReadProfile ReadProfile;		// profile for: 0 read density, 1: read starts, 2:  unique read starts (default is for read density)

char szRsltsFile[_MAX_PATH];	// output results file
int NumInputFiles;			// number of input files
char *pszInFileSpecs[cMaxInFileSpecs];		// input files
char szFeatureFile[_MAX_PATH];	// annotated feature filter file
char szExcludeFile[_MAX_PATH];	// exclude element overlap file

// command line args
struct arg_lit  *help    = arg_lit0("hH","help",                "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"diagnostics log file");

struct arg_int *pmode = arg_int0("m","mode","<int>",		    "input reads loci file format: 0 - Auto, 1 - CSV, 2 - BED, 3 - SAM (default = 0)");
struct arg_int *format = arg_int0("M","format","<int>",		    "output format: 0 - CSV (currently only CSV is supported)");

struct arg_int *readprofile = arg_int0("P","readprofile","<int>","profile for: 0 read density, 1: read starts, 2:  unique read starts (default is for read density)");
struct arg_int *pkwinlen = arg_int0("w","peakwinlen","<int>",	"use this window size when determining transcript peak counts (read counts only, range 10 to 500, default 100bp)");

struct arg_int *profpol = arg_int0("p","profpol","<int>",		"profile overlap policy: 0 - accept all reads, 1 - do not accept reads overlapping other features, 2 - do not accept reads within -K distance of overlapping other features (default: 0)");
struct arg_int *kdist = arg_int0("k","kdist","<int>",			"do not accept reads within this distance of overlapping other features (default = 1000)");

struct arg_int  *strand = arg_int0("s","strand","<int>",		"profile for this strand: 0 - any, 1 - Watson '+', 2 - Crick '-' (default is any)");
struct arg_int *intergeniclen = arg_int0("l","intergeniclen","<int>", "intergenic flank length (default = 1000 for TSS and TSE)");
struct arg_int *intrageniclen = arg_int0("L","intrageniclen","<int>", "intragenic flank length (default = 1000 for TSS and TSE)");
struct arg_int *numbins = arg_int0("n","numbins","<int>",       "profile into this number of bins (default = 100)");
struct arg_int *feature = arg_int0("r","feature","<int>",		"profile for these features: 0 - gene isoforms, 1 - TSS, 2 - TES (default = isoforms)");

struct arg_file *infiles = arg_filen("i","in","<file>",0,cMaxInFileSpecs,"input from these CSV, BED, or SAM reads element loci files (wildcards allowed)");
struct arg_file *featurefile = arg_file1("I","features","<file>",	"annotated genomic feature definition BED file");
struct arg_file *excludefile = arg_file0("x","exclude","<file>","exclude reads from profile if overlaping elements in this BED file");
struct arg_file *outfile = arg_file1("o","out","<file>",		"output regions to this file");

struct arg_end *end = arg_end(20);

void *argtable[] = {help,version,FileLogLevel,LogFile,
					pmode,format,profpol,pkwinlen,kdist,readprofile,strand,feature,intergeniclen,intrageniclen,numbins,infiles,excludefile,featurefile,outfile,
					end};

char **pAllArgs;
int argerrors;
argerrors = CUtility::arg_parsefromfile(argc,(char **)argv,&pAllArgs);
if(argerrors >= 0)
	argerrors = arg_parse(argerrors,pAllArgs,argtable);

    /* special case: '--help' takes precedence over error reporting */
if (help->count > 0)
        {
		printf("\n%s Generate Profiles, Version %s\nOptions ---\n", gszProcName,cpszProgVer);
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

	PMode = (etPMode)(pmode->count ? pmode->ival[0] : 0);
	if(PMode < 0 || PMode >= ePMplaceholder)
		{
		printf("\nError: Input file format '-m%d' specified outside of range %d..%d",PMode,0,(int)ePMplaceholder-1);
		exit(1);
		}

	FMode = (etFMode)(format->count ? format->ival[0] : eFMdefault);
	if(FMode < eFMdefault || FMode >= eFMplaceholder)
		{
		printf("\nError: Output file format mode '-M%d' specified outside of range %d..%d",FMode,eFMdefault,(int)eFMplaceholder-1);
		exit(1);
		}

	PkWinLen = pkwinlen->count ? pkwinlen->ival[0] : 100;
	if(PkWinLen < 10 || PkWinLen > 500)
		{
		printf("\nError: Peak window length '-w%d' specified outside of range 10..500",PkWinLen);
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


	ReadProfile = (etReadProfile)(readprofile->count ? readprofile->ival[0] : eRPdensity);	
	if(ReadProfile < eRPdensity || ReadProfile > eRPplaceholder)
		{
		printf("\nSpecified read profiling type '-P%d' outside of range 0..%d",eRPdensity,eRPplaceholder-1);
		exit(1);
		}

	Feature = (etPFLFeature)(feature->count ? feature->ival[0] : ePFLIsoforms);	
	if(Feature < ePFLIsoforms || Feature > ePFLplaceholder)
		{
		printf("\nSpecified profile centric feature '-f%d' outside of range 0..%d",Feature,ePFLplaceholder-1);
		exit(1);
		}

	NumBins = numbins->count ? numbins->ival[0] : cDfltNumBins;
	if(NumBins < cMinNumBins || NumBins > cMaxNumBins)
		{
		printf("\nError: Number of bins '-n%d' specified outside of range %d..%d",NumBins,cMinNumBins,cMaxNumBins);
		exit(1);
		}
	if(Feature == ePFLIsoforms && NumBins > 200)
		{
		printf("\nError: Number of bins '-n%d' must be <= 200 when profiling isoforms",NumBins);
		exit(1);
		}

	if(Feature > ePFLIsoforms)
		{
		InterGenicLen= intergeniclen->count ? intergeniclen->ival[0] : cDfltFlankLen;
		if(InterGenicLen < cMinFlankLen || InterGenicLen > cMaxFlankLen)
			{
			printf("\nError: Intergenic flank length '-l%d' specified outside of range %d..%d",InterGenicLen,cMinFlankLen,cMaxFlankLen);
			exit(1);
			}

		IntraGenicLen= intrageniclen->count ? intrageniclen->ival[0] : cDfltFlankLen;
		if(IntraGenicLen < cMinFlankLen || IntraGenicLen > cMaxFlankLen)
			{
			printf("\nError: Intragenic flank length '-l%d' specified outside of range %d..%d",IntraGenicLen,cMinFlankLen,cMaxFlankLen);
			exit(1);
			}
		if((InterGenicLen + IntraGenicLen) < 10)
			{
			printf("\nError: Total intragenic and intergenic (%d) flank lengths must total at least 10",IntraGenicLen+InterGenicLen);
			exit(1);
			}

		ProfPol = (etPOPpolicy)(profpol->count ? profpol->ival[0] : ePOPdefault);
		if(ProfPol < ePOPdefault || ProfPol > ePOPplaceholder)
			{
			printf("\nProfile overlap policy '-p%d' outside of range 0..%d",ProfPol,ePOPplaceholder-1);
			exit(1);
			}

		if(ProfPol != ePOPdefault)
			{
			KDist = kdist->count ? kdist->ival[0] :  cDfltKDist;
			if(KDist < cMinKDist || KDist > cMaxKDist)
				{
				printf("\nProfile overlap distance '-k%d' outside of range %d..%d",KDist,cMinKDist,cMaxKDist);
				exit(1);
				}
			}
		else
			KDist = 0;
		if(excludefile->count)
			{
			strncpy(szExcludeFile,excludefile->filename[0],_MAX_PATH);
			szExcludeFile[_MAX_PATH-1] = '\0';
			}
		else
			szExcludeFile[0] = '\0';
		}
	else
		{
		KDist = 0;
		ProfPol = ePOPdefault;
		szExcludeFile[0] = '\0';
		IntraGenicLen = 0;
		InterGenicLen = 0;
		}

	for(NumInputFiles=Idx=0;NumInputFiles < cMaxInFileSpecs && Idx < infiles->count; Idx++)
		{
		pszInFileSpecs[Idx] = NULL;
		if(pszInFileSpecs[NumInputFiles] == NULL)
			pszInFileSpecs[NumInputFiles] = new char [_MAX_PATH];
		strncpy(pszInFileSpecs[NumInputFiles],infiles->filename[Idx],_MAX_PATH);
		pszInFileSpecs[NumInputFiles][_MAX_PATH-1] = '\0';
		CUtility::TrimQuotedWhitespcExtd(pszInFileSpecs[NumInputFiles]);
		if(pszInFileSpecs[NumInputFiles][0] != '\0')
			NumInputFiles++;
		}

	if(!NumInputFiles)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: After removal of whitespace, no input file(s) specified with '-i<filespec>' option)\n");
		exit(1);
		}

	strcpy(szRsltsFile,outfile->filename[0]);

	strcpy(szFeatureFile,featurefile->filename[0]);

			// now that command parameters have been parsed then initialise diagnostics log system
	if(!gDiagnostics.Open(szLogFile,(etDiagLevel)iScreenLogLevel,(etDiagLevel)iFileLogLevel,true))
		{
		printf("\nError: Unable to start diagnostics subsystem.");
		if(szLogFile[0] != '\0')
			printf(" Most likely cause is that logfile '%s' can't be opened/created",szLogFile);
		exit(1);
		}
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Version: %s Processing parameters:",cpszProgVer);
	const char *pszDescr;
	switch(PMode) {
		case ePMauto:
			pszDescr = "Auto-determination of input file formats";
			break;
		case ePMCsv:
			pszDescr = "CSV loci";
			break;
		case ePMBed:			// UCSC BED format
			pszDescr = "UCSC BED loci";
			break;
		case ePMSam:			// SAM format
			pszDescr = "SAM loci";
			break;
		}
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"input file format : '%s'",pszDescr);

	switch(FMode) {
		case eFMdefault:
			pszDescr = "CSV";
			break;
		}
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"output format is : '%s'",pszDescr);

	switch(ReadProfile) {
		case eRPdensity:
			pszDescr = "overlap read density";
			break;
		case eRPstartcnts:				
			pszDescr = "read start counts";
			break;
		case eRPstartuniques:			
			pszDescr = "unique read start loci";
			break;
		}


	if(ReadProfile > eRPdensity)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Peak window length: %d",PkWinLen);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"bin reads by : %s",pszDescr);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"process for this strand only: '%c'",(char)Strand);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Intergenic flank length: %d",InterGenicLen);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Intragenic flank length: %d",IntraGenicLen);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"number of bins: %d",NumBins);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"flanks centered around: %s",Feature2Txt((etPFLFeature)Feature));


	if(Feature != ePFLIsoforms)
		{
		switch(ProfPol) {
			case ePOPdefault:	
				pszDescr = "accept all reads";
				break;
			case ePOPkdist:
				pszDescr = "do not accept reads within -k distance of overlapping other features";
				break;
			default:
				pszDescr = "Unsupported profile policy";
				break;
			}
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Profile overlap policy is : '%s'",pszDescr);
		if(ProfPol == ePOPkdist)
			{
			gDiagnostics.DiagOutMsgOnly(eDLInfo," overlap distance is : %d",KDist);
			if(szExcludeFile[0] != '\0')
				gDiagnostics.DiagOutMsgOnly(eDLInfo,"Exclude reads which also overlap features in this BED file: '%s'",szExcludeFile);
			}
		}

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"feature file: '%s'",szFeatureFile);
	for(Idx = 0; Idx < NumInputFiles; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"input file (%d): '%s'",Idx+1, pszInFileSpecs[Idx]);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Output to file: '%s'",szRsltsFile);

#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	gStopWatch.Start();
	Rslt = Process(PMode,FMode,ReadProfile,PkWinLen,(char)Strand,Feature,InterGenicLen,IntraGenicLen,NumBins,ProfPol,KDist,szExcludeFile,szFeatureFile,NumInputFiles,pszInFileSpecs,szRsltsFile);

	gStopWatch.Stop();
	Rslt = Rslt >=0 ? 0 : 1;
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Exit code: %d Total processing time: %s",Rslt,gStopWatch.Read());
	exit(Rslt);
	}
else
	{
	printf("\n%s Generate Profiles, Version %s\n", gszProcName,cpszProgVer);
	arg_print_errors(stdout,end,gszProcName);
	arg_print_syntax(stdout,argtable,"\nUse '-h' to view option and parameter usage\n");
	exit(1);
	}
}

char *
Feature2Txt(etPFLFeature Feature)
{
switch(Feature) {
	case ePFLIsoforms:
		return((char *)"Gene isoforms");
	case ePFLTSS:
		return((char *)"TSS flanks");
	case ePFLTES:
		return((char *)"TES flanks");
	default:
		break;
	}
return((char *)"Unsupported");
}

char *
TrimWhitespace(char *pTxt)
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

void
Init(void)
{
m_pCSVFile = NULL;
m_pBEDFile = NULL;
m_pExclBEDFile = NULL;

m_hRsltsFile = -1;			// handle for opened results file

for(int ChromIdx = 0; ChromIdx < cMaxChromCov; ChromIdx++)
	{
	m_ChromCnts[ChromIdx].AllocCovCnts = 0;
	m_ChromCnts[ChromIdx].StartOfs = 0;
	m_ChromCnts[ChromIdx].EndOfs = 0;
	m_ChromCnts[ChromIdx].pCovCnts = NULL;
	m_ChromCnts[ChromIdx].szChrom[0] = '\0';
	}
m_NumChromsCov = 0;	
}

void
Reset(void)
{
if(m_hRsltsFile != -1)
	{
	close(m_hRsltsFile);
	m_hRsltsFile = -1;
	}

if(m_pBEDFile != NULL)
	{
	delete m_pBEDFile;
	m_pBEDFile = NULL;
	}
if(m_pExclBEDFile != NULL)
	{
	delete m_pExclBEDFile;
	m_pExclBEDFile = NULL;
	}
if(m_pCSVFile != NULL)
	{
	delete m_pCSVFile;
	m_pCSVFile = NULL;
	}

if(m_NumChromsCov > 0)
	{
	for(int ChromIdx = 0; ChromIdx < m_NumChromsCov; ChromIdx++)
		{
		m_ChromCnts[ChromIdx].StartOfs = 0;
		m_ChromCnts[ChromIdx].EndOfs = 0;
		if(m_ChromCnts[ChromIdx].pCovCnts != NULL)
			{
#ifdef _WIN32
			free(m_ChromCnts[ChromIdx].pCovCnts);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
			if(m_ChromCnts[ChromIdx].pCovCnts != MAP_FAILED)
				munmap(m_ChromCnts[ChromIdx].pCovCnts,m_ChromCnts[ChromIdx].AllocCovCnts * sizeof(UINT32));
#endif
			m_ChromCnts[ChromIdx].pCovCnts = NULL;
			}
		m_ChromCnts[ChromIdx].AllocCovCnts = 0;
		m_ChromCnts[ChromIdx].szChrom[0] = '\0';
		}
	m_NumChromsCov = 0;	
	}
}



int
BuildReadCoverage(etReadProfile ReadProfile,		// profile for: 0 read density, 1: read starts, 2:  unique read starts (default is for read density)
			   char *pszChrom,		// coverage is onto this chrom
		      char Strand,				// read was aligned to this strand
			  int StartOfs,				// coverage start at this offset 
			  int EndOfs,				// and ends at this offset inclusive
			  int Cnt)					// increment coverage by this
{
tsChromCnts *pChrom;
int ChromIdx;
int AllocCovCnts;
UINT32 *pCovCnts;
size_t ReallocTo;

if(pszChrom == NULL || pszChrom[0] == '\0')
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"BuildReadCoverage: No chromosome specified");
	Reset();
	return(eBSFerrChrom);
	}

// arbitary count clamping to be in range 1..10000 inclusive
if(Cnt < 1)
	Cnt = 1;
else
	if(Cnt > 10000)
		Cnt = 10000;

// ensure StartOfs and EndOfs are both >= 0
if(StartOfs < 0)
	StartOfs = 0;
if(EndOfs < 0)
	EndOfs = 0;

// ensure StartOfs <= EndOfs
if(StartOfs > EndOfs)
	{
	int TmpOfs = EndOfs;
	EndOfs = StartOfs;
	StartOfs = TmpOfs;
	}

// check if this is a new chrom or if coverage is onto an existing chrom
pChrom = &m_ChromCnts[0];
ChromIdx = 0;
if(m_NumChromsCov > 0)
	{
	for(ChromIdx = 0; ChromIdx < m_NumChromsCov; ChromIdx++,pChrom++)
		{
		if(!strnicmp(pszChrom,pChrom->szChrom,cMaxDatasetSpeciesChrom))
			break;
		}
	}
if(ChromIdx == m_NumChromsCov)	// if a new or first chrom
	{
	strncpy(pChrom->szChrom,pszChrom,cMaxDatasetSpeciesChrom);
	pChrom->szChrom[cMaxDatasetSpeciesChrom-1] = 0;
	pChrom->StartOfs = StartOfs;
	pChrom->EndOfs = EndOfs;
	AllocCovCnts = EndOfs + cAllocCovCnts;
	ReallocTo =  AllocCovCnts * sizeof(UINT32);
#ifdef _WIN32
	pChrom->pCovCnts = (UINT32 *) malloc(ReallocTo);	// initial and perhaps the only allocation
	if(pChrom->pCovCnts == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"BuildReadCoverage: Memory allocation of %lld bytes - %s",(INT64)ReallocTo,strerror(errno));
		Reset();
		return(eBSFerrMem);
		}
#else
	// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
	pChrom->pCovCnts = (UINT32 *)mmap(NULL,ReallocTo, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
	if(pChrom->pCovCnts == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"BuildReadCoverage: Memory allocation of %lld bytes through mmap()  failed - %s",(INT64)ReallocTo,strerror(errno));
		pChrom->pCovCnts = NULL;
		Reset();
		return(eBSFerrMem);
		}
#endif
	pChrom->AllocCovCnts = AllocCovCnts;
	memset(pChrom->pCovCnts,0,ReallocTo);
	m_NumChromsCov += 1;
	}

// check if chrom coverage cnts needs to be extended
if(EndOfs >= pChrom->AllocCovCnts)
	{
	AllocCovCnts = EndOfs + cAllocCovCnts;
	ReallocTo = AllocCovCnts * sizeof(UINT32);
#ifdef _WIN32
	pCovCnts = (UINT32 *) realloc(pChrom->pCovCnts,ReallocTo);
#else
	pCovCnts = (UINT32 *)mremap(pChrom->pCovCnts,pChrom->AllocCovCnts * sizeof(UINT32),ReallocTo,MREMAP_MAYMOVE);
	if(pCovCnts == MAP_FAILED)
		pCovCnts = NULL;
#endif
	if(pCovCnts == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"BuildReadCoverage: Memory re-allocation to %d bytes - %s",ReallocTo,strerror(errno));
		Reset();
		return(eBSFerrMem);
		}
	pChrom->pCovCnts = pCovCnts;
	memset(&pChrom->pCovCnts[pChrom->AllocCovCnts],0,(AllocCovCnts - pChrom->AllocCovCnts) * sizeof(UINT32));
	pChrom->AllocCovCnts = AllocCovCnts;
	}

if(EndOfs > pChrom->EndOfs)
	pChrom->EndOfs = EndOfs;
if(StartOfs < pChrom->StartOfs)
	pChrom->StartOfs = StartOfs;

pCovCnts = &pChrom->pCovCnts[StartOfs];

if(ReadProfile != eRPdensity)		// if using read starts instead of the default overlap read densities
	{
	if(Strand == '-')
		pCovCnts = &pChrom->pCovCnts[EndOfs];
	if((cMaxAccumCnt - *pCovCnts) > (UINT32)Cnt)
		*pCovCnts += (UINT32)Cnt;
	else
		*pCovCnts = cMaxAccumCnt;
	}
else
	{
	while(StartOfs++ <= EndOfs)
		{
		// clamp accumulated cnts to be no greater than cMaxAccumCnt
		if((cMaxAccumCnt - *pCovCnts) > (UINT32)Cnt)
			*pCovCnts += (UINT32)Cnt;
		else
			*pCovCnts = cMaxAccumCnt;
		pCovCnts += 1;
		}
	}
return(eBSFSuccess);
}



// RetainReadCoverage
// Marks read coverage regions as being regions to be retained by
// setting MSB of coverage counts 
// Later when coverage is being reported, only those counts with
// bit 15 set will be reported on 
int
RetainReadCoverage(char *pszChrom,		// filtering is onto this chrom
			  int StartOfs,				// filtering is to start at this offset 
			  int EndOfs)				// and ends at this offset inclusive
{
tsChromCnts *pChrom;
int ChromIdx;
UINT32 *pCovCnts;

if(pszChrom == NULL || pszChrom[0] == '\0')
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"FilterReadCoverage: No chromosome specified");
	Reset();
	return(eBSFerrChrom);
	}

// ensure StartOfs <= EndOfs
if(StartOfs > EndOfs)
	{
	int TmpOfs = EndOfs;
	EndOfs = StartOfs;
	StartOfs = TmpOfs;
	}

// check that coverage is onto an existing chrom
pChrom = &m_ChromCnts[0];
ChromIdx = 0;
if(m_NumChromsCov > 0)
	{
	for(ChromIdx = 0; ChromIdx < m_NumChromsCov; ChromIdx++,pChrom++)
		{
		if(!strnicmp(pszChrom,pChrom->szChrom,cMaxDatasetSpeciesChrom))
			break;
		}
	}
if(ChromIdx == m_NumChromsCov)	// not an existing chrom so no filtering required
	return(eBSFSuccess);

if(StartOfs >= pChrom->AllocCovCnts || EndOfs < pChrom->StartOfs)
	return(eBSFSuccess);

if(EndOfs >= pChrom->AllocCovCnts)
	EndOfs = pChrom->AllocCovCnts - 1;

if(StartOfs < pChrom->StartOfs)
	pChrom->StartOfs = StartOfs;

pCovCnts = &pChrom->pCovCnts[StartOfs];
while(StartOfs++ <= EndOfs)
	*pCovCnts++ |= cRetainCnt;

return(eBSFSuccess);
}





int
WriteRegions(char *pszSrcFile,				// file from which read loci were sourced
			 etFMode FMode,					// output in this format to m_hRsltsFile
			char *pszTitle,					// CSV species or title used if BED output format
			 char Strand,					// ROI are on this strand 
			  int MinMedianCov,				// minimum median coverage required in reported regions
				int MinRegionLen,			// report regions which are of at least this length
				int MaxGapLen,				// report regions containing no read gaps of <= this length 
				bool bNoFilter)				// true if all counts to be processed
{
int BuffIdx;
char szLineBuff[8096];
tsChromCnts *pChrom;
UINT32 *pCnts;
int Cnts;
UINT32 SumCnts;
int ChromIdx;
int SeqIdx;
int StartOfRegion;
int EndOfRegion;
int CurGapLen;
int NumMedCnts;	
int SubRegionLen;
int TotRegionLen;
int NumRegions;
int Score;

if(FMode == 0) //eFMbed)
	{
	BuffIdx = sprintf(szLineBuff,"track type=bed name=\"%s\" description=\"%s\"\n",pszTitle,pszTitle);
	CUtility::SafeWrite(m_hRsltsFile,szLineBuff,BuffIdx);
	}

BuffIdx = 0;
pChrom = &m_ChromCnts[0];
NumRegions = 0;
for(ChromIdx = 0 ; ChromIdx < m_NumChromsCov; ChromIdx++,pChrom++)
	{
	CurGapLen = 0;				
	NumMedCnts = 0;
	SumCnts = 0;
	SubRegionLen = -1;
	StartOfRegion = -1;
	EndOfRegion = -1;
	
	pCnts = &pChrom->pCovCnts[pChrom->StartOfs];
	for(SeqIdx = pChrom->StartOfs; SeqIdx <= pChrom->EndOfs+1; SeqIdx++,pCnts++)
		{
		if(bNoFilter || *pCnts >= cRetainCnt)
			Cnts = *pCnts & cMaxAccumCnt;
		else
			Cnts = 0;
		if(Cnts == 0 ||				// if in gap with no reads
			SeqIdx > pChrom->EndOfs)    // or past the last read
			{
			if(StartOfRegion < 0)		// if not actually started a region then continue until region or new chrom
				continue;
			if(CurGapLen == 0 || SeqIdx > pChrom->EndOfs)	// if gap just started then check if subregion can be accepted
				{			
				if(NumMedCnts >= (SubRegionLen+1)/2)	// more than 50% of bases in subregion had at least the minimum coverage?
					EndOfRegion = SeqIdx - 1;
				}
			if(CurGapLen++ == MaxGapLen ||	// if gap is too large then terminate region 
				SeqIdx > pChrom->EndOfs)	// or at end of chrom also terminates region
				{
				if(EndOfRegion != -1)
					{
					TotRegionLen = 1 + EndOfRegion - StartOfRegion; // region meets minimum length requirements?
					if(TotRegionLen >= MinRegionLen)
						{
						// output region here!!	
						NumRegions += 1;
						Score = SumCnts/(UINT32)TotRegionLen;
						switch(FMode) {
							case eFMdefault:
								BuffIdx+=sprintf(&szLineBuff[BuffIdx],"%d,\"ROI\",\"%s\",\"%s\",%d,%d,%d,\"%c\"\n",NumRegions,pszTitle,pChrom->szChrom,StartOfRegion,EndOfRegion,TotRegionLen,Strand);
								break;

							default: //case eFMbed:
								BuffIdx+=sprintf(&szLineBuff[BuffIdx],"%s\t%d\t%d\tROI%d\t%d\t%c\n",
									pChrom->szChrom,StartOfRegion,EndOfRegion,NumRegions,Score,Strand != '-' ? '+' : Strand);
								break;
							}
								
						if((BuffIdx + 200)> sizeof(szLineBuff))
							{
							CUtility::SafeWrite(m_hRsltsFile,szLineBuff,BuffIdx);
							BuffIdx = 0;
							}
						}
					}
				StartOfRegion = -1;
				EndOfRegion = -1;
				SumCnts = 0;
				}
			SubRegionLen = -1;
			continue;
			}

		CurGapLen = 0;
		if(StartOfRegion < 0)	// starting a putative region?
			{
			StartOfRegion = SeqIdx;
			SubRegionLen = -1;
			SumCnts = 0;
			}

		if(SubRegionLen < 0)	// starting a subregion?
			{
			NumMedCnts = 0;
			SubRegionLen = 1;
			}
		else
			SubRegionLen += 1;

		if(Cnts >= MinMedianCov)
			NumMedCnts += 1;
		SumCnts += Cnts;
		}
	}
if(BuffIdx)
	CUtility::SafeWrite(m_hRsltsFile,szLineBuff,BuffIdx);

return(NumRegions);
}

int GenCSVdensity(etReadProfile ReadProfile,	// profile for: 0 read density, 1: read starts, 2:  unique read starts (default is for read density)
			     char Strand,					// ROI strand
				 char *pszInFile)				// CSV loci input file
{
int Rslt;
int NumEls;
int NumFields;
char *pszChrom;
int StartLoci;
int EndLoci;
char *pszStrand;

if((m_pCSVFile = new CCSVFile)==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CCSVFile");
	Reset();
	return(eBSFerrObj);
	}

if((Rslt=m_pCSVFile->Open(pszInFile)) !=eBSFSuccess)
	{
	while(m_pCSVFile->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,m_pCSVFile->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open '%s' for processing",pszInFile);
	Reset();
	return(eBSFerrOpnFile);
	}
NumEls = 0;
while((Rslt=m_pCSVFile->NextLine()) > 0)	// onto next line containing fields
	{
	NumFields = m_pCSVFile->GetCurFields();
	if(NumFields < 7)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expected at least 7 fields in '%s', GetCurFields() returned '%d'",pszInFile,NumFields);
		return(eBSFerrFieldCnt);
		}

	if(NumFields >= 8)
		m_pCSVFile->GetText(8,&pszStrand);
	else
		pszStrand = (char *)"+";

	if(Strand != '*')
		{
		if(NumFields < 8)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Strand filtering, expected at least 8 fields in '%s', GetCurFields() returned '%d'",pszInFile,NumFields);
			return(eBSFerrFieldCnt);
			}
		if(*pszStrand != Strand)
			continue;
		}

	NumEls += 1;
	m_pCSVFile->GetText(4,&pszChrom);
	m_pCSVFile->GetInt(5,&StartLoci);
	m_pCSVFile->GetInt(6,&EndLoci);

	if((Rslt=BuildReadCoverage(ReadProfile,pszChrom,*pszStrand,StartLoci,EndLoci,1))!=eBSFSuccess)
		break;
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Iteration of %d elements completed",NumEls);
delete m_pCSVFile;
m_pCSVFile = NULL;
return(Rslt);
}

int GenBEDdensity(etReadProfile ReadProfile,		// profile for: 0 read density, 1: read starts, 2:  unique read starts (default is for read density)
				  char ROIStrand,				// ROI strand
				 char *pszInFile)				// UCSC BED input file
{
int Rslt;
int CurFeatureID;
int StartLoci;
int EndLoci;
int Score;
char szChrom[128];
char szPrevChrom[128];
char szFeatName[128];
char Strand;
int RefID;
int IntergenicStart;
int NumEls;

if((m_pBEDFile = new CBEDfile)==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CBEDfile");
	Reset();
	return(eBSFerrObj);
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Starting to load %s",pszInFile);
if((Rslt=m_pBEDFile->Open(pszInFile,eBTAnyBed)) !=eBSFSuccess)
	{
	while(m_pBEDFile->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,m_pBEDFile->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open '%s' for processing",pszInFile);
	Reset();
	return(eBSFerrOpnFile);
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Load completed, iterating elements...");
if(m_pBEDFile->ContainsGeneDetail())			// returns true if file contains gene detail (utr/cds/intron etc)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"bed file %s contains gene regions",pszInFile);
	Reset();
	return(eBSFerrFeature);
	}

// now iterate over the features
szPrevChrom[0] = '\0';
CurFeatureID = 0;
IntergenicStart = 0;
RefID = 0;
Rslt = eBSFSuccess;
NumEls = 0;
while(Rslt == eBSFSuccess && (CurFeatureID = m_pBEDFile->GetNextFeatureID(CurFeatureID)) > 0)
	{
	NumEls += 1;

	m_pBEDFile->GetFeature(CurFeatureID,	// feature instance identifier
				szFeatName,				// where to return feature name
				szChrom,				// where to return chromosome name
				&StartLoci,				// where to return feature start on chromosome (0..n) 
				&EndLoci,				// where to return feature end on chromosome
 				&Score,					// where to return score
 				&Strand);				// where to return strand

	if(CurFeatureID == 1 || stricmp(szChrom,szPrevChrom))	// if new chromosome then reset IntergenicStart
		{
		strcpy(szPrevChrom,szChrom);
		IntergenicStart = 0;
		}

	if(ROIStrand != '*')
		{
		if(ROIStrand != Strand)
			continue;
		}
	
	Rslt=BuildReadCoverage(ReadProfile,szChrom,Strand,StartLoci,EndLoci,1);
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Iteration of %d elements completed",NumEls);
delete m_pBEDFile;
m_pBEDFile = NULL;
return(Rslt);
}

int GenSAMdensity(etReadProfile ReadProfile,	// profile for: 0 read density, 1: read starts, 2:  unique read starts (default is for read density)
				  char ROIStrand,				// ROI strand
				 char *pszInFile)				// SAM input file
{
teBSFrsltCodes Rslt;
FILE *pSAMStream;
int NumEls;
char szLine[32000];				// buffer input lines
char *pTxt;
char szDescriptor[18];			// parsed out descriptor
int Flags;						// parsed out flags
char szChrom[128];				// parsed out chrom
int StartLoci;					// start loci
char Strand;
char szCigar[128];
char szRNext[128];
int MAPQ;
int PNext;
int TLen;
int LineNum;
int NumFields;

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Starting to load %s",pszInFile);

// open SAM for reading
if(pszInFile == NULL || *pszInFile == '\0')
	return(eBSFerrParams);
if((pSAMStream = fopen(pszInFile,"r"))==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to fopen SAM format file %s error: %s",pszInFile,strerror(errno));
	return(eBSFerrOpnFile);
	}

NumEls = 0;
LineNum = 0;
while(fgets(szLine,sizeof(szLine),pSAMStream)!= NULL)
	{
	LineNum += 1;
	pTxt = TrimWhitespace(szLine);
	if(*pTxt=='\0' || *pTxt=='@')	// simply slough lines which were just whitespace or start with '@'
		continue;

	NumEls += 1;

	// expecting to parse as "%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t", szDescriptor, Flags, m_szSAMTargChromName, StartLoci+1,MAPQ,szCigar,pszRNext,PNext,TLen);
	// interest is in the chromname, startloci, length, and strand (encoded into Flags)
	NumFields = sscanf(szLine,"%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t",szDescriptor, &Flags, szChrom, &StartLoci,&MAPQ,szCigar,szRNext,&PNext,&TLen);
	if(NumFields != 9)	// expected minimum number of fields?
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expected to parse at least 9 fields (%d parsed) from SAM file at line %d",NumFields,LineNum);
		return(eBSFerrLocField);
		}

	NumEls += 1;
	Strand = Flags & 0x010 ? '-' : '+';
	if(ROIStrand != '*')
		{
		if(ROIStrand != Strand)
			continue;
		}

	Rslt=(teBSFrsltCodes)BuildReadCoverage(ReadProfile,szChrom,Strand,StartLoci-1,StartLoci+TLen,1);
	if(Rslt == eBSFSuccess)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddElCore failed");
		break;
		}

	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Iteration of %d elements completed",NumEls);
fclose(pSAMStream);
return(Rslt);
}

int
ReportCenteredAt(etReadProfile ReadProfile,		// profile for: 0 read density, 1: read starts, 2:  unique read starts (default is for read density)
				 int PkWinLen,					// window length over which to calculate the peak counts
				 char *pszIdentifier,	// identifies this instance
				 int NumExons,			// number of exons
				 char Strand,			// on this strand
				 char *pszChrom,		// is on this chrom	
				 UINT32 LociOfInterest,	// loci of interest (TSS or TES)
 				 int NumBins,			// profile counts into this many bins
				 int IntergenicLen,		// profiling counts for this length intergenic
				 int IntragenicLen,		// profiling counts for this length intergenic
				 int StartOfs,			// offset in profile at which count profiling starts, if < 0 then do not profile, set all bins to 0
 				 int StartLoci,			// start count profiling at this loci
				 int EndLoci)			// stop profiling counts at this loci
{
int BinIdx;
int ChromIdx;
tsChromCnts *pChrom;
int TotFlankLen;
UINT32 *pCnts;
UINT32 *pWinCnts;
UINT32 BinCnts[cMaxNumBins];
char szOutBuff[cMaxNumBins*8];
int BuffIdx;
int TotCnts;
int TotBinCnts;
int Cnts;
int PeakWinCnts;
int PeakWinIdx;
int CurWinCnts;

// validate that parameters are reasonable
if(pszIdentifier == NULL || pszIdentifier[0] == '\0' || !(Strand == '+' || Strand == '-'))
	return(-1);
if(LociOfInterest < 0 || NumBins < 1 || IntergenicLen < 0 || IntragenicLen < 0 || (IntergenicLen + IntragenicLen) < 1)
	return(-1);

// check that coverage is onto an existing chrom
pChrom = &m_ChromCnts[0];
ChromIdx = 0;
if(m_NumChromsCov > 0)
	{
	for(ChromIdx = 0; ChromIdx < m_NumChromsCov; ChromIdx++,pChrom++)
		{
		if(!strnicmp(pszChrom,pChrom->szChrom,cMaxDatasetSpeciesChrom))
			break;
		}
	}
if(ChromIdx == m_NumChromsCov)	// not an existing chrom so no reporting required
	return(eBSFSuccess);

TotFlankLen = IntergenicLen + IntragenicLen;
memset(BinCnts,0,NumBins * sizeof(UINT32));
TotCnts = 0;
TotBinCnts = 0;
PeakWinCnts = 0;
PeakWinIdx = 0;
CurWinCnts = 0;
if(StartOfs >= 0)
	{
	if(Strand == '+')
		{
		pCnts = &pChrom->pCovCnts[StartLoci];
		pWinCnts = pCnts;
		PeakWinIdx = 0;
		while(StartLoci <= EndLoci && StartLoci <= pChrom->EndOfs)
			{
			Cnts = *pCnts;
			if(ReadProfile == eRPstartuniques && Cnts > 1)
				Cnts = 1;
			BinIdx = (StartOfs++ * NumBins)/TotFlankLen;
			BinCnts[BinIdx] += Cnts;
			TotBinCnts += Cnts;
			TotCnts += *pCnts++;
			StartLoci += 1;
			CurWinCnts += Cnts;
			if(CurWinCnts > PeakWinCnts)
				PeakWinCnts = CurWinCnts;
			PeakWinIdx += 1;
			if(PeakWinIdx < PkWinLen)
				continue;
			Cnts = *pWinCnts++;
			if(ReadProfile == eRPstartuniques && Cnts > 1)
				Cnts = 1;
			CurWinCnts -= Cnts;
			}
		}
	else
		{
		pCnts = &pChrom->pCovCnts[StartLoci];
		pWinCnts = pCnts;
		PeakWinIdx = 0;
		while(StartLoci >= EndLoci)
			{
			if(StartLoci <= pChrom->EndOfs)
				{
				Cnts = *pCnts;
				if(ReadProfile == eRPstartuniques && Cnts > 1)
					Cnts = 1;
				BinIdx = (StartOfs * NumBins)/TotFlankLen;
				BinCnts[BinIdx] += Cnts;
				TotBinCnts += Cnts;
				TotCnts += *pCnts;
				CurWinCnts += Cnts;
				if(CurWinCnts > PeakWinCnts)
					PeakWinCnts = CurWinCnts;
				PeakWinIdx += 1;
				if(PeakWinIdx >= PkWinLen)
					{
					Cnts = *pWinCnts--;
					if(ReadProfile == eRPstartuniques && Cnts > 1)
						Cnts = 1;
					CurWinCnts -= Cnts;
					}
				}
			pCnts -= 1;
			StartOfs += 1;
			StartLoci -= 1;
			}
		}
	}

BuffIdx = sprintf(szOutBuff,"\"%s\",\"%s\",\"%c\",%d,%d,%d,%d,%d,%d,%d",pszIdentifier,pszChrom,Strand,NumExons,IntergenicLen,IntragenicLen,LociOfInterest,TotCnts,PeakWinCnts,TotBinCnts);

for(BinIdx = 0; BinIdx < NumBins; BinIdx++)
	BuffIdx += sprintf(&szOutBuff[BuffIdx],",%d",BinCnts[BinIdx]);
BuffIdx += sprintf(&szOutBuff[BuffIdx],"\n");
CUtility::SafeWrite(m_hRsltsFile,szOutBuff,BuffIdx);
return(0);
}

typedef struct TAG_sIntron {
	int Start;	// intron start on chrom
	int End;	// intron end on chrom
} tsIntron;

const int cMaxIntrons = 50000;				// allow for upto 50K introns in any gene isoform	
int m_NumIntrons;						// current gene isoform has this many introns
tsIntron m_Introns[cMaxIntrons];

int
ReportIsoform(etReadProfile ReadProfile,		// profile for: 0 read density, 1: read starts, 2:  unique read starts (default is for read density)
				 int PkWinLen,					// window length over which to calculate the peak counts
				 int CurFeatureID,			// identifies isoform to report on
				 int InterGenicLen,		// 5' upstream flank length
				 int IntraGenicLen,		// 3' upstream flank length
				 int NumBins)			// profile into this many bins over length of gene isoform + upstream and downstream
{
int TranscribedLen;
int StartLoci;
int EndLoci;
int LociOfInterest;
char FeatStrand;
int Score;
char szChrom[128];
char szFeatName[128];
int Ofs;
int BinIdx;
int ChromIdx;
tsChromCnts *pChrom;
int Start;
int End;
int CurIntron;
int TransIdx;
UINT32 *pCnts;
UINT32 *pWinCnts;
UINT32 BinCnts[cMaxNumBins];
char szOutBuff[cMaxNumBins*8];
int BuffIdx;
int TotCnts;
int Cnts;
int PeakWinCnts;
int PeakWinIdx;
int CurWinCnts;
int TotBinCnts;
int TotalLen;

m_pBEDFile->GetFeature(CurFeatureID,	// feature instance identifier
			szFeatName,				// where to return feature name
			szChrom,				// where to return chromosome name
			&StartLoci,				// where to return feature start on chromosome (0..n) 
			&EndLoci,				// where to return feature end on chromosome
			&Score,					// where to return score
			&FeatStrand);			// where to return strand



// check that coverage is onto an existing chrom
pChrom = &m_ChromCnts[0];
ChromIdx = 0;
if(m_NumChromsCov > 0)
	{
	for(ChromIdx = 0; ChromIdx < m_NumChromsCov; ChromIdx++,pChrom++)
		{
		if(!strnicmp(szChrom,pChrom->szChrom,cMaxDatasetSpeciesChrom))
			break;
		}
	}
if(ChromIdx == m_NumChromsCov)	// not an existing chrom so no reporting required
	return(eBSFSuccess);

memset(BinCnts,0,NumBins * sizeof(UINT32));

if((m_NumIntrons = m_pBEDFile->GetNumIntrons(CurFeatureID)) > 0)
	{
	for(CurIntron = 0; CurIntron < m_NumIntrons; CurIntron++)
		{
		m_Introns[CurIntron].Start = m_pBEDFile->GetIntronStart(CurFeatureID,CurIntron+1);
		m_Introns[CurIntron].End = m_pBEDFile->GetIntronEnd(CurFeatureID,CurIntron+1);
		}
	}

if(FeatStrand == '+')
	{
	Start = StartLoci - InterGenicLen;
	End = EndLoci + IntraGenicLen;
	LociOfInterest = Start;
	}
else
	{
	Start = StartLoci - IntraGenicLen;
	End = EndLoci + InterGenicLen;
	LociOfInterest = End;
	}
	
TranscribedLen = m_pBEDFile->GetTranscribedLen(CurFeatureID) + InterGenicLen + IntraGenicLen;
TotalLen = 1 + InterGenicLen + IntraGenicLen + End - Start;
TotCnts = 0;
TotBinCnts = 0;
PeakWinCnts = 0;
PeakWinIdx = 0;
CurWinCnts = 0;

if(FeatStrand == '+')
	{
	pCnts = &pChrom->pCovCnts[Start];
	pWinCnts = pCnts;
	TransIdx = 0;
	CurIntron = 0;
	for(Ofs = Start; Ofs <= End; Ofs++, pCnts++)
		{
		if(Ofs < 0)
			continue;

		if(CurIntron < m_NumIntrons)
			{
			if(Ofs >= m_Introns[CurIntron].End)
				{
				CurIntron += 1;
				continue;
				}
			if(Ofs >= m_Introns[CurIntron].Start)
				continue;
			}
		
		Cnts = *pCnts;
		if(ReadProfile == eRPstartuniques && Cnts > 1)
			Cnts = 1;
		BinIdx = (TransIdx++ * NumBins)/TranscribedLen;
		if(Ofs <= pChrom->EndOfs)
			{
			BinCnts[BinIdx] += Cnts;
			TotBinCnts += Cnts;
			TotCnts += *pCnts;

			CurWinCnts += Cnts;
			if(CurWinCnts > PeakWinCnts)
				PeakWinCnts = CurWinCnts;
			PeakWinIdx += 1;
			if(PeakWinIdx < PkWinLen)
				continue;
			Cnts = *pWinCnts++;
			if(ReadProfile == eRPstartuniques && Cnts > 1)
				Cnts = 1;
			CurWinCnts -= Cnts;
			}
		}
	}
else
	{
	pCnts = &pChrom->pCovCnts[End];
	pWinCnts = pCnts;
	TransIdx = 0;
	CurIntron = m_NumIntrons - 1;
	for(Ofs = End; Ofs >= Start; Ofs--, pCnts--)
		{
		if(Ofs < 0)
			break;
		if(CurIntron >= 0)
			{
			if(Ofs <= m_Introns[CurIntron].End && Ofs > m_Introns[CurIntron].Start)
				continue;
			if(Ofs == m_Introns[CurIntron].Start)
				{
				CurIntron -= 1;
				continue;
				}
			}

		Cnts = *pCnts;
		if(ReadProfile == eRPstartuniques && Cnts > 1)
			Cnts = 1;
		BinIdx = (TransIdx++ * NumBins)/TranscribedLen;
		if(Ofs <= pChrom->EndOfs)
			{
			BinCnts[BinIdx] += *pCnts;
			TotBinCnts += Cnts;
			TotCnts += *pCnts;
			CurWinCnts += Cnts;
			if(CurWinCnts > PeakWinCnts)
				PeakWinCnts = CurWinCnts;
			PeakWinIdx += 1;
			if(PeakWinIdx >= PkWinLen)
				{
				Cnts = *pWinCnts--;
				if(ReadProfile == eRPstartuniques && Cnts > 1)
					Cnts = 1;
				CurWinCnts -= Cnts;
				}

			}
		}
	}

BuffIdx = sprintf(szOutBuff,"\"%s\",\"%s\",\"%c\",%d,%d,%d,%d,%d,%d,%d",szFeatName,szChrom,FeatStrand,m_NumIntrons+1,TotalLen,TranscribedLen,LociOfInterest,TotCnts,PeakWinCnts,TotBinCnts);
for(BinIdx = 0; BinIdx < NumBins; BinIdx++)
	BuffIdx += sprintf(&szOutBuff[BuffIdx],",%d",BinCnts[BinIdx]);
BuffIdx += sprintf(&szOutBuff[BuffIdx],"\n");
CUtility::SafeWrite(m_hRsltsFile,szOutBuff,BuffIdx);
return(0);
}

// returns loci at which 1st overlap with feature in BED file is detected
// starting from StartLoci through to EndLoci
// if EndLoci is >= StartLoci then search is towards the 3' end of targeted chrom
// otherwise if EndLoci < StartLoci then search is towards the 5' end of targeted chrom
// If no overlap detected then returns -1
// 
int										// loci at which there is a feature overlap
DistNearFeat(CBEDfile *pBEDFile,		// contains features to be checked against
			 int ChromID,				// nearest feature will be on this chromosome
			 int StartLoci,				// start loci on chrom to check for overlaps
			 int EndLoci)				// end loci on chrom to check
{
int ExclLoci;
ExclLoci = StartLoci;
while(ExclLoci >= 0)
	{
	if(pBEDFile->InAnyFeature(ChromID,ExclLoci,ExclLoci))
		return(ExclLoci);
	if(EndLoci < StartLoci)
		{
		ExclLoci -= 1;
		if(ExclLoci < EndLoci)
			break;
		}
	else
		{
		ExclLoci += 1;
		if(ExclLoci > EndLoci)
			break;
		}
	}
return(-1);
}

int
ProcessFeatures(etReadProfile ReadProfile,		// profile for: 0 read density, 1: read starts, 2:  unique read starts (default is for read density)
				int PkWinLen,					// window length over which to calculate the peak counts
				char StrandOI,					// strand of interest
				etPFLFeature Feature,			// center profiles around this functional feature
				int InterGenicLen,				// intergenic flank length
				int IntraGenicLen,				// intragenic flank length
				int NumBins,					// profile into this many bins upstream and downstream
				etPOPpolicy ProfPol,			// profile overlap policy
				int KDist,						// overlap distance
				char *pszExcludeFile,			// reads exclude features
				char *pszFeatureFile,			// file containing annotation features of interest to user
				char *pszRsltsFile)				// output CSV profiles file
{
int Rslt;
int NumExons;
int FiltCntLoci;
int StartOfs;
int ProfStartLoci;
int ProfEndLoci;

int CurFeatureID;
int StartLoci;
int EndLoci;
int LociOfInterest;
int Score;
int ExclChromID;
int ChromID;
char szChrom[128];
char szPrevChrom[128];
char szFeatName[128];
char Strand;
int RefID;
int IntergenicStart;

char szOutBuff[cMaxNumBins*8];
int BuffIdx;
int BinIdx;

int NumEls;

CBEDfile *pFiltBed;

if((m_pBEDFile = new CBEDfile)==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CBEDfile");
	Reset();
	return(eBSFerrObj);
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Starting to load features from %s",pszFeatureFile);
if((Rslt=m_pBEDFile->Open(pszFeatureFile,eBTAnyBed)) !=eBSFSuccess)
	{
	while(m_pBEDFile->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,m_pBEDFile->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open '%s' for processing",pszFeatureFile);
	Reset();
	return(eBSFerrOpnFile);
	}

pFiltBed = m_pBEDFile;

if(pszExcludeFile[0] != '\0')
	{
	if((m_pExclBEDFile = new CBEDfile)==NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CBEDfile");
		Reset();
		return(eBSFerrObj);
		}

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Starting to load exclusion features from %s",pszExcludeFile);
	if((Rslt=m_pExclBEDFile->Open(pszExcludeFile,eBTAnyBed)) !=eBSFSuccess)
		{
		while(m_pExclBEDFile->NumErrMsgs())
			gDiagnostics.DiagOut(eDLFatal,gszProcName,m_pExclBEDFile->GetErrMsg());
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open '%s' for processing",pszExcludeFile);
		Reset();
		return(eBSFerrOpnFile);
		}
	pFiltBed = m_pExclBEDFile;
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Load completed, iterating features...");

#ifdef _WIN32
if((m_hRsltsFile = open(pszRsltsFile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
if((m_hRsltsFile = open(pszRsltsFile, O_RDWR | O_CREAT | O_TRUNC, S_IREAD | S_IWRITE ))==-1)
#endif
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create results profile file: %s - %s",pszRsltsFile,strerror(errno));
	Reset();
	return(eBSFerrCreateFile);
	}	

BuffIdx = sprintf(szOutBuff,"\"Name\",\"Chrom\",\"Strand\",\"NumExons\",\"DNAlen\",\"RNALen\",\"Loci\",\"TotCnts\",\"PeakWinCnts\",\"TotBinCnts\"");
for(BinIdx = 0; BinIdx < NumBins; BinIdx++)
	BuffIdx += sprintf(&szOutBuff[BuffIdx],",\"Bin%d\"",BinIdx+1);
BuffIdx += sprintf(&szOutBuff[BuffIdx],"\n");
CUtility::SafeWrite(m_hRsltsFile,szOutBuff,BuffIdx);

szPrevChrom[0] = '\0';
CurFeatureID = 0;
ChromID = 0;
ExclChromID = 0;
IntergenicStart = 0;
RefID = 0;
Rslt = eBSFSuccess;
NumEls = 0;
while(Rslt == eBSFSuccess && (CurFeatureID = m_pBEDFile->GetNextFeatureID(CurFeatureID)) > 0)
	{
	NumEls += 1;

	m_pBEDFile->GetFeature(CurFeatureID,	// feature instance identifier
				szFeatName,				// where to return feature name
				szChrom,				// where to return chromosome name
				&StartLoci,				// where to return feature start on chromosome (0..n) 
				&EndLoci,				// where to return feature end on chromosome
 				&Score,					// where to return score
 				&Strand);				// where to return strand

	if(CurFeatureID == 1 || stricmp(szChrom,szPrevChrom))	// if new chromosome then reset IntergenicStart
		{
		strcpy(szPrevChrom,szChrom);
		ChromID = m_pBEDFile->LocateChromIDbyName(szChrom);
		if(m_pExclBEDFile != '\0')
			ExclChromID = m_pExclBEDFile->LocateChromIDbyName(szChrom);
		IntergenicStart = 0;
		}

	if(StrandOI != '*')
		{
		if(StrandOI != Strand)
			continue;
		}

	if(Feature == ePFLIsoforms)
		{
		ReportIsoform(ReadProfile,PkWinLen,CurFeatureID,InterGenicLen,IntraGenicLen,NumBins);
		continue;
		}

	NumExons = m_pBEDFile->GetNumExons(CurFeatureID);
	StartOfs = 0;
	switch(Feature) {
		case ePFLTSS:			// only process TSS flanking regions
			if(Strand == '+')
				{
				LociOfInterest = StartLoci;
				ProfStartLoci = StartLoci - InterGenicLen;
				ProfEndLoci = StartLoci + IntraGenicLen - 1;
				if(ProfStartLoci < 0)
					{
					StartOfs = abs(ProfStartLoci);
					ProfStartLoci = 0;
					}
				if(ProfPol >= ePOPkdist)
					{
					if((FiltCntLoci = DistNearFeat(pFiltBed,ChromID,LociOfInterest-1,ProfStartLoci - KDist))!=-1)
						{
						if(FiltCntLoci + KDist >= LociOfInterest)
							StartOfs = -1;
						else
							{
							StartOfs = KDist + FiltCntLoci - ProfStartLoci;
							ProfStartLoci = FiltCntLoci + KDist;
							}
						}
					}
				}	
			else				// on minus strand
				{
				LociOfInterest = EndLoci;
				ProfStartLoci = EndLoci + InterGenicLen;
				ProfEndLoci = EndLoci - IntraGenicLen + 1;
				if(ProfEndLoci < 0)
					{
					StartOfs = abs(ProfEndLoci);
					ProfEndLoci = 0; 
					}
				if(ProfPol >= ePOPkdist)
					{
					if((FiltCntLoci = DistNearFeat(pFiltBed,ChromID,LociOfInterest+1,ProfStartLoci + KDist))!=-1)
						{
						if(FiltCntLoci - KDist <= LociOfInterest)
							StartOfs = -1;
						else
							{
							StartOfs = ProfStartLoci + KDist - FiltCntLoci;
							ProfStartLoci = FiltCntLoci - KDist;
							}
						}
					}
				}
			break;
		case ePFLTES:			// only process TES flanking regions
			if(Strand == '+')
				{
				LociOfInterest = EndLoci;
				ProfStartLoci = EndLoci - IntraGenicLen + 1;
				ProfEndLoci = LociOfInterest + InterGenicLen;
				if(ProfStartLoci < 0)
					{
					StartOfs = abs(ProfStartLoci);
					ProfStartLoci = 0;
					}
				if(ProfPol >= ePOPkdist)
					{
					if((FiltCntLoci = DistNearFeat(pFiltBed,ChromID,LociOfInterest-1,ProfStartLoci - KDist))!=-1)
						{
						if(FiltCntLoci + KDist >= LociOfInterest)
							StartOfs = -1;
						else
							{
							StartOfs = KDist + FiltCntLoci - ProfStartLoci;
							ProfStartLoci = FiltCntLoci + KDist;
							}
						}
					}
				}
			else				// on minus strand
				{
				LociOfInterest = StartLoci;
				ProfStartLoci = StartLoci + IntraGenicLen - 1;
				ProfEndLoci = StartLoci - InterGenicLen;
				if(ProfEndLoci < 0)
					{
					StartOfs = abs(ProfEndLoci);
					ProfEndLoci = 0;
					}
				if(ProfPol >= ePOPkdist)
					{
					if((FiltCntLoci = DistNearFeat(pFiltBed,ChromID,LociOfInterest+1,ProfStartLoci + KDist))!=-1)
						{
						if(FiltCntLoci - KDist <= LociOfInterest)
							StartOfs = -1;
						else
							{
							StartOfs = ProfStartLoci + KDist - FiltCntLoci;
							ProfStartLoci = FiltCntLoci - KDist;
							}
						}
					}
				}
			break;
		}

	ReportCenteredAt(ReadProfile,PkWinLen,szFeatName,NumExons,Strand,szChrom,LociOfInterest,NumBins,InterGenicLen,IntraGenicLen,StartOfs,ProfStartLoci,ProfEndLoci);
	}

m_pBEDFile->Close();
if(m_hRsltsFile != -1)
	{
	close(m_hRsltsFile);
	m_hRsltsFile = -1;
	}
return(0);
}

int
Process(etPMode PMode,					// processing mode
		etFMode FMode,					// output format - CSV or BED
		etReadProfile ReadProfile,		// profile for: 0 read density, 1: read starts, 2:  unique read starts (default is for read density)
		int PkWinLen,					// window length over which to calculate the peak counts
		char Strand,					// ROI are on this strand
		etPFLFeature Feature,			// center profiles around this functional feature
		int InterGenicLen,					// 5' upstream flank length
		int IntraGenicLen,					// 3' upstream flank length
		int NumBins,					// profile into this many bins upstream and downstream
		etPOPpolicy ProfPol,			// profile overlap policy
		int KDist,						// overlap distance
		char *pszExcludeFile,			// reads exclude features
		char *pszFeatureFile,			// file containing annotation features of interest to user
		int NumInputFiles,				// number of input files
		char **ppszInFileSpecs,			// input files
		char *pszRsltsFile)				// output CSV profilefile
{
int Rslt;
char *pszInFile;
int FileIdx;

Init();

etClassifyFileType FileType;

CSimpleGlob glob(SG_GLOB_FULLSORT);

for(FileIdx = 0; FileIdx < NumInputFiles; FileIdx++)
	{
	glob.Init();

	if(glob.Add(ppszInFileSpecs[FileIdx]) < SG_SUCCESS)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to glob '%s",ppszInFileSpecs[FileIdx]);
		Reset();
		return(eBSFerrOpnFile);
		}

	if(glob.FileCount() <= 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to locate any source raw reads file matching '%s",ppszInFileSpecs[FileIdx]);
		continue;
		}

	for (int FileID = 0; FileID < glob.FileCount(); FileID += 1)
		{
		pszInFile = glob.File(FileID);

		if(PMode == 0)
			FileType = CUtility::ClassifyFileType(pszInFile);
		else
			FileType = (etClassifyFileType)(PMode - 1);

		switch(FileType) {
			case eCFTopenerr:		// unable to open file for reading
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open file: '%s'",pszInFile);
				return(eBSFerrOpnFile);

			case eCFTlenerr:		// file length is insufficent to classify type
				gDiagnostics.DiagOut(eDLInfo,gszProcName,"Unable to classify file type (insuffient data points): '%s'",pszInFile);
				return(eBSFerrFileAccess);

			case eCFTunknown:		// unable to reliably classify
				gDiagnostics.DiagOut(eDLInfo,gszProcName,"Unable to reliably classify file type: '%s'",pszInFile);
				return(eBSFerrFileType);

			case eCFTCSV:			// file has been classified as being CSV
				gDiagnostics.DiagOut(eDLInfo,gszProcName,"Parsing CSV file: '%s'",pszInFile);
				Rslt = GenCSVdensity(ReadProfile,Strand,pszInFile);
				break;

			case eCFTBED:			// file has been classified as being BED
				gDiagnostics.DiagOut(eDLInfo,gszProcName,"Parsing BED file: '%s'",pszInFile);
				Rslt = GenBEDdensity(ReadProfile,Strand,pszInFile);
				break;

			case eCFTSAM:			// file has been classified as being SAM
				gDiagnostics.DiagOut(eDLInfo,gszProcName,"Parsing SAM file: '%s'",pszInFile);
				Rslt = GenSAMdensity(ReadProfile,Strand,pszInFile);
				break;
			}
		if(Rslt != eBSFSuccess)
			{
			Reset();
			return(Rslt);
			}
		}
	}

// now generate the histogram bin profiles for features
if(Rslt == eBSFSuccess)
	Rslt = ProcessFeatures(ReadProfile,PkWinLen,Strand,Feature,InterGenicLen,IntraGenicLen,NumBins,ProfPol,KDist,pszExcludeFile,pszFeatureFile,pszRsltsFile);

Reset();
return(Rslt);
}




