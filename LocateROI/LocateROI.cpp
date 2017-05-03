// LocateROI.cpp : Defines the entry point for the console application
// Locate regions of interest from short read alignments
// Locates regions which have a median read coverage within a specified range, region lengths are in a specified range
// and any inter-read gap is less than a maximum specified

// Version 1.2.0     Added SAM loci file input and autodetermination of input loci file format
// Version 1.2.1     Added SAM/BAM loci file input and autodetermination of input loci file format

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

const char *cpszProgVer = "1.3.2";			// increment with each release

const int cMaxNumAssocFiles = 20;			// alllow for at most this many feature association files 

const int cChromSeqReAlloc = 5000000;	// realloc chrom sequence size
const UINT32 cMaxAccumCnt= 0x07ffffff;  // truncate accumulated counts to be at most this value
const UINT32 cRetainCnt  = 0x08000000;  // accumulated counts >= this value are counts marked as being to be retained

const int cMinMedianCov  = 1;			// minimum user specified region median coverage
const int cDfltMedianCov = 2;			// default minimum region median coverage
const int cMaxMedianCov = 1000;			// maximum user specified minimum median coverage
const int cMinRegionLen  = 50;			// minimum user specified region length
const int cDfltRegionLen = 100;			// default minimum region length
const int cMaxRegionLen = 100000;		// maximum user specified minimum region length

const int cMinGapLen = 0;				// minimum user specified gap without any aligned reads in reported regions
const int cDfltGapLen = 10;				// default user specified gap without any aligned reads in reported regions
const int cMaxGapLen = 100;				// maximum user specified gap without any aligned reads in reported regions

const int cDfltMaxDist = 1000000000;	// used to flag that any feature is at least this distance away
const char *cpszNoFeat = "FxNxL";		// unlocated features are assigned this name

// processing modes
typedef enum TAG_ePMode {		
	ePMdefault,					// Auto determine
	ePMCsv,						// CSV loci
	ePMBed,						// UCSC BED format
	ePMSam,						// SAM format
	ePMplaceholder				// used to set the enumeration range
	} etPMode;

// strand processing modes
typedef enum TAG_eStrandProc {
		eStrandDflt,			// default is to ignore the strand
		eStrandWatson,			// process for Watson
		eStrandCrick,			// process for Crick
		eStrandPlaceholder
} etStrandProc;

// output format modes
typedef enum TAG_eFMode {
	eFMsumCSV,					// default is for summary CSV in which the minimum distance from all BEDs is reported 
	eFMallCSV,					// CSV in which the minimum distance for each BED is reported 
	eFMbed,						// UCSC BED format
	eFMplaceholder				// used to set the enumeration range
	} etFMode;

typedef enum eBEDRegion {
	eMEGRAny = 0,				// process any region
	eMEGRIntergenic,			// only process intergenic
	eMEGRExons,					// only process exons
	eMEGRIntrons,				// only process introns
	eMEGRCDS,					// only process CDSs
	eMEGUTR,					// only process UTRs
	eMEG5UTR,					// only process 5'UTRs
	eMEG3UTR					// only process 3'UTRs
} etBEDRegion;

#pragma pack(1)

int m_NumFeatFiles;						// number of feature file names in m_pszFeatFiles[] following
char *m_pszFeatFiles[cMaxNumAssocFiles];	// feature file names

typedef struct TAG_sROI {
	int RegionID;					// unique region identifier
	int ChromIdx;					// ROI is on this chromosome (m_ChromCnts[])
	int StartOfRegion;				// region starts at this loci
	int EndOfRegion;				// region ends at this loci
	char Strand;					// and is on this strand
	double BPKM;					// bases per thousand per Million aligned bases

	char USFeatStrand;				// nearest upstream feature is on this strand
	int USFeatDist;					// nearest upstream feature distance
	char szUSFeatID[cMaxGeneNameLen];	// nearest upstream feature name
	int USFeatFileID;				// nearest upstream feature file (in m_pszFeatFiles[USFeatFileID-1])

	char DSFeatStrand;				// nearest downstream feature is on this strand
	int DSFeatDist;					// nearest downstream feature distance
	char szDSFeatID[cMaxGeneNameLen];	// nearest downstream feature name
	int DSFeatFileID;				// nearest downstream feature file (in m_pszFeatFiles[DSFeatFileID-1])
} tsROI;
#pragma pack()

const size_t cROIAlloc = 10000;		// incrementally allocate for this many ROIs


int
Process(etPMode PMode,					// processing mode
		etFMode FMode,					// output format - CSV or BED
		char *pszTitle,					// CSV species or title used if BED output format
		char ReadStrand,				// only accept reads from this strand ('*' if either)
		char RestainStrand,				// filter ROI by filter elements on this strand ('*' if either)
		char FeatDistStrand,			// distances to features on this strand ('*' if either)
		etBEDRegion Region,				// functional region of interest
		int Limit,						// limit (0 if no limit) processing to this many reads total
		int MinMedianCov,				// minimum median coverage required in reported regions
		int MinRegionLen,				// report regions which are of at least this length
		int MaxGapLen,					// report regions containing no read gaps of <= this length 
		char *pszRetainFile,			// file containing regions to be retained which of interest to user
		int NumAssocFiles,				// number of association files in pszAssocFiles[]
		char *pszAssocFiles[],			// distance association files, if prefixed with #n then n is the feature tye
		char *pszInFile,				// input CSV or BED file containing read alignment loci
		char *pszRsltsFile);			// output CSV region loci file

bool								// returns false if errors
DistanceToFeatures(int AnnoFileID,				// annotated file identifer
				   char FiltStrand,				// only interested in distances to nearest feature which is on this strand
				   int ROIChromID,				// ROI is on this chromosome
				   char *pszROIChrom,			// ROI chrom name
				    tsROI *pROI);


char *Region2Txt(etBEDRegion Region);

CStopWatch gStopWatch;
CDiagnostics gDiagnostics;				// for writing diagnostics messages to log file
char gszProcName[_MAX_FNAME];			// process name

int m_hRsltsFile;						// handle for opened results file

char m_ReadStrand;						// accept reads on this strand
char m_RestainStrand;					// filter ROI by filter elements on this strand 
char m_FeatDistStrand;					// distances to features on this strand 

CCSVFile *m_pCSVFile;					// used if processing input CSV file for coverage
CBEDfile *m_pBEDFile;					// used if processing input BED files for coverage

CBEDfile *m_pDistBEDFile;				// current distance annotation bed file

const int cMaxChromCov = 1000;					// can handle at most this many chromosomes
const int cAllocCovCnts = 0x07fffff;			// allocate for chrom coverage cnts in this sized increments

typedef struct TAG_sChromCnts {
	char szChrom[cMaxDatasetSpeciesChrom+1];	// coverage is on this chromosome
	int AllocCovCnts;							// allocated to hold cnts for this sized chromosome
	int StartOfs;								// pCovCnts[offset] of first coverage cnt
	int EndOfs;									// pCovCnts[offset] of last coverage cnt
	UINT32 *pCovCnts;							// coverage counts
} tsChromCnts;



int m_NumOfROIs;					// current total number of ROIs in m_pROIs
int m_AllocNumROIs;					// m_pROIs has been allocated to hold at most this number of ROIs
size_t	m_AllocNumROIsMem;			// actual memory size allocated
tsROI *m_pROIs;						// allocated to hold all ROIs
UINT32 m_NumAcceptedReads;			// total number of accepted reads or element loci

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
	return _T("LocateROL");
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
etFMode FMode;				// format output mode
char szTitle[cMaxDatasetSpeciesChrom];	// species if output to CSV or track title if output format is UCSC BED

int ReadStrand;					// accept reads which are on this strand
int RestainStrand;				// ROI filter strand
int FeatDistStrand;				// feature distance strand

etBEDRegion Region;				// process for this functional region only

int MinMedianCov;				// minimum median coverage required in reported regions
int MinRegionLen;				// report regions which are of at least this length
int MaxGapLen;					// report regions containing no read gaps of <= this length 

char szRsltsFile[_MAX_PATH];	// output results file
char szInFile[_MAX_PATH];		// input file 
char szRegionFile[_MAX_PATH];	// annotated regions to retain filter file
int  NumAssocFiles;				// number of association input files
char *pszAssocFiles[cMaxNumAssocFiles];  // input (wildcards allowed) association BED files

int Limit;						// 0 if no limit, otherwise only process for at most this number of bases

// command line args
struct arg_lit  *help    = arg_lit0("hH","help",                "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"diagnostics log file");

struct arg_int *pmode = arg_int0("m","mode","<int>",		    "input reads loci file: 0 - Auto determine, 1 - CSV, 2 - BED, 3 - SAM or BAM (default = 0)");
struct arg_int *format = arg_int0("M","format","<int>",		    "output format: 0 - CSV min ROI distances of all associations, 1 - CSV min ROI distances for each association BED, 2 - UCSC BED (default: 0)");
struct arg_str  *title = arg_str1("t","title","<string>",       "species if CSV output, track title if BED output");

struct arg_int  *readstrand = arg_int0("s","readstrand","<int>", "Accept reads on this strand: 0 - either, 1 - Watson '+', 2 - Crick '-' (default is either)");
struct arg_int  *retainstrand = arg_int0("S","retainstrand","<int>", "Retain ROI filter file by this strand: 0 - either, 1 - Watson '+', 2 - Crick '-' (default is either)");
struct arg_int  *featdiststrand = arg_int0("d","featdiststrand","<int>", "Report nearest feature distances on this strand: 0 - either, 1 - Watson '+', 2 - Crick '-' (default is either)");

struct arg_int *minmediancov = arg_int0("c","minmediancov","<int>","reported regions must have at least this median coverage (default = 2");
struct arg_int *minregionlen = arg_int0("l","minregionlen","<int>","reported regions must be of at least this length (default = 100");
struct arg_int *maxgaplen = arg_int0("g","maxgaplen","<int>",   "reported regions must not contain gaps longer than this length (default = 10nt");
struct arg_int *region = arg_int0("r","genomicregion","<int>",	"Retain annotated regions 0:ALL,1:Intergenic,2:Exons,3:Introns,4:CDS,5:UTRs,6:5'UTR,7:3'UTR (default = ALL)");
struct arg_file *infile = arg_file1("i","in","<file>",			"input from this CSV, BED or SAM/BAM reads loci file");
struct arg_file *filterfile = arg_file0("I","filter","<file>",	"retain read base loci by regions in this BED file");
struct arg_file *assocfiles = arg_filen("a","assoc","<file>",0, cMaxNumAssocFiles,"optionally associate ROIs to nearest features in these annotated feature BED files");
struct arg_file *outfile = arg_file1("o","out","<file>",		"output regions to this file");
struct arg_int *limit = arg_int0("L","limit","<int>",		    "limit number of reads processed whilst debugging (0 == no limit)");


struct arg_end *end = arg_end(20);

void *argtable[] = {help,version,FileLogLevel,LogFile,
					pmode,format,title,readstrand,retainstrand,featdiststrand,region,minmediancov,minregionlen,maxgaplen,infile,filterfile,assocfiles,outfile,limit,
					end};

char **pAllArgs;
int argerrors;
argerrors = CUtility::arg_parsefromfile(argc,(char **)argv,&pAllArgs);
if(argerrors >= 0)
	argerrors = arg_parse(argerrors,pAllArgs,argtable);

    /* special case: '--help' takes precedence over error reporting */
if (help->count > 0)
        {
		printf("\n%s Locate Regions of Interest, Version %s\nOptions ---\n", gszProcName,cpszProgVer);
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

	PMode = (etPMode)(pmode->count ? pmode->ival[0] : 0);
	if(PMode < 0 || PMode >= ePMplaceholder)
		{
		printf("\nError: Processing mode '-m%d' specified outside of range %d..%d",PMode,0,(int)ePMplaceholder-1);
		exit(1);
		}

	FMode = (etFMode)(format->count ? format->ival[0] : eFMsumCSV);
	if(FMode < eFMsumCSV || FMode >= eFMplaceholder)
		{
		printf("\nError: Output format mode '-m%d' specified outside of range %d..%d",FMode,eFMsumCSV,(int)eFMplaceholder-1);
		exit(1);
		}

	strncpy(szTitle,title->sval[0],sizeof(szTitle));
	szTitle[sizeof(szTitle)-1] = '\0';
	CUtility::TrimQuotedWhitespcExtd(szTitle);
	CUtility::ReduceWhitespace(szTitle);

	ReadStrand = readstrand->count ? readstrand->ival[0] : 0;
	if(ReadStrand < eStrandDflt || ReadStrand >= eStrandPlaceholder)
		{
		printf("\nError: Read strand '-s%d' specified outside of range 0..2",ReadStrand);
		exit(1);
		}
	switch(ReadStrand) {
		case 1: ReadStrand = (int)'+'; break;
		case 2: ReadStrand = (int)'-'; break;
		case 0: ReadStrand = (int)'*'; break;
		}

	RestainStrand = retainstrand->count ? retainstrand->ival[0] : 0;
	if(RestainStrand < 0 || RestainStrand > 2)
		{
		printf("\nError: Retain strand '-S%d' specified outside of range 0..2",RestainStrand);
		exit(1);
		}
	switch(RestainStrand) {
		case 1: RestainStrand = (int)'+'; break;
		case 2: RestainStrand = (int)'-'; break;
		case 0: RestainStrand = (int)'*'; break;
		}

	FeatDistStrand = featdiststrand->count ? featdiststrand->ival[0] : 0;
	if(FeatDistStrand < 0 || FeatDistStrand > 2)
		{
		printf("\nError: Filter strand '-S%d' specified outside of range 0..2",FeatDistStrand);
		exit(1);
		}
	switch(FeatDistStrand) {
		case 1: FeatDistStrand = (int)'+'; break;
		case 2: FeatDistStrand = (int)'-'; break;
		case 0: FeatDistStrand = (int)'*'; break;
		}


	Region = (etBEDRegion)(region->count ? region->ival[0] : eMEGRAny);	// default as being any region
	if(Region < eMEGRAny || Region > eMEG3UTR)
		{
		printf("\nSpecified region '-r%d' outside of range 0..%d",Region,eMEG3UTR);
		exit(1);
		}
	if(!filterfile->count && Region != eMEGRAny)
		{
		printf("\nError: Specified region '-r%d' (%s) but no annotated region filter file specified with '-I<file>'",Region,Region2Txt((etBEDRegion)Region));
		exit(1);
		}

	Limit = limit->count ? limit->ival[0] : 0;
	if(Limit < 0)
		Limit = 0;

	MinMedianCov = minmediancov->count ? minmediancov->ival[0] : cDfltMedianCov;
	if(MinMedianCov < cMinMedianCov || MinMedianCov > cMaxMedianCov)
		{
		printf("\nError: Minimum median coverage '-c%d' specified outside of range %d..%d",MinMedianCov,1,cMaxMedianCov);
		exit(1);
		}

	MinRegionLen = minregionlen->count ? minregionlen->ival[0] : cDfltRegionLen;
	if(MinRegionLen < cMinRegionLen || MinRegionLen > cMaxRegionLen)
		{
		printf("\nError: Minimum region length '-l%d' specified outside of range %d..%d",MinRegionLen,cMinRegionLen,cMaxRegionLen);
		exit(1);
		}

	MaxGapLen = maxgaplen->count ? maxgaplen->ival[0] : cDfltGapLen;
	if(MaxGapLen < cMinGapLen || MaxGapLen > cMaxGapLen)
		{
		printf("\nError: Maximum gap length '-g%d' specified outside of range %d..%d",MaxGapLen,cMinGapLen,cMaxGapLen);
		exit(1);
		}

	strcpy(szInFile,infile->filename[0]);
	strcpy(szRsltsFile,outfile->filename[0]);
	if(filterfile->count)
		strcpy(szRegionFile,filterfile->filename[0]);
	else
		szRegionFile[0] = '\0';

	if(!assocfiles->count)
		{
		if(szRegionFile[0] != '\0')
			{
			printf("\nWarning: No features distance association files specified - defaulting to using %s as the  association file",szRegionFile);
			NumAssocFiles = 1;
			pszAssocFiles[0] = new char [_MAX_PATH];
			strcpy(pszAssocFiles[0],szRegionFile);
			}
		else
			{
			NumAssocFiles = 0;
			pszAssocFiles[0] = NULL;
			}
		}
	else
		{
		for(NumAssocFiles=Idx=0;NumAssocFiles < cMaxNumAssocFiles && Idx < assocfiles->count; Idx++)
			{
			pszAssocFiles[Idx] = NULL;
			if(pszAssocFiles[NumAssocFiles] == NULL)
				pszAssocFiles[NumAssocFiles] = new char [_MAX_PATH];
			strncpy(pszAssocFiles[NumAssocFiles],assocfiles->filename[Idx],_MAX_PATH);
			pszAssocFiles[NumAssocFiles][_MAX_PATH-1] = '\0';
			CUtility::TrimQuotedWhitespcExtd(pszAssocFiles[NumAssocFiles]);
			if(pszAssocFiles[NumAssocFiles][0] != '\0')
				NumAssocFiles++;
			}

		if(!NumAssocFiles)
			{
			printf("\nError: After removal of whitespace, no distance association input file(s) specified with '-a<filespec>' option)");
			exit(1);
			}
		}

			// now that command parameters have been parsed then initialise diagnostics log system
	if(!gDiagnostics.Open(szLogFile,(etDiagLevel)iScreenLogLevel,(etDiagLevel)iFileLogLevel,true))
		{
		printf("\nError: Unable to start diagnostics subsystem.");
		if(szLogFile[0] != '\0')
			printf(" Most likely cause is that logfile '%s' can't be opened/created",szLogFile);
		exit(1);
		}

// show user current resource limits
#ifndef _WIN32
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Resources: %s",CUtility::ReportResourceLimits());
#endif

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Version: %s Processing parameters:",cpszProgVer);
	const char *pszDescr;
	switch(PMode) {
		case ePMdefault:
			pszDescr = "Auto determine";
			break;
		case ePMCsv:
			pszDescr = "CSV loci";
			break;
		case ePMBed:			// UCSC BED format
			pszDescr = "UCSC BED loci";
			break;
		case ePMSam:
			pszDescr = "SAM";
			break;
		}
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"processing mode is : '%s' input file format",pszDescr);

	switch(FMode) {
		case eFMsumCSV:
			pszDescr = "CSV minimum ROI distances from all association BED files";
			break;
		case eFMallCSV:
			pszDescr = "CSV minimum ROI distance for each association BED files";
			break;
		case eFMbed:
			pszDescr = "UCSC BED";
			break;
		}
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"output format is : '%s'",pszDescr);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"title: '%s'",szTitle);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"accept reads on this strand only: '%c'",(char)ReadStrand);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"retain by features on this strand only: '%c'",(char)RestainStrand);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"report distances to features on this strand only: '%c'",(char)FeatDistStrand);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"report regions of at least this length: %d",MinRegionLen);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"report regions with at least this coverage (median): %d",MinMedianCov);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"report regions with no gaps longer than: %d",MaxGapLen);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"input %s file: '%s'",pszDescr, szInFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Retain Region: %s",Region2Txt((etBEDRegion)Region));

	if(szRegionFile[0])
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Retain regions from BED file: '%s'",szRegionFile);
	if(NumAssocFiles)
		for(Idx=0;Idx < NumAssocFiles; Idx++)
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"associate to features in this input BED file (%d): '%s'",Idx+1,pszAssocFiles[Idx]);		

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"output Regions Of Interest to file: '%s'",szRsltsFile);
	if(Limit > 0)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"limit processing to first %d reads",Limit);

	#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	gStopWatch.Start();
	Rslt = Process(PMode,FMode,szTitle,(char)ReadStrand,(char)RestainStrand,(char)FeatDistStrand,Region,Limit,MinMedianCov,MinRegionLen,MaxGapLen,szRegionFile,NumAssocFiles,pszAssocFiles,szInFile,szRsltsFile);

	gStopWatch.Stop();
	Rslt = Rslt >=0 ? 0 : 1;
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Exit code: %d Total processing time: %s",Rslt,gStopWatch.Read());
	exit(Rslt);
	}
else
	{
	printf("\n%s Locate Regions of Interest, Version %s\n", gszProcName,cpszProgVer);
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
		return((char *)"All except Intergenic");

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

void
Init(void)
{
m_pCSVFile = NULL;
m_pBEDFile = NULL;
m_pDistBEDFile = NULL;

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
m_NumFeatFiles = 0;
m_AllocNumROIsMem = 0;
m_AllocNumROIs = 0;
m_NumOfROIs = 0;
m_NumAcceptedReads = 0;
m_pROIs = NULL;
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
if(m_pDistBEDFile != NULL)
	{
	delete m_pDistBEDFile;
	m_pDistBEDFile = NULL;
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

if(m_pROIs != NULL)
	{
#ifdef _WIN32
	free(m_pROIs);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pROIs != MAP_FAILED)
		munmap(m_pROIs,m_AllocNumROIsMem);
#endif
	m_pROIs = NULL;
	}

if(m_NumFeatFiles > 0)
	{
	for(int Idx = 0; Idx < m_NumFeatFiles; Idx++)
		{
		if(m_pszFeatFiles[Idx] != NULL)
			{
			delete m_pszFeatFiles[Idx];
			m_pszFeatFiles[Idx] = NULL;
			}
		}
	m_NumFeatFiles = 0;
	}

m_NumAcceptedReads = 0;
}

int
BuildReadCoverage(char *pszChrom,		// coverage is onto this chrom
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
while(StartOfs++ <= EndOfs)
	{
	// clamp accumulated cnts to be no greater than cMaxAccumCnt
	if((cMaxAccumCnt - *pCovCnts) > (UINT32)Cnt)
		*pCovCnts += (UINT32)Cnt;
	else
		*pCovCnts = cMaxAccumCnt;
	pCovCnts += 1;
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
IdentROI(etFMode FMode,					// output in this format to m_hRsltsFile
		  int MinMedianCov,				// minimum median coverage required in reported regions
			int MinRegionLen,			// report regions which are of at least this length
			int MaxGapLen,				// report regions containing no read gaps of <= this length 
			bool bNoFilter)				// true if all counts to be processed
{
tsROI *pROI;
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
UINT64 TotGenomeCnts;
double GenomeCntsPerM;
size_t memreq;

if(m_pROIs == NULL)		// should be the case but who knows :-)
	{
		// initial allocation of memory to hold regions of interest
	memreq = cROIAlloc * 2 * sizeof(tsROI);

#ifdef _WIN32
	m_pROIs = (tsROI *) malloc(memreq);	// initial allocation
	if(m_pROIs == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"IdentROI: Memory allocation of %lld bytes - %s",(INT64)memreq,strerror(errno));
		Reset();
		return(eBSFerrMem);
		}
#else
	// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
	m_pROIs = (tsROI *)mmap(NULL,memreq, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
	if(m_pROIs == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"IdentROI: Memory allocation of %lld bytes through mmap()  failed - %s",(INT64)memreq,strerror(errno));
		m_pROIs = NULL;
		Reset();
		return(eBSFerrMem);
		}
#endif
	m_AllocNumROIsMem = memreq;
	m_NumOfROIs = 0;
	m_AllocNumROIs = cROIAlloc * 2;
	}

TotGenomeCnts = 0;
pChrom = &m_ChromCnts[0];
for(ChromIdx = 0 ; ChromIdx < m_NumChromsCov; ChromIdx++,pChrom++)
	{
	pCnts = &pChrom->pCovCnts[pChrom->StartOfs];
	for(SeqIdx = pChrom->StartOfs; SeqIdx <= pChrom->EndOfs+1; SeqIdx++,pCnts++)
		TotGenomeCnts += *pCnts & cMaxAccumCnt;
	}
GenomeCntsPerM = 1000000.0 / TotGenomeCnts;

pChrom = &m_ChromCnts[0];
NumRegions = 0;
pROI = &m_pROIs[m_NumOfROIs];
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
						// check though that m_pROIs can hold this ROI
						if((m_NumOfROIs + 10) >= m_AllocNumROIs)
							{
							tsROI *pTmpAlloc;
							memreq = m_AllocNumROIsMem + (cROIAlloc * sizeof(tsROI));
#ifdef _WIN32
							pTmpAlloc = (tsROI *) realloc(m_pROIs,memreq);
#else
							pTmpAlloc = (tsROI *)mremap(m_pROIs,m_AllocNumROIsMem,memreq,MREMAP_MAYMOVE);
							if(pTmpAlloc == MAP_FAILED)
								pTmpAlloc = NULL;
#endif
							if(pTmpAlloc == NULL)
								{
								gDiagnostics.DiagOut(eDLFatal,gszProcName,"IdentROI: Memory re-allocation to %d bytes - %s",memreq,strerror(errno));
								Reset();	
								return(eBSFerrMem);
								}
							m_pROIs = pTmpAlloc;
							m_AllocNumROIsMem = memreq;
							m_AllocNumROIs += cROIAlloc;
							pROI = &m_pROIs[m_NumOfROIs];
							}
						NumRegions += 1;
						pROI->RegionID = NumRegions;
						pROI->ChromIdx = ChromIdx+1;
						pROI->StartOfRegion = StartOfRegion;
						pROI->EndOfRegion = EndOfRegion;
						pROI->Strand = m_ReadStrand == '*' ? '+' : m_ReadStrand;
						pROI->BPKM = (SumCnts * 1000.0 / TotRegionLen) * GenomeCntsPerM;
						pROI->szDSFeatID[0] = '\0';
						pROI->DSFeatDist = -1;
						pROI->DSFeatStrand = '*';
						pROI->szUSFeatID[0] = '\0';
						pROI->USFeatDist = -1;
						pROI->USFeatStrand = '*';
						m_NumOfROIs += 1;
						pROI += 1;
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
return(NumRegions);
}

int
WriteRegions(bool bFinal,					// true if last call to this function
			 char *pszAnnoFile,				// file containing annotated features for distance calculations
			 int AnnoFileID,				// annotated file identifer
			 etFMode FMode,					// output in this format to m_hRsltsFile
			char *pszTitle,					// CSV species or title used if BED output format
			 char FeatDistStrand)			// distance to features on this strand 
{
static bool bFirst = true;
int Rslt;
int BuffIdx;
char szLineBuff[8096];
tsChromCnts *pChrom;
int ChromID;
int ROIidx;
tsROI *pROI;

if(m_pDistBEDFile != NULL)
	{
	delete m_pDistBEDFile;
	m_pDistBEDFile = NULL;
	}

if(pszAnnoFile != NULL && pszAnnoFile[0] != '\0')
	{
	if((m_pDistBEDFile = new CBEDfile)==NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CBEDfile");
		Reset();
		return(eBSFerrObj);
		}

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Starting to load %s",pszAnnoFile);
	if((Rslt=m_pDistBEDFile->Open(pszAnnoFile,eBTAnyBed)) !=eBSFSuccess)
		{
		while(m_pDistBEDFile->NumErrMsgs())
			gDiagnostics.DiagOut(eDLFatal,gszProcName,m_pDistBEDFile->GetErrMsg());
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open '%s' for processing",pszAnnoFile);
		Reset();
		return(eBSFerrOpnFile);
		}
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Load completed, iterating elements...");
	}

if(bFirst == true && FMode == eFMbed)
	{
	BuffIdx = sprintf(szLineBuff,"track type=bed name=\"%s\" description=\"%s\"\n",pszTitle,pszTitle);
	CUtility::SafeWrite(m_hRsltsFile,szLineBuff,BuffIdx);
	bFirst = false;
	}

BuffIdx = 0;
pROI = m_pROIs;
for(ROIidx = 0; ROIidx < m_NumOfROIs; ROIidx++,pROI++)
	{
	pChrom = &m_ChromCnts[pROI->ChromIdx-1];
	if(m_pDistBEDFile != NULL)
		ChromID = m_pDistBEDFile->LocateChromIDbyName(pChrom->szChrom);
	else
		ChromID = 0;

	switch(FMode) {
		case eFMsumCSV:
		case eFMallCSV:
			DistanceToFeatures(AnnoFileID,FeatDistStrand,ChromID,pChrom->szChrom,pROI);
			if(FMode == eFMallCSV)
				{
				BuffIdx+=sprintf(&szLineBuff[BuffIdx],"%d,\"ROI\",\"%s\",\"%s\",%d,%d,%d,\"%c\",\"%s\",\"%c\",%d,\"%s\",\"%s\",\"%c\",%d,\"%s\",%1.3f\n",
								pROI->RegionID,pszTitle,pChrom->szChrom,pROI->StartOfRegion,pROI->EndOfRegion,1 + pROI->EndOfRegion - pROI->StartOfRegion,
											pROI->Strand,
								m_pszFeatFiles[AnnoFileID-1],
							    pROI->USFeatStrand,pROI->USFeatDist,pROI->szUSFeatID,
								m_pszFeatFiles[AnnoFileID-1],
								pROI->DSFeatStrand,pROI->DSFeatDist,pROI->szDSFeatID,
								pROI->BPKM);
				}
			break;

		case eFMbed:
			BuffIdx+=sprintf(&szLineBuff[BuffIdx],"%s\t%d\t%d\tROI%d\t%d\t%c\n",
							pChrom->szChrom,pROI->StartOfRegion,pROI->EndOfRegion,pROI->RegionID,(int)pROI->BPKM,pROI->Strand == '*' ? '+' : pROI->Strand);
			break;
		}

									
	if((BuffIdx + 200)> sizeof(szLineBuff))
		{
		CUtility::SafeWrite(m_hRsltsFile,szLineBuff,BuffIdx);
		BuffIdx = 0;
		}
	}

if(FMode == eFMsumCSV && bFinal)
	{
	BuffIdx = 0;
	pROI = m_pROIs;
	for(ROIidx = 0; ROIidx < m_NumOfROIs; ROIidx++,pROI++)
		{
		pChrom = &m_ChromCnts[pROI->ChromIdx-1];
	
		BuffIdx+=sprintf(&szLineBuff[BuffIdx],"%d,\"ROI\",\"%s\",\"%s\",%d,%d,%d,\"%c\",\"%s\",\"%c\",%d,\"%s\",\"%s\",\"%c\",%d,\"%s\",%1.3f\n",
								pROI->RegionID,pszTitle,pChrom->szChrom,pROI->StartOfRegion,pROI->EndOfRegion,1 + pROI->EndOfRegion - pROI->StartOfRegion,
								pROI->Strand,
								m_pszFeatFiles[pROI->USFeatFileID-1],
								pROI->USFeatStrand,pROI->USFeatDist,pROI->szUSFeatID,
								m_pszFeatFiles[pROI->DSFeatFileID-1],
								pROI->DSFeatStrand,pROI->DSFeatDist,pROI->szDSFeatID,pROI->BPKM);
								
		if((BuffIdx + 200) > sizeof(szLineBuff))
			{
			CUtility::SafeWrite(m_hRsltsFile,szLineBuff,BuffIdx);
			BuffIdx = 0;
			}
		}
	}

if(BuffIdx)
	CUtility::SafeWrite(m_hRsltsFile,szLineBuff,BuffIdx);
delete m_pDistBEDFile;
m_pDistBEDFile = NULL;
return(m_NumOfROIs);
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

int											// returns number of elements accepted or error result if < 0
GenCSVWiggle(int Limit,						// limit (0 if no limit) processing to this many elements total
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
int NumFiltStrand;

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
NumFiltStrand = 0;
while((Rslt=m_pCSVFile->NextLine()) > 0)	// onto next line containing fields
	{
	NumFields = m_pCSVFile->GetCurFields();
	if(NumFields < 7)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expected at least 7 fields in '%s', GetCurFields() returned '%d'",pszInFile,NumFields);
		return(eBSFerrFieldCnt);
		}
	NumEls += 1;
	if(Strand != '*')
		{
		if(NumFields < 8)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Strand filtering, expected at least 8 fields in '%s', GetCurFields() returned '%d'",pszInFile,NumFields);
			return(eBSFerrFieldCnt);
			}
		m_pCSVFile->GetText(8,&pszStrand);
		if(*pszStrand != Strand)
			{
			NumFiltStrand += 1;
			continue;
			}
		}

	m_pCSVFile->GetText(4,&pszChrom);
	m_pCSVFile->GetInt(5,&StartLoci);
	m_pCSVFile->GetInt(6,&EndLoci);


	if((Rslt=BuildReadCoverage(pszChrom,StartLoci,EndLoci,1))!=eBSFSuccess)
		break;
	m_NumAcceptedReads += 1;
	if(Limit && (NumEls - NumFiltStrand) > Limit)
		break;
	}
if(Rslt == eBSFSuccess)
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"From %d elements there were %d accepted and %d filtered out because of strand",NumEls,NumEls-NumFiltStrand,NumFiltStrand);
return(Rslt == eBSFSuccess ? (NumEls-NumFiltStrand) : Rslt);
}

int GenBEDWiggle(int Limit,						// limit (0 if no limit) processing to this many bases total
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
int NumFiltStrand;

if((m_pBEDFile = new CBEDfile)==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CBEDfile");
	Reset();
	return(eBSFerrObj);
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Starting to load reads %s",pszInFile);
if((Rslt=m_pBEDFile->Open(pszInFile,eBTAnyBed)) !=eBSFSuccess)
	{
	while(m_pBEDFile->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,m_pBEDFile->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open '%s' for processing",pszInFile);
	Reset();
	return(eBSFerrOpnFile);
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Load completed, iterating reads...");

// now iterate over the reads 
szPrevChrom[0] = '\0';
CurFeatureID = 0;
IntergenicStart = 0;
RefID = 0;
Rslt = eBSFSuccess;
NumEls = 0;
NumFiltStrand = 0;
while(Rslt == eBSFSuccess && (CurFeatureID = m_pBEDFile->GetNextFeatureID(CurFeatureID)) > 0)
	{
	NumEls += 1;
	if(Limit && NumEls > Limit)
		{
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Reached limit of %d reads",NumEls);
		break;
		}

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
			{
			NumFiltStrand += 1;
			continue;
			}
		}
	
	Rslt=BuildReadCoverage(szChrom,StartLoci,EndLoci,1);
	if(Rslt >= eBSFSuccess)
		m_NumAcceptedReads += 1;
	}
if(Rslt == eBSFSuccess)
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"From %d reads there were %d accepted and %d filtered out because of strand",NumEls,NumEls-NumFiltStrand,NumFiltStrand);
return(Rslt == eBSFSuccess ? (NumEls-NumFiltStrand) : Rslt);
}


int GenSAMWiggle(int Limit,						// limit (0 if no limit) processing to this many bases total
				 char ROIStrand,				// ROI strand
				 char *pszInFile)				// SAM input file
{
int Rslt;
int StartLoci;
char szDescriptor[128];			// parsed out descriptor
int Flags;						// parsed out flags
int MAPQ;
int TLen;
char szCigar[128];
char szRNext[128];
int PNext;
char szChrom[128];
char szLine[16000];				// buffer input lines
char Strand;
int NumEls;
int NumFiltStrand;
char *pTxt;
int LineLen;

CSAMfile BAMfile;

// open SAM for reading
if(pszInFile == NULL || *pszInFile == '\0')
	return(eBSFerrParams);
if((Rslt = BAMfile.Open(pszInFile)) != eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"GenSAMWiggle: Unable to load reads from from '%s'",pszInFile);
	return((teBSFrsltCodes)Rslt);
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Starting to load reads %s",pszInFile);

NumEls = 0;
m_NumAcceptedReads = 0;
NumFiltStrand = 0;
Rslt = eBSFSuccess;
while((LineLen = BAMfile.GetNxtSAMline(szLine)) > 0)
	{
	szLine[sizeof(szLine)-1] = '\0';
	pTxt = TrimWhitespace(szLine);
	if(*pTxt=='\0' || *pTxt=='@')	// simply slough lines which were just whitespace or start with '@'
		continue;

	NumEls += 1;
	if(Limit && NumEls > Limit)
		{
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Reached limit of %d reads",NumEls);
		break;
		}

	// expecting to parse as "%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t", szDescriptor, Flags, m_szSAMTargChromName, StartLoci+1,MAPQ,szCigar,pszRNext,PNext,TLen);
	// interest is in the chromname, startloci, and length
	if((Rslt = sscanf(szLine,"%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t",szDescriptor, &Flags, szChrom, &StartLoci,&MAPQ,szCigar,szRNext,&PNext,&TLen)) != 9)
		{
		continue;
		}

	if(Flags & 0x04 || szCigar[0] == '*')	// unmapped?
		continue;

	Strand = Flags & 0x010 ? '-' : '+';
	if(ROIStrand != '*')
		{
		if(ROIStrand != Strand)
			{
			NumFiltStrand += 1;
			continue;
			}
		}

	if(szRNext[0] != '*' && szRNext[0] != '=')
		{
		if(!stricmp(szRNext,szChrom)) // if read ends aligning onto separate chroms or contigs then treat as if SE
			szRNext[0] = '=';
		else
			szRNext[0] = '*';
		}
	
	if(szRNext[0] == '*')		// SE alignment
		{
		TLen = BAMfile.CigarAlignLen(szCigar);
		}
	else                        // else was a PE alignment with both ends on the same chrom or contig   
		{
		// treating the TLen as the length but only wanting to count it once, not twice for each end...
		if(PNext < StartLoci)   // counting the forward TLen only for coverage
			{
			m_NumAcceptedReads += 1;
			continue;
			}
		}

	Rslt=BuildReadCoverage(szChrom,StartLoci-1,StartLoci + TLen - 2,1);
	if(Rslt >= eBSFSuccess)
		m_NumAcceptedReads += 1;
	}
BAMfile.Close();
if(Rslt == eBSFSuccess)
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"From %d reads there were %d accepted and %d filtered out because of strand",NumEls,NumEls-NumFiltStrand,NumFiltStrand);
return(Rslt == eBSFSuccess ? (NumEls-NumFiltStrand) : Rslt);
}

// need function which returns distance to nearest feature...
bool								// returns false if errors
DistanceToFeatures(int AnnoFileID,				// annotated file identifer
				   char FiltStrand,				// only interested in distances to nearest feature which is on this strand
				   int ROIChromID,				// ROI is on this chromosome
				   char *pszROIChrom,			// ROI chrom name
				   tsROI *pROI)
{
int Rslt;
int FeatID;
int FeatStart;
int FeatEnd;
char FeatStrand;
char szFeatChrom[cMaxDatasetSpeciesChrom+1];
char szFeatName[cMaxGeneNameLen+1];
int Nth;

// assume unable to locate a nearby feature
if(pROI->USFeatDist < 0)
	{
	pROI->USFeatStrand =  '*' ? '+' : FiltStrand;
	pROI->USFeatDist = cDfltMaxDist;
	pROI->USFeatFileID = AnnoFileID;
	strcpy(pROI->szUSFeatID,cpszNoFeat);
	}
if(pROI->DSFeatDist < 0)
	{
	pROI->DSFeatStrand = '*' ? '+' : FiltStrand;
	pROI->DSFeatDist = cDfltMaxDist;
	pROI->DSFeatFileID = AnnoFileID;
	strcpy(pROI->szDSFeatID,cpszNoFeat);
	}
if(m_pDistBEDFile == NULL || ROIChromID < 1)
	return(false);

// check if contained or even partially overlapping a feature
Nth = 0;
while((FeatID = m_pDistBEDFile->LocateFeatureIDinRangeOnChrom(ROIChromID,pROI->StartOfRegion,pROI->EndOfRegion,++Nth)) > 0)
	{
	Rslt = m_pDistBEDFile->GetFeature(FeatID,szFeatName,szFeatChrom,&FeatStart,&FeatEnd,NULL,&FeatStrand);
	if(FiltStrand != '*' && FeatStrand != FiltStrand)
		continue;
	pROI->DSFeatFileID = AnnoFileID;
	pROI->USFeatFileID = AnnoFileID;
	pROI->USFeatDist = 0;
	pROI->DSFeatDist = 0;
	pROI->USFeatStrand = FeatStrand;
	pROI->DSFeatStrand = FeatStrand;
	strcpy(pROI->szDSFeatID,szFeatName);
	strcpy(pROI->szUSFeatID,szFeatName);
	return(true);
	}

// now check for feature which is on strand and downstream of the ROI
Rslt = -1;
if((FeatID = m_pDistBEDFile->LocateFeatureAfter(ROIChromID,pROI->EndOfRegion)) > 0)
	{
	// determine feature strand and start loci
	while((Rslt = m_pDistBEDFile->GetFeature(FeatID,szFeatName,szFeatChrom,&FeatStart,NULL,NULL,&FeatStrand))==0)
		{
		// ensure still on same chromosome as ROI
		if(stricmp(szFeatChrom,pszROIChrom))
			{
			Rslt = 0;
			break;
			}
		// feature start must be after ROIend and strand also has to match
		if(FeatStart > pROI->EndOfRegion && (FiltStrand == '*' || FeatStrand == FiltStrand))
			{
			Rslt = 1;
			break;
			}
		FeatID += 1;	// try next feature, feature ID's are ordered by feature chrom.start loci
		}
	if(Rslt == 1)
		{
		if(pROI->DSFeatDist >= (FeatStart - pROI->EndOfRegion))
			{
			pROI->DSFeatFileID = AnnoFileID;
			pROI->DSFeatDist = FeatStart - pROI->EndOfRegion;
			pROI->DSFeatStrand = FeatStrand;
			strcpy(pROI->szDSFeatID,szFeatName);
			}
		}
	}
if(FeatID < 0 || Rslt < 0)
	return(false);

// now check for feature which is upstream of ROI
if((FeatID = m_pDistBEDFile->LocateFeatureBefore(ROIChromID,pROI->StartOfRegion)) > 0)
	{
	// determine feature strand and end loci
	while(FeatID > 0 && (Rslt = m_pDistBEDFile->GetFeature(FeatID,szFeatName,szFeatChrom,NULL,&FeatEnd,NULL,&FeatStrand))==0)
		{
		// ensure still on same chromosome as ROI
		if(stricmp(szFeatChrom,pszROIChrom))
			{
			Rslt = 0;
			break;
			}
		// feature has to end before the ROI starts
		if(FeatEnd < pROI->StartOfRegion && (FiltStrand == '*' || FeatStrand == FiltStrand))
			{
			Rslt = 1;
			break;
			}
		FeatID -= 1;	// try previous feature
		}
	if(FeatID > 0 && Rslt == 1)
		{
		if(pROI->USFeatDist >= (pROI->StartOfRegion - FeatEnd))
			{
			pROI->USFeatFileID = AnnoFileID;
			pROI->USFeatDist = pROI->StartOfRegion - FeatEnd;
			pROI->USFeatStrand = FeatStrand;
			strcpy(pROI->szUSFeatID,szFeatName);
			}
		}
	}
return((FeatID < 0 || Rslt < 0) ? false : true);
}


int 
FilterRegions(etBEDRegion Region,			// which regions are of interest
				 char ROIStrand,			// region on this strand
				 char *pszInFile)			// UCSC BED containing regions
{
int Rslt;
int CurFeatureID;
int NumExons;
int NumIntrons;
int CDSstart;
int CDSend;
int StartLoci;
int EndLoci;
int Score;
char szChrom[128];
char szPrevChrom[128];
char szFeatName[128];
char Strand;
int RefID;
int IntergenicStart;
int Idx;
int NumEls;
int NumFiltStrand;

if(pszInFile[0] == '\0')
	return(eBSFSuccess);

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

if(!m_pBEDFile->ContainsGeneDetail() && Region > eMEGRIntergenic)			// returns true if file contains gene detail (utr/cds/intron etc)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"retention bed file %s does not contain gene regions",pszInFile);
	Reset();
	return(eBSFerrFeature);
	}

// now iterate over the features, filtering as may be appropriate
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Load completed, processing ROI against retention features...");
szPrevChrom[0] = '\0';
CurFeatureID = 0;
IntergenicStart = 0;
RefID = 0;
Rslt = eBSFSuccess;
NumEls = 0;
NumFiltStrand = 0;
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
			{
			NumFiltStrand += 1;
			continue;
			}
		}
	
	if(Region != eMEGRAny)
		{
		NumExons = m_pBEDFile->GetNumExons(CurFeatureID);					// returns number of exons - includes UTRs + CDS
		NumIntrons = m_pBEDFile->GetNumIntrons(CurFeatureID);
		CDSstart = StartLoci + m_pBEDFile->GetCDSStart(CurFeatureID);		// returns relative start offset of CDS - NOTE add to '+' strand gene start, subtract on '-' strand gene start
		CDSend = StartLoci + m_pBEDFile->GetCDSEnd(CurFeatureID);			// returns relative end offset of CDS - NOTE add to '+' strand gene start, subtract on '-' strand gene start
		}
	Rslt = eBSFSuccess;

	switch(Region) {
		case eMEGRAny:			// retain any region
			Rslt=RetainReadCoverage(szChrom,StartLoci,EndLoci);
			continue;

		case eMEGRIntergenic:	// only retain intergenic
			if(IntergenicStart < StartLoci)
				Rslt=RetainReadCoverage(szChrom,IntergenicStart,StartLoci-1);
			if(IntergenicStart <= EndLoci)
				IntergenicStart = EndLoci+1;
			continue;

		case eMEGRExons:		// only retain exons
			for(Idx = 1; Idx <= NumExons; Idx++)
				{
				StartLoci = m_pBEDFile->GetExonStart(CurFeatureID,Idx);
				EndLoci   = m_pBEDFile->GetExonEnd(CurFeatureID,Idx);
				if(StartLoci <= EndLoci)
					Rslt=RetainReadCoverage(szChrom,StartLoci,EndLoci);
				}
			continue;

		case eMEGRIntrons:		// only retain introns
			for(Idx = 1; Idx <= NumIntrons; Idx++)
				{
				StartLoci = m_pBEDFile->GetIntronStart(CurFeatureID,Idx);
				EndLoci = m_pBEDFile->GetIntronEnd(CurFeatureID,Idx);
				if(StartLoci <= EndLoci)
					Rslt=RetainReadCoverage(szChrom,StartLoci,EndLoci);
				}
			continue;

		case eMEGRCDS:			// only retain CDSs
			for(Idx = 1; Idx <= NumExons; Idx++)
				{
				StartLoci = m_pBEDFile->GetExonStart(CurFeatureID,Idx);
				EndLoci   = m_pBEDFile->GetExonEnd(CurFeatureID,Idx);
				if(EndLoci < CDSstart || StartLoci > CDSend)
					continue;
				if(StartLoci < CDSstart)
					StartLoci = CDSstart;
				if(EndLoci > CDSend)
					EndLoci = CDSend;
				if(StartLoci <= EndLoci)
					RetainReadCoverage(szChrom,StartLoci,EndLoci);
				}
			continue;

		case eMEGUTR:			// only process UTRs - single exon may have both 5' and 3' UTRs
			for(Idx = 1; Idx <= NumExons; Idx++)
				{
				StartLoci = m_pBEDFile->GetExonStart(CurFeatureID,Idx);
				EndLoci   = m_pBEDFile->GetExonEnd(CurFeatureID,Idx);
				if(EndLoci <= CDSend && StartLoci >= CDSstart) // is exon CDS only?
					continue;

				// check if 5' UTR 
				if(StartLoci < CDSstart)
					{
					if(EndLoci >= CDSstart)
						Rslt=RetainReadCoverage(szChrom,StartLoci,CDSstart-1);
					else
						Rslt=RetainReadCoverage(szChrom,StartLoci,EndLoci);
					}
					
				// check if 3'UTR
				if(EndLoci > CDSend)
					{
					if(StartLoci <= CDSend)
						StartLoci = CDSend+1;
					Rslt=RetainReadCoverage(szChrom,StartLoci,EndLoci);
					}
				}
			continue;

		case eMEG5UTR:			// only process 5'UTRs - strand sensitive
			for(Idx = 1; Idx <= NumExons; Idx++)
				{
				StartLoci = m_pBEDFile->GetExonStart(CurFeatureID,Idx);
				EndLoci   = m_pBEDFile->GetExonEnd(CurFeatureID,Idx);
				if(EndLoci <= CDSend && StartLoci >= CDSstart) // is exon CDS only?
					continue;

				if(Strand != '-')
					{
					// check if 5' UTR on '+' strand 
					if(StartLoci < CDSstart)
						{
						if(EndLoci >= CDSstart)
							EndLoci = CDSstart - 1;
						}
					}
				else 
					{	
					// check if 5'UTR on '-' strand
					if(EndLoci > CDSend)
						{
						if(StartLoci <= CDSend)
							StartLoci = CDSend+1;
						}
					}
				Rslt=RetainReadCoverage(szChrom,StartLoci,EndLoci);
				}
			continue;

		case eMEG3UTR:			// only process 3'UTRs  - strand sensitive
			for(Idx = 1; Idx <= NumExons; Idx++)
				{
				StartLoci = m_pBEDFile->GetExonStart(CurFeatureID,Idx);
				EndLoci   = m_pBEDFile->GetExonEnd(CurFeatureID,Idx);
				if(EndLoci <= CDSend && StartLoci >= CDSstart) // is exon CDS only?
					continue;

				if(Strand == '-')
					{
					// check if 3' UTR on '-' strand 
					if(StartLoci < CDSstart)
						{
						if(EndLoci >= CDSstart)
							EndLoci = CDSstart - 1;
						}
					}
				else 
					{	
					// check if 3'UTR on '+' strand
					if(EndLoci > CDSend)
						{
						if(StartLoci <= CDSend)
							StartLoci = CDSend+1;
						}
					}
				Rslt=RetainReadCoverage(szChrom,StartLoci,EndLoci);
				}
			continue;
		}
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing of %d retention elements completed of which %d filtered out because of strand",NumEls,NumFiltStrand);
return(Rslt);
}


int
Process(etPMode PMode,					// processing mode
		etFMode FMode,					// output format - CSV or BED
		char *pszTitle,					// CSV species or title used if BED output format
		char ReadStrand,				// only accept reads from this strand ('*' if either)
		char RestainStrand,				// filter ROI by filter elements on this strand ('*' if either)
		char FeatDistStrand,			// distances to features on this strand ('*' if either)
		etBEDRegion Region,				// filter by retaining this functional region
		int Limit,						// limit (0 if no limit) processing to this many bases total
		int MinMedianCov,				// minimin median coverage required in reported regions
		int MinRegionLen,				// report regions which are of at least this length
		int MaxGapLen,					// report regions containing no read gaps of <= this length 
		char *pszRetainFile,			// file containing annotated regions to retain
		int NumAssocFiles,				// number of association files in pszAssocFiles[]
		char *pszAssocFiles[],			// distance association files
		char *pszInFile,				// input CSV or BED file containing read alignment loci
		char *pszRsltsFile)				// output CSV region loci file
{
int Rslt;
bool bNoFilter;
char *pszAssocFile;
char *pszInLociFile;
int Idx;

Init();

CSimpleGlob glob(SG_GLOB_FULLSORT);
if (glob.Add(pszInFile) >= SG_SUCCESS)
	{
	for (int n = 0; n < glob.FileCount(); ++n)
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Will process in this order: %d '%s", n+1,glob.File(n));
	}
else
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to glob '%s",pszInFile);
	Reset();
	return(eBSFerrOpnFile);	
    }



#ifdef _WIN32
if((m_hRsltsFile = open(pszRsltsFile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
if((m_hRsltsFile = open(pszRsltsFile, O_RDWR | O_CREAT | O_TRUNC, S_IREAD | S_IWRITE ))==-1)
#endif
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create results file: %s - %s",pszRsltsFile,strerror(errno));
	Reset();
	return(eBSFerrCreateFile);
	}

m_ReadStrand = ReadStrand;
m_RestainStrand = RestainStrand;
m_FeatDistStrand = FeatDistStrand;

for (int n = 0; n < glob.FileCount(); ++n)
	{
	pszInLociFile = glob.File(n);

	etClassifyFileType FileType;
	if(PMode == ePMdefault)
		FileType = CUtility::ClassifyFileType(pszInLociFile);
	else
		FileType = (etClassifyFileType)(PMode - 1);

	switch(FileType) {
		case eCFTopenerr:		// unable to open file for reading
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open file: '%s'",pszInLociFile);
			return(eBSFerrOpnFile);

		case eCFTlenerr:		// file length is insufficent to classify type
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"Unable to classify file type (insuffient data points): '%s'",pszInLociFile);
			return(eBSFerrFileAccess);

		case eCFTunknown:		// unable to reliably classify
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"Unable to reliably classify file type: '%s'",pszInLociFile);
			return(eBSFerrFileType);

		case eCFTCSV:		// CSV loci processing
			Rslt = GenCSVWiggle(Limit,ReadStrand,pszInLociFile);
			break;

		case eCFTBED:			// UCSC BED processing
			Rslt = GenBEDWiggle(Limit,ReadStrand,pszInLociFile);
			break;

		case eCFTSAM:			// file has been classified as being SAM
			Rslt = GenSAMWiggle(Limit,ReadStrand,pszInLociFile);
			break;

		default:
			break;
		}
	if(Rslt < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Errors processing file: %s",pszInLociFile);
		Reset();
		return(Rslt);
		}

	if(Rslt == 0)
		gDiagnostics.DiagOut(eDLWarn,gszProcName,"No accepted reads in file: %s",pszInLociFile);
	}

if(!m_NumAcceptedReads)
	{
	gDiagnostics.DiagOut(eDLWarn,gszProcName,"Nothing to do, no accepted reads in processed source files");
	Reset();
	return(1);
	}


// now filter by genomic regions
if(pszRetainFile != NULL && pszRetainFile[0] != '\0')
	bNoFilter = false;
else
	bNoFilter = true;

if(Rslt > eBSFSuccess && !bNoFilter)
	Rslt = FilterRegions(Region,RestainStrand,pszRetainFile);

if(Rslt < eBSFSuccess)
	{
	Reset();
	return(Rslt);
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Identification of ROIs started...");
if((Rslt = IdentROI(FMode,			// output in this format to m_hRsltsFile
			MinMedianCov,			// minimum median coverage required in reported regions
			MinRegionLen,			// report regions which are of at least this length
			MaxGapLen,				// report regions containing no read gaps of <= this length 
			bNoFilter)) < 1)		// true if all counts to be processed
	{
	Reset();
	return(Rslt);
	}

char szAssocName[_MAX_PATH];
char szAssocExtn[_MAX_PATH];

if(Rslt >= 1)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Identified %d ROIs",Rslt);	
	if(NumAssocFiles)
		{
		for(Idx = 0; Idx < NumAssocFiles; Idx++)
			{
			pszAssocFile = pszAssocFiles[Idx];
#ifdef _WIN32
			_splitpath(pszAssocFile,NULL,NULL,szAssocName,szAssocExtn);
			strcat(szAssocName,szAssocExtn);
#else
			CUtility::splitpath(pszAssocFile,NULL,szAssocName);
#endif
			m_pszFeatFiles[Idx] = new char [strlen(szAssocName)+1];
			strcpy(m_pszFeatFiles[Idx],szAssocName);
			m_NumFeatFiles += 1;

			gDiagnostics.DiagOut(eDLInfo,gszProcName,"distances to feature processing on %s started...",pszAssocFile);
			Rslt = WriteRegions(((Idx + 1) == NumAssocFiles) ? true : false, pszAssocFile,Idx+1,FMode,pszTitle,FeatDistStrand);
			if(Rslt < 0)
				{
				gDiagnostics.DiagOut(eDLInfo,gszProcName,"Errors whilst identifying ROIs");
				Reset();
				return(Rslt);
				}
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"distances to feature processing on %s completed",pszAssocFile);
			}
		}
	else
		{
		Rslt = WriteRegions(true, NULL,1,FMode,pszTitle,FeatDistStrand);
		if(Rslt < 0)
			{
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"Errors whilst identifying ROIs");
			Reset();
			return(Rslt);
			}
		}
	}
else
	if(Rslt == 0)
		gDiagnostics.DiagOut(eDLWarn,gszProcName,"No ROI were identified");
Reset();
return(Rslt);
}




