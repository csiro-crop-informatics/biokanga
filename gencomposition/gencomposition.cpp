// gencomposition.cpp : Defines the entry point for the console application.


#include "stdafx.h"
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#if _WIN32

#include "../libbiokanga/commhdrs.h"
#else
#include "../libbiokanga/commhdrs.h"
#endif

const unsigned int cProgVer = 106;		// increment with each release

const int cMaxNMerLen = 12;				// maximum length N-Mer (don't change unless code dependencies checked and plenty of memory!)
const int cDfltNMerLen= 5;				// default N-Mer length

const int cMaxIncludeLoci = 10;			// maximum number of include loci filter files
const int cMaxExcludeLoci = 10;			// maximum number of exclude loci filter files

const int cMaxIncludeChroms = 20;		// max number of include chromosomes regular expressions
const int cMaxExcludeChroms = 20;		// max number of exclude chromosomes regular expressions

const int cMaxLengthRange = 750000000;	// maximal element length when sequence is from bioseq file
const int cDfltMinLengthRange = 10;		// default minimum element length

const int cDfltAllocSeqLen = 100000;	// default allocation for holding subsequence from bioseq file

typedef enum TAG_eProcMode {
	eProcModeNMerDistAllSeqs,			// process N-Mer frequency distribution over all sequences
	eProcModeNMerDistPerSeq				// process per sequence N-Mer frequency distribution
} etProcMode;

typedef enum TAG_eRsltsFormat {
	eRsltsMultifasta = 0,				// output sequences as multifasta
	eRsltsFasta							// output sequences concatenated as single fasta record
	} teRsltsFormat;


typedef struct TAG_sFiltState {
	bool bFiltOut;						// true if this element has been marked to be filtered out
	unsigned int fOverLen:1;			// element is overlength
	unsigned int fUnderLen:1;			// element is underlength
	unsigned int fExcludeLoci:1;		// element is on loci to be excluded
	unsigned int fIncludeLoci:1;		// element is not on loci to be included
	unsigned int fChrom:1;				// element is on filtered chrom
	unsigned int fOutRefID:1;			// element filtered RefID
	unsigned int fInRefID:1;			// element filtered RefID
	unsigned int fExcludeRegion:1;		// element filtered because in regions to be filtered out
	unsigned int fIncludeRegion:1;		// element filtered because not in regions to be retained
	} tsFiltState;


typedef struct TAG_sElChrom {
	int ChromID;							// unique identifier for this chromosome
	UINT16 Hash;							// hash over chromosome name
	char szChrom[cMaxDatasetSpeciesChrom];	// chromosome
} tsElChrom;

const int cMaxExcludeHistory = 1000;
typedef struct TAG_sExcludeEl {
	struct TAG_sExcludeEl *pNext;
	struct TAG_sExcludeEl *pPrev;
	char szChrom[cMaxDatasetSpeciesChrom];	// which chromosome
	bool bExclude;		// true if to be excluded, false if not
	} tsExcludeEl;

typedef struct TAG_sProcParams 
	{
	etProcMode ProcMode;		// processing mode

	int MaxNMerLen;				// max length N-Mer to process
	int MinNMerLen;				// min N-Mer to process

	int RegionsIn;				// exclusive regions to be retained
	int RegionsOut;				// regions to be filtered out

	teCSVFormat CSVFormat;		// expected input CSV format
	char *pszInLociFile;		// input CSV loci file
	CCSVFile *pCSV;

	char *pszInSeqFile;			// file containing genome assembly sequences
	CBioSeqFile *pBioseq;		// to hold opened bioseq assembly file containing chromosome sequences
	
	char *pszCurChrom;			// sequence in pSeq is from this chromosome
	int StartLoci;				// sequence starts at this chrom loci
	int EndLoci;				// sequence ends at this chrom loci
	int CurSeqLen;				// current length of sequence loaded into pSeq
	int CurRegion;				// characterised as being in this region
	int AllocdSeqLen;			// currently allocated for subsequences
	etSeqBase *pSeq;			// memory to hold loci subsequence loaded from pBioSeq

	INT64 *pCntStepCnts;		// array of stats counters organised in [Step] order
	int NumCntSteps;				// number of steps in pCntStepCnts
	int NumRegions;					// number of regions per step
	int CntStepOfs[cMaxNMerLen];	// ofs into pCntStepCnts for each NMer length

	teRsltsFormat RsltsFormat;		// output results format
	char *pszRsltsFile;				// output results file to create/write
	int hRsltsFile;					// file handle for results file

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
	int	MinElLen;					// minimum element length required
	int	MaxElLen;					// maximum element length required

	CFilterRefIDs *pFilterInRefIDs;	// used when filtering in by RefID
	CFilterRefIDs *pFilterOutRefIDs;// used when filtering out by RefID

	int gNumExcludeEls;				// current number of elements in gExcludeChroms
	tsExcludeEl *gpMRA;				// pts to most recently accessed or added
	tsExcludeEl *gpLRA;				// pts to least recently accessed
	tsExcludeEl gExcludeChroms[cMaxExcludeHistory];

	CFilterLoci *pExcludeFiltLoci;	// exclude these loci
	CFilterLoci *pIncludeFiltLoci;	// include these loci

	int CachElChromID;				// identifier of last retrieved chromosome from pElChroms
	tsElChrom *pElChroms;			// pts to mem allocd to hold array of chromosome names
	int MaxElChroms;				// max allocd chroms
	int NumElChroms;				// number of chroms

} tsProcParams; 


int Process(etProcMode ProcMode,		// processing mode
			teCSVFormat CSVFormat,		// expected input CSV format
			char *pszInLociFile,		// CSV file containing elements
			char *pszInSeqFile,				// file containing genome assembly sequences
			teRsltsFormat RsltsFormat,	// results format 
			char *pszRsltsFile,				// write results to this file
			int FilterRegionsIn,		// retain any of these (exclusive) regions
			int FilterRegionsOut,		// remove any of these regions
			char *pszFilterOutRefIDFile, // RefIDs to filter out
			char *pszFilterInRefIDFile,	 // RefIDs to filter in unless filtered out
			int MinElLen,				 // elements must be of at least this length
			int MaxElLen,				 // no longer than this length
			int MaxNMerLen,				// max length N-Mer to process
			int MinNMerLen,				// min N-Mer to process
			int NumIncludeChroms,		 // number of chromosome regular expressions to include
			char **ppszIncludeChroms,	 // array of include chromosome regular expressions
			int NumExcludeChroms,		 // number of chromosome expressions to exclude
			char **ppszExcludeChroms,	 // array of exclude chromosome regular expressions
			int NumIncludeLoci,			 // number of include region files 
			char **ppszIncludeLoci,		 // array of include region files
			int NumExcludeLoci,			 // number of exclude region files
			char **ppszExcludeLoci);	 // array of exclude region files



char *RsltsFormat2Text(teRsltsFormat Format);
char *CSVFormat2Text(teCSVFormat Format);


char *ProcMode2Txt(etProcMode ProcMode);
int TrimQuotes(char *pszTxt);
int ParseRegions(char *pszRegionList);
char *Regions2Txt(int Regions);
bool ExcludeThisChrom(char *pszChrom,tsProcParams *pProcParams);
tsExcludeEl *LocateExclude(char *pszChrom,tsProcParams *pParams);
bool AddExcludeHistory(char *pszChrom,bool bExclude,tsProcParams *pParams);
int FilterCSV(tsProcParams *pProcParams);
int ProcessSequence(tsProcParams *pParams);
UINT16 GenNameHash(char *pszName);
tsElChrom *LocateElChrom(char *pszChrom,tsProcParams *pProcParams);
int AddChrom(char *pszChrom,tsProcParams *pParams);
int GenSeqIdx(int SeqLen,etSeqBase *pSeq);
char *StepIdx2Seq(int SeqLen,int SeqIdx);
int OutputResults(tsProcParams *pParams);

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
	return _T("gencomposition");
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
int iProcMode;				// processing mode
int iCSVFormat;				// expected input loci file format
int iRsltsFormat;			// output results in this format
bool bSkipFirst;			// true if first line contains header and should be skipped
int iMinLength;				// core elements must be of at least this length
int iMaxLength;				// and no longer than this length
int iMinNMerLen;			// composition for this N-Mer upto 
int iMaxNMerLen;			// this N-Mer 

int iRegionsIn;				// number of regions to retain
int iRegionsOut;			// number of regions to exclude
char szRegionsIn[128];		// regions to retain
char szRegionsOut[128];		// regions to exclude

int iNumIncludeChroms;		// number of chromosome filtering patterns to retain
char *pszIncludeChroms[cMaxIncludeChroms];
int iNumExcludeChroms;		// number of chromosome filtering patterns to exclude
char *pszExcludeChroms[cMaxExcludeChroms];
int iNumIncludeLoci;			// number of loci files with loci to retain
char *pszIncludeLoci[cMaxIncludeLoci];
int iNumExcludeLoci;			// number of loci files with loci to exclude
char *pszExcludeLoci[cMaxExcludeLoci];

char szInLociFile[_MAX_PATH];	// input element loci from this file
char szInSeqFile[_MAX_PATH];	// input bioseq file containing assembly
char szRsltsFile[_MAX_PATH];	// output loci + sequences or gene identifiers to this file

int LenFileList;	
int Idx;
char szFilterOutRefIDFile[_MAX_PATH];
char szFilterInRefIDFile[_MAX_PATH];


// command line args
struct arg_lit  *help    = arg_lit0("hH","help",                "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_int *ScreenLogLevel=arg_int0("S", "ScreenLogLevel",	"<int>","Level of diagnostics written to logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"diagnostics log file");

struct arg_int  *ProcMode = arg_int0("m","mode","<int>",		"processing mode 0: over all sequences, 1: per sequence");
struct arg_int  *CSVFormat = arg_int0("M","informat","<int>",	"input CSV file type 0:Loci 1:locsfx probe 2:locsfx target");
struct arg_lit  *SkipFirst    = arg_lit0("x","skipfirst",       "skip first line of input CSV - header line");
struct arg_file *InLociFile = arg_file0("i","inloci","<file>",	"element loci CSV file");
struct arg_file *InSeqFile = arg_file1("I","assembly","<file>",	"genome assembly bioseq file");
struct arg_file *RsltsFile = arg_file1("o","output","<file>",	"output file");
struct arg_int  *RsltsFormat = arg_int0("p","outformat","<int>","output as 0: default");
struct arg_int  *MinLength = arg_int0("l","minlength","<int>",	"minimum element length (default 10)");
struct arg_int  *MaxLength = arg_int0("L","maxlength","<int>",	"maximum element length (default 1000000000)");
struct arg_int *MinNMerLen=arg_int0("k", "minnmerlen","<int>",  "process N-Mers from this length up (default is 1)");
struct arg_int *MaxNMerLen=arg_int0("K", "maxnmerlen","<int>",  "process N-Mers upto and including this length (default is 5)");

struct arg_str  *RegionsOut = arg_str0("R","regionsout","<string>","Remove these regions (space or comma delimit), 1: Intergenic, 2: US, 3: 5'UTR, 4: CDS, 5: Intron, 6: 3'UTR, 7: DS, 8: 5'Splice, 9: 3'Splice");
struct arg_str  *RegionsIn = arg_str0("r","regionsin","<string>","Retain these exclusive regions (space or comma delimit), 1: Intergenic, 2: US, 3: 5'UTR, 4: CDS, 5: Intron, 6: 3'UTR, 7: DS");
struct arg_file *FilterOutRefIDFile = arg_file0("X",NULL,"<file>",	"First: filter out with RefIDs from this filter file");
struct arg_file *FilterInRefIDFile = arg_file0("x",NULL,"<file>",	"Second: filter in with RefIDs from this filter file");
struct arg_file *ExcludeLoci = arg_filen("E","exclude","<file>",0,cMaxExcludeLoci,	"First: filter out these Loci");
struct arg_file *IncludeLoci = arg_filen("e","include","<file>",0,cMaxIncludeLoci,	"Second: filter in these Loci");
struct arg_str  *ExcludeChroms = arg_strn("Z","chromexclude","<string>",0,cMaxExcludeChroms,"First: regular expressions defining species.chromosomes to exclude");
struct arg_str  *IncludeChroms = arg_strn("z","chromeinclude","<string>",0,cMaxIncludeChroms,"Second: regular expressions defining species.chromosomes to include");

struct arg_end *end = arg_end(20);

void *argtable[] = {help,version,FileLogLevel,ScreenLogLevel,LogFile,
					ProcMode,CSVFormat,InLociFile,InSeqFile,RsltsFile,RsltsFormat,SkipFirst,MinLength,MaxLength,MinNMerLen,MaxNMerLen,
					RegionsOut,RegionsIn,FilterOutRefIDFile,FilterInRefIDFile,ExcludeLoci,IncludeLoci,ExcludeChroms,IncludeChroms,
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

	iProcMode = ProcMode->count ? ProcMode->ival[0] : eProcModeNMerDistAllSeqs;
	if(iProcMode < eProcModeNMerDistAllSeqs || iProcMode > eProcModeNMerDistPerSeq)
		{
		printf("\nError: Requested processing mode '-x%d' not supported",iProcMode);
		exit(1);
		}

		iMinNMerLen = MinNMerLen->count ? MinNMerLen->ival[0] : 1;
	if(iMinNMerLen < 1 || iMinNMerLen > cMaxNMerLen)
		{
		printf("Error: min N-Mer length requested '-l%d' outside of range %d..%d\n",iMinNMerLen,1,cMaxNMerLen);
		exit(1);
		}
	if(iMinNMerLen > cDfltNMerLen)
		iMaxNMerLen = iMinNMerLen;
	else
		iMaxNMerLen = cDfltNMerLen;

	iMaxNMerLen = MaxNMerLen->count ? MaxNMerLen->ival[0] : iMaxNMerLen;
	if(iMaxNMerLen < iMinNMerLen || iMaxNMerLen > cMaxNMerLen)
		{
		printf("Error: max N-Mer length requested '-n%d' outside of range %d..%d\n",iMaxNMerLen,iMinNMerLen,cMaxNMerLen);
		exit(1);
		}
	iCSVFormat = CSVFormat->count ? CSVFormat->ival[0] : eCSVFdefault;
	if(iCSVFormat < eCSVFdefault || iCSVFormat > eCSVFtarget)
		{
		printf("\nError: expected input CSV format specified '-m%d' must be in range 0..2",iCSVFormat);
		exit(1);
		}

	if(RegionsIn->count)
		{
		strcpy(szRegionsIn,RegionsIn->sval[0]);
		TrimQuotes(szRegionsIn);
		if((iRegionsIn = ParseRegions(szRegionsIn)) < 0)
			{
			printf("Error: unable to parse '-r%s' into regions to retain",szRegionsIn);
			exit(1);
			}
		if(iRegionsIn & (cIntronExonSpliceSite | cExonIntronSpliceSite))
			{
			printf("Error: Can't specify 5'Splice and 3'Splice to be exclusively retained with '-r%s'",szRegionsIn);
			exit(1);
			}
		}
	else
		{
		szRegionsIn[0] = '\0';
		iRegionsIn = 0;
		}

	if(RegionsOut->count)
		{
		strcpy(szRegionsOut,RegionsOut->sval[0]);
		TrimQuotes(szRegionsOut);
		if((iRegionsOut = ParseRegions(szRegionsOut)) < 0)
			{
			printf("Error: unable to parse '-R%s' into regions to remove",szRegionsOut);
			exit(1);
			}
		}
	else
		{
		szRegionsOut[0] = '\0';
		iRegionsOut = 0;
		}


	if(FilterOutRefIDFile->count)
		{
		strncpy(szFilterOutRefIDFile,FilterOutRefIDFile->filename[0],sizeof(szFilterOutRefIDFile));
		szFilterOutRefIDFile[sizeof(szFilterOutRefIDFile)-1] = '\0';
		}
	else
		szFilterOutRefIDFile[0] = '\0';


	if(FilterInRefIDFile->count)
		{
		strncpy(szFilterInRefIDFile,FilterInRefIDFile->filename[0],sizeof(szFilterInRefIDFile));
		szFilterInRefIDFile[sizeof(szFilterInRefIDFile)-1] = '\0';
		}
	else
		szFilterInRefIDFile[0] = '\0';


	iNumIncludeLoci = IncludeLoci->count;
	for(Idx=0;Idx < IncludeLoci->count; Idx++)
		{
		LenFileList = (int)strlen(IncludeLoci->filename[Idx]);
		pszIncludeLoci[Idx] = new char [LenFileList+1];
		strcpy(pszIncludeLoci[Idx],IncludeLoci->filename[Idx]);
		TrimQuotes(pszIncludeLoci[Idx]);
		}
	iNumExcludeLoci = ExcludeLoci->count;
	for(Idx=0;Idx < ExcludeLoci->count; Idx++)
		{
		LenFileList = (int)strlen(ExcludeLoci->filename[Idx]);
		pszExcludeLoci[Idx] = new char [LenFileList+1];
		strcpy(pszExcludeLoci[Idx],ExcludeLoci->filename[Idx]);
		TrimQuotes(pszExcludeLoci[Idx]);
		}

	iNumIncludeChroms = IncludeChroms->count;
	for(Idx=0;Idx < IncludeChroms->count; Idx++)
		{
		LenFileList = (int)strlen(IncludeChroms->sval[Idx]);
		pszIncludeChroms[Idx] = new char [LenFileList+1];
		strcpy(pszIncludeChroms[Idx],IncludeChroms->sval[Idx]);
		TrimQuotes(pszIncludeChroms[Idx]);
		}

	iNumExcludeChroms = ExcludeChroms->count;
	for(Idx=0;Idx < ExcludeChroms->count; Idx++)
		{
		LenFileList = (int)strlen(ExcludeChroms->sval[Idx]);
		pszExcludeChroms[Idx] = new char [LenFileList+1];
		strcpy(pszExcludeChroms[Idx],ExcludeChroms->sval[Idx]);
		TrimQuotes(pszExcludeChroms[Idx]);
		}




	iRsltsFormat = RsltsFormat->count ? RsltsFormat->ival[0] : eRsltsMultifasta;
	if(iRsltsFormat < eRsltsMultifasta || iRsltsFormat > eRsltsFasta)
		{
		printf("\nError: RsltsFormat specified '-p%d' must be in range 0..1",iRsltsFormat);
		exit(1);
		}

	bSkipFirst = SkipFirst->count ? true : false;

	iMinLength = MinLength->count ? MinLength->ival[0] : cDfltMinLengthRange;
	if(iMinLength < 1 || iMinLength > cMaxLengthRange)
		{
		printf("Error: Mininum element length '-l%d' is not in range 1..%d",iMinLength,cMaxLengthRange);
		exit(1);
		}

	iMaxLength = MaxLength->count ? MaxLength->ival[0] : cMaxLengthRange;
	if(iMaxLength < iMinLength || iMaxLength > cMaxLengthRange)
		{
		printf("Error: Maximum element length '-L%d' is not in range %d..%d",iMaxLength,iMinLength,cMaxLengthRange);
		exit(1);
		}

	if(InLociFile->count)
		strncpy(szInLociFile,InLociFile->filename[0],_MAX_PATH);
	else
		szInLociFile[0] = '\0';
	szInLociFile[_MAX_PATH-1] = '\0';
	strncpy(szInSeqFile,InSeqFile->filename[0],_MAX_PATH);
	szInSeqFile[_MAX_PATH-1] = '\0';
	strncpy(szRsltsFile,RsltsFile->filename[0],_MAX_PATH);
	szRsltsFile[_MAX_PATH-1] = '\0';

		// now that command parameters have been parsed then initialise diagnostics log system
	if(!gDiagnostics.Open(szLogFile,(etDiagLevel)iScreenLogLevel,(etDiagLevel)iFileLogLevel,true))
		{
		printf("\nError: Unable to start diagnostics subsystem.");
		if(szLogFile[0] != '\0')
			printf(" Most likely cause is that logfile '%s' can't be opened/created",szLogFile);
		exit(1);
		}
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Version: %d.%2.2d Processing parameters:",cProgVer/100,cProgVer%100);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Processing mode: %s",ProcMode2Txt((etProcMode)iProcMode));
	if(szInLociFile[0] == '\0')
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Input CSV element loci file: None specified");
	else
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Input CSV element loci file: '%s'",szInLociFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Expecting CSV format as: %s",CSVFormat2Text((teCSVFormat)iCSVFormat));
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Input bioseq genome assembly file: '%s'",szInSeqFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Output to file: '%s'",szRsltsFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Generating Results as: %s",RsltsFormat2Text((teRsltsFormat)iRsltsFormat));
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"First line contains header: %s",bSkipFirst ? "yes" : "no");
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Minimum element length: %d",iMinLength);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Maximum element length: %d",iMaxLength);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Minimum N-Mer length: %d",iMinNMerLen);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Maximum N-Mer length: %d",iMaxNMerLen);

	if(iRegionsIn)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Filter in these regions: '%s'",Regions2Txt(iRegionsIn));
	if(iRegionsOut)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Filter out these regions: '%s'",Regions2Txt(iRegionsOut));
	
	if(szFilterOutRefIDFile[0])
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"filter out with RefIDs from this filter file: %s",szFilterOutRefIDFile);
	if(szFilterOutRefIDFile[0])
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"filter in with RefIDs from this filter file: %s",szFilterInRefIDFile);

	for(Idx = 0; Idx < iNumIncludeChroms; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"reg expressions defining chroms to include: '%s'",pszIncludeChroms[Idx]);
	for(Idx = 0; Idx < iNumExcludeChroms; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"reg expressions defining chroms to exclude: '%s'",pszExcludeChroms[Idx]); 

	for(Idx = 0; Idx < iNumIncludeLoci; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"files with loci to include: '%s'",pszIncludeLoci[Idx]);
	for(Idx = 0; Idx < iNumExcludeLoci; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"files with loci to include: '%s'",pszExcludeLoci[Idx]); 


#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	// processing here...
	Rslt = Process((etProcMode)iProcMode,			// processing mode
			(teCSVFormat)iCSVFormat,	// expected input CSV format
			szInLociFile,				// input loci file
			szInSeqFile,				// file containing genome assembly sequences
			(teRsltsFormat)iRsltsFormat,	// results format 
			szRsltsFile,				// write results to this file
			iRegionsIn,					// regions to retain
			iRegionsOut,				// regions to remove
			szFilterOutRefIDFile,		// RefIDs to filter out
			szFilterInRefIDFile,		// RefIDs to filter in unless filtered out
			iMinLength,					// elements must be of at least this length
			iMaxLength,					// no longer than this length
			iMaxNMerLen,				// max length N-Mer to process
			iMinNMerLen,				// min N-Mer to process
			iNumIncludeChroms,			// number of chromosome regular expressions to include
			pszIncludeChroms,			// array of include chromosome regular expressions
			iNumExcludeChroms,			// number of chromosome expressions to exclude
			pszExcludeChroms,			// array of exclude chromosome regular expressions
			iNumIncludeLoci,			// number of include region files 
			pszIncludeLoci,				// array of include region files
			iNumExcludeLoci,			// number of exclude region files
			pszExcludeLoci);			// array of exclude region files

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
}

char *
RsltsFormat2Text(teRsltsFormat Format)
{
switch(Format) {
	case eRsltsFasta:	// all sequences concatenated as single fasta record
		return((char *)"Sequences concatenated into single Fasta record");
	case eRsltsMultifasta:	// sequences as multifasta
		return((char *)"Sequences as multiple fasta records");
	default:
		break;
	}
return((char *)"Unsupported");
}

char *
CSVFormat2Text(teCSVFormat Format)
{
switch(Format) {
	case eCSVFdefault:
		return((char *)"Default 9+ field loci");
	case eCSVFprobe:
		return((char *)"locsfx probe loci");
	case eCSVFtarget:
		return((char *)"locsfx target loci");
	default:
		break;
	}
return((char *)"Unsupported");
}

char *
ProcMode2Txt(etProcMode ProcMode)
{
switch(ProcMode) {
	case eProcModeNMerDistAllSeqs:
		return((char *)"process N-Mer frequency distribution over all sequences");
	case eProcModeNMerDistPerSeq:
		return((char *)"process per sequence N-Mer frequency distribution");

	default:
		break;
	}
return((char *)"Unsupported");
}

// ParseRegions
// Parses space or comma delimited list of regions in which
// 1 == Intergenic, 2 == US, 3 == 5'UTR, 4 == CDS, 5 == Intron, 6 == 3'UTR, 7 == DS, 8 == 5'Splice, 9 == 3'Splice
//
// Returns bitmap of regions or -1 if parse errors
int
ParseRegions(char *pszRegionList)
{
// parse out region list
char Chr;
int Region = 0;
if(pszRegionList == NULL || *pszRegionList == '\0')
	return(-1);

while(Chr = *pszRegionList++) {
	if(isspace(Chr) || Chr == ',')		// accept spaces and commas as separators
		continue;
	switch(Chr) {
		case '1':						// intergenic to be filtered
			Region |= cFeatBitIG;
			break;
		case '2':						// 5'US to be filtered
			Region |= cFeatBitUpstream;
			break;
		case '3':						// 5'UTR to be filtered
			Region |= cFeatBit5UTR;
			break;
		case '4':
			Region |= cFeatBitCDS;		// CDS to be filtered
			break;
		case '5':
			Region |=  cFeatBitIntrons;	// any intronic to be filtered
			break;
		case '6':
			Region |=  cFeatBit3UTR;	// any 3'UTR to be filtered
			break;
		case '7':
			Region |=  cFeatBitDnstream;	// any 3'DS to be filtered 	
			break;
		case '8':
			Region |=  cIntronExonSpliceSite;	// any 5' splice site to be filtered
			break;
		case '9':
			Region |=  cExonIntronSpliceSite; // any 3' splice site to be filtered
			break;
		default:
			return(-1);
		}
	}
return(Region);
}

// Regions2Txt
// Returns textual representation of regions
char *
Regions2Txt(int Regions)
{
static char szRegions[200];
if(!Regions)
	return((char *)"None specified");
if(Regions & cFeatBitIG || Regions == 0)
	strcpy(szRegions,"Intergenic");
else
	szRegions[0] = '\0';
if(Regions & cFeatBitUpstream)
	{
	if(szRegions[0] != '\0')
		strcat(szRegions,", ");
	strcat(szRegions,"5'US");
	}
if(Regions & cFeatBit5UTR)
	{
	if(szRegions[0] != '\0')
		strcat(szRegions,", ");
	strcat(szRegions,"5'UTR");
	}
if(Regions & cFeatBitCDS)
	{
	if(szRegions[0] != '\0')
		strcat(szRegions,", ");
	strcat(szRegions,"CDS");
	}
if(Regions & cFeatBitIntrons)
	{
	if(szRegions[0] != '\0')
		strcat(szRegions,", ");
	strcat(szRegions,"Introns");
	}
if(Regions & cFeatBit3UTR)
	{
	if(szRegions[0] != '\0')
		strcat(szRegions,", ");
	strcat(szRegions,"3'UTR");
	}
if(Regions & cFeatBitDnstream)
	{
	if(szRegions[0] != '\0')
		strcat(szRegions,", ");
	strcat(szRegions,"3'DS");
	}
if(Regions & cIntronExonSpliceSite)
	{
	if(szRegions[0] != '\0')
		strcat(szRegions,", ");
	strcat(szRegions,"3'Splice");
	}
return(szRegions);
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


void
Cleanup(tsProcParams *pProcParams)
{
if(pProcParams->pFilterInRefIDs != NULL)
	{
	delete pProcParams->pFilterInRefIDs;
	pProcParams->pFilterInRefIDs = NULL;
	}

if(pProcParams->pFilterOutRefIDs != NULL)
	{
	delete pProcParams->pFilterOutRefIDs;
	pProcParams->pFilterOutRefIDs = NULL;
	}

if(pProcParams->pIncludeFiltLoci != NULL)
	{
	delete pProcParams->pIncludeFiltLoci;
	pProcParams->pIncludeFiltLoci = NULL;
	}

if(pProcParams->pExcludeFiltLoci != NULL)
	{
	delete pProcParams->pExcludeFiltLoci;
	pProcParams->pExcludeFiltLoci = NULL;
	}

if(pProcParams->pBioseq != NULL)
	{
	delete pProcParams->pBioseq;
	pProcParams->pBioseq = NULL;
	}

if(pProcParams->pSeq != NULL)
	{
	delete pProcParams->pSeq;
	pProcParams->pSeq = NULL;
	}

if(pProcParams->pCSV != NULL)
	{
	delete pProcParams->pCSV;
	pProcParams->pCSV = NULL;
	}

if(pProcParams->hRsltsFile != -1)
	{
	close(pProcParams->hRsltsFile);
	pProcParams->hRsltsFile = -1;
	}

if(pProcParams->pElChroms != NULL)
	{
	delete pProcParams->pElChroms;
	pProcParams->pElChroms = NULL;
	}

if(pProcParams->pCntStepCnts != NULL)
	{
	delete pProcParams->pCntStepCnts;
	pProcParams->pCntStepCnts = NULL;
	}
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
// Locates - starting from the MRA - a tsExcludeEl which matches on szSpecies and szChrom
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

// check if this chromosome are already known to be included/excluded
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

return(AddExcludeHistory(pszChrom,pProcParams->NumIncludeChroms > 0 ? true : false ,pProcParams));
}

int
GenSeqIdx(int SeqLen,etSeqBase *pSeq)
{
int Idx;
int SeqIdx;
int Base;
for(Idx=SeqIdx=0; Idx < SeqLen; Idx++,pSeq++)
	{
	Base = *pSeq & ~cRptMskFlg;
	if(Base > eBaseT)
		return(-1);
	SeqIdx <<= 2;
	SeqIdx |= Base;
	}
return(SeqIdx);
}

char *
StepIdx2Seq(int SeqLen,int SeqIdx)
{
static char szSeqBuff[128];
char *pChr;
int Base;
int Idx;

szSeqBuff[SeqLen] = '\0';
pChr = &szSeqBuff[SeqLen-1];
for(Idx=0; Idx < SeqLen; Idx++,pChr--)
	{
	Base = SeqIdx & 0x03;
	switch(Base) {
		case 0: *pChr = 'a'; break;
		case 1: *pChr = 'c'; break;
		case 2: *pChr = 'g'; break;
		case 3: *pChr = 't'; break;
		}
	SeqIdx >>= 2;
	}
return(szSeqBuff);
}



int Process(etProcMode ProcMode,		// processing mode
			teCSVFormat CSVFormat,		// expected input CSV format
			char *pszInLociFile,		// CSV file containing element loci
			char *pszInSeqFile,			// file containing genome assembly sequences
			teRsltsFormat RsltsFormat,	// results format 
			char *pszRsltsFile,			// write results to this file
			int FilterRegionsIn,		// retain any of these (exclusive) regions
			int FilterRegionsOut,		// remove any of these regions
			char *pszFilterOutRefIDFile, // RefIDs to filter out
			char *pszFilterInRefIDFile,	 // RefIDs to filter in unless filtered out
			int MinElLen,				 // elements must be of at least this length
			int MaxElLen,				 // no longer than this length
			int MaxNMerLen,				// max length N-Mer to process
			int MinNMerLen,				// min N-Mer to process
			int NumIncludeChroms,		 // number of chromosome regular expressions to include
			char **ppszIncludeChroms,	 // array of include chromosome regular expressions
			int NumExcludeChroms,		 // number of chromosome expressions to exclude
			char **ppszExcludeChroms,	 // array of exclude chromosome regular expressions
			int NumIncludeLoci,			 // number of include region files 
			char **ppszIncludeLoci,		 // array of include region files
			int NumExcludeLoci,			 // number of exclude region files
			char **ppszExcludeLoci)		// array of exclude region files
{
int Idx;
int Rslt;
tsProcParams ProcParams;

#ifndef _WIN32
int RegErr;				// regular expression parsing error
char szRegErr[128];			// to hold RegErr as textual representation ==> regerror();
#endif

memset(&ProcParams,0,sizeof(ProcParams));
ProcParams.ProcMode = ProcMode;
ProcParams.CSVFormat = CSVFormat;
ProcParams.RsltsFormat = RsltsFormat;
ProcParams.hRsltsFile = -1;
ProcParams.RegionsIn = FilterRegionsIn;
ProcParams.RegionsOut= FilterRegionsOut;
ProcParams.MaxElLen = MaxElLen;
ProcParams.MinElLen = MinElLen;
ProcParams.MinNMerLen = MinNMerLen;
ProcParams.MaxNMerLen = MaxNMerLen;
ProcParams.pszRsltsFile = pszRsltsFile;
ProcParams.pszInLociFile = pszInLociFile;
ProcParams.pszInSeqFile = pszInSeqFile;

if(NumIncludeLoci)
	{
	if((ProcParams.pIncludeFiltLoci = new CFilterLoci) == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instance CFilterLoci");
		return(-1);
		}
	for(Idx=0;Idx<NumIncludeLoci; Idx++)
		{
		if((Rslt = ProcParams.pIncludeFiltLoci->Load(ppszIncludeLoci[Idx])) < 0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open/load %s",ppszIncludeLoci[Idx]);
			Cleanup(&ProcParams);
			return(Rslt);
			}
		}
	}
if(NumExcludeLoci)
	{
	if((ProcParams.pExcludeFiltLoci = new CFilterLoci) == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instance CFilterLoci");
		Cleanup(&ProcParams);
		return(-1);
		}
	for(Idx=0;Idx<NumExcludeLoci; Idx++)
		{
		if((Rslt = ProcParams.pExcludeFiltLoci->Load(ppszExcludeLoci[Idx]))<0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open/load %s",ppszExcludeLoci[Idx]);
			Cleanup(&ProcParams);
			return(Rslt);
			}
		}
	}

if(pszFilterInRefIDFile != NULL && pszFilterInRefIDFile[0])
	{
	ProcParams.pFilterInRefIDs = new CFilterRefIDs;
	if((Rslt=ProcParams.pFilterInRefIDs->Open(pszFilterInRefIDFile)) < 0)
		{
		while(ProcParams.pFilterInRefIDs->NumErrMsgs())
			gDiagnostics.DiagOut(eDLFatal,gszProcName,ProcParams.pFilterInRefIDs->GetErrMsg());
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open/parse input filter RefID file %s",pszFilterInRefIDFile);
		Cleanup(&ProcParams);
		return(Rslt);
		}
	}

if(pszFilterOutRefIDFile != NULL && pszFilterOutRefIDFile[0])
	{
	ProcParams.pFilterOutRefIDs = new CFilterRefIDs;
	if((Rslt=ProcParams.pFilterOutRefIDs->Open(pszFilterOutRefIDFile)) < 0)
		{
		while(ProcParams.pFilterOutRefIDs->NumErrMsgs())
			gDiagnostics.DiagOut(eDLFatal,gszProcName,ProcParams.pFilterOutRefIDs->GetErrMsg());
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open/parse input filter RefID file %s",pszFilterOutRefIDFile);
		Cleanup(&ProcParams);
		return(Rslt);
		}
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

ProcParams.pCSV = new CCSVFile;
if(ProcParams.pCSV == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CCSVfile");
	Cleanup(&ProcParams);
	return(eBSFerrObj);
	}

if(pszInLociFile[0] != '\0')
	{
	if((Rslt=ProcParams.pCSV->Open(pszInLociFile))!=eBSFSuccess)
		{
		while(ProcParams.pCSV->NumErrMsgs())
			gDiagnostics.DiagOut(eDLFatal,gszProcName,ProcParams.pCSV->GetErrMsg());
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open file: %s",pszInLociFile);
		Cleanup(&ProcParams);
		return(Rslt);
		}
	}


if((ProcParams.pBioseq = new CBioSeqFile) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CBioSeqFile object");
	Cleanup(&ProcParams);
	return(eBSFerrObj);
	}

if((Rslt = ProcParams.pBioseq->Open(pszInSeqFile))!=eBSFSuccess)
	{
	while(ProcParams.pBioseq->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,ProcParams.pBioseq->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open assembly sequence file '%s'",pszInSeqFile);
	Cleanup(&ProcParams);
	return(Rslt);
	}

if((ProcParams.pSeq = new unsigned char [cDfltAllocSeqLen])==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate memory (%d bytes) for bioseq assembly subsequence",cDfltAllocSeqLen);
	Cleanup(&ProcParams);
	return(eBSFerrMem);
	}
ProcParams.AllocdSeqLen = cDfltAllocSeqLen;

int NumCnts = 0;
int RegionCnts = 0;
for(Idx = 0; Idx < MaxNMerLen; Idx++)
	{
	ProcParams.CntStepOfs[Idx] = NumCnts;
	if(!Idx)
		RegionCnts = 4;
	else
		RegionCnts *= 4;
	NumCnts += RegionCnts;
	}
	
if((ProcParams.pCntStepCnts = new INT64[NumCnts])==NULL)
	{
	Cleanup(&ProcParams);
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"nUnable to allocate memory (%d bytes) for holding composition statistics",sizeof(int) * NumCnts);
	return(eBSFerrMem);
	}
memset(ProcParams.pCntStepCnts,0,NumCnts * sizeof(INT64));
ProcParams.NumCntSteps = NumCnts;

if(pszRsltsFile != NULL && pszRsltsFile[0] != '\0')
#ifdef _WIN32
	if((ProcParams.hRsltsFile = open(pszRsltsFile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
	if((ProcParams.hRsltsFile = open(pszRsltsFile, O_RDWR | O_CREAT | O_TRUNC, S_IREAD | S_IWRITE ))==-1)
#endif
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create output results file: %s - %s",pszRsltsFile,strerror(errno));
		Cleanup(&ProcParams);
		return(eBSFerrCreateFile);
		}


Rslt = FilterCSV(&ProcParams);

Cleanup(&ProcParams);
return(Rslt);
}


int 
FilterCSV(tsProcParams *pParams)
{
int Rslt = 0;
int NumFields;
int RefID;
char *pszType;
char *pszSpecies;
char *pszChrom;
tChromID ChromID;
int StartLoci;
int EndLoci;
int Region;
int Len;

int BuffLen;
char szLineBuff[0x07fff];

int NumElsRead;		// number of elements before filtering
int NumElsAccepted;	// number of elements accepted after filtering
int NumElsRejected; // number of elements rejected

tsFiltState FiltState; // filter state for current element

bool bUseCSV;
bUseCSV = pParams->pszInLociFile[0] == '\0' ? false : true;

int NumSeqEntries = pParams->pBioseq->NumEntries();
int EntryID = 1;
char szSeqName[128];
char szSeqTitle[128];

int Idx;

int OverLens = 0;
int UnderLens = 0;
int ExcludeLoci = 0;
int IncludeLoci = 0;
int IncludeRegions = 0;
int ExcludeRegions = 0;
int Chroms = 0;
int Loci= 0;
int OutRefIDs = 0;
int InRefIDs = 0;

NumElsRead =0;		// number of elements before filtering
NumElsAccepted =0;	// number of elements accepted after filtering
NumElsRejected = 0; // number of elements rejected
memset(pParams->pCntStepCnts,0,pParams->NumCntSteps * sizeof(INT64));
while((!bUseCSV && NumSeqEntries) || (bUseCSV && (Rslt=pParams->pCSV->NextLine()) > 0))	// onto next line containing fields
	{
	if(bUseCSV)
		{
		NumFields = pParams->pCSV->GetCurFields();
		if(NumFields < 7)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expected at least 7 fields in '%s', GetCurFields() returned '%d'",pParams->pszInLociFile,NumFields);
			return(eBSFerrFieldCnt);
			}

		NumElsRead += 1;
		memset(&FiltState,0,sizeof(FiltState));
		pParams->pCSV->GetInt(1,&RefID);
		if(pParams->pFilterOutRefIDs != NULL && pParams->pFilterOutRefIDs->Locate(RefID))
			{
			FiltState.fOutRefID = 1;
			FiltState.bFiltOut = true;
			}

		if(pParams->pFilterInRefIDs != NULL && !pParams->pFilterInRefIDs->Locate(RefID))
			{	
			FiltState.fInRefID = 1;
			FiltState.bFiltOut = true;
			}

		pParams->pCSV->GetInt(7,&Len);
		if(Len < pParams->MinElLen)
			{
			FiltState.fUnderLen = 1;
			FiltState.bFiltOut = true;
			}
		else
			if(Len > pParams->MaxElLen)
				{
				FiltState.fOverLen = 1;
				FiltState.bFiltOut = true;
				}

		pParams->pCSV->GetText(2,&pszType);
		pParams->pCSV->GetText(3,&pszSpecies);
		
		pParams->pCSV->GetText(4,&pszChrom);
		if(ExcludeThisChrom(pszChrom,pParams))
			{
			FiltState.fChrom = 1;
			FiltState.bFiltOut = true;
			}

		pParams->pCSV->GetInt(5,&StartLoci);
		pParams->pCSV->GetInt(6,&EndLoci);
		}
	else
		{
		NumElsRead += 1;
		memset(&FiltState,0,sizeof(FiltState));
		Len = pParams->pBioseq->GetDataLen(EntryID);
		pParams->pBioseq->GetName(EntryID,sizeof(szSeqName)-1,szSeqName);
		pszChrom = szSeqName; 
		pszType = (char *)"seq";
		pParams->pBioseq->GetTitle(sizeof(szSeqTitle)-1,szSeqTitle);
		pszSpecies = szSeqTitle;
		StartLoci = 0;
		EndLoci = Len - 1;
		EntryID += 1;
		NumSeqEntries -= 1;
		}

	if(pParams->pExcludeFiltLoci != NULL && pParams->pExcludeFiltLoci->Locate(pszChrom,StartLoci,EndLoci))
		{	
		FiltState.fExcludeLoci = 1;
		FiltState.bFiltOut = true;
		}

	if(pParams->pIncludeFiltLoci != NULL && !pParams->pIncludeFiltLoci->Locate(pszChrom,StartLoci,EndLoci))
		{	
		FiltState.fIncludeLoci = 1;
		FiltState.bFiltOut = true;
		}

	if(bUseCSV && (pParams->RegionsIn || pParams->RegionsOut))
		{
		if(NumFields < 9)				// need at least 9 fields to filter on regions
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expected at least 9 fields for region filtering in '%s', GetCurFields() returned '%d'",pParams->pszInLociFile,NumFields);
			return(eBSFerrFieldCnt);
			}
		pParams->pCSV->GetInt(9,&Region);
		if(pParams->RegionsOut)
			{
			if((Region == 0 && (pParams->RegionsOut & cFeatBitIG)) || Region & pParams->RegionsOut)	
				{
				FiltState.fExcludeRegion = 1;
				FiltState.bFiltOut = true;
				}
			}
		if(pParams->RegionsIn)
			{
			if(Region == 0 && !(pParams->RegionsIn & cFeatBitIG))	
				{
				FiltState.fIncludeRegion = 1;
				FiltState.bFiltOut = true;
				}
			else
				if(Region != 0)
					{
					int Msk = 0x01;
					for(Idx = 0; Idx < 6; Idx++,Msk <<= 1)
						{
						if((pParams->RegionsIn & Msk) && (Region & 0x03f) == Msk)
							break;
						}
					if(Idx == 6)
						{
						FiltState.fIncludeRegion = 1;
						FiltState.bFiltOut = true;
						}
					}
			}
		}

	if(!FiltState.bFiltOut)	// if this element is thus far unfiltered then get it's sequence
		{
		NumElsAccepted += 1;
		pParams->CurSeqLen = 0;
		if(pParams->pSeq == NULL || Len > pParams->AllocdSeqLen)
			{
			if(pParams->pSeq != NULL)
				delete pParams->pSeq;
			if((pParams->pSeq = new unsigned char [Len + cDfltAllocSeqLen])==NULL)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate memory (%d bytes) for bioseq assembly subsequence",Len + cDfltAllocSeqLen);
				return(eBSFerrMem);
				}
			pParams->AllocdSeqLen = Len + cDfltAllocSeqLen;
			}

		if((Rslt= ChromID = pParams->pBioseq->LocateEntryIDbyName(pszChrom))<=0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to locate chrom '%s' in assembly file '%s'",pszChrom,pParams->pszInSeqFile);
			return(eBSFerrChrom);
			}

		if((Rslt=pParams->pBioseq->GetData(ChromID,eSeqBaseType,StartLoci,pParams->pSeq,Len)) != Len)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Loading sequence failed from chrom: %s loci %d len: %d file: '%s'",pszChrom,StartLoci,Len,pParams->pszInSeqFile);
			return(eBSFerrOfs);
			}
		pParams->CurSeqLen = Len;
		pParams->pszCurChrom = pszChrom;
		pParams->StartLoci = StartLoci;
		pParams->EndLoci = EndLoci;
		pParams->CurRegion = Region; 

		if(pParams->ProcMode == eProcModeNMerDistPerSeq) 
			memset(pParams->pCntStepCnts,0,pParams->NumCntSteps * sizeof(INT64));
		if((Rslt=ProcessSequence(pParams)) < 0)
			return(Rslt);
		if(pParams->ProcMode == eProcModeNMerDistPerSeq && 
			(Rslt = OutputResults(pParams)) < 0)
			return(Rslt);
		}
	else
		{
		if(FiltState.fOverLen)
			OverLens+=1;
		if(FiltState.fUnderLen)
			UnderLens+=1;
		if(FiltState.fExcludeLoci)
			ExcludeLoci+=1;
		if(FiltState.fIncludeLoci)
			IncludeLoci+=1;
		if(FiltState.fChrom)
			Chroms+=1;
		if(FiltState.fOutRefID)
			OutRefIDs+=1;
		if(FiltState.fInRefID)
			InRefIDs+=1;
		if(FiltState.fIncludeRegion)
			IncludeRegions += 1;
		if(FiltState.fExcludeRegion)
			ExcludeRegions += 1;
		NumElsRejected += 1;
		}
	}

if(Rslt >= eBSFSuccess && pParams->ProcMode == eProcModeNMerDistAllSeqs) 
	Rslt = OutputResults(pParams);

BuffLen = sprintf(szLineBuff,"Processed: %d, Accepted: %d, Rejected: %d, Overlen: %d, UnderLen: %d, ExcludeLoci: %d, Include Loci: %d, Chroms: %d, OutRefIDs: %d, InRefIDs: %d, IncludeRegions: %d, ExcludeRegions: %d",
				  NumElsRead,NumElsAccepted,NumElsRejected,OverLens,UnderLens, ExcludeLoci,IncludeLoci,Chroms,OutRefIDs,InRefIDs,IncludeRegions,ExcludeRegions); 

gDiagnostics.DiagOut(eDLInfo,gszProcName,szLineBuff);
return(Rslt);
}

int
ProcessSequence(tsProcParams *pParams)
{
INT64 *pCnt;
int Idx;
int SeqIdx;
int CurSeqIdxLen;

for(Idx = 0; Idx < pParams->CurSeqLen; Idx++)
	{
	for(CurSeqIdxLen = pParams->MaxNMerLen; CurSeqIdxLen >= pParams->MinNMerLen; CurSeqIdxLen--)
		{
		if((Idx + CurSeqIdxLen) > pParams->CurSeqLen)
			continue;
		SeqIdx = GenSeqIdx(CurSeqIdxLen,&pParams->pSeq[Idx]);
		if(SeqIdx < 0)		// < 0 generally indicates that sequences contains 'N' or some other none a,c,g,t base
			continue;

		pCnt = pParams->pCntStepCnts;
		pCnt += pParams->CntStepOfs[CurSeqIdxLen-1];
		pCnt[SeqIdx] += 1;
		}
	}
return(eBSFSuccess);
}

int
OutputResults(tsProcParams *pParams)
{
char szLineBuff[cMaxReadLen];
INT64 *pCnt;

INT64 Tot;
int Len;
int NumSteps;
int StepIdx;
int Idx;

NumSteps = 0;
for(Idx = 0; Idx < pParams->MaxNMerLen; Idx++)
	{
	if(Idx == 0)
		NumSteps = 1;
	NumSteps *= 4;

	if(Idx < (pParams->MinNMerLen -1))
		continue;

	pCnt = &pParams->pCntStepCnts[pParams->CntStepOfs[Idx]];
	Tot = 0;
	for(StepIdx = 0; StepIdx < NumSteps; StepIdx++)
		Tot += *pCnt++;
	if(!Tot)
		continue;
	pCnt = &pParams->pCntStepCnts[pParams->CntStepOfs[Idx]];
	for(StepIdx = 0; StepIdx < NumSteps; StepIdx++,pCnt++)
		{
		if(!*pCnt)
			continue;
		Len = sprintf(szLineBuff,"\"%s\",%d,%d,\"%s\",%lld,%.10g",
			pParams->ProcMode == eProcModeNMerDistAllSeqs ? pParams->pszInLociFile : pParams->pszCurChrom,
				Idx+1,StepIdx+1,StepIdx2Seq(Idx+1,StepIdx),*pCnt,*pCnt/(double)Tot);
		Len += sprintf(&szLineBuff[Len],"\n");
		CUtility::SafeWrite(pParams->hRsltsFile,szLineBuff,Len);
		}
	}
return(eBSFSuccess);
}


