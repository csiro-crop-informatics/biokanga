// genzygosity.cpp : Defines the entry point for the console application.
// This process is expected to be manually renamed to be - 'kangaz' - when ready for release
//
// Purpose -
// Loads suffix array for targeted genome of interest and generates a score matrix whereby the scores represent the degree of zygosity
// between all chrom/contigs contained in the targeted gnome
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

#include "genzygosity.h"

const char *cpszProgVer = "0.0.4";		// increment with each release

// processing modes
typedef enum TAG_ePMode {
	ePMdefault,					// default processing mode
	ePMMoreSens,				// more sensitive - slower
	ePMUltraSens,				// ultra sensitive - much slower
	ePMLessSens,				// less sensitive - quicker
	ePMplaceholder				// used to set the enumeration range
	} etPMode;

// output format modes
typedef enum TAG_eFMode {
	eFMdefault,					// default is to only report those with a zygosity above threshold
	eFMall,						// report all
	eFMplaceholder				// used to set the enumeration range
	} etFMode;


#pragma pack(4)
typedef struct TAG_sThreadMatchPars {
	int ThreadIdx;				// index of this thread (1..m_NumThreads)
	int NumIdentNodes;			// number of ident nodes allocd for exlusive use by this thread
	tsIdentNode *pIdentNodes;	// thread to use these nodes

	UINT32 NumSfxEntries;		// number of entries in pLocEntryCnts and pBlockEntryCnts
	UINT32 *pLocEntryCnts;		// allocated to hold temp local counts for each suffix entry
	UINT32 *pEntryCnts;			// allocated to hold counts for each suffix entry

#ifdef _WIN32
	HANDLE threadHandle;			// handle as returned by _beginthreadex()
	unsigned int threadID;			// identifier as set by _beginthreadex()
#else
	int threadRslt;					// result as returned by pthread_create ()
	pthread_t threadID;				// identifier as set by pthread_create ()
#endif
	int CurBlockID;					// current suffix block identifier
 	int Rslt;						// returned result code
	int NumReadsProc;				// returned number of reads processed by this thread instance
} tsThreadMatchPars;


typedef struct TAG_sSubseqsToMatchBlock {
	int EntryID;			// suffix array entry identifier from which subsequences are derived
	UINT32 StartOfs;		// starting offset in entry at which these subsequences start
	int NumSubseqs;			// number of subsequences for processing in this block
} tsSubseqsToMatchBlock;
#pragma pack()


#ifdef _WIN32
HANDLE m_hMtxIterReads;
unsigned __stdcall ThreadedCoredApprox(void * pThreadPars);
#else
pthread_mutex_t m_hMtxIterReads;
void *ThreadedCoredApprox(void * pThreadPars);
#endif

void ResetThreadedIterSubseqs(void);
bool ThreadedIterSubseqs(tsSubseqsToMatchBlock *pRetBlock);		// returns false if no more reads availing for processing by calling thread


int
Process(etPMode PMode,					// processing mode
		int SubseqLen,					// subsequence length
		int MaxNs,					    // allow at most this number of indeterminate eBaseNs in read before deeming as nonalignable
		int MaxSubs,					// maximum number of substitutions allowed
		int MaxMatches,					// 0 if no limit, otherwise if more than this number of matches then filter out this match
		etFMode FMode,					// output format mode
		int NumThreads,					// number of worker threads to use
		double zygosityThreshold,		// only report those with a homozygosity equal or above this threshold
		char *pszRawRsltsFile,			// optionally report all raw to this file
		char *pszRsltsFile,				// where to write thresholded results
		char *pszSfxFile);				// target as suffix array

int LoadReads(char *pszRdsFile);			// load preprocessed reads (genreads output)

int LocateCoredSubseqs();



UINT32 ApproxNumSubseqsProc(void);		// gets an approximation of the number of reads thus far aligned

int
AppendStr(char *pszBuff,	// write to this buffer
		  char LeadSep,		// if != '\0' then prefix with this separator (usually ',' or '\t')
		  char Quote,		// if != '\0' then quote string with this char (usually single or double quote char)
		  char *pStr,		// '\0' terminated string
		  char TrailSep);	// if != '\0' then suffix with this separator (usually ',' or '\t' or '\n')

int
AppendChrs(char *pszBuff,	// write to this buffer
		  char LeadSep,		// if != '\0' then prefix with this separator (usually ',' or '\t')
		  char Quote,		// if != '\0' then quote chars with this char (usually single or double quote char)
		  int NumChrs,		// number of chars to append
		  char *Chrs,		// pts to chars to append
		  char TrailSep);	// if != '\0' then suffix with this separator (usually ',' or '\t' or '\n')

int							// length written
AppendUInt(char *pszBuff,	// write to this buffer
		  char LeadSep,		// if > '\0' then prefix with this separator (usually ',' or '\t')
		  UINT32 Value,
		  char TrailSep);	// if > '\0' then suffix with this separator (usually ',' or '\t' or '\n')



static int SortReadIDs(const void *arg1, const void *arg2);
static int SortHitMatch(const void *arg1, const void *arg2);

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
	return _T("kanga");
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

int iFileLogLevel;			// level of file diagnostics
int iScreenLogLevel;		// level of file diagnostics
char szLogFile[_MAX_PATH];	// write diagnostics to this file
int Rslt = 0;   			// function result code >= 0 represents success, < 0 on failure

etPMode PMode;				// processing mode
etFMode FMode;				// format output mode

int NumberOfProcessors;		// number of installed CPUs
int NumThreads;				// number of threads (0 defaults to number of CPUs)

int SubseqLen;				// subsequence sampling length
int MaxSubs;				// maximum number of substitutions allowed
int MaxNs;				    // allow at most this number of indeterminate eBaseNs in read before deeming as nonalignable
int MaxMatches;				// max allowed matches before filtering out any subsequence

double FiltZygosityThres;				// report zygosity equal or above this threshold
char szRsltsFile[_MAX_PATH];			// zygosity thresholded results to this file
char szRawRsltsFile[_MAX_PATH];			// all raw results to this file
char szTargFile[_MAX_PATH];				// align against this target suffix array genome file

// command line args
struct arg_lit  *help    = arg_lit0("h","help",                 "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to screen and logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"diagnostics log file");

struct arg_int *pmode = arg_int0("m","mode","<int>",		    "alignment processing mode: 0 - standard sensitivity, 1 - more sensitive (slower), 2 - ultra sensitive (slowest), 3 - less sensitive (quicker)");
struct arg_int *format = arg_int0("M","format","<int>",		    "output format: 0 - currently CSV only");

struct arg_dbl *filtzygositythres = arg_dbl0("z","zygosity","<double>","report zygosity proportion equal or above this threshold (defaults to 0.25)");

struct arg_file *sfxfile = arg_file1("i","sfx","<file>",		"process against this suffix array (kangax generated) file");
struct arg_file *rsltsfile = arg_file1("o","out","<file>",		"output results above threshold to this file");
struct arg_file *rawrsltsfile = arg_file0("O","rawrslts","<file>",	"optionally output all raw results to this file");

struct arg_int *subseqlen = arg_int0("l","subseqlen","<int>",	"subsequence sample length (default is 25)");
struct arg_int *maxsubs = arg_int0("s","substitutions","<int>",	"accept up to this number of aligner induced substitutions in alignments (default is 2)");
struct arg_int *maxns = arg_int0("n","maxns","<int>",	        "maximum number of indeterminate 'N's in subsequence before treating as unalignable (default is 1, max 5)");

struct arg_int *maxmatches = arg_int0("x","maxmatches","<int>",	"exclude subsequences having more matches than this limit - 0 if no limit (default is 5000)");

struct arg_int *threads = arg_int0("T","threads","<int>",		"number of processing threads 0..128 (defaults to 0 which sets threads to number of CPU cores)");
struct arg_end *end = arg_end(20);

void *argtable[] = {help,version,FileLogLevel,LogFile,
					pmode,subseqlen,maxsubs,maxns,maxmatches,filtzygositythres,format,sfxfile,rsltsfile,rawrsltsfile,threads,
					end};

char **pAllArgs;
int argerrors;
argerrors = CUtility::arg_parsefromfile(argc,(char **)argv,&pAllArgs);
if(argerrors >= 0)
	argerrors = arg_parse(argerrors,pAllArgs,argtable);

/* special case: '--help' takes precedence over error reporting */
if (help->count > 0)
        {
		printf("\n%s zygosity processor, Version %s\nOptions ---\n", gszProcName,cpszProgVer);
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

	PMode = (etPMode)(pmode->count ? pmode->ival[0] : ePMdefault);
	if(PMode < ePMdefault || PMode >= ePMplaceholder)
		{
		printf("\nError: Processing mode '-m%d' specified outside of range %d..%d",PMode,ePMdefault,(int)ePMplaceholder-1);
		exit(1);
		}

	FMode = (etFMode)(format->count ? format->ival[0] : eFMdefault);

	SubseqLen = subseqlen->count ? subseqlen->ival[0] : cDlftSubseqLen;
	if(SubseqLen < cMinSubseqLen || SubseqLen > cMaxSubseqLen)
		{
		printf("\nError: subsequence sample length  '-l%d' specified outside of range %d..%d",SubseqLen,cMinSubseqLen,cMaxSubseqLen);
		exit(1);
		}

	if(maxsubs->count)
		{
		MaxSubs = maxsubs->ival[0]; 

		// ensure user doesn't specify too many allowed substitutions relative to subsequence length
		int MaxSubs2Allow = min(cMaxAllowedSubs,SubseqLen - 18);
		if(MaxSubs > MaxSubs2Allow)
			{
			printf("\nError: Max allowed substitutions '-s%d' for subsequence length %d specified outside of range %d..%d",MaxSubs,SubseqLen,0,MaxSubs2Allow);
			exit(1);
			}
		}
	else
		MaxSubs = cDfltAllowedSubs;	

	MaxMatches = maxmatches->count ? maxmatches->ival[0] : cDfltMaxMatches;
	if(MaxMatches < 0 || MaxMatches > cMaxMaxMatches)
		{
		printf("\nError: max allowed matches for any subsequence  '-x%d' specified outside of range 0..%d",MaxMatches,cMaxMaxMatches);
		exit(1);
		}

	MaxNs = maxns->count ? maxns->ival[0] : cDfltMaxNs;
	if(MaxNs < 0 || MaxNs > cMaxNs)
		{
		printf("\nError: Allowed number of indeterminate bases in reads '-n%d' specified outside of range 0..%d",MaxNs,cMaxNs);
		exit(1);
		}

	FiltZygosityThres = filtzygositythres->count ? filtzygositythres->dval[0] : cDfltzygosityThreshold;
	if(FiltZygosityThres < 0.01f || FiltZygosityThres > 1.0f)
		{
		printf("\nError: Zygosity filter '-z%f' specified outside of range 0.01 to 1.0",FiltZygosityThres);
		exit(1);
		}

#ifdef _WIN32
	SYSTEM_INFO SystemInfo;
	GetSystemInfo(&SystemInfo);
	NumberOfProcessors = SystemInfo.dwNumberOfProcessors;
#else
	NumberOfProcessors = sysconf(_SC_NPROCESSORS_CONF);
#endif
	if((NumThreads = threads->count ? threads->ival[0] : NumberOfProcessors)==0)
		NumThreads = NumberOfProcessors;
	if(NumThreads < 0)
		{
		printf("\nError: Number of threads '-T%d' specified must be >= 0",NumThreads);
		exit(1);
		}
	if(NumThreads > NumberOfProcessors)				// silently truncate to be at most the actual number of processors
		NumThreads = NumberOfProcessors;

	FMode = (etFMode)(format->count ? format->ival[0] : eFMdefault);
	if(FMode < eFMdefault || FMode >= eFMplaceholder)
		{
		printf("\nError: Output format mode '-m%d' specified outside of range %d..%d",FMode,eFMdefault,(int)eFMplaceholder-1);
		exit(1);
		}

	strcpy(szTargFile,sfxfile->filename[0]);
	strcpy(szRsltsFile,rsltsfile->filename[0]);
	if(rawrsltsfile->count)
		{
		strncpy(szRawRsltsFile,rawrsltsfile->filename[0],_MAX_PATH);
		szRawRsltsFile[_MAX_PATH-1] = '\0';
		}
	else
		szRawRsltsFile[0] = '\0';


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
			pszDescr = "Standard alignment sensitivity";
			break;
		case ePMMoreSens:
			pszDescr = "More sensitive alignment (slower)";
			break;
		case ePMUltraSens:
			pszDescr = "Ultra sensitive alignment (very slow)";
			break;
		case ePMLessSens:
		default:
			pszDescr = "Less sensitive alignment (quicker)";
			break;
		}
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Processing mode is : '%s'",pszDescr);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"subsequence length : %d",SubseqLen);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"maximum aligner induced substitutions : %d",MaxSubs);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"filter out subsequences matching more than this limit : %d",MaxMatches);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"maximum number of indeterminate 'N's : %d",MaxNs);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"report zygosity equal or above : %f%%",FiltZygosityThres);

	switch(FMode) {
		case eFMdefault:
			pszDescr = "CSV only";
			break;
		}
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"output format is : '%s'",pszDescr);


	gDiagnostics.DiagOutMsgOnly(eDLInfo,"input target sequence(s) suffix array file: '%s'",szTargFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"output results file: '%s'",szRsltsFile);
	if(szRawRsltsFile[0] != '\0')
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"output raw results file: '%s'",szRawRsltsFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"number of threads : %d",NumThreads);


#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	gStopWatch.Start();
	Rslt = Process(PMode,SubseqLen,MaxNs,MaxSubs,MaxMatches,FMode,NumThreads,FiltZygosityThres,szRawRsltsFile,szRsltsFile,szTargFile);
	gStopWatch.Stop();
	Rslt = Rslt >=0 ? 0 : 1;
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Exit code: %d Total processing time: %s",Rslt,gStopWatch.Read());
	exit(Rslt);
	}
else
	{
	printf("\n%s the K-mer Adaptive Next Generation Aligner, Version %s\n",gszProcName,cpszProgVer);
	arg_print_errors(stdout,end,gszProcName);
	arg_print_syntax(stdout,argtable,"\nUse '-h' to view option and parameter usage\n");
	exit(1);
	}
return 0;
}


size_t m_BlockTotSeqLen;		// total sequence length in currently loaded suffix block
int	m_MinCoreLen;				// minimum core length allowed
int m_MaxIter;					// max iteration depth when matching cores

int m_PerThreadAllocdIdentNodes;    // each thread can use this many tsIdentNodes
int m_TotAllocdIdentNodes;			// total number of tsIdentNodes allocated
tsIdentNode *m_pAllocsIdentNodes;	// memory allocated to hold tsIdentNodes required by all threads
int m_TotAllocdEntryCnts;			// total number of allocated suffix entry counts
UINT32 *m_pAllocEntryCnts;			// memory allocated to hold each threads suffix entry counts 

UINT32 *m_pEntryCntsMatrix;			// memory allocated to hold accumulated entry counts for all targeted genome entries
UINT32 *m_pEntrySubseqCnts;			// memory allocated to hold accumulated number of entry unique subsequences

UINT32 m_SubseqLen;			// subsequences to be processed are of this length
int m_hInFile;				// input file handle
int m_hRsltsFile;			// results output file handle
int m_hRawRsltsFile;		// raw results output file

char *m_pszRsltsFile;		// results to this file
char *m_pszRawRsltsFile;	// raw results to this file
char *m_pszSfxFile;			// target as suffix array
UINT32 m_NumSfxEntries;		// suffix array contains this many entries (chroms/contigs etc)
UINT32 m_CurSfxEntryProc;	// current suffix entry being processed
UINT32 m_CurSfxEntryProcLen;	// length of current sfx array entry being processed
UINT32 m_CurSfxEntryStartOfs;	// offset into sfx array entry at which next block of subsequences are to be processed from
UINT32 m_TotSfxEntriesLen;	// total length of all sfx array entries

int m_szLineBuffIdx;	// offset into m_pszLineBuff at which to next write
char *m_pszLineBuff;	// allocated to hold output line buffering

CSfxArrayV3 *m_pSfxArray; // suffix array holds genome of interest
char m_szTargSpecies[cMaxDatasetSpeciesChrom+1]; // suffix array was generated over this targeted species

etPMode m_PMode;		// processing mode
etFMode m_FMode;		// output format mode
int m_NumThreads;		// number of worker threads to use
int m_MaxNs;			// max number of indeterminate bases 'N' to acept in read before deeming read as unalignable
int m_MaxSubs;			// maximum number of substitutions allowed
int m_MaxMatches;		// if more than 0 and encountering more matches than this limit then return match count -1
int m_MaxNumSlides;		// max number of times core can be slide to right
double ZygosityThreshold;	// only report those with a homozygosity equal or above this threshold

void
Init(void)
{
m_hInFile = -1;
m_hRsltsFile = -1;
m_hRawRsltsFile = -1;
m_pSfxArray = NULL;
m_pAllocsIdentNodes = NULL;
m_pAllocEntryCnts = NULL;
m_pEntryCntsMatrix = NULL;
m_pEntrySubseqCnts = NULL;
m_szLineBuffIdx = 0;
m_TotAllocdIdentNodes = 0;
m_PerThreadAllocdIdentNodes = 0;
m_TotAllocdEntryCnts = 0;
m_MaxSubs = 0;
m_MaxNs = 0;
m_MaxMatches = 0;
m_MaxIter = 0;
m_MaxNumSlides = 0;
m_pszRsltsFile = NULL;	
m_pszRawRsltsFile = NULL;
ZygosityThreshold = cDfltzygosityThreshold;
}

void
Reset(void)
{
if(m_hInFile != -1)
	{
	close(m_hInFile);
	m_hInFile = -1;
	}
if(m_hRsltsFile != -1)
	{
	close(m_hRsltsFile);
	m_hRsltsFile = -1;
	}
if(m_hRawRsltsFile != -1)
	{
	close(m_hRawRsltsFile);
	m_hRawRsltsFile = -1;
	}

if(m_pszLineBuff != NULL)
	{
	delete m_pszLineBuff;
	m_pszLineBuff = NULL;
	}
if(m_pAllocsIdentNodes != NULL)
	{
	delete m_pAllocsIdentNodes;
	m_pAllocsIdentNodes = NULL;
	}

if(m_pAllocEntryCnts != NULL)
	{
	delete m_pAllocEntryCnts;
	m_pAllocEntryCnts = NULL;
	}

if(m_pSfxArray != NULL)
	{
	delete m_pSfxArray;
	m_pSfxArray = NULL;
	}
if(m_pEntryCntsMatrix != NULL)
	{
	delete m_pEntryCntsMatrix;
	m_pEntryCntsMatrix = NULL;
	}
if(m_pEntrySubseqCnts != NULL)
	{
	delete m_pEntrySubseqCnts;
	m_pEntrySubseqCnts = NULL;
	}
Init();
}



int
CreateOrTruncResultFiles(void)
{

// create/truncate output file
#ifdef _WIN32
m_hRsltsFile = open(m_pszRsltsFile,O_CREATETRUNC );
#else
if((m_hRsltsFile = open(m_pszRsltsFile,O_RDWR | O_CREAT,S_IREAD | S_IWRITE))!=-1)
    if(ftruncate(m_hRsltsFile,0)!=0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate %s - %s",m_pszRsltsFile,strerror(errno));
			Reset();
			return(eBSFerrCreateFile);
			}
#endif

if(m_hRsltsFile < 0)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate output file '%s'",m_pszRsltsFile);
	Reset();
	return(eBSFerrCreateFile);
	}

if(m_pszRawRsltsFile != NULL && m_pszRawRsltsFile[0] != '\0')
	{
#ifdef _WIN32
	m_hRawRsltsFile = open(m_pszRawRsltsFile,O_CREATETRUNC );
#else
	if((m_hRawRsltsFile = open(m_pszRawRsltsFile,O_RDWR | O_CREAT,S_IREAD | S_IWRITE))!=-1)
		if(ftruncate(m_hRawRsltsFile,0)!=0)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate %s - %s",m_pszRawRsltsFile,strerror(errno));
				Reset();
				return(eBSFerrCreateFile);
				}
#endif

	if(m_hRawRsltsFile < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate output file '%s'",m_pszRawRsltsFile);
		Reset();
		return(eBSFerrCreateFile);
		}

	if((m_pszLineBuff = new char [cAllocLineBuffSize])==NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to allocate memory (%d bytes) for output line buffer",cAllocLineBuffSize);
		Reset();
		return(eBSFerrMem);
		}
	m_szLineBuffIdx = 0;
	}
else
	m_hRawRsltsFile = -1;

return(eBSFSuccess);
}

int
Process(etPMode PMode,					// processing mode
		int SubseqLen,					// subsequence length
		int MaxNs,					    // allow at most this number of indeterminate eBaseNs in read before deeming as nonalignable
		int MaxSubs,					// maximum number of substitutions allowed
		int MaxMatches,					// 0 if no limit, otherwise if more than this number of matches then filter out this match
		etFMode FMode,					// output format mode
		int NumThreads,					// number of worker threads to use
		double zygosityThreshold,		// only report those with a homozygosity equal or above this threshold
		char *pszRawRsltsFile,				// optionally report all raw to this file
		char *pszRsltsFile,				// where to write thresholded results
		char *pszSfxFile)				// target as suffix array
{
int Rslt;
int BuffLen = 0;
int BuffOfs = 0;
int SeqIdx;
Init();

m_PMode = PMode;
m_FMode = FMode;
m_NumThreads = NumThreads;
m_MaxSubs = MaxSubs;
m_MaxNs = MaxNs;
m_MaxMatches = MaxMatches;
m_SubseqLen = SubseqLen;
m_pszRsltsFile = pszRsltsFile;
m_pszRawRsltsFile = pszRawRsltsFile;
m_pszSfxFile = pszSfxFile;

// open bioseq file containing suffix array for targeted assembly to align reads against
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Loading suffix array file '%s'", pszSfxFile);
if((m_pSfxArray = new CSfxArrayV3()) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CSfxArrayV2");
	Reset();
	return(eBSFerrObj);
	}
if((Rslt=m_pSfxArray->Open(pszSfxFile,false,false,false))!=eBSFSuccess)
	{
	while(m_pSfxArray->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,m_pSfxArray->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open input bioseq suffix array file '%s'",pszSfxFile);
	Reset();
	return(Rslt);
	}

// report to user some sfx array metadata as conformation the targeted assembly is correct
strcpy(m_szTargSpecies,m_pSfxArray->GetDatasetName());
if(m_pSfxArray->IsSOLiD())
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Input bioseq suffix array file '%s' was generated for SOLiD colorspace, needs to have been basespace",pszSfxFile);
	Reset();
	return(-1);
	}

tsSfxHeaderV3 SfxHeader;
m_pSfxArray->GetSfxHeader(&SfxHeader);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Genome Assembly Name: '%s' Descr: '%s' Title: '%s' Version: %d",
					 m_szTargSpecies,SfxHeader.szDescription,SfxHeader.szTitle,SfxHeader.Version);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Assembly has blocks: %d, max block size: %llu",SfxHeader.NumSfxBlocks,SfxHeader.SfxBlockSize);

m_NumSfxEntries = m_pSfxArray->GetNumEntries();		// ensure that there is at least one entry in suffix array file
if(m_NumSfxEntries < 1)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Can't process input bioseq suffix array file '%s' as no entries",pszSfxFile);
	Reset();
	return(-1);
	}

// check that at least one entry has sequence of a length that one or more subsequences which can be processed from these entries
for(m_CurSfxEntryProc = 1; m_CurSfxEntryProc <= m_NumSfxEntries; m_CurSfxEntryProc++)
	{
	if(m_pSfxArray->GetSeqLen(m_CurSfxEntryProc) >= (UINT32)SubseqLen)
		break;
	}
if( m_CurSfxEntryProc > m_NumSfxEntries)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Input bioseq suffix array file '%s' entry lengths are all shorter than the requested subsequence length %d",pszSfxFile,SubseqLen);
	Reset();
	return(-1);
	}
m_CurSfxEntryProc = 0;

m_TotSfxEntriesLen = (UINT32)m_pSfxArray->GetTotSeqsLen();

// restrict the max core iterations according to the requested sensitivity
int MaxIter;
switch(PMode) {
	case ePMdefault:			// default processing mode
		MaxIter = cDfltSensCoreIters;
		break;
	case ePMMoreSens:			// more sensitive - slower
		MaxIter = cMoreSensCoreIters;
		break;
	case ePMUltraSens:			// ultra sensitive - much slower
		MaxIter = cUltraSensCoreIters;
		break;
	case ePMLessSens:			// less sensitive - quicker
	default:
		MaxIter = cMinSensCoreIters;
	}
m_MaxIter = MaxIter;
m_pSfxArray->SetMaxIter(m_MaxIter);

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing subsequences of length %d...",SubseqLen);

// iterate and allocate for matrix of all suffix array chroms or contigs

// reasonably confident that there will be a resultset generated now that the targeted assembly have been loaded
// so create/trunc required result files
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Creating/truncating result files..");
if((Rslt=CreateOrTruncResultFiles())!=eBSFSuccess)
	{
	Reset();
	return(Rslt);
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Creating/truncating result files completed");

// start processing
Rslt = LocateCoredSubseqs();
if(Rslt<eBSFSuccess)
	{
	Reset();
	return(Rslt);
	}

BuffLen = 0;
BuffOfs = 0;
SeqIdx;

// now start reporting the alignment results....

// now time to write out the read hits
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Reporting of aligned resultset started...");

char szRslts[4096];
int RsltsIdx;
double Zigosity;
char szSrcEntryName[100];
UINT32 SrcEntryID;
char szTargEntryName[100];
UINT32 TargEntryID;
UINT32 *pTargEntryCnts;
UINT32 *pSrcEntryCnts;

pSrcEntryCnts = m_pEntrySubseqCnts;
pTargEntryCnts = m_pEntryCntsMatrix;
m_szLineBuffIdx = 0;
RsltsIdx = 0;

for(SrcEntryID = 1; SrcEntryID <= m_NumSfxEntries; SrcEntryID++, pSrcEntryCnts++)
	{
	m_pSfxArray->GetIdentName(SrcEntryID,sizeof(szSrcEntryName),szSrcEntryName);
	for(TargEntryID = 1; TargEntryID <= m_NumSfxEntries; TargEntryID++, pTargEntryCnts++)
		{
		m_pSfxArray->GetIdentName(TargEntryID,sizeof(szSrcEntryName),szTargEntryName);
		if(*pSrcEntryCnts > 0 && ((Zigosity = ((double)*pTargEntryCnts / *pSrcEntryCnts)) >= zygosityThreshold))
			{
			RsltsIdx += sprintf(&szRslts[RsltsIdx],"\"%s\",%d,\"%s\",%d,%1.6f\n",szSrcEntryName,*pSrcEntryCnts,szTargEntryName,*pTargEntryCnts,Zigosity);
			if((RsltsIdx + 200) > sizeof(szRslts))
				{
				CUtility::SafeWrite(m_hRsltsFile,szRslts,RsltsIdx);
				RsltsIdx = 0;
				}
			}
		if(m_hRawRsltsFile != -1)
			{
			m_szLineBuffIdx += sprintf(&m_pszLineBuff[m_szLineBuffIdx],"\"%s\",%d,\"%s\",%d\n",szSrcEntryName,*pSrcEntryCnts,szTargEntryName,*pTargEntryCnts);
			if((m_szLineBuffIdx + 200) > cAllocLineBuffSize)
				{
				CUtility::SafeWrite(m_hRawRsltsFile,m_pszLineBuff,m_szLineBuffIdx);
				m_szLineBuffIdx = 0;
				}
			}
		}
	}
if(RsltsIdx > 0)
	{
	CUtility::SafeWrite(m_hRsltsFile,szRslts,RsltsIdx);
	RsltsIdx = 0;
	}

if(m_hRawRsltsFile != -1 && m_szLineBuffIdx > 0)
	{
	CUtility::SafeWrite(m_hRawRsltsFile,m_pszLineBuff,m_szLineBuffIdx);
	m_szLineBuffIdx = 0;
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Reporting completed");

Reset();
return(Rslt);
}


int
AppendStr(char *pszBuff,	// write to this buffer
		  char LeadSep,		// if != '\0' then prefix with this separator (usually ',' or '\t')
		  char Quote,		// if != '\0' then quote string with this char (usually single or double quote char)
		  char *pStr,		// '\0' terminated string
		  char TrailSep)	// if != '\0' then suffix with this separator (usually ',' or '\t' or '\n')
{
int Len = 0;
if(LeadSep != '\0')
	{
	*pszBuff++ = LeadSep;
	Len += 1;
	}
if(Quote != '\0')
	{
	*pszBuff++ = Quote;
	Len += 1;
	}
while(*pStr)
	{
	*pszBuff++ = *pStr++;
	Len += 1;
	}
if(Quote != '\0')
	{
	*pszBuff++ = Quote;
	Len += 1;
	}
if(TrailSep != '\0')
	{
	*pszBuff++ = TrailSep;
	Len += 1;
	}
*pszBuff = '\0';
return(Len);
}

int
AppendChrs(char *pszBuff,	// write to this buffer
		  char LeadSep,		// if != '\0' then prefix with this separator (usually ',' or '\t')
		  char Quote,		// if != '\0' then quote chars with this char (usually single or double quote char)
		  int NumChrs,		// number of chars to append
		  char *Chrs,		// pts to chars to append
		  char TrailSep)	// if != '\0' then suffix with this separator (usually ',' or '\t' or '\n')
{
int Len = 0;
if(LeadSep != '\0')
	{
	*pszBuff++ = LeadSep;
	Len += 1;
	}
if(Quote != '\0')
	{
	*pszBuff++ = Quote;
	Len += 1;
	}
if(NumChrs)
	{
	while(NumChrs--)
		{
		*pszBuff++ = *Chrs++;
		Len += 1;
		}
	}
if(Quote != '\0')
	{
	*pszBuff++ = Quote;
	Len += 1;
	}
if(TrailSep != '\0')
	{
	*pszBuff++ = TrailSep;
	Len += 1;
	}
*pszBuff = '\0';
return(Len);
}


// very fast version of uitoa
int							// length written
AppendUInt(char *pszBuff,	// write to this buffer
		  char LeadSep,		// if > '\0' then prefix with this separator (usually ',' or '\t')
		  UINT32 Value,
		  char TrailSep)	// if > '\0' then suffix with this separator (usually ',' or '\t' or '\n')
{
int Len = 0;
char *pChr;
char *pMark;
char Tmp;
if(LeadSep != '\0')
	{
	*pszBuff++ = LeadSep;
	Len += 1;
	}

if(Value)
	{
	pChr = pszBuff;
	while(Value)
		{
		*pChr++ = '0' + (char)(Value%10);
		Value/=10;
		Len += 1;
		}
	pMark = pChr;
	*pChr-- = '\0';
	while(pszBuff < pChr)
		{
		Tmp = *pChr;
		*pChr-- = *pszBuff;
		*pszBuff++ = Tmp;
		}
	pszBuff = pMark;
	}
else
	{
	Len += 1;
	*pszBuff++ = '0';
	}
if(TrailSep)
	{
	Len += 1;
	*pszBuff++ = TrailSep;
	}
*pszBuff = '\0';
return(Len);
}


// LocateCoredSubseqs
// Locates all cored subsequences
int
LocateCoredSubseqs(void)
{
int Rslt;
int CurBlockID;							// current suffix block being processed
tBSFEntryID CurChromID;				    // current suffix array entry being processed
UINT32 TotNumReadsProc;					// total number of reads processed

UINT32 CurReadsProc;
UINT32 PrevReadsProc;

int ThreadIdx;
tsThreadMatchPars WorkerThreads[cMaxWorkerThreads];

m_PerThreadAllocdIdentNodes = cMaxNumIdentNodes;
m_TotAllocdIdentNodes = m_PerThreadAllocdIdentNodes * m_NumThreads;
if((m_pAllocsIdentNodes = new tsIdentNode [m_TotAllocdIdentNodes])==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: unable to allocate memory for %d tsIdentNodes",m_TotAllocdIdentNodes);
	Reset();
	return(eBSFerrMem);
	}

m_TotAllocdEntryCnts = m_NumSfxEntries * m_NumThreads * 2;
if((m_pAllocEntryCnts = new UINT32 [m_TotAllocdEntryCnts])==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: unable to allocate memory for %d entry counts",m_TotAllocdEntryCnts);
	Reset();
	return(eBSFerrMem);
	}


if((m_pEntryCntsMatrix = new UINT32 [m_NumSfxEntries * m_NumSfxEntries])==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: unable to allocate memory for %d matrix entry counts",m_NumSfxEntries * m_NumSfxEntries);
	Reset();
	return(eBSFerrMem);
	}
memset(m_pEntryCntsMatrix,0,m_NumSfxEntries * m_NumSfxEntries * sizeof(UINT32));

if((m_pEntrySubseqCnts = new UINT32 [m_NumSfxEntries])==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: unable to allocate memory for %d entry unique subseq counts",m_NumSfxEntries);
	Reset();
	return(eBSFerrMem);
	}
memset(m_pEntrySubseqCnts,0,m_NumSfxEntries * sizeof(UINT32));

#ifdef _WIN32
if((m_hMtxIterReads = CreateMutex(NULL,false,NULL))==NULL)
#else
if(pthread_mutex_init (&m_hMtxIterReads,NULL)!=0)
#endif
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: unable to create mutex");
	return(eBSFerrInternal);
	}

// load single SfxBlock, expected to contain all chromosomes, and process all reads against that block
CurChromID = 0;
CurBlockID = 1;
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Loading genome assembly suffix array...");
if((Rslt=m_pSfxArray->SetTargBlock(CurBlockID))<0)
	{
	while(m_pSfxArray->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,m_pSfxArray->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: unable to load genome assembly suffix array");
	return(Rslt);
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Genome assembly suffix array loaded");

	// determine minimum core length from targeted sequence length
	// core length is a balance between sensitivity and throughput
	// reducing core size has a relatively minor effect on sensitivity but significantly reduces throughput
	// large genomes require larger cores, more sensitive alignments require smaller cores
m_BlockTotSeqLen = m_pSfxArray->GetTotSeqsLen();

if(m_BlockTotSeqLen < 20000000)				    // covers yeast
	m_MinCoreLen = cMinCoreLen;
else
	if(m_BlockTotSeqLen < 250000000)		    // covers arabidopsis and fly
		m_MinCoreLen = cMinCoreLen+3;
	else
		m_MinCoreLen = cMinCoreLen+5;			// covers the big guys...

switch(m_PMode) {
	case ePMUltraSens:				// ultra sensitive - much slower
		m_MaxNumSlides = 8;			// leave m_MinCoreLen at it's minimum
		break;						
	case ePMMoreSens:				// more sensitive - slower
		m_MinCoreLen += 2;
		m_MaxNumSlides = 6;
		break;						
	case ePMdefault:				// default processing mode
		m_MinCoreLen += 4;			
		m_MaxNumSlides = 4;		
		break;
	case ePMLessSens:				// less sensitive - quicker
		m_MinCoreLen += 6;
		m_MaxNumSlides = 3;
		break;
	default:
		m_MaxNumSlides = 3;
		break;
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Now aligning with minimum core size of %d...\n",m_MinCoreLen);
#ifndef _WIN32
printf("Number of reads processed: %1.9u",0);
#endif
ResetThreadedIterSubseqs();
memset(WorkerThreads,0,sizeof(WorkerThreads));
for(ThreadIdx = 0; ThreadIdx < m_NumThreads; ThreadIdx++)
	{
	WorkerThreads[ThreadIdx].ThreadIdx = ThreadIdx + 1;
	WorkerThreads[ThreadIdx].NumIdentNodes = m_PerThreadAllocdIdentNodes;
	WorkerThreads[ThreadIdx].pIdentNodes = &m_pAllocsIdentNodes[m_PerThreadAllocdIdentNodes * ThreadIdx];
	WorkerThreads[ThreadIdx].CurBlockID = CurBlockID;
	WorkerThreads[ThreadIdx].NumSfxEntries = m_NumSfxEntries;
	WorkerThreads[ThreadIdx].pEntryCnts = &m_pAllocEntryCnts[m_NumSfxEntries * ThreadIdx * 2];
	WorkerThreads[ThreadIdx].pLocEntryCnts = &m_pAllocEntryCnts[(m_NumSfxEntries * ThreadIdx * 2) + m_NumSfxEntries];

#ifdef _WIN32
	WorkerThreads[ThreadIdx].threadHandle = (HANDLE)_beginthreadex(NULL,0x0fffff,ThreadedCoredApprox,&WorkerThreads[ThreadIdx],0,&WorkerThreads[ThreadIdx].threadID);
#else
	WorkerThreads[ThreadIdx].threadRslt =	pthread_create (&WorkerThreads[ThreadIdx].threadID , NULL , ThreadedCoredApprox , &WorkerThreads[ThreadIdx] );
#endif
	}

// wait for all threads to have completed
TotNumReadsProc = 0;
PrevReadsProc = 0;
for(ThreadIdx = 0; ThreadIdx < m_NumThreads; ThreadIdx++)
	{
#ifdef _WIN32
	while(WAIT_TIMEOUT == WaitForSingleObject( WorkerThreads[ThreadIdx].threadHandle, 60000 * 10))
		{
		CurReadsProc =	ApproxNumSubseqsProc();
		if(CurReadsProc > PrevReadsProc)
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"Progress: approx %u subsequences aligned",CurReadsProc);
		PrevReadsProc = CurReadsProc;
		}
	CloseHandle( WorkerThreads[ThreadIdx].threadHandle);
#else
	pthread_join(WorkerThreads[ThreadIdx].threadID,NULL);
#endif
	TotNumReadsProc += WorkerThreads[ThreadIdx].NumReadsProc;
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Progress: approx %u subsequences aligned",ApproxNumSubseqsProc());
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Progress: alignment of subsequences completed");

m_PerThreadAllocdIdentNodes = 0;
m_TotAllocdIdentNodes = 0;
if(m_pAllocsIdentNodes != NULL)
	{
	delete m_pAllocsIdentNodes;
	m_pAllocsIdentNodes = NULL;
	}

#ifdef _WIN32
CloseHandle(m_hMtxIterReads);
#else
pthread_mutex_destroy(&m_hMtxIterReads);
#endif

return(eBSFSuccess);
}

bool
AccumEntryMatchCnts(UINT32 SrcEntryID,UINT32 NumSrcSubseqs,UINT32 *pEntryCnts)
{
UINT32 Idx;
UINT32 MatrixIdx;
UINT32 *pColCnts;

if(m_pEntryCntsMatrix == NULL || m_pEntrySubseqCnts == NULL || SrcEntryID == 0 || SrcEntryID > m_NumSfxEntries)
	return(false);

MatrixIdx = (SrcEntryID - 1) * m_NumSfxEntries;
pColCnts = &m_pEntryCntsMatrix[MatrixIdx];

#ifdef _WIN32
WaitForSingleObject(m_hMtxIterReads,INFINITE);
#else
pthread_mutex_lock(&m_hMtxIterReads);
#endif
for(Idx = 0; Idx < m_NumSfxEntries; Idx++)
	*pColCnts++ += *pEntryCnts++;
m_pEntrySubseqCnts[SrcEntryID-1] += NumSrcSubseqs;
#ifdef _WIN32
ReleaseMutex(m_hMtxIterReads);
#else
pthread_mutex_unlock(&m_hMtxIterReads);
#endif
return(true);
}

#ifdef _WIN32
unsigned __stdcall ThreadedCoredApprox(void * pThreadPars)
#else
void *ThreadedCoredApprox(void * pThreadPars)
#endif
{
tsThreadMatchPars *pPars = (tsThreadMatchPars *)pThreadPars; // makes it easier not having to deal with casts!
etSeqBase Sequence[cMaxFastQSeqLen+1];	// to hold sequence (sans quality scores) for current subsequence

UINT32 SeqIdx;
int NumNs;
etSeqBase *pSeq;

tsSubseqsToMatchBlock SubseqsToMatchBlock;			// block of reads for this thread to process
int HitRslt;

// iterate each subsequence starting from the first
pPars->NumReadsProc = 0;
SubseqsToMatchBlock.NumSubseqs = 0;

UINT32 Idx;
int CurEntryID;					// current entry being processed
UINT32 EntrySourceSfxEntryCnts;	// number of subseqs in current source sfx entry which were unique to that sfx entry
UINT32 *pTmpAllEntryCount;		// temp ptr into pAllEntryCounts[]
UINT32 *pTmplocEntryCount;		// temp ptr into pLocEntryCounts[]
UINT32 MatchLoci;

EntrySourceSfxEntryCnts = 0;
CurEntryID = 0;
while(ThreadedIterSubseqs(&SubseqsToMatchBlock))	// get next block of subsequences to be processed
	{
	if(SubseqsToMatchBlock.EntryID != CurEntryID)
		{
		if(EntrySourceSfxEntryCnts > 0)
			{
			AccumEntryMatchCnts(CurEntryID,EntrySourceSfxEntryCnts,pPars->pEntryCnts);		
			EntrySourceSfxEntryCnts = 0;
			}

		CurEntryID = SubseqsToMatchBlock.EntryID;
		memset(pPars->pEntryCnts,0,sizeof(*pPars->pEntryCnts) * pPars->NumSfxEntries);
		}

	MatchLoci = SubseqsToMatchBlock.StartOfs;
	for(MatchLoci = SubseqsToMatchBlock.StartOfs; MatchLoci < (SubseqsToMatchBlock.StartOfs + SubseqsToMatchBlock.NumSubseqs); MatchLoci++)
		{
		// get subsequence
		m_pSfxArray->GetSeq(SubseqsToMatchBlock.EntryID,MatchLoci,Sequence,m_SubseqLen);
		pPars->NumReadsProc+=1;

		// ensure subsequence doesn't contain too many indeterminate bases
		pSeq = Sequence;
		NumNs = 0;
		for(SeqIdx = 0; SeqIdx < m_SubseqLen; SeqIdx++,pSeq++)
			{
			if((*pSeq = (*pSeq & 0x07)) > eBaseN)
				break;
			*pSeq &= ~cRptMskFlg;
			if(*pSeq == eBaseN)
				{
				if(++NumNs > m_MaxNs)
					break;
				}
			}
		if(SeqIdx != m_SubseqLen) // if too many 'N's then slough this subsequence...
			continue;

		// The window core length is set to be read length / (m_MaxSubs+1) 
		// The window core length is clamped to be at least m_MinCoreLen
		int CoreLen = max(m_MinCoreLen,(int)m_SubseqLen/(m_MaxSubs+1));

		// returns < 0 if errors, 0 if no matches or too many, otherwise the total match instances  
		HitRslt = m_pSfxArray->LocateAllNearMatches(m_MaxSubs,			// max number of mismatches allowed
						 CoreLen,						// core window length 
						 CoreLen,						// core window offset increment (1..n)
						 m_MaxNumSlides,				// max number of times to slide core on each strand
						 m_MaxMatches,				// if more than 0 and encountering more matches than this limit then return match count -1
 					 	 pPars->NumSfxEntries,			// number of entries in pEntryMatch
						 pPars->pLocEntryCnts,			// return number of matches to each entry in this array
						 CurEntryID,					// subsequence was from this entry
						 MatchLoci,						// starting at this offset
						 m_SubseqLen,					// probe length
						 Sequence,						// probe
						 m_MaxIter,						// max allowed iterations per subsegmented sequence when matching that subsegment
						 pPars->NumIdentNodes,			// memory has been allocated by caller for holding upto this many tsIdentNodes
						 pPars->pIdentNodes);			// memory allocated by caller for holding tsIdentNodes

		if(HitRslt < 0)		// < 0 if too many matches or sequence not unique in EntryID 
			continue;

		// Sequence was unique in entry
		if(HitRslt > 0)
			{
			// have matches from subsequence current entry into other chrom/contig entries  so accumulate the counts
			pTmpAllEntryCount = pPars->pEntryCnts;
			pTmplocEntryCount = pPars->pLocEntryCnts;
			for(Idx = 0; Idx < pPars->NumSfxEntries; Idx++,pTmpAllEntryCount++,pTmplocEntryCount++)
				{
				if(*pTmplocEntryCount > 0)
					*pTmpAllEntryCount += 1;
				}
			}
		EntrySourceSfxEntryCnts += 1;
		}
	}
if(CurEntryID > 0 &&  EntrySourceSfxEntryCnts > 0)
	AccumEntryMatchCnts(CurEntryID,EntrySourceSfxEntryCnts,pPars->pEntryCnts);					

pPars->Rslt = 1;
#ifdef _WIN32
_endthreadex(0);
return(eBSFSuccess);
#else
pthread_exit(NULL);
#endif
}



UINT32 m_NumSubseqsProc;	// number of subsequences thus far processed - note this is total subsequences handed out to processing threads
							// and should be treated as a guide only



unsigned long ProcessingStartSecs;

void
ResetThreadedIterSubseqs(void) // must be called by master thread prior to worker threads calling ThreadedIterSubseqs()
{
m_NumSubseqsProc = 0;
m_CurSfxEntryProc = 0;
m_CurSfxEntryProcLen = 0;
m_CurSfxEntryStartOfs = 0;
ProcessingStartSecs = gStopWatch.ReadUSecs();
}

// ApproxNumSubseqsProc
// Returns number of subsequences thus far returned to threads for processing 
// Only termed as approximate because processing may not have yet completed on all these subsequences
UINT32
ApproxNumSubseqsProc(void)
{
UINT32 NumProc;
#ifdef _WIN32
WaitForSingleObject(m_hMtxIterReads,INFINITE);
#else
pthread_mutex_lock(&m_hMtxIterReads);
#endif
NumProc = m_NumSubseqsProc;
#ifdef _WIN32
ReleaseMutex(m_hMtxIterReads);
#else
pthread_mutex_unlock(&m_hMtxIterReads);
#endif
return(NumProc);
}

bool						// returns false if no more reads available for processing by calling thread
ThreadedIterSubseqs(tsSubseqsToMatchBlock *pRetBlock)
{
UINT32 MaxSubSeqs2Proc;
pRetBlock->NumSubseqs = 0;
pRetBlock->EntryID = 0;
pRetBlock->StartOfs = 0;

#ifdef _WIN32
WaitForSingleObject(m_hMtxIterReads,INFINITE);
#else
pthread_mutex_lock(&m_hMtxIterReads);
#endif
if(m_pSfxArray == NULL)
	{
#ifdef _WIN32
	ReleaseMutex(m_hMtxIterReads);
#else
	pthread_mutex_unlock(&m_hMtxIterReads);
#endif
	return(false);
	}

if(m_CurSfxEntryProc == 0 || (m_CurSfxEntryStartOfs + m_SubseqLen) > m_CurSfxEntryProcLen)	
	{
	m_CurSfxEntryProcLen = 0;
	m_CurSfxEntryStartOfs = 0;
	m_CurSfxEntryProc += 1;

	while(m_CurSfxEntryProc <= m_NumSfxEntries && (m_CurSfxEntryProcLen = m_pSfxArray->GetSeqLen(m_CurSfxEntryProc)) > 0)
		{
		if(m_CurSfxEntryProcLen >= m_SubseqLen)
			break;
		}

	if(m_CurSfxEntryProcLen < 1 || m_CurSfxEntryProc > m_NumSfxEntries)
		{
#ifdef _WIN32
		ReleaseMutex(m_hMtxIterReads);
#else
		pthread_mutex_unlock(&m_hMtxIterReads);
#endif
		return(false);
		}
	}

// idea is to maximise the number of threads still processing when most subsequences have been processed so that
// the last thread processing doesn't end up with a large block of subsequences needing lengthly processing
MaxSubSeqs2Proc = min(cMaxSubseqsPerBlock,250 + ((m_TotSfxEntriesLen - m_NumSubseqsProc) / (UINT32)m_NumThreads));

pRetBlock->NumSubseqs = min(MaxSubSeqs2Proc,1 + m_CurSfxEntryProcLen - (m_CurSfxEntryStartOfs + m_SubseqLen));
pRetBlock->EntryID = m_CurSfxEntryProc;
pRetBlock->StartOfs = m_CurSfxEntryStartOfs;
m_CurSfxEntryStartOfs += pRetBlock->NumSubseqs;
m_NumSubseqsProc += pRetBlock->NumSubseqs;
#ifndef _WIN32
unsigned long ElapsedSecs = gStopWatch.ReadUSecs();
if((ElapsedSecs - ProcessingStartSecs) > 30)
	{
	ProcessingStartSecs = ElapsedSecs;
	printf("\b\b\b\b\b\b\b\b\b\b %1.9u",m_NumSubseqsProc);
	fflush(stdout);
	}
#endif

#ifdef _WIN32
ReleaseMutex(m_hMtxIterReads);
#else
pthread_mutex_unlock(&m_hMtxIterReads);
#endif
return(true);
}


