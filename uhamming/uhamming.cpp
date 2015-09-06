// uhamming.cpp : Defines the entry point for the console application.
// generates Hamming edit distances for all sequences of specified length over a target genome
// 1.4.2 increased Hamming max K-mer from previous 500, now 5000 

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

const char *cpszProgVer = "1.5.0";		// increment with each release

const int cMinSeqLen = 20;				// minimum sequence length for Hamming distances
const int cDfltSeqLen = 100;			// default sequence length for Hamming distances
const int cMaxSeqLen = 5000;			// maximum sequence K-mer length for exhustative Hamming distances

const int cMinRHamming = 1;				// restricted hamming lower limit
const int cDfltRHamming = 3;			// restricted hamming default limit
const int cMaxRHamming = 10;			// restricted hamming upper limit
const int cMaxNthRHamming = 100;		// sampling is to sample every Nth K-mer, with Nth is the range 1..cSampRHamming

const int cMinCoreLen = 8;				// restricted hamming minimum core length supported

const int cMaxWorkerThreads = 128;		// limiting max number of threads to this many
const int cMaxNumNodes = 10000;			// allow for upto this many nodes if processing is distributed over multiple nodes

const int cRptBuffAllocsize = 0x0fffff; // reporting buffer allocation size

#pragma pack(4)
typedef struct TAG_sWorkerThreads {
#ifdef _WIN32
	HANDLE threadHandle;			// handle as returned by _beginthreadex()
	unsigned int threadID;			// identifier as set by _beginthreadex()
#else
	int threadRslt;					// result as returned by pthread_create ()
	pthread_t threadID;				// identifier as set by pthread_create ()
#endif
	int Rslt;						// processing result
	int DummyValue;					// currently just a place holder
} tsWorkerThreads;

typedef struct TAG_sThreadParams {
	int ThreadID;		// uniquely (1..N) identifies this set of thread processing parameters
	bool bWatsonOnly;	// true if watson strand only (no crick revcpl processing)
	int SampleN;		// sample (process) every N sweep instance
	UINT32 SSofs;		// start processing with this relative (to first) subsequence (1..m_NumSubSeqs-1)
	UINT32 SeqDelta;	// process every SeqDelta subsequences
	UINT32 NumSeqs;	    // process, at most, this number of subsequences
	UINT64 HamDistOfs;  // offset into m_HamDist at which the edit distances for this processing instance are to be written
	int State;		    // processing state: 0 if waiting to be processed, 1 if currently being processed, 2 if processing completed
} tsThreadParams;
#pragma pack()

#pragma pack(1)
// Hamming specific structures
const int cMaxHammingChroms = 1000;		// can handle at most this many chromosomes with hammings
typedef struct TAG_sHamChrom {
	UINT32 ChromID;						// uniquely identifies this chromosome
	UINT8  szChrom[cMaxDatasetSpeciesChrom];	// chrom name
	UINT32 NumEls;						// number of subsequences with hammings on this chrom
	UINT16 Dists[1];					// array, in ascending loci order, of hamming distances
} tsHamChrom;

typedef struct TAG_sHamHdr {
	UINT8 Magic[4];		        // magic chars 'bham' to identify this file as a biosequence file containing hamming edit distances
	UINT32 Version;				// structure version
	INT32 Len;					// file length, also allocation size required when loading hammings into memory
	UINT16 NumChroms;		    // number of chromosomes with Hammings
	UINT32 ChromOfs[cMaxHammingChroms];	// offsets to each chromosomes respective tsHamChrom
} tsHamHdr;

const int cAllocHamDist = sizeof(tsHamHdr) + sizeof(tsHamChrom) + 10000000; // allocate for hamming distances in this sized chunks

#pragma pack()

// processing modes
typedef enum TAG_ePMode {
	ePMrestrict = 0,			// default processing mode is for restricted Hamming processing on a single node
	ePMnode,					// Single node exhaustive Hamming edit distance processing
	ePMdist,					// distributed processing over multiple nodes
	ePMmerge,					// merge multiple node Hammings
	ePMtrans,					// transform Hamming csv file into quick load binary
	ePMtransCSV,				// transform Hamming quick load binary file into csv
	ePMplaceholder				// used to set the enumeration range
	} etPMode;

// restricted Hamming output file format
typedef enum TAG_eResFormat {
	cRHFcsv = 0,				// CSV with Hamming loci ranges (chrom, startloci, length, Hamming)
	cRHFBed,					// UCSC BED file format
	cRHFWiggle,					// UCSC Wiggle file format
	cRHFplaceholder
} etResFormat;

// restricted Hamming processing sensitivity (determines core search depth)
typedef enum TAG_eSensitivity {
	eSensDefault,						// default processing mode
	eSensMore,							// more sensitive - slower
	eSensUltra,							// ultra sensitive - much slower
	eSensLess,							// less sensitive - quicker
	eSensPlaceholder					// used to set the enumeration range
} etSensitivity;


int
Process(etPMode PMode,			// processing mode
		bool bWatsonOnly,		// true if watson strand only processing
		etSensitivity Sensitivity, // restricted hamming processing sensitivity
		etResFormat ResFormat,	// restricted Hamming output file format
		int RHamm,			    // if > 0 then restricted hammings limit
		UINT32 SweepStart,		// start processing from this sweep instance (1..GenomeLen) or if distributed processing then the total number of nodes
		UINT32 SweepEnd,		// finish processing at this sweep instance (0 for all remaining or SweepStart..GenomeLen) or if distributed processing then the node instance
		int SeqLen,				// Hamming for these length sequences
		int SampleN,			// sample (process) every N sweep instance (or if restricted Hammings then every Nth K-mer) 
		int NumThreads,			// process with upto this number of threads
		char *pszGenomeFile,	// bioseq file containing targeted genome assembly, or if restricted hammings then the assembly suffix array file
		char *pszInSeqFile,		// if restricted Hammings then kmer sequences from this file
		char *pszHammingFile);	// writeout Hamming edit distances into this file

int LoadGenome(char *pszBioSeqFile); // load genome from this file
extern void GHamDistWatson(UINT16 *pHDs,		// where to return Hamming differentials for each subsequence
			  int SubSeqLen,	// generate Hammings edit distances for subsequences of this length
			  UINT32 SSofs,		// offset between subsequences for current pass
			  UINT8 *pGenomeSeq, // genome sequence (concatenated chrom seqs, separated by eBaseEOSs) with final chrom terminated by eBaseEOG, not eBaseEOS
			  UINT32 GenomeLen);	 // total genome length including chrom eBaseEOS markers, but exluding start/end eBaseEOG markers

extern void GHamDistCrick(UINT16 *pHDs,		// where to return Hamming differentials for each subsequence
			  int SubSeqLen,	// generate Hammings edit distances for subsequences of this length
			  UINT32 SSofs,		// offset between subsequences for current pass
			  UINT8 *pGenomeSeq, // genome sequence (concatenated chrom seqs, separated by eBaseEOSs) with final chrom terminated by eBaseEOG, not eBaseEOS
			  UINT32 GenomeLen);	 // total genome length including chrom eBaseEOS markers, but exluding start/end eBaseEOG markers

// use to return the processing parameters for use by calling thread instance
tsThreadParams *					// returns NULL if no more processing parameters available
ThreadedIterParams(void);

int MinHamCnt(int MaxHamCnt,		// only need to cnt upto this many mismatches
			   etSeqBase *pSeq1,	// determine Hamming edit distance between this sequence
		   etSeqBase *pSeq2);		// and this sequence

static int SortSeqs(const void *arg1, const void *arg2);


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
	return _T("ubsalign");
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

etPMode PMode;				// processing mode
etSensitivity Sensitivity;  // restricted hamming processing sensitivity

int NumberOfProcessors;		// number of installed Core CPUs
int NumThreads;				// number of threads (0 defaults to number of Core CPUs)
int SeqLen;					// Hamming edit distances for this length sequences
int SweepStart;			    // process starting from this sweep instance inclusive
int SweepEnd;				// complete processing at this sweep instance inclusive

int Node;					// node instance (1..N) if processing over multiple nodes
int NumNodes;				// total number of nodes (N) if processing over multiple nodes

int SampleN;				// sample every N sweep instances

bool bWatsonOnly;			// true if watson only strand processing - Crick is rather slow...
int CoreLen;				// core length to use when processing restricted maximal Hammings
int RHamm;					// restricted hamming upper limit
etResFormat ResFormat;				// resdtricted Hamming file output format
char szInFile[_MAX_PATH];	// input genome assembly file
char szInSeqFile[_MAX_PATH];	// optional input kmer sequences file
char szOutFile[_MAX_PATH];	// generated Hamming distances file


// command line args
struct arg_lit  *help    = arg_lit0("h","help",                 "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to screen and logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"diagnostics log file");

struct arg_int *pmode = arg_int0("m","mode","<int>",		    "distributed processing mode: 0 - restricted Hammings, 1 - exhaustive single node Hammings, 2 - exhaustive multiple node Hammings, 3 - merge multiple Hamming files, 4 - transform Hamming CSV into quick load binary format, 5  - transform quick load Hamming binary format into CSV (default = 0)");
struct arg_int *sensitivity = arg_int0("s","sensitivity","<int>","restricted Hamming sensitivity: 0 - normal, 1 - high, 2 - ultra, 3 - low (default = 0)");
struct arg_int *resformat = arg_int0("S","resformat","<int>",	"restricted Hamning file output format: 0 - csv, 1 - UCSC BED, 2 - UCSC Wiggle (default = 0)");


struct arg_lit  *crick = arg_lit0("c","strandcrick",            "process Crick in addition to Watson strand - Caution: very slow processing");

struct arg_int *rhamm = arg_int0("r","rhamm","<int>",			"restricted hamming upper limit (1..10, default 3) only applies in mode 0");

struct arg_int *numnodes = arg_int0("n","numnodes","<int>",	    "total number of nodes (2..10000) if processing over multiple nodes");
struct arg_int *node = arg_int0("N","node","<int>",	            "node instance (1..N) if processing over multiple nodes");

struct arg_int *sweepstart = arg_int0("b","sweepstart","<int>",	"process starting from this sweep instance inclusive (default = 1 for 1st)");
struct arg_int *sweepend = arg_int0("B","sweepend","<int>",		"complete processing at this sweep instance inclusive (default = 0 for all remaining, or >= Sweep start)");
struct arg_int *seqlen = arg_int0("K","seqlen","<int>",			"Hamming edit distances for these length k-mer subsequences (range 20..5000, default is 100)");
struct arg_file *infile = arg_file1("i","in","<file>",			"in mode 0 input sfx file, in mode 1 and 2, bioseq genome assembly file or in mode 3 merge from this input Hamming file");
struct arg_file *inseqfile = arg_file0("I","seq","<file>",		"if restricted hamming processing then optional file containing source kmer sequences");

struct arg_int *sample = arg_int0("k","sample","<int>",		    "sample every -S<N> sweep instances (default is 1) useful if only interested in overall distributions\n\t\tin restricted Hamming processing then sample every Nth (max 100) K-mer");

struct arg_file *outfile = arg_file0("o","out","<file>",		"output (merged) Hamming distances to this file");
struct arg_int *threads = arg_int0("T","threads","<int>",		"number of processing threads 0..128 (defaults to 0 which sets threads to number of CPU cores)");
struct arg_end *end = arg_end(20);

void *argtable[] = {help,version,FileLogLevel,LogFile,
					pmode,sensitivity,rhamm,resformat,crick,numnodes,node,sweepstart,sweepend,seqlen,sample,infile,inseqfile,outfile,threads,
					end};

char **pAllArgs;
int argerrors;
argerrors = CUtility::arg_parsefromfile(argc,(char **)argv,&pAllArgs);
if(argerrors >= 0)
	argerrors = arg_parse(argerrors,pAllArgs,argtable);

/* special case: '--help' takes precedence over error reporting */
if (help->count > 0)
        {
		printf("\n%s the K-mer Hamming distance generator, Version %s\nOptions ---\n", gszProcName,cpszProgVer);
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
		printf("\nError: FileLogLevel '-f%d' specified but no logfile '-F<logfile>\n'",FileLogLevel->ival[0]);
		exit(1);
		}

	iScreenLogLevel = iFileLogLevel = FileLogLevel->count ? FileLogLevel->ival[0] : eDLInfo;
	if(iFileLogLevel < eDLNone || iFileLogLevel > eDLDebug)
		{
		printf("\nError: FileLogLevel '-l%d' specified outside of range %d..%d\n",iFileLogLevel,eDLNone,eDLDebug);
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

	// now that log parameters have been parsed then initialise diagnostics log system
	if(!gDiagnostics.Open(szLogFile,(etDiagLevel)iScreenLogLevel,(etDiagLevel)iFileLogLevel,true))
		{
		printf("\nError: Unable to start diagnostics subsystem\n");
		if(szLogFile[0] != '\0')
			printf(" Most likely cause is that logfile '%s' can't be opened/created\n",szLogFile);
		exit(1);
		}

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Version: %s",cpszProgVer);


	PMode = (etPMode)(pmode->count ? pmode->ival[0] : ePMrestrict);
	if(PMode < ePMrestrict || PMode >= ePMplaceholder)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Processing mode '-m%d' specified outside of range %d..%d",PMode,0,(int)ePMplaceholder-1);
		exit(1);
		}


	SeqLen = cDfltSeqLen;
	NumNodes = 1;
	SweepStart = 0;
	SweepEnd = 0;
	NumThreads = 1;
	SampleN = 1;
	CoreLen = 0;
	RHamm = 0;
	ResFormat = cRHFcsv;
	Sensitivity = eSensDefault;
	szInSeqFile[0] = '\0';
	bWatsonOnly = true;

	if(PMode <= ePMdist)
		{
		SeqLen = seqlen->count ? seqlen->ival[0] : cDfltSeqLen;
		if(SeqLen < cMinSeqLen || SeqLen > cMaxSeqLen)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: k-mer sequence length '-K%d' specified outside of range %d..%d",SeqLen,cMinSeqLen,cMaxSeqLen);
			exit(1);
			}
		}
	else
		SeqLen = 0;

	if(PMode == ePMrestrict)
		{
		Sensitivity = (etSensitivity)(sensitivity->count ? sensitivity->ival[0] : eSensDefault);
		if(Sensitivity < eSensDefault || Sensitivity >= eSensPlaceholder)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Restricted hamming sensitivity '-s%d' specified outside of range %d..%d",Sensitivity,eSensDefault,eSensPlaceholder-1);
			exit(1);
			}

		RHamm= rhamm->count ? rhamm->ival[0] : cDfltRHamming;
		if(RHamm < cMinRHamming || RHamm > cMaxRHamming)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Restricted Hamming limit '-r%d' specified outside of range %d..%d",RHamm,cMinRHamming,cMaxRHamming);
			exit(1);
			}



		ResFormat= (etResFormat)(resformat->count ? resformat->ival[0] : cRHFcsv);
		if(ResFormat < cRHFcsv || ResFormat > cRHFWiggle)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Restricted Hamming output file format '-S%d' specified outside of range %d..%d",ResFormat,cRHFcsv,cRHFWiggle);
			exit(1);
			}

		CoreLen = SeqLen/(RHamm+1);
		if(CoreLen < cMinCoreLen)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Restricted hamming limit '-r%d' is incompatible with k-mer sequence length '-k%d'",RHamm,SeqLen);
			exit(1);
			}

		SampleN = sample->count ? sample->ival[0] : 1;
		if(SampleN < 1 || SampleN > cMaxNthRHamming)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Sampling '-S%d' specified outside of range 1..%d",SampleN,cMaxNthRHamming);
			exit(1);
			}

		if(crick->count > 0)
			bWatsonOnly = false;

		if(inseqfile->count)
			{
			strncpy(szInSeqFile,inseqfile->filename[0],_MAX_PATH);
			szInSeqFile[_MAX_PATH-1] = '\0';
			}

		}

	strcpy(szInFile,infile->filename[0]);

	if(PMode <= ePMdist)
		{
#ifdef _WIN32
		SYSTEM_INFO SystemInfo;
		GetSystemInfo(&SystemInfo);
		NumberOfProcessors = SystemInfo.dwNumberOfProcessors;
#else
		NumberOfProcessors = sysconf(_SC_NPROCESSORS_CONF);
#endif
		int MaxAllowedThreads = min(cMaxWorkerThreads,NumberOfProcessors);	// limit to be at most cMaxWorkerThreads
		if((NumThreads = threads->count ? threads->ival[0] : MaxAllowedThreads)==0)
			NumThreads = MaxAllowedThreads;
		if(NumThreads < 0 || NumThreads > MaxAllowedThreads)
			{
			gDiagnostics.DiagOut(eDLWarn,gszProcName,"Warning: Number of threads '-T%d' specified was outside of range %d..%d",NumThreads,1,MaxAllowedThreads);
			gDiagnostics.DiagOut(eDLWarn,gszProcName,"Warning: Defaulting number of threads to %d",MaxAllowedThreads);
			NumThreads = MaxAllowedThreads;
			}
		}


	if(PMode == ePMnode || PMode == ePMdist)
		{
		if(PMode == ePMnode)
			{
			SweepStart = sweepstart->count ? sweepstart->ival[0] : 1;
			if(SweepStart < 0)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Sweep start '-b%d' must be >= 1",SweepStart);
				exit(1);
				}
			if(SweepStart == 0)	// allow a little latitude, treat 0 as being the same as 1...
				SweepStart = 1;
			SweepEnd = sweepend->count ? sweepend->ival[0] : 0;
			if(SweepEnd != 0 && SweepEnd < SweepStart)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Sweep end '-B%d' must be either 0 or >= %d",SweepEnd,SweepStart);
				exit(1);
				}
			}

		if(PMode == ePMdist)
			{
			if(!numnodes->count || !node->count)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: In distributed processing mode both number of nodes '-n<num>' and node instance '-N<node>' must be specified");
				exit(1);
				}
			NumNodes = numnodes->ival[0];
			if(NumNodes < 2 || NumNodes > cMaxNumNodes)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Number of nodes '-n%d' must be in the range 2..%d",NumNodes,cMaxNumNodes);
				exit(1);
				}

			Node = node->ival[0];
			if(Node < 1)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: node instance '-N%d' must be in range 1..%d",Node,NumNodes);
				exit(1);
				}
			}

		SampleN = sample->count ? sample->ival[0] : 1;
		if(SampleN < 1 || SampleN > 10000000)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Sampling '-S%d' specified outside of range 1..10000000",SampleN);
			exit(1);
			}

		if(crick->count > 0)
			bWatsonOnly = false;
		}

	szOutFile[0] = '\0';
	if(SampleN > 1 && outfile->count)
		gDiagnostics.DiagOut(eDLWarn,gszProcName,"Warning: When sampling no output to file is supported");
	else
		{
		if(outfile->count)
			{
			strncpy(szOutFile,outfile->filename[0],_MAX_PATH);
			szOutFile[_MAX_PATH-1] = '\0';
			}
		}
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing parameters:");

	const char *pszDescr;
	switch(PMode) {
		case ePMrestrict:
			pszDescr = "Single node restricted Hamming edit distance processing";
			break;
		case ePMnode:
			pszDescr = "Single node exhaustive Hamming edit distance processing";
			break;
		case ePMdist:
			pszDescr = "Multiple nodes exhaustive Hamming edit distance processing";
			break;
		case ePMmerge:
			pszDescr = "Merge two existing Hamming edit distance files";
			break;
		case ePMtrans:
			pszDescr = "Transform Hamming csv file into quick load binary file";
			break;
		case ePMtransCSV:
			pszDescr = "Transform Hamming quick load binary file into csv file";
			break;
		}
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Processing mode is : '%s'",pszDescr);

	switch(PMode) {
		case ePMnode:
			if(SweepStart == 1 && SweepEnd == 0)
				gDiagnostics.DiagOutMsgOnly(eDLInfo,"Process all sweep instances");
			else
				{
				gDiagnostics.DiagOutMsgOnly(eDLInfo,"Process starting from this sweep instance : %d",SweepStart);
				if(SweepEnd == 0)
					gDiagnostics.DiagOutMsgOnly(eDLInfo,"Process all remaining sweeps");
				else
					gDiagnostics.DiagOutMsgOnly(eDLInfo,"Process until this sweep instance: %d",SweepEnd);
				}
			break;
		case ePMdist:
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Distribute processing over this many nodes: %d",NumNodes);
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Process for node instance: %d",Node);
			SweepStart = NumNodes;
			SweepEnd = Node;
			break;
		}

	switch(PMode) {
		case ePMmerge:
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Merge from this Hamming distance file: '%s'",szInFile);
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Merge into this Hamming distance file: '%s'",szOutFile);
			break;
		case ePMtrans:
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Transform from this Hamming distance CSV file: '%s'",szInFile);
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Transform as quick load binary into this file: '%s'",szOutFile);
			break;
		case ePMtransCSV:
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Transform from this quick load binary file: '%s'",szInFile);
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Transform as CSV into this file: '%s'",szOutFile);
			break;
		case ePMnode:
		case ePMdist:
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Process %s",bWatsonOnly ? "Watson only strand" : "both Watson and Crick strands - Caution: very slow -");
			if(SampleN > 1)
				gDiagnostics.DiagOutMsgOnly(eDLInfo,"Only sample (process) every %d sweep",SampleN);
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Process k-mer subsequences of this length: %d",SeqLen);
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Input bioseq genome assembly file: '%s'",szInFile);
			if(szOutFile[0] != '\0')
				gDiagnostics.DiagOutMsgOnly(eDLInfo,"Output Hamming distance file: '%s'",szOutFile);
			else
				gDiagnostics.DiagOutMsgOnly(eDLInfo,"Summary only, no output Hamming distance file will be generated");
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"number of threads : %d",NumThreads);
			break;

		case ePMrestrict:
			switch(Sensitivity) {
				case eSensDefault:
					pszDescr = "default";
					break;
				case eSensMore:
					pszDescr = "more sensitive - slower";
					break;
				case eSensUltra:
					pszDescr = "ultra sensitive - much slower";
					break;
				case eSensLess:
					pszDescr = "less sensitive - quicker";
					break;
					}
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Sensitivity is : '%s'",pszDescr);
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Process %s",bWatsonOnly ? "Watson only strand" : "both Watson and Crick strands");
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Process k-mer subsequences of this length: %d",SeqLen);
			if(SampleN > 1)
				gDiagnostics.DiagOutMsgOnly(eDLInfo,"Only sample (process) every %d K-mer",SampleN);
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Input genome assembly suffix array file: '%s'",szInFile);
			if(szOutFile[0] != '\0')
				{
				switch(ResFormat) {
					case cRHFcsv:
						pszDescr = "CSV with Hamming loci ranges (chrom, startloci, length, Hamming)";
						break;
					case cRHFBed:
						pszDescr = "UCSC BED";
						break;

					case cRHFWiggle:
						pszDescr = "UCSC Wiggle";
						break;

					}
				gDiagnostics.DiagOutMsgOnly(eDLInfo,"Output Hamming distance file format: %s",pszDescr);
				gDiagnostics.DiagOutMsgOnly(eDLInfo,"Output Hamming distance file: '%s'",szOutFile);
				}
			else
				gDiagnostics.DiagOutMsgOnly(eDLInfo,"Summary only, no output Hamming distance file will be generated");
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Restrict to at most Hamming: %d",RHamm);
			if(szInSeqFile[0] != '\0')
				gDiagnostics.DiagOutMsgOnly(eDLInfo,"Restricted hamming processing source k-mer sequence file: '%s'",szInSeqFile);
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"number of threads : %d",NumThreads);
			break;
		}


#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	gStopWatch.Start();
	Rslt = Process(PMode,bWatsonOnly,Sensitivity,ResFormat,RHamm,(UINT32)SweepStart,(UINT32)SweepEnd,SeqLen,SampleN,NumThreads,szInFile,szInSeqFile,szOutFile);
	gStopWatch.Stop();
	Rslt = Rslt >=0 ? 0 : 1;
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Exit code: %d Total processing time: %s",Rslt,gStopWatch.Read());
	exit(Rslt);
	}
else
	{
	printf("\n%s theK-mer Hamming distance generator, Version %s\n",gszProcName,cpszProgVer);
	arg_print_errors(stdout,end,gszProcName);
	arg_print_syntax(stdout,argtable,"\nUse '-h' to view option and parameter usage\n");
	exit(1);
	}
return 0;
}


typedef struct TAG_sGChrom {
	int ChromSeqID;			// m_pChromSeqs[ChromSeqID-1] of this tsChromSeq instance
	int ChromID;
	char szChromName[cMaxDatasetSpeciesChrom];
	UINT32 Len;				// actual length of this chromosome
	UINT32 SeqOfs;			// offset in m_pGenomeSeq[] at which this chromosome sequence starts
	etSeqBase *pSeq;		// pts to start of this chromosome sequence in m_pGenomeSeq
	UINT32 NumSubSeqs;		// number of subsequences of length m_SeqLen on this chromosome
} tsGChrom;

#ifdef _WIN32
HANDLE m_hMtxThreadParams;
unsigned __stdcall ThreadedGHamDist(void * pThreadPars);
#else
pthread_mutex_t m_hMtxThreadParams;
void *ThreadedGHamDist(void * pThreadPars);
#endif

etPMode m_PMode;				// processing mode
UINT32 m_SubSeqLen;				// subsequence lengths over which to calc Hamming edit distances
int m_NumProcThreads;			// actual number of processing threads
int m_hOutFile;					// output results file handle
char *m_pszOutFile;				// output file name

CBioSeqFile *m_pBioSeqFile;		// genome assembly

char m_szSpecies[cMaxDatasetSpeciesChrom+1];		// species (title) as retrieved from bioseqfile
int m_NumGChroms;				// number of chromosomes loaded
tsGChrom *m_pGChroms;			// pts to gchrom array
UINT8 *m_pGenomeSeq;			// allocated to hold concatenated (starts with eBaseEOG marker) chromosome sequences terminated by final eBaseEOG
INT64 m_AllocGenomeSeq;			// allocation size for m_pGenomeSeq
UINT32 m_GenomeLen;				// total genome length including separator markers but excluding inital and final eBaseEOG
UINT16 *m_pHamDist;				// allocated to hold Hamming edit distances (expected to be m_NumProcThreads * (m_GenomeLen-2)), start/final eBaseEOG have no hammings
INT64 m_AllocHamDist;			// allocation size for m_pHamDist
UINT32 m_NumSubSeqs;			// total number of actual Hamming subsequences in genome

double m_PercentComplete;		// proportion of all hammings generated

tsHamHdr *m_pHamHdr;			// header for binary format hamming edit distances
tsHamChrom *m_pCurHamChrom;		// pts to current chromosome specific binary hammings

tsHamHdr *m_pPregenHamHdr;		// header for loading pregenerated binary format hamming edit distances

tsThreadParams *m_pThreadParams;	// holds initialised phase parameter sets for each thread or processing core

int m_NumThreads;					// number of threads
int m_PerThreadAllocdIdentNodes;    // each thread can use this many tsIdentNodes
int m_TotAllocdIdentNodes;			// total number of tsIdentNodes allocated
tsIdentNode *m_pAllocsIdentNodes;	// memory allocated to hold tsIdentNodes required by all threads
etSensitivity m_Sensitivity;		// restricted hamming processing sensitivity
int m_RHamm;						// if > 0 then restricted hammings limit
int m_SampleN;						// sample (process) every N sweep instance (or if restricted Hammings then every Nth K-mer) 
int m_KMerLen;						// Hammings for these K-mer length sequences
bool m_bWatsonOnly;					// true if watson strand only processing
UINT8 *m_pRHammings;				// to hold restricted hammings
size_t m_TotAllocHammings;			// memory allocation size for holding restricted hammings

CSfxArrayV3 *m_pSfxArray; // suffix array holds genome of interest
char m_szTargSpecies[cMaxDatasetSpeciesChrom+1]; // suffix array was generated over this targeted species

void
DumpHamPair(int Seq1Idx,int Seq2Idx,int PrvHam1,int PrvHam2,int HamCnt,etSeqBase *pSeq1,etSeqBase *pSeq2)
{
char szSeq1[1000];
char szSeq2[1000];
CSeqTrans::MapSeq2Ascii(pSeq1,m_SubSeqLen,szSeq1);
CSeqTrans::MapSeq2Ascii(pSeq2,m_SubSeqLen,szSeq2);
gDiagnostics.DiagOutMsgOnly(eDLInfo,"%1.2d = %1.2d %1.3d %s\n\t     %1.2d %1.3d %s",HamCnt,PrvHam1,Seq1Idx,szSeq1,PrvHam2,Seq2Idx,szSeq2);
}

void
DumpSeq(int SeqIdx)
{
char szSeq1[1000];
CSeqTrans::MapSeq2Ascii(&m_pGenomeSeq[SeqIdx],m_SubSeqLen,szSeq1);
gDiagnostics.DiagOutMsgOnly(eDLInfo,"%1.8d = %s",SeqIdx,szSeq1);
}


void
Reset(bool bSync)	// true if output files to be commited/fsync'd before closing
{
if(m_hOutFile != -1)
	{
	if(bSync)
#ifdef _WIN32
		_commit(m_hOutFile);
#else
		fsync(m_hOutFile);
#endif
	close(m_hOutFile);
	m_hOutFile = -1;
	}
if(m_pBioSeqFile != NULL)
	{
	delete m_pBioSeqFile;
	m_pBioSeqFile = NULL;
	}
if(m_pGChroms != NULL)
	{
	delete m_pGChroms;
	m_pGChroms = NULL;
	}
if(m_pGenomeSeq != NULL)
	{
#ifdef _WIN32
	free(m_pGenomeSeq);
#else
	if(m_pGenomeSeq != MAP_FAILED)
		munmap(m_pGenomeSeq,m_AllocGenomeSeq);
#endif
	m_pGenomeSeq = NULL;
	}
if(m_pHamDist != NULL)
	{
#ifdef _WIN32
	free(m_pHamDist);
#else
	if(m_pHamDist != MAP_FAILED)
		munmap(m_pHamDist,m_AllocHamDist);
#endif
	m_pHamDist = NULL;
	}
if(m_pThreadParams != NULL)
	{
	delete m_pThreadParams;
	m_pThreadParams = NULL;
	}
if(m_pHamHdr != NULL)
	{
	free(m_pHamHdr);					// was allocated with malloc/realloc
	m_pHamHdr = NULL;
	}
if(m_pPregenHamHdr != NULL)
	{
	free(m_pPregenHamHdr);
	m_pPregenHamHdr = NULL;
	}

if(m_pSfxArray != NULL)
	{
	delete m_pSfxArray;
	m_pSfxArray = NULL;
	}

if(m_pRHammings != NULL)
	{
#ifdef _WIN32
	free(m_pRHammings);
#else
	if(m_pRHammings != MAP_FAILED)
		munmap(m_pRHammings,m_TotAllocHammings);
#endif
	m_pRHammings = NULL;
	}

if(m_pAllocsIdentNodes != NULL)
	{
	delete m_pAllocsIdentNodes;
	m_pAllocsIdentNodes = NULL;
	}

m_TotAllocHammings = 0;
m_AllocGenomeSeq = 0;
m_AllocHamDist = 0;
m_pCurHamChrom = NULL;
m_SubSeqLen=0;
m_NumProcThreads=0;
m_szSpecies[0] = '\0';
m_szTargSpecies[0] = '\0';
m_GenomeLen = 0;
m_NumSubSeqs = 0;
m_NumGChroms= 0;
m_NumThreads = 1;
m_TotAllocdIdentNodes = 0;
m_PerThreadAllocdIdentNodes = 0;
m_PercentComplete = 0.0;
}

void
Init(void)
{
m_pBioSeqFile = NULL;
m_pGChroms = NULL;
m_pCurHamChrom = NULL;
m_pHamHdr = NULL;
m_pPregenHamHdr = NULL;
m_pSfxArray = NULL;
m_pRHammings = NULL;
m_hOutFile = -1;
m_pGenomeSeq = NULL;
m_pHamDist = NULL;
m_pThreadParams = NULL;
m_pAllocsIdentNodes = NULL;
m_TotAllocdIdentNodes = 0;
m_PerThreadAllocdIdentNodes = 0;
Reset(false);
}

#ifdef _WIN32
unsigned __stdcall ThreadedGHamDist(void * pThreadPars)
#else
void *ThreadedGHamDist(void * pThreadPars)
#endif
{
tsWorkerThreads *pPars = (tsWorkerThreads *)pThreadPars; // makes it easier not having to deal with casts!

// iterate over thread processing parameter sets until processing completed for this thread
tsThreadParams *pParams;
UINT32 Seq;
UINT32 Idx;
double PercentComplete;
while((pParams = ThreadedIterParams())!=NULL)
	{
	for(Seq = pParams->SSofs, Idx = 0; Idx < pParams->NumSeqs; Idx++, Seq+=pParams->SeqDelta)
		{
		if((Seq + m_SubSeqLen) >= m_GenomeLen)
			break;
		if(pParams->SampleN > 1)	// if sampling....
			{
			if(Idx > 0 && (Idx % pParams->SampleN))
				continue;
			}
		if(!(Idx % 10))
			{
			PercentComplete = (double)((UINT64)Idx*100)/pParams->NumSeqs;
#ifdef _WIN32
			WaitForSingleObject(m_hMtxThreadParams,INFINITE);
#else
			pthread_mutex_lock(&m_hMtxThreadParams);
#endif
			if(PercentComplete > m_PercentComplete)
				PercentComplete = PercentComplete;
#ifdef _WIN32
			ReleaseMutex(m_hMtxThreadParams);
#else
			pthread_mutex_unlock(&m_hMtxThreadParams);
#endif
			}

		// no hammings to be generated for inital/final eBaseEOG markers
		if((Seq + m_SubSeqLen) < m_GenomeLen-1)
			GHamDistWatson(&m_pHamDist[pParams->HamDistOfs], m_SubSeqLen,Seq,&m_pGenomeSeq[1],m_GenomeLen-2);

		if(!pParams->bWatsonOnly && (Seq + m_SubSeqLen) < m_GenomeLen)
			GHamDistCrick(&m_pHamDist[pParams->HamDistOfs], m_SubSeqLen,Seq-1,&m_pGenomeSeq[1],m_GenomeLen-2);
		}
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Thread %d processing subsequence %d completed",pParams->ThreadID,Seq-1);
	}

#ifdef _WIN32
_endthreadex(0);
return(eBSFSuccess);
#else
pthread_exit(NULL);
#endif
}

// LoadHammings
// Loads pregenerated binary format Hammings file into memory
teBSFrsltCodes
LoadPregenHammings(char *pszHammings)
{
int hHamFile;
tsHamHdr HamHdr;

if(m_pPregenHamHdr != NULL)
	{
	free(m_pPregenHamHdr);
	m_pPregenHamHdr = NULL;
	}

#ifdef _WIN32
hHamFile = open(pszHammings, O_READSEQ ); // file access will be sequential..
#else
hHamFile = open64(pszHammings, O_READSEQ ); // file access will be sequential..
#endif
if(hHamFile == -1)					// check if file open succeeded
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open %s - %s",pszHammings,strerror(errno));
	Reset(false);
	return(eBSFerrOpnFile);
	}

if(read(hHamFile,&HamHdr,sizeof(tsHamHdr))!=sizeof(tsHamHdr))
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to read %s - %s",pszHammings,strerror(errno));
	close(hHamFile);
	Reset(false);
	return(eBSFerrParse);
	}
if(HamHdr.Magic[0] != 'b' || HamHdr.Magic[1] != 'h' || HamHdr.Magic[2] != 'a' || HamHdr.Magic[3] != 'm')
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Not a binary format Hamming distance file - %s",pszHammings);
	close(hHamFile);
	Reset(false);
	return(eBSFerrFileType);
	}
if(HamHdr.NumChroms < 1)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Hamming distance file contains no edit distances - %s",pszHammings);
	close(hHamFile);
	Reset(false);
	return(eBSFerrNoEntries);
	}

if((m_pPregenHamHdr = (tsHamHdr *)new UINT8 [HamHdr.Len])==NULL)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Unable to allocate memory (%d bytes) for holding Hamming distances loaded from - %s",HamHdr.Len,pszHammings);
	close(hHamFile);
	Reset(false);
	return(eBSFerrMem);
	}
memcpy(m_pPregenHamHdr,&HamHdr,sizeof(tsHamHdr));
if(read(hHamFile,(UINT8 *)m_pPregenHamHdr+sizeof(tsHamHdr),HamHdr.Len-sizeof(tsHamHdr))!=HamHdr.Len-sizeof(tsHamHdr))
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Unable to read all Hamming edit distances from - %s",pszHammings);
	close(hHamFile);
	Reset(false);
	return(eBSFerrFileAccess);
	}

close(hHamFile);
return(eBSFSuccess);
}

// TransHammingsCSV
// Write out existing binary format Hammings file as CSV file for subsequent post processing
int
TransHammingsCSV(char *pszHamFile,char *pszOutFile)
{
int Rslt;
int ChromIdx;
UINT32 Loci;
tsHamChrom *pCurChrom;
char szBuff[32000];
int BuffOfs;

if((Rslt=LoadPregenHammings(pszHamFile))!=eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Unable to dump Hammings as csv from - %s",pszHamFile);
	Reset(false);
	return(eBSFerrFileAccess);
	}

#ifdef _WIN32
m_hOutFile = open(pszOutFile,O_CREATETRUNC );
#else
if((m_hOutFile = open(pszOutFile,O_RDWR | O_CREAT,S_IREAD | S_IWRITE))!=-1)
    if(ftruncate(m_hOutFile,0)!=0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate %s - %s",pszOutFile,strerror(errno));
			Reset(false);
			return(eBSFerrCreateFile);
			}
#endif

if(m_hOutFile < 0)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate output file '%s'",pszOutFile);
	Reset(false);
	return(eBSFerrCreateFile);
	}

BuffOfs = sprintf(szBuff,"\"Chrom\",\"Loci\",\"Hamming\"\n");
for(ChromIdx = 0; ChromIdx < m_pPregenHamHdr->NumChroms; ChromIdx++)
	{
	pCurChrom = (tsHamChrom *)((UINT8 *)m_pPregenHamHdr + m_pPregenHamHdr->ChromOfs[ChromIdx]);
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Transforming chrom: '%s' with %d Hammings",pCurChrom->szChrom,pCurChrom->NumEls);
	for(Loci = 0; Loci < pCurChrom->NumEls; Loci++)
		{
		BuffOfs += sprintf(&szBuff[BuffOfs],"\"%s\",%d,%d\n",pCurChrom->szChrom,Loci,pCurChrom->Dists[Loci]);
		if(BuffOfs + 100 > sizeof(szBuff))
			{
			CUtility::SafeWrite(m_hOutFile,szBuff,BuffOfs);
			BuffOfs = 0;
			}
		}
	}
if(BuffOfs)
	CUtility::SafeWrite(m_hOutFile,szBuff,BuffOfs);
Reset(true);
return(eBSFSuccess);
}

int		// returned minimum Hamming distance
MinHammingDist(int MinHamming,	// initial minimum Hamming distance
			   int ReadLen,		// read length
			   int ChromID,     // from which chromosome was this read was derived
			   int ReadLoci,	// read start loci
			   etSeqBase *pRead) // read sequence
{
int ChromIdx;
int SeqIdx;
int ChromSeqIdx;
int CurHamming;
int EndIdx;
tsGChrom *pChrom;
etSeqBase *pChromBase;
etSeqBase *pChromSeq;
etSeqBase *pReadBase;

if(MinHamming <= 0)
	return(0);

pChrom = &m_pGChroms[0];
for(ChromIdx = 0; ChromIdx < m_NumGChroms; ChromIdx++, pChrom++)
	{
	pChromSeq = pChrom->pSeq;
	EndIdx = pChrom->Len - ReadLen;
	for(ChromSeqIdx=0;ChromSeqIdx <= EndIdx; ChromSeqIdx++,pChromSeq++)
		{
		if(ChromSeqIdx == ReadLoci && ChromID == pChrom->ChromID)
			continue;

		pChromBase = pChromSeq;
		pReadBase = pRead;
		CurHamming = 0;
		for(SeqIdx=0;SeqIdx < ReadLen; SeqIdx++,pChromBase++,pReadBase++)
			if(*pChromBase != *pReadBase)
				{
				CurHamming++;
				if(CurHamming > MinHamming)
					break;
				}
		if(CurHamming < MinHamming)
			{
			MinHamming = CurHamming;
			if(MinHamming == 0)
				return(MinHamming);
			}
		}
	}
return(MinHamming);
}


typedef struct TAG_sMergeHammings {
	int ChromID;
	int Loci;
	int Dist;		// current minimum edit distances
} sMergeHammings;

int
MergeHammings(char *pszMergeFromFile,char *pszMergeIntoFile)
{
int RsltFrom;
int RsltInto;
bool bCopy;
char szBuff[10000];
int BuffIdx;

CCSVFile *pMergeFrom;
CCSVFile *pMergeInto;
char szTmpFile[_MAX_PATH];
int hOutFile;

// if pszMergeIntoFile is empty or doesn't exist then simply make a copy of pszMergeFromFile
// a) load record from pszMergeFromFile
// b) load record from pszMergeToFile
// c) ensure chrom and loci are same
// d) writeout record with minimum of Hamming distance

if((pMergeFrom = new CCSVFile())==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CCSVfile");
	return(eBSFerrObj);
	}
if((RsltFrom=pMergeFrom->Open(pszMergeFromFile))!=eBSFSuccess)
	{
	while(pMergeFrom->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,pMergeFrom->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open file: %s",pszMergeFromFile);
	delete pMergeFrom;
	return(RsltFrom);
	}

if((pMergeInto = new CCSVFile())==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CCSVfile");
	delete pMergeFrom;
	return(eBSFerrObj);
	}

char *pszFromChrom;
int FromLoci;
int FromHammingDist;
char *pszToChrom;
int ToLoci;
int ToHammingDist;
int NumFields;
int LineNum;
int NumErrs;
int Idx;
int NumSubSeqs;

hOutFile = -1;

// now try opening the merge into Hamming file
// if the open fails then assume file does not exist and the merged into file will be a
// simple copy of the merged from file
if((RsltInto=pMergeInto->Open(pszMergeIntoFile))!=eBSFSuccess)
	{
	strcpy(szTmpFile,pszMergeIntoFile);
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Unable to open '%s' for merging, will create and copy from '%s'",pszMergeIntoFile,pszMergeFromFile);
	RsltInto = eBSFSuccess;
	bCopy = true;
	}
else
	{
	strcpy(szTmpFile,pszMergeIntoFile);
	strcat(szTmpFile,".tmp");
	bCopy = false;
	}

// output results
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Creating%sfile '%s' to hold merged Hammings",bCopy ? " " : " temp ", szTmpFile);

#ifdef _WIN32
hOutFile = open(szTmpFile,O_CREATETRUNC );
#else
if((hOutFile = open(szTmpFile,O_RDWR | O_CREAT,S_IREAD | S_IWRITE))!=-1)
    if(ftruncate(hOutFile,0)!=0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate %s - %s",szTmpFile,strerror(errno));
			delete pMergeInto;
			delete pMergeFrom;
			return(eBSFerrCreateFile);
			}
#endif

if(hOutFile < 0)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate output file '%s'",szTmpFile);
	delete pMergeInto;
	delete pMergeFrom;
	Reset(false);
	return(eBSFerrCreateFile);
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Starting %s from '%s' into '%s' ...",bCopy ? "copy" : "merge", pszMergeFromFile,szTmpFile);

int CntHist[256];					// used to generate distance histogram
memset(CntHist,0,sizeof(CntHist));

LineNum = 0;
NumErrs = 0;
NumSubSeqs = 0;
BuffIdx = sprintf(szBuff,"\"Chrom\",\"Loci\",\"Hamming\"");
while(1)	// onto next line containing fields
	{
	LineNum+=1;
	if(!(LineNum % 10000000))
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing line %d",LineNum);
	if((RsltFrom=pMergeFrom->NextLine()) <= 0)
		break;
	if(!bCopy && (RsltInto=pMergeInto->NextLine()) <= 0)
		break;

	NumFields = pMergeFrom->GetCurFields();
	if(NumFields < 3)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expected at least %d fields at line %d in '%s', only parsed '%d'",3,LineNum,pszMergeFromFile,NumFields);
		RsltFrom = -1;
		break;
		}
	if(!bCopy)
		{
		NumFields = pMergeInto->GetCurFields();
		if(NumFields < 3)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expected at least %d fields at line %d in '%s', only parsed '%d'",3,LineNum,pszMergeIntoFile,NumFields);
			RsltFrom = -1;
			break;
			}
		}

	// simply slough any records which are not in the format of "chrom",loci,distance
	// assume these are descriptor records
	if(!pMergeFrom->GetQuoted(1) || pMergeFrom->GetQuoted(2) || pMergeFrom->GetQuoted(3))
		continue;
	if(!bCopy && (!pMergeInto->GetQuoted(1) || pMergeInto->GetQuoted(2) || pMergeInto->GetQuoted(3)))
		continue;

	pMergeFrom->GetText(1,&pszFromChrom);
	pMergeFrom->GetInt(2,&FromLoci);
	pMergeFrom->GetInt(3,&FromHammingDist);
	if(!bCopy)
		{
		pMergeInto->GetText(1,&pszToChrom);
		pMergeInto->GetInt(2,&ToLoci);
		pMergeInto->GetInt(3,&ToHammingDist);
		// check that the chrom.loci are identical, which they will be if generated by this application
		// because of a bug (now fixed) in the first version of this application then a non-matching last record is simply sloughed
		if(stricmp(pszFromChrom,pszToChrom) || FromLoci != ToLoci)
			{
			if(NumErrs++ > 1)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expected chromosome and loci to be identical at line %d between merge files %s and %s",LineNum,pszMergeIntoFile,pszMergeFromFile);
				RsltFrom = -1;
				break;
				}
			continue;
			}
		}
	else
		ToHammingDist = FromHammingDist;

	NumSubSeqs += 1;
	FromHammingDist = min(FromHammingDist,ToHammingDist);
	CntHist[FromHammingDist] += 1;

	// here need to write out the merged record into the temp file
	BuffIdx += sprintf(&szBuff[BuffIdx],"\n\"%s\",%d,%d",pszFromChrom,FromLoci,FromHammingDist);
	if((BuffIdx + 100) > sizeof(szBuff))
		{
		CUtility::SafeWrite(hOutFile,szBuff,BuffIdx);
		BuffIdx = 0;
		}
	}

if(hOutFile != -1 && (RsltFrom >= 0 && RsltInto >= 0))
	{
	if(BuffIdx)
		CUtility::SafeWrite(hOutFile,szBuff,BuffIdx);
	}
if(hOutFile != -1)
	close(hOutFile);
hOutFile = -1;
pMergeInto->Close();
pMergeFrom->Close();
delete pMergeInto;
delete pMergeFrom;


if(RsltFrom >= 0 && RsltInto >= 0)
	{
	// log the basic count distribution.....
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Distribution:\nEditDist,Freq,Proportion");
	for(Idx = 0; Idx < 50; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"%d,%d,%1.3f",Idx,CntHist[Idx],(CntHist[Idx]*100.0f)/NumSubSeqs);

	if(!bCopy)
		{
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Renaming '%s' to '%s'",szTmpFile,pszMergeIntoFile);
		if(remove(pszMergeIntoFile))
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to delete '%s' prior to renaming",pszMergeIntoFile);
			return(-1);
			}
		if(rename(szTmpFile,pszMergeIntoFile))
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to rename '%s' to '%s'",szTmpFile,pszMergeIntoFile);
			return(-1);
			}
		}
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Merge of %d subsequences completed",NumSubSeqs);
	}

return(min(RsltFrom,RsltInto));
}

teBSFrsltCodes
LoadCSVHammings(char *pszHammings)
{
int Rslt;
int NumFields;
int LineNum;
char *pszChrom;
int Loci;
int Dist;
char szCurChrom[cMaxDatasetSpeciesChrom+1];
int CurChromID;
int AllocLen;

tsHamChrom *pCurHamChrom;
tsHamHdr *pTmpHdr;
CCSVFile *pHammings;

if(m_pHamHdr != NULL)
	free(m_pHamHdr);
m_pCurHamChrom = NULL;
if((m_pHamHdr = (tsHamHdr *)malloc(cAllocHamDist))==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate memory (%d bytes) for Hamming distances",cAllocHamDist);
	return(eBSFerrMem);
	}
AllocLen = cAllocHamDist;
memset(m_pHamHdr,0,sizeof(tsHamHdr));

m_pHamHdr->Magic[0] = 'b';
m_pHamHdr->Magic[1] = 'h';
m_pHamHdr->Magic[2] = 'a';
m_pHamHdr->Magic[3] = 'm';
m_pHamHdr->Version = 1;
m_pHamHdr->Len = sizeof(tsHamHdr);

if((pHammings = new CCSVFile())==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CCSVfile");
	return(eBSFerrObj);
	}
if((Rslt=pHammings->Open(pszHammings))!=eBSFSuccess)
	{
	while(pHammings->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,pHammings->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open file: %s",pszHammings);
	delete pHammings;
	return((teBSFrsltCodes)Rslt);
	}

LineNum = 0;
szCurChrom[0] = '\0';
CurChromID = -1;
while((Rslt=pHammings->NextLine()) > 0)
	{
	LineNum += 1;
	NumFields = pHammings->GetCurFields();
	if(NumFields < 3)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expected at least %d fields at line %d in '%s', only parsed '%d'",3,LineNum,pszHammings,NumFields);
		Rslt = -1;
		break;
		}
	// simply slough any records which are not in the format of "chrom",loci,distance
	// assume these are descriptor records
	if(!pHammings->GetQuoted(1) || pHammings->GetQuoted(2) || pHammings->GetQuoted(3))
		continue;

	pHammings->GetText(1,&pszChrom);
	pHammings->GetInt(2,&Loci);
	pHammings->GetInt(3,&Dist);

	// check if needing to start a new chrom, and if more memory is needed to hold the new chrom...
	if(szCurChrom[0] == '\0' || stricmp(pszChrom,szCurChrom))
		{
		if(m_pHamHdr->NumChroms > 0)
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"chrom %s has %d Hammings",szCurChrom,pCurHamChrom->NumEls);
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Loading Hammings for chrom: %s",pszChrom);
		if((AllocLen - m_pHamHdr->Len) < (sizeof(tsHamChrom) + (cAllocHamDist/16)))	// allow some hammings as well as the tsHamChrom
			{
			if((pTmpHdr = (tsHamHdr *)realloc(m_pHamHdr,AllocLen + cAllocHamDist))==NULL)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to reallocate memory (%d bytes) for Hamming distances",AllocLen + cAllocHamDist);
				return(eBSFerrMem);
				}
			m_pHamHdr = pTmpHdr;
			AllocLen += cAllocHamDist;
			}
		m_pHamHdr->ChromOfs[m_pHamHdr->NumChroms] = m_pHamHdr->Len;
		m_pHamHdr->Len += sizeof(tsHamChrom) - 1;
		pCurHamChrom = (tsHamChrom *)((UINT8 *)m_pHamHdr + m_pHamHdr->ChromOfs[m_pHamHdr->NumChroms++]);
		pCurHamChrom->NumEls = 0;
		pCurHamChrom->ChromID = m_pHamHdr->NumChroms;
		strcpy((char *)pCurHamChrom->szChrom,pszChrom);
		strcpy(szCurChrom,pszChrom);
		CurChromID = m_pHamHdr->NumChroms;
		}

	// check if needing to extend
	if(m_pHamHdr->Len == AllocLen)
		{
		if((pTmpHdr = (tsHamHdr *)realloc(m_pHamHdr,AllocLen + cAllocHamDist))==NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to reallocate memory (%d bytes) for Hamming distances",AllocLen + cAllocHamDist);
			return(eBSFerrMem);
			}
		m_pHamHdr = pTmpHdr;
		AllocLen += cAllocHamDist;
		pCurHamChrom = (tsHamChrom *)((UINT8 *)m_pHamHdr + m_pHamHdr->ChromOfs[m_pHamHdr->NumChroms-1]);
		}

	// check that the hamming loci are monotonically ascending
	if(Loci != pCurHamChrom->NumEls)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Hamming loci are not monotonically ascending on chrom %s, check line %d in '%s'",szCurChrom,LineNum,pszHammings);
		return(eBSFerrParse);
		}

	pCurHamChrom->Dists[Loci] = Dist;
	pCurHamChrom->NumEls += 1;
	m_pHamHdr->Len += 1;
	}
if(m_pHamHdr->NumChroms > 0)
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"chrom %s has %d Hammings",szCurChrom,pCurHamChrom->NumEls);
m_pCurHamChrom = NULL;
delete pHammings;

return((teBSFrsltCodes)Rslt);
}

int
TransHammings(char *pszInFile,char *pszOutFile)
{
int Rslt;

// load the CSV Hammings into memory
if((Rslt=LoadCSVHammings(pszInFile)) < eBSFSuccess)
	return(Rslt);

// output as binary
#ifdef _WIN32
m_hOutFile = open(pszOutFile,O_CREATETRUNC );
#else
if((m_hOutFile = open(pszOutFile,O_RDWR | O_CREAT,S_IREAD | S_IWRITE))!=-1)
    if(ftruncate(m_hOutFile,0)!=0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate %s - %s",pszOutFile,strerror(errno));
			Reset(false);
			return(eBSFerrCreateFile);
			}
#endif

if(m_hOutFile < 0)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate output file '%s'",pszOutFile);
	Reset(false);
	return(eBSFerrCreateFile);
	}

CUtility::SafeWrite(m_hOutFile,m_pHamHdr,m_pHamHdr->Len);

// that's it!
#ifdef _WIN32
	_commit(m_hOutFile);
#else
	fsync(m_hOutFile);
#endif
close(m_hOutFile);
m_hOutFile = -1;
return(eBSFSuccess);
}

// EstL
// Estimate Lx from genome length, number of nodes and current node instance
UINT32
EstL(UINT32 GenomeLen,int NodeInst,int NumNodes)
{
if(NodeInst < 1)
	return(0);
if(NodeInst == NumNodes)
	return(GenomeLen);
double TotWork = ((double)GenomeLen * GenomeLen) + ((double)GenomeLen * GenomeLen)/2;
double TargWork = (TotWork * NodeInst) / NumNodes;
double ActWork;
double L2;
double L2Est;

// make an initial estimate and then iteratively refine this estimate until the estimated work load is the same as the targeted load
L2 = sqrt((2.0f * NodeInst * TotWork)/NumNodes);
if(L2 > GenomeLen)
	L2 = GenomeLen;
do {
	L2Est = L2;
	ActWork = ((UINT64)GenomeLen * L2Est) + (L2Est * L2Est)/2;
	double WrkRatio = TargWork/ActWork;
	L2 = ((L2Est * TargWork)/ActWork);
	}
while(abs((long)(((INT64)L2Est - (INT64)L2))) > 1);
return((UINT32)L2Est);
}

// following are the structures used for restricted hamming processing
#pragma pack(1)
typedef struct TAG_sHammProcState {
	UINT32 NumEntryIDs;		// number of entry identifiers
	UINT32 CurEntryID;		// current entry for which hammings are being generated
	UINT32 CurEntryLen;		// total length of this entries sequence
	UINT32 CurSeqOfs;		// last sequence was processed starting at this chrom offset
	UINT32 CurSeqLen;		// length of sequence which was last processed
	UINT32 NxtSeqOfs;		// start next hamming generation from this offset
	bool bWatsonOnly;		// true if watson (sense) strand only processing
	int RHamm;				// restricted hammings limit
	UINT32 MinCoreDepth;	// initially explore cores to at least this depth and only continue exploring if depth would be no more than MaxCoreDepth
	UINT32 MaxCoreDepth;	// explore cores to at most this depth
	UINT32 KMerLen;			// Hammings for these KMer length sequences
	UINT32 SampleN;			// sample (process) every Nth K-mer 
	UINT32 MaxHammSeqLen;	// generate hammings over this maximal length sequence
	bool bProcSepKMers;		// true if kmer sequences are separate from suffix array sequences
	etSeqBase *pSeq;		// if separate kmer sequences then pts to start of sequence
	UINT32 EntryHammOfs;	// offset into pHammings of where current entry hammings to be written
	UINT8 *pHammings;		// pts to start of returned hammings for all entries
	UINT32 ApproxNumChroms; // approximate number of chromosomes for which hammings have beein generated by all threads
	INT64 ApproxNumHammings; // approximate total number of hammings generated by all threads
	} tsHammProcState;

// restricted hamming processing threads
typedef struct TAG_RHThread {
#ifdef _WIN32
	HANDLE threadHandle;			// handle as returned by _beginthreadex()
	unsigned int threadID;			// identifier as set by _beginthreadex()
#else
	int threadRslt;					// result as returned by pthread_create ()
	pthread_t threadID;				// identifier as set by pthread_create ()
#endif
	int NumIdentNodes;		// number of ident nodes allocd for use by this thread
	tsIdentNode *pIdentNodes;	// thread to use these nodes
	UINT32 State;			// processing state
	tsHammProcState *pHammState; // current hamming processing state common to all threads
	UINT32 CurHammSeqLen;		// current length of sequence being processed by this thread
	UINT32 EntryID;			// if using suffix array sequences then the entry identifer
	UINT32 SeqOfs;			// if using suffix array sequences then the sequence offset
	UINT8 *pSeq;			// if separate kmer sequences then pts to start of sequence
	UINT8 *pHammings;		// Hammings returned
} tsRHThread;

#pragma pack()


#ifdef _WIN32
unsigned __stdcall RestrictedHammingThread(void * pThreadPars)
#else
void *RestrictedHammingThread(void * pThreadPars)
#endif
{
int Rslt;
bool bCompleted;
UINT8 *pHammings;
tsHammProcState *pHammState;
tsGChrom *pChrom;
tsRHThread *pPars = (tsRHThread *)pThreadPars;
pHammings = pPars->pHammings;
pHammState = pPars->pHammState;
bCompleted = false;
Rslt = 0;
do {
#ifdef _WIN32
	WaitForSingleObject(m_hMtxThreadParams,INFINITE);
#else
	pthread_mutex_lock(&m_hMtxThreadParams);
#endif
	if(pHammState->CurEntryID <= pHammState->NumEntryIDs) // can only process at most NumEntryIDs chrom sequences
		{
		if(pHammState->CurEntryID == 0 ||	// check if starting or current chromosomes hammings completed...
			((pHammState->NxtSeqOfs + pHammState->KMerLen) > pHammState->CurEntryLen))
			{
			pHammState->CurEntryID += 1;
			if(pHammState->CurEntryID <= pHammState->NumEntryIDs)	// check if still processing within entries
				{
				pHammState->CurSeqOfs = 0;
				if(pHammState->CurEntryID == 1)
					{
					pHammState->EntryHammOfs = 0;
					pHammState->ApproxNumChroms = 1;
					pHammState->ApproxNumHammings = 0;
					}
				else
					{
					pHammState->EntryHammOfs += pHammState->CurEntryLen;
					pHammState->ApproxNumChroms += 1;
					}

				if(pHammState->bProcSepKMers)	// if k-mers in separate sequence....
					{
					pChrom = &m_pGChroms[pHammState->CurEntryID-1];
					pHammState->CurEntryLen = pChrom->Len;
					pHammState->pSeq = pChrom->pSeq;
					}
				else							// k-mers in suffix array sequence....
					{
					pHammState->CurEntryLen = m_pSfxArray->GetSeqLen(pHammState->CurEntryID);
					pHammState->pSeq = NULL;
					}
				pHammState->CurSeqLen = min(pHammState->MaxHammSeqLen,pHammState->CurEntryLen);
				pHammState->NxtSeqOfs = 1 + pHammState->CurSeqLen - pHammState->KMerLen;
				}
			}
		else		// still processing the current chromosome
			{
			pHammState->CurSeqOfs = pHammState->NxtSeqOfs;
			pHammState->CurSeqLen = min(pHammState->MaxHammSeqLen,pHammState->CurEntryLen - pHammState->CurSeqOfs);
			pHammState->NxtSeqOfs += (1 + pHammState->CurSeqLen - pHammState->KMerLen);
			}
		}

	if(pHammState->CurEntryID <= pHammState->NumEntryIDs)
		{
		if(pHammState->CurEntryLen >= pHammState->KMerLen)
			{
			pPars->EntryID = pHammState->CurEntryID;
			pPars->CurHammSeqLen = pHammState->CurSeqLen;
			if(pHammState->bProcSepKMers)
				pPars->pSeq = &pHammState->pSeq[pHammState->CurSeqOfs];
			pPars->SeqOfs = pHammState->CurSeqOfs;
			pPars->pHammings = &pHammState->pHammings[pHammState->EntryHammOfs + pHammState->CurSeqOfs];
			}
		}
	else
		bCompleted = true;
	if(!bCompleted && pHammState->CurEntryLen >= pHammState->KMerLen)
		pHammState->ApproxNumHammings += (INT64)(1 + pPars->CurHammSeqLen - pHammState->KMerLen);

#ifdef _WIN32
	ReleaseMutex(m_hMtxThreadParams);
#else
	pthread_mutex_unlock(&m_hMtxThreadParams);
#endif

	if(bCompleted)
		{
		Rslt = 0;
		break;
		}

	if(pHammState->CurEntryLen >= pHammState->KMerLen)
		{
			
		if(!pHammState->bProcSepKMers)
			Rslt = m_pSfxArray->LocateSfxHammings(pHammState->RHamm,!pHammState->bWatsonOnly,pHammState->KMerLen,pHammState->SampleN,pPars->CurHammSeqLen,pHammState->MinCoreDepth,pHammState->MaxCoreDepth,pPars->EntryID,pPars->SeqOfs,pPars->pHammings,pPars->NumIdentNodes,pPars->pIdentNodes);
		else
			Rslt = m_pSfxArray->LocateHammings(pHammState->RHamm,!pHammState->bWatsonOnly,pHammState->KMerLen,pHammState->SampleN,pPars->CurHammSeqLen,pHammState->MinCoreDepth,pHammState->MaxCoreDepth,pPars->pSeq,0,0,pPars->pHammings,pPars->NumIdentNodes,pPars->pIdentNodes);
		}
	}
while(!bCompleted && Rslt >= 0);

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Thread %d completed",pPars->threadID);
pPars->State = Rslt;

#ifdef _WIN32
_endthreadex(0);
return(Rslt);
#else
pthread_exit(NULL);
#endif
}

// ReportRestrictedHammingsCSV
int
ReportRestrictedHammingsCSV(int KMerLen,						// Hamming for this length k-mer sequences
					bool bProcSepKMers,							// true if processing probe k-mers from sequences not within target suffix array
					char *pszHammingFile)						// report into this output file
{
char *pszBuffer;
int BuffOfs;
UINT32 CurLoci;
int ChromID;
tsGChrom *pChrom;
UINT8 *pRHammings;
int NumEntries;
int ChromLen;
char szChromName[100];
int CurHammingLen;
int CurHamming;
int CurHammingLoci;
UINT32 RHammingsOfs;

	// output results
#ifdef _WIN32
m_hOutFile = open(pszHammingFile,O_CREATETRUNC );
#else
if((m_hOutFile = open(pszHammingFile,O_RDWR | O_CREAT,S_IREAD | S_IWRITE))!=-1)
	if(ftruncate(m_hOutFile,0)!=0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate %s - %s",pszHammingFile,strerror(errno));
		Reset(false);
		return(eBSFerrCreateFile);
		}
#endif

if(m_hOutFile < 0)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate output file '%s'",pszHammingFile);
	Reset(false);
	return(eBSFerrCreateFile);
	}

if((pszBuffer = new char [cRptBuffAllocsize])==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to allocate memory for report buffering");
	Reset(false);
	return(eBSFerrMem);
	}



BuffOfs = sprintf(pszBuffer,"\"Chrom\",\"StartLoci\",\"Len\",\"Hamming\"\n"); // header line
CUtility::SafeWrite(m_hOutFile,pszBuffer,BuffOfs);
BuffOfs = 0;
CurLoci = 0;

pChrom = &m_pGChroms[0];
pRHammings = m_pRHammings;

if(bProcSepKMers)
	NumEntries = m_NumGChroms;
else
	NumEntries = m_pSfxArray->GetNumEntries();

RHammingsOfs	= 0;
for(ChromID = 0; ChromID < NumEntries; ChromID++)
	{
	CurLoci = 0;
	if(bProcSepKMers)
		{
		pChrom = &m_pGChroms[ChromID];
		ChromLen = pChrom->Len;
		strcpy(szChromName,pChrom->szChromName);
		}
	else
		{
		ChromLen = m_pSfxArray->GetSeqLen(ChromID+1);
		m_pSfxArray->GetIdentName(ChromID+1,sizeof(szChromName)-1,szChromName);
		}
	pRHammings = &m_pRHammings[RHammingsOfs];
	RHammingsOfs += ChromLen;
	if(ChromLen < KMerLen)
		continue;

	CurHammingLen = 0;
	CurHamming = *pRHammings;
	CurHammingLoci = 0;
	for(CurLoci = 1; CurLoci <= (UINT32)(ChromLen - KMerLen); CurLoci++,pRHammings++)
		{
		if(*pRHammings == CurHamming)
			{
			CurHammingLen += 1;
			continue;
			}
		BuffOfs += sprintf(&pszBuffer[BuffOfs],"\"%s\",%d,%d,%d\n",szChromName,CurHammingLoci,CurHammingLen,CurHamming);
		CurHammingLen = 1;
		CurHamming = *pRHammings;
		CurHammingLoci = CurLoci;

		if((BuffOfs + 1000) > cRptBuffAllocsize)
			{
			CUtility::SafeWrite(m_hOutFile,pszBuffer,BuffOfs);
			BuffOfs = 0;
			}
		}
	BuffOfs += sprintf(&pszBuffer[BuffOfs],"\"%s\",%d,%d,%d\n",szChromName,CurHammingLoci,CurHammingLen,CurHamming);
	}

if(BuffOfs)
	{
	CUtility::SafeWrite(m_hOutFile,pszBuffer,BuffOfs);
	BuffOfs = 0;
	}
#ifdef _WIN32
_commit(m_hOutFile);
#else
fsync(m_hOutFile);
#endif
close(m_hOutFile);
m_hOutFile = -1;
delete pszBuffer;
return(eBSFSuccess);
}


// ReportRestrictedHammingsBed
int
ReportRestrictedHammingsBed(etSensitivity Sensitivity, // restricted hamming processing sensitivity
					int RHamm,							// restricted hammings limit
					int KMerLen,						// Hamming for this length k-mer sequences
					bool bWatsonOnly,					// true if watson strand only processing
					bool bProcSepKMers,					// true if processing probe k-mers from sequences not within target suffix array
					char *pszHammingFile)				// report into this output file
{
char szBuffer[16000];
int BuffOfs;
UINT32 CurLoci;
int ChromID;
tsGChrom *pChrom;
UINT8 *pRHammings;
int NumEntries;
int ChromLen;
char szChromName[100];
int CurHammingLen;
int CurHamming;
int CurHammingLoci;
UINT32 RHammingsOfs;

	// output results
#ifdef _WIN32
m_hOutFile = open(pszHammingFile,O_CREATETRUNC );
#else
if((m_hOutFile = open(pszHammingFile,O_RDWR | O_CREAT,S_IREAD | S_IWRITE))!=-1)
	if(ftruncate(m_hOutFile,0)!=0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate %s - %s",pszHammingFile,strerror(errno));
		Reset(false);
		return(eBSFerrCreateFile);
		}
#endif

if(m_hOutFile < 0)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate output file '%s'",pszHammingFile);
	Reset(false);
	return(eBSFerrCreateFile);
	}

BuffOfs = sprintf(szBuffer,"track type=bedGraph name=ResHamming%d_%d description=\"Restricted Hammings for K-mer length %d and Hamming limit %d\"\n",KMerLen,RHamm,KMerLen,RHamm); // header line
CUtility::SafeWrite(m_hOutFile,szBuffer,BuffOfs);
BuffOfs = 0;
CurLoci = 0;

pChrom = &m_pGChroms[0];
pRHammings = m_pRHammings;

if(bProcSepKMers)
	NumEntries = m_NumGChroms;
else
	NumEntries = m_pSfxArray->GetNumEntries();

RHammingsOfs	= 0;
for(ChromID = 0; ChromID < NumEntries; ChromID++)
	{
	CurLoci = 0;
	if(bProcSepKMers)
		{
		pChrom = &m_pGChroms[ChromID];
		ChromLen = pChrom->Len;
		strcpy(szChromName,pChrom->szChromName);
		}
	else
		{
		ChromLen = m_pSfxArray->GetSeqLen(ChromID+1);
		m_pSfxArray->GetIdentName(ChromID+1,sizeof(szChromName)-1,szChromName);
		}
	pRHammings = &m_pRHammings[RHammingsOfs];
	RHammingsOfs += ChromLen;
	if(ChromLen < KMerLen)
		continue;

	CurHammingLen = 1;
	CurHamming = *pRHammings;
	CurHammingLoci = 0;
	for(CurLoci = 1; CurLoci <= (UINT32)(ChromLen - KMerLen); CurLoci++,pRHammings++)
		{
		if(*pRHammings == CurHamming)
			{
			CurHammingLen += 1;
			continue;
			}
		BuffOfs += sprintf(&szBuffer[BuffOfs],"%s\t%d\t%d\t%d\n",szChromName,CurHammingLoci,CurHammingLoci+CurHammingLen-1,CurHamming);
		CurHammingLen = 1;
		CurHamming = *pRHammings;
		CurHammingLoci = CurLoci;

		if((BuffOfs + 200) > sizeof(szBuffer))
			{
			CUtility::SafeWrite(m_hOutFile,szBuffer,BuffOfs);
			BuffOfs = 0;
			}
		}
	BuffOfs += sprintf(&szBuffer[BuffOfs],"%s\t%d\t%d\t%d\n",szChromName,CurHammingLoci,CurHammingLoci+CurHammingLen-1,CurHamming);
	}

if(BuffOfs)
	{
	CUtility::SafeWrite(m_hOutFile,szBuffer,BuffOfs);
	BuffOfs = 0;
	}
#ifdef _WIN32
_commit(m_hOutFile);
#else
fsync(m_hOutFile);
#endif
close(m_hOutFile);
m_hOutFile = -1;
return(eBSFSuccess);
}

// ReportRestrictedHammingsWiggle
int
ReportRestrictedHammingsWiggle(etSensitivity Sensitivity, // restricted hamming processing sensitivity
					int RHamm,							// restricted hammings limit
					int KMerLen,						// Hamming for this length k-mer sequences
					bool bWatsonOnly,					// true if watson strand only processing
					bool bProcSepKMers,					// true if processing probe k-mers from sequences not within target suffix array
					char *pszHammingFile)				// report into this output file
{
char szBuffer[16000];
int BuffOfs;
UINT32 CurLoci;
int ChromID;
tsGChrom *pChrom;
UINT8 *pRHammings;
int NumEntries;
int ChromLen;
char szChromName[100];
UINT32 RHammingsOfs;

	// output results
#ifdef _WIN32
m_hOutFile = open(pszHammingFile,O_CREATETRUNC );
#else
if((m_hOutFile = open(pszHammingFile,O_RDWR | O_CREAT,S_IREAD | S_IWRITE))!=-1)
	if(ftruncate(m_hOutFile,0)!=0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate %s - %s",pszHammingFile,strerror(errno));
		Reset(false);
		return(eBSFerrCreateFile);
		}
#endif

if(m_hOutFile < 0)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate output file '%s'",pszHammingFile);
	Reset(false);
	return(eBSFerrCreateFile);
	}


pChrom = &m_pGChroms[0];
pRHammings = m_pRHammings;

BuffOfs = 0;
CurLoci = 0;

if(bProcSepKMers)
	NumEntries = m_NumGChroms;
else
	NumEntries = m_pSfxArray->GetNumEntries();

BuffOfs += sprintf(&szBuffer[BuffOfs],"track type=wiggle_0 color=50,150,255 autoScale=off maxHeightPixels=128:32:8 name=\"Hammings (%d,%d,%d) - %s\" description=\"Hammings Sensitivity: %d KMerLen: %d RHamm: %d  for %s\"\n",
						Sensitivity,KMerLen,RHamm,pszHammingFile,Sensitivity,KMerLen,RHamm,pszHammingFile);

RHammingsOfs	= 0;
for(ChromID = 0; ChromID < NumEntries; ChromID++)
	{
	CurLoci = 0;
	if(bProcSepKMers)
		{
		pChrom = &m_pGChroms[ChromID];
		ChromLen = pChrom->Len;
		strcpy(szChromName,pChrom->szChromName);
		}
	else
		{
		ChromLen = m_pSfxArray->GetSeqLen(ChromID+1);
		m_pSfxArray->GetIdentName(ChromID+1,sizeof(szChromName)-1,szChromName);
		}
	pRHammings = &m_pRHammings[RHammingsOfs];
	RHammingsOfs += ChromLen;
	if(ChromLen < KMerLen)
		continue;
	BuffOfs += sprintf(&szBuffer[BuffOfs],"fixedStep chrom=%s start=%d step=1 span=1\n",szChromName,1);
	for(CurLoci = 0; CurLoci <= (UINT32)(ChromLen - KMerLen); CurLoci++,pRHammings++)
		{
		BuffOfs += sprintf(&szBuffer[BuffOfs],"%d\n",*pRHammings);
		if((BuffOfs + 500) > sizeof(szBuffer))
			{
			CUtility::SafeWrite(m_hOutFile,szBuffer,BuffOfs);
			BuffOfs = 0;
			}
		}
	}

if(BuffOfs)
	{
	CUtility::SafeWrite(m_hOutFile,szBuffer,BuffOfs);
	BuffOfs = 0;
	}
#ifdef _WIN32
_commit(m_hOutFile);
#else
fsync(m_hOutFile);
#endif
close(m_hOutFile);
m_hOutFile = -1;
return(eBSFSuccess);
}


int
GenRestrictedHammings(etSensitivity Sensitivity, // restricted hamming processing sensitivity
		etResFormat ResFormat,				// restricted Hamming output file format
		bool bWatsonOnly,					// true if watson strand only processing
		int RHamm,							// if > 0 then restricted hammings limit
		int KMerLen,						// Hamming for this length k-mer sequences
		int SampleN,						// sample (process) every Nth K-mer
		int NumThreads,						// process with upto this number of threads
		char *pszSfxFile,					// assembly suffix array file
		char *pszInSeqFile,					// optionally kmer sequences from this file
		char *pszHammingFile)				// write Hamming edit distances into this file
{
int Rslt;
bool bProcSepKMers;
UINT32 Len;
INT64 TotLen;
int ChromID;
int ChromIdx;
tsGChrom *pChrom;
etSeqBase *pSeq;
etSeqBase *pGseq;

m_NumThreads = NumThreads;
m_PerThreadAllocdIdentNodes = cMaxNumIdentNodes;
m_TotAllocdIdentNodes = m_PerThreadAllocdIdentNodes * m_NumThreads;
if((m_pAllocsIdentNodes = new tsIdentNode [m_TotAllocdIdentNodes])==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: unable to allocate memory for %d tsIdentNodes",m_TotAllocdIdentNodes);
	Reset(false);
	return(eBSFerrMem);
	}

// load targeted suffix array assembly
// open bioseq file containing suffix array for targeted assembly to align reads against
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Loading suffix array file '%s'", pszSfxFile);
if((m_pSfxArray = new CSfxArrayV3()) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CSfxArrayV2");
	Reset(false);
	return(eBSFerrObj);
	}
if((Rslt=m_pSfxArray->Open(pszSfxFile,false,false,false))!=eBSFSuccess)
	{
	while(m_pSfxArray->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,m_pSfxArray->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open input bioseq suffix array file '%s'",pszSfxFile);
	Reset(false);
	return(Rslt);
	}

// report to user some sfx array metadata as conformation the targeted assembly is correct
strcpy(m_szTargSpecies,m_pSfxArray->GetDatasetName());

tsSfxHeaderV3 SfxHeader;
m_pSfxArray->GetSfxHeader(&SfxHeader);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Genome Assembly Name: '%s' Descr: '%s' Title: '%s' Version: %d",
					 m_szTargSpecies,SfxHeader.szDescription,SfxHeader.szTitle,SfxHeader.Version);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Assembly has blocks: %d, max block size: %llu",SfxHeader.NumSfxBlocks,SfxHeader.SfxBlockSize);

// load single SfxBlock, expected to contain all chromosomes, and sample K-mers against that block
int CurBlockID = 1;
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Loading genome assembly suffix array...");
if((Rslt=m_pSfxArray->SetTargBlock(CurBlockID))<0)
	{
	while(m_pSfxArray->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,m_pSfxArray->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: unable to load genome assembly suffix array");
	return(Rslt);
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Genome assembly suffix array loaded");

if(pszInSeqFile != NULL && pszInSeqFile[0] != '\0')
	bProcSepKMers = true;
else
	bProcSepKMers = false;

// if specified that kmer sequences are to be separately loaded then load from that file
if(bProcSepKMers)
	{
	if((m_pBioSeqFile = new CBioSeqFile) == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CBioSeqFile object");
		return(eBSFerrObj);
		}
	if((Rslt = m_pBioSeqFile->Open(pszInSeqFile))!=eBSFSuccess)
		{
		while(m_pBioSeqFile->NumErrMsgs())
			gDiagnostics.DiagOut(eDLFatal,gszProcName,m_pBioSeqFile->GetErrMsg());
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open assembly sequence file '%s'",pszInSeqFile);
		Reset(false);
		return(Rslt);
		}

	m_pBioSeqFile->GetTitle(sizeof(m_szSpecies),m_szSpecies);

	m_NumGChroms = m_pBioSeqFile->NumEntries();
	if((m_pGChroms = new tsGChrom [m_NumGChroms])==NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"unable to allocate memory (%d bytes) for %d ChromSeqs",sizeof(tsGChrom) * m_NumGChroms,m_NumGChroms);
		Reset(false);
		return(eBSFerrMem);
		}
	memset(m_pGChroms,0,sizeof(tsGChrom) * m_NumGChroms);

	// determine total sequence length required for allocating to hold genome as one concatenated sequence
	// concatenated sequence starts with eBaseEOG followed be each chromosome sequence separated from the
	// next chromosome sequence with eBaseEOS, and finally the last chromosome is terminated by eBaseEOG, not eBaseEOS
	m_GenomeLen = 1;		// allow for initial eBaseEOG marker
	ChromID = 0;
	ChromIdx = 0;
	while((ChromID = m_pBioSeqFile->Next(ChromID)) > 0) {
		pChrom = &m_pGChroms[ChromID-1];
		pChrom->ChromID = ChromID;
		pChrom->SeqOfs = m_GenomeLen;
		m_pBioSeqFile->GetName(ChromID,sizeof(pChrom->szChromName),pChrom->szChromName);
		Len = m_pBioSeqFile->GetDataLen(ChromID);
		m_GenomeLen += Len + 1;		// allow for a concatenation separator (eBaseEOS) or final (eBaseEOG)
		pChrom->Len = Len;
		pChrom->pSeq = NULL;
		pChrom->ChromSeqID = ++ChromIdx;
		}

	m_AllocGenomeSeq = m_GenomeLen;
	#ifdef _WIN32
	// use malloc instead of new() because of gnu new()/malloc issues - helps portability
	m_pGenomeSeq = (UINT8 *) malloc((size_t)m_AllocGenomeSeq);
	if(m_pGenomeSeq == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadGenome: Memory allocation of %lld bytes failed",(INT64)m_AllocGenomeSeq);
		Reset(false);
		return(eBSFerrMem);
		}
	#else
	// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations so use mmap
	m_pGenomeSeq = (UINT8 *)mmap(NULL,(size_t)m_AllocGenomeSeq, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
	if(m_pGenomeSeq == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadGenome: Memory allocation of %lld bytes through mmap()  failed",(INT64)m_AllocGenomeSeq,strerror(errno));
		m_pGenomeSeq = NULL;
		Reset(false);
		return(eBSFerrMem);
		}
	#endif

	// now initialise m_pGenomeSeq with the concatenated sequences
	ChromID = 0;
	pSeq = m_pGenomeSeq;
	*pSeq++ = eBaseEOG;		// start with this marker so subsequent crick strand processing functions 'know' when all sequence bases processed
	TotLen = 0;
	while((ChromID = m_pBioSeqFile->Next(ChromID)) > 0) {
		pChrom = &m_pGChroms[ChromID-1];
		pChrom->pSeq = pSeq;
		TotLen += pChrom->Len;
		m_pBioSeqFile->GetData(ChromID,eSeqBaseType,0,pSeq,pChrom->Len);
		pGseq = pSeq;
		pSeq += pChrom->Len;
		*pSeq++ = eBaseEOS;		// mark end of current chrom
		while(*pGseq++ != eBaseEOS)		// processing expects repeat masking removed
			pGseq[-1] &= ~ cRptMskFlg;
		}
	pSeq[-1] = eBaseEOG;   // marks prev chromosome sequence as being the last
	delete m_pBioSeqFile;
	m_pBioSeqFile = NULL;

	// determine the total number of subsequences over which Hammings are required
	int ChromLen;
	m_NumSubSeqs = 0;
	m_TotAllocHammings = 0;
	pChrom = m_pGChroms;
	for(ChromIdx = 0; ChromIdx < m_NumGChroms; ChromIdx++, pChrom++)
		{
		m_TotAllocHammings += pChrom->Len;
		ChromLen = pChrom->Len - (m_SubSeqLen-1);
		if(ChromLen < 1)
			{
			pChrom->NumSubSeqs = 0;
			continue;
			}
		pChrom->NumSubSeqs = (UINT32)ChromLen;
		m_NumSubSeqs += (UINT32)ChromLen;
		}
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Targeted suffix genome has %d sequences containing %llu total nucleotides, %u sequences of at least K-mer length %d",
														m_NumGChroms,(UINT64)m_TotAllocHammings,m_NumSubSeqs,m_SubSeqLen);


	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Genome containing %llu total nucleotides loaded with %u subsequences of K-mer length %d...",m_TotAllocHammings,m_NumSubSeqs,m_SubSeqLen);
	}
else
	{
	// determine the total number of subsequences over which Hammings are required
	int EntryID;
	UINT32 ChromLen;
	int NumEntries = m_pSfxArray->GetNumEntries();
	m_NumSubSeqs = 0;
	m_TotAllocHammings = 0;
	for(EntryID = 1; EntryID <= NumEntries; EntryID++)
		{
		ChromLen = m_pSfxArray->GetSeqLen(EntryID);
		m_TotAllocHammings += ChromLen;
		ChromLen -= (m_SubSeqLen-1);
		if(ChromLen < 1)					// don't bother with chromosome too short to contain a subseq
			continue;
		m_NumSubSeqs += (UINT32)ChromLen;
		}
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Targeted suffix genome has %d sequences containing %llu total nucleotides, %u sequences of at least K-mer length %d",
														NumEntries,(UINT64)m_TotAllocHammings,m_NumSubSeqs,m_SubSeqLen);
	}

m_pRHammings = NULL;

#ifdef _WIN32
	// use malloc instead of new() because of gnu new()/malloc issues - helps portability
m_pRHammings = (UINT8 *) malloc((size_t)m_TotAllocHammings);
if(m_pRHammings == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadGenome: Memory allocation of %llu bytes failed",(INT64)m_TotAllocHammings);
	Reset(false);
	return(eBSFerrMem);
	}
#else
	// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations so use mmap
m_pRHammings = (UINT8 *)mmap(NULL,(size_t)m_TotAllocHammings, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
if(m_pRHammings == MAP_FAILED)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadGenome: Memory allocation of %llu bytes through mmap()  failed",(INT64)m_TotAllocHammings,strerror(errno));
	m_pRHammings = NULL;
	Reset(false);
	return(eBSFerrMem);
	}
#endif
memset(m_pRHammings,0x0ff,m_TotAllocHammings);   // 0x0ff used to show that no hamming was generated at m_pRHammings[loci], sequence too short ...

tsHammProcState HammProcState;
tsRHThread RHThreads[cMaxWorkerThreads];
tsWorkerThreads WorkerThreads[cMaxWorkerThreads];
tsRHThread *pThreadParam;
tsWorkerThreads *pWorkerThread;
int ThreadIdx;

int MinCoreDepth;
int MaxCoreDepth;

switch(Sensitivity) {
	case eSensLess:
		MinCoreDepth = 200;
		MaxCoreDepth = 5000;
		break;

	case eSensDefault:
		MinCoreDepth = 1000;
		MaxCoreDepth = 10000;
		break;

	case eSensMore:
		MinCoreDepth = 2000;
		MaxCoreDepth = 20000;
		break;

	case eSensUltra:
		MinCoreDepth = 5000;
		MaxCoreDepth = 50000;
		break;
	}

HammProcState.CurEntryID = 0;
HammProcState.CurEntryLen = 0;
HammProcState.CurSeqLen = 0;
HammProcState.CurSeqOfs = 0;
HammProcState.EntryHammOfs = 0;
HammProcState.NumEntryIDs = bProcSepKMers ? m_NumGChroms : m_pSfxArray->GetNumEntries();
HammProcState.pHammings = m_pRHammings;
HammProcState.pSeq = NULL;
HammProcState.bProcSepKMers = bProcSepKMers;
HammProcState.bWatsonOnly = bWatsonOnly;
HammProcState.MaxHammSeqLen = (int)max(10000,(int)(m_TotAllocHammings/(size_t)10000));
HammProcState.KMerLen = KMerLen;
HammProcState.SampleN = SampleN;
HammProcState.RHamm = RHamm;
HammProcState.ApproxNumChroms = 0;
HammProcState.ApproxNumHammings = 0;
HammProcState.MinCoreDepth = MinCoreDepth;
HammProcState.MaxCoreDepth = MaxCoreDepth;

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Starting restricted Hamming edit distance processing");
pThreadParam = &RHThreads[0];
pWorkerThread = WorkerThreads;
for(ThreadIdx = 0; ThreadIdx < m_NumProcThreads; ThreadIdx++,pThreadParam++,pWorkerThread++)
	{
	pThreadParam->threadID = ThreadIdx+1;
	pThreadParam->EntryID = 0;
	pThreadParam->pHammings = NULL;
	pThreadParam->pHammState = &HammProcState;
	pThreadParam->NumIdentNodes = m_PerThreadAllocdIdentNodes;
	pThreadParam->pIdentNodes = &m_pAllocsIdentNodes[m_PerThreadAllocdIdentNodes * ThreadIdx];
	pThreadParam->State = 0;
#ifdef _WIN32
	pWorkerThread->threadHandle = (HANDLE)_beginthreadex(NULL,0x0fffff,RestrictedHammingThread,pThreadParam,0,&WorkerThreads[ThreadIdx].threadID);
#else
	pWorkerThread->threadRslt =	pthread_create (&pWorkerThread->threadID , NULL , RestrictedHammingThread , pThreadParam );
#endif
	}

// wait for all threads to have completed
UINT32 ApproxNumChroms;
INT64 ApproxNumHammings;
UINT32 PrevApproxNumChroms = 0;
INT64 PrevApproxNumHammings = 0;

// allow threads 10 seconds to startup
#ifdef _WIN32
Sleep(10000);
#else
sleep(10);
#endif

#ifdef _WIN32
WaitForSingleObject(m_hMtxThreadParams,INFINITE);
PrevApproxNumChroms = HammProcState.ApproxNumChroms;
PrevApproxNumHammings = HammProcState.ApproxNumHammings;
ReleaseMutex(m_hMtxThreadParams);
#else
pthread_mutex_lock(&m_hMtxThreadParams);
PrevApproxNumChroms = HammProcState.ApproxNumChroms;
PrevApproxNumHammings = HammProcState.ApproxNumHammings;
pthread_mutex_unlock(&m_hMtxThreadParams);
#endif
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Progress: %llu K-mers processed from %d sequences",PrevApproxNumHammings,PrevApproxNumChroms);

for(ThreadIdx = 0; ThreadIdx < m_NumProcThreads; ThreadIdx++)
	{
#ifdef _WIN32
	while(WAIT_TIMEOUT == WaitForSingleObject( WorkerThreads[ThreadIdx].threadHandle, 60000))
		{
		WaitForSingleObject(m_hMtxThreadParams,INFINITE);
		ApproxNumChroms = HammProcState.ApproxNumChroms;
		ApproxNumHammings = HammProcState.ApproxNumHammings;
		ReleaseMutex(m_hMtxThreadParams);
		if(ApproxNumChroms > PrevApproxNumChroms || ApproxNumHammings > PrevApproxNumHammings)
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"Progress: %llu K-mers processed from %d sequences",ApproxNumHammings,ApproxNumChroms);
		PrevApproxNumChroms = ApproxNumChroms;
		PrevApproxNumHammings = ApproxNumHammings;
		}
	CloseHandle( WorkerThreads[ThreadIdx].threadHandle);
#else
	struct timespec ts;
	int JoinRlt;
	clock_gettime(CLOCK_REALTIME, &ts);
	ts.tv_sec += 60;
	while((JoinRlt = pthread_timedjoin_np(WorkerThreads[ThreadIdx].threadID, NULL, &ts)) != 0)
		{
		pthread_mutex_lock(&m_hMtxThreadParams);
		ApproxNumChroms = HammProcState.ApproxNumChroms;
		ApproxNumHammings = HammProcState.ApproxNumHammings;
		pthread_mutex_unlock(&m_hMtxThreadParams);
		if(ApproxNumChroms > PrevApproxNumChroms || ApproxNumHammings > PrevApproxNumHammings)
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"Progress: %llu K-mers  processed from %d sequences",ApproxNumHammings,ApproxNumChroms);
		PrevApproxNumChroms = ApproxNumChroms;
		PrevApproxNumHammings = ApproxNumHammings;
		ts.tv_sec += 60;
		}
	pthread_join(WorkerThreads[ThreadIdx].threadID,NULL);
#endif
	}

#ifdef _WIN32
CloseHandle(m_hMtxThreadParams);
#else
pthread_mutex_destroy(&m_hMtxThreadParams);
#endif

ApproxNumChroms = HammProcState.ApproxNumChroms;
ApproxNumHammings = HammProcState.ApproxNumHammings;

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Progress: %llu K-mers processed from %d sequences",ApproxNumHammings,ApproxNumChroms);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Completed restricted hammings generation...");

UINT32 CurLoci;
//tsGChrom *pChrom;
UINT8 *pRHammings;
UINT32 SeqIdx;
int NumEntries;
UINT32 ChromLen;
UINT32 RHammingsOfs;

if(pszHammingFile != NULL && pszHammingFile[0] != '\0')
	{
	// write out the Hamming distances as chrom,loci,Hamming - a simple format easily post processed
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Writing Hamming edit distances to file: '%s'",pszHammingFile);

	switch(ResFormat) {
		case cRHFcsv:				// CSV with Hamming loci ranges (chrom, startloci, length, Hamming)
			Rslt = ReportRestrictedHammingsCSV(KMerLen,			// Hamming for this length k-mer sequences
						bProcSepKMers,							// true if processing probe k-mers from sequences not within target suffix array
						pszHammingFile);						// report into this output file
			break;

		case cRHFBed:
			Rslt = ReportRestrictedHammingsBed(Sensitivity,RHamm,KMerLen,bWatsonOnly,bProcSepKMers,pszHammingFile);
			break;

		case cRHFWiggle:
			Rslt = ReportRestrictedHammingsWiggle(Sensitivity,RHamm,KMerLen,bWatsonOnly,bProcSepKMers,pszHammingFile);
			break;

		}
	}

UINT64 CntHist[256];
memset(CntHist,0,sizeof(CntHist));

CurLoci = 0;

if(bProcSepKMers)
	NumEntries = m_NumGChroms;
else
	NumEntries = m_pSfxArray->GetNumEntries();

RHammingsOfs = 0;
UINT32 NumSampled = 0;
UINT32 NumUnsampled = 0;
for(ChromID = 0; ChromID < NumEntries; ChromID++)
	{
	CurLoci = 0;
	if(bProcSepKMers)
		{
		pChrom = &m_pGChroms[ChromID];
		ChromLen = pChrom->Len;
		}
	else
		ChromLen = m_pSfxArray->GetSeqLen(ChromID+1);
	pRHammings = &m_pRHammings[RHammingsOfs];
	RHammingsOfs += ChromLen;
	if(ChromLen < (UINT32)KMerLen)
		continue;

	for(CurLoci = 0; CurLoci <= (UINT32)(ChromLen - (UINT32)KMerLen); CurLoci++,pRHammings++)
		{
		if(*pRHammings < 0x7f)
			{
			NumSampled += 1;
			CntHist[*pRHammings] += 1;
			}
		else
			NumUnsampled += 1;
		}
	}
if(*pRHammings < 0x7f)
	{
	CntHist[*pRHammings] += 1;	
	NumSampled += 1;
	}
else
	NumUnsampled += 1;

// report to log the basic count distribution.....
if(SampleN > 1)
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Sampled %u K-mers, did not sample %u",NumSampled,NumUnsampled);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Distribution:\nEditDist,Freq,Proportion");
for(SeqIdx = 0; SeqIdx <= (UINT32)(RHamm+1); SeqIdx++)
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"%d%c,%llu,%1.3f",SeqIdx,SeqIdx == RHamm+1 ? '+' : ' ',CntHist[SeqIdx],(CntHist[SeqIdx]*100.0f)/NumSampled);

Reset(true);

delete m_pSfxArray;
m_pSfxArray = NULL;
Reset(true);
return(0);
}


int
Process(etPMode PMode,			// processing mode
		bool bWatsonOnly,		// true if watson strand only processing
		etSensitivity Sensitivity, // restricted hamming processing sensitivity
		etResFormat ResFormat,	// restricted Hamming output file format
		int RHamm,			    // if > 0 then restricted hammings limit
		UINT32 SweepStart,		// start processing from this sweep instance (1..GenomeLen) or if distributed processing then the total number of nodes
		UINT32 SweepEnd,		// finish processing at this sweep instance (0 for all remaining or SweepStart..GenomeLen) or if distributed processing then the node instance
		int SeqLen,				// Hamming for this length sequences
		int SampleN,			// sample (process) every N sweep instance (or if restricted Hammings then every Nth K-mer)
		int NumThreads,			// process with upto this number of threads
		char *pszGenomeFile,	// bioseq file containing targeted genome assembly, or if restricted hammings then the assembly suffix array file
		char *pszInSeqFile,		// if restricted Hammings then kmer sequences from this file
		char *pszHammingFile)	// write Hamming edit distances into this file
{
int Rslt;
UINT32 SSeqStart;
UINT32 SSeqEnd;

Init();


#ifdef _WIN32
if((m_hMtxThreadParams = CreateMutex(NULL,false,NULL))==NULL)
#else
if(pthread_mutex_init (&m_hMtxThreadParams,NULL)!=0)
#endif
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: unable to create mutex");
	Reset(false);
	return(eBSFerrInternal);
	}

m_PMode = PMode;
m_SubSeqLen = SeqLen;
m_NumProcThreads = NumThreads;
m_KMerLen = SeqLen;
m_bWatsonOnly = bWatsonOnly;
m_RHamm = RHamm;
m_Sensitivity = Sensitivity;
m_SampleN = SampleN;

if(PMode == ePMrestrict)
	return(GenRestrictedHammings(Sensitivity,ResFormat,bWatsonOnly,RHamm,SeqLen,SampleN,NumThreads,pszGenomeFile,pszInSeqFile,pszHammingFile));

if(PMode == ePMmerge)
	return(MergeHammings(pszGenomeFile,pszHammingFile));

if(PMode == ePMtrans)
	return(TransHammings(pszGenomeFile,pszHammingFile));

if(PMode == ePMtransCSV)
	return(TransHammingsCSV(pszGenomeFile,pszHammingFile));

if((Rslt=LoadGenome(pszGenomeFile))!=eBSFSuccess)
	{
	Reset(false);
	return(Rslt);
	}

if(PMode == ePMdist)
	{
	UINT32 L1;
	UINT32 L2;
	if(!bWatsonOnly)
		{
		L1 = EstL(m_GenomeLen,SweepEnd,SweepStart);
		L2 = EstL(m_GenomeLen,SweepEnd-1,SweepStart);
		}
	else
		{
		UINT64 Area = ((UINT64)m_GenomeLen * (UINT64)m_GenomeLen)/2;
		L1 = (UINT32)sqrt((2.0f * SweepEnd * Area)/SweepStart);
		L2 = (UINT32)sqrt((2.0f * (SweepEnd - 1) * Area)/SweepStart) - 1;
		}
	SSeqStart = m_GenomeLen - L1;
	SSeqEnd = m_GenomeLen - L2;
	// adjust SweepStart (down) and SweepEnd (up) to allow for the presence of interchromosome EOS markers and rounding errors
	if(SSeqStart > (UINT32)(2 * m_NumGChroms))
		SSeqStart -= (2 * m_NumGChroms);
	else
		SSeqStart = 0;
	SSeqEnd += 10 + (2 * m_NumGChroms); // allowance for EOS markers plus safty margin to ensure node hammings do overlap
	if(SSeqEnd > m_GenomeLen)
		SSeqEnd = m_GenomeLen;
	if(SSeqStart > 10)				   // safty margin to ensure node hammings do overlap
		SSeqStart -= 10;
	else
	    SSeqStart = 1;
	}
else	// single node processing
	{
	// silently clamp SweepStart and SweepEnd to be no more than the actual genome length
	// SweepStart must be <= to the genome length and SweepEnd must be 0 or <= genome length
	if(SweepStart > m_GenomeLen)
		SSeqStart = m_GenomeLen;
	else
		SSeqStart = SweepStart;

	if(SweepEnd == 0)
		SSeqEnd = m_GenomeLen;
	else
		if(SweepEnd > m_GenomeLen)
			SSeqEnd = m_GenomeLen;
		else
			SSeqEnd = SweepEnd;
	}

// log the actual sweep start/end
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Node sweep start is %u, and sweep end is %u",SSeqStart,SSeqEnd);

// adjust m_NumProcThreads so as to ensure each processing thread has at least 1 sweep to process
m_NumProcThreads = (int)min(1 + SSeqEnd - SSeqStart,(UINT32)m_NumProcThreads);

if(pszHammingFile != NULL && pszHammingFile[0] != '\0')
	{
	// ensure output results file can be accessed/created
#ifdef _WIN32
	m_hOutFile = open(pszHammingFile,O_CREATETRUNC );
#else
	if((m_hOutFile = open(pszHammingFile,O_RDWR | O_CREAT,S_IREAD | S_IWRITE))!=-1)
		if(ftruncate(m_hOutFile,0)!=0)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate %s - %s",pszHammingFile,strerror(errno));
				Reset(false);
				return(eBSFerrCreateFile);
				}
#endif

	if(m_hOutFile < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate output file '%s'",pszHammingFile);
		Reset(false);
		return(eBSFerrCreateFile);
		}
	close(m_hOutFile);
	m_hOutFile = -1;
	}

m_PercentComplete = 0.0;
tsThreadParams *pThreadParam;
UINT32 BeginSeq;
UINT32 EndSeq;
UINT64 HamDistOfs;
UINT32 NumSubSeqs;
int ThreadIdx;
HamDistOfs = 0;
pThreadParam = m_pThreadParams;

BeginSeq = SSeqStart;
EndSeq = SSeqEnd;
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Starting Hamming edit distance processing for %d subsequences",1 + EndSeq - BeginSeq);
for(ThreadIdx = 0; ThreadIdx < m_NumProcThreads; ThreadIdx++,pThreadParam++)
	{
	pThreadParam->ThreadID = ThreadIdx+1;
	pThreadParam->bWatsonOnly = bWatsonOnly;
	NumSubSeqs = (EndSeq - (SSeqStart - 1))/(m_NumProcThreads - ThreadIdx);
	pThreadParam->SSofs = BeginSeq;
	pThreadParam->SeqDelta = 1;
	pThreadParam->NumSeqs = NumSubSeqs;
	pThreadParam->HamDistOfs = HamDistOfs;
	pThreadParam->SampleN = SampleN;
	EndSeq -= NumSubSeqs;
	BeginSeq += NumSubSeqs;
	HamDistOfs += (m_GenomeLen-2);		// initial/final eBaseEOG do not have hammings
	pThreadParam->State = 0;
	}


gDiagnostics.DiagOut(eDLInfo,gszProcName,"Starting Hamming edit distance processing");
tsWorkerThreads WorkerThreads[cMaxWorkerThreads];

for(ThreadIdx = 0; ThreadIdx < m_NumProcThreads; ThreadIdx++)
	{
	WorkerThreads[ThreadIdx].Rslt = 0;
#ifdef _WIN32
	WorkerThreads[ThreadIdx].threadHandle = (HANDLE)_beginthreadex(NULL,0x0fffff,ThreadedGHamDist,&WorkerThreads[ThreadIdx],0,&WorkerThreads[ThreadIdx].threadID);
#else
	WorkerThreads[ThreadIdx].threadRslt =	pthread_create (&WorkerThreads[ThreadIdx].threadID , NULL , ThreadedGHamDist , &WorkerThreads[ThreadIdx] );
#endif
	}

// allow threads 10 seconds to startup
#ifdef _WIN32
Sleep(10000);
#else
sleep(10);
#endif

double CurPercentComplete;
double PrevCurPercentComplete = -1.0;

#ifdef _WIN32
WaitForSingleObject(m_hMtxThreadParams,INFINITE);
#else
pthread_mutex_lock(&m_hMtxThreadParams);
#endif
PrevCurPercentComplete = m_PercentComplete;
#ifdef _WIN32
ReleaseMutex(m_hMtxThreadParams);
#else
pthread_mutex_unlock(&m_hMtxThreadParams);
#endif

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Progress: %1.2f%% completed",PrevCurPercentComplete);

// wait for all threads to have completed and every 5 min give user an indication of progress
for(ThreadIdx = 0; ThreadIdx < m_NumProcThreads; ThreadIdx++)
	{
#ifdef _WIN32
	while(WAIT_TIMEOUT == WaitForSingleObject( WorkerThreads[ThreadIdx].threadHandle, 60000 * 10))
		{
		WaitForSingleObject(m_hMtxThreadParams,INFINITE);
		CurPercentComplete = m_PercentComplete;
		ReleaseMutex(m_hMtxThreadParams);
		if(CurPercentComplete > PrevCurPercentComplete)
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"Progress: %1.2f%% completed",CurPercentComplete);
		PrevCurPercentComplete = CurPercentComplete;
		}
	CloseHandle( WorkerThreads[ThreadIdx].threadHandle);
#else
	struct timespec ts;
	int JoinRlt;
	clock_gettime(CLOCK_REALTIME, &ts);
	ts.tv_sec += 60;
	while((JoinRlt = pthread_timedjoin_np(WorkerThreads[ThreadIdx].threadID, NULL, &ts)) != 0)
		{
		pthread_mutex_lock(&m_hMtxThreadParams);
		CurPercentComplete = m_PercentComplete;
		pthread_mutex_unlock(&m_hMtxThreadParams);
		if(CurPercentComplete > PrevCurPercentComplete)
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"Progress: %1.2f%% completed",CurPercentComplete);
		PrevCurPercentComplete = CurPercentComplete;
		ts.tv_sec += 60;
		}
	pthread_join(WorkerThreads[ThreadIdx].threadID,NULL);
#endif
	}

#ifdef _WIN32
CloseHandle(m_hMtxThreadParams);
#else
pthread_mutex_destroy(&m_hMtxThreadParams);
#endif
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Progress: %1.2f%% completed",m_PercentComplete);

// now merge the Hammings

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Merging Hamming edit distances...");

int Idx;
UINT32 SeqIdx;
UINT16 *pHamDist1;
UINT16 *pHamDist2;

if(m_NumProcThreads > 1)
	{
	for(Idx=1; Idx < m_NumProcThreads; Idx++)
		{
		pHamDist1 = m_pHamDist;
		pHamDist2 = &m_pHamDist[(INT64)Idx * (m_GenomeLen-2)];		// no hammings for initial and final eBaseEOG markers)
		for(SeqIdx = 0; SeqIdx < (m_GenomeLen-2); SeqIdx++,pHamDist1++,pHamDist2++)
			{
			if(*pHamDist2 < *pHamDist1)
				*pHamDist1 = *pHamDist2;
			}
		}
	}

char szBuffer[16000];
int BuffOfs;
UINT32 CurLoci;
tsGChrom *pChrom;

if(pszHammingFile != NULL && pszHammingFile[0] != '\0')
	{
	// write out the Hamming distances as chrom,loci,Hamming - a simple format easily post processed
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Writing Hamming edit distances to file: '%s'",pszHammingFile);

	// output results
#ifdef _WIN32
	m_hOutFile = open(pszHammingFile,O_CREATETRUNC );
#else
	if((m_hOutFile = open(pszHammingFile,O_RDWR | O_CREAT,S_IREAD | S_IWRITE))!=-1)
		if(ftruncate(m_hOutFile,0)!=0)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate %s - %s",pszHammingFile,strerror(errno));
				Reset(false);
				return(eBSFerrCreateFile);
				}
#endif

	if(m_hOutFile < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate output file '%s'",pszHammingFile);
		Reset(false);
		return(eBSFerrCreateFile);
		}

	pChrom = &m_pGChroms[0];
	pHamDist1 = m_pHamDist;
	BuffOfs = 0;
	CurLoci = 0;
	BuffOfs = sprintf(szBuffer,"%u,%d,%d\n",m_GenomeLen,SSeqStart+1,SSeqEnd);
	for(SeqIdx = 0; SeqIdx < (m_GenomeLen-2); SeqIdx++,pHamDist1++,CurLoci++)
		{
		if(CurLoci >= pChrom->NumSubSeqs)
			{
			if(pChrom->ChromSeqID == m_NumGChroms)
				break;
			pChrom += 1;
			CurLoci = 0;
			SeqIdx += m_SubSeqLen;
			pHamDist1 += m_SubSeqLen;
			}
		if(*pHamDist1 <= SeqLen)	// if Hamming has been generated for this loci...
			{
			BuffOfs += sprintf(&szBuffer[BuffOfs],"\"%s\",%d,%d\n",pChrom->szChromName,CurLoci,*pHamDist1);
			if((BuffOfs + 200) > sizeof(szBuffer))
				{
				CUtility::SafeWrite(m_hOutFile,szBuffer,BuffOfs);
				BuffOfs = 0;
				}
			}
		}
	if(BuffOfs)
		{
		CUtility::SafeWrite(m_hOutFile,szBuffer,BuffOfs);
		BuffOfs = 0;
		}
#ifdef _WIN32
	_commit(m_hOutFile);
#else
	fsync(m_hOutFile);
#endif
	close(m_hOutFile);
	m_hOutFile = -1;
	}

// log basic count distribution.....
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Distribution:\nEditDist,Freq,Proportion");
int CntHist[cMaxSeqLen+1];
memset(CntHist,0,sizeof(CntHist));
pChrom = &m_pGChroms[0];
pHamDist1 = m_pHamDist;
BuffOfs = 0;
CurLoci = 0;
for(SeqIdx = 0; SeqIdx < (m_GenomeLen-2); SeqIdx++,pHamDist1++,CurLoci++)
	{
	if(CurLoci >= pChrom->NumSubSeqs)
		{
		if(pChrom->ChromSeqID == m_NumGChroms)
			break;
		pChrom += 1;
		CurLoci = 0;
		SeqIdx += m_SubSeqLen;
		pHamDist1 += m_SubSeqLen;
		}
	CntHist[*pHamDist1] += 1;
	}

for(SeqIdx = 0; SeqIdx < (UINT32)min(66,m_SubSeqLen); SeqIdx++)
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"%d,%d,%1.3f",SeqIdx,CntHist[SeqIdx],(CntHist[SeqIdx]*100.0f)/m_NumSubSeqs);

Reset(true);
return(eBSFSuccess);
}

int
MinHamCnt(int MaxHamCnt,			// only cnt upto this many mismatches
		   etSeqBase *pSeq1,		// determine Hamming edit distance between this sequence
		   etSeqBase *pSeq2)		// and this sequence
{
int Idx;
int CurCnt=0;
for(Idx = 0; Idx < (int)m_SubSeqLen; Idx++)
	if(*pSeq1++ != *pSeq2++ && CurCnt++ >= MaxHamCnt)
		return(CurCnt);
return(CurCnt);
}

int
LoadGenome(char *pszBioSeqFile) // load genome from this file which must be as a bioseq file
{
int Rslt;
int Len;
INT64 TotLen;
int ChromID;
int ChromIdx;
tsGChrom *pChrom;
etSeqBase *pSeq;
etSeqBase *pGseq;


if((m_pBioSeqFile = new CBioSeqFile) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CBioSeqFile object");
	return(eBSFerrObj);
	}
if((Rslt = m_pBioSeqFile->Open(pszBioSeqFile))!=eBSFSuccess)
	{
	while(m_pBioSeqFile->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,m_pBioSeqFile->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open assembly sequence file '%s'",pszBioSeqFile);
	Reset(false);
	return(Rslt);
	}

m_pBioSeqFile->GetTitle(sizeof(m_szSpecies),m_szSpecies);

m_NumGChroms = m_pBioSeqFile->NumEntries();
if((m_pGChroms = new tsGChrom [m_NumGChroms])==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"unable to allocate memory (%d bytes) for %d ChromSeqs",sizeof(tsGChrom) * m_NumGChroms,m_NumGChroms);
	Reset(false);
	return(eBSFerrMem);
	}
memset(m_pGChroms,0,sizeof(tsGChrom) * m_NumGChroms);

// determine total sequence length required for allocating to hold genome as one concatenated sequence
// concatenated sequence starts with eBaseEOG followed be each chromosome sequence separated from the
// next chromosome sequence with eBaseEOS, and finally the last chromosome is terminated by eBaseEOG, not eBaseEOS
m_GenomeLen = 1;		// allow for initial eBaseEOG marker
ChromID = 0;
ChromIdx = 0;
while((ChromID = m_pBioSeqFile->Next(ChromID)) > 0) {
	pChrom = &m_pGChroms[ChromID-1];
	pChrom->ChromID = ChromID;
	pChrom->SeqOfs = m_GenomeLen;
	m_pBioSeqFile->GetName(ChromID,sizeof(pChrom->szChromName),pChrom->szChromName);
	Len = m_pBioSeqFile->GetDataLen(ChromID);
	m_GenomeLen += Len + 1;		// allow for a concatenation separator (eBaseEOS) or final (eBaseEOG)
	pChrom->Len = Len;
	pChrom->pSeq = NULL;
	pChrom->ChromSeqID = ++ChromIdx;
	}

m_AllocGenomeSeq = m_GenomeLen;
#ifdef _WIN32
// use malloc instead of new() because of gnu new()/malloc issues - helps portability
m_pGenomeSeq = (UINT8 *) malloc((size_t)m_AllocGenomeSeq);
if(m_pGenomeSeq == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadGenome: Memory allocation of %lld bytes failed",(INT64)m_AllocGenomeSeq);
	Reset(false);
	return(eBSFerrMem);
	}
#else
// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations so use mmap
m_pGenomeSeq = (UINT8 *)mmap(NULL,(size_t)m_AllocGenomeSeq, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
if(m_pGenomeSeq == MAP_FAILED)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadGenome: Memory allocation of %lld bytes through mmap()  failed",(INT64)m_AllocGenomeSeq,strerror(errno));
	m_pGenomeSeq = NULL;
	Reset(false);
	return(eBSFerrMem);
	}
#endif

// now initialise m_pGenomeSeq with the concatenated sequences
ChromID = 0;
pSeq = m_pGenomeSeq;
*pSeq++ = eBaseEOG;		// start with this marker so subsequent crick strand processing functions 'know' when all sequence bases processed
TotLen = 0;
while((ChromID = m_pBioSeqFile->Next(ChromID)) > 0) {
	pChrom = &m_pGChroms[ChromID-1];
	pChrom->pSeq = pSeq;
	TotLen += pChrom->Len;
	m_pBioSeqFile->GetData(ChromID,eSeqBaseType,0,pSeq,pChrom->Len);
	pGseq = pSeq;
	pSeq += pChrom->Len;
	*pSeq++ = eBaseEOS;		// mark end of current chrom
	while(*pGseq++ != eBaseEOS)		// processing expects repeat masking removed
		pGseq[-1] &= ~ cRptMskFlg;
	}
pSeq[-1] = eBaseEOG;   // marks prev chromosome sequence as being the last
delete m_pBioSeqFile;
m_pBioSeqFile = NULL;

// determine the total number of subsequences over which Hammings are required
m_NumSubSeqs = 0;
pChrom = m_pGChroms;
int ChromLen;
for(ChromIdx = 0; ChromIdx < m_NumGChroms; ChromIdx++, pChrom++)
	{
	ChromLen = pChrom->Len + 1 - m_SubSeqLen;
	if(ChromLen < 1) // don't bother with chromosome too short to contain a subseq
		{
		pChrom->NumSubSeqs = 0;
		continue;
		}
	pChrom->NumSubSeqs = (UINT32)ChromLen;
	m_NumSubSeqs += ChromLen;
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Genome containing %llu total nucleotides loaded with %u subsequences of length %d...",TotLen,m_NumSubSeqs,m_SubSeqLen);

// now allocate and initialise the Hamming edit distance array
// note that each processing thread will have it's own region of the array so that serialisation of r/w's is not required
m_AllocHamDist = ((INT64)(m_GenomeLen-2) * m_NumProcThreads * sizeof(UINT16)); // no hammings for start/end eBaseEOGs
#ifdef _WIN32
m_pHamDist = (UINT16 *) malloc((size_t)m_AllocHamDist);
if(m_pHamDist == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadGenome: Memory allocation of %lld bytes failed",(INT64)m_AllocHamDist);
	Reset(false);
	return(eBSFerrMem);
	}
#else
// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
m_pHamDist = (UINT16 *)mmap(NULL,(size_t)m_AllocHamDist, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
if(m_pHamDist == MAP_FAILED)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadGenome: Memory allocation of %lld bytes through mmap()  failed",(INT64)m_AllocHamDist,strerror(errno));
	m_pHamDist = NULL;
	Reset(false);
	return(eBSFerrMem);
	}
#endif

memset(m_pHamDist,m_SubSeqLen+1,(size_t)m_AllocHamDist);	// m_SubSeqLen+1 used as marker to show no hammings written to this loci

if((m_pThreadParams = new tsThreadParams [m_NumProcThreads]) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"unable to allocate memory (%d bytes) for thread parameter sets",sizeof(tsThreadParams) * m_NumProcThreads);
	Reset(false);
	return(eBSFerrMem);
	}
memset(m_pThreadParams,0,sizeof(tsThreadParams) * m_NumProcThreads);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Hamming distances initialised...");

return(eBSFSuccess);
}

// ThreadedIterParams
// use to return the processing parameters set for use calling thread instance
tsThreadParams *						// returns NULL if no more processing parameters available
ThreadedIterParams(void)
{
tsThreadParams *pThreadParam;
int Idx;
if((pThreadParam = m_pThreadParams) == NULL || m_NumProcThreads == 0)
	return(NULL);

#ifdef _WIN32
WaitForSingleObject(m_hMtxThreadParams,INFINITE);
#else
pthread_mutex_lock(&m_hMtxThreadParams);
#endif
for(Idx = 0; Idx < m_NumProcThreads; Idx++,pThreadParam++)
	{
	if(pThreadParam->State == 0)
		{
		pThreadParam->State = 1;
#ifdef _WIN32
		ReleaseMutex(m_hMtxThreadParams);
#else
		pthread_mutex_unlock(&m_hMtxThreadParams);
#endif
		return(pThreadParam);
		}
	}
#ifdef _WIN32
ReleaseMutex(m_hMtxThreadParams);
#else
pthread_mutex_unlock(&m_hMtxThreadParams);
#endif
return(NULL);
}
