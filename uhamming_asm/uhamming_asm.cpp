// uhamming.cpp : Defines the entry point for the console application.
// generates Hamming edit distances for all sequences of specified length over a target genome
// low level hamming generation is written in assembler on the postulate that this would reduce the computational cost
// as yet though, benchmarking is showing no real throughput advantage relative to the same functions written in 'C'
//
// EXPERIMENTAL!

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

const unsigned int cProgVer = 112;		// increment with each release

const int cMinSeqLen = 5;				// minimum sequence length for Hamming distances
const int cDfltSeqLen = 36;				// default sequence length for Hamming distances
const int cMaxSeqLen = 500;				// maximum sequence length for Hamming distances

const int cMaxWorkerThreads = 128;			// limiting max number of threads to this many
const int cMaxNumNodes = 10000;			// allow for upto this many nodes if processing is distributed over multiple nodes

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
	UINT32 HamDistOfs;  // offset into m_HamDist at which the edit distances for this processing instance are to be written
	int State;		    // processing state: 0 if waiting to be processed, 1 if currently being processed, 2 if processing completed
} tsThreadParams;

#pragma pack(1)
// Hamming specific structures
const int cMaxHammingChroms = 200;	// can handle at most this many chromosomes with hammings
typedef struct TAG_sHamChrom {
	UINT32 ChromID;					// uniquely identifies this chromosome
	UINT8  szChrom[cMaxDatasetSpeciesChrom];	// chrom name
	UINT32 NumEls;					// number of subsequences with hammings on this chrom
	UINT8 Dists[1];					// array, in ascending loci order, of hamming distances
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
	ePMdefault,					// default processing mode
	ePMdist,					// distributed processing over multiple nodes
	ePMmerge,					// merge multiple node Hammings
	ePMtrans,					// transform Hamming csv file into quick load binary
	ePMplaceholder				// used to set the enumeration range
	} etPMode;


int
Process(etPMode PMode,			// processing mode
		bool bWatsonOnly,		// true if watson strand only processing
		UINT32 SweepStart,		// start processing from this sweep instance (1..GenomeLen) or if distributed processing then the total number of nodes
		UINT32 SweepEnd,		// finish processing at this sweep instance (0 for all remaining or SweepStart..GenomeLen) or if distributed processing then the node instance
		int SeqLen,				// Hamming for these length sequences
		int SampleN,			// sample (process) every N sweep instance
		int NumThreads,			// process with upto this number of threads
		char *pszGenomeFile,	// bioseq file containing targeted genome assembly
		char *pszHammingFile);	// writeout Hamming edit distances into this file

int LoadGenome(char *pszBioSeqFile); // load genome from this file

extern void GHamDistWatson(UINT8 *pHDs,		// where to return Hamming differentials for each subsequence 
			  int SubSeqLen,	// generate Hammings edit distances for subsequences of this length
			  UINT32 SSofs,		// offset between subsequences for current pass
			  UINT8 *pGenomeSeq, // genome sequence (concatenated chrom seqs, separated by eBaseEOSs) with final chrom terminated by eBaseEOG, not eBaseEOS
			  UINT32 GenomeLen);	 // total genome length including chrom eBaseEOS markers, but exluding start/end eBaseEOG markers

extern void GHamDistCrick(UINT8 *pHDs,		// where to return Hamming differentials for each subsequence 
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

int NumberOfProcessors;		// number of installed Core CPUs
int NumThreads;				// number of threads (0 defaults to number of Core CPUs)
int SeqLen;					// Hamming edit distances for this length sequences
int SweepStart;			    // process starting from this sweep instance inclusive
int SweepEnd;				// complete processing at this sweep instance inclusive

int Node;					// node instance (1..N) if processing over multiple nodes
int NumNodes;				// total number of nodes (N) if processing over multiple nodes

int SampleN;				// sample every N sweep instances

bool bWatsonOnly;			// true if watson only strand processing - Crick is rather slow...

char szInFile[_MAX_PATH];	// input genome assembly file
char szOutFile[_MAX_PATH];	// generated Hamming distances file


// command line args
struct arg_lit  *help    = arg_lit0("hH","help",                "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_int *ScreenLogLevel=arg_int0("S", "ScreenLogLevel",	"<int>","Level of diagnostics written to screen 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"diagnostics log file");

struct arg_int *pmode = arg_int0("m","mode","<int>",		    "distributed processing mode: 0 - single node, 1 - multiple nodes, 2 - merge multiple Hamming files, 3 - transform Hamming CSV into quick load binary format  (default = 0)");
struct arg_lit  *crick = arg_lit0("c","crick",                  "process Crick in addition to Watson strand - Caution: very slow processing");


struct arg_int *numnodes = arg_int0("n","numnodes","<int>",	    "total number of nodes (2..N) if processing over multiple nodes");
struct arg_int *node = arg_int0("N","node","<int>",	            "node instance (1..N) if processing over multiple nodes");

struct arg_int *sweepstart = arg_int0("b","sweepstart","<int>",	"process starting from this sweep instance inclusive (default = 1 for 1st)");
struct arg_int *sweepend = arg_int0("B","sweepend","<int>",		"complete processing at this sweep instance inclusive (default = 0 for all remaining, or >= Sweep start)");
struct arg_int *seqlen = arg_int0("s","seqlen","<int>",			"Hamming edit distances for these length subsequences (default is 36)");
struct arg_file *infile = arg_file1("i","in","<file>",			"in mode 0 and 1, input from this bioseq genome assembly file or in mode 2 merge from this input Hamming file");

struct arg_int *sample = arg_int0("k","sample","<int>",		    "sample every -S<N> sweep instances (default is 1) useful if only interested in overall distributions");

struct arg_file *outfile = arg_file0("o","out","<file>",		"output (merged) Hamming distances to this file");
struct arg_int *threads = arg_int0("T","threads","<int>",		"number of processing threads 0..128 (defaults to 0 which sets threads to number of CPU cores)");
struct arg_end *end = arg_end(20);

void *argtable[] = {help,version,FileLogLevel,ScreenLogLevel,LogFile,
					pmode,crick,numnodes,node,sweepstart,sweepend,seqlen,sample,infile,outfile,threads,
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

	PMode = (etPMode)(pmode->count ? pmode->ival[0] : 0);
	if(PMode < 0 || PMode >= ePMplaceholder)
		{
		printf("\nError: Processing mode '-m%d' specified outside of range %d..%d",PMode,0,(int)ePMplaceholder-1);
		exit(1);
		}


	SeqLen = cDfltSeqLen;
	NumNodes = 1;
	SweepStart = 0;
	SweepEnd = 0;
	NumThreads = 1;
	SampleN = 1;
	bWatsonOnly = true;
	if(PMode == ePMdefault || PMode == ePMdist)
		{
#ifdef _WIN32
		SYSTEM_INFO SystemInfo;
		GetSystemInfo(&SystemInfo);
		NumberOfProcessors = SystemInfo.dwNumberOfProcessors;
#else
		NumberOfProcessors = sysconf(_SC_NPROCESSORS_CONF);
#endif
		if((NumThreads = threads->count ? threads->ival[0] : NumberOfProcessors)==0)
			NumThreads = NumberOfProcessors;
		if(NumThreads < 0 || NumThreads > NumberOfProcessors)
			{
			printf("\nError: Number of threads '-T%d' specified was outside of range %d..%d",PMode,0,NumberOfProcessors);
			exit(1);
			}

		if(PMode == ePMdefault)
			{
			SweepStart = sweepstart->count ? sweepstart->ival[0] : 1;
			if(SweepStart < 0)
				{
				printf("\nError: Sweep start '-b%d' must be >= 1",SweepStart);
				exit(1);
				}
			if(SweepStart == 0)	// allow a little latitude, treat 0 as being the same as 1...
				SweepStart = 1;
			SweepEnd = sweepend->count ? sweepend->ival[0] : 0;
			if(SweepEnd != 0 && SweepEnd < SweepStart)
				{
				printf("\nError: Sweep end '-B%d' must be either 0 or >= %d",SweepEnd,SweepStart);
				exit(1);
				}
			}
		if(PMode == ePMdist)
			{
			if(!numnodes->count || !node->count)
				{
				printf("\nError: In distributed processing mode both number of nodes '-n<num>' and node instance '-N<node>' must be specified");
				exit(1);
				}
			NumNodes = numnodes->ival[0];
			if(NumNodes < 2 || NumNodes > cMaxNumNodes)
				{
				printf("\nError: Number of nodes '-n%d' must be in the range 2..%d",NumNodes,cMaxNumNodes);
				exit(1);
				}

			Node = node->ival[0];
			if(Node < 1)
				{
				printf("\nError: node instance '-N%d' must be in range 1..%d",NumNodes);
				exit(1);
				}
			}

		SeqLen = seqlen->count ? seqlen->ival[0] : cDfltSeqLen;
		if(SeqLen < cMinSeqLen || SeqLen > cMaxSeqLen)
			{
			printf("\nError: Core length '-c%d' specified outside of range %d..%d",SeqLen,cMinSeqLen,cMaxSeqLen);
			exit(1);
			}

		SampleN = sample->count ? sample->ival[0] : 1;
		if(SampleN < 1 || SampleN > 100000000)
			{
			printf("\nError: Sampling '-S%d' specified outside of range 1..100000000",SampleN);
			exit(1);
			}
		if(crick->count > 0)
			bWatsonOnly = false;
		}

	strcpy(szInFile,infile->filename[0]);
	if(SampleN > 1)
		szOutFile[0] = '\0';
	else
		{
		if(outfile->count)
			{
			strncpy(szOutFile,outfile->filename[0],_MAX_PATH);
			szOutFile[_MAX_PATH-1] = '\0';
			}
		else
			{
			printf("\nError: No output results file specified with '-o<file' option");
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

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Version: %d.%2.2d Processing parameters:",cProgVer/100,cProgVer%100);

	const char *pszDescr;
	switch(PMode) {
		case ePMdefault:
			pszDescr = "Single node Hamming edit distance processing";
			break;
		case ePMdist:
			pszDescr = "Multiple nodes Hamming edit distance processing";
			break;
		case ePMmerge:
			pszDescr = "Merge two existing Hamming edit distance files";
			break;
		case ePMtrans:
			pszDescr = "Transform Hamming csv file into quick load binary file";
			break;
		default:
			pszDescr = "Non-standard processing";
			break;
		}
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Processing mode is : '%s'",pszDescr);
	switch(PMode) {
		case ePMdefault:
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
		default:
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Process %s",bWatsonOnly ? "Watson only strand" : "both Watson and Crick strands - Caution: very slow -");
			if(SampleN > 1)
				gDiagnostics.DiagOutMsgOnly(eDLInfo,"Only sample (process) every %d sweep",SampleN);
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Process subsequences of this length : %d",SeqLen);
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Input bioseq genome assembly file: '%s'",szInFile);
			if(SampleN == 1)
				gDiagnostics.DiagOutMsgOnly(eDLInfo,"Output Hamming distance file: '%s'",szOutFile);
			else
				gDiagnostics.DiagOutMsgOnly(eDLInfo,"Sampling so no output Hamming distance file will be generated");
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"number of threads : %d",NumThreads);
			break;
		}


#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	gStopWatch.Start();
	Rslt = Process(PMode,bWatsonOnly,(UINT32)SweepStart,(UINT32)SweepEnd,SeqLen,SampleN,NumThreads,szInFile,szOutFile);
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


typedef struct TAG_sGChrom {
	int ChromSeqID;			// m_pChromSeqs[ChromSeqID-1] of this tsChromSeq instance
	int ChromID;
	char szChromName[cMaxDatasetSpeciesChrom];
	int Len;				// actual length
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
int m_SubSeqLen;				// subsequence lengths over which to calc Hamming edit distances
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
UINT8 *m_pHamDist;				// allocated to hold Hamming edit distances (expected to be m_NumProcThreads * (m_GenomeLen-2)), start/final eBaseEOG have no hammings
INT64 m_AllocHamDist;			// allocation size for m_pHamDist
UINT32 m_NumSubSeqs;			// total number of actual Hamming subsequences in genome

tsHamHdr *m_pHamHdr;			// header for binary format hamming edit distances 
tsHamChrom *m_pCurHamChrom;		// pts to current chromosome specific binary hammings

tsThreadParams *m_pThreadParams;	// holds initialised phase parameter sets for each thread or processing core

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
Reset(void)
{
if(m_hOutFile != -1)
	{
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
m_AllocGenomeSeq = 0;
m_AllocHamDist = 0;
m_pCurHamChrom = NULL;
m_SubSeqLen=0;					
m_NumProcThreads=0;			
m_szSpecies[0] = '\0';
m_GenomeLen = 0;
m_NumSubSeqs = 0;
m_NumGChroms= 0;
m_hOutFile = -1;
}

void
Init(void)
{
m_pBioSeqFile = NULL;
m_pGChroms = NULL;
m_pCurHamChrom = NULL;
m_pHamHdr = NULL;

m_hOutFile = -1;
m_pGenomeSeq = NULL;		// allocated to hold concatenated chromosome sequences with eBaseEOS/eBaseEOG markers
m_pHamDist = NULL;			// allocated to hold Hamming edit distances
m_pThreadParams = NULL;	    // holds initialised thread parameter sets
Reset();
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
while((pParams = ThreadedIterParams())!=NULL)
	{
	for(Seq = pParams->SSofs, Idx = 0; Idx <= pParams->NumSeqs; Idx++, Seq+=pParams->SeqDelta)
		{
		if((Seq + m_SubSeqLen) >= m_GenomeLen)
			break;
		if(!(Idx % 1000000))
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"Thread %d Processing subsequence %u",pParams->ThreadID,Seq);
		if(pParams->SampleN > 1)	// if sampling....
			{
			if(Idx > 0 && (Idx % pParams->SampleN))
				continue;
			}
		// no hammings to be generated for inital/final eBaseEOG markers
		if((Seq + m_SubSeqLen) < m_GenomeLen-1)
			GHamDistWatson(&m_pHamDist[pParams->HamDistOfs], m_SubSeqLen,Seq,&m_pGenomeSeq[1],m_GenomeLen-2);
		if(!pParams->bWatsonOnly && (Seq + m_SubSeqLen) <= (m_GenomeLen-1))
			GHamDistCrick(&m_pHamDist[pParams->HamDistOfs], m_SubSeqLen,Seq-1,&m_pGenomeSeq[1],m_GenomeLen-2);
		}
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Thread %d Processing completed",pParams->ThreadID);
	}

#ifdef _WIN32
_endthreadex(0);
return(eBSFSuccess);
#else
pthread_exit(NULL);
#endif
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
				CurHamming++;
					break;
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
	Reset();
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
		// because of a bug in the first version of this application then a non-matching last record is simply sloughed
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
LoadHammings(char *pszHammings)
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
if((Rslt=LoadHammings(pszInFile)) < eBSFSuccess)
	return(Rslt);

// output as binary
#ifdef _WIN32
m_hOutFile = open(pszOutFile,O_CREATETRUNC );
#else
if((m_hOutFile = open(pszOutFile,O_RDWR | O_CREAT,S_IREAD | S_IWRITE))!=-1)
    if(ftruncate(m_hOutFile,0)!=0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate %s - %s",pszOutFile,strerror(errno));
			return(eBSFerrCreateFile);
			}
#endif

if(m_hOutFile < 0)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate output file '%s'",pszOutFile);
	Reset();
	return(eBSFerrCreateFile);
	}

CUtility::SafeWrite(m_hOutFile,m_pHamHdr,m_pHamHdr->Len);

// that's it!
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

int
Process(etPMode PMode,			// processing mode
		bool bWatsonOnly,		// true if watson strand only processing
		UINT32 SweepStart,		// start processing from this sweep instance (1..GenomeLen) or if distributed processing then the total number of nodes
		UINT32 SweepEnd,		// finish processing at this sweep instance (0 for all remaining or SweepStart..GenomeLen) or if distributed processing then the node instance
		int SeqLen,				// Hamming for this length sequences
		int SampleN,			// sample (process) every N sweep instance
		int NumThreads,			// process with upto this number of threads
		char *pszGenomeFile,	// bioseq file containing targeted genome assembly
		char *pszHammingFile)	// write Hamming edit distances into this file
{
int Rslt;
UINT32 SSeqStart;
UINT32 SSeqEnd;
tsThreadParams *pThreadParam;
UINT32 BeginSeq;
UINT32 EndSeq;
UINT32 HamDistOfs;
UINT32 NumSubSeqs;
int ThreadIdx;
char szBuffer[16000];
int BuffOfs;
UINT32 CurLoci;
tsGChrom *pChrom;
int Idx;
UINT32 SeqIdx;
UINT8 *pHamDist1;
UINT8 *pHamDist2;

Init();

m_PMode = PMode;
m_SubSeqLen = SeqLen;
m_NumProcThreads = NumThreads;

if(PMode == ePMmerge)
	return(MergeHammings(pszGenomeFile,pszHammingFile));

if(PMode == ePMtrans)
	return(TransHammings(pszGenomeFile,pszHammingFile));

if((Rslt=LoadGenome(pszGenomeFile))!=eBSFSuccess)
	{
	Reset();
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
	}
else	// single node processing
	{
	// silently clamp SweepStart and SweepEnd to be no more than the actual genome length
	// SweepStart must be <= to the genome length and SweepEnd must be 0 or <= genome length
	if(SweepStart > m_GenomeLen)
		SweepStart = m_GenomeLen;
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

// adjust SweepStart (down) and SweepEnd (up) to allow for the presence of interchromosome EOS markers and rounding errors
if(SSeqStart > (UINT32)(2 * m_NumGChroms))
	SSeqStart -= (2 * m_NumGChroms);
else
	SSeqStart = 0;
SSeqEnd += (2 * m_NumGChroms);
if(SSeqEnd > m_GenomeLen)
	SSeqEnd = m_GenomeLen;

// log the actual sweep start/end
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Node sweep start is %u, and sweep end is %u",SSeqStart+1,SSeqEnd);

// adjust m_NumProcThreads so as to ensure each processing thread has at least 1 sweep to process
m_NumProcThreads = (int)min(1 + SSeqEnd - SSeqStart,(UINT32)m_NumProcThreads);

if(SampleN == 1)
{
// ensure output results file can be accessed/created
#ifdef _WIN32
	m_hOutFile = open(pszHammingFile,O_CREATETRUNC );
#else
	if((m_hOutFile = open(pszHammingFile,O_RDWR | O_CREAT,S_IREAD | S_IWRITE))!=-1)
		if(ftruncate(m_hOutFile,0)!=0)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate %s - %s",pszHammingFile,strerror(errno));
				Reset();
				return(eBSFerrCreateFile);
				}
#endif

	if(m_hOutFile < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate output file '%s'",pszHammingFile);
		Reset();
		return(eBSFerrCreateFile);
		}
	close(m_hOutFile);
	m_hOutFile = -1;
}



HamDistOfs = 0;
pThreadParam = m_pThreadParams;

BeginSeq = SSeqStart;
EndSeq = SSeqEnd;
for(ThreadIdx = 0; ThreadIdx < m_NumProcThreads; ThreadIdx++,pThreadParam++)
	{
	pThreadParam->ThreadID = ThreadIdx+1;
	pThreadParam->bWatsonOnly = bWatsonOnly;
	NumSubSeqs = (EndSeq - SSeqStart)/(m_NumProcThreads - ThreadIdx);
	pThreadParam->SSofs = ++BeginSeq;
	pThreadParam->SeqDelta = m_NumProcThreads;
	pThreadParam->NumSeqs = NumSubSeqs;
	pThreadParam->HamDistOfs = HamDistOfs;
	pThreadParam->SampleN = SampleN;
	EndSeq -= NumSubSeqs;
	HamDistOfs += (m_GenomeLen-2);		// initial/final eBaseEOG do not have hammings
	pThreadParam->State = 0;
	}

#ifdef _WIN32
if((m_hMtxThreadParams = CreateMutex(NULL,false,NULL))==NULL)
#else
if(pthread_mutex_init (&m_hMtxThreadParams,NULL)!=0)
#endif
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: unable to create mutex");
	Reset();
	return(eBSFerrInternal);
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

// wait for all threads to have completed
for(ThreadIdx = 0; ThreadIdx < m_NumProcThreads; ThreadIdx++)
	{
#ifdef _WIN32
	WaitForSingleObject( WorkerThreads[ThreadIdx].threadHandle, INFINITE );
	CloseHandle( WorkerThreads[ThreadIdx].threadHandle);
#else
	pthread_join(WorkerThreads[ThreadIdx].threadID,NULL);
#endif
	}

#ifdef _WIN32
CloseHandle(m_hMtxThreadParams);
#else
pthread_mutex_destroy(&m_hMtxThreadParams);
#endif

// now merge the Hammings
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Merging Hamming edit distances...");

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

if(SampleN == 1)
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
				return(eBSFerrCreateFile);
				}
#endif

	if(m_hOutFile < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate output file '%s'",pszHammingFile);
		Reset();
		return(eBSFerrCreateFile);
		}

	pChrom = &m_pGChroms[0];
	pHamDist1 = m_pHamDist;
	BuffOfs = 0;
	CurLoci = 0;
	BuffOfs = sprintf(szBuffer,"%lld,%d,%d\n",m_GenomeLen,SSeqStart+1,SSeqEnd);
	for(SeqIdx = 0; SeqIdx < (m_GenomeLen-2); SeqIdx++,pHamDist1++,CurLoci++)
		{
		if(CurLoci >= pChrom->NumSubSeqs)
			{
			if(pChrom->ChromSeqID == m_NumGChroms)
				break;
			pChrom += 1;
			CurLoci = 0;
			SeqIdx += (2*m_SubSeqLen)-1;
			pHamDist1 += (2*m_SubSeqLen)-1;
			}
		BuffOfs += sprintf(&szBuffer[BuffOfs],"\"%s\",%d,%d\n",pChrom->szChromName,CurLoci,*pHamDist1);
		if((BuffOfs + 200) > sizeof(szBuffer))
			{
			CUtility::SafeWrite(m_hOutFile,szBuffer,BuffOfs);
			BuffOfs = 0;
			}
		}
	if(BuffOfs)
		{
		CUtility::SafeWrite(m_hOutFile,szBuffer,BuffOfs);
		BuffOfs = 0;
		}
	close(m_hOutFile);
	m_hOutFile = -1;
	}

// log basic count distribution.....
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Distribution:\nEditDist,Freq,Proportion");
int CntHist[256];
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
		SeqIdx += (2*m_SubSeqLen)-1;
		pHamDist1 += (2*m_SubSeqLen)-1;
		}
	CntHist[*pHamDist1] += 1;
	}

for(SeqIdx = 0; SeqIdx < (UINT32)min(66,m_SubSeqLen); SeqIdx++)
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"%d,%d,%1.3f",SeqIdx,CntHist[SeqIdx],(CntHist[SeqIdx]*100.0f)/m_NumSubSeqs);

Reset();
return(eBSFSuccess);
}

int
MinHamCnt(int MaxHamCnt,			// only cnt upto this many mismatches
		   etSeqBase *pSeq1,		// determine Hamming edit distance between this sequence
		   etSeqBase *pSeq2)		// and this sequence
{
int Idx;
int CurCnt=0;
for(Idx = 0; Idx < m_SubSeqLen; Idx++)
	if(*pSeq1++ != *pSeq2++ && CurCnt++ >= MaxHamCnt)
		return(CurCnt);
return(CurCnt); 
}

int
LoadGenome(char *pszBioSeqFile) // load genome from this file
{
int Rslt;
int Len;
INT64 TotLen;
int ChromID;
int ChromIdx;
UINT32 NumSubSeqs;
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
	Reset();
	return(Rslt);
	}

m_pBioSeqFile->GetTitle(sizeof(m_szSpecies),m_szSpecies);

m_NumGChroms = m_pBioSeqFile->NumEntries();
if((m_pGChroms = new tsGChrom [m_NumGChroms])==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"unable to allocate memory (%d bytes) for %d ChromSeqs",sizeof(tsGChrom) * m_NumGChroms,m_NumGChroms);
	Reset();
	return(eBSFerrMem);
	}
memset(m_pGChroms,0,sizeof(tsGChrom) * m_NumGChroms);

// determine total sequence length required for allocating to hold genome as one concatenated sequence
// concatenated sequence starts with eBaseEOG followed be each chromosome sequence separated from the 
// next chromosome sequence by m_SubSeqLen eBaseEOSs, and finally the last chromosome is terminated by eBaseEOG, not eBaseEOS
m_GenomeLen = 1;		// allow for initial eBaseEOG marker
ChromID = 0;
ChromIdx = 0;
while((ChromID = m_pBioSeqFile->Next(ChromID)) > 0) {
	pChrom = &m_pGChroms[ChromID-1];
	pChrom->ChromID = ChromID;
	pChrom->SeqOfs = m_GenomeLen;
	m_pBioSeqFile->GetName(ChromID,sizeof(pChrom->szChromName),pChrom->szChromName);
	Len = m_pBioSeqFile->GetDataLen(ChromID);
	m_GenomeLen += Len + m_SubSeqLen;		// allow for a concatenation separators (eBaseEOS)
	pChrom->Len = Len;
	pChrom->pSeq = NULL;
	pChrom->ChromSeqID = ++ChromIdx;
	}
m_GenomeLen -= (m_SubSeqLen - 1);
m_AllocGenomeSeq = m_GenomeLen + sizeof(UINT64); // allocate for a little more so later could perhaps use int64 loads without memeory access errors 
#ifdef _WIN32
// use malloc instead of new() because of gnu new()/malloc issues - helps portability
m_pGenomeSeq = (UINT8 *) malloc((size_t)m_AllocGenomeSeq);	
if(m_pGenomeSeq == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadGenome: Memory allocation of %lld bytes failed",(INT64)m_AllocGenomeSeq);
	Reset();
	return(eBSFerrMem);
	}
#else
// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations so use mmap
m_pGenomeSeq = (UINT8 *)mmap(NULL,(size_t)m_AllocGenomeSeq, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
if(m_pGenomeSeq == MAP_FAILED)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadGenome: Memory allocation of %lld bytes through mmap()  failed",(INT64)m_AllocGenomeSeq,strerror(errno));
	m_pGenomeSeq = NULL;
	Reset();
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
	if(--ChromIdx)
		{
		memset(pSeq,eBaseEOS,m_SubSeqLen); // mark end of current chrom
		pSeq += m_SubSeqLen;
		}
	else
		*pSeq = eBaseEOS;
	while(*pGseq++ != eBaseEOS)		// processing expects repeat masking removed
		pGseq[-1] &= ~ cRptMskFlg;
	}
*pSeq = eBaseEOG;   // marks prev chromosome sequence as being the last
delete m_pBioSeqFile;
m_pBioSeqFile = NULL;

// determine the total number of subsequences over which Hammings are required
m_NumSubSeqs = 0;
pChrom = m_pGChroms;
for(ChromIdx = 0; ChromIdx < m_NumGChroms; ChromIdx++, pChrom++)
	{
	if((NumSubSeqs = (UINT32)(pChrom->Len + 1 - m_SubSeqLen)) < 1) // don't bother with chromosome too short to contain a subseq
		{
		pChrom->NumSubSeqs = 0;
		continue;
		}
	pChrom->NumSubSeqs = NumSubSeqs;
	m_NumSubSeqs += NumSubSeqs;
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Genome containing %llu total nucleotides loaded with %llu subsequences of length %d...",TotLen,(size_t)m_NumSubSeqs,m_SubSeqLen);

// now allocate and initialise the Hamming edit distance array
// note that each processing thread will have it's own region of the array so that serialisation of r/w's is not required
m_AllocHamDist = ((INT64)(m_GenomeLen-2) * m_NumProcThreads); // no hammings for start/end eBaseEOGs
#ifdef _WIN32
m_pHamDist = (UINT8 *) malloc((size_t)m_AllocHamDist);	
if(m_pHamDist == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadGenome: Memory allocation of %lld bytes failed",(INT64)m_AllocHamDist);
	Reset();
	return(eBSFerrMem);
	}
#else
// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
m_pHamDist = (UINT8 *)mmap(NULL,(size_t)m_AllocHamDist, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
if(m_pHamDist == MAP_FAILED)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadGenome: Memory allocation of %lld bytes through mmap()  failed",(INT64)m_AllocHamDist,strerror(errno));
	m_pHamDist = NULL;
	Reset();
	return(eBSFerrMem);
	}
#endif

memset(m_pHamDist,m_SubSeqLen+1,(size_t)m_AllocHamDist);	// m_SubSeqLen+1 used as marker to show no hammings written to this loci

if((m_pThreadParams = new tsThreadParams [m_NumProcThreads]) == NULL)	 
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"unable to allocate memory (%d bytes) for thread parameter sets",sizeof(tsThreadParams) * m_NumProcThreads);
	Reset();
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
