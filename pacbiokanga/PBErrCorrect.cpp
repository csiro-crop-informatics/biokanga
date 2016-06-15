// Copyright 2013 CSIRO  ( http://www.csiro.au/ ) 
//   Licensed under the Apache License, Version 2.0 (the "License");
//   you may not use this file except in compliance with the License.
//   You may obtain a copy of the License at
//       http://www.apache.org/licenses/LICENSE-2.0
//
//   Unless required by applicable law or agreed to in writing, software
//   distributed under the License is distributed on an "AS IS" BASIS,
//   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//   See the License for the specific language governing permissions and
//   limitations under the License
//   Please contact stuart.stephen@csiro.au for support or 
//   to submit modifications to this source
#include "stdafx.h"
// Supporting at most this many concurrent TCP sessions between service requester (this server) and all service providers
// this could be increased to an internally restricted maximum of 511 but there could then be a significant performance throughput degradation
// because of the use of select() instead of ePoll or derivatives.
// It is essential that there is consistency between the various source files in the cMaxConcurrentSessions value
#ifdef WIN32
#define cMaxConcurrentSessions 100                     // limit to less than 511
#else
#define cMaxConcurrentSessions 100
#endif

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#if _WIN32
#include <process.h>
#include "../libbiokanga/commhdrs.h"
// default on Windows is for a maximum of 64 sockets in a FD_SET select() descriptor set
// increase this limit to the maximum number of supported TCP sessions for all service provider sessions
// plus 1 for the listening socket and 1 for the control socket

#if (cMaxConcurrentSessions + 2 > 64)
#define FD_SETSIZE  (cMaxConcurrentSessions + 2)
#endif
#include <winsock2.h>
#include <ws2tcpip.h>
#else
#include <sys/time.h>
#include <sys/resource.h>
#include <sys/mman.h>
#include <pthread.h>
#include <sys/ioctl.h>
#include <sys/socket.h>
#include <netdb.h>
typedef struct sockaddr_storage SOCKADDR_STORAGE;
#include "../libbiokanga/commhdrs.h"
#endif

#include "pacbiokanga.h"

#include "../libbiokanga/bgzf.h"
#include "SSW.h"
#include "pacbiocommon.h"
#include "SeqStore.h"
#include "AssembGraph.h"
#include "./BKScommon.h"
#include "./BKSRequester.h"
#include "PBErrCorrect.h"

int
ProcPacBioErrCorrect(etPBPMode PMode,		// processing mode
		char *pszHostName,			// listening on this host name or IPv4/IPv5 address for connections by service providers 
		char *pszServiceName,		// Listen on this service name or port for for connections by service providers
		int MaxRMI,					// max number of RMI service provider instances supported
		int MaxNonRMI,				// max number of non-RMI SW threads
		int SampleInRate,			// sample input sequences at this rate per 100
		int SampleAcceptRate,		// sample accepted input sequences which are at least MinPBSeqLen bp at this rate per 1000 accepted
		int FiltMinHomoLen,			// filter PacBio reads for homopolymer runs >= this length (0 to disable filtering) 
		bool bSenseOnlyOvlps,		// process for sense only overlaps (default is for sense/sense and sense/antisense overlaps)
		int DeltaCoreOfs,			// offset by this many bp the core windows of coreSeqLen along the probe sequence when checking for overlaps
		int MaxSeedCoreDepth,		// only further process a seed core if there are no more than this number of matching cores in all targeted sequences
		int MinSeedCoreLen,			// use seed cores of this length when identifying putative overlapping scaffold sequences
		int MinNumSeedCores,        // require at least this many seed cores between overlapping scaffold sequences
		int SWMatchScore,			// score for matching bases (0..50)
		int SWMismatchPenalty,		// mismatch penalty (-50..0)
		int SWGapOpenPenalty,		// gap opening penalty (-50..0)
		int SWGapExtnPenalty,		// gap extension penalty (-50..0)
		int SWProgExtnPenaltyLen,	// progressive gap scoring then only apply gap extension score if gap at least this length (0..63) - use if aligning PacBio
		int TransciptomeLens,		// 0 if disabled, processing transcript reads, putatively overlapping reads must have length differential no more than this percentage and overlaps to be nearly full length
		int MinPBSeqLen,			// only accepting PacBio reads of at least this length (defaults to 10Kbp) and if
		int MaxPBSeqLen,			// no more than this length (defaults to 35Kbp)
		int MinPBSeqOverlap,		// any overlap of a PacBio onto a target PacBio must be of at least this many bp to be considered for contributing towards error correction (defaults to 5Kbp) 
		int MaxArtefactDev,			// classify overlaps as artefactual if sliding window of 500bp over any overlap deviates by more than this percentage from the overlap mean
		int MinHCSeqLen,			// only accepting hiconfidence reads of at least this length (defaults to 1Kbp)
		int MinHCSeqOverlap,		// any overlap of a hiconfidence read onto a target PacBio read must be of at least this many bp to be considered for contributing towards error correction (defaults to 1Kbp)
		int HCRelWeighting,         // hiconfidence read overlaps are usually weighted higher than normal lesser confidence read overlaps when calling consensus bases 
		int MinErrCorrectLen,		// error corrected sequences must be at least this minimum length
		int MinConcScore,			// error corrected sequences trimmed until mean 50bp concensus score is at least this threshold
		int NumPacBioFiles,			// number of input pacbio file specs
		char *pszPacBioFiles[],		// input pacbio files
		int NumHiConfFiles,			// number of input hiconfidence file specs
		char *pszHiConfFiles[],		// input hiconfidence files		
	    char *pszErrCorFile,		// name of file into which write error corrected sequences
		char *pszMultiAlignFile,	// name of file into which write multiple alignments
		char *pszChkPtsFile,        // name of file used for checkpointing in case resume processing is required
		int NumThreads);			// maximum number of worker threads to use


#ifdef _WIN32
int ProcErrCorrect(int argc, char* argv[])
{
// determine my process name
_splitpath(argv[0],NULL,NULL,gszProcName,NULL);
#else
int
ProcErrCorrect(int argc, char** argv)
{
// determine my process name
CUtility::splitpath((char *)argv[0],NULL,gszProcName);
#endif

int iFileLogLevel;			// level of file diagnostics
int iScreenLogLevel;		// level of file diagnostics
char szLogFile[_MAX_PATH];	// write diagnostics to this file
int Idx;                    // file index 
int Rslt = 0;   			// function result code >= 0 represents success, < 0 on failure

int PMode;					// processing mode
bool bSenseOnlyOvlps;		// process for sense only overlaps (default is for sense/sense and sense/antisense overlaps)
int FiltMinHomoLen;			// filter PacBio reads for homopolymer runs >= this length (0 to disable filtering) 
int MinSeedCoreLen;			// use seed cores of this length when identifying putative overlapping sequences
int MinNumSeedCores;        // require at least this many seed cores between overlapping sequences before attempting SW
int DeltaCoreOfs;			// offset by this many bp the core windows of coreSeqLen along the probe sequence when checking for overlaps
int MaxSeedCoreDepth;		// only further process a seed core if there are no more than this number of matching cores in all targeted sequences

int SWMatchScore;			// score for matching bases (0..50)
int SWMismatchPenalty;		// mismatch penalty (0..50)
int SWGapOpenPenalty;		// gap opening penalty (0..50)
int SWGapExtnPenalty;		// gap extension penalty (0..50)
int SWProgExtnPenaltyLen;	// progressive gap scoring then only apply gap extension score if gap at least this length (0..63) - use if aligning PacBio

int TranscriptomeLens;		// transcriptome assembly - overlap read lengths must match within this percentage and overlaps be full length

int MinPBSeqLen;			// only accepting PacBio reads of at least this length (defaults to 15Kbp) and if
int MaxPBSeqLen;			// no more than this length (defaults to 35Kbp)
int MinPBSeqOverlap;		// any overlap of a PacBio onto a target PacBio must be of at least this many bp to be considered for contributing towards error correction (defaults to 5Kbp) 
int MinHCSeqLen;			// only accepting hiconfidence reads of at least this length (defaults to 1Kbp)
int MinHCSeqOverlap;		// any overlap of a hiconfidence read onto a target PacBio read must be of at least this many bp to be considered for contributing towards error correction (defaults to 1Kbp) 
int HCRelWeighting;         // hiconfidence read overlaps are usually weighted higher than normal lesser confidence read overlaps when calling consensus bases 
int MinConcScore;			// error corrected sequences trimmed until mean 50bp concensus score is at least this threshold
int MinErrCorrectLen;		// error corrected and sequence trimmed sequences must be of at least this length
int MaxArtefactDev;			// classify overlaps as artefactual if sliding window of 500bp over any overlap deviates by more than this percentage from the overlap mean

CSimpleGlob glob(SG_GLOB_FULLSORT);
char szTempFilePath[_MAX_PATH+1];
char *pszInfile;
int FileID;
int NumPacBioFiles;						// number of input pacbio file specs
char *pszPacBioFiles[cMaxInFileSpecs+1];	// input pacbio files
int NumHiConfFiles;						// number of input hiconfidence file specs
char *pszHiConfFiles[cMaxInFileSpecs+1];	// input hiconfidence files

char szOutFile[_MAX_PATH+1];				// where to write (ePBPMConsensus or ePBPMErrCorrect) error corrected sequences or (ePBPMOverlapDetail) overlap detail
char szOutMAFile[_MAX_PATH+1];				// where to write optional multialignments
char szChkPtsFile[_MAX_PATH+1];				// used for error correcting checkpointing

int SampleInRate;							// sample input sequences at this rate per 100
int SampleAcceptRate;						// sample accepted input sequences at this rate per 1000

int NumberOfProcessors;						// number of installed CPUs
int NumThreads;								// number of threads (0 defaults to number of CPUs)

int MaxNonRMI;								// max number of non-RMI SW threads
int MaxRMI;									// max number of RMI service provider instances supported
char szHostName[cMaxHostNameLen];			// listening on this host name or IPv4/IPv5 address for connections by service providers 
char szServiceName[cMaxServiceNameLen];		// Listen on this service name or port for for connections by service providers

char szSQLiteDatabase[_MAX_PATH];	// results summaries to this SQLite file
char szExperimentName[cMaxDatasetSpeciesChrom+1];			// experiment name
char szExperimentDescr[1000];		// describes experiment

//
struct arg_lit  *help    = arg_lit0("h","help",                 "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to screen and logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"diagnostics log file");

struct arg_int *pmode = arg_int0("m","pmode","<int>",					"processing mode - 0 error correct, 1 consensus sequence only, 2 scaffold overlap detail only, 3 consolidate transcripts into representative");

struct arg_lit  *senseonlyovlps    = arg_lit0("a","senseonlyovlps",    "process for sense only overlaps (default is for sense and antisense overlaps)");
struct arg_int  *transcriptomelens = arg_int0("t","transcriptome","<int>",	   "transcriptome assembly - overlap read lengths must match within this percentage and overlaps be full length");


struct arg_int *minpbseqlen = arg_int0("l","minpbseqlen","<int>",		"minimum individual PacBio sequence length to error correct (default 10000, range 500 to 100000)");
struct arg_int *maxpbseqlen = arg_int0("L", "maxpbseqlen", "<int>",		"maximum individual PacBio sequence length (default 35000, minimum minpbseqlen)");

struct arg_int *minpbseqovl = arg_int0("b","minpbseqovl","<int>",		"minimum PacBio overlap onto PacBio length required (default 5000, range 500 to 100000)");

struct arg_int *minhcseqlen = arg_int0("p","minhcseqlen","<int>",		"minimum individual high confidence sequence length (default 1000, range 250 to 100000)");
struct arg_int *minhcseqovl = arg_int0("P","minhcseqovl","<int>",		"minimum high confidence sequence overlap onto PacBio length required (default 500, range 250 to 100000)");
struct arg_int *hcrelweighting = arg_int0("r","hcrelweighting","<int>",	"high confidence sequence relative weighting when consensus base calling (default 3, range 1 to 10)");

struct arg_int *minfilthomolen = arg_int0("H","minfilthomolen","<int>",			"filtering for near homopolymer runs of at least this length (default 16, 0 to disable, accepted range 12 to 25)");


struct arg_int *minseedcorelen = arg_int0("c","seedcorelen","<int>",			"use seed cores of this length when identifying putative overlapping sequences (default 14, range 12 to 50)");
struct arg_int *minseedcores = arg_int0("C","minseedcores","<int>",				"require at least this many accepted seed cores between overlapping sequences to use SW (default 10, range 1 to 50)");

struct arg_int *deltacoreofs = arg_int0("d","deltacoreofs","<int>",				"offset cores (default 2, range 1 to 10)");
struct arg_int *maxcoredepth = arg_int0("D","maxcoredepth","<int>",				"explore cores of less than this maximum depth (default 15000, range 1000 to 20000)");

struct arg_int *matchscore = arg_int0("x","matchscore","<int>",					"SW score for matching bases (default 3, range 1 to 50)");
struct arg_int *mismatchpenalty = arg_int0("X","mismatchpenalty","<int>",		"SW mismatch penalty (default 7, range 1 to 50)");
struct arg_int *gapopenpenalty = arg_int0("y","gapopenpenalty","<int>",			"SW gap opening penalty (default 4, range 1 to 50)");
struct arg_int *gapextnpenalty = arg_int0("Y","gapextnpenalty","<int>",			"SW gap extension penalty (default 1, range 1 to 50)");
struct arg_int *progextnpenaltylen = arg_int0("z","progextnpenaltylen","<int>",	"SW gap extension penalty only applied for gaps of at least this number of bases (default 2, range 1 to 63)");

struct arg_int *minconcscore = arg_int0("s","minconcscore","<int>",			     "error corrected sequences trimmed until mean 50bp concensus score is at least this threshold (default 3, range 0 to 9)");
struct arg_int *minerrcorrectlen = arg_int0("S","minerrcorrectlen","<int>",		 "error corrected and trimmed sequences must be at least this minimum length (default 5000, range 500 to 20000)");
struct arg_int *maxartefactdev = arg_int0("A","artefactdev","<int>",			 "classify overlaps as artefactual if 500bp window score deviates by more than this percentage from complete overlap mean (0 to disable, range 1 to 70)");

struct arg_file *hiconffiles = arg_filen("I","hiconffile","<file>",0,cMaxInFileSpecs,		"optional, names of input files containing higher confidence reads or sequences to be used in error correcton of PacBio reads (wildcards allowed)");
struct arg_file *pacbiofiles = arg_filen("i","pacbiofile","<file>",1,cMaxInFileSpecs,		"names of input files containing PacBio sequences to be error corrected (wildcards allowed)");
struct arg_file *outfile = arg_file1("o","out","<file>",									"output error corrected PacBio reads to this file");
struct arg_file *mafile = arg_file0("O","mafile","<file>",						"optional, output multialignments to this file, caution can grow very large");
struct arg_file *scaffovrlapsfile = arg_file0("e","scaffovrlapsfile","<file>",	"optional, output scaffolding overlap detail to this file");

struct arg_int *sampleinrate = arg_int0("R","sampleinrate","<int>",					"accept input sequences at this rate per 100 (default 100, range 1 to 100)");
struct arg_int *sampleacceptrate = arg_int0("Z","sampleacceptrate","<int>",			"sample accepted input sequences at this rate per 1000 (default 1000, range 1 to 1000)");

struct arg_int *threads = arg_int0("T","threads","<int>",						"number of processing threads 0..128 (defaults to 0 which sets threads to number of CPU cores)");

struct arg_int *maxnonrmi = arg_int0("N","maxnonrmi","<int>",					"if RMI SW processing then limit non-RMI to this maximum number of SW threads (defaults to threads, range 0 ... threads)");
struct arg_int *maxrmi = arg_int0("n","maxrmi","<int>",							"maximum number of RMI SW instances supported (defaults to 0 which sets max RMI instances to 4 x number of threads, range 16 to 500)");
struct arg_str  *rmihost = arg_str0("u", "rmihost", "<string>",					"listening on this host name or IPv4/IPv5 address for connections by SW service providers (default 127.0.0.1)");
struct arg_str  *rmiservice = arg_str0("U", "rmiservice", "<string>",			"Listen on this service name or port for connections by SW service providers (default 43123)");

struct arg_file *summrslts = arg_file0("q","sumrslts","<file>",				"Output results summary to this SQLite3 database file");
struct arg_str *experimentname = arg_str0("w","experimentname","<str>",		"experiment name SQLite3 database file");
struct arg_str *experimentdescr = arg_str0("W","experimentdescr","<str>",	"experiment description SQLite3 database file");

struct arg_end *end = arg_end(200);

void *argtable[] = {help,version,FileLogLevel,LogFile,
					pmode,rmihost,rmiservice,maxnonrmi,maxrmi,minfilthomolen,senseonlyovlps,minseedcorelen,minseedcores,deltacoreofs,maxcoredepth,
					matchscore,mismatchpenalty,gapopenpenalty,gapextnpenalty,progextnpenaltylen,
					transcriptomelens,minpbseqlen,maxpbseqlen,minpbseqovl,minhcseqlen,minhcseqovl,hcrelweighting,minconcscore,minerrcorrectlen,maxartefactdev,
					summrslts,pacbiofiles,hiconffiles,experimentname,experimentdescr,
					outfile,mafile,scaffovrlapsfile,sampleinrate,sampleacceptrate,threads,
					end};

char **pAllArgs;
int argerrors;
argerrors = CUtility::arg_parsefromfile(argc,(char **)argv,&pAllArgs);
if(argerrors >= 0)
	argerrors = arg_parse(argerrors,pAllArgs,argtable);

/* special case: '--help' takes precedence over error reporting */
if (help->count > 0)
        {
		printf("\n%s %s %s, Version %s\nOptions ---\n", gszProcName,gpszSubProcess->pszName,gpszSubProcess->pszFullDescr,cpszProgVer);
        arg_print_syntax(stdout,argtable,"\n");
        arg_print_glossary(stdout,argtable,"  %-25s %s\n");
		printf("\nNote: Parameters can be entered into a parameter file, one parameter per line.");
		printf("\n      To invoke this parameter file then precede it's name with '@'");
		printf("\n      e.g. %s %s @myparams.txt\n",gszProcName,gpszSubProcess->pszName);
		printf("\nPlease report any issues regarding usage of %s to stuart.stephen@csiro.au\n\n",gszProcName);
		return(1);
        }

    /* special case: '--version' takes precedence error reporting */
if (version->count > 0)
        {
		printf("\n%s %s Version %s\n",gszProcName,gpszSubProcess->pszName,cpszProgVer);
		return(1);
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

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Subprocess %s Version %s starting",gpszSubProcess->pszName,cpszProgVer);
	gExperimentID = 0;
	gProcessID = 0;
	gProcessingID = 0;
	szSQLiteDatabase[0] = '\0';
	szExperimentName[0] = '\0';
	szExperimentDescr[0] = '\0';


	if(experimentname->count)
		{
		strncpy(szExperimentName,experimentname->sval[0],sizeof(szExperimentName));
		szExperimentName[sizeof(szExperimentName)-1] = '\0';
		CUtility::TrimQuotedWhitespcExtd(szExperimentName);
		CUtility::ReduceWhitespace(szExperimentName);
		}
	else
		szExperimentName[0] = '\0';

	gExperimentID = 0;
	gProcessID = 0;
	gProcessingID = 0;
	szSQLiteDatabase[0] = '\0';
	szExperimentDescr[0] = '\0';

	if(summrslts->count)
		{
		strncpy(szSQLiteDatabase,summrslts->filename[0],sizeof(szSQLiteDatabase)-1);
		szSQLiteDatabase[sizeof(szSQLiteDatabase)-1] = '\0';
		CUtility::TrimQuotedWhitespcExtd(szSQLiteDatabase);
		if(strlen(szSQLiteDatabase) < 1)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: After removal of whitespace, no SQLite database specified with '-q<filespec>' option");
			return(1);
			}

		if(strlen(szExperimentName) < 1)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: After removal of whitespace, no SQLite experiment name specified with '-w<str>' option");
			return(1);
			}
		if(experimentdescr->count)
			{
			strncpy(szExperimentDescr,experimentdescr->sval[0],sizeof(szExperimentDescr)-1);
			szExperimentDescr[sizeof(szExperimentDescr)-1] = '\0';
			CUtility::TrimQuotedWhitespcExtd(szExperimentDescr);
			}
		if(strlen(szExperimentDescr) < 1)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: After removal of whitespace, no SQLite experiment description specified with '-W<str>' option");
			return(1);
			}

		gExperimentID = gSQLiteSummaries.StartExperiment(szSQLiteDatabase,false,true,szExperimentName,szExperimentName,szExperimentDescr);
		if(gExperimentID < 1)
			return(1);
		gProcessID = gSQLiteSummaries.AddProcess((char *)gpszSubProcess->pszName,(char *)gpszSubProcess->pszName,(char *)gpszSubProcess->pszFullDescr);
		if(gProcessID < 1)
			return(1);
		gProcessingID = gSQLiteSummaries.StartProcessing(gExperimentID,gProcessID,(char *)cpszProgVer);
		if(gProcessingID < 1)
			return(1);
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Initialised SQLite database '%s' for results summary collection",szSQLiteDatabase);
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"SQLite database experiment identifier for '%s' is %d",szExperimentName,gExperimentID);
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"SQLite database process identifier for '%s' is %d",(char *)gpszSubProcess->pszName,gProcessID);
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"SQLite database processing instance identifier is %d",gProcessingID);
		}
	else
		{
		szSQLiteDatabase[0] = '\0';
		szExperimentDescr[0] = '\0';
		}

	PMode = (etPBPMode)(pmode->count ? pmode->ival[0] : (int)ePBPMErrCorrect);
	if(PMode < ePBPMErrCorrect || PMode > ePBMConsolidate)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Processing mode '-m%d' must be in range 0..%d",PMode,ePBMConsolidate);
		return(1);
		}

	FiltMinHomoLen = 0;
	bSenseOnlyOvlps = false;
	MinSeedCoreLen = cDfltSeedCoreLen;
	MinNumSeedCores = cDfltNumSeedCores;
	DeltaCoreOfs = cDfltDeltaCoreOfs;
	MaxSeedCoreDepth = cDfltMaxSeedCoreDepth;
	SWMatchScore = cDfltSWMatchScore;
	SWMismatchPenalty = abs(cDfltSWMismatchPenalty);
	SWGapOpenPenalty = abs(cDfltSWGapOpenPenalty);
	SWGapExtnPenalty = abs(cDfltSWGapExtnPenalty);
	SWProgExtnPenaltyLen = cDfltSWProgExtnLen;
	TranscriptomeLens = 0;
	MinPBSeqLen = cDfltMinPBSeqLen;
	MaxPBSeqLen = cDfltMaxPBSeqLen;
	MinPBSeqOverlap = cDfltMinErrCorrectLen; 
	MinHCSeqLen = cDfltMinHCSeqLen;
	MinHCSeqOverlap = cDfltMinHCSeqOverlap;
	HCRelWeighting = cDfltHCRelWeighting;
	MaxArtefactDev = cDfltMaxArtefactDev;
	NumPacBioFiles = 0;
	NumHiConfFiles = 0;
	szOutMAFile[0] = '\0';

	szServiceName[0] = '\0';
	szHostName[0] = '\0';
	MaxRMI = 0;
	MaxNonRMI = 0;
	if(!(PMode == ePBMConsolidate || PMode == ePBPMConsensus) && rmihost->count)
		{
		strncpy(szHostName,rmihost->sval[0],sizeof(szHostName)-1);
		szHostName[sizeof(szHostName) - 1] = '\0';
		CUtility::TrimQuotedWhitespcExtd(szHostName);
		}

	szServiceName[0] = '\0';
	if (!(PMode == ePBMConsolidate || PMode == ePBPMConsensus) && rmiservice->count)
		{
		strncpy(szServiceName, rmiservice->sval[0], sizeof(szServiceName) - 1);
		szServiceName[sizeof(szServiceName) - 1] = '\0';
		CUtility::TrimQuotedWhitespcExtd(szServiceName);
		}

	if(szHostName[0] == '\0' && szServiceName[0] != '\0')
		strcpy(szHostName, "127.0.0.1");
	if(szHostName[0] != '\0' && szServiceName[0] == '\0')
		strcpy(szServiceName, "43123");

	SampleInRate = sampleinrate->count ? sampleinrate->ival[0] : 100;
	if(SampleInRate < 1 || SampleInRate > 100)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Input sequence acceptance sampling rate per 100 '-R%d' must be in range 1..100",SampleInRate);
		return(1);
		}

	SampleAcceptRate = sampleacceptrate->count ? sampleacceptrate->ival[0] : 1000;
	if(SampleAcceptRate < 1 || SampleAcceptRate > 1000)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Sample accepted sequence rate per 1000 '-Z%d' must be in range 1..1000",SampleAcceptRate);
		return(1);
		}

	if(PMode == ePBPMErrCorrect || PMode == ePBMConsolidate)
		{
		if(PMode == ePBPMErrCorrect)
			{
			// simply clamping FiltMinHomoLen to be 0 or in the range of 12..25 with a default of 16
			FiltMinHomoLen = minfilthomolen->count ? minfilthomolen->ival[0] : 16;
			if(FiltMinHomoLen <= 0)
				FiltMinHomoLen = 0;
			else
				if(FiltMinHomoLen < 12)
					FiltMinHomoLen = 12;
			if(FiltMinHomoLen > 25)
				FiltMinHomoLen = 25;	
			}
		else
			FiltMinHomoLen = 0;
		if(PMode == ePBPMErrCorrect)
			MinSeedCoreLen = minseedcorelen->count ? minseedcorelen->ival[0] : cDfltSeedCoreLen;
		else
			MinSeedCoreLen = minseedcorelen->count ? minseedcorelen->ival[0] : cDfltConsolidateSeedCoreLen;
		if(MinSeedCoreLen < cMinSeedCoreLen || MinSeedCoreLen > cMaxSeedCoreLen)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: seed core length '-c%d' must be in range %d..%d",MinSeedCoreLen,cMinSeedCoreLen,cMaxSeedCoreLen);
			return(1);
			}
		if(PMode == ePBPMErrCorrect)
			MinNumSeedCores = minseedcores->count ? minseedcores->ival[0] : cDfltNumSeedCores;
		else
			MinNumSeedCores = minseedcores->count ? minseedcores->ival[0] : cDfltConsolidateNumSeedCores;
		if(MinNumSeedCores < cMinNumSeedCores || MinNumSeedCores > cMaxNumSeedCores)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Minimum number of accepted seed cores '-C%d' must be in range %d..%d",MinNumSeedCores,cMinNumSeedCores,cMaxNumSeedCores);
			return(1);
			}

		if(PMode == ePBPMErrCorrect)		
			DeltaCoreOfs = deltacoreofs->count ? deltacoreofs->ival[0] : cDfltDeltaCoreOfs;
		else
			DeltaCoreOfs = deltacoreofs->count ? deltacoreofs->ival[0] : cDfltConsolidateDeltaCoreOfs;
		if(DeltaCoreOfs < 1 || DeltaCoreOfs > 10)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: offset seed cores '-d%d' must be in range 1..10",DeltaCoreOfs);
			return(1);
			}
		if(PMode == ePBPMErrCorrect)
			MaxSeedCoreDepth = maxcoredepth->count ? maxcoredepth->ival[0] : cDfltMaxSeedCoreDepth;
		else
			MaxSeedCoreDepth = maxcoredepth->count ? maxcoredepth->ival[0] : cDfltMaxConsolidateSeedCoreDepth;
		if(MaxSeedCoreDepth < 1000 || MaxSeedCoreDepth > 100000)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: maximum depth to explore seed cores '-D%d' must be in range 1000..100000",MaxSeedCoreDepth);
			return(1);
			}

		if(PMode == ePBPMErrCorrect)
			SWMatchScore = matchscore->count ? matchscore->ival[0] : cDfltSWMatchScore;
		else
			SWMatchScore = matchscore->count ? matchscore->ival[0] : cDfltConsolidateSWMatchScore;
		if(SWMatchScore < 1 || SWMatchScore > cMaxAllowedSWScore)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: SW Match score '-x%d' must be in range 1..%d",SWMatchScore,cMaxAllowedSWScore);
			return(1);
			}
		if(PMode == ePBPMErrCorrect)
			SWMismatchPenalty = mismatchpenalty->count ? mismatchpenalty->ival[0] : abs(cDfltSWMismatchPenalty);
		else
			SWMismatchPenalty = mismatchpenalty->count ? mismatchpenalty->ival[0] : abs(cDfltConsolidateSWMismatchPenalty);
		if(SWMismatchPenalty < 1 || SWMismatchPenalty > cMaxAllowedSWScore)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: SW Mismatch penalty '-X%d' must be in range 1..%d",SWMismatchPenalty,cMaxAllowedSWScore);
			return(1);
			}

		if(PMode == ePBPMErrCorrect)
			SWGapOpenPenalty = gapopenpenalty->count ? gapopenpenalty->ival[0] : abs(cDfltSWGapOpenPenalty);
		else
			SWGapOpenPenalty = gapopenpenalty->count ? gapopenpenalty->ival[0] : abs(cDfltConsolidateSWGapOpenPenalty);
		if(SWGapOpenPenalty < 1 || SWGapOpenPenalty > cMaxAllowedSWScore)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: SW Gap open penalty '-y%d' must be in range 1..%d",SWGapOpenPenalty,cMaxAllowedSWScore);
			return(1);
			}

		if(PMode == ePBPMErrCorrect)
			SWGapExtnPenalty = gapextnpenalty->count ? gapextnpenalty->ival[0] : abs(cDfltSWGapExtnPenalty);
		else
			SWGapExtnPenalty = gapextnpenalty->count ? gapextnpenalty->ival[0] : abs(cDfltConsolidateSWGapExtnPenalty);
		if(SWGapExtnPenalty < 1 || SWGapExtnPenalty > cMaxAllowedSWScore)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: SW Gap extension penalty '-Y%d' must be in range 1..%d",SWGapExtnPenalty,cMaxAllowedSWScore);
			return(1);
			}
		if(PMode == ePBPMErrCorrect)
			SWProgExtnPenaltyLen = progextnpenaltylen->count ? progextnpenaltylen->ival[0] : cDfltSWProgExtnLen;
		else
			SWProgExtnPenaltyLen = progextnpenaltylen->count ? progextnpenaltylen->ival[0] : cDfltConsolidateSWProgExtnLen;
		if(SWProgExtnPenaltyLen < 1 || SWProgExtnPenaltyLen > 63)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: SW Apply gap extension progress penalty for gaps at least '-z%d' must be in range 1..%d",SWProgExtnPenaltyLen,63);
			return(1);
			}

		if(PMode == ePBPMErrCorrect)
			TranscriptomeLens = transcriptomelens->count ? transcriptomelens->ival[0] : 0;
		else
			{
			TranscriptomeLens = transcriptomelens->count ? transcriptomelens->ival[0] : 5;
			if(TranscriptomeLens == 0)
				TranscriptomeLens = 1;
			}
		if(TranscriptomeLens < 0 || TranscriptomeLens > 15)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Transcript overlap read length differential percentage '-t%d' must be either 0 to disable or in range 1..15",TranscriptomeLens);
			return(1);
			}

		if(PMode == ePBPMErrCorrect)
			MinPBSeqLen = minpbseqlen->count ? minpbseqlen->ival[0] : cDfltMinPBSeqLen;
		else
			MinPBSeqLen = minpbseqlen->count ? minpbseqlen->ival[0] : cDfltMinConsolidatePBSeqLen;
		if(MinPBSeqLen < cMinPBSeqLen || MinPBSeqLen > cMaxMinPBSeqLen)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Minimum accepted PacBio length '-l%d' for error correction must be in range %d..%dbp",MinPBSeqLen,cMinPBSeqLen,cMaxMinPBSeqLen);
			return(1);
			}

		MaxPBSeqLen = maxpbseqlen->count ? maxpbseqlen->ival[0] : max(MinPBSeqLen,cDfltMaxPBSeqLen);
		if (MaxPBSeqLen < MinPBSeqLen || MaxPBSeqLen > cMaxMaxPBSeqLen)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: Maximum accepted PacBio length '-L%d' must be in range %d..%dbp", MaxPBSeqLen,MinPBSeqLen, cMaxMaxPBSeqLen);
			return(1);
			}


		if(PMode == ePBPMErrCorrect)
			MinPBSeqOverlap = minpbseqovl->count ? minpbseqovl->ival[0] : cDfltMinErrCorrectLen;
		else
			MinPBSeqOverlap = minpbseqovl->count ? minpbseqovl->ival[0] : cDfltMinConsolidatePBSeqLen;

		if(MinPBSeqOverlap < cMinPBSeqLen || MinPBSeqOverlap > cMaxMinPBSeqLen)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Minimum PacBio overlap required length '-b%d' must be in range %d..%dbp",MinPBSeqOverlap,cMinPBSeqLen,cMaxMinPBSeqLen);
			return(1);
			}

		if(PMode == ePBPMErrCorrect)
			MaxArtefactDev = maxartefactdev->count ? maxartefactdev->ival[0] : cDfltMaxArtefactDev;
		else
			MaxArtefactDev = maxartefactdev->count ? maxartefactdev->ival[0] : cDfltMaxConsolidateArtefactDev;

		if(MaxArtefactDev <= 0 || MaxArtefactDev > cMaxMaxArtefactDev)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Max overlap artefactual deviation '-A%d' must be either 0 or in range %d..%d",MaxArtefactDev,cMinMaxArtefactDev,cMaxMaxArtefactDev);
			return(1);
			}

		if(PMode == ePBPMErrCorrect)
			MinHCSeqLen = minhcseqlen->count ? minhcseqlen->ival[0] : 1000;
		else
			MinHCSeqLen = minhcseqlen->count ? minhcseqlen->ival[0] : cDfltMinConsolidatePBSeqLen;

		if(MinHCSeqLen < cMinPBSeqLen || MinHCSeqLen > cMaxMinPBSeqLen)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Minimum accepted high confidence length '-l%d' must be in range %d..%dbp",MinHCSeqLen,cMinPBSeqLen,cMaxMinPBSeqLen);
			return(1);
			}

		MinHCSeqOverlap = minhcseqovl->count ? minhcseqovl->ival[0] : 500;
		if(MinHCSeqOverlap < cMinPBSeqLen || MinHCSeqOverlap > cMaxMinPBSeqLen)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Minimum high confidence overlap required length '-L%d' must be in range %d..%dbp",MinHCSeqOverlap,cMinPBSeqLen,cMaxMinPBSeqLen);
			return(1);
			}

		HCRelWeighting = hcrelweighting->count ? hcrelweighting->ival[0] : 3;
		if(HCRelWeighting < 1 || HCRelWeighting > 10)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: High confidence relative weighting '-r%d' must be in range %d..%dbp",HCRelWeighting,1,10);
			return(1);
			}

		// widcards are allowed and filepaths will be expanded
		pszPacBioFiles[0] = NULL;
		for(NumPacBioFiles=Idx=0;NumPacBioFiles < cMaxInFileSpecs && Idx < pacbiofiles->count; Idx++)
			{
			strncpy(szTempFilePath,pacbiofiles->filename[Idx],_MAX_PATH);
			szTempFilePath[_MAX_PATH-1] = '\0';
			CUtility::TrimQuotedWhitespcExtd(szTempFilePath);
			if(szTempFilePath[0] == '\0')
				continue;
			glob.Init();
			if(glob.Add(szTempFilePath) < SG_SUCCESS)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to glob '%s",szTempFilePath);
				exit(1);
				}
			if(glob.FileCount() <= 0)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to locate any PacBio file matching '%s', checking for other files",szTempFilePath);
				continue;
				}
			for (FileID = 0; Rslt >= eBSFSuccess &&  FileID < glob.FileCount(); FileID++)
				{
				pszInfile = glob.File(FileID);
				if(pszPacBioFiles[NumPacBioFiles] == NULL)
					pszPacBioFiles[NumPacBioFiles] = new char [_MAX_PATH + 1];
				strncpy(pszPacBioFiles[NumPacBioFiles],pszInfile,_MAX_PATH);
				pszPacBioFiles[NumPacBioFiles][_MAX_PATH-1] = '\0';
				NumPacBioFiles += 1;
				if(NumPacBioFiles == cMaxInFileSpecs)
					{
					gDiagnostics.DiagOut(eDLWarn,gszProcName,"At limit (%d) of input PacBio files accepted, ignoring any additional",cMaxInFileSpecs);
					break;
					}
				else
					pszPacBioFiles[NumPacBioFiles] = NULL;
				}
			}

		if(!NumPacBioFiles)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: After removal of whitespace, no input PacBio file(s) specified with '-i<filespec>' option)\n");
			exit(1);
			}

		pszHiConfFiles[0] = NULL;
		for(NumHiConfFiles=Idx=0;NumHiConfFiles < cMaxInFileSpecs && Idx < hiconffiles->count; Idx++)
			{
			strncpy(szTempFilePath,hiconffiles->filename[Idx],_MAX_PATH);
			szTempFilePath[_MAX_PATH-1] = '\0';
			CUtility::TrimQuotedWhitespcExtd(szTempFilePath);
			if(szTempFilePath[0] == '\0')
				continue;
			glob.Init();
			if(glob.Add(szTempFilePath) < SG_SUCCESS)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to glob '%s",szTempFilePath);
				exit(1);
				}
			if(glob.FileCount() <= 0)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to locate any high confidence file matching '%s', checking for other files",szTempFilePath);
				continue;
				}
			for (FileID = 0; Rslt >= eBSFSuccess &&  FileID < glob.FileCount(); FileID++)
				{
				pszInfile = glob.File(FileID);
				if(pszHiConfFiles[NumHiConfFiles] == NULL)
					pszHiConfFiles[NumHiConfFiles] = new char [_MAX_PATH + 1];
				strncpy(pszHiConfFiles[NumHiConfFiles],pszInfile,_MAX_PATH);
				pszHiConfFiles[NumHiConfFiles][_MAX_PATH-1] = '\0';
				NumHiConfFiles += 1;
				if(NumHiConfFiles == cMaxInFileSpecs)
					{
					gDiagnostics.DiagOut(eDLWarn,gszProcName,"At limit (%d) of high confidence input files accepted, ignoring any additional",cMaxInFileSpecs);
					break;
					}
				else
					pszHiConfFiles[NumHiConfFiles] = NULL;
				}
			}


		if(!NumHiConfFiles)
			{
			NumHiConfFiles = 0;
			pszHiConfFiles[0] = NULL;
			}


		if(!(outfile->count || mafile->count))
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: No output files specified with either/both '-o<filespec>' or '-O<filespec>' option)\n");
			exit(1);
			}

		if(outfile->count)
			{
			strncpy(szOutFile,outfile->filename[0],sizeof(szOutFile));
			szOutFile[sizeof(szOutFile)-1] = '\0';
			CUtility::TrimQuotedWhitespcExtd(szOutFile);
			if(szOutFile[0] == '\0')
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: After removal of whitespace, no output file specified with '-o<filespec>' option)\n");
				exit(1);
				}
			}
		else
			szOutFile[0] = '\0';

		if(mafile->count)
			{
			strncpy(szOutMAFile,mafile->filename[0],sizeof(szOutFile));
			szOutMAFile[sizeof(szOutMAFile)-1] = '\0';
			CUtility::TrimQuotedWhitespcExtd(szOutMAFile);
			if(szOutMAFile[0] == '\0')
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: After removal of whitespace, no multialignment file specified with '-O<mafile>' option)\n");
				exit(1);
				}
			}
		else
			szOutMAFile[0] = '\0';

		if(szOutMAFile[0] != '\0' && szOutFile[0] != '\0' && !strcmp(szOutMAFile,szOutFile))
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Ensure different names used for the multialignment and the error corrected output files\n");
			exit(1);
			}
		}

	if(PMode != ePBPMOverlapDetail)
		{
		if(PMode == ePBMConsolidate)
			MinConcScore = minconcscore->count ? minconcscore->ival[0] : 0;
		else
			MinConcScore = minconcscore->count ? minconcscore->ival[0] : 3;
		if(MinConcScore < 0 || MinConcScore > 9)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: error corrected sequences trimming threshold score '-s%d' must be in range 0..9",MinConcScore);
			return(1);
			}

		if(PMode == ePBMConsolidate)
			MinErrCorrectLen = minerrcorrectlen->count ? minerrcorrectlen->ival[0] : cDfltMinConsolidatePBSeqLen;
		else
			MinErrCorrectLen = minerrcorrectlen->count ? minerrcorrectlen->ival[0] : cDfltMinErrCorrectLen;

		if(MinErrCorrectLen < cMinPBSeqLen || MinErrCorrectLen > cMaxMinPBSeqLen)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Error corrected and trimmed sequence minimum length '-S%d' must be in range %d..%dbp",MinErrCorrectLen,cMinPBSeqLen,cMaxMinPBSeqLen);
			return(1);
			}
		}

	if(PMode == ePBPMOverlapDetail)
		{
		MaxArtefactDev = maxartefactdev->count ? maxartefactdev->ival[0] : cDfltScaffMaxArtefactDev;
		if(MaxArtefactDev <= 0 || MaxArtefactDev > cMaxMaxArtefactDev)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Max overlap artefactual deviation '-A%d' must be either 0 or in range %d..%d",MaxArtefactDev,cMinMaxArtefactDev,cMaxMaxArtefactDev);
			return(1);
			}
		if(MaxArtefactDev < 0)
			MaxArtefactDev = 0;

		MinSeedCoreLen = minseedcorelen->count ? minseedcorelen->ival[0] : cDfltScaffSeedCoreLen;
		if(MinSeedCoreLen < cMinSeedCoreLen || MinSeedCoreLen > cMaxSeedCoreLen)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: seed core length '-c%d' must be in range %d..%d",MinSeedCoreLen,cMinSeedCoreLen,cMaxSeedCoreLen);
			return(1);
			}

		MinNumSeedCores = minseedcores->count ? minseedcores->ival[0] : cDfltScaffSeedCores;
		if(MinNumSeedCores < cMinNumSeedCores || MinNumSeedCores > cMaxNumSeedCores)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Minimum number of accepted seed cores '-C%d' must be in range %d..%d",MinNumSeedCores,cMinNumSeedCores,cMaxNumSeedCores);
			return(1);
			}

		DeltaCoreOfs = deltacoreofs->count ? deltacoreofs->ival[0] : 10;
		if(DeltaCoreOfs < 1 || DeltaCoreOfs > 10)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: offset seed cores '-d%d' must be in range 1..10",DeltaCoreOfs);
			return(1);
			}

		MaxSeedCoreDepth = maxcoredepth->count ? maxcoredepth->ival[0] : cDfltMaxSeedCoreDepth;
		if(MaxSeedCoreDepth < 1000 || MaxSeedCoreDepth > 100000)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: maximum depth to explore seed cores '-D%d' must be in range 1000..100000",MaxSeedCoreDepth);
			return(1);
			}

		SWMatchScore = matchscore->count ? matchscore->ival[0] : 1;
		if(SWMatchScore < 1 || SWMatchScore > cMaxAllowedSWScore)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: SW Match score '-x%d' must be in range 1..%d",SWMatchScore,cMaxAllowedSWScore);
			return(1);
			}
		SWMismatchPenalty = mismatchpenalty->count ? mismatchpenalty->ival[0] : 10;
		if(SWMismatchPenalty < 1 || SWMismatchPenalty > cMaxAllowedSWScore)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: SW Mismatch penalty '-X%d' must be in range 1..%d",SWMismatchPenalty,cMaxAllowedSWScore);
			return(1);
			}
		SWGapOpenPenalty = gapopenpenalty->count ? gapopenpenalty->ival[0] : 12;
		if(SWGapOpenPenalty < 1 || SWGapOpenPenalty > cMaxAllowedSWScore)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: SW Gap open penalty '-y%d' must be in range 1..%d",SWGapOpenPenalty,cMaxAllowedSWScore);
			return(1);
			}
		SWGapExtnPenalty = gapextnpenalty->count ? gapextnpenalty->ival[0] : 6;
		if(SWGapExtnPenalty < 1 || SWGapExtnPenalty > cMaxAllowedSWScore)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: SW Gap extension penalty '-Y%d' must be in range 1..%d",SWGapExtnPenalty,cMaxAllowedSWScore);
			return(1);
			}
		SWProgExtnPenaltyLen = progextnpenaltylen->count ? progextnpenaltylen->ival[0] : 1;
		if(SWProgExtnPenaltyLen < 1 || SWProgExtnPenaltyLen > 63)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: SW Apply gap extension progress penalty for gaps at least '-z%d' must be in range 1..%d",SWProgExtnPenaltyLen,63);
			return(1);
			}

		MinPBSeqLen = minpbseqlen->count ? minpbseqlen->ival[0] : PMode == ePBPMOverlapDetail ? cDfltMinErrCorrectLen : cDfltMinPBSeqLen;
		if(MinPBSeqLen < cMinPBSeqLen || MinPBSeqLen > cMaxMinPBSeqLen)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Minimum accepted error corrected sequence length '-l%d' must be in range %d..%dbp",MinPBSeqLen,cMinPBSeqLen,cMaxMinPBSeqLen);
			return(1);
			}

		MaxPBSeqLen = maxpbseqlen->count ? maxpbseqlen->ival[0] : cDfltMaxPBSeqLen;
		if (MaxPBSeqLen < MinPBSeqLen || MaxPBSeqLen > cMaxMaxPBSeqLen)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: Maximum accepted PacBio length '-L%d' must be in range %d..%dbp", MaxPBSeqLen, MinPBSeqLen, cMaxMaxPBSeqLen);
			return(1);
			}

		MinPBSeqOverlap = minpbseqovl->count ? minpbseqovl->ival[0] : min(cDfltMinErrCorrectLen,MinPBSeqLen);
		if(MinPBSeqOverlap < cMinPBSeqLen || MinPBSeqOverlap > cMaxMinPBSeqLen)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Minimum PacBio overlap required length '-L%d' must be in range %d..%dbp",MinPBSeqOverlap,cMinPBSeqLen,cMaxMinPBSeqLen);
			return(1);
			}


		// widcards are allowed and filepaths will be expanded
		pszPacBioFiles[0] = NULL;
		for(NumPacBioFiles=Idx=0;NumPacBioFiles < cMaxInFileSpecs && Idx < pacbiofiles->count; Idx++)
			{
			strncpy(szTempFilePath,pacbiofiles->filename[Idx],_MAX_PATH);
			szTempFilePath[_MAX_PATH-1] = '\0';
			CUtility::TrimQuotedWhitespcExtd(szTempFilePath);
			if(szTempFilePath[0] == '\0')
				continue;
			glob.Init();
			if(glob.Add(szTempFilePath) < SG_SUCCESS)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to glob '%s",szTempFilePath);
				exit(1);
				}
			if(glob.FileCount() <= 0)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to locate any PacBio file matching '%s', checking for other files",szTempFilePath);
				continue;
				}
			for (FileID = 0; Rslt >= eBSFSuccess &&  FileID < glob.FileCount(); FileID++)
				{
				pszInfile = glob.File(FileID);
				if(pszPacBioFiles[NumPacBioFiles] == NULL)
					pszPacBioFiles[NumPacBioFiles] = new char [_MAX_PATH + 1];
				strncpy(pszPacBioFiles[NumPacBioFiles],pszInfile,_MAX_PATH);
				pszPacBioFiles[NumPacBioFiles][_MAX_PATH-1] = '\0';
				NumPacBioFiles += 1;
				if(NumPacBioFiles == cMaxInFileSpecs)
					{
					gDiagnostics.DiagOut(eDLWarn,gszProcName,"At limit (%d) of input PacBio files accepted, ignoring any additional",cMaxInFileSpecs);
					break;
					}
				else
					pszPacBioFiles[NumPacBioFiles] = NULL;
				}
			}

		if(!NumPacBioFiles)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: After removal of whitespace, no input error corrected PacBio file(s) specified with '-i<filespec>' option)\n");
			exit(1);
			}

		NumHiConfFiles = 0;
		pszHiConfFiles[0] = NULL;

		if(!outfile->count)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: No output files specified with '-o<filespec>' option)\n");
			exit(1);
			}

		strncpy(szOutFile,outfile->filename[0],sizeof(szOutFile));
		szOutFile[sizeof(szOutFile)-1] = '\0';
		CUtility::TrimQuotedWhitespcExtd(szOutFile);
		if(szOutFile[0] == '\0')
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: After removal of whitespace, no output file specified with '-o<filespec>' option)\n");
			exit(1);
			}
	
		szOutMAFile[0] = '\0';
		}

	if(PMode == ePBMConsolidate || TranscriptomeLens > 0)
		bSenseOnlyOvlps = false;
	else
		bSenseOnlyOvlps = senseonlyovlps->count ? true : false;

	if(PMode == ePBPMConsensus)
		{
		if(pacbiofiles->count > 1)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: In processing mode 1 only one input mutltialignment file can be specified with '-i<inputfile>'");
			return(1);
			}
		else
			if(pacbiofiles->count == 0)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: In processing mode 1 one input mutltialignment file must be specified with '-i<inputfile>'");
				return(1);
				}

		pszPacBioFiles[0] = new char [_MAX_PATH];
		strncpy(pszPacBioFiles[NumPacBioFiles],pacbiofiles->filename[0],_MAX_PATH);
		pszPacBioFiles[0][_MAX_PATH-1] = '\0';
		CUtility::TrimQuotedWhitespcExtd(pszPacBioFiles[0]);
		if(pszPacBioFiles[0] == '\0')
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: After removal of whitespace, no input mutltialignment file must be specified with '-i<inputfile>'\n");
			exit(1);
			}
		NumPacBioFiles = 1;
		if(outfile->count == 0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: In processing mode 1 one output consensus sequence file must be specified with '-o<inputfile>'");
			return(1);
			}

		strncpy(szOutFile,outfile->filename[0],sizeof(szOutFile));
		szOutFile[sizeof(szOutFile)-1] = '\0';
		CUtility::TrimQuotedWhitespcExtd(szOutFile);
		if(szOutFile[0] == '\0')
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: After removal of whitespace, no output file specified with '-o<filespec>' option\n");
			exit(1);
			}
		}

// show user current resource limits
#ifndef _WIN32
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Resources: %s",CUtility::ReportResourceLimits());
#endif


#ifdef _WIN32
	SYSTEM_INFO SystemInfo;
	GetSystemInfo(&SystemInfo);
	NumberOfProcessors = SystemInfo.dwNumberOfProcessors;
#else
	NumberOfProcessors = sysconf(_SC_NPROCESSORS_ONLN);
#endif
	int MaxAllowedThreads = min(cMaxWorkerThreads,NumberOfProcessors);	// limit to be at most cMaxWorkerThreads
	if((NumThreads = threads->count ? threads->ival[0] : MaxAllowedThreads)==0)
		NumThreads = MaxAllowedThreads;
	if(NumThreads < 1 || NumThreads > MaxAllowedThreads)
		{
		gDiagnostics.DiagOut(eDLWarn,gszProcName,"Warning: Defaulting number of threads to %d",MaxAllowedThreads);
		NumThreads = MaxAllowedThreads;
		}

	if(szHostName[0] != '\0')
		{
		MaxRMI = maxrmi->count ? maxrmi->ival[0] : (NumThreads * cRMIThreadsPerCore);
		if(MaxRMI == 0)
			MaxRMI = (NumThreads * cRMIThreadsPerCore);
		if(MaxRMI < cRMIThreadsPerCore || MaxRMI > 500)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Number of RMI SW service instances '-n%d' must be in range %d..500",MaxRMI,cRMIThreadsPerCore);
			return(1);
			}
		MaxNonRMI = maxnonrmi->count ? maxnonrmi->ival[0] : NumThreads;
		if(MaxNonRMI > NumThreads)		// silently clamping to no more than the max number of threads requested
			MaxNonRMI = NumThreads;
		if(MaxNonRMI < 0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Number of non-RMI SW threads '-N%d' must be in range 0..%d",MaxNonRMI,NumThreads);
			return(1);
			}
		}

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing parameters:");

	char *pszMode;
	switch(PMode) {
		case ePBPMErrCorrect:		// error correction
			pszMode = (char *)"Error correct PacBio and generate consensus sequences";
			break;
		case ePBPMConsensus:		// generate consensus  from previously generated multiple alignments
			pszMode = (char *)"Generate consensus sequence from previously generated multiple alignments";
			break;
		case ePBPMOverlapDetail:		// Generate overlap detail from previously generated consensus sequences
			pszMode = (char *)"Generate overlap detail from previously generated consensus sequences";
			break;

		case ePBMConsolidate:		// consolidate one or more error corrected sequences into single representative sequence (transcript)
			pszMode = (char *)"Consolidate one or more error corrected sequences into a single representative sequence";
			break;

		}

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Processing mode: '%s'",pszMode);

	if(PMode == ePBPMErrCorrect || PMode == ePBMConsolidate)
		{
		strncpy(szChkPtsFile,szOutFile,sizeof(szChkPtsFile)-5);
		strcat(szChkPtsFile,".chk");
		}
	else
		szChkPtsFile[0] = '\0';
	
	if(szHostName[0] != '\0')
		{
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Listening on this host name for RMI SW service connections: '%s'",szHostName);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Listening on this service/port for RMI SW service connections: '%s'",szServiceName);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Allowing for a maximum of this many RMI SW service instances: '%d'",MaxRMI);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Allowing a maximum of this many non-RMI SW threads: '%d'",MaxNonRMI);
		}


	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Sampling input sequences for acceptance rate per 100: %d",SampleInRate);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Sampling accepted sequence rate per 1000: %d",SampleAcceptRate);


	if(PMode != ePBPMConsensus)
		{
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Overlap processing: '%s'",bSenseOnlyOvlps ? "Sense only" : "Sense and antisense");
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Use seed cores of this length when identifying putative overlapping sequences: %dbp",MinSeedCoreLen);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Require at least this many seed cores between overlapping sequences: %d",MinNumSeedCores);

		if(PMode == ePBPMErrCorrect)
			{
			if(FiltMinHomoLen != 0)
				gDiagnostics.DiagOutMsgOnly(eDLInfo,"Filtering PacBio reads for near homopolymer runs which are at least this length: %dbp",FiltMinHomoLen);
			else
				gDiagnostics.DiagOutMsgOnly(eDLInfo,"Filtering PacBio reads for near homopolymer runs which are at least this length: No filtering");
			}

		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Offset cores by this many bp: %d",DeltaCoreOfs);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Maximum seed core depth: %d",MaxSeedCoreDepth);

		gDiagnostics.DiagOutMsgOnly(eDLInfo,"SW score for matching bases: %d",SWMatchScore);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"SW mismatch penalty: %d",SWMismatchPenalty);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"SW gap opening penalty: %d",SWGapOpenPenalty);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"SW gap extension penalty: %d",SWGapExtnPenalty);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"SW gap extension penalty only applied for gaps of at least this size: %d",SWProgExtnPenaltyLen);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"classify overlaps as artefactual if sliding window of 500bp over any overlap deviates by more than this percentage: %d",MaxArtefactDev);
		if(PMode == ePBPMErrCorrect || PMode == ePBMConsolidate)
			{
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Minimum PacBio sequence length for error correction: %dbp",MinPBSeqLen);
			gDiagnostics.DiagOutMsgOnly(eDLInfo, "Maximum PacBio sequence length: %dbp", MaxPBSeqLen);
			if(TranscriptomeLens > 0)
				{
				gDiagnostics.DiagOutMsgOnly(eDLInfo, "Error correcting transcriptome reads, putative overlapping reads must have max length differential: %d%%", TranscriptomeLens);
				gDiagnostics.DiagOutMsgOnly(eDLInfo,"Minimum PacBio overlap required for error correction contribution: must be near full length");
				}
			else
				gDiagnostics.DiagOutMsgOnly(eDLInfo,"Minimum PacBio overlap required for error correction contribution: %d",MinPBSeqOverlap);	
			}
		else
			{
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Minimum error corrected sequence length: %dbp",MinPBSeqLen);
			gDiagnostics.DiagOutMsgOnly(eDLInfo, "Maximum error corrected sequence length: %dbp", MaxPBSeqLen);
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Minimum overlap required: %d",MinPBSeqOverlap);
			}


		if(NumHiConfFiles)
			{
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Minimum high confidence sequence length: %dbp",MinHCSeqLen);
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Minimum high confidence sequence onto PacBio overlap required for error correction contribution: %d",MinHCSeqOverlap);
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"High confidence sequence relative weighting: %d",HCRelWeighting);

			}

		for(Idx=0; Idx < NumPacBioFiles; Idx++)
			{
			if(PMode == ePBPMErrCorrect || PMode == ePBMConsolidate)
				gDiagnostics.DiagOutMsgOnly(eDLInfo,"Input PacBio sequences file spec: '%s'",pszPacBioFiles[Idx]);
			else
				gDiagnostics.DiagOutMsgOnly(eDLInfo,"Input consensus sequences file spec: '%s'",pszPacBioFiles[Idx]);
			}

		if(PMode != ePBPMOverlapDetail)
			{
			if(NumHiConfFiles)
				{
				for(Idx=0; Idx < NumHiConfFiles; Idx++)
					gDiagnostics.DiagOutMsgOnly(eDLInfo,"Input hiconfidence sequences file spec: '%s'",pszHiConfFiles[Idx]);
				}
			else
				gDiagnostics.DiagOutMsgOnly(eDLInfo,"Input hiconfidence sequences file spec: None Specified");

			if(szOutMAFile[0] != '\0')
				gDiagnostics.DiagOutMsgOnly(eDLInfo,"Output multialignment file: '%s'",szOutMAFile);
			}
		}

	if(PMode != ePBPMOverlapDetail)
		{
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Trimming error corrected PacBio sequences until mean 50bp score at least: %d",MinConcScore);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Error corrected and trimmed PacBio sequences must be at least this long: %d",MinErrCorrectLen);
		}

	if(PMode == ePBPMConsensus)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Input multialignment file spec: '%s'",pszPacBioFiles[ 0]);

	if(PMode != ePBPMOverlapDetail)
		{
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Output error corrected sequences file: '%s'",szOutFile);
		if(PMode == ePBPMErrCorrect || PMode == ePBMConsolidate)
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Checkpoint error corrected sequences file: '%s'",szChkPtsFile);
		}
	else
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Output error corrected sequences overlap detail file: '%s'",szOutFile);

	if(szExperimentName[0] != '\0')
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"This processing reference: %s",szExperimentName);
	
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"number of threads : %d",NumThreads);

	if(gExperimentID > 0)
		{
		int ParamID;
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szLogFile),"log",szLogFile);

		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(PMode),"pmode",&PMode);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(SampleInRate),"sampleinrate",&SampleInRate);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(SampleAcceptRate),"sampleacceptrate",&SampleAcceptRate);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(FiltMinHomoLen),"minfilthomolen",&FiltMinHomoLen);

		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(bSenseOnlyOvlps),"senseonlyovlps",&bSenseOnlyOvlps);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(MinSeedCoreLen),"seedcorelen",&MinSeedCoreLen);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(MinNumSeedCores),"minseedcores",&MinNumSeedCores);

		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(MaxSeedCoreDepth),"maxcoredepth",&MaxSeedCoreDepth);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(MinSeedCoreLen),"seedcorelen",&MinSeedCoreLen);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(MinNumSeedCores),"minseedcores",&MinNumSeedCores);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(DeltaCoreOfs),"deltacoreofs",&DeltaCoreOfs);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(SWMatchScore),"matchscore",&SWMatchScore);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(SWMismatchPenalty),"mismatchpenalty",&SWMismatchPenalty);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(SWGapOpenPenalty),"gapopenpenalty",&SWGapOpenPenalty);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(SWGapExtnPenalty),"gapextnpenalty",&SWGapExtnPenalty);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(SWProgExtnPenaltyLen),"progextnpenaltylen",&SWProgExtnPenaltyLen);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(MinNumSeedCores),"minseedcores",&MinNumSeedCores);


		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(TranscriptomeLens),"transcriptomelens",&TranscriptomeLens);		
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(MinPBSeqLen),"minpbseqlen",&MinPBSeqLen);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID, ePTInt32, (int)sizeof(MaxPBSeqLen), "maxpbseqlen", &MaxPBSeqLen);

		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(MinPBSeqOverlap),"minpbseqovl",&MinPBSeqOverlap);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(MaxArtefactDev),"artefactdev",&MaxArtefactDev);

		if(NumHiConfFiles)
			{
			ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(MinHCSeqLen),"minhcseqlen",&MinHCSeqLen);
			ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(MinHCSeqOverlap),"minhcseqovl",&MinHCSeqOverlap);
			ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(HCRelWeighting),"hcrelweighting",&HCRelWeighting);
			}

		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(MinConcScore),"minconcscore",&MinConcScore);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(MinErrCorrectLen),"minerrcorrectlen",&MinErrCorrectLen);


		for(Idx=0; Idx < NumPacBioFiles; Idx++)
			ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(pszPacBioFiles[Idx]),"pacbiofile",pszPacBioFiles[Idx]);

		if(NumHiConfFiles)
				{
				for(Idx=0; Idx < NumHiConfFiles; Idx++)
					ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(pszHiConfFiles[Idx]),"hiconffile",pszHiConfFiles[Idx]);
				}

		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szOutFile),"out",szOutFile);
		if(PMode == ePBPMErrCorrect || PMode == ePBMConsolidate)
			gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szChkPtsFile),"out",szChkPtsFile);

		if(szHostName[0] != '\0')
			{
			ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szHostName),"rmihostname",szHostName);
			ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szServiceName),"rmiservice",szServiceName);
			ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(MaxRMI),"maxrmi",&MaxRMI);
			ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(MaxNonRMI),"maxnonrmi",&MaxNonRMI);
			}

		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(NumThreads),"threads",&NumThreads);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(NumberOfProcessors),"cpus",&NumberOfProcessors);

		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szSQLiteDatabase),"sumrslts",szSQLiteDatabase);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szExperimentName),"experimentname",szExperimentName);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szExperimentDescr),"experimentdescr",szExperimentDescr);
		}


#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	gStopWatch.Start();
	Rslt = ProcPacBioErrCorrect((etPBPMode)PMode,szHostName,szServiceName,MaxRMI,MaxNonRMI,SampleInRate,SampleAcceptRate,FiltMinHomoLen,bSenseOnlyOvlps,DeltaCoreOfs,MaxSeedCoreDepth,MinSeedCoreLen,MinNumSeedCores,SWMatchScore,-1 * SWMismatchPenalty,-1 * SWGapOpenPenalty,-1 * SWGapExtnPenalty,SWProgExtnPenaltyLen,
								TranscriptomeLens,MinPBSeqLen, MaxPBSeqLen,MinPBSeqOverlap,MaxArtefactDev,MinHCSeqLen,MinHCSeqOverlap,HCRelWeighting,MinErrCorrectLen,MinConcScore,
								NumPacBioFiles,pszPacBioFiles,NumHiConfFiles,pszHiConfFiles,szOutFile,szOutMAFile,szChkPtsFile,NumThreads);
	Rslt = Rslt >=0 ? 0 : 1;
	if(gExperimentID > 0)
		{
		if(gProcessingID)
			gSQLiteSummaries.EndProcessing(gProcessingID,Rslt);
		gSQLiteSummaries.EndExperiment(gExperimentID);
		}
	gStopWatch.Stop();

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Exit code: %d Total processing time: %s",Rslt,gStopWatch.Read());
	exit(Rslt);
	}
else
	{
    printf("\n%s %s %s, Version %s\n", gszProcName,gpszSubProcess->pszName,gpszSubProcess->pszFullDescr,cpszProgVer);
	arg_print_errors(stdout,end,gszProcName);
	arg_print_syntax(stdout,argtable,"\nUse '-h' to view option and parameter usage\n");
	exit(1);
	}
return 0;
}


int
ProcPacBioErrCorrect(etPBPMode PMode,		// processing mode
		char *pszHostName,			// listening on this host name or IPv4/IPv5 address for connections by service providers 
		char *pszServiceName,		// Listen on this service name or port for for connections by service providers
		int MaxRMI,					// max number of RMI service provider instances supported
		int MaxNonRMI,				// max number of non-RMI SW threads
		int SampleInRate,			// sample input sequences at this rate per 100
		int SampleAcceptRate,		// sample accepted input sequences at this rate per 1000
		int FiltMinHomoLen,			// filter PacBio reads for homopolymer runs >= this length (0 to disable filtering) 		
		bool bSenseOnlyOvlps,		// process for sense only overlaps (default is for sense/sense and sense/antisense overlaps)
		int DeltaCoreOfs,			// offset by this many bp the core windows of coreSeqLen along the probe sequence when checking for overlaps
		int MaxSeedCoreDepth,		// only further process a seed core if there are no more than this number of matching cores in all targeted sequences
		int MinSeedCoreLen,			// use seed cores of this length when identifying putative overlapping scaffold sequences
		int MinNumSeedCores,        // require at least this many seed cores between overlapping scaffold sequences
		int SWMatchScore,			// score for matching bases (0..50)
		int SWMismatchPenalty,		// mismatch penalty (-50..0)
		int SWGapOpenPenalty,		// gap opening penalty (-50..0)
		int SWGapExtnPenalty,		// gap extension penalty (-50..0)
		int SWProgExtnPenaltyLen,	// progressive gap scoring then only apply gap extension score if gap at least this length (0..63) - use if aligning PacBio
		int TranscriptomeLens,		// 0 if disabled, processing transcript reads, putatively overlapping reads must have length differential no more than this percentage and overlaps to be nearly full length
		int MinPBSeqLen,			// only accepting PacBio reads for error correction of at least this length (defaults to 10Kbp) and if
        int MaxPBSeqLen,			// no more than this length (defaults to 35Kbp)
		int MinPBSeqOverlap,		// any overlap of a PacBio onto a target PacBio must be of at least this many bp to be considered for contributing towards error correction (defaults to 5Kbp) 
		int MaxArtefactDev,			// classify overlaps as artefactual if sliding window of 500bp over any overlap deviates by more than this percentage from the overlap mean
		int MinHCSeqLen,			// only accepting hiconfidence reads of at least this length (defaults to 1Kbp)
		int MinHCSeqOverlap,		// any overlap of a hiconfidence read onto a target PacBio read must be of at least this many bp to be considered for contributing towards error correction (defaults to 1Kbp) 
		int HCRelWeighting,         // hiconfidence read overlaps are usually weighted higher than normal lesser confidence read overlaps when calling consensus bases 
	    int MinErrCorrectLen,		// error corrected and trimmed sequences must be at least this minimum length
		int MinConcScore,			// error corrected sequences trimmed until mean 50bp concensus score is at least this threshold
		int NumPacBioFiles,			// number of input pacbio file specs
		char *pszPacBioFiles[],		// input pacbio files
		int NumHiConfFiles,			// number of input hiconfidence file specs
		char *pszHiConfFiles[],		// input hiconfidence files		
	    char *pszOutFile,			// where to write error corrected sequences
		char *pszOutMAFile,			// where to write multiple alignments
		char *pszChkPtsFile,        // name of file used for checkpointing in case resume processing is required
		int NumThreads)				// maximum number of worker threads to use
{
int Rslt;
CPBErrCorrect *pPBErrCorrect;

if((pPBErrCorrect = new CPBErrCorrect)==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: Unable to instantiate PBErrCorrect");
	return(eBSFerrObj);
	}

Rslt = pPBErrCorrect->Process(PMode,pszHostName,pszServiceName,MaxRMI,MaxNonRMI,SampleInRate,SampleAcceptRate,FiltMinHomoLen,bSenseOnlyOvlps,DeltaCoreOfs,MaxSeedCoreDepth,MinSeedCoreLen,MinNumSeedCores,SWMatchScore,SWMismatchPenalty,SWGapOpenPenalty,SWGapExtnPenalty,SWProgExtnPenaltyLen,
								TranscriptomeLens,MinPBSeqLen, MaxPBSeqLen, MinPBSeqOverlap,MaxArtefactDev,MinHCSeqLen,MinHCSeqOverlap,HCRelWeighting,MinErrCorrectLen,MinConcScore,
								NumPacBioFiles,pszPacBioFiles,NumHiConfFiles,pszHiConfFiles,pszOutFile,pszOutMAFile,pszChkPtsFile,NumThreads);
delete pPBErrCorrect;
return(Rslt);
}

CPBErrCorrect::CPBErrCorrect() // relies on base classes constructors
{
m_pSfxArray = NULL;
m_pPBScaffNodes = NULL;
m_pMapEntryID2NodeIDs = NULL;
m_pRequester = NULL;
m_bMutexesCreated = false;
m_hErrCorFile = -1;
m_hChkPtsFile = -1;
m_hMultiAlignFile = -1;

Init();
}

CPBErrCorrect::~CPBErrCorrect() // relies on base classes destructor
{
Reset();
}


void
CPBErrCorrect::Init(void)
{
if(m_pSfxArray != NULL)
	{
	delete m_pSfxArray;
	m_pSfxArray = NULL;
	}

if(m_pPBScaffNodes != NULL)
	{
	delete m_pPBScaffNodes;
	m_pPBScaffNodes = NULL;
	}
if(m_pMapEntryID2NodeIDs != NULL)
	{
	delete m_pMapEntryID2NodeIDs;
	m_pMapEntryID2NodeIDs = NULL;
	}

if(m_hErrCorFile != -1)
	{
#ifdef _WIN32
	_commit(m_hErrCorFile);
#else
	fsync(m_hErrCorFile);
#endif
	close(m_hErrCorFile);
	m_hErrCorFile = -1;
	}

if(m_hChkPtsFile != -1)
	{
#ifdef _WIN32
	_commit(m_hChkPtsFile);
#else
	fsync(m_hChkPtsFile);
#endif
	close(m_hChkPtsFile);
	m_hChkPtsFile = -1;
	}

if(m_hMultiAlignFile != -1)
	{
#ifdef _WIN32
	_commit(m_hMultiAlignFile);
#else
	fsync(m_hMultiAlignFile);
#endif
	close(m_hMultiAlignFile);
	m_hMultiAlignFile = -1;
	}

if(m_pRequester != NULL)
	{
	m_pRequester->Terminate(120);
	delete m_pRequester;
	m_pRequester = NULL;
	}

m_bRMI = false;
m_szRMIHostName[0] = '\0';
m_szRMIServiceName[0] = '\0';

m_szErrCorFile[0] = '\0';
m_szChkPtsFile[0] = '\0';
m_szMultiAlignFile[0] = '\0';

m_NumPBScaffNodes = 0;
m_AllocdPBScaffNodes = 0;
m_MaxPBSeqLen = 0;

m_NumOverlapProcessed = 0;
m_ProvOverlapping = 0;
m_ProvOverlapped = 0;
m_ProvContained = 0;
m_ProvArtefact = 0;
m_ProvSWchecked = 0;
m_MultiAlignFileUnsyncedSize = 0;
m_ErrCorFileUnsyncedSize = 0;

m_bSenseOnlyOvlps = false;

m_MinErrCorrectLen = cDfltMinErrCorrectLen;
m_MinConcScore = 3;

m_PMode = ePBPMErrCorrect;

m_szScaffLineBuff[0] = '\0';
m_ScaffLineBuffIdx = 0;

m_OverlapFloat = cDfltMaxOverlapFloat;
m_TranscriptomeLens = 0;
m_MinPBSeqLen = cDfltMinPBSeqLen;
m_MaxPBRdSeqLen = cDfltMaxPBSeqLen;
m_MinPBSeqOverlap = cDfltMinErrCorrectLen;
m_MaxArtefactDev = cDfltMaxArtefactDev;
m_MinHCSeqLen = cDfltMinHCSeqLen;
m_MinHCSeqOverlap = cDfltMinHCSeqOverlap;
m_HCRelWeighting = cDfltHCRelWeighting;

m_DeltaCoreOfs = cDfltDeltaCoreOfs;
m_MaxSeedCoreDepth = cDfltMaxSeedCoreDepth;
m_MinSeedCoreLen = cDfltSeedCoreLen;
m_MinNumSeedCores = cDfltNumSeedCores;

m_SWMatchScore = cDfltSWMatchScore;
m_SWMismatchPenalty = cDfltSWMismatchPenalty;	
m_SWGapOpenPenalty = cDfltSWGapOpenPenalty;
m_SWGapExtnPenalty = cDfltSWGapExtnPenalty;
m_SWProgExtnPenaltyLen = cDfltSWProgExtnLen;	

memset(m_ExactKmerDists,0,sizeof(m_ExactKmerDists));
m_TotAlignSeqLen = 0;
m_TotAlignSeqs = 0;
 
m_NumPacBioFiles = 0;
m_NumHiConfFiles = 0;

memset(m_szPacBioFiles,0,sizeof(m_szPacBioFiles));
memset(m_szHiConfFiles,0,sizeof(m_szHiConfFiles));

m_MaxRMIInstances = 0;

m_ProcessStatsThen = 0;
m_NumCPUCores = 0;

if(m_bMutexesCreated)
	DeleteMutexes();
m_bMutexesCreated = false; 
}

void
CPBErrCorrect::Reset(void)			// reset state back to that immediately following instantiation
{
Init();
}


int
CPBErrCorrect::CreateMutexes(void)
{
if(m_bMutexesCreated)
	return(eBSFSuccess);

m_CASSerialise = 0;
m_CASLock = 0;
m_CASThreadPBErrCorrect = 0;
m_bMutexesCreated = true;
return(eBSFSuccess);
}

void
CPBErrCorrect::DeleteMutexes(void)
{
if(!m_bMutexesCreated)
	return;
m_bMutexesCreated = false;
}


void
CPBErrCorrect::AcquireCASSerialise(void)
{
int SpinCnt = 10;
int BackoffMS = 1;

#ifdef _WIN32
while(InterlockedCompareExchange(&m_CASSerialise,1,0)!=0)
	{
	if(SpinCnt -= 1)
		continue;
	CUtility::SleepMillisecs(BackoffMS);
	SpinCnt = 10;
	if(BackoffMS < 50)
		BackoffMS += 1;
	else
		BackoffMS = 1 + (rand() % 31);
	}
#else
while(__sync_val_compare_and_swap(&m_CASSerialise,0,1)!=0)
	{
	if(SpinCnt -= 1)
		continue;
	CUtility::SleepMillisecs(BackoffMS);
	SpinCnt = 10;
	if(BackoffMS < 50)
		BackoffMS += 1;
	else
		BackoffMS = 1 + (rand() % 31);
	}
#endif
}

void
CPBErrCorrect::ReleaseCASSerialise(void)
{
#ifdef _WIN32
InterlockedCompareExchange(&m_CASSerialise,0,1);
#else
__sync_val_compare_and_swap(&m_CASSerialise,1,0);
#endif
}


void
CPBErrCorrect::AcquireCASLock(void)
{
int SpinCnt = 10;
int BackoffMS = 1;

#ifdef _WIN32
while(InterlockedCompareExchange(&m_CASLock,1,0)!=0)
	{
	if(SpinCnt -= 1)
		continue;
	CUtility::SleepMillisecs(BackoffMS);
	SpinCnt = 10;
	if(BackoffMS < 50)
		BackoffMS += 1;
	else
		BackoffMS = 1 + (rand() % 31);
	}
#else
while(__sync_val_compare_and_swap(&m_CASLock,0,1)!=0)
	{
	if(SpinCnt -= 1)
		continue;
	CUtility::SleepMillisecs(BackoffMS);
	SpinCnt = 10;
	if(BackoffMS < 50)
		BackoffMS += 1;
	else
		BackoffMS = 1 + (rand() % 31);
	}
#endif
}

void
CPBErrCorrect::ReleaseCASLock(void)
{
#ifdef _WIN32
InterlockedCompareExchange(&m_CASLock,0,1);
#else
__sync_val_compare_and_swap(&m_CASLock,1,0);
#endif
}


void
CPBErrCorrect::AcquireCASThreadPBErrCorrect(void)
{
int SpinCnt = 10;
int BackoffMS = 1;

#ifdef _WIN32
while(InterlockedCompareExchange(&m_CASThreadPBErrCorrect,1,0)!=0)
	{
	if(SpinCnt -= 1)
		continue;
	CUtility::SleepMillisecs(BackoffMS);
	SpinCnt = 10;
	if(BackoffMS < 50)
		BackoffMS += 1;
	else
		BackoffMS = 1 + (rand() % 31);
	}
#else
while(__sync_val_compare_and_swap(&m_CASThreadPBErrCorrect,0,1)!=0)
	{
	if(SpinCnt -= 1)
		continue;
	CUtility::SleepMillisecs(BackoffMS);
	SpinCnt = 10;
	if(BackoffMS < 50)
		BackoffMS += 1;
	else
		BackoffMS = 1 + (rand() % 31);
	}
#endif
}

void
CPBErrCorrect::ReleaseCASThreadPBErrCorrect(void)
{
#ifdef _WIN32
InterlockedCompareExchange(&m_CASThreadPBErrCorrect,0,1);
#else
__sync_val_compare_and_swap(&m_CASThreadPBErrCorrect,1,0);
#endif
}

// ProcessBioseqFile
// Process input biosequence file into suffix file
int
CPBErrCorrect::ProcessBioseqFile(int MinSeqLen,		// only accept for indexing sequences of at least this length
								 int MaxSeqLen,					// and which are no longer than this length
				 char *pszFile,						// file containing sequences to be parsed and indexed
				int Flags)							// default is for flags = cFlgLCSeq used with PacBio read sequences
{
CBioSeqFile BioSeqFile;
etSeqBase *pSeqBuff = NULL;
UINT32 AllocLen = 0;
char szSource[cBSFSourceSize];
char szDescription[cBSFDescriptionSize];
UINT32 SeqLen;
int Rslt;
UINT32 SeqID;
UINT32 NumSampled;
UINT32 NumSeqsUnderlength;
UINT32 NumSeqsOverlength;
UINT32 NumSeqsAccepted;
size_t TotAcceptedLen;
tBSFEntryID CurEntryID;

if((Rslt=BioSeqFile.Open(pszFile,cBSFTypeSeq,false))!=eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessBioseqFile - error %d %s",Rslt,BioSeqFile.GetErrMsg());
	return(Rslt);
	}

CurEntryID = 0;
SeqID = 0;
NumSeqsUnderlength = 0;
NumSeqsOverlength = 0;
NumSampled = 0;
NumSeqsAccepted = 0;
TotAcceptedLen = 0;
while((Rslt = CurEntryID = BioSeqFile.Next(CurEntryID)) > eBSFSuccess)
	{
	SeqID += 1;
	BioSeqFile.GetNameDescription(CurEntryID,cBSFSourceSize-1,(char *)&szSource,
											cBSFDescriptionSize-1,(char *)&szDescription);
	SeqLen = BioSeqFile.GetDataLen(CurEntryID);
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing %s|%s",szSource,szDescription);
	if((NumSeqsAccepted + NumSeqsUnderlength) > ((SeqID * (UINT64)m_SampleInRate) / 100))
		continue;
	NumSampled += 1;
	if(SeqLen < (UINT32)MinSeqLen)		// only accept for indexing sequences of at least this length)
		{
		NumSeqsUnderlength += 1;
		continue;
		}

	if (SeqLen > (UINT32)MaxSeqLen)		// only accept for indexing sequences of no more than this length)
	{
		NumSeqsOverlength += 1;
		continue;
	}


	if(AllocLen < (SeqLen + 1))
		{
		if(pSeqBuff != NULL)
			delete pSeqBuff;
		AllocLen = SeqLen;
		pSeqBuff = new unsigned char [AllocLen];
		if(pSeqBuff == NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessBioseqFile - Unable to alloc %u bytes memory for pSeqBuff",AllocLen);
			Rslt = eBSFerrMem;
			break;
			}
		}

	if((Rslt = BioSeqFile.GetData(CurEntryID,eSeqBaseType,0,pSeqBuff,SeqLen)) != SeqLen)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessBioseqFile - error %d %s",Rslt,BioSeqFile.GetErrMsg());
		break;
		}



	// remove any repeat masking flags so that sorts can actually sort
	// if run of more than 25 Ns and at least 5 Ns to end of buffer then randomly mutate
	// every 13th N
	//	e.g <25Ns>r<12Ns>r<12Ns> where r is a pseudorandom base
	etSeqBase *pMskBase = pSeqBuff;
	int SeqNs = 0;
	for(UINT32 MskIdx = 0; MskIdx < SeqLen; MskIdx++,pMskBase++)
		{
		*pMskBase &= ~cRptMskFlg;
		if(*pMskBase == eBaseN && (MskIdx+5) < SeqLen)
			{
			if(++SeqNs > 25 &&
				pMskBase[1] == eBaseN &&
				pMskBase[2] == eBaseN &&
				pMskBase[3] == eBaseN &&
				pMskBase[4] == eBaseN)
				{
				if(!(SeqNs % 13))	// mutate every 13th
					*pMskBase = rand() % 4;
				}
			}
		else
			SeqNs = 0;
		}

	if((Rslt=m_pSfxArray->AddEntry(szSource,pSeqBuff,SeqLen,Flags)) < eBSFSuccess)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessBioseqFile - error %d %s",Rslt,m_pSfxArray->GetErrMsg());
		break;
		}
	NumSeqsAccepted += 1;
	TotAcceptedLen += SeqLen;
	}
if(Rslt == eBSFerrEntry)
	Rslt = eBSFSuccess;
if(pSeqBuff != NULL)
	delete pSeqBuff;
BioSeqFile.Close();
gDiagnostics.DiagOut(eDLInfo,gszProcName,"ProcessFastaFile - %d parsed, %d sampled, %d accepted, %dbp mean length, %d sequences not accepted for indexing as length under %dbp, %d sequences length over %dbp ",
					SeqID,NumSampled,NumSeqsAccepted,NumSeqsAccepted == 0 ? 0 :(int)(TotAcceptedLen/NumSeqsAccepted),NumSeqsUnderlength,MinSeqLen, NumSeqsOverlength, MaxSeqLen);
return(Rslt);
}

// ProcessFastaFile
// Parse input fasta format file into a biosequence suffix array file
// The sequence can be optionally processed for near homopolymer runs and these runs removed from the sequence

int
CPBErrCorrect::ProcessFastaFile(int MinSeqLen,		// only accept for indexing sequences of at least this length
					int MaxSeqLen,					// and which are no longer than this length
				char *pszFile,						// file containing sequences
				int Flags)							// default is for flags = cFlgLCSeq used with PacBio read sequences
{
CFasta Fasta;
unsigned char *pSeqBuff;
unsigned char *pMskBase;
UINT32 MskIdx;
size_t BuffOfs;
size_t AllocdBuffSize;
size_t AvailBuffSize;
char szName[cBSFSourceSize];
char szDescription[cBSFDescriptionSize];
UINT32 SeqLen;
int Descrlen;
bool bFirstEntry;
bool bEntryCreated;
int Rslt;
int SeqID;
int NumSampled;
int NumSeqsAccepted;
size_t TotAcceptedLen;
UINT32 NumSeqsUnderlength;
UINT32 NumSeqsOverlength;

if((Rslt=Fasta.Open(pszFile,true))!=eBSFSuccess)
	{
	if(Rslt != eBSFerrNotFasta)
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessFastaFile: Unable to open '%s' [%s] %s",pszFile,Fasta.ErrText((teBSFrsltCodes)Rslt),Fasta.GetErrMsg());
	return(Rslt);
	}

AllocdBuffSize = (size_t)cMaxAllocBuffChunk * 16;
// note malloc is used as can then simply realloc to expand as may later be required
if((pSeqBuff = (unsigned char *)malloc(AllocdBuffSize)) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessFastaFile:- Unable to allocate memory (%u bytes) for sequence buffer",(UINT32)AllocdBuffSize);
	Fasta.Close();
	return(eBSFerrMem);
	}
AvailBuffSize = AllocdBuffSize;

gDiagnostics.DiagOut(eDLInfo,gszProcName,"ProcessFastaFile:- Adding %s..",pszFile);

bFirstEntry = true;
bEntryCreated = false;
SeqID = 0;
NumSampled = 0;
BuffOfs = 0;
NumSeqsUnderlength = 0;
NumSeqsOverlength = 0;
NumSeqsAccepted = 0;
TotAcceptedLen = 0;
while((Rslt = SeqLen = Fasta.ReadSequence(&pSeqBuff[BuffOfs],(int)min(AvailBuffSize,(size_t)cMaxAllocBuffChunk),true,false)) > eBSFSuccess)
	{
	if(SeqLen == eBSFFastaDescr)		// just read a descriptor line
		{
		SeqID++;
		if(bEntryCreated)				// add any previous entry
			{
			if((NumSeqsAccepted + NumSeqsUnderlength) <= ((SeqID * (UINT64)m_SampleInRate) / 100))
				{
				NumSampled += 1;
				if(BuffOfs < (size_t)MinSeqLen || BuffOfs > (size_t)MaxSeqLen)
					{
					if(BuffOfs < (size_t)MinSeqLen)
						NumSeqsUnderlength += 1;
					else
						NumSeqsOverlength += 1;
					}
				else
					{
					if((Rslt=m_pSfxArray->AddEntry(szName,pSeqBuff,(UINT32)BuffOfs,Flags)) < eBSFSuccess)
						{
						gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessFastaFile - error %d %s",Rslt,m_pSfxArray->GetErrMsg());
						break;
						}
					else
						{
						NumSeqsAccepted += 1;
						TotAcceptedLen += BuffOfs;
						}
					}
				}
			Rslt = eBSFSuccess;
			}
		Descrlen = Fasta.ReadDescriptor(szDescription,cBSFDescriptionSize);
		// An assumption - will one day bite real hard - is that the
		// fasta descriptor line starts with some form of unique identifier.
		// Use this identifier as the entry name.
		if(sscanf(szDescription," %s[ ,]",szName)!=1)
			sprintf(szName,"%s.%d",pszFile,++SeqID);

		bFirstEntry = false;
		bEntryCreated = true;
		BuffOfs = 0;
		continue;
		}
	else
		if(bFirstEntry)	// if there was no descriptor then dummy up one...
			{
			SeqID++;
			sprintf(szName,"%s.%d",pszFile,SeqID);
			strcpy(szDescription,"No Description provided");
			bFirstEntry = false;
			bEntryCreated = true;
			}

	// remove any repeat masking flags so that sorts can actually sort
	// if run of more than 25 Ns and at least 5 Ns to end of buffer then randomly mutate
	// every 13th N
	//	e.g <25Ns>r<12Ns>r<12Ns> where r is a pseudorandom base
	pMskBase = &pSeqBuff[BuffOfs];
	int SeqNs = 0;
	for(MskIdx = 0; MskIdx < SeqLen; MskIdx++,pMskBase++)
		{
		*pMskBase &= ~cRptMskFlg;
		if(*pMskBase == eBaseN && (MskIdx+5) < SeqLen)
			{
			if(++SeqNs > 25 &&
				pMskBase[1] == eBaseN &&
				pMskBase[2] == eBaseN &&
				pMskBase[3] == eBaseN &&
				pMskBase[4] == eBaseN)
				{
				if(!(SeqNs % 13))	// mutate every 13th
					*pMskBase = rand() % 4;
				}
			}
		else
			SeqNs = 0;
		}


	if(m_FiltMinHomoLen > 0)
		{
	// checking for near homopolymer runs and marking these bases for deletion
		UINT8 *pHomoBase;
		UINT8 *pDelHomoBase;
		UINT8 HomoBase;
		UINT8 CurHomoBase;
		UINT32 TotHomoLen;
		UINT32 HomoIdx;
		UINT32 StartHomoIdx;
		bool bHomoRuns;
		TotHomoLen = 0;
		StartHomoIdx = 0;
		bHomoRuns = false;
		pHomoBase = &pSeqBuff[BuffOfs];
		for(HomoIdx = 0; HomoIdx < SeqLen; HomoIdx++,pHomoBase++)
			{
			HomoBase = *pHomoBase & 0x03;
			if(TotHomoLen == 0)
				{
				StartHomoIdx = HomoIdx;
				CurHomoBase = HomoBase;
				TotHomoLen = 1;
				}
			else
				{
				if(HomoBase == CurHomoBase)
					TotHomoLen += 1;
				else
					{
					if(TotHomoLen >= 12)				// requiring an initial seed K-mer of at least 12 bases before accepting as possibly a homopolymer
						{
						// although not matching the current homopolymer base if the next three bases would match then accept as still being part of the homopolymer
						if((HomoIdx + 4) < SeqLen)
							{ 
							if(((pHomoBase[1] & 0x03) == CurHomoBase) && ((pHomoBase[2] & 0x03) == CurHomoBase) && ((pHomoBase[3] & 0x03) == CurHomoBase))
								{
								TotHomoLen += 1;
								continue;
								}
							}

						if(TotHomoLen >= (UINT32)m_FiltMinHomoLen)				// accepting as homopolymer if at least m_FiltMinHomoLen long
							{
							pDelHomoBase = &pSeqBuff[StartHomoIdx+6];	// retaining the first 6 bases as these started the homopolymer run
							TotHomoLen -= 6;							// and retaining the last 6 bases as these terminated the homopolymer run
							while(TotHomoLen--) 
								*pDelHomoBase++ = eBaseInDel;			// marking base for subsequent deletion 
							bHomoRuns = true;
							}
						}
					TotHomoLen = 0;
					StartHomoIdx = 0;
					}
				}
 			}
		if(TotHomoLen >= (UINT32)m_FiltMinHomoLen)			// accepting as homopolymer if at least m_FiltMinHomoLen long
			{
			pDelHomoBase = &pSeqBuff[StartHomoIdx+6];	// retaining the first 5 bases as these started the homopolymer run
			TotHomoLen -= 6;					    // and retaining the last 5 bases as these terminated the homopolymer run
			while(TotHomoLen--) 
				*pDelHomoBase++ = eBaseInDel;      // marking base for subsequent deletion 
			bHomoRuns = true;
			}

		if(bHomoRuns)
			{
			TotHomoLen = 0;
			pDelHomoBase = pHomoBase = &pSeqBuff[BuffOfs];
			for(HomoIdx = 0; HomoIdx < SeqLen; HomoIdx++,pHomoBase++)
				{
				HomoBase = *pHomoBase & 0x07;
				if(HomoBase == eBaseInDel)
					continue;
				*pDelHomoBase++ = *pHomoBase;
				TotHomoLen += 1;
				}
			SeqLen = TotHomoLen;
			}
		}

	BuffOfs += SeqLen;
	AvailBuffSize -= SeqLen;
	if(AvailBuffSize < (size_t)(cMaxAllocBuffChunk / 8))
		{
		size_t NewSize = (size_t)cMaxAllocBuffChunk + AllocdBuffSize;
		unsigned char *pTmp;
		if((pTmp = (unsigned char *)realloc(pSeqBuff,NewSize))==NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessFastaFile:- Unable to reallocate memory (%u bytes) for sequence buffer",(UINT32)NewSize);
			return(eBSFerrMem);
			}
		pSeqBuff = pTmp;
		AllocdBuffSize = NewSize;
		AvailBuffSize = AllocdBuffSize - BuffOfs;
		}
	}
if(Rslt < eBSFSuccess && Rslt != eBSErrSession)
	{
	while(Fasta.NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,Fasta.GetErrMsg());
	}


if(Rslt >= eBSFSuccess && bEntryCreated && BuffOfs > 0)			// close entry
	{
	if((NumSeqsAccepted + NumSeqsUnderlength) <= ((SeqID * (UINT64)m_SampleInRate) / 100))
		{
		NumSampled += 1;
		if(BuffOfs < (size_t)MinSeqLen || BuffOfs > (size_t)MaxSeqLen)
			{
			if (BuffOfs < (size_t)MinSeqLen)
				NumSeqsUnderlength += 1;
			else
				NumSeqsOverlength += 1;
			Rslt = eBSFSuccess;
			}
		else
			{
			if((Rslt=m_pSfxArray->AddEntry(szName,pSeqBuff,(UINT32)BuffOfs,Flags)) < eBSFSuccess)
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessFastaFile - error %d %s",Rslt,m_pSfxArray->GetErrMsg());
			else
				{
				Rslt = eBSFSuccess;
				NumSeqsAccepted += 1;
				TotAcceptedLen += BuffOfs;
				}
			}
		}
	}
if(pSeqBuff != NULL)
	free(pSeqBuff);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"ProcessFastaFile - %d parsed, %d sampled, %d accepted, %dbp mean length, %d sequences not accepted for indexing as length under %dbp, %d as length over %dbp",
					SeqID,NumSampled,NumSeqsAccepted,NumSeqsAccepted == 0 ? 0 :(int)(TotAcceptedLen/NumSeqsAccepted),NumSeqsUnderlength,MinSeqLen, NumSeqsOverlength, MaxSeqLen);
return(Rslt);
}

INT64
CPBErrCorrect::EstSumSeqLens(int NumTargFiles,char **pszTargFiles)		// guestimate total sequence length by simply summing the lengths of each file - likely to grossly over estimate
{
INT64 SumFileSizes;
int Idx;
int NumGlobs;

CSimpleGlob glob(SG_GLOB_FULLSORT);

	// determine crude estimate of sum of sequence lengths from the sum of file sizes
for(Idx = 0; Idx < NumTargFiles; Idx++)
	if (glob.Add(pszTargFiles[Idx]) < SG_SUCCESS)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to glob '%s",pszTargFiles[Idx]);
		Reset();
		return((INT64)eBSFerrOpnFile);
   		}


SumFileSizes = 0;
for (NumGlobs = 0; NumGlobs < glob.FileCount(); NumGlobs += 1)
	{
#ifdef _WIN32
	struct _stat64 st;
	if(!_stat64(glob.File(NumGlobs),&st))
#else
	struct stat64 st;
	if(!stat64(glob.File(NumGlobs),&st))
#endif
		SumFileSizes += (INT64)st.st_size;
	}
if(NumGlobs == 0)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to locate any input files");
	Reset();
	return(eBSFerrOpnFile);
	}

return(SumFileSizes);
}


int
CPBErrCorrect::GenConsensusFromMAF(int MinErrCorrectLen,		// error corrected sequences must be at least this minimum length
			 int MinConcScore,			// error corrected sequences trimmed until mean 50bp concensus score is at least this threshold
			char *pszErrCorFile,		// name of file into which write error corrected sequences
			char *pszMultiAlignFile)	// name of file containing multiple alignments to process
{
int Rslt;
CSSW SSW;
Rslt = SSW.GenConsensusFromMAF(MinErrCorrectLen,MinConcScore,pszErrCorFile,pszMultiAlignFile);
SSW.Reset();
return(Rslt);
}

int 
CPBErrCorrect::LoadSeqs(int MinSeqLen,					 // only accept for indexing sequences of at least this length
						int MaxSeqLen,                  // and no longer than this length (bp)
				int NumTargFiles,char **pszTargFiles,		// parse, and index sequences in these files into in memory suffix array; file expected to contain either fasta or fastq sequences
				int Flags)			// which by default are low confidence PacBio read sequences
{
int Rslt;
int Idx;
int NumGlobs;
INT64 SumFileSizes;
UINT32 DupEntries[cMaxDupEntries];
int NumDupEntries;
char szDupEntry[100];


CSimpleGlob glob(SG_GLOB_FULLSORT);

for(Idx = 0; Idx < NumTargFiles; Idx++)
	if (glob.Add(pszTargFiles[Idx]) < SG_SUCCESS)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to glob '%s",pszTargFiles[Idx]);
		m_pSfxArray->Close(false);
		Reset();
		return(eBSFerrOpnFile);
   		}


SumFileSizes = 0;
for (NumGlobs = 0; NumGlobs < glob.FileCount(); NumGlobs += 1)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Will process in this order: %d '%s", NumGlobs+1,glob.File(NumGlobs));
#ifdef _WIN32
	struct _stat64 st;
	if(!_stat64(glob.File(NumGlobs),&st))
#else
	struct stat64 st;
	if(!stat64(glob.File(NumGlobs),&st))
#endif
		SumFileSizes += (INT64)st.st_size;
	}
if(NumGlobs == 0)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to locate any input files");
	m_pSfxArray->Close(false);
	Reset();
	return(eBSFerrOpnFile);
	}

Rslt = eBSFSuccess;

for (int n = 0; Rslt >= eBSFSuccess &&  n < glob.FileCount(); ++n)
	{
		// try opening as a fasta file, if that fails then try as a bioseq
	Rslt = ProcessFastaFile(MinSeqLen, MaxSeqLen, glob.File(n),Flags);
	if(Rslt == eBSFerrNotFasta)
		Rslt = ProcessBioseqFile(MinSeqLen, MaxSeqLen,glob.File(n),Flags);
	if(Rslt < eBSFSuccess)
		{
		m_pSfxArray->Close(false);
		Reset();
		return(Rslt);
		}

		// check for duplicate entry names
	if((NumDupEntries = m_pSfxArray->ChkDupEntries(cMaxDupEntries,&DupEntries[0])) > 0)
		{
		while(NumDupEntries--)
			{
			m_pSfxArray->GetIdentName(DupEntries[NumDupEntries],sizeof(szDupEntry)-1,szDupEntry); // get sequence name for specified entry identifier
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"CreateBioseqSuffixFile, duplicate sequence entry name '%s' in file '%s'",szDupEntry,glob.File(n));
			}
		m_pSfxArray->Close(false);
		Reset();
		return(eBSFerrFastqSeqID);
		}
	}
return(eBSFSuccess);

}

int
CPBErrCorrect::Process(etPBPMode PMode,		// processing mode
		char *pszHostName,			// listening on this host name or IPv4/IPv5 address for connections by service providers 
		char *pszServiceName,			// Listen on this service name or port for for connections by service providers
		int MaxRMI,					// max number of RMI service provider instances supported
		int MaxNonRMI,				// max number of non-RMI SW threads supported
		int SampleInRate,			// sample input sequences at this rate per 100 (1..100)
		int SampleAcceptRate,		// sample accepted input sequences at this rate per 1000 (1..1000)
		int FiltMinHomoLen,			// filter PacBio reads for homopolymer runs >= this length (0 to disable filtering) 
		bool bSenseOnlyOvlps,		// process for sense only overlaps (default is for sense/sense and sense/antisense overlaps)
		int DeltaCoreOfs,			// offset by this many bp the core windows of coreSeqLen along the probe sequence when checking for overlaps
		int MaxSeedCoreDepth,		// only further process a seed core if there are no more than this number of matching cores in all targeted sequences
		int MinSeedCoreLen,			// use seed cores of this length when identifying putative overlapping scaffold sequences
		int MinNumSeedCores,        // require at least this many seed cores between overlapping scaffold sequences
		int SWMatchScore,			// score for matching bases (0..50)
		int SWMismatchPenalty,		// mismatch penalty (-50..0)
		int SWGapOpenPenalty,		// gap opening penalty (-50..0)
		int SWGapExtnPenalty,		// gap extension penalty (-50..0)
		int SWProgExtnPenaltyLen,	// progressive gap scoring then only apply gap extension score if gap at least this length (0..63) - use if aligning PacBio
		int TranscriptomeLens,		// 0 if disabled, processing transcript reads, putatively overlapping reads must have length differential no more than this percentage and overlaps to be nearly full length
		int MinPBSeqLen,			// only accepting PacBio reads to be error corrected of at least this length (defaults to 10Kbp) and if 
		int MaxPBSeqLen,			// no more than this length (defaults to 35Kbp)
		int MinPBSeqOverlap,		// any overlap of a PacBio onto a target PacBio must be of at least this many bp to be considered for contributing towards error correction (defaults to 5Kbp) 
		int MaxArtefactDev,			// classify overlaps as artefactual if sliding window of 500bp over any overlap deviates by more than this percentage from the overlap mean
		int MinHCSeqLen,			// only accepting hiconfidence reads of at least this length (defaults to 1Kbp)
		int MinHCSeqOverlap,		// any overlap of a hiconfidence read onto a target PacBio read must be of at least this many bp to be considered for contributing towards error correction (defaults to 1Kbp) 
		int HCRelWeighting,             // hiconfidence read overlaps are usually weighted higher than normal lesser confidence read overlaps when calling consensus bases 
	    int MinErrCorrectLen,		// error corrected and trimmed sequences must be at least this minimum length
		int MinConcScore,			// error corrected sequences trimmed until mean 50bp concensus score is at least this threshold
		int NumPacBioFiles,			// number of input pacbio file specs
		char *pszPacBioFiles[],		// input pacbio files
		int NumHiConfFiles,			// number of input hiconfidence file specs
		char *pszHiConfFiles[],		// input hiconfidence files		
	    char *pszErrCorFile,		// name of file into which write error corrected sequences
		char *pszMultiAlignFile,	// name of file into which write multiple alignments
		char *pszChkPtsFile,        // name of file used for checkpointing in case resume processing is required
		int NumThreads)				// maximum number of worker threads to use
{
int Rslt = eBSFSuccess;
int Idx;
UINT32 NumTargSeqs;
UINT32 CurNodeID;
UINT32 MaxSeqLen;
tsPBEScaffNode *pCurPBScaffNode;

Reset();

CreateMutexes();

bool bRslt;
if(pszHostName[0] != '\0')
	{
	m_bRMI = true;
	strncpy(m_szRMIHostName,pszHostName,sizeof(m_szRMIHostName));
	m_szRMIHostName[sizeof(m_szRMIHostName)-1] = '\0';
	strncpy(m_szRMIServiceName,pszServiceName,sizeof(m_szRMIServiceName));
	m_szRMIServiceName[sizeof(m_szRMIServiceName)-1] = '\0';
    m_MaxRMIInstances = MaxRMI;
	m_MaxNonRMIThreads = MaxNonRMI;
	if((m_pRequester = new CBKSRequester)==NULL)
		{
		Reset();
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: Unable to instantiate CBKSRequester");
		return(eBSFerrObj);
		}

	if((Rslt = m_pRequester->Initialise(pszHostName, pszServiceName))!=eBSFSuccess)
		{
		Reset();
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: Unable to initialise CBKSRequester instance");
		return(cBSFNWSProtErr);
		}

	if((bRslt = m_pRequester->Run())!=true)
		{
		Reset();
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: Unable to run CBKSRequester");
		return(cBSFNWSProtErr);
		}
	}
else
	{
	m_bRMI = false;
	m_MaxNonRMIThreads = NumThreads;
	m_MaxRMIInstances = 0;
	m_pRequester = NULL;
	m_szRMIHostName[0] = '\0';
	m_szRMIServiceName[0] = '\0';
	}

m_SWMatchScore = SWMatchScore;
m_SWMismatchPenalty = SWMismatchPenalty;
m_SWGapOpenPenalty = SWGapOpenPenalty;
m_SWGapExtnPenalty = SWGapExtnPenalty;
m_SWProgExtnPenaltyLen = SWProgExtnPenaltyLen;

m_PMode = PMode;
m_DeltaCoreOfs = DeltaCoreOfs;
m_MaxSeedCoreDepth = MaxSeedCoreDepth;
m_MinSeedCoreLen = MinSeedCoreLen;
m_MinNumSeedCores = MinNumSeedCores;

m_MinPBSeqLen = MinPBSeqLen;	
m_MaxPBRdSeqLen = MaxPBSeqLen;

m_MinPBSeqOverlap = MinPBSeqOverlap; 
m_MaxArtefactDev = MaxArtefactDev;
m_MinHCSeqLen = MinHCSeqLen;
m_MinHCSeqOverlap = MinHCSeqOverlap;
m_HCRelWeighting = HCRelWeighting;  

m_MinErrCorrectLen = MinErrCorrectLen;
m_MinConcScore = MinConcScore;
m_SampleInRate = SampleInRate;
m_SampleAcceptRate = SampleAcceptRate;
m_bSenseOnlyOvlps = bSenseOnlyOvlps;
m_FiltMinHomoLen = 0;
if(FiltMinHomoLen > 0)
	FiltMinHomoLen = max(FiltMinHomoLen,12);

m_TranscriptomeLens = TranscriptomeLens;

switch(PMode) {
	case ePBPMOverlapDetail:
		m_OverlapFloat = cDfltScaffMaxOverlapFloat;
		break;

	case ePBMConsolidate:
		m_OverlapFloat = cDfltConsolidateMaxOverlapFloat;
		break;

	case ePBPMConsensus:
		Rslt = GenConsensusFromMAF(MinErrCorrectLen,MinConcScore,pszErrCorFile,pszPacBioFiles[0]);
		Reset();
		return(Rslt);

	default:  // error correcting PacBio reads
		m_OverlapFloat = cDfltMaxOverlapFloat;
		m_FiltMinHomoLen = FiltMinHomoLen;
		break;
	}

m_NumPacBioFiles = NumPacBioFiles;
memset(m_szPacBioFiles,0,sizeof(m_szPacBioFiles));
if(NumPacBioFiles > 0)
	{
	for(Idx = 0; Idx < NumPacBioFiles; Idx++)
		strcpy(m_szPacBioFiles[Idx],pszPacBioFiles[Idx]);
	}
m_NumHiConfFiles = NumHiConfFiles;
memset(m_szHiConfFiles,0,sizeof(m_szHiConfFiles));
if(NumHiConfFiles > 0)
	{
	for(Idx = 0; Idx < NumHiConfFiles; Idx++)
		strcpy(m_szHiConfFiles[Idx],pszHiConfFiles[Idx]);
	}	

if(pszErrCorFile != NULL && pszErrCorFile[0] != '\0')
	{
	strncpy(m_szErrCorFile,pszErrCorFile,sizeof(m_szErrCorFile));
	m_szErrCorFile[sizeof(m_szErrCorFile)-1] = '\0';
	}
else
	m_szErrCorFile[0] = '\0';
m_hErrCorFile = -1;
m_ErrCorFileUnsyncedSize = 0;

if(pszChkPtsFile != NULL && pszChkPtsFile[0] != '\0')
	{
	strncpy(m_szChkPtsFile,pszChkPtsFile,sizeof(m_szChkPtsFile));
	m_szChkPtsFile[sizeof(m_szChkPtsFile)-1] = '\0';
	}
else
	m_szChkPtsFile[0] = '\0';
m_hChkPtsFile = -1;

if(pszMultiAlignFile != NULL && pszMultiAlignFile[0] != '\0')
	{
	strncpy(m_szMultiAlignFile,pszMultiAlignFile,sizeof(m_szMultiAlignFile));
	m_szMultiAlignFile[sizeof(m_szMultiAlignFile)-1] = '\0';
	}	
else
	m_szMultiAlignFile[0] = '\0';
m_hMultiAlignFile = -1;
m_MultiAlignFileUnsyncedSize = 0;

INT64 SumFileSizes;

SumFileSizes = EstSumSeqLens(NumPacBioFiles,pszPacBioFiles);	// guestimate sum of all sequences; will return -1 if errors
if(SumFileSizes == -1)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: Unable to access the PacBio files");
	return(eBSFerrObj);
	}

if(NumHiConfFiles && pszHiConfFiles != NULL)
	{
	INT64 HiConfFileSizes;
	HiConfFileSizes = EstSumSeqLens(NumHiConfFiles,pszHiConfFiles);	// guestimate sum of all sequences; will return eBSFerrOpnFile if errors
	if(HiConfFileSizes == (INT64)eBSFerrOpnFile)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: Unable to access the high confidence files");
		return(eBSFerrObj);
		}
	SumFileSizes += HiConfFileSizes;
	}

if(PMode == ePBMConsolidate)
	m_NumOvlpCores = 1;
else
	m_NumOvlpCores = NumThreads;	
if(m_pSfxArray != NULL)
	delete m_pSfxArray;
if((m_pSfxArray = new CSfxArrayV3) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadTargetSeqs: Unable to instantiate instance of CSfxArrayV3");
	return(eBSFerrObj);
	}
m_pSfxArray->Reset(false);
m_pSfxArray->SetMaxQSortThreads(min(NumThreads,32));


m_pSfxArray->SetMaxBaseCmpLen(200);     // matches are only down at relatively small K-mer cores so no point in sorting the reads out to their full length 

Rslt=m_pSfxArray->Open(false,false);
if(Rslt !=eBSFSuccess)
	{
	while(m_pSfxArray->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,m_pSfxArray->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"CreateBioseqSuffixFile, unable to create in-memory suffix array - error %d %s",Rslt,m_pSfxArray->GetErrMsg());
	Reset();
	return(Rslt);
	}

if((Rslt=m_pSfxArray->SetDescription((char *)"inmem")) != eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"CreateBioseqSuffixFile: Unable to set description 'inmem'");
	Reset();
	return(Rslt);
	}
if((Rslt=m_pSfxArray->SetTitle((char *)"inmem")) != eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"CreateBioseqSuffixFile: Unable to set title 'inmem'");
	Reset();
	return(Rslt);
	}

if((Rslt = m_pSfxArray->SetDatasetName((char *)"inmem")) != eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"CreateBioseqSuffixFile: Unable to set dataset name 'inmem'");
	Reset();
	return(Rslt);
	}
m_pSfxArray->SetInitalSfxAllocEls(SumFileSizes);	// just a hint which is used for initial allocations by suffix processing


if((Rslt = LoadSeqs(m_MinPBSeqOverlap,m_MaxPBRdSeqLen, NumPacBioFiles,pszPacBioFiles,cFlgLCSeq)) < eBSFSuccess)
	{
	Reset();
	return(Rslt);
	}
	// get number of sequences BioSeqs accepted for error correction
if((int)(NumTargSeqs = m_pSfxArray->GetNumEntries()) < 1)
	NumTargSeqs = 0;
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Loaded and accepted for processing a total of %d PacBio sequences",NumTargSeqs);

if(NumTargSeqs < 3)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Need at least 3 PacBio sequences for error correction accepted from file(s)");
	m_pSfxArray->Close(false);
	Reset();
	return(eBSFerrNoEntries);
	}

m_NumAcceptTargSeqs = NumTargSeqs;

int NumHCSeqs = 0;
if(NumHiConfFiles && pszHiConfFiles != NULL)			// load any optional high confidence sequences requested by user
	{
	m_FiltMinHomoLen = 0;
	if((Rslt = LoadSeqs(m_MinHCSeqLen, m_MaxPBRdSeqLen, NumHiConfFiles,pszHiConfFiles,cFlgHCSeq)) < eBSFSuccess)
		{
		m_pSfxArray->Close(false);
		Reset();
		return(Rslt);
		}
		// get number of high confidence sequences accepted
	NumHCSeqs = m_pSfxArray->GetNumEntries() - NumTargSeqs;
	if(NumHCSeqs < 1)
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"No high confidence sequences for error correction accepted from file(s)");
	else
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Loaded and accepted for processing a total of %d high confidence sequences",NumHCSeqs);
	}
NumTargSeqs = m_pSfxArray->GetNumEntries();
if(NumHCSeqs > 0)
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Loaded and accepted for processing a total of %d sequences including high confidence sequences",NumTargSeqs);

gDiagnostics.DiagOut(eDLInfo,gszProcName,"CreateBioseqSuffixFile: sorting suffix array...");
m_pSfxArray->Finalise();
gDiagnostics.DiagOut(eDLInfo,gszProcName,"CreateBioseqSuffixFile: sorting completed");

if(MinSeedCoreLen <= cMaxKmerLen)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Initialising for over occurring ( > %d) K-mers of length %d",m_MaxSeedCoreDepth,MinSeedCoreLen);
	if((Rslt = m_pSfxArray->InitOverOccKMers((int)MinSeedCoreLen,m_MaxSeedCoreDepth))!=eBSFSuccess)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Failed to initialise for over occurring K-mers");
		Reset();
		return(Rslt);
		}
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Initialised for over occurring K-mers");
	}

if((m_pPBScaffNodes = new tsPBEScaffNode [NumTargSeqs + 1]) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate for %d nodes",NumTargSeqs);
	Reset();
	return(eBSFerrMem);
	}
memset(m_pPBScaffNodes,0,sizeof(tsPBEScaffNode) * (NumTargSeqs+1));
m_AllocdPBScaffNodes = NumTargSeqs;
if((m_pMapEntryID2NodeIDs = new UINT32 [NumTargSeqs + 1]) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate for %d mapping entries nodes",NumTargSeqs);
	Reset();
	return(eBSFerrMem);
	}
memset(m_pMapEntryID2NodeIDs,0,sizeof(UINT32) * (NumTargSeqs+1));
m_NumPBScaffNodes = 0;

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Initialising for %d sequences",NumTargSeqs);	
MaxSeqLen = 0;
UINT32 NumSloughed;
UINT32 NumAccepted;
double PropAccepted;
double ExpPropAccept;

pCurPBScaffNode = m_pPBScaffNodes;
NumSloughed = 0;
NumAccepted = 0;
ExpPropAccept = (double)m_SampleAcceptRate / 1000.0;
for(CurNodeID = 1; CurNodeID <= NumTargSeqs; CurNodeID++,pCurPBScaffNode++)
	{
	UINT16 SeqFlags;
	pCurPBScaffNode->SeqLen = m_pSfxArray->GetSeqLen(CurNodeID);
	pCurPBScaffNode->EntryID = CurNodeID;
	SeqFlags = m_pSfxArray->GetIdentFlags(CurNodeID);
	pCurPBScaffNode->flgHCseq = SeqFlags & cFlgHCSeq ? 1 : 0;
	pCurPBScaffNode->flgUnderlength = pCurPBScaffNode->SeqLen < m_MinPBSeqLen ? 1 : 0;
	pCurPBScaffNode->flgSlough = 0;
	if(m_SampleAcceptRate < 1000 && !pCurPBScaffNode->flgHCseq && !pCurPBScaffNode->flgUnderlength)
		{
		if(NumAccepted > 0)
			PropAccepted = (double)NumAccepted / (double)(NumAccepted + NumSloughed);
		else
			PropAccepted = 0;
		if(PropAccepted < ExpPropAccept)
			NumAccepted += 1;
		else
			{
			NumSloughed += 1;
			pCurPBScaffNode->flgSlough = 1;
			}
		}

	if(MaxSeqLen == 0 || pCurPBScaffNode->SeqLen > (UINT32)MaxSeqLen)
		MaxSeqLen = pCurPBScaffNode->SeqLen;
	}
if(m_SampleAcceptRate < 1000)
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Sampling accepted read rate: %d per 1000, out of %d accepted reads for error correction will be processing %d",m_SampleAcceptRate,NumAccepted+NumSloughed, NumAccepted);

m_NumPBScaffNodes = NumTargSeqs;
m_MaxPBSeqLen = MaxSeqLen;

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Sorting %d sequences (max length %d) ...",NumTargSeqs, m_MaxPBSeqLen);

if(m_NumPBScaffNodes > 1)
	{
	// sort scaffold nodes by sequence length descending
	m_mtqsort.SetMaxThreads(min(NumThreads,16));
	m_mtqsort.qsort(m_pPBScaffNodes,m_NumPBScaffNodes,sizeof(tsPBEScaffNode),SortLenDescending);
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Sorting %d sequences completed",NumTargSeqs);
pCurPBScaffNode = m_pPBScaffNodes;
for(CurNodeID = 1; CurNodeID <= NumTargSeqs; CurNodeID++,pCurPBScaffNode++)
	{
	pCurPBScaffNode->NodeID = CurNodeID;
	m_pMapEntryID2NodeIDs[pCurPBScaffNode->EntryID-1] = CurNodeID;
	}
m_LowestCpltdProcNodeID = 0;
m_hChkPtsFile = -1;
m_hErrCorFile = -1;

// if checkpoints (syncpoints) points processing
// checkpoints file doesn't exist or contains too few entries ( less than 1000) to be worth the effort of parsing etc., then truncate and error correct reads from the beginning again
if(m_szChkPtsFile[0] != '\0' && m_szErrCorFile[0] != '\0')
	{
	tsECChkPt CurChkPt;
	UINT32 NodeID;
	UINT32 NumChkPtsExptd;
	UINT32 NumChkPtsAccepted;
	int ChkPtsStatRslt;
	off_t ChkPtFileSize;
	int ECStatRslt;
	off_t ECFileSize;

#ifdef _WIN32
	struct _stat64 ChkPtsStat;
	ChkPtsStatRslt = _stat64(m_szChkPtsFile, &ChkPtsStat);		// checking if exists, is a file, and if likely to contain sufficient entries (at least 5) to actually resume or if easier to restart
	struct _stat64 ECStat;
	ECStatRslt = _stat64(m_szErrCorFile, &ECStat);		// checking if exists, is a file, and if likely to contain sufficient entries (at least 5) to actually resume or if easier to restart
#else
	struct stat64 ChkPtsStat;
	ChkPtsStatRslt = stat64(m_szChkPtsFile, &ChkPtsStat);		// checking if exists, is a file, and if likely to contain sufficient entries (at least 5) to actually resume or if easier to restart
	struct stat64 ECStat;
	ECStatRslt = stat64(m_szErrCorFile, &ECStat);		// checking if exists, is a file, and if likely to contain sufficient entries (at least 5) to actually resume or if easier to restart
#endif
	// both checkpoint and error corrected reads files must exist and not be empty
	if(ChkPtsStatRslt >= 0 && ECStatRslt >= 0)
		{
		if(!(ChkPtsStat.st_mode & S_IFREG))
			{
			ChkPtsStatRslt = -1;
			ECStatRslt = -1;
			gDiagnostics.DiagOut(eDLWarn,gszProcName,"Checkpoint exists but not a regular file, restart error correction and checkpointing from 1st read");
			}
		else
			{
			ChkPtFileSize = (off_t)ChkPtsStat.st_size;
			if(ChkPtFileSize < sizeof(tsECChkPt) * 5)
				{
				ChkPtsStatRslt = -1;
				ECStatRslt =-1;
				gDiagnostics.DiagOut(eDLWarn,gszProcName,"Checkpoint exists but less than 5 checkpoints, restart error correction and checkpointing from 1st read");
				}
			}
		}
	else
		{
		ChkPtsStatRslt = -1;
		ECStatRslt =-1;
		}

	if(ECStatRslt >= 0)
		{
		if(!(ECStat.st_mode & S_IFREG))
			{
			ECStatRslt = -1;
			ChkPtsStatRslt = -1;
			gDiagnostics.DiagOut(eDLWarn,gszProcName,"Error corrected reads file exists but not a regular file, restart error correction and checkpointing from 1st read");
			}
		else
			{
			ECFileSize = (off_t)ECStat.st_size;
			if(ECFileSize < 100)			// purely arbitrary!
				{
				ECStatRslt = -1;
				ChkPtsStatRslt = -1;
				gDiagnostics.DiagOut(eDLWarn,gszProcName,"Checkpoint exists but empty or too small, restart error correction and checkpointing from 1st read");
				}
			}
		}


	if(ChkPtsStatRslt >= 0)
		{
#ifdef _WIN32
        m_hChkPtsFile = open(m_szChkPtsFile, _O_BINARY | _O_RDWR | _O_SEQUENTIAL, _S_IREAD | _S_IWRITE);
#else
        m_hChkPtsFile = open64(m_szChkPtsFile,O_RDWR,S_IREAD | S_IWRITE);
#endif
		if(m_hChkPtsFile < 0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to open existing checkpoint file '%s'",m_szChkPtsFile);
			Reset();
			return(eBSFerrOpnFile);
			}

#ifdef _WIN32
		m_hErrCorFile = open(m_szErrCorFile,_O_BINARY | _O_RDWR | _O_SEQUENTIAL, _S_IREAD | _S_IWRITE);
#else
		m_hErrCorFile = open(m_szErrCorFile,O_RDWR,S_IREAD | S_IWRITE);
#endif
		if(m_hErrCorFile < 0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to open existing error corrected reads file '%s'",m_szErrCorFile);
			Reset();
			return(eBSFerrOpnFile);
			}
		}

	if(ChkPtsStatRslt >= 0)
		{
		bool bRestartChkPts = false;
		UINT32 NumOverlapProcessed = 0;
		INT64 MaxECFileSize = 0;
		INT64 MaxECFileSizeAccepted = 0;
		// iterate over all checkpoint entries
		// stop iterating if either an inconsistency between checkpointed entry and local, or if offset is after error corrected file size
		NumChkPtsExptd = (UINT32)(ChkPtFileSize / sizeof(tsECChkPt));
		NumChkPtsAccepted = 0;
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Checkpoint file contains approximately %d entries, now processing and validating these entries",NumChkPtsExptd);
		while(NumChkPtsExptd && read(m_hChkPtsFile,&CurChkPt,sizeof(tsECChkPt)) == sizeof(tsECChkPt))
			{
			NodeID = CurChkPt.NodeID;
			if(NodeID == 0 || NodeID > m_NumPBScaffNodes)	// any inconsistencies are treated as if checkpoint file was originally empty
				{
				if(NumChkPtsAccepted >= 200 && NumChkPtsExptd < 100)		// if accepted at least 200 check points and less than 100 more check points remaining when an inconsistency is identified then accept checkpoints up to and including last accepted
					{
					_lseeki64(m_hChkPtsFile,(off_t)((UINT64)NumChkPtsAccepted *  sizeof(tsECChkPt)),SEEK_SET);
					break;					
					}
				bRestartChkPts = true;
				gDiagnostics.DiagOut(eDLWarn,gszProcName,"Checkpoint file node inconsistency at node %u, restart error correction and checkpointing from 1st read",NodeID);
				break;
				}

			pCurPBScaffNode = &m_pPBScaffNodes[CurChkPt.NodeID-1];
			if(pCurPBScaffNode->EntryID != CurChkPt.EntryID ||
				pCurPBScaffNode->flgUnderlength != CurChkPt.flgUnderlength ||
				pCurPBScaffNode->SeqLen != CurChkPt.SeqLen ||
				pCurPBScaffNode->flgHCseq != CurChkPt.flgHCseq)
				{
				if(NumChkPtsAccepted >= 200 && NumChkPtsExptd < 100)		// if accepted at least 200 check points and less than 100 more check points remaining when an inconsistency is identified then accept checkpoints up to and including last accepted
					{
					_lseeki64(m_hChkPtsFile,(off_t)((UINT64)NumChkPtsAccepted *  sizeof(tsECChkPt)),SEEK_SET);
					break;					
					}
				bRestartChkPts = true;
				gDiagnostics.DiagOut(eDLWarn,gszProcName,"Checkpoint file node inconsistency at node %u with read length and/or flags, restart error correction and checkpointing from 1st read",NodeID);
				break;
				}

			if(CurChkPt.ECFileOfs > 0)  // will be < 0 if no ecread generated for this sequence
				{
				if(CurChkPt.ECFileOfs > ECFileSize)
					{
					_lseeki64(m_hChkPtsFile,(off_t)((UINT64)NumChkPtsAccepted *  sizeof(tsECChkPt)),SEEK_SET);
					break;
					}
				if(CurChkPt.ECFileOfs > MaxECFileSize)
					{
					if((MaxECFileSize = (INT64)_lseeki64(m_hErrCorFile,(off_t)CurChkPt.ECFileOfs,SEEK_SET))!=CurChkPt.ECFileOfs)
						{
						if(NumChkPtsAccepted >= 200 && NumChkPtsExptd < 100)		// if accepted at least 200 check points and less than 100 more check points remaining when an inconsistency is identified then accept checkpoints up to and including last accepted
							{
							_lseeki64(m_hChkPtsFile,(off_t)((UINT64)NumChkPtsAccepted *  sizeof(tsECChkPt)),SEEK_SET);
							_lseeki64(m_hErrCorFile,(off_t)MaxECFileSizeAccepted,SEEK_SET);
							break;					
							}
						bRestartChkPts = true;
						gDiagnostics.DiagOut(eDLWarn,gszProcName,"Checkpoint file node inconsistency at node %u, unable to access error corrected reads file at offset %lld, restart error correction and checkpointing from 1st read",NodeID,CurChkPt.ECFileOfs);
						break;
						}
					}
				}


			pCurPBScaffNode->flgCpltdProc = 1;
			pCurPBScaffNode->flgContained = CurChkPt.flgContained;
			pCurPBScaffNode->flgContains = CurChkPt.flgContains;
			if(CurChkPt.NodeID > NumOverlapProcessed)
				NumOverlapProcessed = CurChkPt.NodeID;
			NumChkPtsExptd -= 1;
			NumChkPtsAccepted += 1;
			MaxECFileSizeAccepted = MaxECFileSize;

			pCurPBScaffNode = &m_pPBScaffNodes[m_LowestCpltdProcNodeID];
			for(CurNodeID = m_LowestCpltdProcNodeID + 1; CurNodeID <= m_NumPBScaffNodes; CurNodeID++, pCurPBScaffNode++)
				{
				if(pCurPBScaffNode->flgCpltdProc == 1 ||	// completed processing
						pCurPBScaffNode->flgHCseq	  ||	// not error correcting already assumed to be high confidence sequences 
						pCurPBScaffNode->flgUnderlength == 1 || // must be of a user specified minimum length 
						pCurPBScaffNode->SeqLen < m_MinPBSeqLen)
						continue;
				 m_LowestCpltdProcNodeID = CurNodeID - 1;
				break;
				 }
			if(CurNodeID == m_NumPBScaffNodes+1)
				m_LowestCpltdProcNodeID = CurNodeID - 1;
			}

		if(bRestartChkPts)
			{
			close(m_hChkPtsFile);
			close(m_hErrCorFile);
			m_hChkPtsFile = -1;
			m_hErrCorFile = -1;
			ChkPtsStatRslt = -1;
			ECStatRslt = -1;
			m_LowestCpltdProcNodeID = 0;

			// restore nodes back to all unprocessed
			pCurPBScaffNode = m_pPBScaffNodes;
			for(NumChkPtsExptd = 0; NumChkPtsExptd < m_NumPBScaffNodes; NumChkPtsExptd += 1, pCurPBScaffNode++)
				{
				pCurPBScaffNode->flgCpltdProc = 0;
				pCurPBScaffNode->flgContained = 0;
				pCurPBScaffNode->flgContains = 0;
				}
			}
		else
			{
			if(m_LowestCpltdProcNodeID == NumOverlapProcessed && NumOverlapProcessed == m_NumPBScaffNodes)
				{
				gDiagnostics.DiagOut(eDLWarn,gszProcName,"Completed checkpoint file processing, nothing to do; already completed error correction of %u reads", NumOverlapProcessed);
				close(m_hChkPtsFile);
				close(m_hErrCorFile);
				m_hChkPtsFile = -1;
				m_hErrCorFile = -1;
				Reset();
				return(eBSFSuccess);
				}
			gDiagnostics.DiagOut(eDLWarn,gszProcName,"Completed checkpoint file processing, resuming error correction from entry %u with highest checkpointed entry %u", m_LowestCpltdProcNodeID, NumOverlapProcessed);
			}
		}

	if(ChkPtsStatRslt < 0)
		{
#ifdef _WIN32
		m_hChkPtsFile = open(m_szChkPtsFile,( O_WRONLY | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC),(_S_IREAD | _S_IWRITE));
#else
		if((m_hChkPtsFile = open(m_szChkPtsFile,O_WRONLY | O_CREAT,S_IREAD | S_IWRITE))!=-1)
			if(ftruncate(m_hChkPtsFile,0)!=0)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate %s - %s",m_szChkPtsFile,strerror(errno));
				Reset();
				return(eBSFerrCreateFile);
				}
#endif
		if(m_hChkPtsFile < 0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate checkpoints file '%s'",m_szChkPtsFile);
			Reset();
			return(eBSFerrCreateFile);
			}
		}
	
	}

if(m_szErrCorFile[0] != '\0' && m_hErrCorFile == -1)
	{
#ifdef _WIN32
	m_hErrCorFile = open(m_szErrCorFile,( O_WRONLY | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC),(_S_IREAD | _S_IWRITE));
#else
	if((m_hErrCorFile = open(m_szErrCorFile,O_WRONLY | O_CREAT,S_IREAD | S_IWRITE))!=-1)
		if(ftruncate(m_hErrCorFile,0)!=0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate %s - %s",m_szErrCorFile,strerror(errno));
			Reset();
			return(eBSFerrCreateFile);
			}
#endif
	if(m_hErrCorFile < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate output file '%s'",m_szErrCorFile);
		Reset();
		return(eBSFerrCreateFile);
		}
	}

if(m_PMode == ePBPMOverlapDetail)
	{
	m_ScaffLineBuffIdx=sprintf(m_szScaffLineBuff,"\"Class\",\"ProbeID\",\"ProbDescr\",\"TargID\",\"TargDescr\",\"SeedHits\",\"ProbeSense\",\"TargSense\",\"ProbeLen\",\"TargLen\",\"ProbeAlignLength\",\"TargAlignLength\",\"PeakScore\",\"FinalScore\",\"NumAlignedBases\",\"NumExactBases\",\"NumProbeInserts\",\"NumProbeInsertBases\",\"NumTargInserts\",\"NumTargInsertBases\",\"ProbeStartOfs\",\"TargStartOfs\",\"ProbeEndOfs\",\"TargEndOfs\",\"ProbeOfs5\",\"TargOfs5\",\"ProbeOfs3\",\"TargOfs3\"");
	CUtility::SafeWrite(m_hErrCorFile,m_szScaffLineBuff,m_ScaffLineBuffIdx);
	m_ScaffLineBuffIdx = 0;
	}

if(m_szMultiAlignFile[0] != '\0')
	{
#ifdef _WIN32
	m_hMultiAlignFile = open(m_szMultiAlignFile,( O_WRONLY | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC),(_S_IREAD | _S_IWRITE));
#else
	if((m_hMultiAlignFile = open(m_szMultiAlignFile,O_WRONLY | O_CREAT,S_IREAD | S_IWRITE))!=-1)
		if(ftruncate(m_hMultiAlignFile,0)!=0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate %s - %s",m_szMultiAlignFile,strerror(errno));
			Reset();
			return(eBSFerrCreateFile);
			}
#endif
	if(m_hMultiAlignFile < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate output file '%s'",m_szMultiAlignFile);
		Reset();
		return(eBSFerrCreateFile);
		}
	}
else
	m_hMultiAlignFile = -1;

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Identifying sequence overlaps ...");

UINT32 MaxWorkerThreads;
UINT32 ReqRMICores;
if(m_bRMI)
	{
	MaxWorkerThreads = m_MaxRMIInstances;
	ReqRMICores = (m_MaxRMIInstances+cRMIThreadsPerCore-1)/cRMIThreadsPerCore;
	if(ReqRMICores < m_NumOvlpCores)
		MaxWorkerThreads += min(m_NumOvlpCores - ReqRMICores, m_MaxNonRMIThreads);
	}
else
	MaxWorkerThreads = m_NumOvlpCores;
Rslt = IdentifySequenceOverlaps(MaxSeqLen,m_NumOvlpCores,MaxWorkerThreads);
#ifdef _DEBUG
#ifdef _WIN32
_ASSERTE( _CrtCheckMemory());
#endif
#endif
Reset();
return(Rslt);
}

#ifdef _WIN32
unsigned __stdcall PBErrCorrectThread(void * pThreadPars)
#else
void *PBErrCorrectThread(void * pThreadPars)
#endif
{
int Rslt;
tsThreadPBErrCorrect *pPars = (tsThreadPBErrCorrect *)pThreadPars;			// makes it easier not having to deal with casts!
CPBErrCorrect *pPBErrCorrect = (CPBErrCorrect *)pPars->pThis;

Rslt = pPBErrCorrect->ThreadPBErrCorrect(pPars);
pPars->Rslt = Rslt;
#ifdef _WIN32
_endthreadex(0)                                                                             ;
return(eBSFSuccess);
#else
pthread_exit(NULL);
#endif
}


int
CPBErrCorrect::UpdateProcessStats(void)	// determine current CPU utilisation by this process and numbers of commited and uncommited service provider classses
{
UINT32 NumClasses;
UINT32 NumCommitedClasses;
UINT32 NumUncommitedClasses;
UINT32 CurNodeID;
int DeltaThreads;					// > 0 to enable increase in non-RMI active threads, < 0 to enable decrease in active non-RMI threads 
tsPBEScaffNode *pCurPBScaffNode;
time_t Now;

UINT32 TargNonRMIThreads;
UINT32 TargRMIThreads;

if(m_ProcessStatsThen == 0)
	{
	Now = m_ProcessStatsThen = time(NULL);
#ifdef WIN32
	GetSystemTimes(&m_PrevIdleTime,&m_PrevKernelTime,&m_PrevUserTime);
#endif
	}
else
	Now = time(NULL);

DeltaThreads = 0;
NumClasses = 0;
NumCommitedClasses = 0;
NumUncommitedClasses = 0;

if(m_bRMI && (Now - m_ProcessStatsThen) >= 60)
	{
	NumClasses = m_pRequester->GetNumClassInstances(eBKSPTSmithWaterman,&NumCommitedClasses,&NumUncommitedClasses);
	UINT32 PercentIdle;

#ifdef WIN32
	UINT64 CurIdleNanoSecs;    
	UINT64 PrevIdleNanoSecs;
	UINT64 Norm1SecIdleNanoSecs;  // idle nanosecs per second over the preceding 60 secs

	FILETIME CurIdleTime;
	FILETIME CurKernelTime;
	FILETIME CurUserTime;
	GetSystemTimes(&CurIdleTime,&CurKernelTime,&CurUserTime);
	CurIdleNanoSecs = (UINT64)CurIdleTime.dwLowDateTime;
	CurIdleNanoSecs |= (UINT64)CurIdleTime.dwHighDateTime << 32;
	PrevIdleNanoSecs = (UINT64)m_PrevIdleTime.dwLowDateTime;
	PrevIdleNanoSecs |= (UINT64)m_PrevIdleTime.dwHighDateTime << 32;
	m_PrevIdleTime = CurIdleTime;
	Norm1SecIdleNanoSecs = (CurIdleNanoSecs - PrevIdleNanoSecs) / (UINT64)(Now - m_ProcessStatsThen);
	PercentIdle = min(100,(UINT32)(Norm1SecIdleNanoSecs / 1000000));
#else
	double Loadavgs[3];
	getloadavg(Loadavgs,3);

	PercentIdle = (UINT32)max(0.0,100.0 - (100.0 * (Loadavgs[0] / m_NumOvlpCores)));

#endif
	if(PercentIdle < 3)
		DeltaThreads = max(-3,-1 * (int)(3 - PercentIdle));
	else
		if(PercentIdle > 5)
			DeltaThreads = min(5,(int)(PercentIdle - 5));

	AcquireCASSerialise();

	m_RMINumCommitedClasses = NumCommitedClasses;
	m_RMINumUncommitedClasses = NumUncommitedClasses;

	TargRMIThreads = min(NumClasses,m_MaxRMIInstances);	
	TargNonRMIThreads = m_CurActiveNonRMIThreads;
	if(DeltaThreads < 0)		// needing to reduce non-RMI threads
		TargNonRMIThreads = (UINT32)max(0,(int)m_CurActiveNonRMIThreads + DeltaThreads);
	else
		if(DeltaThreads > 0)		// wanting to increase non-RMI threads
			TargNonRMIThreads = (UINT32)min(m_MaxNonRMIThreads,m_CurActiveNonRMIThreads + DeltaThreads); 
 
	if(TargNonRMIThreads < m_CurActiveNonRMIThreads)
		m_ReduceNonRMIThreads = min(m_CurActiveNonRMIThreads - TargNonRMIThreads,cRMIThreadsPerCore);
	else
		m_ReduceNonRMIThreads = 0;
	m_CurAllowedActiveOvlpThreads = min(TargRMIThreads + TargNonRMIThreads,m_MaxActiveOvlpThreads);
	}
else
	AcquireCASSerialise();


if((Now - m_ProcessStatsThen) >= 60)
	{
	if(m_PMode == ePBPMErrCorrect)
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Progress: %u processed, SW aligned: %u, Overlapping: %u, Overlapped: %u, Contained: %u, Artifact: %u",
							m_NumOverlapProcessed,m_ProvSWchecked,m_ProvOverlapping,m_ProvOverlapped,m_ProvContained,m_ProvArtefact);
	else
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Progress: %u processed",m_NumOverlapProcessed);
	if(m_bRMI)
		{
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Progress: Current number of RMI SW threads: %d, number of non_RMI SW threads: %d, allowing up to %d active SW threads",
							m_CurActiveRMIThreads,m_CurActiveNonRMIThreads,m_CurAllowedActiveOvlpThreads);
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Progress: %u RMI classes instantiated, %u RMI classes uninstantiated",m_RMINumCommitedClasses,m_RMINumUncommitedClasses);
		}
	m_ProcessStatsThen += 60;
	}
ReleaseCASSerialise();
AcquireCASLock();
pCurPBScaffNode = &m_pPBScaffNodes[m_LowestCpltdProcNodeID];
for(CurNodeID = m_LowestCpltdProcNodeID + 1; CurNodeID <= m_NumPBScaffNodes; CurNodeID++, pCurPBScaffNode++)
	{
	if(pCurPBScaffNode->flgCpltdProc == 1 ||	// completed processing
			pCurPBScaffNode->flgHCseq	  ||	// not error correcting already assumed to be high confidence sequences 
			pCurPBScaffNode->flgUnderlength == 1 || // must be of a user specified minimum length 
			pCurPBScaffNode->SeqLen < m_MinPBSeqLen)
			continue;
	 m_LowestCpltdProcNodeID = CurNodeID - 1;
	break;
	 }
if(CurNodeID == m_NumPBScaffNodes+1)
	m_LowestCpltdProcNodeID = CurNodeID - 1;
ReleaseCASLock();
return(0);
}




int
CPBErrCorrect::IdentifySequenceOverlaps(int MaxSeqLen,		// max length sequence to be overlapped
									int NumOvlpCores,		// targeting to maximise usage of this many cores
									int NumOvlpThreads,		// identify all read overlaps using at most this this many threads 
								    teBKSPType RMIBKSPType,	// workers are to request this service type
									UINT32 RMIBufferSize,	// each worker thread default allocates to process up to this much buffered data
									UINT32 RMIParamDataSize,// each worker thread default allocates to process up to this much parameter data
									UINT32 RMIReqDataSize,	// each worker thread default allocates to process up to this much request data
									UINT32 RMIRespDataSize)	// each worker thread default allocates to process up to this much response data
{
tsThreadPBErrCorrect *pThreadPutOvlps;
int ThreadIdx;
tsThreadPBErrCorrect *pThreadPar;
UINT32 NumCommitedClasses;
UINT32 NumUncommitedClasses;
int NumClasses;
memset(m_ExactKmerDists,0,sizeof(m_ExactKmerDists));
m_TotAlignSeqLen = 0;
m_TotAlignSeqs = 0;

pThreadPutOvlps = new tsThreadPBErrCorrect [NumOvlpThreads];

pThreadPar = pThreadPutOvlps;
for(ThreadIdx = 0; ThreadIdx < NumOvlpThreads; ThreadIdx++,pThreadPar++)
	{
	memset(pThreadPar,0,sizeof(tsThreadPBErrCorrect));
	pThreadPar->AllocdCoreHits = cAllocdNumCoreHits;
	pThreadPar->AllocdCoreHitsSize = cAllocdNumCoreHits * sizeof(tsPBECoreHit);
#ifdef _WIN32
	pThreadPar->pCoreHits = (tsPBECoreHit *)malloc(pThreadPar->AllocdCoreHitsSize);	
	if(pThreadPar->pCoreHits == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"IdentifySequenceOverlaps: Core hits memory allocation of %llu bytes - %s",pThreadPar->AllocdCoreHitsSize,strerror(errno));
		break;
		}
#else
	if((pThreadPar->pCoreHits = (tsPBECoreHit *)mmap(NULL,pThreadPar->AllocdCoreHitsSize, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0)) == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"IdentifySequenceOverlaps: Core hits memory allocation of %llu bytes - %s",pThreadPar->AllocdCoreHitsSize,strerror(errno));
		break;
		}
#endif

	pThreadPar->AllocdProbeSeqSize = max(cAllocdQuerySeqLen,MaxSeqLen);
	if((pThreadPar->pProbeSeq = new etSeqBase [pThreadPar->AllocdProbeSeqSize + 10])==NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"IdentifySequenceOverlaps: Core hits probe sequence memory allocation of %d bytes - %s",pThreadPar->AllocdProbeSeqSize,strerror(errno));
		break;
		}

	pThreadPar->AllocdTargSeqSize = pThreadPar->AllocdProbeSeqSize;
	if((pThreadPar->pTargSeq = new etSeqBase [pThreadPar->AllocdTargSeqSize + 10])==NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"IdentifySequenceOverlaps: Core hits target sequence memory allocation of %d bytes - %s",pThreadPar->AllocdTargSeqSize,strerror(errno));
		break;
		}

	if(m_szErrCorFile[0] != '\0')
		{
		pThreadPar->AllocdErrCorLineBuff = cMaxPacBioErrCorLen;
		if((pThreadPar->pszErrCorLineBuff = new char [cMaxPacBioErrCorLen]) == NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate %d memory for error corection buffering",cMaxPacBioErrCorLen);
			break;
			}
		}

	if(m_szMultiAlignFile[0] != '\0')
		{
		pThreadPar->AllocdMultiAlignLineBuff = cMaxPacBioMAFLen;
		if((pThreadPar->pszMultiAlignLineBuff = new char [cMaxPacBioMAFLen]) == NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate %d memory for multialignment buffering",cMaxPacBioMAFLen);
			break;
			}
		}

	if((pThreadPar->pmtqsort = new CMTqsort) == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"IdentifySequenceOverlaps: Core hits instantiation of CMTqsort failed");
		break;
		}
	pThreadPar->pmtqsort->SetMaxThreads(4);
	pThreadPar->bRevCpl = false;
	pThreadPar->MaxSeedCoreDepth = m_MaxSeedCoreDepth;
	pThreadPar->DeltaCoreOfs = m_DeltaCoreOfs;
	pThreadPar->CoreSeqLen = m_MinSeedCoreLen;
	pThreadPar->MinNumCores = m_MinNumSeedCores;
	pThreadPar->MaxAcceptHitsPerSeedCore = cDfltMaxAcceptHitsPerSeedCore;
	pThreadPar->MinPBSeqLen = m_MinPBSeqLen;
	pThreadPar->MinOverlapLen = m_MinPBSeqOverlap;

// RMI support
	pThreadPar->bRMI = m_bRMI;
	if(m_bRMI)
		{
		pThreadPar->ServiceType = RMIBKSPType;
		pThreadPar->pRequester = m_pRequester;	
		pThreadPar->RMIReqDataSize = RMIReqDataSize;		// pRMIReqData allocation size
		if((pThreadPar->pRMIReqData = (UINT8 *)malloc(RMIReqDataSize))==NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate %d memory for RMIReqData buffering",RMIReqDataSize);
			break;
			}

		pThreadPar->RMIParamDataSize = RMIParamDataSize;	// pRMIParamData allocation size
		if((pThreadPar->pRMIParamData = (UINT8 *)malloc(RMIParamDataSize))==NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate %d memory for RMIParamData buffering",RMIParamDataSize);
			break;
			}

		pThreadPar->RMIRespDataSize = RMIRespDataSize;		// pRMIRespData allocation size
		if((pThreadPar->pRMIRespData = (UINT8 *)malloc(RMIRespDataSize))==NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate %d memory for RMIRespDatabuffering",RMIRespDataSize);
			break;
			}
		pThreadPar->RMIBufferSize = RMIBufferSize;			// RMIBuffer allocation size
		if((pThreadPar->pRMIBuffer = (UINT8 *)malloc(RMIBufferSize))==NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate %d memory for RMIBuffer buffering",RMIBufferSize);
			break;
			}
		}
	}

if(ThreadIdx != NumOvlpThreads)	// any errors whilst allocating memory for core hits?
	{
	do {
		if(pThreadPar->pmtqsort != NULL)
			delete pThreadPar->pmtqsort;

		if(pThreadPar->pszMultiAlignLineBuff != NULL)
			delete pThreadPar->pszMultiAlignLineBuff;

		if(pThreadPar->pszErrCorLineBuff != NULL)
			delete pThreadPar->pszErrCorLineBuff;

		if(pThreadPar->pCoreHits != NULL)
			{
#ifdef _WIN32
			free(pThreadPar->pCoreHits);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
			if(pThreadPar->pCoreHits != MAP_FAILED)
				munmap(pThreadPar->pCoreHits,pThreadPar->AllocdCoreHitsSize);
#endif	
			}
		if(pThreadPar->pProbeSeq != NULL)
			delete pThreadPar->pProbeSeq;

		if(pThreadPar->pRMIReqData != NULL)
			free(pThreadPar->pRMIReqData);
		if(pThreadPar->pRMIParamData != NULL)
			free(pThreadPar->pRMIParamData);
		if(pThreadPar->pRMIRespData != NULL)
			free(pThreadPar->pRMIRespData);
		if(pThreadPar->pRMIBuffer != NULL)
			free(pThreadPar->pRMIBuffer);

		pThreadPar -= 1;
		ThreadIdx -= 1;
		}
	while(ThreadIdx >= 0);
	delete pThreadPutOvlps;
	Reset();
	return((INT64)eBSFerrMem);
	}

if(m_bRMI)
	{
	NumClasses = m_pRequester->GetNumClassInstances(eBKSPTSmithWaterman,&NumCommitedClasses,&NumUncommitedClasses);
	m_RMINumCommitedClasses = NumCommitedClasses;
	m_RMINumUncommitedClasses = NumUncommitedClasses;
	}
else
	{
	NumClasses= 0;
	m_RMINumCommitedClasses = 0;
	m_RMINumUncommitedClasses = 0;
	}

m_NumOvlpCores = NumOvlpCores;
m_ReduceNonRMIThreads = 0;
m_MaxActiveOvlpThreads = NumOvlpThreads;
m_CurAllowedActiveOvlpThreads = m_MaxNonRMIThreads;
m_CurActiveRMIThreads = 0;
m_CurActiveNonRMIThreads = 0;
UpdateProcessStats();


pThreadPar = pThreadPutOvlps;
for (ThreadIdx = 1; ThreadIdx <= NumOvlpThreads; ThreadIdx++, pThreadPar++)
	{
	pThreadPar->ThreadIdx = ThreadIdx;
	pThreadPar->pThis = this;
#ifdef _WIN32
	pThreadPar->threadHandle = (HANDLE)_beginthreadex(NULL, 0x0fffff, PBErrCorrectThread, pThreadPar, 0, &pThreadPar->threadID);
#else
	pThreadPar->threadRslt = pthread_create(&pThreadPar->threadID, NULL, PBErrCorrectThread, pThreadPar);
#endif
	}

// allow threads 30 seconds to startup
#ifdef _WIN32
Sleep(30000);
#else
sleep(30);
#endif
pThreadPar = pThreadPutOvlps;
for (ThreadIdx = 0; ThreadIdx < NumOvlpThreads; ThreadIdx++, pThreadPar++)
	{
#ifdef _WIN32
	while (WAIT_TIMEOUT == WaitForSingleObject(pThreadPar->threadHandle, 30000))
		{
		UpdateProcessStats();
		};
	CloseHandle(pThreadPar->threadHandle);
#else
	struct timespec ts;
	int JoinRlt;
	clock_gettime(CLOCK_REALTIME, &ts);
	ts.tv_sec += 30;
	while ((JoinRlt = pthread_timedjoin_np(pThreadPar->threadID, NULL, &ts)) != 0)
		{
		UpdateProcessStats();
		ts.tv_sec += 30;
		}
#endif
	}

pThreadPar = pThreadPutOvlps;
for(ThreadIdx = 0; ThreadIdx < NumOvlpThreads; ThreadIdx++,pThreadPar++)
	{
	if(pThreadPar->pmtqsort != NULL)
		delete pThreadPar->pmtqsort;

	if(pThreadPar->pCoreHits != NULL)
		{
#ifdef _WIN32
		free(pThreadPar->pCoreHits);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
		if(pThreadPar->pCoreHits != MAP_FAILED)
			munmap(pThreadPar->pCoreHits,pThreadPar->AllocdCoreHitsSize);
#endif	
		pThreadPar->pCoreHits = NULL;
		}
	if(pThreadPar->pProbeSeq != NULL)
		delete pThreadPar->pProbeSeq;
	if(pThreadPar->pTargSeq != NULL)
		delete pThreadPar->pTargSeq;
	if(pThreadPar->pszMultiAlignLineBuff != NULL)
		delete pThreadPar->pszMultiAlignLineBuff;
	if(pThreadPar->pszErrCorLineBuff != NULL)
		delete pThreadPar->pszErrCorLineBuff;

	if(!m_bRMI && pThreadPar->pSW != NULL)
		delete pThreadPar->pSW;

	if(pThreadPar->pRMIReqData != NULL)
		free(pThreadPar->pRMIReqData);
	if(pThreadPar->pRMIParamData != NULL)
		free(pThreadPar->pRMIParamData);
	if(pThreadPar->pRMIRespData != NULL)
		free(pThreadPar->pRMIRespData);
	if(pThreadPar->pRMIBuffer != NULL)
		free(pThreadPar->pRMIBuffer);
	}

delete pThreadPutOvlps;

if(m_PMode == ePBPMErrCorrect)
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Completed: %u processed, SW aligned: %u, Overlapping: %u, Overlapped: %u, Contained: %u, Artefact: %u",
							m_NumOverlapProcessed,m_ProvSWchecked,m_ProvOverlapping,m_ProvOverlapped,m_ProvContained,m_ProvArtefact);
else
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Completed: %u processed",m_NumOverlapProcessed);

if(m_hMultiAlignFile != -1)
	{
#ifdef _WIN32
	_commit(m_hMultiAlignFile);
#else
	fsync(m_hMultiAlignFile);
#endif
	close(m_hMultiAlignFile);
	m_hMultiAlignFile = -1;
	}

if(m_hErrCorFile != -1)
	{
#ifdef _WIN32
	_commit(m_hErrCorFile);
#else
	fsync(m_hErrCorFile);
#endif
	close(m_hErrCorFile);
	m_hErrCorFile = -1;
	}
return(0);
}

int
CPBErrCorrect::ThreadPBErrCorrect(tsThreadPBErrCorrect *pThreadPar)
{
UINT32 AdjOverlapFloat;
UINT32 CurTargCoreHitCnts;
tsPBEScaffNode *pTargNode;

UINT32 TargSeqLen;
UINT32 ProbeSeqLen;
UINT32 LongSAligns;
UINT32 LongAAligns;
tsSSWCell *pPeakMatchesCell;
tsSSWCell PeakMatchesCell;
#ifdef _PEAKSCOREACCEPT
tsSSWCell PeakScoreCell;
#endif
bool bTargSense;
UINT32 Idx;
UINT32 CurSummaryHitCnts;
UINT32 LowestSummaryHitCnts;
sPBECoreHitCnts *pLowestSummaryHitCnts;
UINT32 HitIdx;
tsPBECoreHit *pCoreHit;
tsPBECoreHit *pNxtCoreHit;
UINT32 CurTargHitOfs;
UINT32 CurProbeHitOfs;
UINT32 CurSEntryIDHits;
UINT32 CurAEntryIDHits;
UINT32	CurSTargStartOfs;
UINT32	CurSTargEndOfs;
UINT32	CurATargStartOfs;
UINT32	CurATargEndOfs;
UINT32	CurSProbeStartOfs;
UINT32	CurSProbeEndOfs;
UINT32	CurAProbeStartOfs;
UINT32	CurAProbeEndOfs;
UINT32 ProvOverlapping;
UINT32 ProvOverlapped;
UINT32 ProvContained;
UINT32 ProvArtefact;
UINT32 ProvSWchecked;
UINT32 MinOverlapLen;
UINT32 MinTranscriptOverlapLen;
sPBECoreHitCnts *pSummaryCnts;
UINT32 CurNodeID;
UINT32 LowestCpltdProcNodeID;
int NumInMultiAlignment;

char szTargSeqName[cMaxDatasetSpeciesChrom];
char szProbeSeqName[cMaxDatasetSpeciesChrom];

UINT64 ClassInstanceID;
UINT32 ClassMethodID;
bool bRMIInitialised;
bool bRMIRslt;
int iRMIRslt;
bool bNonRMIRslt;
int iNonRMIRslt;

tsECChkPt CurChkPt;

tsPBEScaffNode *pCurPBScaffNode;
UINT32 RMINumUncommitedClasses;

pThreadPar->pSW = NULL;
pCurPBScaffNode = NULL;
ClassInstanceID = 0;

///////////////////////////////////////////////////////////////////////////////////////////////////////////
RMIRestartThread:							// RMI threads with errors detected are restarted from here with a goto!
///////////////////////////////////////////////////////////////////////////////////////////////////////////
if(ClassInstanceID != 0)   // non-zero if must have been a RMI SW thread which is restarting 
	{
	RMI_delete(pThreadPar,cRMI_SecsTimeout,ClassInstanceID);
	if(pCurPBScaffNode != NULL)
		{
		AcquireCASLock();
		pCurPBScaffNode->flgCurProc = 0;
		pCurPBScaffNode->flgCpltdProc = 0;
		ReleaseCASLock();
		}
	AcquireCASSerialise();
	if(m_CurActiveRMIThreads > 0)
		m_CurActiveRMIThreads -= 1;
	ReleaseCASSerialise();
	}

if(pThreadPar->pSW != NULL)  // non-NULL if must have been a non-RMI SW thread which is restarting
	{
	delete pThreadPar->pSW ;
	pThreadPar->pSW = NULL;
	pCurPBScaffNode = NULL;
	AcquireCASSerialise();
	if(m_CurActiveNonRMIThreads > 0)
		m_CurActiveNonRMIThreads -= 1;
	ReleaseCASSerialise();
	}

ClassInstanceID = 0;
ClassMethodID = 0;
bRMIInitialised = false;
pThreadPar->bRMI = false;
 
// if debugging and only interested in sessions with a specific identifier then set this to the session identifier of interest!!
UINT32 ReqSessionID = 0;

while(1)
	{
	AcquireCASThreadPBErrCorrect();
	AcquireCASLock();
	LowestCpltdProcNodeID = m_LowestCpltdProcNodeID;
	if(LowestCpltdProcNodeID == m_NumPBScaffNodes)
		{
		ReleaseCASLock();
		ReleaseCASThreadPBErrCorrect();
		goto CompletedNodeProcessing;
		}
	ReleaseCASLock();
	AcquireCASSerialise();
	if((m_CurActiveRMIThreads + m_CurActiveNonRMIThreads) >= m_CurAllowedActiveOvlpThreads)
		{
		ReleaseCASSerialise();
		ReleaseCASThreadPBErrCorrect();
		CUtility::SleepMillisecs(15000);
		continue;
		}
	RMINumUncommitedClasses = m_RMINumUncommitedClasses;
	ReleaseCASSerialise();

	if(m_bRMI && RMINumUncommitedClasses > 0)   // RMI class instance threads for SW have higher priority than non-RMI SW processing threads
		{
		if((ClassInstanceID = RMI_new(pThreadPar,cRMI_SecsTimeout,ReqSessionID))!=0)
			{
			pThreadPar->bRMI = true;
			AcquireCASSerialise();
			m_CurActiveRMIThreads += 1;
			ReleaseCASSerialise();
			ReleaseCASThreadPBErrCorrect();
			break;
			}
		ReleaseCASThreadPBErrCorrect();
		CUtility::SleepMillisecs(10000);
		continue;
		}
	else											// not using RMI classes or if using RMI then there are no RMI uncommitted class instances remaining
		{
		AcquireCASSerialise();
		if(m_CurActiveNonRMIThreads >= (UINT32)m_MaxNonRMIThreads || m_ReduceNonRMIThreads > 0)
			{
			ReleaseCASSerialise();
			ReleaseCASThreadPBErrCorrect();
			CUtility::SleepMillisecs(15000);
			continue;
			}
		m_CurActiveNonRMIThreads += 1;
		ReleaseCASSerialise();
		ReleaseCASThreadPBErrCorrect();
		break;
		}
	}

NumInMultiAlignment = 0;
AdjOverlapFloat = m_OverlapFloat + pThreadPar->CoreSeqLen + 120;

for(CurNodeID = (LowestCpltdProcNodeID+1); CurNodeID <= m_NumPBScaffNodes; CurNodeID++)
	{
	AcquireCASSerialise();				// check if needing to reduce core loading
	if(pThreadPar->bRMI == false && m_ReduceNonRMIThreads > 0)
		{
		m_ReduceNonRMIThreads -= 1;
		ReleaseCASSerialise();
		goto RMIRestartThread;
		}
	ReleaseCASSerialise();

	pThreadPar->NumCoreHits = 0;
	pCurPBScaffNode = &m_pPBScaffNodes[CurNodeID-1];
	AcquireCASLock();                    // check if sequence is to be skipped
	if(pCurPBScaffNode->flgCurProc == 1 ||			// another thread already processing this sequence? 
			pCurPBScaffNode->flgCpltdProc == 1 ||	// completed processing
			pCurPBScaffNode->flgHCseq ||			// not error correcting already assumed to be high confidence sequences 
		    pCurPBScaffNode->flgUnderlength == 1 || // must be of a user specified minimum length 
			pCurPBScaffNode->SeqLen < (UINT32)pThreadPar->MinPBSeqLen) 
		{
		ReleaseCASLock();
		continue;
		}

	if(m_PMode != ePBMConsolidate)
		pCurPBScaffNode->flgCurProc = 1;
	pCurPBScaffNode->flgCpltdProc = 0;
	ReleaseCASLock();

	ProvOverlapping = 0;
	ProvOverlapped = 0;
	ProvContained = 0;
	ProvArtefact = 0;
	ProvSWchecked = 0;

	CurChkPt.ECFileOfs = -1;
	CurChkPt.EntryID = pCurPBScaffNode->EntryID;
	CurChkPt.flgContained = pCurPBScaffNode->flgContained;
	CurChkPt.flgContains = pCurPBScaffNode->flgContains;
	CurChkPt.flgCpltdProc = 0;
	CurChkPt.flgHCseq = pCurPBScaffNode->flgHCseq;
	CurChkPt.flgUnderlength = pCurPBScaffNode->flgUnderlength;
	CurChkPt.NodeID = pCurPBScaffNode->NodeID;
	CurChkPt.SeqLen = pCurPBScaffNode->SeqLen;

	if(pCurPBScaffNode->flgSlough)
		{
		pCurPBScaffNode->flgCpltdProc = 1;
		goto SloughedSeq;
		}


	if(m_TranscriptomeLens > 0)
		MinTranscriptOverlapLen =  (pCurPBScaffNode->SeqLen * (100 - m_TranscriptomeLens)) / 100;
	else
		MinTranscriptOverlapLen = 0;

	pThreadPar->bRevCpl = false;
	IdentifyCoreHits(CurNodeID,pThreadPar);

	if(!m_bSenseOnlyOvlps)
		{
		pThreadPar->bRevCpl = true;
		IdentifyCoreHits(CurNodeID,pThreadPar);
		}

	pThreadPar->NumTargCoreHitCnts = 0;
	memset(pThreadPar->TargCoreHitCnts,0,sizeof(pThreadPar->TargCoreHitCnts));

	if(pThreadPar->NumCoreHits >= pThreadPar->MinNumCores)
		{
			// resort core hits by TargNodeID.TargOfs.ProbeNodeID.ProbeOfs ascending
		pThreadPar->pmtqsort->qsort(pThreadPar->pCoreHits,pThreadPar->NumCoreHits,sizeof(tsPBECoreHit),SortCoreHitsByTargProbeOfs);
		pCoreHit = pThreadPar->pCoreHits;
		for(HitIdx = 0; HitIdx < pThreadPar->NumCoreHits; HitIdx++, pCoreHit++)
			pCoreHit->flgMulti = 0;

		// with large target sequences then can have many artifactual core hits
		// process and mark these probable artifact hits by identifying the most spatially related cluster of hits; hits outside of
		// the most spatially related cluster are marked as being artefactual 

		UINT32 CurTargSeqID;
		UINT32 MaxWinSize;
		UINT32 RelWinSize;
		bool bFirstHitNewTargSeq;
		UINT32 PrevAcceptedProbeOfs;
		UINT32 PrevAcceptedTargOfs;
		UINT32 MaxNoHitGapLen;
		UINT32 NumNoHitGaps;
		UINT32 SumNoHitGapLens;
		tsPBECoreHit *pFirstCoreHit;
		tsPBECoreHit *pMaxCoreHit;

		if(m_PMode != ePBMConsolidate)
			MaxNoHitGapLen = 1000;                 // expecting cores to be distributed such that there shouldn't be many separated by more than this intercore bp gap
		else
			MaxNoHitGapLen = 500;
		CurTargSeqID = 0;
		pFirstCoreHit = NULL;
		ProbeSeqLen = pCurPBScaffNode->SeqLen;
		MaxWinSize = (ProbeSeqLen * 115) / 100;
		pCoreHit = pThreadPar->pCoreHits;
		pMaxCoreHit = NULL;
		for (HitIdx = 0; HitIdx < pThreadPar->NumCoreHits; HitIdx++, pCoreHit++)
			{
			if (CurTargSeqID == 0)    // 0 if 1st hit about to be processed for a new target sequence
				{
				pMaxCoreHit = NULL;
				pFirstCoreHit = pCoreHit;
				CurTargSeqID = pCoreHit->TargNodeID;
				TargSeqLen = m_pPBScaffNodes[CurTargSeqID-1].SeqLen;
				}

			// if just checked last core hit for the current target ...
			if (HitIdx + 1 == pThreadPar->NumCoreHits || pCoreHit[1].TargNodeID != CurTargSeqID)
				{
				while (pFirstCoreHit <= pCoreHit)
					{
					pFirstCoreHit->WinHits = 0;
					pFirstCoreHit->flgClustered = 0;
					pNxtCoreHit = pFirstCoreHit;
					RelWinSize = 0;
					PrevAcceptedProbeOfs = 0;
					PrevAcceptedTargOfs = 0;
					SumNoHitGapLens = 0;
					NumNoHitGaps = 0;

					if(pFirstCoreHit->ProbeOfs > MaxNoHitGapLen && pFirstCoreHit->TargOfs > MaxNoHitGapLen)
						{
						pFirstCoreHit += 1;
						continue;
						}
					
					while (pNxtCoreHit <= pCoreHit)
						{
						RelWinSize = pNxtCoreHit->TargOfs - pFirstCoreHit->TargOfs;
						if (RelWinSize > MaxWinSize)
							break;
						if (pNxtCoreHit->flgRevCpl == pFirstCoreHit->flgRevCpl)	// only interested in hits which are same sense as the first hit in window
							{
							if(PrevAcceptedProbeOfs == 0 || pNxtCoreHit->ProbeOfs >= PrevAcceptedProbeOfs + pThreadPar->CoreSeqLen) // a single matched seed core extension may have resulted in multiple hits, reduce counts by requiring a differential of at least m_SeedCoreLen
								{
								if(pFirstCoreHit->WinHits > 0 && (pNxtCoreHit->ProbeOfs - PrevAcceptedProbeOfs) >= MaxNoHitGapLen)
									{
									NumNoHitGaps += 1;
									SumNoHitGapLens += pNxtCoreHit->ProbeOfs - PrevAcceptedProbeOfs;
									}
								PrevAcceptedProbeOfs = pNxtCoreHit->ProbeOfs;
								PrevAcceptedTargOfs = pNxtCoreHit->TargOfs;
								pFirstCoreHit->WinHits += 1;
								}
							}
						pNxtCoreHit += 1;
						}

					if(pFirstCoreHit->WinHits >= pThreadPar->MinNumCores && ((PrevAcceptedProbeOfs + MaxNoHitGapLen) > ProbeSeqLen || (PrevAcceptedTargOfs + MaxNoHitGapLen) > TargSeqLen))
						{
						UINT32 PutativeOverlapLen = 1 + PrevAcceptedProbeOfs - pFirstCoreHit->ProbeOfs;
						if(((SumNoHitGapLens * 100) / PutativeOverlapLen) <= 20)		// only accepting if no more than 20% of alignment sums to gaps > MaxNoHitGapLen
							{
							if (pMaxCoreHit == NULL || pFirstCoreHit->WinHits > pMaxCoreHit->WinHits)
								pMaxCoreHit = pFirstCoreHit;
							}
						}
					pFirstCoreHit += 1;
					}

				if (pMaxCoreHit != NULL)
					{
					RelWinSize = 0;
					pNxtCoreHit = pMaxCoreHit;
					while (pNxtCoreHit <= pCoreHit)
						{
						RelWinSize = pNxtCoreHit->TargOfs - pMaxCoreHit->TargOfs;
						if (RelWinSize > MaxWinSize)
							break;
						if (pNxtCoreHit->flgRevCpl == pMaxCoreHit->flgRevCpl)	// only interested in hits which are same sense as the first hit in window
							pNxtCoreHit->flgClustered = 1;
						pNxtCoreHit->flgMulti = 0;
						pNxtCoreHit += 1;
						}
					}
				CurTargSeqID = 0;  // looking for cores on a new target
				}
			}

		// iterate and count hits for each TargNodeID whilst recording the loci of the first and last hit so can determine if overlap is a likely artefact
		CurSEntryIDHits = 0;
		CurAEntryIDHits = 0;
		CurTargSeqID = 0;
		CurTargHitOfs = 0;
		CurProbeHitOfs = 0;
		CurSTargStartOfs = 0;
		CurSTargEndOfs = 0;
		CurATargStartOfs = 0;
		CurATargEndOfs = 0;
		CurSProbeStartOfs = 0;
		CurSProbeEndOfs = 0;
		CurAProbeStartOfs = 0;
		CurAProbeEndOfs = 0;

		bFirstHitNewTargSeq = false;
		pCoreHit = pThreadPar->pCoreHits;
		for(HitIdx = 0; HitIdx < pThreadPar->NumCoreHits; HitIdx++,pCoreHit++)
			{
			if(CurTargSeqID == 0)    // 0 if 1st hit about to be processed for a new target sequence
				{
				bFirstHitNewTargSeq = true;
				CurTargSeqID = pCoreHit->TargNodeID;
				CurSEntryIDHits = 0;
				CurAEntryIDHits = 0;
				CurSTargStartOfs = 0;
				CurSTargEndOfs = 0;
				CurATargStartOfs = 0;
				CurATargEndOfs = 0;
				CurSProbeStartOfs = 0;
				CurSProbeEndOfs = 0;
				CurAProbeStartOfs = 0;
				CurAProbeEndOfs = 0;
				}

			if(pCoreHit->flgClustered && pCoreHit->TargNodeID == CurTargSeqID) // same target sequence so check for starting/ending offsets and accumulate hit counts 
				{
				CurTargHitOfs = pCoreHit->TargOfs;
				CurProbeHitOfs = pCoreHit->ProbeOfs;
				if(pCoreHit->flgRevCpl == 0)
					{
					if(bFirstHitNewTargSeq == true || CurTargHitOfs < CurSTargStartOfs)
						CurSTargStartOfs = CurTargHitOfs;
					if(CurTargHitOfs > CurSTargEndOfs)
						CurSTargEndOfs = CurTargHitOfs;	
					if(bFirstHitNewTargSeq == true || CurProbeHitOfs < CurSProbeStartOfs)
						CurSProbeStartOfs = CurProbeHitOfs;
					if(CurProbeHitOfs > CurSProbeEndOfs)
						CurSProbeEndOfs = CurProbeHitOfs;
					}
				else
					{
					if(bFirstHitNewTargSeq == true || CurTargHitOfs < CurATargStartOfs)
						CurATargStartOfs = CurTargHitOfs;
					if(CurTargHitOfs > CurATargEndOfs)
						CurATargEndOfs = CurTargHitOfs;	
					if(bFirstHitNewTargSeq == true || CurProbeHitOfs < CurAProbeStartOfs)
						CurAProbeStartOfs = CurProbeHitOfs;
					if(CurProbeHitOfs > CurAProbeEndOfs)
						CurAProbeEndOfs = CurProbeHitOfs;
					}
				bFirstHitNewTargSeq = false;
				if(pCoreHit->flgMulti != 1)
					{
					if(pCoreHit->flgRevCpl == 0)
						CurSEntryIDHits += 1;
					else
						CurAEntryIDHits += 1;
					}
				}

			// if just processed last core hit for the current target ...
			if(HitIdx + 1 == pThreadPar->NumCoreHits || pCoreHit[1].TargNodeID != CurTargSeqID)
				{
				// checking here that the first and last hit are consistent with either a completely contained or overlapped
				TargSeqLen = m_pPBScaffNodes[pCoreHit->TargNodeID-1].SeqLen;
				if(CurSEntryIDHits >= pThreadPar->MinNumCores) 
					{
					if((CurSProbeStartOfs >= m_OverlapFloat &&  CurSTargStartOfs >= m_OverlapFloat) ||
							((TargSeqLen - CurSTargEndOfs) >= AdjOverlapFloat && (ProbeSeqLen - CurSProbeEndOfs) >= AdjOverlapFloat))
						CurSEntryIDHits = 0;
					}
				else
					CurSEntryIDHits = 0;

				if(CurAEntryIDHits >= pThreadPar->MinNumCores)
					{
					if((CurAProbeStartOfs >= m_OverlapFloat && CurATargStartOfs >= m_OverlapFloat) ||
						((TargSeqLen - CurATargEndOfs) >= AdjOverlapFloat && (ProbeSeqLen - CurAProbeEndOfs) >= AdjOverlapFloat))
						CurAEntryIDHits = 0;
					}
				else
					CurAEntryIDHits = 0;

				if(CurSEntryIDHits >= pThreadPar->MinNumCores || CurAEntryIDHits >= pThreadPar->MinNumCores)
					{
					if(pThreadPar->NumTargCoreHitCnts == cSummaryTargCoreHitCnts)
						{
						LowestSummaryHitCnts = 0;
						pLowestSummaryHitCnts = NULL;
						pSummaryCnts = pThreadPar->TargCoreHitCnts; 
						for(Idx = 0; Idx < cSummaryTargCoreHitCnts; Idx++, pSummaryCnts++)
							{
							CurSummaryHitCnts = pSummaryCnts->NumSHits + pSummaryCnts->NumAHits;
							if(LowestSummaryHitCnts == 0 || CurSummaryHitCnts < LowestSummaryHitCnts)
								{
								LowestSummaryHitCnts = CurSummaryHitCnts;
								pLowestSummaryHitCnts = pSummaryCnts;
								}
							}

						if((CurSEntryIDHits + CurAEntryIDHits) <= LowestSummaryHitCnts)
							{
							CurTargSeqID = 0; 
							continue;
							}
						pSummaryCnts = pLowestSummaryHitCnts;
						}
					else
						pSummaryCnts = &pThreadPar->TargCoreHitCnts[pThreadPar->NumTargCoreHitCnts++];
					pSummaryCnts->TargNodeID = CurTargSeqID;
					pSummaryCnts->STargStartOfs = CurSTargStartOfs;
					pSummaryCnts->STargEndOfs = CurSTargEndOfs;
					pSummaryCnts->ATargStartOfs = CurATargStartOfs;
					pSummaryCnts->ATargEndOfs = CurATargEndOfs;
					pSummaryCnts->SProbeStartOfs = CurSProbeStartOfs;
					pSummaryCnts->SProbeEndOfs = CurSProbeEndOfs;
					pSummaryCnts->AProbeStartOfs = CurAProbeStartOfs;
					pSummaryCnts->AProbeEndOfs = CurAProbeEndOfs;
					pSummaryCnts->NumSHits = CurSEntryIDHits;
					pSummaryCnts->NumAHits = CurAEntryIDHits;
					pSummaryCnts->flgProbeHCseq = m_pPBScaffNodes[pCoreHit->ProbeNodeID-1].flgHCseq;
					pSummaryCnts->flgTargHCseq = m_pPBScaffNodes[pCoreHit->TargNodeID-1].flgHCseq;
					}
				CurTargSeqID = 0; 
				}
			}
		}


	// can't process, SW over all would be too resource intensive, all targets which meet the minimum number of core hits requested so choose the top cMaxProbeSWs as ranked by the number of core hits
	if(pThreadPar->NumTargCoreHitCnts > 1)
		{
		pThreadPar->pmtqsort->qsort(pThreadPar->TargCoreHitCnts,pThreadPar->NumTargCoreHitCnts,sizeof(sPBECoreHitCnts),SortCoreHitsDescending);
		if(m_PMode == ePBMConsolidate)
			{
			if(pThreadPar->NumTargCoreHitCnts > cMaxConsolidateProbeSWs)		// when consolidating (usually when generating consensus transcripts) then allow for large depth even though at most cMaxProbeSWs will be used to generate the consensus bases
				pThreadPar->NumTargCoreHitCnts = cMaxConsolidateProbeSWs;
			}
		else
			if(pThreadPar->NumTargCoreHitCnts > cMaxProbeSWs)		// if not consolidating then clamp to no more than this many SW alignments
				pThreadPar->NumTargCoreHitCnts = cMaxProbeSWs;
		}

	NumInMultiAlignment = 0;
	if(pThreadPar->NumTargCoreHitCnts > 0)
		{
		LongSAligns = 0;
		LongAAligns = 0;

		if(!pThreadPar->bRMI && pThreadPar->pSW == NULL)
			{
			bNonRMIRslt = false;
			AcquireCASSerialise();
			if((pThreadPar->pSW = new CSSW) != NULL)
				bNonRMIRslt = true;
			if(bNonRMIRslt == true)
				bNonRMIRslt = pThreadPar->pSW->SetScores(m_SWMatchScore,m_SWMismatchPenalty,m_SWGapOpenPenalty,m_SWGapExtnPenalty,m_SWProgExtnPenaltyLen,min(63,m_SWProgExtnPenaltyLen+3),cAnchorLen);
			if(bNonRMIRslt == true)
				bNonRMIRslt = pThreadPar->pSW->SetCPScores(m_SWMatchScore, m_SWMismatchPenalty, m_SWGapOpenPenalty, m_SWGapExtnPenalty);
			if(bNonRMIRslt == true)
				bNonRMIRslt = pThreadPar->pSW->SetMaxInitiatePathOfs(cDfltMaxOverlapFloat);
			if(bNonRMIRslt == true)
				bNonRMIRslt = pThreadPar->pSW->PreAllocMaxTargLen(m_MaxPBSeqLen+100, m_PMode == ePBPMConsensus ? 0 : m_MaxPBSeqLen+100);
			ReleaseCASSerialise();
			if(bNonRMIRslt == false)
				goto RMIRestartThread;
			}
		if(pThreadPar->bRMI && !bRMIInitialised)
			{
			bRMIRslt = RMI_SetScores(pThreadPar,cRMI_SecsTimeout,ClassInstanceID,m_SWMatchScore,m_SWMismatchPenalty,m_SWGapOpenPenalty,m_SWGapExtnPenalty,m_SWProgExtnPenaltyLen,min(63,m_SWProgExtnPenaltyLen+3),cAnchorLen);
			if(bRMIRslt != false)
				bRMIRslt = RMI_SetCPScores(pThreadPar,cRMI_SecsTimeout,ClassInstanceID,m_SWMatchScore, m_SWMismatchPenalty, m_SWGapOpenPenalty, m_SWGapExtnPenalty);
			if(bRMIRslt != false)
				bRMIRslt = RMI_SetMaxInitiatePathOfs(pThreadPar,cRMI_SecsTimeout,ClassInstanceID,cDfltMaxOverlapFloat);
			if(bRMIRslt != false)
				bRMIRslt = RMI_PreAllocMaxTargLen(pThreadPar,cRMI_SecsTimeout,ClassInstanceID,m_MaxPBSeqLen, m_PMode == ePBPMConsensus ? 0 : m_MaxPBSeqLen);
			if(bRMIRslt == false)
				goto RMIRestartThread;
			bRMIInitialised = true;
			}		

		if((m_PMode == ePBPMErrCorrect || m_PMode == ePBMConsolidate) && pThreadPar->NumTargCoreHitCnts >= 2)
			{
			if(!pThreadPar->bRMI)
				{
				iNonRMIRslt = pThreadPar->pSW->StartMultiAlignments(pCurPBScaffNode->SeqLen,pThreadPar->pProbeSeq,pThreadPar->NumTargCoreHitCnts,pCurPBScaffNode->flgHCseq == 1 ? (0x80 | m_HCRelWeighting) : cLCWeightingFactor);
				if(iNonRMIRslt < 0)
					goto RMIRestartThread;
				}
			else
				{
				iRMIRslt = RMI_StartMultiAlignments(pThreadPar,cRMI_SecsTimeout,ClassInstanceID,pCurPBScaffNode->SeqLen,pThreadPar->pProbeSeq,pThreadPar->NumTargCoreHitCnts,pCurPBScaffNode->flgHCseq == 1 ? (0x80 | m_HCRelWeighting)  : cLCWeightingFactor);
				if(iRMIRslt < 0)
					goto RMIRestartThread;
				}
			}

		if(pCurPBScaffNode->SeqLen > m_MaxPBSeqLen)
			gDiagnostics.DiagOut(eDLWarn,gszProcName,"####### Probe length of %dbp is longer than expected max length of %dbp #######",
						pCurPBScaffNode->SeqLen,m_MaxPBSeqLen);
		m_pSfxArray->GetIdentName(pCurPBScaffNode->EntryID,sizeof(szProbeSeqName),szProbeSeqName);
		if(!pThreadPar->bRMI)
			{
			bNonRMIRslt = pThreadPar->pSW->SetProbe(pCurPBScaffNode->SeqLen,pThreadPar->pProbeSeq);
			if(bNonRMIRslt == false)
				goto RMIRestartThread;
			}
		else
			{
			bRMIRslt = RMI_SetProbe(pThreadPar,cRMI_SecsTimeout,ClassInstanceID,pCurPBScaffNode->SeqLen,pThreadPar->pProbeSeq);
			if(bRMIRslt == false)
				goto RMIRestartThread;
			}

		// iterate over all putative targets and SW these
		UINT64 ReqSummCoverage;
		UINT64 CurSummCoverage;

		ReqSummCoverage = (UINT64)pCurPBScaffNode->SeqLen * cReqConsensusCoverage;
		CurSummCoverage = 0;

		pSummaryCnts = &pThreadPar->TargCoreHitCnts[0];
		NumInMultiAlignment = 0;
		for(CurTargCoreHitCnts = 0; CurTargCoreHitCnts < pThreadPar->NumTargCoreHitCnts && CurSummCoverage < ReqSummCoverage; CurTargCoreHitCnts++,pSummaryCnts++)
			{
			if(pSummaryCnts->NumSHits < pThreadPar->MinNumCores && pSummaryCnts->NumAHits < pThreadPar->MinNumCores)
				{
				pSummaryCnts->NumSHits = 0;
				pSummaryCnts->NumAHits = 0;
				continue;
				}
			bTargSense = pSummaryCnts->NumSHits >= pSummaryCnts->NumAHits ? true :  false;
			pTargNode = &m_pPBScaffNodes[pSummaryCnts->TargNodeID-1];

			if(m_PMode == ePBMConsolidate && pTargNode->flgCpltdProc == 1)
				{
				pSummaryCnts->NumSHits = 0;
				pSummaryCnts->NumAHits = 0;
				continue;
				}
			if(MinTranscriptOverlapLen == 0)
				{
				if(pTargNode->flgHCseq == 1)
					MinOverlapLen = m_MinHCSeqOverlap;
				else
					MinOverlapLen = pThreadPar->MinOverlapLen;
				}
			else
				MinOverlapLen = MinTranscriptOverlapLen;


			TargSeqLen = pTargNode->SeqLen; 
			if(TargSeqLen + 10 > (UINT32)pThreadPar->AllocdTargSeqSize)
				{
				delete pThreadPar->pTargSeq;
				pThreadPar->AllocdTargSeqSize = (TargSeqLen * 150) / 100;
				pThreadPar->pTargSeq = new etSeqBase [pThreadPar->AllocdTargSeqSize];
				}
			m_pSfxArray->GetSeq(pTargNode->EntryID,0,pThreadPar->pTargSeq,TargSeqLen);
			if(!bTargSense)
				CSeqTrans::ReverseComplement(TargSeqLen,pThreadPar->pTargSeq);
			pThreadPar->pTargSeq[TargSeqLen] = eBaseEOS;

			if(TargSeqLen > m_MaxPBSeqLen)
				gDiagnostics.DiagOut(eDLWarn,gszProcName,"####### Target length of %dbp is longer than expected max length of %dbp #######",
							TargSeqLen,m_MaxPBSeqLen);


			// start combined from here!
			tsCombinedTargAlignPars CombinedTargAlignPars;
			tsCombinedTargAlignRet TargAlignRet;
			memset(&CombinedTargAlignPars,0,sizeof(tsCombinedTargAlignPars));
			memset(&TargAlignRet,0,sizeof(tsCombinedTargAlignRet));
			CombinedTargAlignPars.PMode = m_PMode;
			CombinedTargAlignPars.NumTargSeqs = pThreadPar->NumTargCoreHitCnts;
			CombinedTargAlignPars.MinOverlapLen = MinOverlapLen;
			CombinedTargAlignPars.MaxOverlapLen = m_PMode == ePBPMConsensus ? 0 : m_MaxPBSeqLen;
			CombinedTargAlignPars.ProbeSeqLen = pCurPBScaffNode->SeqLen;
			CombinedTargAlignPars.TargSeqLen = TargSeqLen;
			CombinedTargAlignPars.pTargSeq = pThreadPar->pTargSeq;
			CombinedTargAlignPars.OverlapFloat = m_OverlapFloat;
			CombinedTargAlignPars.MaxArtefactDev = m_MaxArtefactDev;
			CombinedTargAlignPars.TargFlags = pSummaryCnts->flgTargHCseq == 1 ? (0x80 | m_HCRelWeighting) : cLCWeightingFactor;

			// restrict the range over which the SW will be processed to that of the overlap +/- m_OverlapFloat

			if(bTargSense)
				{
				if(pSummaryCnts->SProbeStartOfs < m_OverlapFloat)
					pSummaryCnts->SProbeStartOfs = 0;
				else
					pSummaryCnts->SProbeStartOfs -= m_OverlapFloat;
				if(pSummaryCnts->SProbeEndOfs + AdjOverlapFloat >= pCurPBScaffNode->SeqLen)
					pSummaryCnts->SProbeEndOfs = pCurPBScaffNode->SeqLen - 1;
				else
					pSummaryCnts->SProbeEndOfs += AdjOverlapFloat;
				if(pSummaryCnts->STargStartOfs < m_OverlapFloat)
					pSummaryCnts->STargStartOfs = 0;
				else
					pSummaryCnts->STargStartOfs -= m_OverlapFloat;
				if(pSummaryCnts->STargEndOfs + AdjOverlapFloat >= TargSeqLen)
					pSummaryCnts->STargEndOfs = TargSeqLen - 1;
				else
					pSummaryCnts->STargEndOfs += AdjOverlapFloat;

				CombinedTargAlignPars.ProbeStartRelOfs = pSummaryCnts->SProbeStartOfs;
				CombinedTargAlignPars.TargStartRelOfs = pSummaryCnts->STargStartOfs;
				CombinedTargAlignPars.ProbeRelLen = pSummaryCnts->SProbeEndOfs + 1 - pSummaryCnts->SProbeStartOfs;
				CombinedTargAlignPars.TargRelLen = pSummaryCnts->STargEndOfs + 1 - pSummaryCnts->STargStartOfs;
				}
			else
				{
				UINT32 Xchg;
				Xchg = pSummaryCnts->AProbeStartOfs;
				pSummaryCnts->AProbeStartOfs = pCurPBScaffNode->SeqLen - (pSummaryCnts->AProbeEndOfs + 1);
				pSummaryCnts->AProbeEndOfs = pCurPBScaffNode->SeqLen - (Xchg + 1);
				Xchg = pSummaryCnts->ATargStartOfs;
				pSummaryCnts->ATargStartOfs = TargSeqLen - (pSummaryCnts->ATargEndOfs + 1);
				pSummaryCnts->ATargEndOfs = TargSeqLen - (Xchg + 1);

				if(pSummaryCnts->AProbeStartOfs < m_OverlapFloat)
					pSummaryCnts->AProbeStartOfs = 0;
				else
					pSummaryCnts->AProbeStartOfs -= m_OverlapFloat;
				if(pSummaryCnts->AProbeEndOfs + AdjOverlapFloat >= pCurPBScaffNode->SeqLen)
					pSummaryCnts->AProbeEndOfs = pCurPBScaffNode->SeqLen - 1;
				else
					pSummaryCnts->AProbeEndOfs += AdjOverlapFloat;
				if(pSummaryCnts->ATargStartOfs < m_OverlapFloat)
					pSummaryCnts->ATargStartOfs = 0;
				else
					pSummaryCnts->ATargStartOfs -= m_OverlapFloat;
				if(pSummaryCnts->ATargEndOfs + AdjOverlapFloat >= TargSeqLen)
					pSummaryCnts->ATargEndOfs = TargSeqLen - 1;
				else
					pSummaryCnts->ATargEndOfs += AdjOverlapFloat;

				CombinedTargAlignPars.ProbeStartRelOfs = pSummaryCnts->AProbeStartOfs;
				CombinedTargAlignPars.TargStartRelOfs = pSummaryCnts->ATargStartOfs;
				CombinedTargAlignPars.ProbeRelLen = pSummaryCnts->AProbeEndOfs + 1 - pSummaryCnts->AProbeStartOfs;
				CombinedTargAlignPars.TargRelLen = pSummaryCnts->ATargEndOfs + 1 - pSummaryCnts->ATargStartOfs;
				}

			if(CombinedTargAlignPars.ProbeRelLen < MinOverlapLen || CombinedTargAlignPars.TargRelLen < MinOverlapLen)
				continue;


			if(!pThreadPar->bRMI)
				{
				if(!pThreadPar->pSW->CombinedTargAlign(&CombinedTargAlignPars, &TargAlignRet))
					goto CompletedNodeProcessing;
				if(TargAlignRet.ErrRslt != eBSFSuccess || TargAlignRet.ProcPhase < 2 || TargAlignRet.ProcPhase == 3)
					goto CompletedNodeProcessing;
				}
			else
				{
				if(!RMI_CombinedTargAlign(pThreadPar,cRMI_AlignSecsTimeout,ClassInstanceID,&CombinedTargAlignPars, &TargAlignRet))
					goto RMIRestartThread;
				if(TargAlignRet.ErrRslt != eBSFSuccess || TargAlignRet.ProcPhase < 2 || TargAlignRet.ProcPhase == 3)
					goto RMIRestartThread;
				}

			PeakMatchesCell = TargAlignRet.PeakMatchesCell;
			pPeakMatchesCell = &PeakMatchesCell;

			ProvSWchecked += 1;
			if(TargAlignRet.Flags & 0x02)		// set if alignment classified as an artifact
				ProvArtefact += 1;
			
			if(TargAlignRet.ProcPhase == 4)
				{
				if(TargAlignRet.Flags & 0x08)		// set if alignment was accepted and added as a multialignment
					{
					NumInMultiAlignment += 1;
					if(m_PMode == ePBMConsolidate)
						{
						pTargNode->flgCpltdProc = 1;
						CurSummCoverage = 0;
						}
					else
						{
						UINT64 OvlpLen =  (UINT64)(TargAlignRet.ProbeAlignLength + TargAlignRet.TargAlignLength + 1) / 2;
						if(pSummaryCnts->flgTargHCseq == 1)
							CurSummCoverage += OvlpLen * 3 / 2;  // if alignment was onto a high confidence sequence then less coverage depth required for consensus confidence 
						else
							CurSummCoverage += OvlpLen;
						}

					if(TargAlignRet.Flags & 0x04)		// set if alignment classified as contained
						ProvContained += 1;

					if(TargAlignRet.Flags & 0x01)		// set if alignment was classified as overlapping 
						ProvOverlapping += 1;
					ProvOverlapped += 1;
					}

				if(m_PMode == ePBPMOverlapDetail && m_hErrCorFile != -1)
					{
					m_pSfxArray->GetIdentName(pTargNode->EntryID,sizeof(szTargSeqName)-1,szTargSeqName);
					AcquireCASSerialise();
					if(m_ScaffLineBuffIdx > (sizeof(m_szScaffLineBuff) - 1000))
						{
						if(!CUtility::SafeWrite(m_hErrCorFile,m_szScaffLineBuff,m_ScaffLineBuffIdx))
							{
							gDiagnostics.DiagOut(eDLFatal,gszProcName,"SafeWrite() failed writing %u chars to overlap details file",m_ScaffLineBuffIdx);
							ReleaseCASSerialise();
							goto CompletedNodeProcessing;
							}
						m_ScaffLineBuffIdx = 0;
						}

					m_ScaffLineBuffIdx += sprintf(&m_szScaffLineBuff[m_ScaffLineBuffIdx], "\n%d,%d,\"%s\",%d,\"%s\",%d,\"%c\",\"%c\",%1.5d,%1.5d,%1.5d,%1.5d,%1.5d,%1.5d,%1.5d,%1.5d,%1.5d,%1.5d,%1.5d,%1.5d,%1.5d,%1.5d,%1.5d,%1.5d,%1.5d,%1.5d,%1.5d,%1.5d",
														TargAlignRet.Class,pCurPBScaffNode->EntryID,szProbeSeqName,pTargNode->EntryID,szTargSeqName,
														bTargSense ? pSummaryCnts->NumSHits : pSummaryCnts->NumAHits, 'S', bTargSense ? 'S' : 'A',
														pCurPBScaffNode->SeqLen,TargSeqLen,
														TargAlignRet.ProbeAlignLength, TargAlignRet.TargAlignLength, PeakMatchesCell.PeakScore, PeakMatchesCell.CurScore, PeakMatchesCell.NumMatches,
														PeakMatchesCell.NumExacts, PeakMatchesCell.NumGapsIns, PeakMatchesCell.NumBasesIns, PeakMatchesCell.NumGapsDel, PeakMatchesCell.NumBasesDel, PeakMatchesCell.StartPOfs,
														PeakMatchesCell.StartTOfs, PeakMatchesCell.EndPOfs, PeakMatchesCell.EndTOfs, PeakMatchesCell.PFirstAnchorStartOfs, PeakMatchesCell.TFirstAnchorStartOfs, PeakMatchesCell.PLastAnchorEndOfs, PeakMatchesCell.TLastAnchorEndOfs);

					ReleaseCASSerialise();
					}
				}
			}
		}

	if((m_PMode == ePBPMErrCorrect || m_PMode == ePBMConsolidate))
		{
		if((pThreadPar->pSW != NULL || pThreadPar->bRMI) && NumInMultiAlignment >= 1 &&  (m_hMultiAlignFile != -1 || m_hErrCorFile != -1))
			{
			if(!pThreadPar->bRMI)
				pThreadPar->pSW->GenMultialignConcensus();
			else
				{
				iRMIRslt = RMI_GenMultialignConcensus(pThreadPar,cRMI_SecsTimeout,ClassInstanceID);
				if(iRMIRslt < 0)
					goto RMIRestartThread;
				}

			if(m_hErrCorFile != -1)
				{
				if(!pThreadPar->bRMI)
					pThreadPar->ErrCorBuffIdx=pThreadPar->pSW->MAlignCols2fasta(pCurPBScaffNode->EntryID,m_MinConcScore,m_MinErrCorrectLen,pThreadPar->AllocdErrCorLineBuff,pThreadPar->pszErrCorLineBuff);
				else
					{
					pThreadPar->ErrCorBuffIdx = iRMIRslt = RMI_MAlignCols2fasta(pThreadPar,cRMI_SecsTimeout,ClassInstanceID,pCurPBScaffNode->EntryID,m_MinConcScore,m_MinErrCorrectLen,pThreadPar->AllocdErrCorLineBuff,pThreadPar->pszErrCorLineBuff);
					if(iRMIRslt < 0)
						goto RMIRestartThread;
					}
				if(pThreadPar->ErrCorBuffIdx > 0)
					{
					AcquireCASSerialise();
					if(!CUtility::SafeWrite(m_hErrCorFile,pThreadPar->pszErrCorLineBuff,pThreadPar->ErrCorBuffIdx))
						{
						gDiagnostics.DiagOut(eDLFatal,gszProcName,"SafeWrite() failed writing %u chars to error corrected reads file",pThreadPar->ErrCorBuffIdx);
						ReleaseCASSerialise();
						goto CompletedNodeProcessing;
						}
					m_ErrCorFileUnsyncedSize += pThreadPar->ErrCorBuffIdx;
					if(m_ErrCorFileUnsyncedSize > 50000)
						{
	#ifdef _WIN32
						_commit(m_hErrCorFile);
	#else
						fsync(m_hErrCorFile);
	#endif
						m_ErrCorFileUnsyncedSize = 0;
						}

					if(m_hChkPtsFile != -1)
						{
						CurChkPt.ECFileOfs = (INT64)_lseeki64(m_hErrCorFile,(off_t)0,SEEK_CUR);
						CurChkPt.flgContained = pCurPBScaffNode->flgContained;
						CurChkPt.flgContains = pCurPBScaffNode->flgContains;
						CurChkPt.flgCpltdProc = 1;
						if(!CUtility::SafeWrite(m_hChkPtsFile,&CurChkPt,sizeof(tsECChkPt)))
							{
							gDiagnostics.DiagOut(eDLFatal,gszProcName,"SafeWrite() failed writing to checkpoint file");
							ReleaseCASSerialise();
							goto CompletedNodeProcessing;
							}
	#ifdef _WIN32
						_commit(m_hChkPtsFile);
	#else
						fsync(m_hChkPtsFile);
	#endif
						}
					ReleaseCASSerialise();
					pThreadPar->ErrCorBuffIdx = 0;
					}
				}

			if(m_hMultiAlignFile != -1)
				{
				if(!pThreadPar->bRMI)
					pThreadPar->MultiAlignBuffIdx=pThreadPar->pSW->MAlignCols2MFA(pCurPBScaffNode->EntryID,pThreadPar->AllocdMultiAlignLineBuff,pThreadPar->pszMultiAlignLineBuff);
				else
					{
					pThreadPar->MultiAlignBuffIdx = iRMIRslt = RMI_MAlignCols2MFA(pThreadPar,cRMI_SecsTimeout,ClassInstanceID,pCurPBScaffNode->EntryID,pThreadPar->AllocdMultiAlignLineBuff,pThreadPar->pszMultiAlignLineBuff);
					if(iRMIRslt < 0)
						goto RMIRestartThread;
					}
				if(pThreadPar->MultiAlignBuffIdx > 0)
					{
					pThreadPar->pszMultiAlignLineBuff[pThreadPar->MultiAlignBuffIdx++] = '\n';
					AcquireCASSerialise();
					if(!CUtility::SafeWrite(m_hMultiAlignFile,pThreadPar->pszMultiAlignLineBuff,pThreadPar->MultiAlignBuffIdx))
						{
						gDiagnostics.DiagOut(eDLFatal,gszProcName,"SafeWrite() failed writing %u chars to multialignment file",pThreadPar->MultiAlignBuffIdx);
						ReleaseCASSerialise();
						goto CompletedNodeProcessing;
						}
					m_MultiAlignFileUnsyncedSize += pThreadPar->MultiAlignBuffIdx;
					if(m_MultiAlignFileUnsyncedSize > 100000)
						{
	#ifdef _WIN32
						_commit(m_hMultiAlignFile);
	#else
						fsync(m_hMultiAlignFile);
	#endif
						m_MultiAlignFileUnsyncedSize = 0;
						}
					ReleaseCASSerialise();
					pThreadPar->MultiAlignBuffIdx = 0;
					}
				}
			}
		else
			if(m_PMode == ePBMConsolidate &&  m_hErrCorFile != -1)
				{
				// asciify the consensus and report to m_hErrCorFile
				UINT32 BaseIdx = 0;
				int BasesLine = 0;

				pThreadPar->ErrCorBuffIdx += sprintf(&pThreadPar->pszErrCorLineBuff[pThreadPar->ErrCorBuffIdx],">ecseq%u_1 %d|%d\n",pCurPBScaffNode->NodeID,pCurPBScaffNode->SeqLen,1);
				while(BaseIdx < pCurPBScaffNode->SeqLen)
					{
					BasesLine = pCurPBScaffNode->SeqLen - BaseIdx > 80 ? 80 : pCurPBScaffNode->SeqLen - BaseIdx;
					CSeqTrans::MapSeq2Ascii(&pThreadPar->pProbeSeq[BaseIdx],BasesLine,&pThreadPar->pszErrCorLineBuff[pThreadPar->ErrCorBuffIdx]);
					BaseIdx += BasesLine;
					pThreadPar->ErrCorBuffIdx += BasesLine;
					pThreadPar->ErrCorBuffIdx += sprintf(&pThreadPar->pszErrCorLineBuff[pThreadPar->ErrCorBuffIdx],"\n");
					if(pThreadPar->ErrCorBuffIdx + 100 > pThreadPar->AllocdErrCorLineBuff)
						{
						if(!CUtility::SafeWrite(m_hErrCorFile,pThreadPar->pszErrCorLineBuff,pThreadPar->ErrCorBuffIdx))
							{
							gDiagnostics.DiagOut(eDLFatal,gszProcName,"SafeWrite() failed writing %u chars to error corrected reads file",pThreadPar->ErrCorBuffIdx);
							ReleaseCASSerialise();
							goto CompletedNodeProcessing;
							}		
						m_ErrCorFileUnsyncedSize += pThreadPar->ErrCorBuffIdx;				
						pThreadPar->ErrCorBuffIdx = 0;
						}
					}
				if(pThreadPar->ErrCorBuffIdx > 0)
					{
					AcquireCASSerialise();
					if(!CUtility::SafeWrite(m_hErrCorFile,pThreadPar->pszErrCorLineBuff,pThreadPar->ErrCorBuffIdx))
						{
						gDiagnostics.DiagOut(eDLFatal,gszProcName,"SafeWrite() failed writing %u chars to error corrected reads file",pThreadPar->ErrCorBuffIdx);
						ReleaseCASSerialise();
						goto CompletedNodeProcessing;
						}
					m_ErrCorFileUnsyncedSize += pThreadPar->ErrCorBuffIdx;
					pThreadPar->ErrCorBuffIdx = 0;
					}

				if(m_ErrCorFileUnsyncedSize > 50000)
					{
		#ifdef _WIN32
					_commit(m_hErrCorFile);
		#else
					fsync(m_hErrCorFile);
		#endif
					m_ErrCorFileUnsyncedSize = 0;
					}

				if(m_hChkPtsFile != -1)
					{
					CurChkPt.ECFileOfs = (INT64)_lseeki64(m_hErrCorFile,(off_t)0,SEEK_CUR);
					CurChkPt.flgContained = pCurPBScaffNode->flgContained;
					CurChkPt.flgContains = pCurPBScaffNode->flgContains;
					CurChkPt.flgCpltdProc = 1;
					if(!CUtility::SafeWrite(m_hChkPtsFile,&CurChkPt,sizeof(tsECChkPt)))
						{
						gDiagnostics.DiagOut(eDLFatal,gszProcName,"SafeWrite() failed writing to checkpoint file");
						ReleaseCASSerialise();
						goto CompletedNodeProcessing;
						}
		#ifdef _WIN32
					_commit(m_hChkPtsFile);
		#else
					fsync(m_hChkPtsFile);
		#endif
					}
				ReleaseCASSerialise();
				}
		if(m_PMode == ePBMConsolidate)
			pCurPBScaffNode->flgCpltdProc = 1;
		}

SloughedSeq:     // branch to here for target sequences marked as not to be processed
	AcquireCASSerialise();
	if(ProvOverlapping > 0)
		m_ProvOverlapping += 1;
	if(ProvOverlapped > 0)
		m_ProvOverlapped += ProvOverlapped;
	if(ProvContained > 0)
		m_ProvContained += ProvContained;
	if(ProvArtefact > 0)
		m_ProvArtefact += ProvArtefact;
	if(ProvSWchecked > 0)
		m_ProvSWchecked += ProvSWchecked;
	m_NumOverlapProcessed += 1;
	if(m_PMode == ePBMConsolidate)
		m_NumOverlapProcessed += NumInMultiAlignment;

	if(m_hChkPtsFile != -1 && CurChkPt.flgCpltdProc == 0)
		{
		CurChkPt.ECFileOfs = -1;
		CurChkPt.flgContained = pCurPBScaffNode->flgContained;
		CurChkPt.flgContains = pCurPBScaffNode->flgContains;
		CurChkPt.flgCpltdProc = 1;
		if(!CUtility::SafeWrite(m_hChkPtsFile,&CurChkPt,sizeof(tsECChkPt)))
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"SafeWrite() failed writing to checkpoint file");
			ReleaseCASSerialise();
			goto CompletedNodeProcessing;
			}
#ifdef _WIN32
		_commit(m_hChkPtsFile);
#else
		fsync(m_hChkPtsFile);
#endif
		}

	ReleaseCASSerialise();
	AcquireCASLock();
	pCurPBScaffNode->flgCpltdProc = 1;
	ReleaseCASLock();
	ProvOverlapping = 0;
	ProvOverlapped = 0;
	ProvContained = 0;
	ProvArtefact = 0;
	ProvSWchecked = 0;
	NumInMultiAlignment = 0;
	pThreadPar->NumTargCoreHitCnts = 0;
	pThreadPar->NumCoreHits = 0;
	}

CompletedNodeProcessing:     // when no more nodes requiring processing then goto is used to branch here for thread cleanup
AcquireCASSerialise();
if((m_PMode == ePBPMErrCorrect  || m_PMode == ePBMConsolidate) && m_hErrCorFile != -1 && pThreadPar->ErrCorBuffIdx > 0)
	{
	CUtility::SafeWrite(m_hErrCorFile,pThreadPar->pszErrCorLineBuff,pThreadPar->ErrCorBuffIdx);
	pThreadPar->ErrCorBuffIdx = 0;
#ifdef _WIN32
	_commit(m_hErrCorFile);
#else
	fsync(m_hErrCorFile);
#endif
	m_ErrCorFileUnsyncedSize = 0;
	}
if(m_hMultiAlignFile != -1 && pThreadPar->MultiAlignBuffIdx > 0)
	{
	CUtility::SafeWrite(m_hMultiAlignFile,pThreadPar->pszMultiAlignLineBuff,pThreadPar->MultiAlignBuffIdx);
	pThreadPar->MultiAlignBuffIdx = 0;
#ifdef _WIN32
	_commit(m_hMultiAlignFile);
#else
	fsync(m_hMultiAlignFile);
#endif
	m_MultiAlignFileUnsyncedSize = 0;
	}
if(m_PMode == ePBPMOverlapDetail && m_ScaffLineBuffIdx > 0)
	{
	CUtility::SafeWrite(m_hErrCorFile,m_szScaffLineBuff,m_ScaffLineBuffIdx);
	m_ScaffLineBuffIdx = 0;
#ifdef _WIN32
	_commit(m_hErrCorFile);
#else
	fsync(m_hErrCorFile);
	m_ErrCorFileUnsyncedSize = 0;
#endif
	}
ReleaseCASSerialise();

if(pThreadPar->bRMI && ClassInstanceID != 0)
	{
	RMI_delete(pThreadPar,cRMI_SecsTimeout,ClassInstanceID);
	ClassInstanceID = 0;
	}
if(pThreadPar->pSW != NULL)
	{
	delete pThreadPar->pSW;
	pThreadPar->pSW = NULL;
	}
return(0);
}


int					// returns 0 if core overlapped (uses a non-exhaustive search) a previously added core, index 1..N of just added core hit or -1 if errors
CPBErrCorrect::AddCoreHit(UINT32 ProbeNodeID,			// core hit was from this probe scaffold node 
			   bool bRevCpl,					// true if core sequence was revcpl'd before matching
			   UINT32 ProbeOfs,                 // hit started at this probe offset
			   UINT32 TargNodeID,               // probe core matched onto this target scaffold node
			   UINT32 TargOfs,                  // probe core matched starting at this target loci
			   UINT32 HitLen,					// hit was of this length
               tsThreadPBErrCorrect *pPars)			// thread specific
{
UINT32 ExpdHitLen;
UINT32 NumHits2Chk;
UINT32 HitsChkd;
tsPBECoreHit *pCurCoreHit;

if((pPars->NumCoreHits + 5) > pPars->AllocdCoreHits)	// need to realloc memory to hold additional cores?
	{
		// realloc memory with a 25% increase over previous allocation 
	int coresreq;
	size_t memreq;
	void *pAllocd;
	coresreq = (int)(((INT64)pPars->AllocdCoreHits * 125) / (INT64)100);
	memreq = coresreq * sizeof(tsPBECoreHit);

#ifdef _WIN32
		pAllocd = realloc(pPars->pCoreHits,memreq);
#else
		pAllocd = mremap(pPars->pCoreHits,pPars->AllocdCoreHitsSize,memreq,MREMAP_MAYMOVE);
		if(pAllocd == MAP_FAILED)
			pAllocd = NULL;
#endif
		if(pAllocd == NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"SavePartialSeqs: Memory re-allocation to %d bytes - %s",memreq,strerror(errno));
			Reset();
			return(eBSFerrMem);
			}

		pPars->pCoreHits = (tsPBECoreHit *)pAllocd;
		pPars->AllocdCoreHitsSize = memreq;
		pPars->AllocdCoreHits = coresreq; 
		}

// non-exhaustive check to see if core is overlapping a previously added core
// search back over at most 2000 recently added cores for overlaps
ExpdHitLen = min(50,HitLen * 3); // allowing for some float on the overlap loci and length
if(pPars->NumCoreHits > 0)
	{
	NumHits2Chk = min(2000,pPars->NumCoreHits);
	pCurCoreHit=&pPars->pCoreHits[pPars->NumCoreHits-1];
	for(HitsChkd = 0; HitsChkd < NumHits2Chk; HitsChkd+=1, pCurCoreHit-=1)
		{
		if(pCurCoreHit->TargNodeID == TargNodeID &&
		  pCurCoreHit->ProbeNodeID == ProbeNodeID &&
			pCurCoreHit->flgRevCpl == (bRevCpl ? 1 : 0) &&
			(pCurCoreHit->ProbeOfs >= (ProbeOfs < ExpdHitLen ? 0 : ProbeOfs - ExpdHitLen) && pCurCoreHit->ProbeOfs <= (ProbeOfs + ExpdHitLen)) &&
			(pCurCoreHit->TargOfs >= (TargOfs < ExpdHitLen ? 0 : TargOfs - ExpdHitLen) && pCurCoreHit->TargOfs <= (TargOfs + ExpdHitLen)))
			return(0);
		}	
	}
 
		
pCurCoreHit = &pPars->pCoreHits[pPars->NumCoreHits++];

pCurCoreHit->ProbeNodeID = ProbeNodeID;
pCurCoreHit->flgRevCpl = bRevCpl ? 1 : 0;
pCurCoreHit->flgMulti = 0;
pCurCoreHit->ProbeOfs = ProbeOfs;
pCurCoreHit->TargNodeID = TargNodeID;
pCurCoreHit->HitLen = HitLen;
pCurCoreHit->TargOfs = TargOfs;
memset(&pCurCoreHit[1],0,sizeof(tsPBECoreHit));	// ensuring that used cores are always terminated with a marker end of cores initialised to 0
return(pPars->NumCoreHits);
}



// MapEntryID2NodeID
// Given a suffix array entry identifier returns the corresponding node identifier
UINT32										// returned tsPBScaffNode node identifier
CPBErrCorrect::MapEntryID2NodeID(UINT32 EntryID)			// suffix array entry identifier
{
if(EntryID == 0 || EntryID > m_NumPBScaffNodes || m_pMapEntryID2NodeIDs == NULL)
	return(0);
return(m_pMapEntryID2NodeIDs[EntryID-1]);
}



int
CPBErrCorrect::IdentifyCoreHits(UINT32 ProbeNodeID,	// identify all overlaps of this probe sequence ProbeNodeID onto target sequences
				tsThreadPBErrCorrect *pPars)		// thread specific
{
INT64 PrevHitIdx;
INT64 NextHitIdx;
UINT32 HitEntryID;
UINT32 HitLoci;
UINT32 HitsThisCore;
UINT32 HighHitsThisCore;
UINT32 TotHitsAllCores;
UINT32 HitSeqLen;
UINT32 TargOvlpLenDiffbp;
UINT32 ProbeOfs;
UINT32 LastProbeOfs;

UINT32 ChkOvrLapCoreProbeOfs;
UINT32 LastCoreProbeOfs;
int ChkOvrLapCoreStartIdx;

etSeqBase *pCoreSeq;
tsPBEScaffNode *pProbeNode;
tsPBEScaffNode *pTargNode;
etSeqBase *pHomo;
int HomoIdx;
int HomoBaseCnts[4];
int MaxAcceptHomoCnt;

tsQualTarg QualTargs[cMaxQualTargs];
UINT32 QualCoreLen;
int NumQualSeqs;


if(ProbeNodeID < 1 || ProbeNodeID > m_NumPBScaffNodes)
	return(eBSFerrParams);
pProbeNode = &m_pPBScaffNodes[ProbeNodeID-1];

if(pPars->MinOverlapLen > pProbeNode->SeqLen)
	return(0);

if(pPars->pProbeSeq == NULL || pPars->AllocdProbeSeqSize < (pProbeNode->SeqLen + 10))
	{
	if(pPars->pProbeSeq != NULL)
		delete pPars->pProbeSeq;
	pPars->AllocdProbeSeqSize = pProbeNode->SeqLen;
	pPars->pProbeSeq = new etSeqBase [pProbeNode->SeqLen + 10];
	}
m_pSfxArray->GetSeq(pProbeNode->EntryID,0,pPars->pProbeSeq,pProbeNode->SeqLen);
pPars->pProbeSeq[pProbeNode->SeqLen] = eBaseEOS;

if(pPars->bRevCpl)
	CSeqTrans::ReverseComplement(pProbeNode->SeqLen,pPars->pProbeSeq);
pCoreSeq = pPars->pProbeSeq;
ChkOvrLapCoreProbeOfs = 0;
ChkOvrLapCoreStartIdx = 0;
LastCoreProbeOfs = 0;
PrevHitIdx = 0;
HitsThisCore = 0;
HighHitsThisCore = 0;
TotHitsAllCores = 0;
LastProbeOfs = 1 + pProbeNode->SeqLen - pPars->CoreSeqLen;
if(pPars->CoreSeqLen < cMaxPacBioSeedExtn)
	LastProbeOfs -= 120;

TargOvlpLenDiffbp = m_TranscriptomeLens == 0 ? 0 : ((pProbeNode->SeqLen * m_TranscriptomeLens) / 100);
memset(QualTargs,0,sizeof(QualTargs));
QualCoreLen = pPars->CoreSeqLen + 5;
NumQualSeqs = m_pSfxArray->PreQualTargs(pProbeNode->EntryID,pProbeNode->SeqLen,pPars->pProbeSeq, QualCoreLen,(int)pPars->DeltaCoreOfs, TargOvlpLenDiffbp, cMaxQualTargs,QualTargs);

if(NumQualSeqs > 0)
	{
	pCoreSeq = pPars->pProbeSeq;
	MaxAcceptHomoCnt = (pPars->CoreSeqLen * cQualCoreHomopolymer) / 100; // if any core contains more than cQualCoreHomopolymer% of the same base then treat as being a near homopolymer core and slough
	for(ProbeOfs = 0; ProbeOfs < LastProbeOfs; ProbeOfs+=pPars->DeltaCoreOfs,pCoreSeq+=pPars->DeltaCoreOfs)
		{
		// with PacBio reads most homopolymer runs are actual artefact inserts so don't bother processing these homopolymer runs for cores
		if(MaxAcceptHomoCnt > 0)
			{
			HomoBaseCnts[0] = HomoBaseCnts[1] = HomoBaseCnts[2] = HomoBaseCnts[3] = 0;
			for(pHomo = pCoreSeq, HomoIdx = 0; HomoIdx < (int)pPars->CoreSeqLen; HomoIdx+=1, pHomo += 1)
				HomoBaseCnts[*pHomo & 0x03] += 1;
			if(HomoBaseCnts[0] > MaxAcceptHomoCnt || HomoBaseCnts[1] > MaxAcceptHomoCnt || HomoBaseCnts[2] > MaxAcceptHomoCnt || HomoBaseCnts[3] > MaxAcceptHomoCnt)
				continue;
			}


		PrevHitIdx = 0;
		HitsThisCore = 0;

		while((NextHitIdx = m_pSfxArray->IteratePacBio(pCoreSeq,pProbeNode->SeqLen - ProbeOfs,pPars->CoreSeqLen,pProbeNode->EntryID,pPars->MinPBSeqLen,PrevHitIdx,&HitEntryID,&HitLoci,NumQualSeqs,QualTargs)) > 0)
			{
			PrevHitIdx = NextHitIdx;
			pTargNode = &m_pPBScaffNodes[MapEntryID2NodeID(HitEntryID)-1];
			HitSeqLen = pTargNode->SeqLen; 

			if(TargOvlpLenDiffbp > 0 && ((pProbeNode->SeqLen + TargOvlpLenDiffbp) < HitSeqLen  || pProbeNode->SeqLen > (HitSeqLen + TargOvlpLenDiffbp)))
				continue;

			if((pPars->bSelfHits ? pTargNode->NodeID != pProbeNode->NodeID : pTargNode->NodeID == pProbeNode->NodeID) || 
								HitSeqLen < (UINT32)pPars->MinOverlapLen)		// not interested if target sequence length less than min overlap length required
				continue;

  			if(AddCoreHit(ProbeNodeID,pPars->bRevCpl,ProbeOfs,pTargNode->NodeID,HitLoci,pPars->CoreSeqLen,pPars)>0)
				{
				HitsThisCore += 1;
				if(HitsThisCore > pPars->MaxAcceptHitsPerSeedCore)
					break;
				}
			}
		if(HitsThisCore)	// if at least one hit from this core
			{
			if(HitsThisCore > HighHitsThisCore)
				HighHitsThisCore = HitsThisCore;
			TotHitsAllCores += HitsThisCore;
			}
		}
	}
if(pPars->bRevCpl)
	CSeqTrans::ReverseComplement(pProbeNode->SeqLen,pPars->pProbeSeq);

return(pPars->NumCoreHits);
}

int											// marshaled parameter required this many bytes
CPBErrCorrect::MarshalReq(UINT8 *pInto,					// marshal into this list
				teRMIParamType Type,			// parameter type
				void *pValue,				// parameter value
				UINT32 ValLen)				// length of parameter ptd to by pValue, only used if parameter type is pUint8

{
switch(Type) {
	case eRMIPTBool:	// boolean
		*pInto++ = (UINT8)eRMIPTBool;
		*pInto = *(bool *)pValue == true ? 1 : 0;
		return(sizeof(UINT8) + 1);

	case eRMIPTInt8:		// 8bit signed int
		*pInto++ = (UINT8)eRMIPTInt8;
		*pInto = *(INT8 *)pValue;
		return(sizeof(INT8) + 1);

	case eRMIPTUint8:       // 8bit  unsigned int
		*pInto++ = (UINT8)eRMIPTUint8;
		*pInto = *(UINT8 *)pValue;
		return(sizeof(UINT8) + 1);

	case eRMIPTInt32:		// 32bit signed int
		*pInto++ = (UINT8)eRMIPTInt32;
		*(INT32 *)pInto = *(INT32 *)pValue;
		return(sizeof(INT32) + 1);

	case eRMIPTUint32:		// 32bit unsigned int
		*pInto++ = (UINT8)eRMIPTInt32;
		*(UINT32 *)pInto = *(UINT32 *)pValue;
		return(sizeof(UINT32) + 1);

	case eRMIPTInt64:		// 64bit signed int
		*pInto++ = (UINT8)eRMIPTInt64;
		*(INT64 *)pInto = *(INT64 *)pValue;
		return(sizeof(INT64) + 1);

	case eRMIPTUint64:		// 64bit unsigned int
		*pInto++ = (UINT8)eRMIPTUint64;
		*(UINT64 *)pInto = *(UINT64 *)pValue;
		return(sizeof(UINT64) + 1);

	case eRMIPTDouble:		// floating point double
		*pInto++ = (UINT8)eRMIPTDouble;
		*(double *)pInto = *(double *)pValue;
		return(sizeof(double) + 1);

	case eRMIPTVarUint8:		// variable length
		*pInto++ = (UINT8)eRMIPTVarUint8;
		*(UINT32 *)pInto = ValLen;
		if(ValLen > 0 && pValue != NULL)
			{
			pInto += sizeof(UINT32);
			memcpy(pInto,pValue,ValLen);
			}
		return(sizeof(UINT32) + ValLen + 1);

	default:
		break;
	}
	
return(0);
}

int
CPBErrCorrect::UnmarshalResp(UINT32 DataLen,
				UINT8 *pFrom,		// unmarshal from this marshalled parameter list
				void *pValue)
{
UINT32 ValLen;
switch(*pFrom++) {
	case eRMIPTBool:	// boolean
		*(bool *)pValue = *pFrom == 0 ? false : true;
		return(sizeof(UINT8) + 1);

	case eRMIPTInt8:		// 8bit signed int
		*(INT8 *)pValue = *(INT8 *)pFrom;
		return(sizeof(INT8) + 1);

	case eRMIPTUint8:       // 8bit  unsigned int
		*(UINT8 *)pValue = *(UINT8 *)pFrom;
		return(sizeof(UINT8) + 1);

	case eRMIPTInt32:		// 32bit signed int
		*(INT32 *)pValue = *(INT32 *)pFrom;
		return(sizeof(INT32) + 1);

	case eRMIPTUint32:		// 32bit unsigned int
		*(UINT32 *)pValue = *(UINT32 *)pFrom;
		return(sizeof(UINT32) + 1);

	case eRMIPTInt64:		// 64bit signed int
		*(INT64 *)pValue = *(INT64 *)pFrom;
		return(sizeof(INT64) + 1);

	case eRMIPTUint64:		// 64bit unsigned int
		*(UINT64 *)pValue = *(UINT64 *)pFrom;
		return(sizeof(UINT64) + 1);

	case eRMIPTDouble:		// floating point double
		*(double *)pValue = *(double *)pFrom;
		return(sizeof(double) + 1);

	case eRMIPTVarUint8:		// variable length
		ValLen = *(UINT32 *)pFrom;
		if(ValLen > 0)
			{
			pFrom += sizeof(UINT32);
			*(void **)pValue = pFrom;
			return(sizeof(UINT32) + ValLen + 1);
			}
		*(void **)pValue = (void *)NULL;
		return(sizeof(UINT32) + 1);

	default:
		break;
	}
return(0);
}


UINT64		// returned class instance identifier
CPBErrCorrect::RMI_new(tsThreadPBErrCorrect *pThreadPar, UINT32 Timeout,	// allow at most Timeout seconds for class instantiation
						UINT32 SessionID)									// requesting class instance on this specific session, 0 if on least loaded session
{
int Rslt;
UINT32 JobRslt;
tJobIDEx JobID;
UINT64 ClassInstanceID;
UINT32 ClassMethodID;
UINT32	MaxResponseSize = 0;
time_t Then;
time_t Now;
UINT32 SleepTime;

SleepTime = 50;
Then = time(NULL);
while((Rslt = pThreadPar->pRequester->AddJobRequest(&JobID,pThreadPar->ServiceType,0,eSWMConstruct,0,NULL,0,NULL, SessionID))==0)
	{
	Now = time(NULL);
	if((Now - Then) > Timeout)
		return((UINT64)0);
	CUtility::SleepMillisecs(SleepTime);
	if(SleepTime < 1000)
		SleepTime += 50;
	}
if(Rslt < 1)
	return(0);

SleepTime = 50;
while((Rslt = pThreadPar->pRequester->GetJobResponse(JobID,&ClassInstanceID,&ClassMethodID,&JobRslt,&MaxResponseSize,NULL))==0)
	{
	Now = time(NULL);
	if((Now - Then) > Timeout)
		return((UINT64)0);
	CUtility::SleepMillisecs(SleepTime);
	if(SleepTime < 1000)
		SleepTime += 50;
	}
if(Rslt < 1)		// JobID not accepted, service provider session terminated??
	return((UINT64)0);	
return(ClassInstanceID);
}

volatile UINT32 gDeleteFailed = 0;
void
CPBErrCorrect::RMI_delete(tsThreadPBErrCorrect *pThreadPar, UINT32 Timeout, UINT64 ClassInstanceID) // allow at most Timeout seconds for class deletion
{
int Rslt;
UINT32 JobRslt;
tJobIDEx JobID;
UINT32 ClassMethodID;
UINT32	MaxResponseSize = 0;
time_t Then;
time_t Now;
UINT32 SleepTime;

SleepTime = 50;
Then = time(NULL);
while((Rslt = pThreadPar->pRequester->AddJobRequest(&JobID,pThreadPar->ServiceType,ClassInstanceID,eSWMDestruct,0,NULL,0,NULL))==0)
	{
	Now = time(NULL);
	if((Now - Then) > Timeout)
		return;
	CUtility::SleepMillisecs(SleepTime);
	if(SleepTime < 1000)
		SleepTime += 50;
	}
if(Rslt < 1)		// JobID not accepted, service provider session terminated??
	return;

SleepTime = 50;
while((Rslt = pThreadPar->pRequester->GetJobResponse(JobID,&ClassInstanceID,&ClassMethodID,&JobRslt,&MaxResponseSize,NULL))==0)
	{
	Now = time(NULL);
	if((Now - Then) > Timeout)
		return;
	CUtility::SleepMillisecs(SleepTime);
	if(SleepTime < 1000)
		SleepTime += 50;	
	}
if(Rslt < 1)
#ifdef WIN32
	InterlockedIncrement(&gDeleteFailed);
#else
	__sync_fetch_and_add(&gDeleteFailed,1);
#endif
return;
}

bool 
CPBErrCorrect::RMI_SetScores(tsThreadPBErrCorrect *pThreadPar, UINT32 Timeout,  UINT64 ClassInstanceID,
				int MatchScore,									// score for match
				int MismatchPenalty,							// penalty for mismatch
				int GapOpenPenalty,								// penalty for opening a gap
				int GapExtnPenalty,								// penalty if extending already opened gap
				int DlyGapExtn,									// delayed gap penalties, only apply gap extension penalty if gap at least this length
				int ProgPenaliseGapExtn,						// if non-zero then progressively increment gap extension penalty for gaps of at least this length, 0 to disable, used for PacBio style error profiles
				int AnchorLen)								// identified first and last anchors in alignment to be of at least this length
{
int Rslt;
int RespDataOfs;
UINT32 JobRslt;
tJobIDEx JobID;
UINT32 ClassMethodID;
UINT32	MaxResponseSize = 0;
time_t Then;
time_t Now;
UINT32 SleepTime;

SleepTime = 50;
Then = time(NULL);

RespDataOfs = MarshalReq(pThreadPar->pRMIReqData,eRMIPTInt32,&MatchScore,sizeof(MatchScore));
RespDataOfs += MarshalReq(&pThreadPar->pRMIReqData[RespDataOfs],eRMIPTInt32,&MismatchPenalty,sizeof(MismatchPenalty));
RespDataOfs += MarshalReq(&pThreadPar->pRMIReqData[RespDataOfs],eRMIPTInt32,&GapOpenPenalty,sizeof(GapOpenPenalty));
RespDataOfs += MarshalReq(&pThreadPar->pRMIReqData[RespDataOfs],eRMIPTInt32,&GapExtnPenalty,sizeof(GapExtnPenalty));
RespDataOfs += MarshalReq(&pThreadPar->pRMIReqData[RespDataOfs],eRMIPTInt32,&DlyGapExtn,sizeof(DlyGapExtn));
RespDataOfs += MarshalReq(&pThreadPar->pRMIReqData[RespDataOfs],eRMIPTInt32,&ProgPenaliseGapExtn,sizeof(ProgPenaliseGapExtn));
RespDataOfs += MarshalReq(&pThreadPar->pRMIReqData[RespDataOfs],eRMIPTInt32,&AnchorLen,sizeof(AnchorLen));
while((Rslt = pThreadPar->pRequester->AddJobRequest(&JobID,pThreadPar->ServiceType,ClassInstanceID,eSWMSetScores,0,NULL,RespDataOfs,pThreadPar->pRMIReqData))==0)
	{
	Now = time(NULL);
	if((Now - Then) > Timeout)
		return(false);
	CUtility::SleepMillisecs(SleepTime);
	if(SleepTime < 1000)
		SleepTime += 50;
	}
if(Rslt < 1)		// JobID not accepted, service provider session terminated??
	return(false);

SleepTime = 50;
while((Rslt = pThreadPar->pRequester->GetJobResponse(JobID,&ClassInstanceID,&ClassMethodID,&JobRslt,&MaxResponseSize,NULL))==0)
	{
	Now = time(NULL);
	if((Now - Then) > Timeout)
		return(false);
		CUtility::SleepMillisecs(SleepTime);
	if(SleepTime < 1000)
		SleepTime += 50;
	}
return((JobRslt == 0 || Rslt < 1) ? false : true);
}

bool 
CPBErrCorrect::RMI_SetCPScores(tsThreadPBErrCorrect *pThreadPar, UINT32 Timeout,  UINT64 ClassInstanceID,
				int MatchScore,		// ClassifyPath() score for match
				int MismatchPenalty,	// ClassifyPath() penalty for mismatch
				int GapOpenPenalty,	// ClassifyPath() penalty for opening a gap
				int GapExtnPenalty)	// ClassifyPath() penalty if extending already opened gap
{
int Rslt;
int RespDataOfs;
UINT32 JobRslt;
tJobIDEx JobID;
UINT32 ClassMethodID;
UINT32	MaxResponseSize = 0;
time_t Then;
time_t Now;
UINT32 SleepTime;

SleepTime = 50;
Then = time(NULL);

RespDataOfs = MarshalReq(pThreadPar->pRMIReqData,eRMIPTInt32,&MatchScore,sizeof(MatchScore));
RespDataOfs += MarshalReq(&pThreadPar->pRMIReqData[RespDataOfs],eRMIPTInt32,&MismatchPenalty,sizeof(MismatchPenalty));
RespDataOfs += MarshalReq(&pThreadPar->pRMIReqData[RespDataOfs],eRMIPTInt32,&GapOpenPenalty,sizeof(GapOpenPenalty));
RespDataOfs += MarshalReq(&pThreadPar->pRMIReqData[RespDataOfs],eRMIPTInt32,&GapExtnPenalty,sizeof(GapExtnPenalty));
while((Rslt = pThreadPar->pRequester->AddJobRequest(&JobID,pThreadPar->ServiceType,ClassInstanceID,eSWMSetCPScores,0,NULL,RespDataOfs,pThreadPar->pRMIReqData))==0)
	{
	Now = time(NULL);
	if((Now - Then) > Timeout)
		return(false);
		CUtility::SleepMillisecs(SleepTime);
	if(SleepTime < 1000)
		SleepTime += 50;
	}
if(Rslt < 1)		// JobID not accepted, no service provider available, service provider session terminated??
	return(false);

SleepTime = 50;
while((Rslt = pThreadPar->pRequester->GetJobResponse(JobID,&ClassInstanceID,&ClassMethodID,&JobRslt,&MaxResponseSize,NULL))==0)
	{
	Now = time(NULL);
	if((Now - Then) > Timeout)
		return(false);
		CUtility::SleepMillisecs(SleepTime);
	if(SleepTime < 1000)
		SleepTime += 50;
	}
return((JobRslt == 0 || Rslt < 1) ? false : true);
}


bool 
CPBErrCorrect::RMI_SetMaxInitiatePathOfs(tsThreadPBErrCorrect *pThreadPar, UINT32 Timeout,  UINT64 ClassInstanceID,
			int MaxInitiatePathOfs)	// require SW paths to have started within this many bp (0 to disable) on either the probe or target - effectively anchoring the SW 
{
int Rslt;
int RespDataOfs;
UINT32 JobRslt;
tJobIDEx JobID;
UINT32 ClassMethodID;
UINT32	MaxResponseSize = 0;
time_t Then;
time_t Now;
UINT32 SleepTime;

SleepTime = 50;
Then = time(NULL);

RespDataOfs = MarshalReq(pThreadPar->pRMIReqData,eRMIPTInt32,&MaxInitiatePathOfs,sizeof(MaxInitiatePathOfs));
while((Rslt = pThreadPar->pRequester->AddJobRequest(&JobID,pThreadPar->ServiceType,ClassInstanceID,eSWMSetMaxInitiatePathOfs,0,NULL,RespDataOfs,pThreadPar->pRMIReqData))==0)
	{
	Now = time(NULL);
	if((Now - Then) > Timeout)
		return(false);
		CUtility::SleepMillisecs(SleepTime);
	if(SleepTime < 1000)
		SleepTime += 50;
	}
if(Rslt < 1)		// JobID not accepted, no service provider available, service provider session terminated??
	return(false);

SleepTime = 50;
while((Rslt = pThreadPar->pRequester->GetJobResponse(JobID,&ClassInstanceID,&ClassMethodID,&JobRslt,&MaxResponseSize,NULL))==0)
	{
	Now = time(NULL);
	if((Now - Then) > Timeout)
		return(false);
	CUtility::SleepMillisecs(SleepTime);
	if(SleepTime < 1000)
		SleepTime += 50;
	}
return((JobRslt == 0 || Rslt < 1) ? false : true);
}

bool
CPBErrCorrect::RMI_PreAllocMaxTargLen( tsThreadPBErrCorrect *pThreadPar, UINT32 Timeout,  UINT64 ClassInstanceID,
							 UINT32 MaxTargLen,					// preallocate to process targets of this maximal length
							 UINT32 MaxOverlapLen)			// allocating tracebacks for this maximal expected overlap, 0 if no tracebacks required
{
int Rslt;
int RespDataOfs;
UINT32 JobRslt;
tJobIDEx JobID;
UINT32 ClassMethodID;
UINT32	MaxResponseSize = 0;
time_t Then;
time_t Now;
UINT32 SleepTime;

SleepTime = 50;
Then = time(NULL);

RespDataOfs = MarshalReq(pThreadPar->pRMIReqData,eRMIPTUint32,&MaxTargLen,sizeof(MaxTargLen));
RespDataOfs += MarshalReq(&pThreadPar->pRMIReqData[RespDataOfs],eRMIPTUint32,&MaxOverlapLen,sizeof(MaxOverlapLen));
while((Rslt = pThreadPar->pRequester->AddJobRequest(&JobID,pThreadPar->ServiceType,ClassInstanceID,eSWMPreAllocMaxTargLen,0,NULL,RespDataOfs,pThreadPar->pRMIReqData))==0)
	{
	Now = time(NULL);
	if((Now - Then) > Timeout)
		return(false);
		CUtility::SleepMillisecs(SleepTime);
	if(SleepTime < 1000)
		SleepTime += 50;
	}
if(Rslt < 1)		// JobID not accepted, no service provider available, service provider session terminated??
	return(false);

SleepTime = 50;
while((Rslt = pThreadPar->pRequester->GetJobResponse(JobID,&ClassInstanceID,&ClassMethodID,&JobRslt,&MaxResponseSize,NULL))==0)
	{
	Now = time(NULL);
	if((Now - Then) > Timeout)
		return(false);
		CUtility::SleepMillisecs(SleepTime);
	if(SleepTime < 1000)
		SleepTime += 50;
	}
return((JobRslt == 0 || Rslt < 1) ? false : true);
}

int
CPBErrCorrect::RMI_StartMultiAlignments(tsThreadPBErrCorrect *pThreadPar, UINT32 Timeout,  UINT64 ClassInstanceID,
					int SeqLen,						// probe sequence is this length
					etSeqBase *pProbeSeq,			// probe sequence 
					int Alignments,					// number of pairwise alignments to allocate for
					UINT8 Flags)					// bit 0 set true if probe sequence loaded as a high confidence sequence
{
int Rslt;
int RespDataOfs;
UINT32 JobRslt;
tJobIDEx JobID;
UINT32 ClassMethodID;
UINT32	MaxResponseSize = 0;
time_t Then;
time_t Now;
UINT32 SleepTime;

SleepTime = 50;
Then = time(NULL);

RespDataOfs = MarshalReq(pThreadPar->pRMIReqData,eRMIPTInt32,&SeqLen,sizeof(SeqLen));
RespDataOfs += MarshalReq(&pThreadPar->pRMIReqData[RespDataOfs],eRMIPTVarUint8,pProbeSeq,SeqLen);
RespDataOfs += MarshalReq(&pThreadPar->pRMIReqData[RespDataOfs],eRMIPTInt32,&Alignments,sizeof(Alignments));
RespDataOfs += MarshalReq(&pThreadPar->pRMIReqData[RespDataOfs],eRMIPTUint8,&Flags,sizeof(Flags));
while((Rslt = pThreadPar->pRequester->AddJobRequest(&JobID,pThreadPar->ServiceType,ClassInstanceID,eSWMStartMultiAlignments,0,NULL,RespDataOfs,pThreadPar->pRMIReqData))==0)
	{
	Now = time(NULL);
	if((Now - Then) > Timeout)
		return(-1);
	CUtility::SleepMillisecs(SleepTime);
	if(SleepTime < 1000)
		SleepTime += 50;
	}
if(Rslt < 1)		// JobID not accepted, no service provider available, service provider session terminated??
	return(-1);

SleepTime = 50;
while((Rslt = pThreadPar->pRequester->GetJobResponse(JobID,&ClassInstanceID,&ClassMethodID,&JobRslt,&MaxResponseSize,NULL))==0)
	{
	Now = time(NULL);
	if((Now - Then) > Timeout)
		return(-1);
	CUtility::SleepMillisecs(SleepTime);
	if(SleepTime < 1000)
		SleepTime += 50;
	}
return(Rslt < 1 ? -1 : (int)JobRslt);
}

bool 
CPBErrCorrect::RMI_SetProbe(tsThreadPBErrCorrect *pThreadPar, UINT32 Timeout,  UINT64 ClassInstanceID, 
					UINT32 Len,etSeqBase *pSeq)					// set probe sequence to use in subsequent alignments
{
int Rslt;
int RespDataOfs;
UINT32 JobRslt;
tJobIDEx JobID;
UINT32 ClassMethodID;
UINT32	MaxResponseSize = 0;
time_t Then;
time_t Now;
UINT32 SleepTime;

SleepTime = 50;
Then = time(NULL);
RespDataOfs = MarshalReq(pThreadPar->pRMIReqData,eRMIPTInt32,&Len,sizeof(Len));
RespDataOfs += MarshalReq(&pThreadPar->pRMIReqData[RespDataOfs],eRMIPTVarUint8,pSeq,Len);
while((Rslt = pThreadPar->pRequester->AddJobRequest(&JobID,pThreadPar->ServiceType,ClassInstanceID,eSWMSetProbe,0,NULL,RespDataOfs,pThreadPar->pRMIReqData))==0)
	{
	Now = time(NULL);
	if((Now - Then) > Timeout)
		return(false);
	CUtility::SleepMillisecs(SleepTime);
	if(SleepTime < 1000)
		SleepTime += 50;
	}
if(Rslt < 1)		// JobID not accepted, no service provider available, service provider session terminated??
	return(false);

SleepTime = 50;
while((Rslt = pThreadPar->pRequester->GetJobResponse(JobID,&ClassInstanceID,&ClassMethodID,&JobRslt,&MaxResponseSize,NULL))==0)
	{
	Now = time(NULL);
	if((Now - Then) > Timeout)
		return(false);
	CUtility::SleepMillisecs(SleepTime);
	if(SleepTime < 1000)
		SleepTime += 50;
	}
return((JobRslt == 0 || Rslt < 1) ? false : true);
}

bool 
CPBErrCorrect::RMI_SetTarg(tsThreadPBErrCorrect *pThreadPar, UINT32 Timeout,  UINT64 ClassInstanceID,  
			UINT32 Len,etSeqBase *pSeq)					// set target sequence to use in subsequent alignments
{
int Rslt;
int RespDataOfs;
UINT32 JobRslt;
tJobIDEx JobID;
UINT32 ClassMethodID;
UINT32	MaxResponseSize = 0;
time_t Then;
time_t Now;
UINT32 SleepTime;

SleepTime = 50;
Then = time(NULL);
RespDataOfs = MarshalReq(pThreadPar->pRMIReqData,eRMIPTInt32,&Len,sizeof(Len));
RespDataOfs += MarshalReq(&pThreadPar->pRMIReqData[RespDataOfs],eRMIPTVarUint8,pSeq,Len);
while((Rslt = pThreadPar->pRequester->AddJobRequest(&JobID,pThreadPar->ServiceType,ClassInstanceID,eSWMSetTarg,0,NULL,RespDataOfs,pThreadPar->pRMIReqData))==0)
	{
	Now = time(NULL);
	if((Now - Then) > Timeout)
		return(false);
	CUtility::SleepMillisecs(SleepTime);
	if(SleepTime < 1000)
		SleepTime += 50;
	}
if(Rslt < 1)		// JobID not accepted, no service provider available, service provider session terminated??
	return(false);

SleepTime = 50;
while((Rslt = pThreadPar->pRequester->GetJobResponse(JobID,&ClassInstanceID,&ClassMethodID,&JobRslt,&MaxResponseSize,NULL))==0)
	{
	Now = time(NULL);
	if((Now - Then) > Timeout)
		return(false);
	CUtility::SleepMillisecs(20);
	}
return((JobRslt == 0 || Rslt < 1) ? false : true);
}

int 
CPBErrCorrect::RMI_SetAlignRange(tsThreadPBErrCorrect *pThreadPar, UINT32 Timeout,  UINT64 ClassInstanceID,  
				UINT32 m_ProbeStartRelOfs,	// when aligning then start SW from this probe sequence relative offset
				UINT32 m_TargStartRelOfs, 	// and SW starting from this target sequence relative offset
				UINT32 m_ProbeRelLen,	// and SW with this probe relative length starting from m_ProbeStartRelOfs - if 0 then until end of probe sequence
				UINT32 m_TargRelLen)	// and SW with this target relative length starting from m_TargStartRelOfs - if 0 then until end of target sequence
{
int Rslt;
int RespDataOfs;
UINT32 JobRslt;
tJobIDEx JobID;
UINT32 ClassMethodID;
UINT32	MaxResponseSize = 0;
time_t Then;
time_t Now;
UINT32 SleepTime;

SleepTime = 50;
Then = time(NULL);
RespDataOfs = MarshalReq(pThreadPar->pRMIReqData,eRMIPTUint32,&m_ProbeStartRelOfs,sizeof(m_ProbeStartRelOfs));
RespDataOfs += MarshalReq(&pThreadPar->pRMIReqData[RespDataOfs],eRMIPTUint32,&m_TargStartRelOfs,sizeof(m_TargStartRelOfs));
RespDataOfs += MarshalReq(&pThreadPar->pRMIReqData[RespDataOfs],eRMIPTUint32,&m_ProbeRelLen,sizeof(m_ProbeRelLen));
RespDataOfs += MarshalReq(&pThreadPar->pRMIReqData[RespDataOfs],eRMIPTUint32,&m_TargRelLen,sizeof(m_TargRelLen));
while((Rslt = pThreadPar->pRequester->AddJobRequest(&JobID,pThreadPar->ServiceType,ClassInstanceID,eSWMSetAlignRange,0,NULL,RespDataOfs,pThreadPar->pRMIReqData))==0)
	{
	Now = time(NULL);
	if((Now - Then) > Timeout)
		return(-1);
	CUtility::SleepMillisecs(SleepTime);
	if(SleepTime < 1000)
		SleepTime += 50;
	}
if(Rslt < 1)		// JobID not accepted, no service provider available, service provider session terminated??
	return(-1);

SleepTime = 50;
while((Rslt = pThreadPar->pRequester->GetJobResponse(JobID,&ClassInstanceID,&ClassMethodID,&JobRslt,&MaxResponseSize,NULL))==0)
	{
	Now = time(NULL);
	if((Now - Then) > Timeout)
		return(-1);
	CUtility::SleepMillisecs(SleepTime);
	if(SleepTime < 1000)
		SleepTime += 50;
	}
return(Rslt < 1 ? -1 : (int)JobRslt);
}


tsSSWCell *									// smith-waterman style local alignment, returns highest accumulated exact matches scoring cell
CPBErrCorrect::RMI_Align(tsThreadPBErrCorrect *pThreadPar, UINT32 Timeout,  UINT64 ClassInstanceID, 
				tsSSWCell *pPeakScoreCell,	// optionally also return conventional peak scoring cell
				UINT32 MaxOverlapLen)		// process tracebacks for this maximal expected overlap, 0 if no tracebacks required
{
int Rslt;
int RespDataOfs;
UINT32 JobRslt;
tJobIDEx JobID;
UINT32 ClassMethodID;
UINT32	MaxResponseSize = (sizeof(tsSSWCell) + 20) * 2;
time_t Then;
time_t Now;
UINT32 SleepTime;
bool bPeakScoreCell;
tsSSWCell *pCell;

SleepTime = 50;
Then = time(NULL);
if(pPeakScoreCell != NULL)
	bPeakScoreCell = true;
else
	bPeakScoreCell = false;

memset(&pThreadPar->RMIHighScoreCell,0,sizeof(pThreadPar->RMIHighScoreCell));

RespDataOfs = MarshalReq(pThreadPar->pRMIReqData,eRMIPTBool,&bPeakScoreCell,sizeof(bPeakScoreCell));
RespDataOfs += MarshalReq(&pThreadPar->pRMIReqData[RespDataOfs],eRMIPTUint32,&MaxOverlapLen,sizeof(MaxOverlapLen));
while((Rslt = pThreadPar->pRequester->AddJobRequest(&JobID,pThreadPar->ServiceType,ClassInstanceID,eSWMAlign,0,NULL,RespDataOfs,pThreadPar->pRMIReqData))==0)
	{
	Now = time(NULL);
	if((Now - Then) > Timeout)
		return(NULL);
	CUtility::SleepMillisecs(SleepTime);
	if(SleepTime < 1000)
		SleepTime += 50;
	}
if(Rslt < 1)		// JobID not accepted, no service provider available, service provider session terminated??
	return(NULL);

SleepTime = 50;
while((Rslt = pThreadPar->pRequester->GetJobResponse(JobID,&ClassInstanceID,&ClassMethodID,&JobRslt,&MaxResponseSize,pThreadPar->pRMIRespData))==0)
	{
	Now = time(NULL);
	if((Now - Then) > Timeout)
		return(NULL);
	CUtility::SleepMillisecs(SleepTime);
	if(SleepTime < 1000)
		SleepTime += 50;
	}
if(Rslt < 1 || JobRslt != 0 || MaxResponseSize < sizeof(tsSSWCell))
	return(NULL);

RespDataOfs = UnmarshalResp(sizeof(pThreadPar->RMIHighScoreCell),pThreadPar->pRMIRespData,&pCell);
if(RespDataOfs < sizeof(tsSSWCell))
	return(NULL);
if(pCell != NULL && RespDataOfs >= sizeof(tsSSWCell))
	memcpy(&pThreadPar->RMIHighScoreCell,pCell,sizeof(tsSSWCell));
if(bPeakScoreCell)
	{
	RespDataOfs += UnmarshalResp(sizeof(tsSSWCell),&pThreadPar->pRMIRespData[RespDataOfs],&pCell);
	*pPeakScoreCell = *pCell;
	}
return(&pThreadPar->RMIHighScoreCell);
}

      // methods which combines the functionality of SetTarg, SetAlignRange, Align, ClassifyPath, TracebacksToAlignOps, and AddMultiAlignment into a single method 
int //  -3: timeout waiting for job to complete, -2: parameter errors, -1: class instance no longer exists, 0: currently no available service instance 1: if job accepted and processed
CPBErrCorrect::RMI_CombinedTargAlign(tsThreadPBErrCorrect *pThreadPar, UINT32 Timeout,  UINT64 ClassInstanceID,
						tsCombinedTargAlignPars *pAlignPars, // input alignment parameters
						tsCombinedTargAlignRet *pAlignRet)		// returned alignment results
{
int Rslt;
int RespDataOfs;
UINT32 JobRslt;
tJobIDEx JobID;
UINT32 ClassMethodID;
UINT32	MaxResponseSize = (sizeof(tsCombinedTargAlignRet) + 20) * 2;
tsCombinedTargAlignRet *pTmpAlignRet;
time_t Then;
time_t Now;
UINT32 SleepTime;

SleepTime = 50;
memset(pAlignRet,0,sizeof(tsCombinedTargAlignRet));
Then = time(NULL);
RespDataOfs = MarshalReq(pThreadPar->pRMIReqData,eRMIPTVarUint8,pAlignPars,sizeof(tsCombinedTargAlignPars));
RespDataOfs += MarshalReq(&pThreadPar->pRMIReqData[RespDataOfs],eRMIPTVarUint8,pAlignPars->pTargSeq,pAlignPars->TargSeqLen);
while((Rslt = pThreadPar->pRequester->AddJobRequest(&JobID,pThreadPar->ServiceType,ClassInstanceID,eSWMCombinedTargAlign,0,NULL,RespDataOfs,pThreadPar->pRMIReqData))==0)
	{
	Now = time(NULL);
	if((Now - Then) > Timeout)
		return(-3);
	CUtility::SleepMillisecs(SleepTime);
	if(SleepTime < 1000)
		SleepTime += 50;
	}
if(Rslt < 1)		// -2: parameter errors, -1: class instance no longer exists, 0: currently no available service instance 1: if job accepted
	return(Rslt);

SleepTime = 50;
while((Rslt = pThreadPar->pRequester->GetJobResponse(JobID,&ClassInstanceID,&ClassMethodID,&JobRslt,&MaxResponseSize,pThreadPar->pRMIRespData))==0)
	{
	Now = time(NULL);
	if((Now - Then) > Timeout)
		return(-3);
	CUtility::SleepMillisecs(SleepTime);
	if(SleepTime < 1000)
		SleepTime += 50;
	}
if(Rslt < 1 || MaxResponseSize < sizeof(tsCombinedTargAlignRet))
	return(-2);

pTmpAlignRet = NULL;
RespDataOfs = UnmarshalResp(sizeof(tsCombinedTargAlignRet),pThreadPar->pRMIRespData,&pTmpAlignRet);
if(pTmpAlignRet == NULL || RespDataOfs < sizeof(tsCombinedTargAlignRet))
	return(-2);
*pAlignRet = *pTmpAlignRet;
return(1);
}


int // -3: timeout waiting for job to complete, -2: parameter errors, -1: class instance no longer exists, 0: currently no available service instance 1: if job accepted and processed
CPBErrCorrect::RMI_CombinedTargAlign(tsThreadPBErrCorrect *pThreadPar, UINT32 Timeout,  UINT64 ClassInstanceID,
						UINT8 PMode,              // processing mode: 0 error correct , 1 generate consensus from previously generated multiple alignments, 2  generate overlap detail from previously generated consensus sequences
						UINT32 NumTargSeqs,			// current probe putatively overlaying this many targets
						UINT32 ProbeSeqLen,         // use this probe sequence length 
						UINT8 TargFlags,		    // bit 7 set if target loaded as a high confidence sequence, bits 0..3 is weighting factor to apply when generating consensus bases
						UINT32 TargSeqLen,          // target sequence length
						etSeqBase *pTargSeq,        // target sequence
						UINT32 ProbeStartRelOfs,	// when aligning then start SW from this probe sequence relative offset
						UINT32 TargStartRelOfs, 	// and SW starting from this target sequence relative offset
						UINT32 ProbeRelLen,		    // and SW with this probe relative length starting from m_ProbeStartRelOfs - if 0 then until end of probe sequence
						UINT32 TargRelLen,		    // and SW with this target relative length starting from m_TargStartRelOfs - if 0 then until end of target sequence
						UINT32 OverlapFloat,		// allowing up to this much float on overlaps to account for the PacBio error profile
						UINT32 MaxArtefactDev,		// classify overlaps as artefactual if sliding window of 1Kbp over any overlap deviates by more than this percentage from the overlap mean
						UINT32 MinOverlapLen,       // minimum accepted overlap length
						UINT32 MaxOverlapLen,      // max expected overlap length
						UINT8 *pRetClass,			// returned overlap classification
						tsSSWCell *pRetPeakMatchesCell, // returned peak matches cell
						UINT32 *pRetProbeAlignLength, // probe alignment length
						UINT32 *pRetTargAlignLength, // target alignment length
						bool *pRetbProvOverlapping,  // probe overlapping target
						bool *pRetbProvArtefact,	// set true if overlap classified as artefact
						bool *pRetbProvContained,	// probe was contained
						bool *pRetbAddedMultiAlignment) // added as a multialignment
{
int Rslt;
tsCombinedTargAlignPars AlignPars;
tsCombinedTargAlignRet AlignRet;
memset(&AlignRet,0,sizeof(tsCombinedTargAlignRet));
*pRetClass = 0;
memset(pRetPeakMatchesCell,0,sizeof(tsSSWCell));
*pRetProbeAlignLength = 0;
*pRetTargAlignLength = 0;
*pRetbProvOverlapping = false;
*pRetbProvArtefact = false;
*pRetbProvContained = false;
*pRetbAddedMultiAlignment = false;

AlignPars.PMode = PMode;
AlignPars.NumTargSeqs = NumTargSeqs;
AlignPars.ProbeSeqLen = ProbeSeqLen;
AlignPars.TargFlags = TargFlags;
AlignPars.TargSeqLen = TargSeqLen;
AlignPars.pTargSeq = pTargSeq;
AlignPars.ProbeStartRelOfs = ProbeStartRelOfs;
AlignPars.TargStartRelOfs = TargStartRelOfs;
AlignPars.ProbeRelLen = ProbeRelLen;
AlignPars.TargRelLen = TargRelLen;
AlignPars.OverlapFloat = OverlapFloat;
AlignPars.MaxArtefactDev = MaxArtefactDev;
AlignPars.MinOverlapLen = MinOverlapLen;
AlignPars.MaxOverlapLen = MaxOverlapLen;

if((Rslt = RMI_CombinedTargAlign(pThreadPar,Timeout,ClassInstanceID,&AlignPars,&AlignRet)) > 0)
	{
	*pRetClass = AlignRet.Class;
	*pRetPeakMatchesCell = AlignRet.PeakMatchesCell;
	*pRetProbeAlignLength = AlignRet.ProbeAlignLength; 
	*pRetTargAlignLength = AlignRet.TargAlignLength;
	*pRetbProvOverlapping = AlignRet.Flags & 0x01 ? true : false;
	*pRetbProvArtefact = AlignRet.Flags & 0x02 ? true : false;
	*pRetbProvContained = AlignRet.Flags & 0x04 ? true : false;
	*pRetbAddedMultiAlignment = AlignRet.Flags & 0x08 ? true : false; 
	}
return(Rslt);
}

int												// attempting to determine if path is artfact resulting from aligning to a paralogous fragment
CPBErrCorrect::RMI_ClassifyPath(tsThreadPBErrCorrect *pThreadPar, UINT32 Timeout,  UINT64 ClassInstanceID,
			int MaxArtefactDev,			// classify path as artefactual if sliding window of 500bp over any overlap deviates by more than this percentage from the overlap mean
			UINT32 ProbeStartOfs,			// alignment starts at this probe sequence offset (1..n)
			UINT32 ProbeEndOfs,				// alignment ends at this probe sequence offset
			UINT32 TargStartOfs,			// alignment starts at this target sequence offset (1..n)
			UINT32 TargEndOfs)				// alignment ends at this target sequence offset
{
int Rslt;
int RespDataOfs;
UINT32 JobRslt;
tJobIDEx JobID;
UINT32 ClassMethodID;
UINT32	MaxResponseSize = 0;
time_t Then;
time_t Now;
UINT32 SleepTime;

SleepTime = 50;
Then = time(NULL);
RespDataOfs = MarshalReq(pThreadPar->pRMIReqData,eRMIPTInt32,&MaxArtefactDev,sizeof(MaxArtefactDev));
RespDataOfs += MarshalReq(&pThreadPar->pRMIReqData[RespDataOfs],eRMIPTUint32,&ProbeStartOfs,sizeof(ProbeStartOfs));
RespDataOfs += MarshalReq(&pThreadPar->pRMIReqData[RespDataOfs],eRMIPTUint32,&ProbeEndOfs,sizeof(ProbeEndOfs));
RespDataOfs += MarshalReq(&pThreadPar->pRMIReqData[RespDataOfs],eRMIPTUint32,&TargStartOfs,sizeof(TargStartOfs));
RespDataOfs += MarshalReq(&pThreadPar->pRMIReqData[RespDataOfs],eRMIPTUint32,&TargEndOfs,sizeof(TargEndOfs));
while((Rslt = pThreadPar->pRequester->AddJobRequest(&JobID,pThreadPar->ServiceType,ClassInstanceID,eSWMClassifyPath,0,NULL,RespDataOfs,pThreadPar->pRMIReqData))==0)
	{
	Now = time(NULL);
	if((Now - Then) > Timeout)
		return(-1);
	CUtility::SleepMillisecs(SleepTime);
	if(SleepTime < 1000)
		SleepTime += 50;
	}
if(Rslt < 1)		// JobID not accepted, no service provider available, service provider session terminated??
	return(-1);

SleepTime = 50;
while((Rslt = pThreadPar->pRequester->GetJobResponse(JobID,&ClassInstanceID,&ClassMethodID,&JobRslt,&MaxResponseSize,pThreadPar->pRMIRespData))==0)
	{
	Now = time(NULL);
	if((Now - Then) > Timeout)
		return(-1);
	CUtility::SleepMillisecs(SleepTime);
	if(SleepTime < 1000)
		SleepTime += 50;
	}
return(Rslt < 1 ? -1 : (int)JobRslt);
}

int												// number of alignment ops generated
CPBErrCorrect::RMI_TracebacksToAlignOps(tsThreadPBErrCorrect *pThreadPar, UINT32 Timeout,  UINT64 ClassInstanceID,
			UINT32 ProbeStartOfs,			// alignment starts at this probe sequence offset (1..n)
			UINT32 ProbeEndOfs,				// alignment ends at this probe sequence offset
			UINT32 TargStartOfs,			// alignment starts at this target sequence offset (1..n)
			UINT32 TargEndOfs,				// alignment ends at this target sequence offset
			tMAOp **ppAlignOps)				// optionally return ptr to alignment operations
{
int Rslt;
int RespDataOfs;
UINT32 JobRslt;
tJobIDEx JobID;
UINT32 ClassMethodID;
UINT32	MaxResponseSize;
time_t Then;
time_t Now;
UINT32 SleepTime;
bool bAlignOps;

SleepTime = 50;
if(ppAlignOps != NULL)
	{
	MaxResponseSize = 1000000;
	bAlignOps = true;
	}
else
	{
	MaxResponseSize = 0;
	bAlignOps = false;
	}
Then = time(NULL);
RespDataOfs = MarshalReq(pThreadPar->pRMIReqData,eRMIPTUint32,&ProbeStartOfs,sizeof(ProbeStartOfs));
RespDataOfs += MarshalReq(&pThreadPar->pRMIReqData[RespDataOfs],eRMIPTUint32,&ProbeEndOfs,sizeof(ProbeEndOfs));
RespDataOfs += MarshalReq(&pThreadPar->pRMIReqData[RespDataOfs],eRMIPTUint32,&TargStartOfs,sizeof(TargStartOfs));
RespDataOfs += MarshalReq(&pThreadPar->pRMIReqData[RespDataOfs],eRMIPTUint32,&TargEndOfs,sizeof(TargEndOfs));
RespDataOfs += MarshalReq(&pThreadPar->pRMIReqData[RespDataOfs],eRMIPTBool,&bAlignOps,sizeof(bAlignOps));
while((Rslt = pThreadPar->pRequester->AddJobRequest(&JobID,pThreadPar->ServiceType,ClassInstanceID,eSWMTracebacksToAlignOps,0,NULL,RespDataOfs,pThreadPar->pRMIReqData))==0)
	{
	Now = time(NULL);
	if((Now - Then) > Timeout)
		return(-1);
	CUtility::SleepMillisecs(SleepTime);
	if(SleepTime < 1000)
		SleepTime += 50;
	}
if(Rslt < 1)		// JobID not accepted, no service provider available, service provider session terminated??
	return(-1);

SleepTime = 50;
while((Rslt = pThreadPar->pRequester->GetJobResponse(JobID,&ClassInstanceID,&ClassMethodID,&JobRslt,&MaxResponseSize,pThreadPar->pRMIRespData))==0)
	{
	Now = time(NULL);
	if((Now - Then) > Timeout)
		return(-1);
	CUtility::SleepMillisecs(SleepTime);
	if(SleepTime < 1000)
		SleepTime += 50;
	}
if(Rslt < 1 || (bAlignOps && MaxResponseSize < 1))
	return(-1);
if(bAlignOps)
	{
	RespDataOfs = UnmarshalResp(sizeof(pThreadPar->RMIHighScoreCell),pThreadPar->pRMIRespData,pThreadPar->pRMIBuffer);
	*ppAlignOps = (tMAOp *)pThreadPar->pRMIBuffer;
	}
return((int)JobRslt);
}

int
CPBErrCorrect::RMI_AddMultiAlignment(tsThreadPBErrCorrect *pThreadPar, UINT32 Timeout,  UINT64 ClassInstanceID,
				UINT32 ProbeStartOfs,			// alignment starts at this probe sequence offset (1..n)
				UINT32 ProbeEndOfs,				// alignment ends at this probe sequence offset inclusive
				UINT32 TargStartOfs,			// alignment starts at this target sequence offset (1..n)
				UINT32 TargEndOfs,				// alignment ends at this target sequence offset inclusive
				UINT32 TargSeqLen,				// target sequence length
				etSeqBase *pTargSeq,			// alignment target sequence
				UINT8 Flags)                    // bit 7 set if target loaded as a high confidence sequence, bits 0..3 is weighting factor to apply when generating consensus bases
{
int Rslt;
int RespDataOfs;
UINT32 JobRslt;
tJobIDEx JobID;
UINT32 ClassMethodID;
UINT32	MaxResponseSize = 0;
time_t Then;
time_t Now;
UINT32 SleepTime;

SleepTime = 50;
Then = time(NULL);
RespDataOfs = MarshalReq(pThreadPar->pRMIReqData,eRMIPTUint32,&ProbeStartOfs,sizeof(ProbeStartOfs));
RespDataOfs += MarshalReq(&pThreadPar->pRMIReqData[RespDataOfs],eRMIPTUint32,&ProbeEndOfs,sizeof(ProbeEndOfs));
RespDataOfs += MarshalReq(&pThreadPar->pRMIReqData[RespDataOfs],eRMIPTUint32,&TargStartOfs,sizeof(TargStartOfs));
RespDataOfs += MarshalReq(&pThreadPar->pRMIReqData[RespDataOfs],eRMIPTUint32,&TargEndOfs,sizeof(TargEndOfs));
RespDataOfs += MarshalReq(&pThreadPar->pRMIReqData[RespDataOfs],eRMIPTUint32,&TargEndOfs,sizeof(TargSeqLen));
RespDataOfs += MarshalReq(&pThreadPar->pRMIReqData[RespDataOfs],eRMIPTVarUint8,pTargSeq,TargSeqLen);
RespDataOfs += MarshalReq(&pThreadPar->pRMIReqData[RespDataOfs],eRMIPTUint8,&Flags,sizeof(Flags));
while((Rslt = pThreadPar->pRequester->AddJobRequest(&JobID,pThreadPar->ServiceType,ClassInstanceID,eSWMAddMultiAlignment,0,NULL,RespDataOfs,pThreadPar->pRMIReqData))==0)
	{
	Now = time(NULL);
	if((Now - Then) > Timeout)
		return(-1);
	CUtility::SleepMillisecs(SleepTime);
	if(SleepTime < 1000)
		SleepTime += 50;
	}
if(Rslt < 1)		// JobID not accepted, no service provider available, service provider session terminated??
	return(-1);

SleepTime = 50;
while((Rslt = pThreadPar->pRequester->GetJobResponse(JobID,&ClassInstanceID,&ClassMethodID,&JobRslt,&MaxResponseSize,pThreadPar->pRMIRespData))==0)
	{
	Now = time(NULL);
	if((Now - Then) > Timeout)
		return(-1);
	CUtility::SleepMillisecs(SleepTime);
	if(SleepTime < 1000)
		SleepTime += 50;
	}
return(Rslt < 1 ? -1 : (int)JobRslt);
}

int
CPBErrCorrect::RMI_GenMultialignConcensus(tsThreadPBErrCorrect *pThreadPar, UINT32 Timeout,  UINT64 ClassInstanceID)
{
int Rslt;
UINT32 JobRslt;
tJobIDEx JobID;
UINT32 ClassMethodID;
UINT32	MaxResponseSize = 0;
time_t Then;
time_t Now;
UINT32 SleepTime;

SleepTime = 50;
Then = time(NULL);
while((Rslt = pThreadPar->pRequester->AddJobRequest(&JobID,pThreadPar->ServiceType,ClassInstanceID,eSWMGenMultialignConcensus,0,NULL,0,NULL))==0)
	{
	Now = time(NULL);
	if((Now - Then) > Timeout)
		return(-1);
	CUtility::SleepMillisecs(SleepTime);
	if(SleepTime < 1000)
		SleepTime += 50;
	}
if(Rslt < 1)		// JobID not accepted, no service provider available, service provider session terminated??
	return(-1);

SleepTime = 50;
while((Rslt = pThreadPar->pRequester->GetJobResponse(JobID,&ClassInstanceID,&ClassMethodID,&JobRslt,&MaxResponseSize,pThreadPar->pRMIRespData))==0)
	{
	Now = time(NULL);
	if((Now - Then) > Timeout)
		return(-1);
	CUtility::SleepMillisecs(SleepTime);
	if(SleepTime < 1000)
		SleepTime += 50;
	}
return(Rslt < 1 ? -1 : (int)JobRslt);
}

int      // total number of returned chars in pszBuffer for the textual representation of the error corrected consensus sequence (could be multiple consensus sequences)
CPBErrCorrect::RMI_MAlignCols2fasta(tsThreadPBErrCorrect *pThreadPar, UINT32 Timeout,  UINT64 ClassInstanceID,
			UINT32 ProbeID,	// identifies sequence which was used as the probe when determining the multialignments
			int MinConf,				// sequence bases averaged over 100bp must be of at least this confidence (0..9)
			int MinLen,					// and sequence lengths must be of at least this length 
			UINT32 BuffSize,			// buffer allocated to hold at most this many chars
			char *pszBuffer)			// output error corrected sequences to this buffer
{
int Rslt;
int RespDataOfs;
UINT32 JobRslt;
tJobIDEx JobID;
UINT32 ClassMethodID;
UINT32	MaxResponseSize;
char *pszAlignRow;
time_t Then;
time_t Now;
UINT32 SleepTime;

SleepTime = 50;
if(BuffSize > pThreadPar->RMIBufferSize)
	BuffSize = pThreadPar->RMIBufferSize;
MaxResponseSize = BuffSize;
Then = time(NULL);
RespDataOfs = MarshalReq(pThreadPar->pRMIReqData,eRMIPTUint32,&ProbeID,sizeof(ProbeID));
RespDataOfs += MarshalReq(&pThreadPar->pRMIReqData[RespDataOfs],eRMIPTInt32,&MinConf,sizeof(MinConf));
RespDataOfs += MarshalReq(&pThreadPar->pRMIReqData[RespDataOfs],eRMIPTInt32,&MinLen,sizeof(MinLen));
RespDataOfs += MarshalReq(&pThreadPar->pRMIReqData[RespDataOfs],eRMIPTUint32,&BuffSize,sizeof(BuffSize));

while((Rslt = pThreadPar->pRequester->AddJobRequest(&JobID,pThreadPar->ServiceType,ClassInstanceID,eSWMMAlignCols2fasta,0,NULL,RespDataOfs,pThreadPar->pRMIReqData))==0)
	{
	Now = time(NULL);
	if((Now - Then) > Timeout)
		return(-1);
	CUtility::SleepMillisecs(SleepTime);
	if(SleepTime < 1000)
		SleepTime += 50;
	}
if(Rslt < 1)		// JobID not accepted, no service provider available, service provider session terminated??
	return(-1);

SleepTime = 50;
while((Rslt = pThreadPar->pRequester->GetJobResponse(JobID,&ClassInstanceID,&ClassMethodID,&JobRslt,&MaxResponseSize,pThreadPar->pRMIRespData))==0)
	{
	Now = time(NULL);
	if((Now - Then) > Timeout)
		return(-1);
	CUtility::SleepMillisecs(SleepTime);
	if(SleepTime < 1000)
		SleepTime += 50;
	}
if(Rslt < 1)
	return(-1);
if(JobRslt > 0)
	{
	RespDataOfs = UnmarshalResp(JobRslt,pThreadPar->pRMIRespData,(UINT8 *)&pszAlignRow);
	strncpy(pszBuffer,pszAlignRow,JobRslt);
	pszBuffer[JobRslt] = '\0';
	}
return(JobRslt);
}

int      // total number of returned chars in pszBuffer for the textual representation of the multialignment or -1 if job not accepted 
CPBErrCorrect::RMI_MAlignCols2MFA(tsThreadPBErrCorrect *pThreadPar, UINT32 Timeout,  UINT64 ClassInstanceID,
					UINT32 ProbeID,		// identifies sequence which was used as the probe when determining the multialignments
				UINT32 BuffSize,	// buffer allocated to hold at most this many chars
				char *pszBuffer)	// output multialignment textual representation to this buffer
{
int Rslt;
int RespDataOfs;
UINT32 JobRslt;
tJobIDEx JobID;
UINT32 ClassMethodID;
UINT32	MaxResponseSize;
char *pszAlignRow;
time_t Then;
time_t Now;
UINT32 SleepTime;

SleepTime = 50;
if(BuffSize > pThreadPar->RMIBufferSize)
	{
	BuffSize = pThreadPar->RMIBufferSize;
	}
MaxResponseSize = BuffSize;
Then = time(NULL);
RespDataOfs = MarshalReq(pThreadPar->pRMIReqData,eRMIPTUint32,&ProbeID,sizeof(ProbeID));
RespDataOfs += MarshalReq(&pThreadPar->pRMIReqData[RespDataOfs],eRMIPTUint32,&BuffSize,sizeof(BuffSize));

while((Rslt = pThreadPar->pRequester->AddJobRequest(&JobID,pThreadPar->ServiceType,ClassInstanceID,eSWMMAlignCols2MFA,0,NULL,RespDataOfs,pThreadPar->pRMIReqData))==0)
	{
	Now = time(NULL);
	if((Now - Then) > Timeout)
		return(-1);
	CUtility::SleepMillisecs(SleepTime);
	if(SleepTime < 1000)
		SleepTime += 50;
	}
if(Rslt < 1)		// JobID not accepted, no service provider available, service provider session terminated??
	return(-1);

SleepTime = 50;
while((Rslt = pThreadPar->pRequester->GetJobResponse(JobID,&ClassInstanceID,&ClassMethodID,&JobRslt,&MaxResponseSize,pThreadPar->pRMIRespData))==0)
	{
	Now = time(NULL);
	if((Now - Then) > Timeout)
		return(-1);
	CUtility::SleepMillisecs(SleepTime);
	if(SleepTime < 1000)
		SleepTime += 50;
	}
if(Rslt < 1)
	return(-1);
if(JobRslt > 0)
	{
	RespDataOfs = UnmarshalResp(JobRslt,pThreadPar->pRMIRespData,(UINT8 *)&pszAlignRow);
	strncpy(pszBuffer,pszAlignRow,JobRslt);
	pszBuffer[JobRslt] = '\0';
	}
return(JobRslt);
}

// SortLenDescending
// Sort scaffolding nodes by length descending with entry identifiers as tie breaker
int
CPBErrCorrect::SortLenDescending(const void *arg1, const void *arg2)
{
tsPBEScaffNode *pEl1 = (tsPBEScaffNode *)arg1;
tsPBEScaffNode *pEl2 = (tsPBEScaffNode *)arg2;

if(pEl1->SeqLen < pEl2->SeqLen)
	return(1);
if(pEl1->SeqLen > pEl2->SeqLen)
	return(-1);
if(pEl1->EntryID < pEl2->EntryID)	
	return(-1);
if(pEl1->EntryID > pEl2->EntryID)
	return(1);
return(0);
}


// SortCoreHitsByTargNodeID
// Sort core hits by ProbeNodeID.TargNodeID.TargOfs.ProbeOfs.flgRevCpl ascending
int
CPBErrCorrect::SortCoreHitsByTargProbeOfs(const void *arg1, const void *arg2)
{
tsPBECoreHit *pEl1 = (tsPBECoreHit *)arg1;
tsPBECoreHit *pEl2 = (tsPBECoreHit *)arg2;

if(pEl1->ProbeNodeID < pEl2->ProbeNodeID)
	return(-1);
if(pEl1->ProbeNodeID > pEl2->ProbeNodeID)
	return(1);
if(pEl1->TargNodeID < pEl2->TargNodeID)
	return(-1);
if(pEl1->TargNodeID > pEl2->TargNodeID)
	return(1);
if(pEl1->TargOfs < pEl2->TargOfs)	
	return(-1);
if(pEl1->TargOfs > pEl2->TargOfs)
	return(1);
if(pEl1->ProbeOfs < pEl2->ProbeOfs)	
	return(-1);
if(pEl1->ProbeOfs > pEl2->ProbeOfs)
	return(1);
if(pEl1->flgRevCpl != pEl2->flgRevCpl)
	return(1);
return(0);
}

// SortCoreHitsByProbeNodeID
// Sort core hits by ProbeNodeID.TargNodeID.ProbeOfs.TargOfs.flgRevCpl ascending
int
CPBErrCorrect::SortCoreHitsByProbeTargOfs(const void *arg1, const void *arg2)
{
tsPBECoreHit *pEl1 = (tsPBECoreHit *)arg1;
tsPBECoreHit *pEl2 = (tsPBECoreHit *)arg2;

if(pEl1->ProbeNodeID < pEl2->ProbeNodeID)
	return(-1);
if(pEl1->ProbeNodeID > pEl2->ProbeNodeID)
	return(1);
if(pEl1->TargNodeID < pEl2->TargNodeID)
	return(-1);
if(pEl1->TargNodeID > pEl2->TargNodeID)
	return(1);
if(pEl1->ProbeOfs < pEl2->ProbeOfs)	
	return(-1);
if(pEl1->ProbeOfs > pEl2->ProbeOfs)
	return(1);
if(pEl1->TargOfs < pEl2->TargOfs)	
	return(-1);
if(pEl1->TargOfs > pEl2->TargOfs)
	return(1);
if(pEl1->flgRevCpl != pEl2->flgRevCpl)
	return(1);
return(0);
}

// SortCoreHitsDescending
// Sort target core hits by number of hits descending
int
CPBErrCorrect::SortCoreHitsDescending(const void *arg1, const void *arg2)
{
sPBECoreHitCnts *pEl1 = (sPBECoreHitCnts *)arg1;
sPBECoreHitCnts *pEl2 = (sPBECoreHitCnts *)arg2;

UINT32 El1NumHits;
UINT32 El2NumHits;
El1NumHits = max(pEl1->NumSHits,pEl1->NumAHits);
El2NumHits = max(pEl2->NumSHits,pEl2->NumAHits);
if(El1NumHits > El2NumHits)
	return(-1);
if(El1NumHits < El2NumHits)
	return(1);
return(0);
}





