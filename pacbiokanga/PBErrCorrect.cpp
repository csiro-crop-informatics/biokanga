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

#include "pacbiokanga.h"

#include "../libbiokanga/bgzf.h"
#include "SSW.h"
#include "pacbiocommon.h"
#include "SeqStore.h"
#include "AssembGraph.h"
#include "PBErrCorrect.h"

int
ProcPacBioErrCorrect(etPBPMode PMode,		// processing mode
		int SampleRate,				// sample input sequences at this rate (1..100)
		int DeltaCoreOfs,			// offset by this many bp the core windows of coreSeqLen along the probe sequence when checking for overlaps
		int MaxSeedCoreDepth,		// only further process a seed core if there are no more than this number of matching cores in all targeted sequences
		int MinSeedCoreLen,			// use seed cores of this length when identifying putative overlapping scaffold sequences
		int MinNumSeedCores,        // require at least this many seed cores between overlapping scaffold sequences
		int SWMatchScore,			// score for matching bases (0..50)
		int SWMismatchPenalty,		// mismatch penalty (-50..0)
		int SWGapOpenPenalty,		// gap opening penalty (-50..0)
		int SWGapExtnPenalty,		// gap extension penalty (-50..0)
		int SWProgExtnPenaltyLen,	// progressive gap scoring then only apply gap extension score if gap at least this length (0..63) - use if aligning PacBio
		int MinPBSeqLen,			// only accepting PacBio reads of at least this length (defaults to 15Kbp)
		int MinPBSeqOverlap,		// any overlap of a PacBio onto a target PacBio must be of at least this many bp to be considered for contributing towards error correction (defaults to 5Kbp) 
		int MaxArtefactDev,			// classify overlaps as artefactual if sliding window of 1Kbp over any overlap deviates by more than this percentage from the overlap mean
		int MinHCSeqLen,			// only accepting hiconfidence reads of at least this length (defaults to 1Kbp)
		int MinHCSeqOverlap,		// any overlap of a hiconfidence read onto a target PacBio read must be of at least this many bp to be considered for contributing towards error correction (defaults to 1Kbp) 
		int MinErrCorrectLen,		// error corrected sequences must be at least this minimum length
		int MinConcScore,			// error corrected sequences trimmed until mean 100bp concensus score is at least this threshold
		int NumPacBioFiles,			// number of input pacbio file specs
		char *pszPacBioFiles[],		// input pacbio files
		int NumHiConfFiles,			// number of input hiconfidence file specs
		char *pszHiConfFiles[],		// input hiconfidence files		
	    char *pszErrCorFile,		// name of file into which write error corrected sequences
		char *pszMultiAlignFile,	// name of file into which write multiple alignments
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
int MinSeedCoreLen;			// use seed cores of this length when identifying putative overlapping sequences
int MinNumSeedCores;        // require at least this many seed cores between overlapping sequences before attempting SW
int DeltaCoreOfs;			// offset by this many bp the core windows of coreSeqLen along the probe sequence when checking for overlaps
int MaxSeedCoreDepth;		// only further process a seed core if there are no more than this number of matching cores in all targeted sequences

int SWMatchScore;			// score for matching bases (0..50)
int SWMismatchPenalty;		// mismatch penalty (0..50)
int SWGapOpenPenalty;		// gap opening penalty (0..50)
int SWGapExtnPenalty;		// gap extension penalty (0..50)
int SWProgExtnPenaltyLen;	// progressive gap scoring then only apply gap extension score if gap at least this length (0..63) - use if aligning PacBio

int MinPBSeqLen;			// only accepting PacBio reads of at least this length (defaults to 15Kbp)
int MinPBSeqOverlap;		// any overlap of a PacBio onto a target PacBio must be of at least this many bp to be considered for contributing towards error correction (defaults to 5Kbp) 
int MinHCSeqLen;			// only accepting hiconfidence reads of at least this length (defaults to 1Kbp)
int MinHCSeqOverlap;		// any overlap of a hiconfidence read onto a target PacBio read must be of at least this many bp to be considered for contributing towards error correction (defaults to 1Kbp) 
int MinConcScore;			// error corrected sequences trimmed until mean 50bp concensus score is at least this threshold
int MinErrCorrectLen;		// error corrected and sequence trimmed sequences must be of at least this length
int MaxArtefactDev;			// classify overlaps as artefactual if sliding window of 1000bp over any overlap deviates by more than this percentage from the overlap mean

int NumPacBioFiles;						// number of input pacbio file specs
char *pszPacBioFiles[cMaxInFileSpecs];	// input pacbio files
int NumHiConfFiles;						// number of input hiconfidence file specs
char *pszHiConfFiles[cMaxInFileSpecs];	// input hiconfidence files

char szOutFile[_MAX_PATH];				// where to write (ePBPMConsensus or ePBPMErrCorrect) error corrected sequences or (ePBPMOverlapDetail) overlap detail
char szOutMAFile[_MAX_PATH];			// where to write optional multialignments

int SampleRate;				// sample input sequences at this rate
int NumberOfProcessors;		// number of installed CPUs
int NumThreads;				// number of threads (0 defaults to number of CPUs)

char szSQLiteDatabase[_MAX_PATH];	// results summaries to this SQLite file
char szExperimentName[cMaxDatasetSpeciesChrom+1];			// experiment name
char szExperimentDescr[1000];		// describes experiment

//
struct arg_lit  *help    = arg_lit0("h","help",                 "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to screen and logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"diagnostics log file");

struct arg_int *pmode = arg_int0("m","pmode","<int>",				"processing mode - 0 error correct, 1 consensus sequence only, 2 scaffold overlap detail only");
struct arg_int *minpbseqlen = arg_int0("l","minpbseqlen","<int>",		"minimum individual PacBio sequence length (default 15000, range 500 to 100000)");
struct arg_int *minpbseqovl = arg_int0("L","minpbseqovl","<int>",		"minimum PacBio overlap onto PacBio length required (default 5000, range 500 to 100000)");

struct arg_int *minhcseqlen = arg_int0("p","minhcseqlen","<int>",		"minimum individual high confidence sequence length (default 1000, range 250 to 100000)");
struct arg_int *minhcseqovl = arg_int0("P","minhcseqovl","<int>",		"minimum high confidence sequence overlap onto PacBio length required (default 500, range 250 to 100000)");

struct arg_int *minseedcorelen = arg_int0("c","seedcorelen","<int>",			"use seed cores of this length when identifying putative overlapping sequences (default 12, range 10 to 50)");
struct arg_int *minseedcores = arg_int0("C","minseedcores","<int>",			"require at least this many accepted seed cores between overlapping sequences to use SW (default 10, range 1 to 50)");

struct arg_int *deltacoreofs = arg_int0("d","deltacoreofs","<int>",				"offset cores (default 2, range 1 to 10)");
struct arg_int *maxcoredepth = arg_int0("D","maxcoredepth","<int>",				"explore cores of less than this maximum depth (default 10000, range 1000 to 20000)");

struct arg_int *matchscore = arg_int0("x","matchscore","<int>",					"SW score for matching bases (default 3, range 1 to 50)");
struct arg_int *mismatchpenalty = arg_int0("X","mismatchpenalty","<int>",		"SW mismatch penalty (default 7, range 1 to 50)");
struct arg_int *gapopenpenalty = arg_int0("y","gapopenpenalty","<int>",			"SW gap opening penalty (default 4, range 1 to 50)");
struct arg_int *gapextnpenalty = arg_int0("Y","gapextnpenalty","<int>",			"SW gap extension penalty (default 1, range 1 to 50)");
struct arg_int *progextnpenaltylen = arg_int0("z","progextnpenaltylen","<int>",	"SW gap extension penalty only applied for gaps of at least this number of bases (default 2, range 1 to 63)");

struct arg_int *minconcscore = arg_int0("s","minconcscore","<int>",			     "error corrected sequences trimmed until mean 50bp concensus score is at least this threshold (default 3, range 0 to 9)");
struct arg_int *minerrcorrectlen = arg_int0("S","minerrcorrectlen","<int>",		 "error corrected and trimmed sequences must be at least this minimum length (default 5000, range 500 to 20000)");
struct arg_int *maxartefactdev = arg_int0("a","artefactdev","<int>",			 "classify overlaps as artefactual if 1Kbp window score deviates by more than this percentage from complete overlap mean (0 to disable, range 1 to 25)");

struct arg_file *hiconffiles = arg_filen("I","hiconffile","<file>",0,cMaxInFileSpecs,		"optional, names of input files containing higher confidence reads or sequences to be used in error correcton of PacBio reads");
struct arg_file *pacbiofiles = arg_filen("i","pacbiofile","<file>",1,cMaxInFileSpecs,		"names of input files containing PacBio sequences to be error corrected");
struct arg_file *outfile = arg_file1("o","out","<file>",									"output error corrected PacBio reads to this file");
struct arg_file *mafile = arg_file0("O","mafile","<file>",						"optional, output multialignments to this file, caution can grow very large");
struct arg_file *scaffovrlapsfile = arg_file0("e","scaffovrlapsfile","<file>",	"optional, output scaffolding overlap detail to this file");

struct arg_int *samplerate = arg_int0("R","samplerate","<int>",					"sample input sequences at this rate (default 100, range 1 to 100)");

struct arg_int *threads = arg_int0("T","threads","<int>",						"number of processing threads 0..128 (defaults to 0 which sets threads to number of CPU cores)");

struct arg_file *summrslts = arg_file0("q","sumrslts","<file>",				"Output results summary to this SQLite3 database file");
struct arg_str *experimentname = arg_str0("w","experimentname","<str>",		"experiment name SQLite3 database file");
struct arg_str *experimentdescr = arg_str0("W","experimentdescr","<str>",	"experiment description SQLite3 database file");

struct arg_end *end = arg_end(200);

void *argtable[] = {help,version,FileLogLevel,LogFile,
					pmode,minseedcorelen,minseedcores,deltacoreofs,maxcoredepth,
					matchscore,mismatchpenalty,gapopenpenalty,gapextnpenalty,progextnpenaltylen,
					minpbseqlen,minpbseqovl,minhcseqlen,minhcseqovl,minconcscore,minerrcorrectlen,maxartefactdev,
					summrslts,pacbiofiles,hiconffiles,experimentname,experimentdescr,
					outfile,mafile,scaffovrlapsfile,samplerate,threads,
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
	if(PMode < ePBPMErrCorrect || PMode > ePBPMOverlapDetail)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Processing mode '-m%d' must be in range 0..%d",PMode,ePBPMOverlapDetail);
		return(1);
		}

	MinSeedCoreLen = cDfltSeedCoreLen;
	MinNumSeedCores = cDfltNumSeedCores;
	DeltaCoreOfs = cDfltDeltaCoreOfs;
	MaxSeedCoreDepth = cDfltMaxSeedCoreDepth;
	SWMatchScore = cDfltSWMatchScore;
	SWMismatchPenalty = abs(cDfltSWMismatchPenalty);
	SWGapOpenPenalty = abs(cDfltSWGapOpenPenalty);
	SWGapExtnPenalty = abs(cDfltSWGapExtnPenalty);
	SWProgExtnPenaltyLen = cDfltSWProgExtnLen;
	MinPBSeqLen = cDfltMinPBSeqLen;
	MinPBSeqOverlap = cDfltMinErrCorrectLen; 
	MinHCSeqLen = 1000;
	MinHCSeqOverlap = 500;
	MaxArtefactDev = cDfltMaxArtefactDev;
	NumPacBioFiles = 0;
	NumHiConfFiles = 0;
	szOutMAFile[0] = '\0';

	SampleRate = samplerate->count ? samplerate->ival[0] : 100;
	if(SampleRate < 1 || SampleRate > 100)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Input sequence sampling rate '-R%d' must be in range 1..100",SampleRate);
		return(1);
		}

	if(PMode == ePBPMErrCorrect)
		{
		MinSeedCoreLen = minseedcorelen->count ? minseedcorelen->ival[0] : cDfltSeedCoreLen;
		if(MinSeedCoreLen < cMinSeedCoreLen || MinSeedCoreLen > cMaxSeedCoreLen)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: seed core length '-c%d' must be in range %d..%d",MinSeedCoreLen,cMinSeedCoreLen,cMaxSeedCoreLen);
			return(1);
			}

		MinNumSeedCores = minseedcores->count ? minseedcores->ival[0] : cDfltNumSeedCores;
		if(MinNumSeedCores < cMinNumSeedCores || MinNumSeedCores > cMaxNumSeedCores)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Minimum number of accepted seed cores '-C%d' must be in range %d..%d",MinNumSeedCores,cMinNumSeedCores,cMaxNumSeedCores);
			return(1);
			}

		DeltaCoreOfs = deltacoreofs->count ? deltacoreofs->ival[0] : cDfltDeltaCoreOfs;
		if(DeltaCoreOfs < 1 || DeltaCoreOfs > 10)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: offset seed cores '-d%d' must be in range 1..10",DeltaCoreOfs);
			return(1);
			}

		MaxSeedCoreDepth = maxcoredepth->count ? maxcoredepth->ival[0] : cDfltMaxSeedCoreDepth;
		if(MaxSeedCoreDepth < 1000 || MaxSeedCoreDepth > 20000)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: maximum depth to explore seed cores '-D%d' must be in range 1000..25000",MaxSeedCoreDepth);
			return(1);
			}

		SWMatchScore = matchscore->count ? matchscore->ival[0] : cDfltSWMatchScore;
		if(SWMatchScore < 1 || SWMatchScore > cMaxAllowedSWScore)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: SW Match score '-x%d' must be in range 1..%d",SWMatchScore,cMaxAllowedSWScore);
			return(1);
			}
		SWMismatchPenalty = mismatchpenalty->count ? mismatchpenalty->ival[0] : abs(cDfltSWMismatchPenalty);
		if(SWMismatchPenalty < 1 || SWMismatchPenalty > cMaxAllowedSWScore)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: SW Mismatch penalty '-X%d' must be in range 1..%d",SWMismatchPenalty,cMaxAllowedSWScore);
			return(1);
			}
		SWGapOpenPenalty = gapopenpenalty->count ? gapopenpenalty->ival[0] : abs(cDfltSWGapOpenPenalty);
		if(SWGapOpenPenalty < 1 || SWGapOpenPenalty > cMaxAllowedSWScore)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: SW Gap open penalty '-y%d' must be in range 1..%d",SWGapOpenPenalty,cMaxAllowedSWScore);
			return(1);
			}
		SWGapExtnPenalty = gapextnpenalty->count ? gapextnpenalty->ival[0] : abs(cDfltSWGapExtnPenalty);
		if(SWGapExtnPenalty < 1 || SWGapExtnPenalty > cMaxAllowedSWScore)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: SW Gap extension penalty '-Y%d' must be in range 1..%d",SWGapExtnPenalty,cMaxAllowedSWScore);
			return(1);
			}
		SWProgExtnPenaltyLen = progextnpenaltylen->count ? progextnpenaltylen->ival[0] : cDfltSWProgExtnLen;
		if(SWProgExtnPenaltyLen < 1 || SWProgExtnPenaltyLen > 63)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: SW Apply gap extension progress penalty for gaps at least '-z%d' must be in range 1..%d",SWProgExtnPenaltyLen,63);
			return(1);
			}

		MinPBSeqLen = minpbseqlen->count ? minpbseqlen->ival[0] : cDfltMinPBSeqLen;
		if(MinPBSeqLen < cMinPBSeqLen || MinPBSeqLen > cMaxMinPBSeqLen)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Minimum accepted PacBio length '-l%d' must be in range %d..%dbp",MinPBSeqLen,cMinPBSeqLen,cMaxMinPBSeqLen);
			return(1);
			}

		MinPBSeqOverlap = minpbseqovl->count ? minpbseqovl->ival[0] : cDfltMinErrCorrectLen;
		if(MinPBSeqOverlap < cMinPBSeqLen || MinPBSeqOverlap > cMaxMinPBSeqLen)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Minimum PacBio overlap required length '-L%d' must be in range %d..%dbp",MinPBSeqOverlap,cMinPBSeqLen,cMaxMinPBSeqLen);
			return(1);
			}


		MaxArtefactDev = maxartefactdev->count ? maxartefactdev->ival[0] : cDfltMaxArtefactDev;
		if(MaxArtefactDev <= 0 || MaxArtefactDev > cMaxMaxArtefactDev)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Max overlap artefactual deviation '-a%d' must be either 0 or in range %d..%d",MaxArtefactDev,cMinMaxArtefactDev,cMaxMaxArtefactDev);
			return(1);
			}
		if(MaxArtefactDev < 0)
			MaxArtefactDev = 0;

		MinHCSeqLen = minhcseqlen->count ? minhcseqlen->ival[0] : 1000;
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

		for(NumPacBioFiles=Idx=0;NumPacBioFiles < cMaxInFileSpecs && Idx < pacbiofiles->count; Idx++)
			{
			pszPacBioFiles[Idx] = NULL;
			if(pszPacBioFiles[NumPacBioFiles] == NULL)
				pszPacBioFiles[NumPacBioFiles] = new char [_MAX_PATH];
			strncpy(pszPacBioFiles[NumPacBioFiles],pacbiofiles->filename[Idx],_MAX_PATH);
			pszPacBioFiles[NumPacBioFiles][_MAX_PATH-1] = '\0';
			CUtility::TrimQuotedWhitespcExtd(pszPacBioFiles[NumPacBioFiles]);
			if(pszPacBioFiles[NumPacBioFiles][0] != '\0')
				NumPacBioFiles++;
			}

		if(!NumPacBioFiles)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: After removal of whitespace, no input PacBio file(s) specified with '-i<filespec>' option)\n");
			exit(1);
			}

		if(hiconffiles->count)
			{
			for(NumHiConfFiles=Idx=0;NumHiConfFiles < cMaxInFileSpecs && Idx < hiconffiles->count; Idx++)
				{
				pszHiConfFiles[Idx] = NULL;
				if(pszHiConfFiles[NumHiConfFiles] == NULL)
					pszHiConfFiles[NumHiConfFiles] = new char [_MAX_PATH];
				strncpy(pszHiConfFiles[NumHiConfFiles],hiconffiles->filename[Idx],_MAX_PATH);
				pszHiConfFiles[NumHiConfFiles][_MAX_PATH-1] = '\0';
				CUtility::TrimQuotedWhitespcExtd(pszHiConfFiles[NumHiConfFiles]);
				if(pszHiConfFiles[NumHiConfFiles][0] != '\0')
					NumHiConfFiles++;
				}

			if(!NumHiConfFiles)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: After removal of whitespace, no input high confidence reads file(s) specified with '-I<filespec>' option)\n");
				exit(1);
				}
			}
		else
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
		MinConcScore = minconcscore->count ? minconcscore->ival[0] : 3;
		if(MinConcScore < 0 || MinConcScore > 9)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: error corrected sequences trimming threshold score '-s%d' must be in range 0..9",MinConcScore);
			return(1);
			}

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
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Max overlap artefactual deviation '-a%d' must be either 0 or in range %d..%d",MaxArtefactDev,cMinMaxArtefactDev,cMaxMaxArtefactDev);
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

		MinNumSeedCores = minseedcores->count ? minseedcores->ival[0] : 10;
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
		if(MaxSeedCoreDepth < 1000 || MaxSeedCoreDepth > 25000)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: maximum depth to explore seed cores '-D%d' must be in range 1000..25000",MaxSeedCoreDepth);
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

		MinPBSeqOverlap = minpbseqovl->count ? minpbseqovl->ival[0] : min(cDfltMinErrCorrectLen,MinPBSeqLen);
		if(MinPBSeqOverlap < cMinPBSeqLen || MinPBSeqOverlap > cMaxMinPBSeqLen)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Minimum PacBio overlap required length '-L%d' must be in range %d..%dbp",MinPBSeqOverlap,cMinPBSeqLen,cMaxMinPBSeqLen);
			return(1);
			}

		for(NumPacBioFiles=Idx=0;NumPacBioFiles < cMaxInFileSpecs && Idx < pacbiofiles->count; Idx++)
			{
			pszPacBioFiles[Idx] = NULL;
			if(pszPacBioFiles[NumPacBioFiles] == NULL)
				pszPacBioFiles[NumPacBioFiles] = new char [_MAX_PATH];
			strncpy(pszPacBioFiles[NumPacBioFiles],pacbiofiles->filename[Idx],_MAX_PATH);
			pszPacBioFiles[NumPacBioFiles][_MAX_PATH-1] = '\0';
			CUtility::TrimQuotedWhitespcExtd(pszPacBioFiles[NumPacBioFiles]);
			if(pszPacBioFiles[NumPacBioFiles][0] != '\0')
				NumPacBioFiles++;
			}

		if(!NumPacBioFiles)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: After removal of whitespace, no input error corrected file(s) specified with '-i<filespec>' option)\n");
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
		}

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Processing mode: '%s'",pszMode);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Sampling input sequence rate: %d",SampleRate);

	if(PMode != ePBPMConsensus)
		{
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Use seed cores of this length when identifying putative overlapping sequences: %dbp",MinSeedCoreLen);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Require at least this many seed cores between overlapping sequences: %d",MinNumSeedCores);

		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Offset cores by this many bp: %d",DeltaCoreOfs);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Maximum seed core depth: %d",MaxSeedCoreDepth);

		gDiagnostics.DiagOutMsgOnly(eDLInfo,"SW score for matching bases: %d",SWMatchScore);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"SW mismatch penalty: %d",SWMismatchPenalty);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"SW gap opening penalty: %d",SWGapOpenPenalty);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"SW gap extension penalty: %d",SWGapExtnPenalty);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"SW gap extension penalty only applied for gaps of at least this size: %d",SWProgExtnPenaltyLen);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"classify overlaps as artefactual if sliding window of 1Kbp over any overlap deviates by more than this percentage: %d",MaxArtefactDev);
		if(PMode == ePBPMErrCorrect)
			{
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Minimum PacBio sequence length: %dbp",MinPBSeqLen);
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Minimum PacBio overlap required for error correction contribution: %d",MinPBSeqOverlap);
			}
		else
			{
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Minimum error corrected sequence length: %dbp",MinPBSeqLen);
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Minimum overlap required: %d",MinPBSeqOverlap);
			}


		if(NumHiConfFiles)
			{
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Minimum high confidence sequence length: %dbp",MinHCSeqLen);
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Minimum high confidence sequence onto PacBio overlap required for error correction contribution: %d",MinHCSeqOverlap);
			}

		for(Idx=0; Idx < NumPacBioFiles; Idx++)
			{
			if(PMode == ePBPMErrCorrect)
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
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Trimming error corrected PacBio sequences until mean 100bp score at least: %d",MinConcScore);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Error corrected and trimmed PacBio sequences must be at least this long: %d",MinErrCorrectLen);
		}

	if(PMode == ePBPMConsensus)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Input multialignment file spec: '%s'",pszPacBioFiles[ 0]);

	if(PMode != ePBPMOverlapDetail)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Output error corrected sequences file: '%s'",szOutFile);
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
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(SampleRate),"samplerate",&SampleRate);
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

		
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(MinPBSeqLen),"minpbseqlen",&MinPBSeqLen);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(MinPBSeqOverlap),"minpbseqovl",&MinPBSeqOverlap);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(MaxArtefactDev),"artefactdev",&MaxArtefactDev);

		if(NumHiConfFiles)
			{
			ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(MinHCSeqLen),"minhcseqlen",&MinHCSeqLen);
			ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(MinHCSeqOverlap),"minhcseqovl",&MinHCSeqOverlap);
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
	Rslt = ProcPacBioErrCorrect((etPBPMode)PMode,SampleRate,DeltaCoreOfs,MaxSeedCoreDepth,MinSeedCoreLen,MinNumSeedCores,SWMatchScore,-1 * SWMismatchPenalty,-1 * SWGapOpenPenalty,-1 * SWGapExtnPenalty,SWProgExtnPenaltyLen,
								MinPBSeqLen,MinPBSeqOverlap,MaxArtefactDev,MinHCSeqLen,MinHCSeqOverlap,MinErrCorrectLen,MinConcScore,
								NumPacBioFiles,pszPacBioFiles,NumHiConfFiles,pszHiConfFiles,szOutFile,szOutMAFile,NumThreads);
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
		int SampleRate,				// sample input sequences at this rate (1..100)
		int DeltaCoreOfs,			// offset by this many bp the core windows of coreSeqLen along the probe sequence when checking for overlaps
		int MaxSeedCoreDepth,		// only further process a seed core if there are no more than this number of matching cores in all targeted sequences
		int MinSeedCoreLen,			// use seed cores of this length when identifying putative overlapping scaffold sequences
		int MinNumSeedCores,        // require at least this many seed cores between overlapping scaffold sequences
		int SWMatchScore,			// score for matching bases (0..50)
		int SWMismatchPenalty,		// mismatch penalty (-50..0)
		int SWGapOpenPenalty,		// gap opening penalty (-50..0)
		int SWGapExtnPenalty,		// gap extension penalty (-50..0)
		int SWProgExtnPenaltyLen,	// progressive gap scoring then only apply gap extension score if gap at least this length (0..63) - use if aligning PacBio
		int MinPBSeqLen,			// only accepting PacBio reads of at least this length (defaults to 5Kbp)
		int MinPBSeqOverlap,		// any overlap of a PacBio onto a target PacBio must be of at least this many bp to be considered for contributing towards error correction (defaults to 5Kbp) 
		int MaxArtefactDev,			// classify overlaps as artefactual if sliding window of 1Kbp over any overlap deviates by more than this percentage from the overlap mean
		int MinHCSeqLen,			// only accepting hiconfidence reads of at least this length (defaults to 1Kbp)
		int MinHCSeqOverlap,		// any overlap of a hiconfidence read onto a target PacBio read must be of at least this many bp to be considered for contributing towards error correction (defaults to 1Kbp) 
		int MinErrCorrectLen,		// error corrected and trimmed sequences must be at least this minimum length
		int MinConcScore,			// error corrected sequences trimmed until mean 100bp concensus score is at least this threshold
		int NumPacBioFiles,			// number of input pacbio file specs
		char *pszPacBioFiles[],		// input pacbio files
		int NumHiConfFiles,			// number of input hiconfidence file specs
		char *pszHiConfFiles[],		// input hiconfidence files		
	    char *pszOutFile,			// where to write error corrected sequences
		char *pszOutMAFile,			// where to write multiple alignments
		int NumThreads)				// maximum number of worker threads to use
{
int Rslt;
CPBErrCorrect *pPBErrCorrect;

if((pPBErrCorrect = new CPBErrCorrect)==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: Unable to instantiate PBErrCorrect");
	return(eBSFerrObj);
	}

Rslt = pPBErrCorrect->Process(PMode,SampleRate,DeltaCoreOfs,MaxSeedCoreDepth,MinSeedCoreLen,MinNumSeedCores,SWMatchScore,SWMismatchPenalty,SWGapOpenPenalty,SWGapExtnPenalty,SWProgExtnPenaltyLen,
								MinPBSeqLen,MinPBSeqOverlap,MaxArtefactDev,MinHCSeqLen,MinHCSeqOverlap,MinErrCorrectLen,MinConcScore,
								NumPacBioFiles,pszPacBioFiles,NumHiConfFiles,pszHiConfFiles,pszOutFile,pszOutMAFile,NumThreads);
delete pPBErrCorrect;
return(Rslt);
}

CPBErrCorrect::CPBErrCorrect() // relies on base classes constructors
{
m_pSfxArray = NULL;
m_pPBScaffNodes = NULL;
m_pMapEntryID2NodeIDs = NULL;
m_bMutexesCreated = false;
m_hErrCorFile = -1;
m_hMultiAlignFile = -1;
Init();
}

CPBErrCorrect::~CPBErrCorrect() // relies on base classes destructors
{
Reset(false);
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

m_NumPBScaffNodes = 0;
m_AllocdPBScaffNodes = 0;
m_MaxPBSeqLen = 0;

m_NumOverlapProcessed = 0;
m_ProvOverlapping = 0;
m_ProvOverlapped = 0;
m_ProvContained = 0;
m_ProvArtefact = 0;
m_ProvSWchecked = 0;

m_MinErrCorrectLen = cDfltMinErrCorrectLen;
m_MinConcScore = 2;

m_PMode = ePBPMErrCorrect;

m_szScaffLineBuff[0] = '\0';
m_ScaffLineBuffIdx = 0;

m_OverlapFloat = cDfltMaxOverlapFloat;
m_MinPBSeqLen = cDfltMinPBSeqLen;	
m_MinPBSeqOverlap = cDfltMinErrCorrectLen;
m_MaxArtefactDev = cDfltMaxArtefactDev;
m_MinHCSeqLen = cDfltMinHCSeqLen;
m_MinHCSeqOverlap = cDfltMinHCSeqOverlap;

m_DeltaCoreOfs = cDfltDeltaCoreOfs;
m_MaxSeedCoreDepth = cDfltMaxSeedCoreDepth;
m_MinSeedCoreLen = cDfltSeedCoreLen;
m_MinNumSeedCores = cDfltNumSeedCores;

m_SWMatchScore = cDfltSWMatchScore;
m_SWMismatchPenalty = cDfltSWMismatchPenalty;	
m_SWGapOpenPenalty = cDfltSWGapOpenPenalty;
m_SWGapExtnPenalty = cDfltSWGapExtnPenalty;
m_SWProgExtnPenaltyLen = cDfltSWProgExtnLen;	

m_BinClusterSize = cDfltBinClusterSize;
m_MinPropBinned = cDfltMinPropBinned;

memset(m_ExactKmerDists,0,sizeof(m_ExactKmerDists));
m_TotAlignSeqLen = 0;
m_TotAlignSeqs = 0;
 
m_NumPacBioFiles = 0;
m_NumHiConfFiles = 0;

memset(m_szPacBioFiles,0,sizeof(m_szPacBioFiles));
memset(m_szHiConfFiles,0,sizeof(m_szHiConfFiles));

m_NumThreads = 0;
if(m_bMutexesCreated)
	DeleteMutexes();
m_bMutexesCreated = false; 
}

void
CPBErrCorrect::Reset(bool bSync)			// if bSync true then fsync before closing output file handles
{

Init();
}


int
CPBErrCorrect::CreateMutexes(void)
{
if(m_bMutexesCreated)
	return(eBSFSuccess);

#ifdef _WIN32
InitializeSRWLock(&m_hRwLock);
#else
if(pthread_rwlock_init (&m_hRwLock,NULL)!=0)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: unable to create rwlock");
	return(eBSFerrInternal);
	}
#endif

#ifdef _WIN32
if((m_hMtxIterReads = CreateMutex(NULL,false,NULL))==NULL)
	{
#else
if(pthread_mutex_init (&m_hMtxIterReads,NULL)!=0)
	{
	pthread_rwlock_destroy(&m_hRwLock);
#endif
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: unable to create mutex");
	return(eBSFerrInternal);
	}

#ifdef _WIN32
if((m_hMtxMHReads = CreateMutex(NULL,false,NULL))==NULL)
	{
#else
if(pthread_mutex_init (&m_hMtxMHReads,NULL)!=0)
	{
#endif
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: unable to create mutex");
#ifdef _WIN32
	CloseHandle(m_hMtxIterReads);
#else
	pthread_rwlock_destroy(&m_hRwLock);
	pthread_mutex_destroy(&m_hMtxIterReads);
#endif
	return(eBSFerrInternal);
	}

m_bMutexesCreated = true;
return(eBSFSuccess);
}

void
CPBErrCorrect::DeleteMutexes(void)
{
if(!m_bMutexesCreated)
	return;
#ifdef _WIN32
CloseHandle(m_hMtxIterReads);
CloseHandle(m_hMtxMHReads);
#else
pthread_mutex_destroy(&m_hMtxIterReads);
pthread_mutex_destroy(&m_hMtxMHReads);
pthread_rwlock_destroy(&m_hRwLock);
#endif
m_bMutexesCreated = false;
}

void
CPBErrCorrect::AcquireSerialise(void)
{
#ifdef _WIN32
WaitForSingleObject(m_hMtxIterReads,INFINITE);
#else
pthread_mutex_lock(&m_hMtxIterReads);
#endif
}

void
CPBErrCorrect::ReleaseSerialise(void)
{
#ifdef _WIN32
ReleaseMutex(m_hMtxIterReads);
#else
pthread_mutex_unlock(&m_hMtxIterReads);
#endif
}

void
CPBErrCorrect::AcquireSerialiseMH(void)
{
#ifdef _WIN32
WaitForSingleObject(m_hMtxMHReads,INFINITE);
#else
pthread_mutex_lock(&m_hMtxMHReads);
#endif
}

void
CPBErrCorrect::ReleaseSerialiseMH(void)
{
#ifdef _WIN32
ReleaseMutex(m_hMtxMHReads);
#else
pthread_mutex_unlock(&m_hMtxMHReads);
#endif
}

void
CPBErrCorrect::AcquireLock(bool bExclusive)
{
#ifdef _WIN32
if(bExclusive)
	AcquireSRWLockExclusive(&m_hRwLock);
else
	AcquireSRWLockShared(&m_hRwLock);
#else
if(bExclusive)
	pthread_rwlock_wrlock(&m_hRwLock);
else
	pthread_rwlock_rdlock(&m_hRwLock);
#endif
}

void
CPBErrCorrect::ReleaseLock(bool bExclusive)
{
#ifdef _WIN32
if(bExclusive)
	ReleaseSRWLockExclusive(&m_hRwLock);
else
	ReleaseSRWLockShared(&m_hRwLock);
#else
pthread_rwlock_unlock(&m_hRwLock);
#endif
}


// ProcessBioseqFile
// Process input biosequence file into suffix file
int
CPBErrCorrect::ProcessBioseqFile(int MinSeqLen,		// only accept for indexing sequences of at least this length
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
	if((NumSeqsAccepted + NumSeqsUnderlength) > ((SeqID * (UINT64)m_SampleRate) / 100))
		continue;
	NumSampled += 1;
	if(SeqLen < (UINT32)MinSeqLen)		// only accept for indexing sequences of at least this length)
		{
		NumSeqsUnderlength += 1;
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
gDiagnostics.DiagOut(eDLInfo,gszProcName,"ProcessFastaFile - %d parsed, %d sampled, %d accepted, %dbp mean length, %d sequences not accepted for indexing as length under %dbp ",
					SeqID,NumSampled,NumSeqsAccepted,NumSeqsAccepted == 0 ? 0 :(int)(TotAcceptedLen/NumSeqsAccepted),NumSeqsUnderlength,MinSeqLen);
return(Rslt);
}

// ProcessFastaFile
// Parse input fasta format file into a biosequence suffix array file
int
CPBErrCorrect::ProcessFastaFile(int MinSeqLen,		// only accept for indexing sequences of at least this length
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
NumSeqsAccepted = 0;
TotAcceptedLen = 0;
while((Rslt = SeqLen = Fasta.ReadSequence(&pSeqBuff[BuffOfs],(int)min(AvailBuffSize,(size_t)cMaxAllocBuffChunk),true,false)) > eBSFSuccess)
	{
	if(SeqLen == eBSFFastaDescr)		// just read a descriptor line
		{
		SeqID++;
		if(bEntryCreated)				// add any previous entry
			{
			if((NumSeqsAccepted + NumSeqsUnderlength) <= ((SeqID * (UINT64)m_SampleRate) / 100))
				{
				NumSampled += 1;
				if(BuffOfs < (size_t)MinSeqLen)
					NumSeqsUnderlength += 1;
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

if(Rslt >= eBSFSuccess && bEntryCreated && BuffOfs > 0)			// close entry
	{
	if((NumSeqsAccepted + NumSeqsUnderlength) <= ((SeqID * (UINT64)m_SampleRate) / 100))
		{
		NumSampled += 1;
		if(BuffOfs < (size_t)MinSeqLen)
			{
			NumSeqsUnderlength += 1;
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
gDiagnostics.DiagOut(eDLInfo,gszProcName,"ProcessFastaFile - %d parsed, %d sampled, %d accepted, %dbp mean length, %d sequences not accepted for indexing as length under %dbp ",
					SeqID,NumSampled,NumSeqsAccepted,NumSeqsAccepted == 0 ? 0 :(int)(TotAcceptedLen/NumSeqsAccepted),NumSeqsUnderlength,MinSeqLen);
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
		Reset(false);
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
	Reset(false);
	return(eBSFerrOpnFile);
	}

return(SumFileSizes);
}


int
CPBErrCorrect::GenConsensusFromMAF(int MinErrCorrectLen,		// error corrected sequences must be at least this minimum length
			 int MinConcScore,			// error corrected sequences trimmed until mean 100bp concensus score is at least this threshold
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
CPBErrCorrect::LoadSeqs(int MinSeqLen,int NumTargFiles,char **pszTargFiles,		// parse, and index sequences in these files into in memory suffix array; file expected to contain either fasta or fastq sequences
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
		Reset(false);
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
	Reset(false);
	return(eBSFerrOpnFile);
	}

Rslt = eBSFSuccess;

for (int n = 0; Rslt >= eBSFSuccess &&  n < glob.FileCount(); ++n)
	{
		// try opening as a fasta file, if that fails then try as a bioseq
	Rslt = ProcessFastaFile(MinSeqLen,glob.File(n),Flags);
	if(Rslt == eBSFerrNotFasta)
		Rslt = ProcessBioseqFile(MinSeqLen,glob.File(n),Flags);
	if(Rslt < eBSFSuccess)
		{
		m_pSfxArray->Close(false);
		Reset(false);
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
		Reset(false);
		return(eBSFerrFastqSeqID);
		}
	}
return(eBSFSuccess);

}

int
CPBErrCorrect::Process(etPBPMode PMode,		// processing mode
		int SampleRate,				// sample input sequences at this rate (1..100)
		int DeltaCoreOfs,			// offset by this many bp the core windows of coreSeqLen along the probe sequence when checking for overlaps
		int MaxSeedCoreDepth,		// only further process a seed core if there are no more than this number of matching cores in all targeted sequences
		int MinSeedCoreLen,			// use seed cores of this length when identifying putative overlapping scaffold sequences
		int MinNumSeedCores,        // require at least this many seed cores between overlapping scaffold sequences
		int SWMatchScore,			// score for matching bases (0..50)
		int SWMismatchPenalty,		// mismatch penalty (-50..0)
		int SWGapOpenPenalty,		// gap opening penalty (-50..0)
		int SWGapExtnPenalty,		// gap extension penalty (-50..0)
		int SWProgExtnPenaltyLen,	// progressive gap scoring then only apply gap extension score if gap at least this length (0..63) - use if aligning PacBio
		int MinPBSeqLen,			// only accepting PacBio reads of at least this length (defaults to 5Kbp)
		int MinPBSeqOverlap,		// any overlap of a PacBio onto a target PacBio must be of at least this many bp to be considered for contributing towards error correction (defaults to 5Kbp) 
		int MaxArtefactDev,			// classify overlaps as artefactual if sliding window of 1Kbp over any overlap deviates by more than this percentage from the overlap mean
		int MinHCSeqLen,			// only accepting hiconfidence reads of at least this length (defaults to 1Kbp)
		int MinHCSeqOverlap,		// any overlap of a hiconfidence read onto a target PacBio read must be of at least this many bp to be considered for contributing towards error correction (defaults to 1Kbp) 
		int MinErrCorrectLen,		// error corrected and trimmed sequences must be at least this minimum length
		int MinConcScore,			// error corrected sequences trimmed until mean 100bp concensus score is at least this threshold
		int NumPacBioFiles,			// number of input pacbio file specs
		char *pszPacBioFiles[],		// input pacbio files
		int NumHiConfFiles,			// number of input hiconfidence file specs
		char *pszHiConfFiles[],		// input hiconfidence files		
	    char *pszErrCorFile,		// name of file into which write error corrected sequences
		char *pszMultiAlignFile,	// name of file into which write multiple alignments
		int NumThreads)				// maximum number of worker threads to use
{
int Rslt = eBSFSuccess;
int Idx;
UINT32 NumTargSeqs;
UINT32 CurNodeID;
UINT32 MaxSeqLen;
tsPBScaffNode *pCurPBScaffNode;

Reset(false);

CreateMutexes();

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
m_MinPBSeqOverlap = MinPBSeqOverlap; 
m_MaxArtefactDev = MaxArtefactDev;
m_MinHCSeqLen = MinHCSeqLen;
m_MinHCSeqOverlap = MinHCSeqOverlap;

m_MinErrCorrectLen = MinErrCorrectLen;
m_MinConcScore = MinConcScore;
m_SampleRate = SampleRate;


switch(PMode) {
	case ePBPMOverlapDetail:
		m_MinPropBinned = cDfltScaffMinPropBinned;
		m_OverlapFloat = cDfltScaffMaxOverlapFloat;
		break;

	case ePBPMConsensus:
		Rslt = GenConsensusFromMAF(MinErrCorrectLen,MinConcScore,pszErrCorFile,pszPacBioFiles[0]);
		Reset(false);
		return(Rslt);

	default:
		m_MinPropBinned = cDfltMinPropBinned;
		m_OverlapFloat = cDfltMaxOverlapFloat;
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

if(pszMultiAlignFile != NULL && pszMultiAlignFile[0] != '\0')
	{
	strncpy(m_szMultiAlignFile,pszMultiAlignFile,sizeof(m_szMultiAlignFile));
	m_szMultiAlignFile[sizeof(m_szMultiAlignFile)-1] = '\0';
	}	
else
	m_szMultiAlignFile[0] = '\0';

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

m_NumThreads = NumThreads;	
if(m_pSfxArray != NULL)
	delete m_pSfxArray;
if((m_pSfxArray = new CSfxArrayV3) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadTargetSeqs: Unable to instantiate instance of CSfxArrayV3");
	return(eBSFerrObj);
	}
m_pSfxArray->Reset(false);
m_pSfxArray->SetMaxQSortThreads(m_NumThreads);

Rslt=m_pSfxArray->Open(false,false);
if(Rslt !=eBSFSuccess)
	{
	while(m_pSfxArray->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,m_pSfxArray->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"CreateBioseqSuffixFile, unable to create in-memory suffix array - error %d %s",Rslt,m_pSfxArray->GetErrMsg());
	Reset(false);
	return(Rslt);
	}

if((Rslt=m_pSfxArray->SetDescription((char *)"inmem")) != eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"CreateBioseqSuffixFile: Unable to set description 'inmem'");
	Reset(false);
	return(Rslt);
	}
if((Rslt=m_pSfxArray->SetTitle((char *)"inmem")) != eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"CreateBioseqSuffixFile: Unable to set title 'inmem'");
	Reset(false);
	return(Rslt);
	}

if((Rslt = m_pSfxArray->SetDatasetName((char *)"inmem")) != eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"CreateBioseqSuffixFile: Unable to set dataset name 'inmem'");
	Reset(false);
	return(Rslt);
	}
m_pSfxArray->SetInitalSfxAllocEls(SumFileSizes);	// just a hint which is used for initial allocations by suffix processing


if((Rslt = LoadSeqs(m_MinPBSeqLen,NumPacBioFiles,pszPacBioFiles,cFlgLCSeq)) < eBSFSuccess)
	{
	Reset(false);
	return(Rslt);
	}
	// get number of sequences BioSeqs accepted for error correction
if((NumTargSeqs = m_pSfxArray->GetNumEntries()) < 1)
	NumTargSeqs = 0;
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Loaded and accepted for processing a total of %d PacBio sequences",NumTargSeqs);

if(NumTargSeqs < 3)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Need at least 3 PacBio sequences for errror correction accepted from file(s)");
	m_pSfxArray->Close(false);
	Reset(false);
	return(eBSFerrNoEntries);
	}

int NumHCSeqs = 0;
if(NumHiConfFiles && pszHiConfFiles != NULL)			// load any optional high confidence sequences requested by user
	{
	if((Rslt = LoadSeqs(m_MinHCSeqLen,NumHiConfFiles,pszHiConfFiles,cFlgHCSeq)) < eBSFSuccess)
		{
		m_pSfxArray->Close(false);
		Reset(false);
		return(Rslt);
		}
		// get number of high confidence sequences accepted
	NumHCSeqs = m_pSfxArray->GetNumEntries() - NumTargSeqs;
	if(NumHCSeqs < 1)
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"No high confidence sequences for errror correction accepted from file(s)");
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
		Reset(false);
		return(Rslt);
		}
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Initialised for over occurring K-mers");
	}

if((m_pPBScaffNodes = new tsPBScaffNode [NumTargSeqs + 1]) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate for %d nodes",NumTargSeqs);
	Reset(false);
	return(eBSFerrMem);
	}
memset(m_pPBScaffNodes,0,sizeof(tsPBScaffNode) * (NumTargSeqs+1));
m_AllocdPBScaffNodes = NumTargSeqs;
if((m_pMapEntryID2NodeIDs = new UINT32 [NumTargSeqs + 1]) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate for %d mapping entries nodes",NumTargSeqs);
	Reset(false);
	return(eBSFerrMem);
	}
memset(m_pMapEntryID2NodeIDs,0,sizeof(UINT32) * (NumTargSeqs+1));
m_NumPBScaffNodes = 0;

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Initialising for %d sequences",NumTargSeqs);	
MaxSeqLen = 0;
pCurPBScaffNode = m_pPBScaffNodes;
for(CurNodeID = 1; CurNodeID <= NumTargSeqs; CurNodeID++,pCurPBScaffNode++)
	{
	UINT16 SeqFlags;
	pCurPBScaffNode->SeqLen = m_pSfxArray->GetSeqLen(CurNodeID);
	pCurPBScaffNode->EntryID = CurNodeID;
	SeqFlags = m_pSfxArray->GetIdentFlags(CurNodeID);
	pCurPBScaffNode->flgHCseq = SeqFlags & cFlgHCSeq ? 1 : 0;
	pCurPBScaffNode->flgUnderlength = pCurPBScaffNode->SeqLen < m_MinPBSeqLen ? 1 : 0;
	if(MaxSeqLen == 0 || pCurPBScaffNode->SeqLen > (UINT32)MaxSeqLen)
		MaxSeqLen = pCurPBScaffNode->SeqLen;
	}
m_NumPBScaffNodes = NumTargSeqs;
m_MaxPBSeqLen = MaxSeqLen;

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Sorting %d sequences ...",NumTargSeqs);

if(m_NumPBScaffNodes > 1)
	{
	// sort scaffold nodes by sequence length descending
	m_mtqsort.SetMaxThreads(NumThreads);
	m_mtqsort.qsort(m_pPBScaffNodes,m_NumPBScaffNodes,sizeof(tsPBScaffNode),SortLenDescending);
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Sorting %d sequences completed",NumTargSeqs);
pCurPBScaffNode = m_pPBScaffNodes;
for(CurNodeID = 1; CurNodeID <= NumTargSeqs; CurNodeID++,pCurPBScaffNode++)
	{
	pCurPBScaffNode->NodeID = CurNodeID;
	m_pMapEntryID2NodeIDs[pCurPBScaffNode->EntryID-1] = CurNodeID;
	}

if(m_szErrCorFile[0] != '\0')
	{
#ifdef _WIN32
	m_hErrCorFile = open(m_szErrCorFile,( O_WRONLY | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC),(_S_IREAD | _S_IWRITE));
#else
	if((m_hErrCorFile = open(m_szErrCorFile,O_WRONLY | O_CREAT,S_IREAD | S_IWRITE))!=-1)
		if(ftruncate(m_hErrCorFile,0)!=0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate %s - %s",m_szErrCorFile,strerror(errno));
			Reset(false);
			return(eBSFerrCreateFile);
			}
#endif
	if(m_hErrCorFile < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate output file '%s'",m_szErrCorFile);
		Reset(false);
		return(eBSFerrCreateFile);
		}
	}
else
	m_hErrCorFile = -1;

if(m_PMode == ePBPMOverlapDetail)
	{
	m_ScaffLineBuffIdx=sprintf(m_szScaffLineBuff,"\"Class\",\"ProbeID\",\"ProbDescr\",\"TargID\",\"TargDescr\",\"SeedHits\",\"ProbeSense\",\"TargSense\",\"ProbeLen\",\"TargLen\",\"ProbeAlignLength\",\"TargAlignLength\",\"PeakScore\",\"FinalScore\",\"NumAlignedBases\",\"NumExactBases\",\"NumProbeInserts\",\"NumProbeInsertBases\",\"NumTargInserts\",\"NumTargInsertBases\",\"ProbeStartOfs\",\"TargStartOfs\",\"ProbeEndOfs\",\"TargEndOfs\",\"ProbeOfs5\",\"TargOfs5\",\"ProbeOfs3\",\"TargOfs3\"");
#ifdef _PEAKSCOREACCEPT_
	m_ScaffLineBuffIdx+=sprintf(&m_szScaffLineBuff[m_ScaffLineBuffIdx],",\"PSAlignLength\",\"PSPeakScore\",\"PSNumAlignedBases\",\"PSNumExactBases\",\"PSNumProbeInserts\",\"PSNumProbeInsertBases\",\"PSNumTargInserts\",\"PSNumTargInsertBases\",\"PSProbeStartOfs\",\"PSTargStartOfs\",\"ProbeEndOfs\",\"TargEndOfs\",\"PSProbeOfs5\",\"PSTargOfs5\",\"PSProbeOfs3\",\"PSTargOfs3\"");
#endif
#ifdef _EXACTMATCHLENDIST_
	for(int ExactLenIdx = 1; ExactLenIdx <= cExactMatchLenDist; ExactLenIdx++)
		m_ScaffLineBuffIdx+=sprintf(&m_szScaffLineBuff[m_ScaffLineBuffIdx],",\"ExactLen:%d\"",ExactLenIdx);
#endif
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
			Reset(false);
			return(eBSFerrCreateFile);
			}
#endif
	if(m_hMultiAlignFile < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate output file '%s'",m_szMultiAlignFile);
		Reset(false);
		return(eBSFerrCreateFile);
		}
	}
else
	m_hMultiAlignFile = -1;

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Identifying sequence overlaps ...");
Rslt = IdentifySequenceOverlaps(MaxSeqLen,NumThreads);
Reset(false);
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
_endthreadex(0);
return(eBSFSuccess);
#else
pthread_exit(NULL);
#endif
}

int
CPBErrCorrect::IdentifySequenceOverlaps(int MaxSeqLen,		// max length sequence to be overlapped
									int NumOvlpThreads)		// identify all read overlaps using this many threads
{
tsThreadPBErrCorrect *pThreadPutOvlps;
int ThreadIdx;
tsThreadPBErrCorrect *pThreadPar;

memset(m_ExactKmerDists,0,sizeof(m_ExactKmerDists));
m_TotAlignSeqLen = 0;
m_TotAlignSeqs = 0;

pThreadPutOvlps = new tsThreadPBErrCorrect [NumOvlpThreads];

pThreadPar = pThreadPutOvlps;
for(ThreadIdx = 0; ThreadIdx < NumOvlpThreads; ThreadIdx++,pThreadPar++)
	{
	memset(pThreadPar,0,sizeof(tsThreadPBErrCorrect));
	pThreadPar->AllocdCoreHits = cAllocdNumCoreHits;
	pThreadPar->AllocdCoreHitsSize = cAllocdNumCoreHits * sizeof(tsCoreHit);
#ifdef _WIN32
	pThreadPar->pCoreHits = (tsCoreHit *)malloc(pThreadPar->AllocdCoreHitsSize);	
	if(pThreadPar->pCoreHits == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"IdentifySequenceOverlaps: Core hits memory allocation of %llu bytes - %s",pThreadPar->AllocdCoreHitsSize,strerror(errno));
		break;
		}
#else
	if((pThreadPar->pCoreHits = (tsCoreHit *)mmap(NULL,pThreadPar->AllocdCoreHitsSize, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0)) == MAP_FAILED)
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
	pThreadPar->MinPropBinned = m_MinPropBinned;
	pThreadPar->MaxAcceptHitsPerSeedCore = cDfltMaxAcceptHitsPerSeedCore;
	pThreadPar->MinPBSeqLen = m_MinPBSeqLen;
	pThreadPar->MinOverlapLen = m_MinPBSeqOverlap;
	}

if(ThreadIdx != NumOvlpThreads)	// any errors whilst allocating memory for core hits?
	{
	do {
		if(pThreadPar->pmtqsort != NULL)
			delete pThreadPar->pmtqsort;

		if(pThreadPar->pAntisenseKmers != NULL)
			delete pThreadPar->pAntisenseKmers;

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
		if(pThreadPar->pSW != NULL)
			delete pThreadPar->pSW;
		pThreadPar -= 1;
		ThreadIdx -= 1;
		}
	while(ThreadIdx >= 0);
	delete pThreadPutOvlps;
	Reset(false);
	return((INT64)eBSFerrMem);
	}

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

// allow threads a few seconds to startup
#ifdef _WIN32
Sleep(5000);
#else
sleep(5);
#endif
pThreadPar = pThreadPutOvlps;
for (ThreadIdx = 0; ThreadIdx < NumOvlpThreads; ThreadIdx++, pThreadPar++)
	{
#ifdef _WIN32
	while (WAIT_TIMEOUT == WaitForSingleObject(pThreadPar->threadHandle, 60000))
		{
		AcquireSerialise();
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Progress: %u processed, SW aligned: %u, Overlapping: %u, Overlapped: %u, Contained: %u, Artefact: %u",
							m_NumOverlapProcessed,m_ProvSWchecked,m_ProvOverlapping,m_ProvOverlapped,m_ProvContained,m_ProvArtefact);
		ReleaseSerialise();
		};
	CloseHandle(pThreadPar->threadHandle);
#else
	struct timespec ts;
	int JoinRlt;
	clock_gettime(CLOCK_REALTIME, &ts);
	ts.tv_sec += 60;
	while ((JoinRlt = pthread_timedjoin_np(pThreadPar->threadID, NULL, &ts)) != 0)
		{
		AcquireSerialise();
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Progress: %u processed, SW aligned: %u, Overlapping: %u, Overlapped: %u, Contained: %u, Artefact: %u",
							m_NumOverlapProcessed,m_ProvSWchecked,m_ProvOverlapping,m_ProvOverlapped,m_ProvContained,m_ProvArtefact);
		ReleaseSerialise();
		ts.tv_sec += 60;
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
	}

delete pThreadPutOvlps;
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Completed: %u processed, SW aligned: %u, Overlapping: %u, Overlapped: %u, Contained: %u, Artefact: %u",
							m_NumOverlapProcessed,m_ProvSWchecked,m_ProvOverlapping,m_ProvOverlapped,m_ProvContained,m_ProvArtefact);

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
UINT32  OverlapSLen;
UINT32  OverlapALen;
UINT32  PropSBinsOverlap;
UINT32  PropABinsOverlap;
UINT32 TargLen;
UINT32 CurTargCoreHitCnts;
tsPBScaffNode *pTargNode;
UINT32 ProbeAlignLength;
UINT32 TargAlignLength;
UINT32 TargSeqLen;
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
sPBCoreHitCnts *pLowestSummaryHitCnts;
UINT32 HitIdx;
int NumHitsFlgMulti;
tsCoreHit *pCoreHit;
UINT32 CurTargNodeID;
UINT32 CurProbeNodeID;
UINT32 CurProbeOfs;
UINT32 CurTargOfs;
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
tsPBScaffNode *pHitScaffNode;
sPBCoreHitCnts *pSummaryCnts;
UINT32 CurNodeID;
int NumInMultiAlignment;
int Class;
char szTargSeqName[cMaxDatasetSpeciesChrom];
char szProbeSeqName[cMaxDatasetSpeciesChrom];

tsPBScaffNode *pCurPBScaffNode;
pThreadPar->pSW = NULL;
NumInMultiAlignment = 0;
for(CurNodeID = 1; CurNodeID <= m_NumPBScaffNodes; CurNodeID++)
	{
	pThreadPar->NumCoreHits = 0;
	pCurPBScaffNode = &m_pPBScaffNodes[CurNodeID-1];
	AcquireLock(false);
	if(pCurPBScaffNode->flgCurProc == 1 ||    // another thread already processing this sequence? 
			pCurPBScaffNode->flgHCseq ||	  // not error correcting already assumed to be high confidence sequences 
		    pCurPBScaffNode->flgUnderlength == 1 || // must be of a user specified minimum length 
			pCurPBScaffNode->SeqLen < (UINT32)pThreadPar->MinPBSeqLen) 
       	{
		ReleaseLock(false);
		continue;
		}
	pCurPBScaffNode->flgCurProc = 1;
	ReleaseLock(false);
	pThreadPar->bRevCpl = false;
	IdentifyCoreHits(CurNodeID,pThreadPar);

	pThreadPar->bRevCpl = true;
	IdentifyCoreHits(CurNodeID,pThreadPar);

	ProvOverlapping = 0;
	ProvOverlapped = 0;
	ProvContained = 0;
	ProvArtefact = 0;
	ProvSWchecked = 0;


	pThreadPar->NumTargCoreHitCnts = 0;
	memset(pThreadPar->TargCoreHitCnts,0,sizeof(pThreadPar->TargCoreHitCnts));
	if(pThreadPar->NumCoreHits)
		{
		if(pThreadPar->NumCoreHits > 1)
			{
			// sort core hits by ProbeNodeID.ProbeOfs.TargNodeID.TargOfs ascending so multiple hits onto a target from single probe offset can be detected
			pThreadPar->pmtqsort->qsort(pThreadPar->pCoreHits,pThreadPar->NumCoreHits,sizeof(tsCoreHit),SortCoreHitsByProbeTargOfs);
			CurProbeNodeID = 0;
			CurProbeOfs = 0;
			CurTargNodeID = 0;
			CurTargOfs = 0;
			NumHitsFlgMulti = 0;
			pCoreHit = pThreadPar->pCoreHits;
			for(HitIdx = 0; HitIdx < pThreadPar->NumCoreHits; HitIdx++,pCoreHit++)
				{
				if(pCoreHit->ProbeNodeID == CurProbeNodeID && 
					CurTargNodeID == pCoreHit->TargNodeID && 
					(pCoreHit->ProbeOfs >= CurProbeOfs && (pCoreHit->ProbeOfs <= (CurProbeOfs + m_BinClusterSize)))) 
					{
					NumHitsFlgMulti += 1;
					pCoreHit->flgMulti = 1;
					}				
				else
					{
					CurProbeNodeID = pCoreHit->ProbeNodeID;
					CurProbeOfs = pCoreHit->ProbeOfs;
					CurTargNodeID = pCoreHit->TargNodeID;
					CurTargOfs = pCoreHit->TargOfs;
					pCoreHit->flgMulti = 0;
					}    
				}

			// resort core hits by TargNodeID.TargOfs.ProbeNodeID.ProbeOfs ascending
			pThreadPar->pmtqsort->qsort(pThreadPar->pCoreHits,pThreadPar->NumCoreHits,sizeof(tsCoreHit),SortCoreHitsByTargProbeOfs);
			CurTargNodeID = 0;
			CurTargOfs = 0;
			pCoreHit = pThreadPar->pCoreHits;
			for(HitIdx = 0; HitIdx < pThreadPar->NumCoreHits; HitIdx++,pCoreHit++)
				{
				if(pCoreHit->flgMulti)
					continue;
				if(pCoreHit->ProbeNodeID == CurProbeNodeID && 
					CurTargNodeID == pCoreHit->TargNodeID && 
					(pCoreHit->TargOfs >= CurTargOfs && (pCoreHit->TargOfs <= (CurTargOfs + m_BinClusterSize)))) 
					{
					NumHitsFlgMulti += 1;
					pCoreHit->flgMulti = 1;
					}				
				else
					{
					CurProbeNodeID = pCoreHit->ProbeNodeID;
					CurProbeOfs = pCoreHit->ProbeOfs;
					CurTargNodeID = pCoreHit->TargNodeID;
					CurTargOfs = pCoreHit->TargOfs;
					}    
				}
			}

		// iterate and count hits for each TargNodeID whilst recording the loci of the first and last hit so can determine if overlap is a likely artefact
		CurSEntryIDHits = 0;
		CurAEntryIDHits = 0;
		CurTargNodeID = 0;
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

		pCoreHit = pThreadPar->pCoreHits;
		for(HitIdx = 0; HitIdx < pThreadPar->NumCoreHits; HitIdx++,pCoreHit++)
			{
			pHitScaffNode = &m_pPBScaffNodes[pCoreHit->ProbeNodeID-1];
			if(CurTargNodeID == 0)    // 0 if 1st hit about to be processed for a new target sequence
				{
				CurTargNodeID = pCoreHit->TargNodeID;
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

			if(pCoreHit->TargNodeID == CurTargNodeID) // same target sequence so check for starting/ending offsets and accumulate hit counts 
				{
				CurTargHitOfs = pCoreHit->TargOfs;
				CurProbeHitOfs = pCoreHit->ProbeOfs;
				if(pCoreHit->flgRevCpl == 0)
					{
					if(CurSTargEndOfs == 0 || CurTargHitOfs < CurSTargStartOfs)
						CurSTargStartOfs = CurTargHitOfs;
					if(CurTargHitOfs > CurSTargEndOfs)
						CurSTargEndOfs = CurTargHitOfs;	
					if(CurSProbeEndOfs == 0 || CurProbeHitOfs < CurSProbeStartOfs)
						CurSProbeStartOfs = CurProbeHitOfs;
					if(CurProbeHitOfs > CurSProbeEndOfs)
						CurSProbeEndOfs = CurProbeHitOfs;
					}
				else
					{
					if(CurATargEndOfs == 0 || CurTargHitOfs < CurATargStartOfs)
						CurATargStartOfs = CurTargHitOfs;
					if(CurTargHitOfs > CurATargEndOfs)
						CurATargEndOfs = CurTargHitOfs;	
					if(CurAProbeEndOfs == 0 || CurProbeHitOfs < CurAProbeStartOfs)
						CurAProbeStartOfs = CurProbeHitOfs;
					if(CurProbeHitOfs > CurAProbeEndOfs)
						CurAProbeEndOfs = CurProbeHitOfs;
					}

				if(pCoreHit->flgMulti != 1)
					{
					if(pCoreHit->flgRevCpl == 0)
						CurSEntryIDHits += 1;
					else
						CurAEntryIDHits += 1;
					}
				}

			// if just processed last core hit for the current target ...
			if(HitIdx + 1 == pThreadPar->NumCoreHits || pCoreHit[1].TargNodeID != CurTargNodeID)
				{
				// checking here if the overlap is very likely to be artefact
				// requiring that at least MinPropBinned of the bins along the putative alignment length contain hits
				// and that the first and last hit are consistent with either a completely contained or overlapped
				if(CurSEntryIDHits > 1 && ((OverlapSLen = CurSProbeEndOfs - CurSProbeStartOfs) > m_BinClusterSize))
					PropSBinsOverlap = (m_BinClusterSize * CurSEntryIDHits * 100)/OverlapSLen;
				else
					PropSBinsOverlap = 0;

				if(CurAEntryIDHits > 1 && ((OverlapALen = CurAProbeEndOfs - CurAProbeStartOfs) > m_BinClusterSize))
					PropABinsOverlap = (m_BinClusterSize * CurAEntryIDHits * 100)/OverlapALen;
				else
					PropABinsOverlap = 0;
				TargLen = m_pPBScaffNodes[pCoreHit->TargNodeID-1].SeqLen;
				if(PropSBinsOverlap >= pThreadPar->MinPropBinned) 
					{
					if((CurSProbeStartOfs >= m_OverlapFloat &&  CurSTargStartOfs >= m_OverlapFloat) ||
						((TargLen - CurSTargEndOfs) >= m_OverlapFloat && (pCurPBScaffNode->SeqLen - CurSProbeEndOfs) >= m_OverlapFloat))
						PropSBinsOverlap = 0;
					}
				if(PropABinsOverlap >= pThreadPar->MinPropBinned)
					{
					if((CurAProbeStartOfs >= m_OverlapFloat && CurATargStartOfs >= m_OverlapFloat) ||
						((TargLen - CurATargEndOfs) >= m_OverlapFloat && (pCurPBScaffNode->SeqLen - CurAProbeEndOfs) >= m_OverlapFloat))
						PropABinsOverlap = 0;
					}

				if((PropSBinsOverlap >= pThreadPar->MinPropBinned && CurSEntryIDHits >= pThreadPar->MinNumCores) || (PropABinsOverlap >= pThreadPar->MinPropBinned && CurAEntryIDHits >= pThreadPar->MinNumCores))
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
							CurTargNodeID = 0; 
							continue;
							}
						pSummaryCnts = pLowestSummaryHitCnts;
						}
					else
						pSummaryCnts = &pThreadPar->TargCoreHitCnts[pThreadPar->NumTargCoreHitCnts++];
					pSummaryCnts->TargNodeID = CurTargNodeID;
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
					}
				CurTargNodeID = 0; 
				}
			}
		}
 
	// can't process, SW over all would be too resource intensive, all targets which meet the minimum number of core hits requested so choose the top cMaxProbeSWs as ranked by the number of core hits
	if(pThreadPar->NumTargCoreHitCnts > 1)
		{
		pThreadPar->pmtqsort->qsort(pThreadPar->TargCoreHitCnts,pThreadPar->NumTargCoreHitCnts,sizeof(sPBCoreHitCnts),SortCoreHitsDescending);
		if(pThreadPar->NumTargCoreHitCnts > cMaxProbeSWs)		// clamp to no more than this many SW alignments
			pThreadPar->NumTargCoreHitCnts = cMaxProbeSWs;
		}

	NumInMultiAlignment = 0;
	if(pThreadPar->NumTargCoreHitCnts > 0)
		{
		LongSAligns = 0;
		LongAAligns = 0;
		if(pThreadPar->pSW == NULL)
			{
			AcquireSerialise();
			pThreadPar->pSW = new CSSW;
			pThreadPar->pSW->SetScores(m_SWMatchScore,m_SWMismatchPenalty,m_SWGapOpenPenalty,m_SWGapExtnPenalty,m_SWProgExtnPenaltyLen,min(63,m_SWProgExtnPenaltyLen+3),cAnchorLen);
			pThreadPar->pSW->PreAllocMaxTargLen(100000,m_PMode == ePBPMConsensus ? 0 : m_MaxPBSeqLen);
			ReleaseSerialise();
			}

		if(m_PMode == ePBPMErrCorrect && pThreadPar->NumTargCoreHitCnts >= 2)
			pThreadPar->pSW->StartMultiAlignments(pCurPBScaffNode->SeqLen,pThreadPar->pProbeSeq,pThreadPar->NumTargCoreHitCnts);

		m_pSfxArray->GetIdentName(pCurPBScaffNode->EntryID,sizeof(szProbeSeqName),szProbeSeqName);
		pThreadPar->pSW->SetProbe(pCurPBScaffNode->SeqLen,pThreadPar->pProbeSeq);
		pSummaryCnts = &pThreadPar->TargCoreHitCnts[0];
		NumInMultiAlignment = 0;
		for(CurTargCoreHitCnts = 0; CurTargCoreHitCnts < pThreadPar->NumTargCoreHitCnts; CurTargCoreHitCnts++,pSummaryCnts++)
			{
			if(pSummaryCnts->NumSHits < pThreadPar->MinNumCores && pSummaryCnts->NumAHits < pThreadPar->MinNumCores)
				continue;
			bTargSense = pSummaryCnts->NumSHits >= pSummaryCnts->NumAHits ? true :  false;
			pTargNode = &m_pPBScaffNodes[pSummaryCnts->TargNodeID-1];
			if(pTargNode->flgHCseq == 1)
				MinOverlapLen = m_MinHCSeqOverlap;
			else
				MinOverlapLen = pThreadPar->MinOverlapLen;

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
			pThreadPar->pSW->SetTarg(TargSeqLen,pThreadPar->pTargSeq);
			
			// restrict the range over which the SW will be processed to that of the overlap +/- m_OverlapFloat
			int Rslt;
			if(bTargSense)
				{
				if(pSummaryCnts->SProbeStartOfs < m_OverlapFloat)
					pSummaryCnts->SProbeStartOfs = 0;
				else
					pSummaryCnts->SProbeStartOfs -= m_OverlapFloat;
				if(pSummaryCnts->SProbeEndOfs + m_OverlapFloat >= pCurPBScaffNode->SeqLen)
					pSummaryCnts->SProbeEndOfs = pCurPBScaffNode->SeqLen - 1;
				else
					pSummaryCnts->SProbeEndOfs += m_OverlapFloat;
				if(pSummaryCnts->STargStartOfs < m_OverlapFloat)
					pSummaryCnts->STargStartOfs = 0;
				else
					pSummaryCnts->STargStartOfs -= m_OverlapFloat;
				if(pSummaryCnts->STargEndOfs + m_OverlapFloat >= TargSeqLen)
					pSummaryCnts->STargEndOfs = TargSeqLen - 1;
				else
					pSummaryCnts->STargEndOfs += m_OverlapFloat;
				Rslt = pThreadPar->pSW->SetAlignRange(pSummaryCnts->SProbeStartOfs,pSummaryCnts->STargStartOfs,
											pSummaryCnts->SProbeEndOfs + 1 - pSummaryCnts->SProbeStartOfs,pSummaryCnts->STargEndOfs + 1 - pSummaryCnts->STargStartOfs);
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
				if(pSummaryCnts->AProbeEndOfs + m_OverlapFloat >= pCurPBScaffNode->SeqLen)
					pSummaryCnts->AProbeEndOfs = pCurPBScaffNode->SeqLen - 1;
				else
					pSummaryCnts->AProbeEndOfs += m_OverlapFloat;
				if(pSummaryCnts->ATargStartOfs < m_OverlapFloat)
					pSummaryCnts->ATargStartOfs = 0;
				else
					pSummaryCnts->ATargStartOfs -= m_OverlapFloat;
				if(pSummaryCnts->ATargEndOfs + m_OverlapFloat >= TargSeqLen)
					pSummaryCnts->ATargEndOfs = TargSeqLen - 1;
				else
					pSummaryCnts->ATargEndOfs += m_OverlapFloat;

				Rslt = pThreadPar->pSW->SetAlignRange(pSummaryCnts->AProbeStartOfs,pSummaryCnts->ATargStartOfs,
											pSummaryCnts->AProbeEndOfs + 1 - pSummaryCnts->AProbeStartOfs,pSummaryCnts->ATargEndOfs + 1 - pSummaryCnts->ATargStartOfs);
				}
#ifdef _PEAKSCOREACCEPT_
			pPeakMatchesCell = pThreadPar->pSW->Align(&PeakScoreCell,m_PMode == ePBPMConsensus ? 0 : m_MaxPBSeqLen);
#else
			pPeakMatchesCell = pThreadPar->pSW->Align(NULL,m_PMode == ePBPMConsensus ? 0 : m_MaxPBSeqLen);
#endif
			ProvSWchecked += 1;
			if(pPeakMatchesCell != NULL && pPeakMatchesCell->NumMatches >= MinOverlapLen)
				{
				ProvOverlapping += 1;
				PeakMatchesCell = *pPeakMatchesCell;
				ProbeAlignLength = PeakMatchesCell.EndPOfs - PeakMatchesCell.StartPOfs + 1;
				TargAlignLength = PeakMatchesCell.EndTOfs - PeakMatchesCell.StartTOfs + 1;

				// characterise the overlapped target
				// eOLCOverlapping if probe accepted as overlapping, either 5' or 3'
				// eOLCcontaining if both ends of target completely contained within probe
                // eOLCartefact if target is only partially contained
				int PathClass;
				Class = (int)eOLCOverlapping;		
				if((PeakMatchesCell.StartTOfs >= m_OverlapFloat &&  PeakMatchesCell.StartPOfs >= m_OverlapFloat) ||
					 ((TargSeqLen - PeakMatchesCell.EndTOfs) >= m_OverlapFloat && (pCurPBScaffNode->SeqLen - PeakMatchesCell.EndPOfs) >= m_OverlapFloat))
					{
					Class = (int)eOLCartefact;
					ProvArtefact += 1;
					}

				if(m_PMode != ePBPMConsensus && Class == eOLCOverlapping && (PathClass = pThreadPar->pSW->ClassifyPath(m_MaxArtefactDev,
																					PeakMatchesCell.PFirstAnchorStartOfs,PeakMatchesCell.PLastAnchorEndOfs,
																					PeakMatchesCell.TFirstAnchorStartOfs,PeakMatchesCell.TLastAnchorEndOfs)) > 0)
					{
					Class = (int)eOLCartefact;
					ProvArtefact += 1;
					}

				if(Class == eOLCOverlapping)
					{
					if(PeakMatchesCell.StartTOfs < m_OverlapFloat && (TargSeqLen - PeakMatchesCell.EndTOfs) < m_OverlapFloat) // is target completely contained by probe?
						Class = (int)eOLCcontains;
					else
						if(PeakMatchesCell.StartPOfs < m_OverlapFloat && (pCurPBScaffNode->SeqLen - PeakMatchesCell.EndPOfs) < m_OverlapFloat) // or is probe completely contained within target?
							Class = (int)eOLCcontained;
					if(Class != eOLCOverlapping)
						{
						AcquireLock(true);
						if(Class == (int)eOLCcontains)
							pTargNode->flgContains = 1;
						else
							pTargNode->flgContained = 1;
						ReleaseLock(true);	
						ProvContained += 1;
						}
					}
				
				if(m_PMode == ePBPMErrCorrect && Class != (int)eOLCartefact && pThreadPar->NumTargCoreHitCnts >= 2)
					{
#undef _CHECKKMERDIST_
#ifdef _CHECKKMERDIST_	
// define _CHECKKMERDIST_ if wanting to see the exact K-mer length distributions when debugging
					UINT32 CurExactKmerDists[100];
					pThreadPar->pSW->PathKmerCnts(70,CurExactKmerDists,PeakMatchesCell.PFirstAnchorStartOfs,PeakMatchesCell.PLastAnchorEndOfs,
																									PeakMatchesCell.TFirstAnchorStartOfs,PeakMatchesCell.TLastAnchorEndOfs);	
				
					char szKMerDists[2000];
					int KMerDistsIdx;
					AcquireSerialise();
					for(int TTx = 0; TTx < 70; TTx += 1)
						m_ExactKmerDists[TTx] += CurExactKmerDists[TTx];
					m_TotAlignSeqLen += (PeakMatchesCell.PLastAnchorEndOfs - PeakMatchesCell.PFirstAnchorStartOfs);	
					m_TotAlignSeqs += 1;
					if(m_TotAlignSeqs % 100 == 99)
						{
						KMerDistsIdx = sprintf(szKMerDists,"ExactKMerDists over %d seq pairs with mean length %dbp ",m_TotAlignSeqs,m_TotAlignSeqLen/m_TotAlignSeqs);
						for(int TTx = 0; TTx < 70; TTx += 1)
							KMerDistsIdx += sprintf(&szKMerDists[KMerDistsIdx],",%d",m_ExactKmerDists[TTx]);
						gDiagnostics.DiagOut(eDLInfo,gszProcName,szKMerDists);
						}
					ReleaseSerialise();
#endif
					pThreadPar->pSW->TracebacksToAlignOps(PeakMatchesCell.PFirstAnchorStartOfs,PeakMatchesCell.PLastAnchorEndOfs,
															PeakMatchesCell.TFirstAnchorStartOfs,PeakMatchesCell.TLastAnchorEndOfs);
															

					pThreadPar->pSW->AddMultiAlignment(PeakMatchesCell.PFirstAnchorStartOfs,PeakMatchesCell.PLastAnchorEndOfs,
															PeakMatchesCell.TFirstAnchorStartOfs,PeakMatchesCell.TLastAnchorEndOfs,pThreadPar->pTargSeq);
					NumInMultiAlignment += 1;
					}

				if(m_PMode == ePBPMOverlapDetail && m_hErrCorFile != -1)
					{
					m_pSfxArray->GetIdentName(pTargNode->EntryID,sizeof(szTargSeqName)-1,szTargSeqName);
					AcquireSerialise();
					if(m_ScaffLineBuffIdx > (sizeof(m_szScaffLineBuff) - 1000))
						{
						CUtility::SafeWrite(m_hErrCorFile,m_szScaffLineBuff,m_ScaffLineBuffIdx);
						m_ScaffLineBuffIdx = 0;
						}

					m_ScaffLineBuffIdx += sprintf(&m_szScaffLineBuff[m_ScaffLineBuffIdx], "\n%d,%d,\"%s\",%d,\"%s\",%d,\"%c\",\"%c\",%1.5d,%1.5d,%1.5d,%1.5d,%1.5d,%1.5d,%1.5d,%1.5d,%1.5d,%1.5d,%1.5d,%1.5d,%1.5d,%1.5d,%1.5d,%1.5d,%1.5d,%1.5d,%1.5d,%1.5d",
							Class,pCurPBScaffNode->EntryID,szProbeSeqName,pTargNode->EntryID,szTargSeqName,
																			bTargSense ? pSummaryCnts->NumSHits : pSummaryCnts->NumAHits, 'S', bTargSense ? 'S' : 'A',
																			pCurPBScaffNode->SeqLen,TargSeqLen,
							ProbeAlignLength, 
							TargAlignLength,
							PeakMatchesCell.PeakScore,
							PeakMatchesCell.CurScore,
							PeakMatchesCell.NumMatches,
							PeakMatchesCell.NumExacts,
							PeakMatchesCell.NumGapsIns,
							PeakMatchesCell.NumBasesIns,
							PeakMatchesCell.NumGapsDel,
							PeakMatchesCell.NumBasesDel,
							PeakMatchesCell.StartPOfs,
							PeakMatchesCell.StartTOfs,
							PeakMatchesCell.EndPOfs,
							PeakMatchesCell.EndTOfs,
							PeakMatchesCell.PFirstAnchorStartOfs,
							PeakMatchesCell.TFirstAnchorStartOfs,
							PeakMatchesCell.PLastAnchorEndOfs,
							PeakMatchesCell.TLastAnchorEndOfs);
#ifdef _PEAKSCOREACCEPT_
					m_ScaffLineBuffIdx += sprintf(&m_szScaffLineBuff[m_ScaffLineBuffIdx], ",%1.5d,%1.5d,%1.5d,%1.5d,%1.5d,%1.5d,%1.5d,%1.5d,%1.5d,%1.5d,%1.5d,%1.5d,%1.5d,%1.5d,%1.5d,%1.5d",
							PeakScoreCell.NumMatches + PeakScoreCell.NumBasesIns - PeakScoreCell.NumBasesDel, 
							PeakScoreCell.PeakScore,
							PeakScoreCell.NumMatches,
							PeakScoreCell.NumExacts,
							PeakScoreCell.NumGapsIns,
							PeakScoreCell.NumBasesIns,
							PeakScoreCell.NumGapsDel,
							PeakScoreCell.NumBasesDel,
							PeakScoreCell.StartPOfs,
							PeakScoreCell.StartTOfs,
							PeakMatchesCell.EndPOfs,
							PeakMatchesCell.EndTOfs,
							PeakScoreCell.PFirstAnchorStartOfs,
							PeakScoreCell.TFirstAnchorStartOfs,
							PeakScoreCell.PLastAnchorEndOfs,
							PeakScoreCell.TLastAnchorEndOfs);
#endif
#ifdef _EXACTMATCHLENDIST_
					for(int ExactLenIdx = 0; ExactLenIdx < cExactMatchLenDist; ExactLenIdx++)
						m_ScaffLineBuffIdx+=sprintf(&m_szScaffLineBuff[m_ScaffLineBuffIdx],",%d",PeakMatchesCell.ExactMatchLenDist[ExactLenIdx]);
#endif
					ReleaseSerialise();
					}
				ProvOverlapped += 1;
				}
			}
		}

	if(m_PMode == ePBPMErrCorrect && pThreadPar->pSW != NULL && NumInMultiAlignment >= 1 &&  (m_hMultiAlignFile != -1 || m_hErrCorFile != -1))
		{
		pThreadPar->pSW->GenMultialignConcensus();
		AcquireSerialise();
		if(m_hErrCorFile != -1)
			{
			if((pThreadPar->ErrCorBuffIdx=pThreadPar->pSW->MAlignCols2fasta(pCurPBScaffNode->EntryID,m_MinConcScore,m_MinErrCorrectLen,pThreadPar->AllocdErrCorLineBuff,pThreadPar->pszErrCorLineBuff)) > 0)
				{
				CUtility::SafeWrite(m_hErrCorFile,pThreadPar->pszErrCorLineBuff,pThreadPar->ErrCorBuffIdx);
				pThreadPar->ErrCorBuffIdx = 0;
				}
			}

		if(m_hMultiAlignFile != -1)
			{
			if((pThreadPar->MultiAlignBuffIdx=pThreadPar->pSW->MAlignCols2MFA(pCurPBScaffNode->EntryID,pThreadPar->AllocdMultiAlignLineBuff,pThreadPar->pszMultiAlignLineBuff)) > 0)
				{
				pThreadPar->pszMultiAlignLineBuff[pThreadPar->MultiAlignBuffIdx++] = '\n';
				CUtility::SafeWrite(m_hMultiAlignFile,pThreadPar->pszMultiAlignLineBuff,pThreadPar->MultiAlignBuffIdx);
				pThreadPar->MultiAlignBuffIdx = 0;
				}
			}
		ReleaseSerialise();
		}
	AcquireSerialise();
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
	ReleaseSerialise();
	ProvOverlapping = 0;
	ProvOverlapped = 0;
	ProvContained = 0;
	ProvArtefact = 0;
	ProvSWchecked = 0;
	NumInMultiAlignment = 0;
	pThreadPar->NumTargCoreHitCnts = 0;
	pThreadPar->NumCoreHits = 0;
	}
AcquireSerialise();
if(m_PMode == ePBPMErrCorrect && m_hErrCorFile != -1 && pThreadPar->ErrCorBuffIdx > 0)
	{
	CUtility::SafeWrite(m_hErrCorFile,pThreadPar->pszErrCorLineBuff,pThreadPar->ErrCorBuffIdx);
	pThreadPar->ErrCorBuffIdx = 0;
	}
if(m_hMultiAlignFile != -1 && pThreadPar->MultiAlignBuffIdx > 0)
	{
	CUtility::SafeWrite(m_hMultiAlignFile,pThreadPar->pszMultiAlignLineBuff,pThreadPar->MultiAlignBuffIdx);
	pThreadPar->MultiAlignBuffIdx = 0;
	}
if(m_PMode == ePBPMOverlapDetail && m_ScaffLineBuffIdx > 0)
	{
	CUtility::SafeWrite(m_hErrCorFile,m_szScaffLineBuff,m_ScaffLineBuffIdx);
	m_ScaffLineBuffIdx = 0;
	}
ReleaseSerialise();

if(pThreadPar->pSW != NULL)
	{
	delete pThreadPar->pSW;
	pThreadPar->pSW = NULL;
	}
return(0);
}


int					// returns index 1..N of just added core hit or -1 if errors
CPBErrCorrect::AddCoreHit(UINT32 ProbeNodeID,			// core hit was from this probe scaffold node 
			   bool bRevCpl,					// true if core sequence was revcpl'd before matching
			   UINT32 ProbeOfs,                 // hit started at this probe offset
			   UINT32 TargNodeID,               // probe core matched onto this target scaffold node
			   UINT32 TargOfs,                  // probe core matched starting at this target loci
			   UINT32 HitLen,					// hit was of this length
               tsThreadPBErrCorrect *pPars)			// thread specific
{
tsCoreHit *pCoreHit;

if((pPars->NumCoreHits + 5) > pPars->AllocdCoreHits)	// need to realloc memory to hold additional cores?
	{
		// realloc memory with a 25% increase over previous allocation 
	int coresreq;
	size_t memreq;
	void *pAllocd;
	coresreq = (int)(((INT64)pPars->AllocdCoreHits * 125) / (INT64)100);
	memreq = coresreq * sizeof(tsCoreHit);

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
			Reset(false);
			return(eBSFerrMem);
			}

		pPars->pCoreHits = (tsCoreHit *)pAllocd;
		pPars->AllocdCoreHitsSize = memreq;
		pPars->AllocdCoreHits = coresreq; 
		}
		
pCoreHit = &pPars->pCoreHits[pPars->NumCoreHits++];

pCoreHit->ProbeNodeID = ProbeNodeID;
pCoreHit->flgRevCpl = bRevCpl ? 1 : 0;
pCoreHit->flgMulti = 0;
pCoreHit->ProbeOfs = ProbeOfs;
pCoreHit->TargNodeID = TargNodeID;
pCoreHit->HitLen = HitLen;
pCoreHit->TargOfs = TargOfs;
memset(&pCoreHit[1],0,sizeof(tsCoreHit));	// ensuring that used cores are always terminated with a marker end of cores initialised to 0
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
UINT32 ProbeOfs;
UINT32 LastProbeOfs;

UINT32 ChkOvrLapCoreProbeOfs;
UINT32 LastCoreProbeOfs;
int ChkOvrLapCoreStartIdx;

etSeqBase *pCoreSeq;
tsPBScaffNode *pProbeNode;
tsPBScaffNode *pTargNode;
etSeqBase *pHomo;
int HomoIdx;
int HomoBaseCnts[4];
int MaxAcceptHomoCnt;

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
MaxAcceptHomoCnt = (pPars->CoreSeqLen * cQualCoreHomopolymer) / 100; // if any core contains more than cQualCoreHomopolymer% of the same base then treat as being a near homopolymer core and slough
LastProbeOfs = 1 + pProbeNode->SeqLen - pPars->CoreSeqLen;
if(pPars->CoreSeqLen < cMaxPacBioSeedExtn)
	LastProbeOfs -= 120;


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

    while((NextHitIdx = m_pSfxArray->IteratePacBio(pCoreSeq,pProbeNode->SeqLen - ProbeOfs,pPars->CoreSeqLen,pProbeNode->EntryID,pPars->MinPBSeqLen,PrevHitIdx,&HitEntryID,&HitLoci)) > 0)
		{
		PrevHitIdx = NextHitIdx;
		pTargNode = &m_pPBScaffNodes[MapEntryID2NodeID(HitEntryID)-1];
		HitSeqLen = pTargNode->SeqLen; 
		AcquireLock(false);
		if((pPars->bSelfHits ? pTargNode->NodeID != pProbeNode->NodeID : pTargNode->NodeID == pProbeNode->NodeID) || 
				 pTargNode->flgUnderlength == 1 ||	// not interested in selfhits or underlength targets
							HitSeqLen < (UINT32)pPars->MinPBSeqLen)		// not interested if target sequence length less than min sequence length to be processed
			{
			ReleaseLock(false);
			continue;
			}
		ReleaseLock(false);

  		AddCoreHit(ProbeNodeID,pPars->bRevCpl,ProbeOfs,pTargNode->NodeID,HitLoci,pPars->CoreSeqLen,pPars);
		HitsThisCore += 1;
		if(HitsThisCore > pPars->MaxAcceptHitsPerSeedCore)
			break;

		}
	if(HitsThisCore)	// if at least one hit from this core
		{
		if(HitsThisCore > HighHitsThisCore)
			HighHitsThisCore = HitsThisCore;
		TotHitsAllCores += HitsThisCore;
		}
	}
if(pPars->bRevCpl)
	CSeqTrans::ReverseComplement(pProbeNode->SeqLen,pPars->pProbeSeq);

return(pPars->NumCoreHits);
}


// SortLenDescending
// Sort scaffolding nodes by length descending with entry identifiers as tie breaker
int
CPBErrCorrect::SortLenDescending(const void *arg1, const void *arg2)
{
tsPBScaffNode *pEl1 = (tsPBScaffNode *)arg1;
tsPBScaffNode *pEl2 = (tsPBScaffNode *)arg2;

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
tsCoreHit *pEl1 = (tsCoreHit *)arg1;
tsCoreHit *pEl2 = (tsCoreHit *)arg2;

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
tsCoreHit *pEl1 = (tsCoreHit *)arg1;
tsCoreHit *pEl2 = (tsCoreHit *)arg2;

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
sPBCoreHitCnts *pEl1 = (sPBCoreHitCnts *)arg1;
sPBCoreHitCnts *pEl2 = (sPBCoreHitCnts *)arg2;

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





