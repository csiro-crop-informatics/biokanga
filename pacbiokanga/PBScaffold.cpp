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
#include "../libbiokanga/SSW.h"
#include "pacbiocommon.h"
#include "AssembGraph.h"
#include "PBScaffold.h"


int
ProcPacBioScaffolds(etPBPMode PMode,		// processing mode
		int DeltaCoreOfs,			// offset by this many bp the core windows of coreSeqLen along the probe sequence when checking for overlaps
		int MaxSeedCoreDepth,		// only further process a seed core if there are no more than this number of matching cores in all targeted sequences
		int MinSeedCoreLen,			// use seed cores of this length when identifying putative overlapping scaffold sequences
		int MinNumSeedCores,        // require at least this many seed cores between overlapping scaffold sequences
		int SWMatchScore,			// score for matching bases (0..50)
		int SWMismatchPenalty,		// mismatch penalty (-50..0)
		int SWGapOpenPenalty,		// gap opening penalty (-50..0)
		int SWGapExtnPenalty,		// gap extension penalty (-50..0)
		int SWProgExtnPenaltyLen,	// progressive gap scoring then only apply gap extension score if gap at least this length (0..63) - use if aligning PacBio
		int MinScaffoldLen,			// individual scaffold sequences must be of at least this length (defaults to 5Kbp)
		int MinScaffOverlap,		// pairs of targeted scaffold sequences must overlap by at least this many bp to be considered for merging into a longer scaffold sequence (defaults to 5Kbp) 
		char *pszPacBioOvlps,		// pregenerated PacBio sequence overlap loci details
		char *pszPacBioSfxFile,		// suffix indexed PacBio sequences
		int NumPacBioFiles,			// number of input pacbio file specs
		char *pszPacBioFiles[],		// input pacbio files
	    char *pszOutFile,			// where to write merged scaffolded sequences
		int NumThreads);			// maximum number of worker threads to use


#ifdef _WIN32
int ProcScaffold(int argc, char* argv[])
{
// determine my process name
_splitpath(argv[0],NULL,NULL,gszProcName,NULL);
#else
int
ProcScaffold(int argc, char** argv)
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
int MinSeedCoreLen;			// use seed cores of this length when identifying putative overlapping scaffold sequences
int MinNumSeedCores;        // require at least this many seed cores between overlapping scaffold sequences

int DeltaCoreOfs;			// offset by this many bp the core windows of coreSeqLen along the probe sequence when checking for overlaps
int MaxSeedCoreDepth;		// only further process a seed core if there are no more than this number of matching cores in all targeted sequences


int SWMatchScore;			// score for matching bases (0..50)
int SWMismatchPenalty;		// mismatch penalty (0..50)
int SWGapOpenPenalty;		// gap opening penalty (0..50)
int SWGapExtnPenalty;		// gap extension penalty (0..50)
int SWProgExtnPenaltyLen;	// progressive gap scoring then only apply gap extension score if gap at least this length (0..63) - use if aligning PacBio

int MinScaffoldLen;			// individual target scaffold sequences must be of at least this length (defaults to 5Kbp)
int MinScaffOverlap;		// pairs of targeted scaffold sequences must overlap by at least this many bp to be considered for merging into a longer scaffold sequence (defaults to 5Kbp) 

char szPacBioSfxFile[_MAX_PATH];		// suffix indexed PacBio sequences
char szPacBioOvlps[_MAX_PATH];			// pregenerated PacBio sequence overlap loci details
int NumPacBioFiles;						// number of input pacbio file specs
char *pszPacBioFiles[cMaxInFileSpecs];	// input pacbio files

char szOutFile[_MAX_PATH];				// where to write merged scaffolded sequences

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

struct arg_int *pmode = arg_int0("m","pmode","<int>",				"processing mode - 0 generate scaffolds from error corrected sequences, 1 use pregenerated overlap loci details file ");
struct arg_int *minscafflen = arg_int0("l","minscafflen","<int>",	"minimum individual sequence length (default 5000, range 500 to 100000)");
struct arg_int *minscaffovl = arg_int0("L","minscaffovl","<int>",	"minimum scaffold overlap required to merge into single scaffold (default 5000, range 250 to 100000)");

struct arg_int *minseedcorelen = arg_int0("c","seedcorelen","<int>",			"use seed cores of this length when identifying putative overlapping scaffold sequences (default 12, range 10 to 25)");
struct arg_int *minnumseedcores = arg_int0("C","minseedcores","<int>",			"require at least this many accepted seed cores between overlapping scaffold sequences (default 10, range 1 to 50)");

struct arg_int *deltacoreofs = arg_int0("d","deltacoreofs","<int>",							"offset cores (default 1, range 1 to 10)");
struct arg_int *maxcoredepth = arg_int0("D","maxcoredepth","<int>",							"explore cores of less than this maximum depth (default 10000, range 1000 to 20000)");

struct arg_int *matchscore = arg_int0("x","matchscore","<int>",								"SW score for matching bases (default 3, range 1 to 50)");
struct arg_int *mismatchpenalty = arg_int0("X","mismatchpenalty","<int>",					"SW mismatch penalty (default 7, range 1 to 50)");
struct arg_int *gapopenpenalty = arg_int0("y","gapopenpenalty","<int>",						"SW gap opening penalty (default 4, range 1 to 50)");
struct arg_int *gapextnpenalty = arg_int0("Y","gapextnpenalty","<int>",						"SW gap extension penalty (default 1, range 1 to 50)");
struct arg_int *progextnpenaltylen = arg_int0("z","progextnpenaltylen","<int>",				"SW gap extension penalty only applied for gaps of at least this number of bases (default 2, range 1 to 63)");
struct arg_file *pacbiosovlps = arg_file0("p","pacbiosovlps","<file>",						"input file containing pregenerated PacBio overlap loci detail");

struct arg_file *pacbiosfxfile = arg_file0("I","pacbiosfx","<file>",						"input file containing suffix indexed PacBio sequences");
struct arg_file *pacbiofiles = arg_filen("i","pacbiofile","<file>",0,cMaxInFileSpecs,		"names of input files containing PacBio sequences to be used for scaffolding");
struct arg_file *outfile = arg_file1("o","out","<file>",			"output merged scaffold sequences to this file");

struct arg_int *threads = arg_int0("T","threads","<int>",		"number of processing threads 0..128 (defaults to 0 which sets threads to number of CPU cores)");

struct arg_file *summrslts = arg_file0("q","sumrslts","<file>",				"Output results summary to this SQLite3 database file");
struct arg_str *experimentname = arg_str0("w","experimentname","<str>",		"experiment name SQLite3 database file");
struct arg_str *experimentdescr = arg_str0("W","experimentdescr","<str>",	"experiment description SQLite3 database file");

struct arg_end *end = arg_end(200);

void *argtable[] = {help,version,FileLogLevel,LogFile,
					pmode,minseedcorelen,minnumseedcores,deltacoreofs,maxcoredepth,
					matchscore,mismatchpenalty,gapopenpenalty,gapextnpenalty,progextnpenaltylen,
					minscafflen,minscaffovl,
					summrslts,pacbiosovlps,pacbiosfxfile,pacbiofiles,experimentname,experimentdescr,
					outfile,threads,
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

	PMode = (etPBPMode)(pmode->count ? pmode->ival[0] : (int)ePBPMOverlaps);
	if(PMode < ePBPMOverlaps || PMode > ePBPMScaffold)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Processing mode '-m%d' must be in range 0..%d",PMode,ePBPMScaffold);
		return(1);
		}

	MinSeedCoreLen = minseedcorelen->count ? minseedcorelen->ival[0] : cDfltSeedCoreLen;
	if(MinSeedCoreLen < cMinSeedCoreLen || MinSeedCoreLen > cMaxSeedCoreLen)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: seed core length '-c%d' must be in range %d..%d",MinSeedCoreLen,cMinSeedCoreLen,cMaxSeedCoreLen);
		return(1);
		}

	MinNumSeedCores = minnumseedcores->count ? minnumseedcores->ival[0] : cDfltNumSeedCores;
	if(MinNumSeedCores < cMinNumSeedCores || MinNumSeedCores > cMaxNumSeedCores)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Minimum number of accepted seed cores '-C%d' must be in range %d..%d",MinNumSeedCores,cMinNumSeedCores,cMaxNumSeedCores);
		return(1);
		}

	DeltaCoreOfs = deltacoreofs->count ? deltacoreofs->ival[0] : 1;
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

	MinScaffoldLen = minscafflen->count ? minscafflen->ival[0] : cDfltMinScaffSeqLen;
	if(MinScaffoldLen < cMinScaffSeqLen || MinScaffoldLen > cMaxMinScaffSeqLen)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Minimum accepted scaffold length '-l%d' must be in range %d..%dbp",MinScaffoldLen,cMinScaffSeqLen,cMaxMinScaffSeqLen);
		return(1);
		}

	MinScaffOverlap = minscaffovl->count ? minscaffovl->ival[0] : cDfltMinOverlapLen;
	if(MinScaffOverlap < cMinOverlapLen || MinScaffOverlap > cMaxMinScaffSeqLen)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Minimum scaffold overlap required length '-L%d' must be in range %d..%dbp",MinScaffOverlap,MinScaffOverlap,cMaxMinScaffSeqLen);
		return(1);
		}

	if(pacbiosfxfile->count && pacbiofiles->count)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: input PacBio suffix file or PacBio sequences are exclusive, specify one only");
		return(1);
		}

	if(!(pacbiosfxfile->count || pacbiofiles->count))
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: either PacBio suffix file or PacBio sequences file must be specified, specify one only");
		return(1);
		}

	if(pacbiosfxfile->count)
		{
		strncpy(szPacBioSfxFile,pacbiosfxfile->filename[0],sizeof(szPacBioSfxFile));
		szPacBioSfxFile[sizeof(szPacBioSfxFile)-1] = '\0';
		CUtility::TrimQuotedWhitespcExtd(szPacBioSfxFile);
		if(szPacBioSfxFile[0] == '\0')
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: After removal of whitespace, no input suffix indexed PacBio file specified with '-I<filespec>' option)\n");
			exit(1);
			}
		}
	else
		szPacBioSfxFile[0] = '\0';

	if(pacbiofiles->count)
		{
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
		}
	else
		{
		NumPacBioFiles = 0;
		pszPacBioFiles[0] = NULL;
		}

	szPacBioOvlps[0] = '\0';
	if(PMode == ePBPMScaffold)
		{
		if(pacbiosovlps->count)
			{
			strncpy(szPacBioOvlps,pacbiosovlps->filename[0],sizeof(szOutFile));
			szPacBioOvlps[sizeof(szPacBioOvlps)-1] = '\0';
			CUtility::TrimQuotedWhitespcExtd(szPacBioOvlps);
			}
		if(szPacBioOvlps[0] == '\0')
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: After removal of whitespace, no input pregenerated PacBio sequence overlap loci detail file specified with '-p<filespec>' option)\n");
			exit(1);
			}
		}

	strncpy(szOutFile,outfile->filename[0],sizeof(szOutFile));
	szOutFile[sizeof(szOutFile)-1] = '\0';
	CUtility::TrimQuotedWhitespcExtd(szOutFile);
	if(szOutFile[0] == '\0')
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: After removal of whitespace, no output file specified with '-o<filespec>' option)\n");
		exit(1);
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
		case ePBPMOverlaps:		// overlap discovery
			pszMode = (char *)"Identify overlaps for subsequent scaffolding";
			break;
		case ePBPMScaffold:		// scaffolding
			pszMode = (char *)"Using pregenerated overlaps for scaffolding";
			break;

			}

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Processing mode: '%s'",pszMode);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Use seed cores of this length when identifying putative overlapping scaffold sequences: %dbp",MinSeedCoreLen);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Require at least this many seed cores between overlapping scaffold sequences: %d",MinNumSeedCores);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Offset cores by this many bp: %d",DeltaCoreOfs);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Maximum seed core depth: %d",MaxSeedCoreDepth);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"SW score for matching bases: %d",SWMatchScore);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"SW mismatch penalty: %d",SWMismatchPenalty);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"SW gap opening penalty: %d",SWGapOpenPenalty);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"SW gap extension penalty: %d",SWGapExtnPenalty);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"SW gap extension penalty only applied for gaps of at least this size: %d",SWProgExtnPenaltyLen);


	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Minimum individual scaffolding sequence length: %dbp",MinScaffoldLen);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Minimum scaffold overlap required to merge into single scaffold: %d",MinScaffOverlap);

	if(PMode == ePBPMScaffold)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Input file containing pregenerated PacBio overlap loci detail: '%s'",szPacBioOvlps);

	if(NumPacBioFiles)
		{
		for(Idx=0; Idx < NumPacBioFiles; Idx++)
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Input PacBio sequences file spec: '%s'",pszPacBioFiles[Idx]);
		}
	else
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Input PacBio suffix array file: '%s'",szPacBioSfxFile);

	if(PMode == ePBPMScaffold)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Output merged scaffold sequences file: '%s'",szOutFile);
	else
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Output PacBio overlap loci detail file: '%s'",szOutFile);

	if(szExperimentName[0] != '\0')
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"This processing reference: %s",szExperimentName);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"number of threads : %d",NumThreads);

	if(gExperimentID > 0)
		{
		int ParamID;
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szLogFile),"log",szLogFile);

		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(PMode),"pmode",&PMode);

		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(MinSeedCoreLen),"seedcorelen",&MinSeedCoreLen);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(MinNumSeedCores),"minseedcores",&MinNumSeedCores);

		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(DeltaCoreOfs),"deltacoreofs",&DeltaCoreOfs);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(MaxSeedCoreDepth),"maxcoredepth",&MaxSeedCoreDepth);


		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(MinScaffoldLen),"minscafflen",&MinScaffoldLen);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(MinScaffOverlap),"minscaffovl",&MinScaffOverlap);
		if(NumPacBioFiles)
			{
			for(Idx=0; Idx < NumPacBioFiles; Idx++)
				ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(pszPacBioFiles[Idx]),"pacbiofile",pszPacBioFiles[Idx]);
			}
		else
			ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szPacBioSfxFile),"pacbiosfx",szPacBioSfxFile);

		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szOutFile),"out",szOutFile);

		if(PMode == ePBPMScaffold)
			ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szPacBioOvlps),"pacbiosovlps",szPacBioOvlps);
		
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
	Rslt = ProcPacBioScaffolds((etPBPMode)PMode,DeltaCoreOfs,MaxSeedCoreDepth,MinSeedCoreLen,MinNumSeedCores,SWMatchScore,-1 * SWMismatchPenalty,-1 * SWGapOpenPenalty,-1 * SWGapExtnPenalty,SWProgExtnPenaltyLen,MinScaffoldLen,MinScaffOverlap,
						szPacBioOvlps,szPacBioSfxFile,NumPacBioFiles,pszPacBioFiles,szOutFile,NumThreads);
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
ProcPacBioScaffolds(etPBPMode PMode,		// processing mode
		int DeltaCoreOfs,			// offset by this many bp the core windows of coreSeqLen along the probe sequence when checking for overlaps
		int MaxSeedCoreDepth,		// only further process a seed core if there are no more than this number of matching cores in all targeted sequences
		int MinSeedCoreLen,			// use seed cores of this length when identifying putative overlapping scaffold sequences
		int MinNumSeedCores,        // require at least this many seed cores between overlapping scaffold sequences
		int SWMatchScore,			// score for matching bases (0..50)
		int SWMismatchPenalty,		// mismatch penalty (-50..0)
		int SWGapOpenPenalty,		// gap opening penalty (-50..0)
		int SWGapExtnPenalty,		// gap extension penalty (-50..0)
		int SWProgExtnPenaltyLen,	// only apply gap extension score if gap at least this length (0..63) - use if aligning PacBio
		int MinScaffoldLen,			// individual target scaffold sequences must be of at least this length (defaults to 5Kbp)
		int MinScaffOverlap,		// pairs of targeted scaffold sequences must overlap by at least this many bp to be considered for merging into a longer scaffold sequence (defaults to 5Kbp) 
		char *pszPacBioOvlps,		// pregenerated PacBio sequence overlap loci details
		char *pszPacBioSfxFile,		// suffix indexed PacBio sequences
		int NumPacBioFiles,			// number of input pacbio file specs
		char *pszPacBioFiles[],		// input pacbio files
		char *pszOutFile,			// where to write merged scaffolded sequences
		int NumThreads)				// maximum number of worker threads to use
{
int Rslt;
CPBScaffold *pPacBioer;

if((pPacBioer = new CPBScaffold)==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: Unable to instantiate CPBScaffold");
	return(eBSFerrObj);
	}

Rslt = pPacBioer->Process(PMode,DeltaCoreOfs,MaxSeedCoreDepth,MinSeedCoreLen,MinNumSeedCores,SWMatchScore,SWMismatchPenalty,SWGapOpenPenalty,SWGapExtnPenalty,SWProgExtnPenaltyLen,MinScaffoldLen,MinScaffOverlap,
								pszPacBioOvlps,pszPacBioSfxFile,NumPacBioFiles,pszPacBioFiles,pszOutFile,NumThreads);
delete pPacBioer;
return(Rslt);
}

CPBScaffold::CPBScaffold() // relies on base classes constructors
{
m_pSfxArray = NULL;
m_pAssembGraph = NULL;
m_pPBScaffNodes = NULL;
m_pMapEntryID2NodeIDs = NULL;
m_bMutexesCreated = false;
m_hScaffDist = -1;
Init();
}

CPBScaffold::~CPBScaffold() // relies on base classes destructors
{
Reset(false);
}


void
CPBScaffold::Init(void)
{
if(m_pSfxArray != NULL)
	{
	delete m_pSfxArray;
	m_pSfxArray = NULL;
	}

if(m_pAssembGraph != NULL)
	{
	delete m_pAssembGraph;
	m_pAssembGraph = NULL;
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

if(m_hScaffDist != -1)
	{
#ifdef _WIN32
	_commit(m_hScaffDist);
#else
	fsync(m_hScaffDist);
#endif
	close(m_hScaffDist);
	m_hScaffDist = -1;
	}

m_NumPBScaffNodes = 0;
m_AllocdPBScaffNodes = 0;

m_NumOverlapProcessed = 0;
m_ProvOverlapping = 0;
m_ProvOverlapped = 0;
m_ProvContained = 0;
m_ProvArtefact = 0;
m_ProvSWchecked = 0;

m_PMode = ePBPMOverlaps;

m_OverlapFloat = cDfltMaxOverlapFloat;
m_MinScaffSeqLen = cDfltMinScaffSeqLen;	
m_MinScaffOverlap = cDfltMinOverlapLen;

m_DeltaCoreOfs = cDfltDeltaCoreOfs;
m_MaxSeedCoreDepth = cDfltMaxSeedCoreDepth;
m_MinSeedCoreLen = cDfltSeedCoreLen;
m_MinNumSeedCores = cDfltNumSeedCores;

m_SWMatchScore = cDfltSWMatchScore;
m_SWMismatchPenalty = cDfltSWMismatchPenalty;	
m_SWGapOpenPenalty = cDfltSWGapOpenPenalty;
m_SWGapExtnPenalty = cDfltSWGapExtnPenalty;
m_SWProgExtnPenaltyLen = cDfltSWProgExtnLen;	
 
m_NumPacBioFiles = 0;
m_szPacBioSfxFile[0] = '\0';
memset(m_szPacBioFiles,0,sizeof(m_szPacBioFiles));

m_szOutFile[0] = '\0';	

m_NumThreads = 0;
if(m_bMutexesCreated)
	DeleteMutexes();
m_bMutexesCreated = false; 
}

void
CPBScaffold::Reset(bool bSync)			// if bSync true then fsync before closing output file handles
{

Init();
}


int
CPBScaffold::CreateMutexes(void)
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
CPBScaffold::DeleteMutexes(void)
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
CPBScaffold::AcquireSerialise(void)
{
#ifdef _WIN32
WaitForSingleObject(m_hMtxIterReads,INFINITE);
#else
pthread_mutex_lock(&m_hMtxIterReads);
#endif
}

void
CPBScaffold::ReleaseSerialise(void)
{
#ifdef _WIN32
ReleaseMutex(m_hMtxIterReads);
#else
pthread_mutex_unlock(&m_hMtxIterReads);
#endif
}

void
CPBScaffold::AcquireSerialiseMH(void)
{
#ifdef _WIN32
WaitForSingleObject(m_hMtxMHReads,INFINITE);
#else
pthread_mutex_lock(&m_hMtxMHReads);
#endif
}

void
CPBScaffold::ReleaseSerialiseMH(void)
{
#ifdef _WIN32
ReleaseMutex(m_hMtxMHReads);
#else
pthread_mutex_unlock(&m_hMtxMHReads);
#endif
}

void
CPBScaffold::AcquireLock(bool bExclusive)
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
CPBScaffold::ReleaseLock(bool bExclusive)
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


UINT32												//  returns number of overlaps loaded and accepted, if > cMaxValidID then cast to teBSFrsltCodes for actual error 
CPBScaffold::LoadPacBioOvlps(char *pszPacBioOvlps,  // load from this pregenerated PacBio sequence overlap loci CSV file
							 bool bValidateOnly)	// if true then simply parse and validate that the CSV format is as expected
{
int Rslt;
int NumFields;
UINT32 NumElsRead;
UINT32 NumEdges;

char *pszTmp;
int Class;
int ProbeID;
char szProbeDescr[cMaxDescrIDLen+1];
int TargID;
char szTargDescr[cMaxDescrIDLen+1];
int SeedHits;
char Sense;
int ProbeLen;
int TargLen;
int	ProbeAlignLength;
int	TargAlignLength;
int	PeakScore;
int FinalScore;
int	NumAlignedBases;
int	NumExactBases;
int	NumProbeInserts;
int	NumProbeInsertBases;
int	NumTargInserts;
int	NumTargInsertBases;
int	ProbeStartOfs;
int	TargStartOfs;
int	ProbeEndOfs;
int	TargEndOfs;
int	ProbeOfs5;
int	TargOfs5;
int	ProbeOfs3;
int	TargOfs3;

if(!bValidateOnly && m_pAssembGraph == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expected m_pAssembGraph to have been instantiated");
	return(eBSFerrObj);
	}

CCSVFile *pCSV = new CCSVFile;
if(pCSV == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CCSVfile");
	return((UINT32)eBSFerrObj);
	}

if((Rslt=pCSV->Open(pszPacBioOvlps))!=eBSFSuccess)
	{
	while(pCSV->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,pCSV->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open file: %s",pszPacBioOvlps);
	delete pCSV;
	return((UINT32)Rslt);
	}
NumElsRead = 0;
NumEdges = 0;
while((Rslt=pCSV->NextLine()) > 0)			// onto next line containing fields
	{
	NumFields = pCSV->GetCurFields();		// PacBio scaffold overlap files contain 27 fields
	if(NumFields != 27)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expected PacBio overlaps CSV file number of fields to be 27 in '%s', GetCurFields() returned '%d'",pszPacBioOvlps,NumFields);
		delete pCSV;
		return((UINT32)eBSFerrFieldCnt);
		}
	if(!NumElsRead && pCSV->IsLikelyHeaderLine())
		continue;
	NumElsRead += 1;

	if((Rslt=pCSV->GetInt(1,&Class)) < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error parsing Class at line %d in file '%s'",NumElsRead,pszPacBioOvlps);
		delete pCSV;
		return((UINT32)Rslt);
		}
	if((Rslt=pCSV->GetInt(2,&ProbeID)) < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error parsing ProbeID at line %d in file '%s'",NumElsRead,pszPacBioOvlps);
		delete pCSV;
		return((UINT32)Rslt);
		};
	if((Rslt=pCSV->GetText(3,&pszTmp)) < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error parsing ProbeDescr at line %d in file '%s'",NumElsRead,pszPacBioOvlps);
		delete pCSV;
		return((UINT32)Rslt);
		};
	strncpy(szProbeDescr,pszTmp,sizeof(szProbeDescr));
	szProbeDescr[sizeof(szProbeDescr)-1] = '\0';
	if((Rslt=pCSV->GetInt(4,&TargID)) < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error parsing TargID at line %d in file '%s'",NumElsRead,pszPacBioOvlps);
		delete pCSV;
		return((UINT32)Rslt);
		};
	if((Rslt=pCSV->GetText(5,&pszTmp)) < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error parsing TargDescr at line %d in file '%s'",NumElsRead,pszPacBioOvlps);
		delete pCSV;
		return((UINT32)Rslt);
		};
	strncpy(szTargDescr,pszTmp,sizeof(szTargDescr));
	szTargDescr[sizeof(szTargDescr)-1] = '\0';
	if((Rslt=pCSV->GetInt(6,&SeedHits)) < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error parsing SeedHits at line %d in file '%s'",NumElsRead,pszPacBioOvlps);
		delete pCSV;
		return((UINT32)Rslt);
		};
	if((Rslt=pCSV->GetText(7,&pszTmp)) < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error parsing Sense at line %d in file '%s'",NumElsRead,pszPacBioOvlps);
		delete pCSV;
		return((UINT32)Rslt);
		};
	Sense = pszTmp[0];
	switch(Sense) {
		case 's': case 'S': break;
		case 'a': case 'A': break;
		default:
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error parsing Sense at line %d in file '%s'",NumElsRead,pszPacBioOvlps);
			delete pCSV;
			return((UINT32)eBSFerrParse);
		};
	if((Rslt=pCSV->GetInt(8,&ProbeLen)) < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error parsing ProbeLen at line %d in file '%s'",NumElsRead,pszPacBioOvlps);
		delete pCSV;
		return((UINT32)Rslt);
		};
	if((Rslt=pCSV->GetInt(9,&TargLen)) < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error parsing TargLen at line %d in file '%s'",NumElsRead,pszPacBioOvlps);
		delete pCSV;
		return((UINT32)Rslt);
		};
	if((Rslt=pCSV->GetInt(10,&ProbeAlignLength)) < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error parsing ProbeAlignLength at line %d in file '%s'",NumElsRead,pszPacBioOvlps);
		delete pCSV;
		return((UINT32)Rslt);
		};
	if((Rslt=pCSV->GetInt(11,&TargAlignLength)) < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error parsing TargAlignLength at line %d in file '%s'",NumElsRead,pszPacBioOvlps);
		delete pCSV;
		return((UINT32)Rslt);
		};
	if((Rslt=pCSV->GetInt(12,&PeakScore)) < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error parsing PeakScore at line %d in file '%s'",NumElsRead,pszPacBioOvlps);
		delete pCSV;
		return((UINT32)Rslt);
		};
	if((Rslt=pCSV->GetInt(13,&FinalScore)) < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error parsing FinalScore at line %d in file '%s'",NumElsRead,pszPacBioOvlps);
		delete pCSV;
		return((UINT32)Rslt);
		};
	if((Rslt=pCSV->GetInt(14,&NumAlignedBases)) < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error parsing NumAlignedBases at line %d in file '%s'",NumElsRead,pszPacBioOvlps);
		delete pCSV;
		return((UINT32)Rslt);
		};
	if((Rslt=pCSV->GetInt(15,&NumExactBases)) < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error parsing NumExactBases at line %d in file '%s'",NumElsRead,pszPacBioOvlps);
		delete pCSV;
		return((UINT32)Rslt);
		};
	if((Rslt=pCSV->GetInt(16,&NumProbeInserts)) < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error parsing NumProbeInserts at line %d in file '%s'",NumElsRead,pszPacBioOvlps);
		delete pCSV;
		return((UINT32)Rslt);
		};
	if((Rslt=pCSV->GetInt(17,&NumProbeInsertBases)) < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error parsing NumProbeInsertBases at line %d in file '%s'",NumElsRead,pszPacBioOvlps);
		delete pCSV;
		return((UINT32)Rslt);
		};
	if((Rslt=pCSV->GetInt(18,&NumTargInserts)) < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error parsing NumTargInserts at line %d in file '%s'",NumElsRead,pszPacBioOvlps);
		delete pCSV;
		return((UINT32)Rslt);
		};
	if((Rslt=pCSV->GetInt(19,&NumTargInsertBases)) < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error parsing NumTargInsertBases at line %d in file '%s'",NumElsRead,pszPacBioOvlps);
		delete pCSV;
		return((UINT32)Rslt);
		};
	if((Rslt=pCSV->GetInt(20,&ProbeStartOfs)) < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error parsing ProbeStartOfs at line %d in file '%s'",NumElsRead,pszPacBioOvlps);
		delete pCSV;
		return((UINT32)Rslt);
		};

	if((Rslt=pCSV->GetInt(21,&TargStartOfs)) < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error parsing TargStartOfs at line %d in file '%s'",NumElsRead,pszPacBioOvlps);
		delete pCSV;
		return((UINT32)Rslt);
		};
	if((Rslt=pCSV->GetInt(22,&ProbeEndOfs)) < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error parsing ProbeEndOfs at line %d in file '%s'",NumElsRead,pszPacBioOvlps);
		delete pCSV;
		return((UINT32)Rslt);
		};
	if((Rslt=pCSV->GetInt(23,&TargEndOfs)) < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error parsing TargEndOfs at line %d in file '%s'",NumElsRead,pszPacBioOvlps);
		delete pCSV;
		return((UINT32)Rslt);
		};
	if((Rslt=pCSV->GetInt(24,&ProbeOfs5)) < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error parsing ProbeOfs5 at line %d in file '%s'",NumElsRead,pszPacBioOvlps);
		delete pCSV;
		return((UINT32)Rslt);
		};
	if((Rslt=pCSV->GetInt(25,&TargOfs5)) < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error parsing TargOfs5 at line %d in file '%s'",NumElsRead,pszPacBioOvlps);
		delete pCSV;
		return((UINT32)Rslt);
		};
	if((Rslt=pCSV->GetInt(26,&ProbeOfs3)) < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error parsing ProbeOfs3 at line %d in file '%s'",NumElsRead,pszPacBioOvlps);
		delete pCSV;
		return((UINT32)Rslt);
		};
	if((Rslt=pCSV->GetInt(27,&TargOfs3)) < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error parsing TargOfs3 at line %d in file '%s'",NumElsRead,pszPacBioOvlps);
		delete pCSV;
		return((UINT32)Rslt);
		};

	// check for consistency in the PacBio read overlaps before edge processing
	// only accepted as consistent if sequence lengths are in the range of 1000..1000000, alignments are at least 100bp, starts+100 < ends, start/ends within read sequence lengths, and SW scores less than 100000000
	// note that all offsets are 1 based 
	if(ProbeID < 1 || TargID < 1 || ProbeID == TargID ||  // not expecting probe to be self overlapping!
	   PeakScore < 0 || PeakScore > 10000000 ||  FinalScore < 0 || FinalScore > PeakScore ||         // SW scores in the millions surely is an error????
	   Class < (int)eOLCOverlapping || Class > (int)eOLCartefact ||	// overlap classification expected to be in this range	
	   ProbeLen < 1000 || ProbeLen > 100000000 || TargLen < 1000 || TargLen > 100000000 ||    // allowing for PacBio reads of 1Kbp upto 100Mbp!	
	   ProbeStartOfs < 1 || (ProbeStartOfs + 100) > ProbeLen ||  // expecting alignments to be at least 100bp and alignment start/end within the read length
	   ProbeEndOfs < (ProbeStartOfs + 99) || ProbeEndOfs > ProbeLen || 
	   TargStartOfs < 1 || (TargStartOfs + 100) > TargLen ||  
	   TargEndOfs < (TargStartOfs + 99) || TargEndOfs > TargLen ||
	   !(Sense == 'S' || Sense == 's' || Sense == 'A' || Sense == 'a')) 
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error, inconsistencies in values parsed from line %d in file '%s'",NumElsRead,pszPacBioOvlps);
		delete pCSV;
		return((UINT32)eBSFerrParse);
		}
	
	if(bValidateOnly || Class > (int)eOLCartefact)
		{
		if(Class <=  (int)eOLCartefact)
			NumEdges += 1;
		continue;
		}

	UINT32 Score;	// score is dependent on: PeakScore, NumExactBases and ProbeAlignLength
	Score = (UINT32)(((INT64)PeakScore * NumExactBases) / ProbeAlignLength);
	if(Score > 0x07fff)
		Score = 0x07fff;

	// note that AddEdge() is expecting offsets to be 0 based hence substracting 1
	if((NumEdges = m_pAssembGraph->AddEdge(ProbeID,TargID,Score,ProbeStartOfs-1,ProbeEndOfs-1,TargStartOfs-1,TargEndOfs-1,(eOverlapClass)Class,Sense == 's' || Sense == 'S' ? false : true)) > cMaxValidID)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddEdge failed at line %d in file '%s'",NumElsRead,pszPacBioOvlps);
		delete pCSV;
		return(NumEdges);
		}
	}

delete pCSV;
return(NumEdges);
}

int
CPBScaffold::LoadTargetSeqs(char *pszTargFile)		// load pre-generated suffix array from pszTargFile
{
int Rslt;
INT64 FileSize;
char szTargSpecies[120];

// check if pszTargFile exists, is of a minimum length, and can be read
#ifdef _WIN32
struct _stat64 st;
if(!_stat64(pszTargFile,&st))
#else
struct stat64 st;
if(!stat64(pszTargFile,&st))
#endif
	FileSize = (INT64)st.st_size;
else
	FileSize = 0;

if(FileSize == 0)		// 0 if file not readable or if 0 length
	{
	gDiagnostics.DiagOut(eDLWarn,gszProcName,"LoadTargetSeqs: Unable to access file '%s', does file exist and not 0 length, or is it readable",pszTargFile);
	return(eBSFerrFileAccess);
	}

if(FileSize < cMinScaffSeqLen)		// arbitary minimum file size...
	{
	gDiagnostics.DiagOut(eDLWarn,gszProcName,"LoadTargetSeqs: file exists but is only %lld long",pszTargFile,FileSize);
	return(eBSFerrFileAccess);
	}

// first try loading as previously generated suffix index over targeted sequences
if((Rslt = m_pSfxArray->Open(pszTargFile)) < eBSFSuccess)
	{
	while(m_pSfxArray->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,m_pSfxArray->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadTargetSeqs: Unable to load input bioseq suffix array file '%s'",pszTargFile);
	delete m_pSfxArray;
	m_pSfxArray = NULL;
	return(Rslt);
	}

// report to user some sfx array metadata as conformation the targeted assembly is correct
strcpy(szTargSpecies,m_pSfxArray->GetDatasetName());
tsSfxHeaderV3 SfxHeader;
m_pSfxArray->GetSfxHeader(&SfxHeader);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"LoadTargetSeqs:  Name: '%s' Descr: '%s' Title: '%s' Version: %d",
					 szTargSpecies,SfxHeader.szDescription,SfxHeader.szTitle,SfxHeader.Version);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Loading target scaffolding sequences suffix array...");
if((Rslt=m_pSfxArray->SetTargBlock(1))<0)
	{
	while(m_pSfxArray->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,m_pSfxArray->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: unable to load target scaffolding sequences suffix array");
	return(Rslt);
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"target scaffolding sequences suffix array loaded");

return(eBSFSuccess);
}

// ProcessBioseqFile
// Process input biosequence file into suffix file
int
CPBScaffold::ProcessBioseqFile(int MinSeqLen,		// only accept for indexing sequences of at least this length
				 char *pszFile)						// file containing sequences to be parsed and indexed
{
CBioSeqFile BioSeqFile;
etSeqBase *pSeqBuff = NULL;
UINT32 AllocLen = 0;
char szSource[cBSFSourceSize];
char szDescription[cBSFDescriptionSize];
UINT32 SeqLen;
int Rslt;
UINT32 NumSeqsUnderlength;
tBSFEntryID CurEntryID;

if((Rslt=BioSeqFile.Open(pszFile,cBSFTypeSeq,false))!=eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessBioseqFile - error %d %s",Rslt,BioSeqFile.GetErrMsg());
	return(Rslt);
	}

CurEntryID = 0;
NumSeqsUnderlength = 0;
while((Rslt = CurEntryID = BioSeqFile.Next(CurEntryID)) > eBSFSuccess)
	{
	BioSeqFile.GetNameDescription(CurEntryID,cBSFSourceSize-1,(char *)&szSource,
											cBSFDescriptionSize-1,(char *)&szDescription);
	SeqLen = BioSeqFile.GetDataLen(CurEntryID);
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing %s|%s",szSource,szDescription);

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

	if((Rslt=m_pSfxArray->AddEntry(szSource,pSeqBuff,SeqLen)) < eBSFSuccess)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessBioseqFile - error %d %s",Rslt,m_pSfxArray->GetErrMsg());
		break;
		}
	}
if(Rslt == eBSFerrEntry)
	Rslt = eBSFSuccess;
if(pSeqBuff != NULL)
	delete pSeqBuff;
BioSeqFile.Close();
if(NumSeqsUnderlength > 0)
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"ProcessBioseqFile - %u sequences not accepted for indexing as length under %dbp ",NumSeqsUnderlength,MinSeqLen);
return(Rslt);
}

// ProcessFastaFile
// Parse input fasta format file into a biosequence suffix array file
int
CPBScaffold::ProcessFastaFile(int MinSeqLen,		// only accept for indexing sequences of at least this length
				char *pszFile)						// file containing sequences
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
			if(BuffOfs < (size_t)MinSeqLen)
				NumSeqsUnderlength += 1;
			else
				if((Rslt=m_pSfxArray->AddEntry(szName,pSeqBuff,(UINT32)BuffOfs)) < eBSFSuccess)
					{
					gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessFastaFile - error %d %s",Rslt,m_pSfxArray->GetErrMsg());
					break;
					}
				else
					{
					NumSeqsAccepted += 1;
					TotAcceptedLen += BuffOfs;
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
	if(BuffOfs < (size_t)MinSeqLen)
		{
		NumSeqsUnderlength += 1;
		Rslt = eBSFSuccess;
		}
	else
		{
		if((Rslt=m_pSfxArray->AddEntry(szName,pSeqBuff,(UINT32)BuffOfs)) < eBSFSuccess)
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessFastaFile - error %d %s",Rslt,m_pSfxArray->GetErrMsg());
		else
			{
			Rslt = eBSFSuccess;
			NumSeqsAccepted += 1;
			TotAcceptedLen += BuffOfs;
			}
		}
	}
if(pSeqBuff != NULL)
	free(pSeqBuff);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"ProcessFastaFile - %d parsed, %d accepted, %dbp mean length, %d sequences not accepted for indexing as length under %dbp ",
					SeqID,NumSeqsAccepted,NumSeqsAccepted == 0 ? 0 :(int)(TotAcceptedLen/NumSeqsAccepted),NumSeqsUnderlength,MinSeqLen);
return(Rslt);
}


int 
CPBScaffold::LoadTargetSeqs(int MinSeqLen,int NumTargFiles,char **pszTargFiles)		// parse, and index sequences in these files into in memory suffix array; file expected to contain either fasta or fastq sequences
{
int Rslt;
int Idx;
int NumGlobs;
INT64 SumFileSizes;
UINT32 DupEntries[cMaxDupEntries];
int NumDupEntries;
char szDupEntry[100];


m_pSfxArray->SetMaxQSortThreads(m_NumThreads);

CSimpleGlob glob(SG_GLOB_FULLSORT);

	// determine crude estimate of total genome size from the sum of file sizes
for(Idx = 0; Idx < NumTargFiles; Idx++)
	if (glob.Add(pszTargFiles[Idx]) < SG_SUCCESS)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to glob '%s",pszTargFiles[Idx]);
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
	Reset(false);
	return(eBSFerrOpnFile);
	}


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

Rslt = eBSFSuccess;

for (int n = 0; Rslt >= eBSFSuccess &&  n < glob.FileCount(); ++n)
	{
		// try opening as a fasta file, if that fails then try as a bioseq
	Rslt = ProcessFastaFile(MinSeqLen,glob.File(n));
	if(Rslt == eBSFerrNotFasta)
		Rslt = ProcessBioseqFile(MinSeqLen,glob.File(n));
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
		return(-1);
		}
	}
if(m_pSfxArray->GetNumEntries() < 3)		// need at least 3 sequences if attempting to error correct
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"CreateBioseqSuffixFile: Insufficent ( < 3) sequences accepted for processing");
	m_pSfxArray->Close(false);
	Reset(false);
	return(eBSFerrNoEntries);
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"CreateBioseqSuffixFile: sorting suffix array...");
m_pSfxArray->Finalise();
gDiagnostics.DiagOut(eDLInfo,gszProcName,"CreateBioseqSuffixFile: sorting completed");
return(Rslt);

}

int
CPBScaffold::Process(etPBPMode PMode,	// processing mode
		UINT32 DeltaCoreOfs,		// offset by this many bp the core windows of coreSeqLen along the probe sequence when checking for overlaps
		UINT32 MaxSeedCoreDepth,    // only further process a seed core if there are no more than this number of matching cores in all targeted sequences
		UINT32 MinSeedCoreLen,		// use seed cores of this length when identifying putative overlapping scaffold sequences
		UINT32 MinNumSeedCores,		// require at least this many seed cores between overlapping scaffold sequences
		int SWMatchScore,			// score for matching bases (0..50)
		int SWMismatchPenalty,		// mismatch penalty (-50..0)
		int SWGapOpenPenalty,		// gap opening penalty (-50..0)
		int SWGapExtnPenalty,		// gap extension penalty (-50..0)
		int SWProgExtnPenaltyLen,	// only apply gap extension score if gap at least this length (0..63) - use if aligning PacBio
		UINT32 MinScaffSeqLen,		// individual target scaffold sequences must be of at least this length (defaults to 5Kbp)
		UINT32 MinScaffOverlap,		// pairs of targeted scaffold sequences must overlap by at least this many bp to be considered for merging into a longer scaffold sequence (defaults to 5Kbp) 
		char *pszPacBioOvlps,		// pregenerated PacBio sequence overlap loci details
		char *pszPacBioSfxFile,		// suffix indexed PacBio sequences
		int NumPacBioFiles,			// number of input pacbio file specs
		char *pszPacBioFiles[],		// input pacbio files
		char *pszOutFile,			// where to write merged scaffolded sequences
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

m_MinScaffSeqLen = MinScaffSeqLen;	
m_MinScaffOverlap = MinScaffOverlap; 

if(PMode == ePBPMScaffold)
	{
	if((UINT32)(Rslt = (int)LoadPacBioOvlps(pszPacBioOvlps, true)) > cMaxValidID)
		{
		Reset(false);
		return(Rslt);
		}
	if(Rslt = 0)
		{
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"LoadPacBioOvlps: Nothing to do, no overlaps to process");
		Reset(false);
		return(Rslt);
		}
	}

if(pszPacBioSfxFile != NULL && pszPacBioSfxFile[0] != '\0')
	{
	strncpy(m_szPacBioSfxFile,pszPacBioSfxFile,sizeof(m_szPacBioSfxFile));
	m_szPacBioSfxFile[sizeof(m_szPacBioSfxFile)-1] = '\0';	
	}
else
	m_szPacBioSfxFile[0] = '\0';	
m_NumPacBioFiles = NumPacBioFiles;
memset(m_szPacBioFiles,0,sizeof(m_szPacBioFiles));
if(NumPacBioFiles > 0)
	{
	for(Idx = 0; Idx < NumPacBioFiles; Idx++)
		strcpy(m_szPacBioFiles[Idx],pszPacBioFiles[Idx]);
	}

strncpy(m_szOutFile,pszOutFile,sizeof(m_szOutFile));
m_szOutFile[sizeof(m_szOutFile)-1] = '\0';	

strcpy(m_szOutScaffFile,pszOutFile);
if(PMode == ePBPMScaffold)
	strcat(m_szOutScaffFile,"scaffolded.fasta");
else
	strcat(m_szOutScaffFile,"scaffold.csv");

m_NumThreads = NumThreads;	
if(m_pSfxArray != NULL)
	delete m_pSfxArray;
if((m_pSfxArray = new CSfxArrayV3) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadTargetSeqs: Unable to instantiate instance of CSfxArrayV3");
	return(eBSFerrObj);
	}
m_pSfxArray->Reset(false);

if(m_pAssembGraph != NULL)
	{
	delete m_pAssembGraph;
	m_pAssembGraph = NULL;
	}

if(PMode == ePBPMScaffold)
	{
	if((m_pAssembGraph = new CAssembGraph) == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadTargetSeqs: Unable to instantiate instance of CAssembGraph");
		Reset(false);
		return(eBSFerrObj);
		}
	if((Rslt=m_pAssembGraph->Init(min(4,NumThreads)))!=eBSFSuccess)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadTargetSeqs: Unable to initialise CAssembGraph");
		Reset(false);
		return(eBSFerrObj);
		} 
	}

if(pszPacBioSfxFile != NULL && pszPacBioSfxFile[0] != '\0')
	{
	if((Rslt = LoadTargetSeqs(pszPacBioSfxFile)) < eBSFSuccess)
		{
		Reset(false);
		return(Rslt);
		}

	// get number of sequences loaded as targets to be used for scaffolding
	NumTargSeqs = m_pSfxArray->GetNumEntries();
	if(NumTargSeqs < 1)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"No targeted scaffold sequences in file '%s'",pszPacBioSfxFile);
		Reset(false);
		return(eBSFerrNoEntries);
		}
	}
else
	{
	if((Rslt = LoadTargetSeqs(MinScaffSeqLen,NumPacBioFiles,pszPacBioFiles)) < eBSFSuccess)
		{
		Reset(false);
		return(Rslt);
		}
	// get number of sequences loaded as targets to be used for scaffolding
	NumTargSeqs = m_pSfxArray->GetNumEntries();
	if(NumTargSeqs < 1)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"No targeted scaffold sequences in file(s)");
		Reset(false);
		return(eBSFerrNoEntries);
		}
	}

if(PMode == ePBPMOverlaps)
	{
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
	}

if((m_pPBScaffNodes = new tsPBScaffNode [NumTargSeqs + 1]) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate for %d scaffold nodes",NumTargSeqs);
	Reset(false);
	return(eBSFerrMem);
	}
memset(m_pPBScaffNodes,0,sizeof(tsPBScaffNode) * (NumTargSeqs+1));
m_AllocdPBScaffNodes = NumTargSeqs;
if((m_pMapEntryID2NodeIDs = new UINT32 [NumTargSeqs + 1]) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate for %d scaffold nodes",NumTargSeqs);
	Reset(false);
	return(eBSFerrMem);
	}
memset(m_pMapEntryID2NodeIDs,0,sizeof(UINT32) * (NumTargSeqs+1));
m_NumPBScaffNodes = 0;

// initialise scaffold nodes
if(PMode == ePBPMScaffold)
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Initialising %d scaffold nodes and graph vertices",NumTargSeqs);
else
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Initialising %d scaffold nodes",NumTargSeqs);	
MaxSeqLen = 0;
pCurPBScaffNode = m_pPBScaffNodes;
for(CurNodeID = 1; CurNodeID <= NumTargSeqs; CurNodeID++,pCurPBScaffNode++)
	{
	pCurPBScaffNode->SeqLen = m_pSfxArray->GetSeqLen(CurNodeID);
	pCurPBScaffNode->EntryID = CurNodeID;
	pCurPBScaffNode->flgUnderlength = pCurPBScaffNode->SeqLen < m_MinScaffSeqLen ? 1 : 0;
	if(MaxSeqLen == 0 || pCurPBScaffNode->SeqLen > (UINT32)MaxSeqLen)
		MaxSeqLen = pCurPBScaffNode->SeqLen;

	if(PMode == ePBPMScaffold)
		{
		if((pCurPBScaffNode->VertexID = m_pAssembGraph->AddVertex(pCurPBScaffNode->SeqLen,CurNodeID)) > cMaxValidID)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddVertex failed");
			Reset(false);
			return((int)pCurPBScaffNode->VertexID);
			}
		}
	}
m_NumPBScaffNodes = NumTargSeqs;
if(PMode == ePBPMScaffold)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Finalising %d graph vertices",NumTargSeqs);
	m_pAssembGraph->FinaliseVertices();
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Finalised graph vertices");
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Sorting %d scaffold nodes",NumTargSeqs);

if(m_NumPBScaffNodes > 1)
	{
	// sort scaffold nodes by sequence length descending
	m_mtqsort.SetMaxThreads(NumThreads);
	m_mtqsort.qsort(m_pPBScaffNodes,m_NumPBScaffNodes,sizeof(tsPBScaffNode),SortLenDescending);
	}

pCurPBScaffNode = m_pPBScaffNodes;
for(CurNodeID = 1; CurNodeID <= NumTargSeqs; CurNodeID++,pCurPBScaffNode++)
	{
	pCurPBScaffNode->NodeID = CurNodeID;
	m_pMapEntryID2NodeIDs[pCurPBScaffNode->EntryID-1] = CurNodeID;
	}

if(PMode == ePBPMOverlaps)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Identifying scaffold overlaps ...");
	Rslt = IdentifyScaffoldOverlaps(MaxSeqLen,NumThreads);
	Reset(false);
	return(Rslt);
	}	

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Loading previously generated PacBio sequence overlaps from file '%s' ...",pszPacBioOvlps);
if((UINT32)(Rslt = LoadPacBioOvlps(pszPacBioOvlps)) > cMaxValidID)
	return(Rslt);
if(Rslt == 0)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"No scaffold overlaps");
	return(Rslt);										
	}

if((UINT32)(Rslt = m_pAssembGraph->FinaliseEdges())  > cMaxValidID)
	{
	Reset(false);
	return(Rslt);
	}
if(Rslt == 0)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"FinaliseEdges: No scaffold overlaps");
	return(Rslt);										
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Accepted %u overlaps from file '%s' ...",(UINT32)Rslt,pszPacBioOvlps);


UINT32 NumDiscComponents = 0;
if((NumDiscComponents = m_pAssembGraph->IdentifyDiscComponents()) > cMaxValidID)
	{
	Reset(false);
	return((int)NumDiscComponents);
	}
m_pAssembGraph->WriteScaffoldSeqs(pszOutFile);

Reset(false);
return(Rslt);
}

#ifdef _WIN32
unsigned __stdcall IdentScaffOverlapsThread(void * pThreadPars)
#else
void *IdentScaffOverlapsThread(void * pThreadPars)
#endif
{
int Rslt;
tsThreadPBScaffold *pPars = (tsThreadPBScaffold *)pThreadPars;			// makes it easier not having to deal with casts!
CPBScaffold *pPBScaffold = (CPBScaffold *)pPars->pThis;

Rslt = pPBScaffold->ThreadIdentScaffOverlaps(pPars);
pPars->Rslt = Rslt;
#ifdef _WIN32
_endthreadex(0);
return(eBSFSuccess);
#else
pthread_exit(NULL);
#endif
}

int
CPBScaffold::IdentifyScaffoldOverlaps(int MaxSeqLen,		// max length sequence to be overlapped
									int NumOvlpThreads)		// identify all scaffold overlaps using this many threads
{
tsThreadPBScaffold *pThreadPutOvlps;
int ThreadIdx;
tsThreadPBScaffold *pThreadPar;

pThreadPutOvlps = new tsThreadPBScaffold [NumOvlpThreads];

pThreadPar = pThreadPutOvlps;
for(ThreadIdx = 0; ThreadIdx < NumOvlpThreads; ThreadIdx++,pThreadPar++)
	{
	memset(pThreadPar,0,sizeof(tsThreadPBScaffold));
	pThreadPar->AllocdCoreHits = cAllocdNumCoreHits;
	pThreadPar->AllocdCoreHitsSize = cAllocdNumCoreHits * sizeof(tsCoreHit);
#ifdef _WIN32
	pThreadPar->pCoreHits = (tsCoreHit *)malloc(pThreadPar->AllocdCoreHitsSize);	
	if(pThreadPar->pCoreHits == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: Core hits memory allocation of %llu bytes - %s",pThreadPar->AllocdCoreHitsSize,strerror(errno));
		break;
		}
#else
	if((pThreadPar->pCoreHits = (tsCoreHit *)mmap(NULL,pThreadPar->AllocdCoreHitsSize, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0)) == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: Core hits memory allocation of %llu bytes - %s",pThreadPar->AllocdCoreHitsSize,strerror(errno));
		break;
		}
#endif

	pThreadPar->AllocdProbeSeqSize = max(cAllocdQuerySeqLen,MaxSeqLen);
	if((pThreadPar->pProbeSeq = new etSeqBase [pThreadPar->AllocdProbeSeqSize + 10])==NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: Core hits probe sequence memory allocation of %d bytes - %s",pThreadPar->AllocdProbeSeqSize,strerror(errno));
		break;
		}

	pThreadPar->AllocdTargSeqSize = pThreadPar->AllocdProbeSeqSize;
	if((pThreadPar->pTargSeq = new etSeqBase [pThreadPar->AllocdTargSeqSize + 10])==NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: Core hits target sequence memory allocation of %d bytes - %s",pThreadPar->AllocdTargSeqSize,strerror(errno));
		break;
		}

	if((pThreadPar->pmtqsort = new CMTqsort) == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: Core hits instantiation of CMTqsort failed");
		break;
		}
	pThreadPar->pmtqsort->SetMaxThreads(4);
	pThreadPar->bRevCpl = false;
	pThreadPar->MaxSeedCoreDepth = m_MaxSeedCoreDepth;
	pThreadPar->DeltaCoreOfs = m_DeltaCoreOfs;
	pThreadPar->CoreSeqLen = m_MinSeedCoreLen;
	pThreadPar->MinNumCores = m_MinNumSeedCores;
	pThreadPar->MaxAcceptHitsPerSeedCore = cDfltMaxAcceptHitsPerSeedCore;
	pThreadPar->MinScaffSeqLen = m_MinScaffSeqLen;
	pThreadPar->MinOverlapLen = m_MinScaffOverlap;
	}

if(ThreadIdx != NumOvlpThreads)	// any errors whilst allocating memory for core hits?
	{
	do {
		if(pThreadPar->pmtqsort != NULL)
			delete pThreadPar->pmtqsort;

		if(pThreadPar->pAntisenseKmers != NULL)
			delete pThreadPar->pAntisenseKmers;

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

#ifdef _WIN32
m_hScaffDist = open(m_szOutScaffFile,( O_WRONLY | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC),(_S_IREAD | _S_IWRITE));
#else
if((m_hScaffDist = open(m_szOutScaffFile,O_WRONLY | O_CREAT,S_IREAD | S_IWRITE))!=-1)
	if(ftruncate(m_hScaffDist,0)!=0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate %s - %s",m_szOutScaffFile,strerror(errno));
		Reset(false);
		return(eBSFerrCreateFile);
		}
#endif
if(m_hScaffDist < 0)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate output file '%s'",m_szOutScaffFile);
	Reset(false);
	return(eBSFerrCreateFile);
	}

#ifdef USETHISCODE
if(m_hScaffDist != -1)
	{
	m_ScaffLineBuffIdx=sprintf(m_szScaffLineBuff,"\"Class\",\"ProbeID\",\"ProbDescr\",\"TargID\",\"TargDescr\",\"SeedHits\",\"Sense\",\"ProbeLen\",\"TargLen\",\"ProbeAlignLength\",\"TargAlignLength\",\"PeakScore\",\"FinalScore\",\"NumAlignedBases\",\"NumExactBases\",\"NumProbeInserts\",\"NumProbeInsertBases\",\"NumTargInserts\",\"NumTargInsertBases\",\"ProbeStartOfs\",\"TargStartOfs\",\"ProbeEndOfs\",\"TargEndOfs\",\"ProbeOfs5\",\"TargOfs5\",\"ProbeOfs3\",\"TargOfs3\"");
#ifdef _PEAKSCOREACCEPT_
	m_ScaffLineBuffIdx+=sprintf(&m_szScaffLineBuff[m_ScaffLineBuffIdx],",\"PSAlignLength\",\"PSPeakScore\",\"PSNumAlignedBases\",\"PSNumExactBases\",\"PSNumProbeInserts\",\"PSNumProbeInsertBases\",\"PSNumTargInserts\",\"PSNumTargInsertBases\",\"PSProbeStartOfs\",\"PSTargStartOfs\",\"ProbeEndOfs\",\"TargEndOfs\",\"PSProbeOfs5\",\"PSTargOfs5\",\"PSProbeOfs3\",\"PSTargOfs3\"");
#endif
#ifdef _EXACTMATCHLENDIST_
	for(int ExactLenIdx = 1; ExactLenIdx <= cExactMatchLenDist; ExactLenIdx++)
		m_ScaffLineBuffIdx+=sprintf(&m_szScaffLineBuff[m_ScaffLineBuffIdx],",\"ExactLen:%d\"",ExactLenIdx);
#endif
	CUtility::SafeWrite(m_hScaffDist,m_szScaffLineBuff,m_ScaffLineBuffIdx);
	m_ScaffLineBuffIdx = 0;
	}
#endif

pThreadPar = pThreadPutOvlps;
for (ThreadIdx = 1; ThreadIdx <= NumOvlpThreads; ThreadIdx++, pThreadPar++)
	{
	pThreadPar->ThreadIdx = ThreadIdx;
	pThreadPar->pThis = this;
#ifdef _WIN32
	pThreadPar->threadHandle = (HANDLE)_beginthreadex(NULL, 0x0fffff, IdentScaffOverlapsThread, pThreadPar, 0, &pThreadPar->threadID);
#else
	pThreadPar->threadRslt = pthread_create(&pThreadPar->threadID, NULL, IdentScaffOverlapsThread, pThreadPar);
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
	}

delete pThreadPutOvlps;
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Completed: %u processed, SW aligned: %u, Overlapping: %u, Overlapped: %u, Contained: %u, Artefact: %u",
							m_NumOverlapProcessed,m_ProvSWchecked,m_ProvOverlapping,m_ProvOverlapped,m_ProvContained,m_ProvArtefact);

if(m_hScaffDist != -1)
	{
	if(m_ScaffLineBuffIdx > 0)
		{
		CUtility::SafeWrite(m_hScaffDist,m_szScaffLineBuff,m_ScaffLineBuffIdx);
		m_ScaffLineBuffIdx = 0;
		}
#ifdef _WIN32
	_commit(m_hScaffDist);
#else
	fsync(m_hScaffDist);
#endif
	close(m_hScaffDist);
	m_hScaffDist = -1;
	}

return(0);
}

// expectations -
// core hits in pThreadPar have bee sorted ProbeID.TargID.TargOfs.ProbeOfs ascending
int 
CPBScaffold::ClusterSpatial(tsThreadPBScaffold *pThreadPar, 
			   UINT32 ProbeLen,					// probe from which cores were used to align against targets was this length
			   tsCoreHitsClusters *pClusters,	// returned clusters
			   UINT32 MinClusterHits,				// clusters must contain at least this number of consistency checked hits
			   UINT32 MinClusterLen)				// clusters must be at least this length
{
UINT32 MaxClusterLen;
UINT32 NumClustHits;
UINT32 SumClustHitLens;

tsCoreHit *pClustCoreHit;


UINT32 CurTargNodeID;
UINT32 CurProbeNodeID;
tsCoreHit *pCurCoreHit;

UINT32 ClustProbeOfs;
UINT32 IntraClustProbeOfs;
UINT32 ClustTargOfs;
UINT32 ClustLen;

double ClustScore;
double LowestClustScore;
tsCluster *pCluster;
int ClusterIdx;
int SelClustIdx;

UINT32 DeltaProbeOfs;
UINT32 DeltaTargOfs;

if( pClusters != NULL)
	memset(pClusters,0,sizeof(tsCoreHitsClusters));

if(pThreadPar == NULL || pThreadPar->pCoreHits == NULL || MinClusterHits < 2 || pThreadPar->NumCoreHits < MinClusterHits || pClusters == NULL || MinClusterLen < 50 || ProbeLen < MinClusterLen)
	return(0);

// sort the hits by Probe.Target.TargOfs.ProbeOfs ascending so can process per target
pThreadPar->pmtqsort->qsort(pThreadPar->pCoreHits,pThreadPar->NumCoreHits,sizeof(tsCoreHit),SortCoreHitsByTargProbeOfs);

MaxClusterLen = (ProbeLen * 110) / 100;			// allowing 10% for PacBio insertions
pCurCoreHit = pThreadPar->pCoreHits;
CurProbeNodeID = pCurCoreHit->ProbeNodeID;
pClusters->ProbeID = CurProbeNodeID; 
pClusters->ProbeSeqLen = m_pPBScaffNodes[CurProbeNodeID-1].SeqLen;

while(pCurCoreHit->ProbeNodeID != 0)
	{
	CurTargNodeID = pCurCoreHit->TargNodeID;

	while(pCurCoreHit->ProbeNodeID != 0 && pCurCoreHit->ProbeNodeID == CurProbeNodeID && pCurCoreHit->TargNodeID == CurTargNodeID)
		{
		if(pThreadPar->bRevCpl != pCurCoreHit->flgRevCpl ? true : false)
			{
			pCurCoreHit += 1;
			continue;
			}
		ClustProbeOfs = pCurCoreHit->ProbeOfs;	// note starting offsets as these are used intra cluster to determine consistency
		ClustTargOfs = pCurCoreHit->TargOfs;
		IntraClustProbeOfs = ClustProbeOfs;
		NumClustHits = 1;
		ClustLen = SumClustHitLens = pCurCoreHit->HitLen;
		pClustCoreHit = pCurCoreHit + 1;
		while(pClustCoreHit->ProbeNodeID == CurProbeNodeID && pClustCoreHit->TargNodeID == CurTargNodeID)
			{
			if(pThreadPar->bRevCpl != pClustCoreHit->flgRevCpl ? true : false)
				{
				pClustCoreHit += 1;
				continue;
				}

			if(pClustCoreHit->TargOfs >= ClustTargOfs + MaxClusterLen)	// target hits must be within the probe length + 10% to allow for the PacBio insertions
				break;
			if(pClustCoreHit->ProbeOfs <= IntraClustProbeOfs) // looking for next ascending probe core which hits onto this target
				{
				pClustCoreHit += 1;
				continue;
				}
			DeltaProbeOfs = pClustCoreHit->ProbeOfs - ClustProbeOfs;
			DeltaTargOfs = pClustCoreHit->TargOfs - ClustTargOfs;
			if(DeltaProbeOfs > max(25,(DeltaTargOfs * 110)/100) || DeltaTargOfs > max(25,(DeltaProbeOfs * 110)/100))
				{
				pClustCoreHit += 1;
				continue;		
				}
			// accepting this core hit as being consistent between probe and target
			NumClustHits += 1;
			SumClustHitLens += pClustCoreHit->HitLen;
			IntraClustProbeOfs = pClustCoreHit->ProbeOfs;
			ClustLen = pClustCoreHit->HitLen + IntraClustProbeOfs - ClustProbeOfs; 
			pClustCoreHit += 1;
			}
		pCurCoreHit += 1;
		if(ClustLen < MinClusterLen || ClustLen > MaxClusterLen || NumClustHits < MinClusterHits)
			continue;

		// if a cluster then is this cluster in the top cMaxClusters?
		// clusters are scored as a function of the SumClustLens and the ratio of the number of cores in cluster to the cluster length
        // e.g (SumClustLens * NumClustHits) / ClustLen
		ClustScore = ((double)SumClustHitLens * NumClustHits) / ClustLen;
		pCluster = &pClusters->Clusters[0];
		SelClustIdx = -1;
		if(pClusters->NumClusters)
			{
			for(ClusterIdx = 0; ClusterIdx < pClusters->NumClusters; ClusterIdx++,pCluster++)
				{
				if(pCluster->TargNodeID == CurTargNodeID)
					{
					if(pCluster->ClustScore < ClustScore)
						SelClustIdx = ClusterIdx;
					else
						SelClustIdx = -2;
					break;
					}
				}
							
			if(SelClustIdx == -1)
				{
				if(pClusters->NumClusters < cMaxClusters)
					{
					pClusters->NumClusters += 1;
					SelClustIdx = ClusterIdx;
					}
				else
					{
					LowestClustScore = ClustScore;
					pCluster = &pClusters->Clusters[0];
					for(ClusterIdx = 0; ClusterIdx < pClusters->NumClusters; ClusterIdx++,pCluster++)
						if(pCluster->ClustScore < LowestClustScore)
							{
							LowestClustScore = pCluster->ClustScore;
							SelClustIdx = ClusterIdx;
							}
					}
				}
			}
		else
			{
			SelClustIdx = 0;
			pClusters->NumClusters = 1;
			}
		if(SelClustIdx >= 0)
			{
			pCluster = pCluster = &pClusters->Clusters[SelClustIdx];
			pCluster->TargNodeID = CurTargNodeID;	
			pCluster->TargSeqLen = m_pPBScaffNodes[CurTargNodeID-1].SeqLen;
			pCluster->ClustProbeOfs = ClustProbeOfs;		
			pCluster->ClustTargOfs = ClustTargOfs;
			pCluster->ClustLen = ClustLen;
			pCluster->NumClustHits = NumClustHits;
			pCluster->SumClustHitLens = SumClustHitLens;
			pCluster->ClustScore = ClustScore;
			}
		}
	}

if(pClusters->NumClusters && m_hScaffDist != -1)
	{
	AcquireSerialise();
	if(m_ScaffLineBuffIdx > sizeof(m_szScaffLineBuff) - 1000)
		{
		CUtility::SafeWrite(m_hScaffDist,m_szScaffLineBuff,m_ScaffLineBuffIdx);
		m_ScaffLineBuffIdx = 0;
		}

	pCluster = pClusters->Clusters; 
	for(ClusterIdx = 0; ClusterIdx < pClusters->NumClusters; ClusterIdx++,pCluster++)
		{
		m_ScaffLineBuffIdx += sprintf(&m_szScaffLineBuff[m_ScaffLineBuffIdx], "\n\"%c\",%1.3f,%d,%d,%d,%d,%d,%d,%d,%d,%d",
							pThreadPar->bRevCpl ? 'A' : 'S',pCluster->ClustScore,pClusters->ProbeID,pClusters->ProbeSeqLen,
							pCluster->TargNodeID,pCluster->TargSeqLen,pCluster->ClustProbeOfs,pCluster->ClustTargOfs,pCluster->ClustLen,pCluster->NumClustHits,pCluster->SumClustHitLens);
		}
	ReleaseSerialise();
	}

return(pClusters->NumClusters);
}


int
CPBScaffold::ThreadIdentScaffOverlaps(tsThreadPBScaffold *pThreadPar)
{
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
UINT32 CurHitLoci;
UINT32 CurSEntryIDHits;
UINT32 CurAEntryIDHits;
UINT32 ProvOverlapping;
UINT32 ProvOverlapped;
UINT32 ProvContained;
UINT32 ProvArtefact;
UINT32 ProvSWchecked;

tsPBScaffNode *pHitScaffNode;
sPBCoreHitCnts *pSummaryCnts;
UINT32 CurNodeID;
int Class;
char szTargSeqName[cMaxDatasetSpeciesChrom];
char szProbeSeqName[cMaxDatasetSpeciesChrom];

tsPBScaffNode *pCurPBScaffNode;
pThreadPar->pSW = NULL;

for(CurNodeID = 1; CurNodeID <= m_NumPBScaffNodes; CurNodeID++)
	{
	pThreadPar->NumCoreHits = 0;
	pCurPBScaffNode = &m_pPBScaffNodes[CurNodeID-1];
	AcquireLock(false);
	// slough if already processed, perhaps by another thread, or if underlength
	if(pCurPBScaffNode->flgCurProc == 1 || pCurPBScaffNode->flgContained == 1|| pCurPBScaffNode->flgUnderlength == 1 || pCurPBScaffNode->SeqLen < ((pThreadPar->MinOverlapLen * 3)/2)) 
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
					(pCoreHit->ProbeOfs >= CurProbeOfs && (pCoreHit->ProbeOfs <= (CurProbeOfs + 50)))) // note that probe seeds are clustered into bins of 50bp
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
					(pCoreHit->TargOfs >= CurTargOfs && (pCoreHit->TargOfs <= (CurTargOfs + 50)))) // note that targ seeds are clustered into bins of 50bp
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

		// iterate and count hits for each TargNodeID
		CurSEntryIDHits = 0;
		CurAEntryIDHits = 0;
		CurTargNodeID = 0;
		CurHitLoci = 0;
		pCoreHit = pThreadPar->pCoreHits;
		for(HitIdx = 0; HitIdx < pThreadPar->NumCoreHits; HitIdx++,pCoreHit++)
			{
			pHitScaffNode = &m_pPBScaffNodes[pCoreHit->ProbeNodeID-1];
			if(CurTargNodeID == 0 || pCoreHit->TargNodeID == CurTargNodeID)
				{
				if(CurTargNodeID == 0)
					{
					CurTargNodeID = pCoreHit->TargNodeID;
					CurSEntryIDHits = 0;
					CurAEntryIDHits = 0;
					}
				if(pCoreHit->flgMulti == 1)
					continue;
				CurHitLoci = pCoreHit->TargOfs;
				if(pCoreHit->flgRevCpl == 0)
					CurSEntryIDHits += 1;
				else
					CurAEntryIDHits += 1;
				continue;
				}
			else   // else cores on a different target 
				{
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
							continue;
						pSummaryCnts = pLowestSummaryHitCnts;
						}
					else
						pSummaryCnts = &pThreadPar->TargCoreHitCnts[pThreadPar->NumTargCoreHitCnts++];
					pSummaryCnts->TargNodeID = CurTargNodeID;
					pSummaryCnts->NumSHits = CurSEntryIDHits;
					pSummaryCnts->NumAHits = CurAEntryIDHits;
					}

				CurHitLoci = pCoreHit->TargOfs;
				CurTargNodeID = pCoreHit->TargNodeID;
				if(pCoreHit->flgRevCpl == 0)
					{
					CurSEntryIDHits = 1;
					CurAEntryIDHits = 0;
					}
				else
					{
					CurSEntryIDHits = 0;
					CurAEntryIDHits = 1;
					}
				}
			}
		if(CurTargNodeID > 0 && (CurSEntryIDHits >= pThreadPar->MinNumCores || CurAEntryIDHits >= pThreadPar->MinNumCores))
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

				if((CurSEntryIDHits + CurAEntryIDHits) > LowestSummaryHitCnts)
					pSummaryCnts = pLowestSummaryHitCnts;
				else
					pSummaryCnts = NULL;
				}
			else
				pSummaryCnts = &pThreadPar->TargCoreHitCnts[pThreadPar->NumTargCoreHitCnts++];
			if(pSummaryCnts != NULL)
				{
				pSummaryCnts->TargNodeID = CurTargNodeID;
				pSummaryCnts->NumSHits = CurSEntryIDHits;
				pSummaryCnts->NumAHits = CurAEntryIDHits;
				}
			}
		}

	// can't process, SW over all would be too resource intensive, all targets which meet the minimum number of core hits requested so choose the top 30 as ranked by the number of core hits
	if(pThreadPar->NumTargCoreHitCnts > 1)
		{
		pThreadPar->pmtqsort->qsort(pThreadPar->TargCoreHitCnts,pThreadPar->NumTargCoreHitCnts,sizeof(sPBCoreHitCnts),SortCoreHitsDescending);
		if(pThreadPar->NumTargCoreHitCnts > 30)
			pThreadPar->NumTargCoreHitCnts = 30;
		}

	if(pThreadPar->NumTargCoreHitCnts > 0)
		{
		UINT32 CurTargCoreHitCnts;
		tsPBScaffNode *pTargNode;
		UINT32 ProbeAlignLength;
		UINT32 TargAlignLength;
		UINT32 TargSeqLen;
		UINT32 LongSAligns = 0;
		UINT32 LongAAligns = 0;


		tsSSWCell *pPeakMatchesCell;
		tsSSWCell PeakMatchesCell;
#ifdef _PEAKSCOREACCEPT
		tsSSWCell PeakScoreCell;
#endif
		if(pThreadPar->pSW == NULL)
			{
			AcquireSerialise();
			pThreadPar->pSW = new CSSW;
			pThreadPar->pSW->SetScores(m_SWMatchScore,m_SWMismatchPenalty,m_SWGapOpenPenalty,m_SWGapExtnPenalty,m_SWProgExtnPenaltyLen,min(63,m_SWProgExtnPenaltyLen+3),cAnchorLen);
			pThreadPar->pSW->PreAllocMaxTargLen(50000);
			ReleaseSerialise();
			}

		if(pThreadPar->NumTargCoreHitCnts >= 2)
			pThreadPar->pSW->StartMultiAlignments(pCurPBScaffNode->SeqLen,pThreadPar->pProbeSeq,pThreadPar->NumTargCoreHitCnts);

		m_pSfxArray->GetIdentName(pCurPBScaffNode->EntryID,sizeof(szProbeSeqName),szProbeSeqName);
		pThreadPar->pSW->SetProbe(pCurPBScaffNode->SeqLen,pThreadPar->pProbeSeq);
		pSummaryCnts = &pThreadPar->TargCoreHitCnts[0];
		for(CurTargCoreHitCnts = 0; CurTargCoreHitCnts < pThreadPar->NumTargCoreHitCnts; CurTargCoreHitCnts++,pSummaryCnts++)
			{
			if(pSummaryCnts->NumSHits < pThreadPar->MinNumCores && pSummaryCnts->NumAHits < pThreadPar->MinNumCores)
				continue;
			bool bTargSense = pSummaryCnts->NumSHits >= pSummaryCnts->NumAHits ? true :  false;
			pTargNode = &m_pPBScaffNodes[pSummaryCnts->TargNodeID-1];
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
#ifdef _PEAKSCOREACCEPT_
			pPeakMatchesCell = pThreadPar->pSW->Align(&PeakScoreCell);
#else
			pPeakMatchesCell = pThreadPar->pSW->Align();
#endif
			ProvSWchecked += 1;
			if(pPeakMatchesCell != NULL && pPeakMatchesCell->NumMatches >= pThreadPar->MinOverlapLen)
				{
				ProvOverlapping += 1;
				PeakMatchesCell = *pPeakMatchesCell;
				ProbeAlignLength = PeakMatchesCell.EndPOfs - PeakMatchesCell.StartPOfs + 1;
				TargAlignLength = PeakMatchesCell.EndTOfs - PeakMatchesCell.StartTOfs + 1;

				// characterise the overlapped target
				// eOLCOverlapping if probe accepted as overlapping, either 5' or 3'
				// eOLCcontaining if both ends of target completely contained within probe
                // eOLCartefact if target is only partially contained
				Class = (int)eOLCOverlapping;		
				if(PeakMatchesCell.StartTOfs < m_OverlapFloat && (TargSeqLen - PeakMatchesCell.EndTOfs) < m_OverlapFloat)
					Class = (int)eOLCcontaining;
				if((PeakMatchesCell.StartTOfs >= m_OverlapFloat &&  PeakMatchesCell.StartPOfs >= m_OverlapFloat) ||
					 ((TargSeqLen - PeakMatchesCell.EndTOfs) >= m_OverlapFloat && (pCurPBScaffNode->SeqLen - PeakMatchesCell.EndPOfs) >= m_OverlapFloat))
					{
					Class = (int)eOLCartefact;
					ProvArtefact += 1;
					}

				if(Class == (int)eOLCcontaining)
					{
					AcquireLock(true);
					pTargNode->flgContained = 1;
					ReleaseLock(true);	
					ProvContained += 1;
					}
				
				if(pThreadPar->NumTargCoreHitCnts >= 2)
					{
					pThreadPar->pSW->TracebacksToAlignOps(PeakMatchesCell.PFirstAnchorStartOfs,PeakMatchesCell.PLastAnchorEndOfs,
															PeakMatchesCell.TFirstAnchorStartOfs,PeakMatchesCell.TLastAnchorEndOfs);
															

					pThreadPar->pSW->AddMultiAlignment(PeakMatchesCell.PFirstAnchorStartOfs,PeakMatchesCell.PLastAnchorEndOfs,
															PeakMatchesCell.TFirstAnchorStartOfs,PeakMatchesCell.TLastAnchorEndOfs,pThreadPar->pTargSeq);
					}

				if(m_hScaffDist != -1)
					{
					m_pSfxArray->GetIdentName(pTargNode->EntryID,sizeof(szTargSeqName)-1,szTargSeqName);
					AcquireSerialise();
					if(m_ScaffLineBuffIdx > (sizeof(m_szScaffLineBuff) - 1000))
						{
						CUtility::SafeWrite(m_hScaffDist,m_szScaffLineBuff,m_ScaffLineBuffIdx);
						m_ScaffLineBuffIdx = 0;
						}

					m_ScaffLineBuffIdx += sprintf(&m_szScaffLineBuff[m_ScaffLineBuffIdx], "\n%d,%d,\"%s\",%d,\"%s\",%d,\"%c\",%1.5d,%1.5d,%1.5d,%1.5d,%1.5d,%1.5d,%1.5d,%1.5d,%1.5d,%1.5d,%1.5d,%1.5d,%1.5d,%1.5d,%1.5d,%1.5d,%1.5d,%1.5d,%1.5d,%1.5d",
							Class,pCurPBScaffNode->EntryID,szProbeSeqName,pTargNode->EntryID,szTargSeqName,
																			bTargSense ? pSummaryCnts->NumSHits : pSummaryCnts->NumAHits, bTargSense ? 'S' : 'A',
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
					ProvOverlapped += 1;
					}
				}
			}
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
	pThreadPar->NumTargCoreHitCnts = 0;
	pThreadPar->NumCoreHits = 0;
	}
AcquireSerialise();
if(m_hScaffDist != -1 && m_ScaffLineBuffIdx > 0)
	{
	CUtility::SafeWrite(m_hScaffDist,m_szScaffLineBuff,m_ScaffLineBuffIdx);
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
CPBScaffold::AddCoreHit(UINT32 ProbeNodeID,			// core hit was from this probe scaffold node 
			   bool bRevCpl,					// true if core sequence was revcpl'd before matching
			   UINT32 ProbeOfs,                 // hit started at this probe offset
			   UINT32 TargNodeID,               // probe core matched onto this target scaffold node
			   UINT32 TargOfs,                  // probe core matched starting at this target loci
			   UINT32 HitLen,					// hit was of this length
               tsThreadPBScaffold *pPars)			// thread specific
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
CPBScaffold::MapEntryID2NodeID(UINT32 EntryID)			// suffix array entry identifier
{
if(EntryID == 0 || EntryID > m_NumPBScaffNodes || m_pMapEntryID2NodeIDs == NULL)
	return(0);
return(m_pMapEntryID2NodeIDs[EntryID-1]);
}



int
CPBScaffold::IdentifyCoreHits(UINT32 ProbeNodeID,	// identify all overlaps of this probe sequence ProbeNodeID onto target sequences
				tsThreadPBScaffold *pPars)		// thread specific
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

for(ProbeOfs = 0; ProbeOfs < (pProbeNode->SeqLen - (pPars->CoreSeqLen + 120)); ProbeOfs+=pPars->DeltaCoreOfs,pCoreSeq+=pPars->DeltaCoreOfs)
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

    while((NextHitIdx = m_pSfxArray->IteratePacBio(pCoreSeq,pProbeNode->SeqLen - ProbeOfs,pPars->CoreSeqLen,pProbeNode->EntryID,pPars->MinScaffSeqLen,PrevHitIdx,&HitEntryID,&HitLoci)) > 0)
		{
		PrevHitIdx = NextHitIdx;
		pTargNode = &m_pPBScaffNodes[MapEntryID2NodeID(HitEntryID)-1];
		HitSeqLen = pTargNode->SeqLen; 
		AcquireLock(false);

		if(pTargNode->flgContained == 1 || 
				(pPars->bSelfHits ? pTargNode->NodeID != pProbeNode->NodeID : pTargNode->NodeID == pProbeNode->NodeID) || 
				 pTargNode->flgUnderlength == 1 ||	// not interested in selfhits or underlength targets
							HitSeqLen < (UINT32)pPars->MinScaffSeqLen ||		// not interested if target sequence length less than min sequence length to be processed
							HitSeqLen > pProbeNode->SeqLen)					// not interested if target sequence length longer than probe sequence length as target already prcessed for overlays onto shorter incl this current probe
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
CPBScaffold::SortLenDescending(const void *arg1, const void *arg2)
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
CPBScaffold::SortCoreHitsByTargProbeOfs(const void *arg1, const void *arg2)
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
CPBScaffold::SortCoreHitsByProbeTargOfs(const void *arg1, const void *arg2)
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
CPBScaffold::SortCoreHitsDescending(const void *arg1, const void *arg2)
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





