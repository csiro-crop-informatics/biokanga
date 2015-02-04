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

// Loads NGS reads and processes same for quality scores, base compositions, kmers, etc., an extension of the functionality provided by fastqc 
// 3.1.2 - increased ngsqc default K-mer max from 5 to 6 so as to explore the hexamer primer artefact conjecture

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

#include "./biokanga.h"
#include "ReadStats.h"

#include "../libbiokanga/bgzf.h"
#include "../libbiokanga/SAMfile.h"
#include "../libBKPLPlot/BKPLPlot.h"

int
ProcessReadsetDist(etRSDMode PMode,				// processing mode; default is for independent readsets
				bool bStrand,					// true if read strand specific distributions
				int Trim5,						// trim this number of bases from 5' end of reads when loading the reads
				int Trim3,						// trim this number of bases from 3' end of reads when loading the reads
				int MaxKMerLen,					// processing is for upto this KMer length inclusive
				int KMerCCC,					// concordance correlation coefficient measure KMer length
				int MaxContamSubRate,			// max allowed contamimant substitution rate (bases per 25bp of contaminant overlap, 1st 15bp of overlap no subs allowed)
				int MinContamLen,				// accept contaminant overlaps if overlap at least this many bases 
				int ReqMaxDupSeeds,				// requested to sample for this many duplicate seeds
				int MinPhredScore,				// only accept reads for duplicate and KMer processing if mean Phred score is at least this threshold 
				int NumThreads,					// number of worker threads to use
				bool bAffinity,					// thread to core affinity
				int NumPE1InputFiles,			// number of PE1 input files
				char *pszInPE1files[],			// input PE1 5' read files
				int NumPE2InputFiles,			// number of PE2 input files
				char *pszInPE2files[],		    // input PE2 3' read files
				char *pszContaminantFile,		// contaminants fasta file
				char *pszOutDistFile,			// where to write distributions CSV file
				char *pszOutHTMLFile);			//  where to write distributions HTML5 file



#ifdef _WIN32
int ReadsetDists(int argc, char* argv[])
{
	// determine my process name
	_splitpath(argv[0], NULL, NULL, gszProcName, NULL);
#else
int
ReadsetDists(int argc, char** argv)
{
	// determine my process name
CUtility::splitpath((char *)argv[0], NULL, gszProcName);
#endif
int iScreenLogLevel;		// level of screen diagnostics
int iFileLogLevel;			// level of file diagnostics
char szLogFile[_MAX_PATH];	// write diagnostics to this file
int Rslt = 0;   			// function result code >= 0 represents success, < 0 on failure
int Idx;

int PMode;					// processing mode
bool bStrand;				// true if read strand specific distributions
int Trim5;					// trim this number of bases from 5' end of reads when loading the reads
int Trim3;					// trim this number of bases from 3' end of reads when loading the reads
int ReqMaxDupSeeds;			// requested to sample for this many duplicate seeds
int MaxKMerLen;				// processing is for upto this KMer length inclusive
int KMerCCC;				// concordance correlation coefficient measure KMer length
int MaxContamSubRate;		// max allowed contamimant substitution rate (bases per 25bp of contaminant overlap, 1st 15bp of overlap no subs allowed)
int MinContamLen;			// accept contaminant overlaps if overlap at least this many bases 
int MinPhredScore;			// only accept reads for duplicate and KMer processing if mean Phred score is at least this threshold 

int NumberOfProcessors;		// number of installed CPUs
int NumThreads;				// number of threads (0 defaults to number of CPUs)
bool bAffinity;				// thread to core affinity

int NumPE1InputFiles;			// number of PE1 input files
char *pszInPE1files[cRSDMaxInFileSpecs];  // input PE1 5' read files
int NumPE2InputFiles;			// number of PE2 input files
char *pszInPE2files[cRSDMaxInFileSpecs];  // input PE2 3' read files

char szOutFile[_MAX_PATH];		// readset distributions to this CSV file
char szOutHTMLFile[_MAX_PATH];	// readset distributions to this file as HTML5

char szContaminantFile[_MAX_PATH];		// contaminants fasta file

char szSQLiteDatabase[_MAX_PATH];	// results summaries to this SQLite file
char szExperimentName[cMaxDatasetSpeciesChrom + 1];			// experiment name
char szExperimentDescr[1000];		// describes experiment

// command line args
struct arg_lit  *help = arg_lit0("h", "help", "print this help and exit");
struct arg_lit  *version = arg_lit0("v", "version,ver", "print version information and exit");
struct arg_int *FileLogLevel = arg_int0("f", "FileLogLevel", "<int>", "Level of diagnostics written to screen and logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F", "log", "<file>", "diagnostics log file");

struct arg_int *pmode = arg_int0("m", "mode", "<int>", "processing mode: 0 - independent processing of single/paired readsets, 1 - pooled processing of single/paired readsets");
struct arg_lit  *strand = arg_lit0("S", "strand", "strand specific processing");
struct arg_int *reqmaxdupseeds = arg_int0("s", "seeds", "<int>", "Use at most this number of seed reads for duplicate processing (default 5000000 range 100000..25000000)");
struct arg_int *trim5 = arg_int0("y", "trim5", "<int>", "trim this number of bases from 5' end of reads when loading raw reads (default is 0)");
struct arg_int *trim3 = arg_int0("Y", "trim3", "<int>", "trim this number of bases from 3' end of reads when loading raw reads (default is 0)");

struct arg_int *minphredscore = arg_int0("p", "minphred", "<int>", "only accept reads for duplicate and KMer processing if mean Phred score is at least this threshold (default 0 to ignore, range 10..40)");

struct arg_int *maxkmerlen = arg_int0("k", "maxkmerlen", "<int>", "maximum K-Mer length processing (default is 6, range 3..12)");
struct arg_int *kmerccc = arg_int0("K", "kmerccc", "<int>", "concordance correlation coefficient measure KMer length (default is 6, range 1..maxkmerlen)");

struct arg_int *maxcontamsubrate = arg_int0("z", "maxcontamsubrate", "<int>", "max allowed contaminant substitution rate (bases per 25bp of contaminant overlap, 1st 15bp of overlap no subs allowed) (default is 1, range 0..3)");
struct arg_int *mincontamlen = arg_int0("Z", "mincontamlen", "<int>", "accept contaminant overlaps if overlap of at least this many bases (default is 5, range 1..100)");

struct arg_file *inpe1files = arg_filen("i", "inpe1", "<file>", 1, cRSDMaxInFileSpecs, "Load single ended, or 5' if paired end, reads from fasta or fastq file(s); if single ended then wildcards allowed");
struct arg_file *inpe2files = arg_filen("u", "inpe2", "<file>", 0, cRSDMaxInFileSpecs, "Load 3' if paired end reads from fasta or fastq file(s)");
struct arg_file *contaminantfile = arg_file0("c", "contaminants", "<file>", "Putative contaminant sequences fasta file");
struct arg_file *outfile = arg_file0("o", "out", "<file>", "Output distributions to this CSV file");
struct arg_file *outhtmlfile = arg_file0("O", "out", "<file>", "Output distributions to this HTML5 file");

struct arg_int *threads = arg_int0("T", "threads", "<int>", "number of processing threads 0..64 (defaults to 0 which sets threads to number of cores)");

struct arg_file *summrslts = arg_file0("q", "sumrslts", "<file>", "Output results summary to this SQLite3 database file");
struct arg_str *experimentname = arg_str0("w", "experimentname", "<str>", "experiment name SQLite3 database file");
struct arg_str *experimentdescr = arg_str0("W", "experimentdescr", "<str>", "experiment description SQLite3 database file");

struct arg_end *end = arg_end(200);

void *argtable[] = { help, version, FileLogLevel, LogFile,
	pmode, strand, trim5, trim3,maxkmerlen,kmerccc, reqmaxdupseeds,minphredscore, maxcontamsubrate,mincontamlen,contaminantfile, inpe1files, inpe2files, outfile, // outhtmlfile,
	summrslts, experimentname, experimentdescr,
	threads,
	end };

char **pAllArgs;
int argerrors;
argerrors = CUtility::arg_parsefromfile(argc, (char **)argv, &pAllArgs);
if (argerrors >= 0)
argerrors = arg_parse(argerrors, pAllArgs, argtable);

/* special case: '--help' takes precedence over error reporting */
if (help->count > 0)
	{
	printf("\n%s %s %s, Version %s\nOptions ---\n", gszProcName, gpszSubProcess->pszName, gpszSubProcess->pszFullDescr, cpszProgVer);
	arg_print_syntax(stdout, argtable, "\n");
	arg_print_glossary(stdout, argtable, "  %-25s %s\n");
	printf("\nNote: Parameters can be entered into a parameter file, one parameter per line.");
	printf("\n      To invoke this parameter file then precede it's name with '@'");
	printf("\n      e.g. %s %s @myparams.txt\n", gszProcName, gpszSubProcess->pszName);
	printf("\nPlease report any issues regarding usage of %s to stuart.stephen@csiro.au\n\n", gszProcName);
	return(1);
	}

/* special case: '--version' takes precedence error reporting */
if (version->count > 0)
	{
	printf("\n%s %s Version %s\n", gszProcName, gpszSubProcess->pszName, cpszProgVer);
	return(1);
	}

if (!argerrors)
	{
	if (FileLogLevel->count && !LogFile->count)
		{
		printf("\nError: FileLogLevel '-f%d' specified but no logfile '-F<logfile>\n'", FileLogLevel->ival[0]);
		return(1);
		}

	iScreenLogLevel = iFileLogLevel = FileLogLevel->count ? FileLogLevel->ival[0] : eDLInfo;
	if (iFileLogLevel < eDLNone || iFileLogLevel > eDLDebug)
		{
		printf("\nError: FileLogLevel '-l%d' specified outside of range %d..%d\n", iFileLogLevel, eDLNone, eDLDebug);
		return(1);
		}

	if (LogFile->count)
		{
		strncpy(szLogFile, LogFile->filename[0], _MAX_PATH);
		szLogFile[_MAX_PATH - 1] = '\0';
		}
	else
		{
		iFileLogLevel = eDLNone;
		szLogFile[0] = '\0';
		}

	// now that log parameters have been parsed then initialise diagnostics log system
	if (!gDiagnostics.Open(szLogFile, (etDiagLevel)iScreenLogLevel, (etDiagLevel)iFileLogLevel, true))
		{
		printf("\nError: Unable to start diagnostics subsystem\n");
		if (szLogFile[0] != '\0')
			printf(" Most likely cause is that logfile '%s' can't be opened/created\n", szLogFile);
		return(1);
		}
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Subprocess %s Version %s starting", gpszSubProcess->pszName, cpszProgVer);
	gExperimentID = 0;
	gProcessID = 0;
	gProcessingID = 0;
	szSQLiteDatabase[0] = '\0';
	szExperimentName[0] = '\0';
	szExperimentDescr[0] = '\0';

	if (experimentname->count)
		{
		strncpy(szExperimentName, experimentname->sval[0], sizeof(szExperimentName));
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
	if (summrslts->count)
		{
		strncpy(szSQLiteDatabase, summrslts->filename[0], sizeof(szSQLiteDatabase)-1);
		szSQLiteDatabase[sizeof(szSQLiteDatabase)-1] = '\0';
		CUtility::TrimQuotedWhitespcExtd(szSQLiteDatabase);
		if (strlen(szSQLiteDatabase) < 1)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: After removal of whitespace, no SQLite database specified with '-q<filespec>' option");
			return(1);
			}

		if (strlen(szExperimentName) < 1)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: After removal of whitespace, no SQLite experiment name specified with '-w<str>' option");
			return(1);
			}
		if (experimentdescr->count)
			{
			strncpy(szExperimentDescr, experimentdescr->sval[0], sizeof(szExperimentDescr)-1);
			szExperimentDescr[sizeof(szExperimentDescr)-1] = '\0';
			CUtility::TrimQuotedWhitespcExtd(szExperimentDescr);
			}
		if (strlen(szExperimentDescr) < 1)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: After removal of whitespace, no SQLite experiment description specified with '-W<str>' option");
			return(1);
			}

		gExperimentID = gSQLiteSummaries.StartExperiment(szSQLiteDatabase, false, true, szExperimentName, szExperimentName, szExperimentDescr);
		if (gExperimentID < 1)
			return(1);
		gProcessID = gSQLiteSummaries.AddProcess((char *)gpszSubProcess->pszName, (char *)gpszSubProcess->pszName, (char *)gpszSubProcess->pszFullDescr);
		if (gProcessID < 1)
			return(1);
		gProcessingID = gSQLiteSummaries.StartProcessing(gExperimentID, gProcessID, (char *)cpszProgVer);
		if (gProcessingID < 1)
			return(1);
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "Initialised SQLite database '%s' for results summary collection", szSQLiteDatabase);
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "SQLite database experiment identifier for '%s' is %d", szExperimentName, gExperimentID);
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "SQLite database process identifier for '%s' is %d", (char *)gpszSubProcess->pszName, gProcessID);
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "SQLite database processing instance identifier is %d", gProcessingID);
		}
	else
		{
		szSQLiteDatabase[0] = '\0';
		szExperimentDescr[0] = '\0';
		}


	PMode = pmode->count ? pmode->ival[0] : (int)eRSDindependent;
	if (PMode < eRSDindependent || PMode > eRSDpooled)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: Processing mode '-m%d' must be in range 0..%d", PMode, eRSDpooled);
		return(1);
		}

	Trim5 = trim5->count ? trim5->ival[0] : 0;
	if (Trim5 < 0 || Trim5 > 50)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: Trim 5' raw reads '-y%d' specified outside of range %d..%d\n", Trim5, 0, 50);
		exit(1);
		}
	Trim3 = trim3->count ? trim3->ival[0] : 0;
	if (Trim3 < 0 || Trim3 > 50)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: Trim 3' raw reads '-y%d' specified outside of range %d..%d\n", Trim3, 0, 50);
		exit(1);
		}

	MaxKMerLen = maxkmerlen->count ? maxkmerlen->ival[0] : 6;
	if (MaxKMerLen < 3 || MaxKMerLen > cMaxKMerLen)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: Max K-Mer length '-k%d' specified outside of range %d..%d\n", MaxKMerLen, 3, cMaxKMerLen);
		exit(1);
		}

	MinPhredScore = minphredscore->count ? minphredscore->ival[0] : 0;
	if (MinPhredScore != 0 && (MinPhredScore < cMinBasePhredScore || MinPhredScore > 40))	// minimum Phred can be either 0, or in the range 10 to 40
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: Minimum Phred score '-p%d' specified outside of 0 (to ignore) or %d..40\n",cMinBasePhredScore, MinPhredScore);
		exit(1);
		}

	KMerCCC = kmerccc->count ? kmerccc->ival[0] : min(6,MaxKMerLen);
	if (KMerCCC < 1 || KMerCCC > MaxKMerLen)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: concordance correlation coefficient measure KMer length '-K%d' specified outside of range %d..%d\n", KMerCCC, 1, MaxKMerLen);
		exit(1);
		}

	MaxContamSubRate = maxcontamsubrate->count ? maxcontamsubrate->ival[0] : 1;
	if (MaxContamSubRate < 0 || MaxContamSubRate > 3)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: Max contaminant rate '-z%d' specified outside of range 0..3\n", MaxContamSubRate);
		exit(1);
		}

	MinContamLen = mincontamlen->count ? mincontamlen->ival[0] : 5;
	if (MinContamLen < 1 || MinContamLen > 100)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: Min contaminant overlap length '-Z%d' specified outside of range 1..100\n", MinContamLen);
		exit(1);
		}

	ReqMaxDupSeeds = reqmaxdupseeds->count ? reqmaxdupseeds->ival[0] : cDfltDupSeeds;
	if (ReqMaxDupSeeds < cMinDupSeeds || ReqMaxDupSeeds > cMaxDupSeeds)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: Max seed reads for duplicates '-s%d' must be in range %d..%d", ReqMaxDupSeeds, cMinDupSeeds, cMaxDupSeeds);
		return(1);
		}
	bStrand = strand->count ? true : false;


	if (!inpe1files->count)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: No input file(s) specified with with '-i<filespec>' option)");
		return(1);
		}

	for (NumPE1InputFiles = Idx = 0; NumPE1InputFiles < cRSDMaxInFileSpecs && Idx < inpe1files->count; Idx++)
		{
		pszInPE1files[Idx] = NULL;
		if (pszInPE1files[NumPE1InputFiles] == NULL)
			pszInPE1files[NumPE1InputFiles] = new char[_MAX_PATH];
		strncpy(pszInPE1files[NumPE1InputFiles], inpe1files->filename[Idx], _MAX_PATH);
		pszInPE1files[NumPE1InputFiles][_MAX_PATH - 1] = '\0';
		CUtility::TrimQuotedWhitespcExtd(pszInPE1files[NumPE1InputFiles]);
		if (pszInPE1files[NumPE1InputFiles][0] != '\0')
			NumPE1InputFiles++;
		}

	if (!NumPE1InputFiles)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: After removal of whitespace, no input file(s) specified with '-i<filespec>' option");
		return(1);
		}

	NumPE2InputFiles = 0;
	if (NumPE1InputFiles > 0 && inpe2files->count)
		{
		for (NumPE2InputFiles = Idx = 0; NumPE2InputFiles < cRSDMaxInFileSpecs && Idx < inpe2files->count; Idx++)
			{
			pszInPE2files[Idx] = NULL;
			if (pszInPE2files[NumPE2InputFiles] == NULL)
				pszInPE2files[NumPE2InputFiles] = new char[_MAX_PATH];
			strncpy(pszInPE2files[NumPE2InputFiles], inpe2files->filename[Idx], _MAX_PATH);
			pszInPE2files[NumPE2InputFiles][_MAX_PATH - 1] = '\0';
			CUtility::TrimQuotedWhitespcExtd(pszInPE2files[NumPE2InputFiles]);
			if (pszInPE2files[NumPE2InputFiles][0] != '\0')
				NumPE2InputFiles++;
			}	
		if (!NumPE2InputFiles)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: After removal of whitespace, no PE2 3' input file(s) specified with '-u<filespec>' option");
			return(1);
			}

		if (NumPE1InputFiles != NumPE2InputFiles)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: Paired end processing, expected number of PE1 files (%d) to equal number of PE2 files (%d)", NumPE1InputFiles, NumPE2InputFiles);
			return(1);
			}
		}

	if (outfile->count)
		{
		strncpy(szOutFile, outfile->filename[0], _MAX_PATH);
		szOutFile[_MAX_PATH - 1] = '\0';
		CUtility::TrimQuotedWhitespcExtd(szOutFile);
		}
	else
		szOutFile[0] = '\0';
	
//	if(outhtmlfile->count)
//		{
//		strncpy(szOutHTMLFile, outhtmlfile->filename[0], _MAX_PATH);
//		szOutHTMLFile[_MAX_PATH - 1] = '\0';
//		CUtility::TrimQuotedWhitespcExtd(szOutHTMLFile);
//		}
//	else
		szOutHTMLFile[0] = '\0';

	if (szOutFile[0] == '\0' /* && szOutHTMLFile[0] == '\0' */)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: Expected output file '-o<outfile>' to be specified");
		return(1);
		}

	if (contaminantfile->count)
		{
		strncpy(szContaminantFile, contaminantfile->filename[0], _MAX_PATH);
		szContaminantFile[_MAX_PATH - 1] = '\0';
		CUtility::TrimQuotedWhitespcExtd(szContaminantFile);
		}
	else
		{
		szContaminantFile[0] = '\0';
		MaxContamSubRate = 0;
		MinContamLen = 0;
		}

#ifdef _WIN32
	SYSTEM_INFO SystemInfo;
	GetSystemInfo(&SystemInfo);
	NumberOfProcessors = SystemInfo.dwNumberOfProcessors;
	bAffinity = false;
#else
	bAffinity = false;
	NumberOfProcessors = sysconf(_SC_NPROCESSORS_ONLN);
#endif
	int MaxAllowedThreads = min(cRSDMaxWorkerThreads, NumberOfProcessors);	// limit to be at most cMaxWorkerThreads
	if ((NumThreads = threads->count ? threads->ival[0] : MaxAllowedThreads) == 0)
		NumThreads = MaxAllowedThreads;
	if (NumThreads < 0 || NumThreads > MaxAllowedThreads)
		{
		gDiagnostics.DiagOut(eDLWarn, gszProcName, "Warning: Number of threads '-T%d' specified was outside of range %d..%d", NumThreads, 1, MaxAllowedThreads);
		gDiagnostics.DiagOut(eDLWarn, gszProcName, "Warning: Defaulting number of threads to %d", MaxAllowedThreads);
		NumThreads = MaxAllowedThreads;
		}

	if (bAffinity)
		{
		// try and set thread affinities
#ifndef _WIN32
		int AffinityIdx;
		int CoresPerCPU;
		int ReqAffinity;
		int ReqCores = NumThreads;
		pthread_t Me = pthread_self();
		cpu_set_t CurCpuSet;
		cpu_set_t NewCpuSet;
		CPU_ZERO(&NewCpuSet);
		pthread_getaffinity_np(Me, sizeof(CurCpuSet), &CurCpuSet);		// get current affinity set
		CoresPerCPU = 0;
		for (AffinityIdx = 0; AffinityIdx < CPU_SETSIZE; AffinityIdx++)	// determine number of CPUs to which affinity is enabled
			{
			if (CPU_ISSET(AffinityIdx, &CurCpuSet))
				CoresPerCPU += 1;
			}
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "Info: %d CPUs (%d cores) available", CoresPerCPU, NumberOfProcessors);
		// number of CPUs with affinity known, estimate cores per CPU

		CoresPerCPU = NumberOfProcessors / CoresPerCPU;

		// now cores per CPU has been guestimated then try and set affinity for CPU's
		ReqAffinity = 0;
		for (AffinityIdx = 0; AffinityIdx < CPU_SETSIZE; AffinityIdx++)
			{
			if (CPU_ISSET(AffinityIdx, &CurCpuSet))
				{
				CPU_SET(AffinityIdx, &NewCpuSet);
				ReqCores -= CoresPerCPU;
				ReqAffinity += 1;
				if (ReqCores < 1)
					break;
				}
			}
		pthread_setaffinity_np(Me, sizeof(NewCpuSet), &NewCpuSet);
#endif
		}

	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Processing parameters:");
	const char *pszDescr;
	switch(PMode)
		{
		case eRSDindependent:
			pszDescr = "Independent processing of single or paired end readsets";
			break;
		case eRSDpooled:
			pszDescr = "Pooled processing of single or paired end readsets";
			break;
		}

	gDiagnostics.DiagOutMsgOnly(eDLInfo, "Processing mode is : '%s'", pszDescr);

	gDiagnostics.DiagOutMsgOnly(eDLInfo, "Max number of duplicate seed reads : %d", ReqMaxDupSeeds);

	gDiagnostics.DiagOutMsgOnly(eDLInfo, "Min Phred threshold score for duplicate and K-mer distributions : %d", MinPhredScore);
	

	gDiagnostics.DiagOutMsgOnly(eDLInfo, "trim 5' ends raw reads by : %d", Trim5);
	gDiagnostics.DiagOutMsgOnly(eDLInfo, "trim 3' ends raw reads by : %d", Trim3);

	gDiagnostics.DiagOutMsgOnly(eDLInfo, "processing is for upto this KMer length inclusive : %d", MaxKMerLen);
	gDiagnostics.DiagOutMsgOnly(eDLInfo, "concordance correlation coefficient measure KMer length : %d", KMerCCC);


	gDiagnostics.DiagOutMsgOnly(eDLInfo, "Strand specific processing : '%s'", bStrand ? "Yes" : "No");

	if (szContaminantFile[0] != '\0')
		{
		gDiagnostics.DiagOutMsgOnly(eDLInfo, "Max contaminant rate : %d", MaxContamSubRate);
		gDiagnostics.DiagOutMsgOnly(eDLInfo, "Min contaminant overlap length : %d", MinContamLen);
		gDiagnostics.DiagOutMsgOnly(eDLInfo, "Contaminant sequences file: %s", szContaminantFile);
		}

	if (NumPE2InputFiles)
		{
		for (Idx = 0; Idx < NumPE1InputFiles; Idx++)
			gDiagnostics.DiagOutMsgOnly(eDLInfo, "Input PE 5' reads file (%d) : '%s'", Idx + 1, pszInPE1files[Idx]);

		for (Idx = 0; Idx < NumPE2InputFiles; Idx++)
			gDiagnostics.DiagOutMsgOnly(eDLInfo, "Input PE 3' reads file (%d) : '%s'", Idx + 1, pszInPE2files[Idx]);
		}
	else
		for (Idx = 0; Idx < NumPE1InputFiles; Idx++)
			gDiagnostics.DiagOutMsgOnly(eDLInfo, "Input single ended reads file (%d) : '%s'", Idx + 1, pszInPE1files[Idx]);

	if (szOutHTMLFile[0] != '\0')
			gDiagnostics.DiagOutMsgOnly(eDLInfo, "Output HTML5 distributions file: %s", szOutHTMLFile);

	if (szOutFile[0] != '\0')
		gDiagnostics.DiagOutMsgOnly(eDLInfo, "Output CSV distributions file: %s", szOutFile);

	if (szExperimentName[0] != '\0')
		gDiagnostics.DiagOutMsgOnly(eDLInfo, "This processing reference: %s", szExperimentName);

	gDiagnostics.DiagOutMsgOnly(eDLInfo, "Number of threads : %d", NumThreads);

	if (gExperimentID > 0)
		{
		int ParamID;
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID, ePTText, (int)strlen(szLogFile), "log", szLogFile);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID, ePTInt32, sizeof(PMode), "mode", &PMode);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID, ePTBool, sizeof(PMode), "strand", &bStrand);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID, ePTInt32, sizeof(ReqMaxDupSeeds), "maxdupseeds", &ReqMaxDupSeeds);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID, ePTInt32, (int)sizeof(Trim5), "trim5", &Trim5);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID, ePTInt32, (int)sizeof(Trim3), "trim3", &Trim3);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID, ePTInt32, (int)sizeof(MaxKMerLen), "maxkmerlen", &MaxKMerLen);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID, ePTInt32, (int)sizeof(MinPhredScore), "minphredscore", &MinPhredScore);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID, ePTInt32, (int)sizeof(KMerCCC), "kmerccc", &KMerCCC);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID, ePTInt32, sizeof(NumPE1InputFiles), "NumPE1InputFiles", &NumPE1InputFiles);
		for (Idx = 0; Idx < NumPE1InputFiles; Idx++)
			ParamID = gSQLiteSummaries.AddParameter(gProcessingID, ePTText, (int)strlen(pszInPE1files[Idx]), "inpe1", pszInPE1files[Idx]);

		ParamID = gSQLiteSummaries.AddParameter(gProcessingID, ePTInt32, sizeof(NumPE2InputFiles), "NumPE2InputFiles", &NumPE2InputFiles);
		for (Idx = 0; Idx < NumPE2InputFiles; Idx++)
			ParamID = gSQLiteSummaries.AddParameter(gProcessingID, ePTText, (int)strlen(pszInPE2files[Idx]), "inpe2", pszInPE2files[Idx]);

		if (szOutFile[0] != '\0')
			ParamID = gSQLiteSummaries.AddParameter(gProcessingID, ePTText, (int)strlen(szOutFile), "out", szOutFile);

		if (szOutHTMLFile[0] != '\0')
			ParamID = gSQLiteSummaries.AddParameter(gProcessingID, ePTText, (int)strlen(szOutHTMLFile), "outhtmlfile", szOutHTMLFile);

		if (szContaminantFile[0] != '\0')
			{
			ParamID = gSQLiteSummaries.AddParameter(gProcessingID, ePTInt32, (int)sizeof(MaxContamSubRate), "maxcontamsubrate", &MaxContamSubRate);
			ParamID = gSQLiteSummaries.AddParameter(gProcessingID, ePTInt32, (int)sizeof(MinContamLen), "mincontamlen", &MinContamLen);
			ParamID = gSQLiteSummaries.AddParameter(gProcessingID, ePTText, (int)strlen(szContaminantFile), "contaminantfile", szContaminantFile);
			}			

		ParamID = gSQLiteSummaries.AddParameter(gProcessingID, ePTInt32, sizeof(NumThreads), "threads", &NumThreads);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID, ePTInt32, sizeof(NumberOfProcessors), "cpus", &NumberOfProcessors);

		ParamID = gSQLiteSummaries.AddParameter(gProcessingID, ePTText, (int)strlen(szSQLiteDatabase), "sumrslts", szSQLiteDatabase);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID, ePTText, (int)strlen(szExperimentName), "experimentname", szExperimentName);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID, ePTText, (int)strlen(szExperimentDescr), "experimentdescr", szExperimentDescr);
		}

#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	gStopWatch.Start();

	Rslt = ProcessReadsetDist((etRSDMode)PMode,			// processing mode; eRSDindependent or eRSDpooled
							  bStrand,					// true if read strand specific distributions
							  Trim5,					// trim this number of bases from 5' end of reads when loading the reads
							  Trim3,					// trim this number of bases from 3' end of reads when loading the reads
			  				  MaxKMerLen,				// processing is for upto this KMer length inclusive
							  KMerCCC,					// concordance correlation coefficient measure KMer length
							  MaxContamSubRate,			// max allowed contamimant substitution rate (bases per 25bp of contaminant overlap, 1st 15bp of overlap no subs allowed)
							  MinContamLen,				// accept contaminant overlaps if overlap at least this many bases 
							  ReqMaxDupSeeds,			// requested to sample for this many duplicate seeds and off target alignments
							  MinPhredScore,			// only accept reads for duplicate and KMer processing if their mean Phred score is at least this threshold 
							  NumThreads,				// number of worker threads to use
							  bAffinity,				// thread to core affinity
							  NumPE1InputFiles,			// number of PE1 input files
							  pszInPE1files,			// input PE1 5' read files
							  NumPE2InputFiles,			// number of PE2 input files
							  pszInPE2files,		    // input PE2 3' read files
							  szContaminantFile,		// contaminants fasta file
							  szOutFile,				// where to write distributions CSV file
							  szOutHTMLFile);			//  where to write distributions HTML5 file

	Rslt = Rslt >= 0 ? 0 : 1;

	if (gExperimentID > 0)
		{
		if (gProcessingID)
			gSQLiteSummaries.EndProcessing(gProcessingID, Rslt);
		gSQLiteSummaries.EndExperiment(gExperimentID);
		}
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Exit code: %d Total processing time: %s", Rslt, gStopWatch.Read());
	return(Rslt);
	}
else
	{
	printf("\n%s %s %s, Version %s\n", gszProcName, gpszSubProcess->pszName, gpszSubProcess->pszFullDescr, cpszProgVer);
	arg_print_errors(stdout, end, gszProcName);
	arg_print_syntax(stdout, argtable, "\nUse '-h' to view option and parameter usage\n");
	return(1);
	}
return 0;
}

// Thread startup
#ifdef _WIN32
unsigned __stdcall IndependentProcNGSQC(void * pThreadPars)
#else
void *IndependentProcNGSQC(void * pThreadPars)
#endif
{
CReadStats *pReadStats;
int Rslt = 0;
tsThreadIndependentNGSQCPars *pPars = (tsThreadIndependentNGSQCPars *)pThreadPars; // makes it easier not having to deal with casts!
if((pReadStats = new CReadStats()) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate instance of CReadStats");
	pPars->Rslt = eBSFerrObj;
#ifdef _WIN32
	_endthreadex(0);
	return(eBSFSuccess);
#else
	pthread_exit(NULL);
#endif
	}

Rslt = pReadStats->ProcessReadsetDist(pPars->PMode,			// processing mode; eRSDindependent or eRSDpooled
						pPars->ProcessingID,							// processing instance identifier, used if processing eRSDindependent to identify output file instances
						pPars->bStrand,						// true if read strand specific distributions
						pPars->Trim5,						// trim this number of bases from 5' end of reads when loading the reads
						pPars->Trim3,						// trim this number of bases from 3' end of reads when loading the reads
						pPars->MaxKMerLen,					// processing is for upto this KMer length inclusive
						pPars->KMerCCC,						// concordance correlation coefficient measure KMer length
						pPars->MaxContamSubRate,			// max allowed contamimant substitution rate (bases per 25bp of contaminant overlap, 1st 15bp of overlap no subs allowed)
						pPars->MinContamLen,				// accept contaminant overlaps if overlap at least this many bases 
						pPars->ReqMaxDupSeeds,				// requested to sample this many reads for duplicates and off target alignments
						pPars->MinPhredScore,				// only accept reads for duplicate and KMer processing if their mean Phred score is at least this threshold
						pPars->NumThreads,					// number of worker threads to use
						pPars->bAffinity,					// thread to core affinity
						pPars->NumPE1InputFiles,			// number of PE1 input files
						pPars->ppszInPE1files,				// input PE1 5' read files
						pPars->NumPE2InputFiles,			// number of PE2 input files
						pPars->ppszInPE2files,				// input PE2 3' read files
						pPars->pszContaminantFile,			// contaminants fasta file
						pPars->pszOutDistFile,				// where to write distributions CSV file
						pPars->pszOutHTMLFile);				// where to write distributions HTML5 file

delete pReadStats;
pPars->Rslt = Rslt;
#ifdef _WIN32
_endthreadex(0);
return(eBSFSuccess);
#else
pthread_exit(NULL);
#endif
}

int
ProcessReadsetDist(etRSDMode PMode,					// processing mode; eRSDindependent or eRSDpooled
					bool bStrand,					// true if read strand specific distributions
					int Trim5,						// trim this number of bases from 5' end of reads when loading the reads
					int Trim3,						// trim this number of bases from 3' end of reads when loading the reads
					int MaxKMerLen,					// processing is for upto this KMer length inclusive
					int KMerCCC,					// concordance correlation coefficient measure KMer length
					int MaxContamSubRate,			// max allowed contamimant substitution rate (bases per 25bp of contaminant overlap, 1st 15bp of overlap no subs allowed)
					int MinContamLen,				// accept contaminant overlaps if overlap at least this many bases 
					int ReqMaxDupSeeds,				// requested to sample for this many duplicate seeds and off target alignments
					int MinPhredScore,				// only accept reads for duplicate and KMer processing if their mean Phred score is at least this threshold 
					int NumThreads,					// number of worker threads to use
					bool bAffinity,					// thread to core affinity
					int NumPE1InputFiles,			// number of PE1 input files
					char *pszInPE1files[],			// input PE1 5' read files
					int NumPE2InputFiles,			// number of PE2 input files
					char *pszInPE2files[],		    // input PE2 3' read files
					char *pszcontaminantFile,		// contaminants fasta file
					char *pszOutDistFile,			// where to write distributions CSV file
					char *pszOutHTMLFile)			//  where to write distributions HTML5 file
{
int Rslt;
int NumErrs;
NumErrs = 0;

// ensure low level PLPlot calls will be serialised - I'm not convinced calls into PLPlot are thread safe
CBKPLPlot *pBKPLPlot;
if((pBKPLPlot = new CBKPLPlot()) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate instance of CBKPLPlot");
	return(-1);
	}
pBKPLPlot->CreateSerialise();								

if(PMode == eRSDpooled || NumPE1InputFiles == 1) 			// if pooling all input readsets or just 1 readset
	{
	CReadStats *pReadStats;
	if((pReadStats = new CReadStats()) == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate instance of CReadStats");
		delete pBKPLPlot;
		return(-1);
		}

	Rslt = pReadStats->ProcessReadsetDist(PMode,		// processing mode; eRSDindependent or eRSDpooled
						1,							// processing instance identifier, used if processing eRSDindependent to identify output file instances
						bStrand,					// true if read strand specific distributions
						Trim5,						// trim this number of bases from 5' end of reads when loading the reads
						Trim3,						// trim this number of bases from 3' end of reads when loading the reads
						MaxKMerLen,					// processing is for upto this KMer length inclusive
						KMerCCC,					// concordance correlation coefficient measure KMer length
						MaxContamSubRate,			// max allowed contamimant substitution rate (bases per 25bp of contaminant overlap, 1st 15bp of overlap no subs allowed)
						MinContamLen,				// accept contaminant overlaps if overlap at least this many bases 
						ReqMaxDupSeeds,				// requested to sample this many reads for duplicates and off target alignments
						MinPhredScore,				// only accept reads for duplicate and KMer processing if their mean Phred score is at least this threshold
						NumThreads,					// number of worker threads to use
						bAffinity,					// thread to core affinity
						NumPE1InputFiles,			// number of PE1 input files
						pszInPE1files,				// input PE1 5' read files
						NumPE2InputFiles,			// number of PE2 input files
						pszInPE2files,				// input PE2 3' read files
						pszcontaminantFile,			// contaminants fasta file
						pszOutDistFile,				// where to write distributions CSV file
						pszOutHTMLFile);			//  where to write distributions HTML5 file
	if(Rslt < 0)
		NumErrs = 1;
	delete pReadStats;
	}
else
	{
	tsThreadIndependentNGSQCPars *pNGSQCThreads;
	tsThreadIndependentNGSQCPars *pCurNGSQCThread;
	tsThreadIndependentNGSQCPars *pNxtNGSQCThread;
	char szFName[_MAX_PATH];
	char *pszAllocOutStatsFiles;
	char *pszOutStatsFile;

	int ProcessingID;
	int NextToProcessID;
	int NumCompletedProc;
	Rslt = 0;

	if(NumThreads > NumPE1InputFiles)
		NumThreads = NumPE1InputFiles;

	pNGSQCThreads = new tsThreadIndependentNGSQCPars[NumPE1InputFiles];
	memset(pNGSQCThreads,0,sizeof(tsThreadIndependentNGSQCPars) * NumPE1InputFiles);
	pszAllocOutStatsFiles = new char [_MAX_PATH * 2 * NumPE1InputFiles];
	memset(pszAllocOutStatsFiles,0,_MAX_PATH * 2 * NumPE1InputFiles);
	pCurNGSQCThread = pNGSQCThreads;
	pszOutStatsFile = pszAllocOutStatsFiles;
	for(ProcessingID = 1; ProcessingID <= NumPE1InputFiles; ProcessingID++,pCurNGSQCThread++,pszOutStatsFile += _MAX_PATH * 2)
		{
		strcpy(pszOutStatsFile,pszOutDistFile);
		int PathLen = (int)strlen(pszOutStatsFile);
		if(!(pszOutStatsFile[PathLen-1] == '.' || pszOutStatsFile[PathLen-1] == '_'))
			{
			pszOutStatsFile[PathLen++] = '_';
			pszOutStatsFile[PathLen] = '\0';
			}

		CUtility::splitpath(pszInPE1files[ProcessingID-1], NULL, szFName);
		sprintf(&pszOutStatsFile[PathLen],"%d_%s",ProcessingID,szFName);
		pszOutStatsFile[_MAX_PATH-1] = '\0';
		pCurNGSQCThread->bProcCompleted = false;
		pCurNGSQCThread->PMode = PMode;
		pCurNGSQCThread->ProcessingID = ProcessingID;
		pCurNGSQCThread->bStrand = bStrand;
		pCurNGSQCThread->Trim5 = Trim5;
		pCurNGSQCThread->Trim3 = Trim3;
		pCurNGSQCThread->MaxKMerLen = MaxKMerLen;
		pCurNGSQCThread->KMerCCC = KMerCCC;
		pCurNGSQCThread->MaxContamSubRate = MaxContamSubRate;
		pCurNGSQCThread->MinContamLen = MinContamLen;
		pCurNGSQCThread->ReqMaxDupSeeds = ReqMaxDupSeeds;
		pCurNGSQCThread->MinPhredScore = MinPhredScore;
		pCurNGSQCThread->NumThreads = 1;
		pCurNGSQCThread->bAffinity = bAffinity;
		pCurNGSQCThread->NumPE1InputFiles = 1;
		pCurNGSQCThread->ppszInPE1files = &pszInPE1files[ProcessingID-1];
		pCurNGSQCThread->NumPE2InputFiles = NumPE2InputFiles == 0 ? 0 : 1;
		pCurNGSQCThread->ppszInPE2files = NumPE2InputFiles == 0 ? pszInPE2files : &pszInPE2files[ProcessingID-1];
		pCurNGSQCThread->pszContaminantFile = pszcontaminantFile;
		pCurNGSQCThread->pszOutDistFile = pszOutStatsFile;
		pCurNGSQCThread->pszOutHTMLFile = pszOutHTMLFile;
		}

	pCurNGSQCThread = pNGSQCThreads;
	for (ProcessingID = 1; ProcessingID <= NumThreads; ProcessingID++, pCurNGSQCThread++)
		{
		pCurNGSQCThread->ThreadIdx = ProcessingID;
#ifdef _WIN32
		pCurNGSQCThread->threadHandle = (HANDLE)_beginthreadex(NULL, 0x0fffff, IndependentProcNGSQC, pCurNGSQCThread, 0, &pCurNGSQCThread->threadID);
#else
		pCurNGSQCThread->threadRslt = pthread_create(&pCurNGSQCThread->threadID, NULL, IndependentProcNGSQC, pCurNGSQCThread);
#endif
		}
	pNxtNGSQCThread = pCurNGSQCThread;
	NextToProcessID = ProcessingID;
	NumCompletedProc = 0;

	do {
#ifdef _WIN32
		Sleep(5000);
#else
		sleep(5);
#endif
		pCurNGSQCThread = pNGSQCThreads;
		for (ProcessingID = 1; ProcessingID <= NumPE1InputFiles; ProcessingID++, pCurNGSQCThread++)
			{
			if(pCurNGSQCThread->ThreadIdx == 0 || pCurNGSQCThread->bProcCompleted)
				continue;
#ifdef _WIN32
			if(WAIT_TIMEOUT != WaitForSingleObject(pCurNGSQCThread->threadHandle, 2000)) 
				{
				CloseHandle(pCurNGSQCThread->threadHandle);
				if(pCurNGSQCThread->Rslt < 0)
					{
					NumErrs += 1;
					Rslt = -1;
					}
				pCurNGSQCThread->bProcCompleted = true;
				NumCompletedProc += 1;
				if(NextToProcessID <= NumPE1InputFiles)
					{
					pNxtNGSQCThread->ThreadIdx = NextToProcessID++;
					pNxtNGSQCThread->threadHandle = (HANDLE)_beginthreadex(NULL, 0x0fffff, IndependentProcNGSQC, pNxtNGSQCThread, 0, &pNxtNGSQCThread->threadID);
					pNxtNGSQCThread += 1;
					}
				}
#else
			struct timespec ts;
			int JoinRlt;
			clock_gettime(CLOCK_REALTIME, &ts);
			ts.tv_sec += 2;
			if ((JoinRlt = pthread_timedjoin_np(pCurNGSQCThread->threadID, NULL, &ts)) == 0)
				{
				if(pCurNGSQCThread->Rslt < 0)
					{
					NumErrs += 1;
					Rslt = -1;
					}
				pCurNGSQCThread->bProcCompleted = true;
				NumCompletedProc += 1;
				if(NextToProcessID <= NumPE1InputFiles)
					{
					pNxtNGSQCThread->ThreadIdx = NextToProcessID++;
					pNxtNGSQCThread->threadRslt = pthread_create(&pNxtNGSQCThread->threadID, NULL, IndependentProcNGSQC, pNxtNGSQCThread);
					pNxtNGSQCThread += 1;
					}
				}
#endif
			}
		}
	while(NumCompletedProc < NumPE1InputFiles);
	delete pNGSQCThreads;
	delete pszAllocOutStatsFiles;
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"There were %d readsets with processing errors",NumErrs);
delete pBKPLPlot;
return(Rslt);
}

CReadStats::CReadStats()
{
m_pSeqHashes = NULL;
m_pSampledSeqs = NULL;
m_pContaminates = NULL;
m_pBaseNs = NULL; 
m_pScores = NULL; 
m_pKMerCnts = NULL;
m_bMutexesCreated = false;
m_hContamRptFile = -1;
m_hKMerDistRptFile = -1;
m_hPearsonDistRptFile = -1;
m_hQScoreDistRptFile = -1;
m_hErrFreeReadDistRptFile = -1;
m_hDuplicatesDistRptFile = -1;
m_hReadLenDistRptFile = -1;
Init();
}


CReadStats::~CReadStats()
{
Reset();
}

void 
CReadStats::Init(void)		// initialisation
{
m_PMode = eRSDindependent;
m_bStrand = false;
m_ReqMaxDupSeeds = 0;
m_ActMaxDupSeeds = 0;
m_NumThreads = 1;
m_bAffinity = false;
m_NumPE1InputFiles = 0;
m_ppszInPE1files = NULL;	
m_NumPE2InputFiles =0;
m_ppszInPE2files = NULL;
m_pszContaminantFile = NULL;
m_pszOutDistFile =NULL;	
m_pszOutHTMLFile =NULL;

m_NumChkdPE1ContamHits = 0;
m_NumPE1ContamHits = 0;
m_NumChkdPE2ContamHits = 0;
m_NumPE2ContamHits = 0;

m_MaxContamSubRate = cDfltContamSubRate;
m_MinContamLen = cDfltMinContamLen; 
m_hContamRptFile = -1;
m_hKMerDistRptFile = -1;
m_hPearsonDistRptFile = -1;
m_hQScoreDistRptFile = -1;
m_hErrFreeReadDistRptFile = -1;
m_hDuplicatesDistRptFile = -1;
m_hReadLenDistRptFile = -1;

m_szContamRptFile[0] = 0;
m_szKMerDistRptFile[0] = 0;
m_szPearsonDistRptFile[0] = 0;
m_szQScoreDistRptFile[0] = 0;
m_szDuplicatesDistRptFile[0] = '\0';
m_szReadLenDistRptFile[0] = '\0';

m_EstPE1MeanReadLen = 0;
m_EstPE2MeanReadLen = 0;

m_MinReadLen = 0;
m_MaxReadLen = 0;
memset(m_ReadLenDist,0,sizeof(m_ReadLenDist));
memset(m_ProbNoReadErrDist,0,sizeof(m_ProbNoReadErrDist));

m_NumSeqHashes = 0;
m_AllocdSampledSeqMem = 0;
m_AllocdSampledSeqWrds = 0;
m_UsedSampledSeqWrds = 0;

m_AllocdMaxReadLen = 0;
m_EstMaxSeqLen = 0;

m_bPEProc = false;
m_NumInFiles = 0;
m_NumPE1InFiles = 0;
m_NumPE2InFiles = 0;

m_bTerminate = false;
memset(m_InReadsFiles, 0, sizeof(m_InReadsFiles));
}

void 
CReadStats::Reset(void)		// reset state back to that immediately following initialisation
{
if(m_hContamRptFile != -1)
	{
#ifdef _WIN32
	_commit(m_hContamRptFile);
#else
	fsync(m_hContamRptFile);
#endif
	close(m_hContamRptFile);
	m_hContamRptFile = -1;
	}

if(m_hReadLenDistRptFile != -1)
	{
#ifdef _WIN32
	_commit(m_hReadLenDistRptFile);
#else
	fsync(m_hReadLenDistRptFile);
#endif
	close(m_hReadLenDistRptFile);
	m_hReadLenDistRptFile = -1;
	}

if(m_hKMerDistRptFile != -1)
	{
#ifdef _WIN32
	_commit(m_hKMerDistRptFile);
#else
	fsync(m_hKMerDistRptFile);
#endif
	close(m_hKMerDistRptFile);
	m_hKMerDistRptFile = -1;
	}

if(m_hPearsonDistRptFile != -1)
	{
#ifdef _WIN32
	_commit(m_hPearsonDistRptFile);
#else
	fsync(m_hPearsonDistRptFile);
#endif
	close(m_hPearsonDistRptFile);
	m_hPearsonDistRptFile = -1;
	}

if(m_hQScoreDistRptFile != -1)
	{
#ifdef _WIN32
	_commit(m_hQScoreDistRptFile);
#else
	fsync(m_hQScoreDistRptFile);
#endif
	close(m_hQScoreDistRptFile);
	m_hQScoreDistRptFile = -1;
	}

if(m_hErrFreeReadDistRptFile != -1)
	{
#ifdef _WIN32
	_commit(m_hErrFreeReadDistRptFile);
#else
	fsync(m_hErrFreeReadDistRptFile);
#endif
	close(m_hErrFreeReadDistRptFile);
	m_hErrFreeReadDistRptFile = -1;
	}

if(m_hDuplicatesDistRptFile != -1)
	{
#ifdef _WIN32
	_commit(m_hDuplicatesDistRptFile);
#else
	fsync(m_hDuplicatesDistRptFile);
#endif
	close(m_hDuplicatesDistRptFile);
	m_hDuplicatesDistRptFile = -1;
	}

if(m_pSeqHashes != NULL)
	{
	delete m_pSeqHashes;
	m_pSeqHashes = NULL;
	}
if(m_pSampledSeqs != NULL)
	{
#ifdef _WIN32
	free(m_pSampledSeqs);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pSampledSeqs != MAP_FAILED)
		munmap(m_pSampledSeqs,m_AllocdSampledSeqMem);
#endif
	m_pSampledSeqs = NULL;
	}

if(m_pKMerCnts != NULL)
	{
#ifdef _WIN32
	free(m_pKMerCnts);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pKMerCnts != MAP_FAILED)
		munmap(m_pKMerCnts,m_AllocdKMerCntsMem);
#endif
	m_pKMerCnts = NULL;
	}

if(m_pContaminates != NULL)
	{
	delete m_pContaminates;
	m_pContaminates = NULL;
	}

if(m_pBaseNs != NULL)
	{
	delete m_pBaseNs;
	m_pBaseNs = NULL;
	}
if(m_pScores != NULL)
	{
	delete m_pScores;
	m_pScores = NULL;
	}
m_AllocdSampledSeqMem = 0;
m_AllocdKMerCntsMem = 0;
DeleteMutexes();
Init();
}

// Thread startup
#ifdef _WIN32
unsigned __stdcall ThreadedNGSQC(void * pThreadPars)
#else
void *ThreadedNGSQC(void * pThreadPars)
#endif
{
	int Rslt = 0;
	tsThreadNGSQCPars *pPars = (tsThreadNGSQCPars *)pThreadPars; // makes it easier not having to deal with casts!
	CReadStats *pThis = (CReadStats *)pPars->pThis;
	Rslt = pThis->ProcNGSQC(pPars);
	pPars->Rslt = Rslt;
#ifdef _WIN32
	_endthreadex(0);
	return(eBSFSuccess);
#else
	pthread_exit(NULL);
#endif
}



int				// <0 if errors, 0 if no matches, >0 at least one contaminant sequence overlapping
CReadStats::LocateContaminentMatch(int SeqLen,			// targ sequence is of this length
								etSeqBase *pSeq,		// target sequence	
								bool bPE2)				// false if processing SE/PE1 read, true if PE2	
{
int Rslt;
if(m_pContaminates == NULL)
	return(0);

if(!bPE2)
	{
	if((Rslt = m_pContaminates->MatchContaminants(eAOF5PE1Targ,m_MaxContamSubRate,m_MinContamLen,SeqLen,pSeq))==0)
		Rslt = m_pContaminates->MatchContaminants(eAOF3PE1Targ,m_MaxContamSubRate,m_MinContamLen,SeqLen,pSeq);
	}
else
	{
	if((Rslt = m_pContaminates->MatchContaminants(eAOF5PE2Targ,m_MaxContamSubRate,m_MinContamLen,SeqLen,pSeq))==0)
		Rslt = m_pContaminates->MatchContaminants(eAOF3PE2Targ,m_MaxContamSubRate,m_MinContamLen,SeqLen,pSeq);
	}
return(Rslt);
}

int
CReadStats::ProcessReadsetDist(etRSDMode PMode,		// processing mode; eRSDindependent or eRSDpooled
				    int ProcessingID,				// processing instance identifier, used if processing eRSDindependent to identify output file instances 
					bool bStrand,					// true if read strand specific distributions
					int Trim5,						// trim this number of bases from 5' end of reads when loading the reads
					int Trim3,						// trim this number of bases from 3' end of reads when loading the reads
					int MaxKMerLen,					// processing is for upto this KMer length inclusive
					int KMerCCC,					// concordance correlation coefficient measure KMer length
					int MaxContamSubRate,			// max allowed contamimant substitution rate (bases per 25bp of contaminant overlap, 1st 15bp of overlap no subs allowed)
					int MinContamLen,				// accept contaminant overlaps if overlap at least this many bases 
					int ReqMaxDupSeeds,				// requested to sample for this many duplicate seeds
					int MinPhredScore,				// only accept reads for duplicate and KMer processing if their mean Phred score is at least this threshold
					int NumThreads,					// number of worker threads to use
					bool bAffinity,					// thread to core affinity
					int NumPE1InputFiles,			// number of PE1 input files
					char *pszInPE1files[],			// input PE1 5' read files
					int NumPE2InputFiles,			// number of PE2 input files
					char *pszInPE2files[],		    // input PE2 3' read files
					char *pszContaminantFile,		// contaminants fasta file
					char *pszOutDistFile,			// where to write distributions CSV file
					char *pszOutHTMLFile)			//  where to write distributions HTML5 file
{
int Rslt;
int Idx;
char *pszInFile;
int SeriesID;
UINT32 EstNumReads;
INT32 EstSeqLen;
INT32 EstMaxSeqLen;
UINT64 EstTotNumPE1Reads;
UINT64 EstTotNumPE2Reads;
UINT64 EstTotPE1SeqLen;
UINT64 EstTotPE2SeqLen;

bool bIsSAM;			// true if current file being processed is SAM or BAM
CFasta FastaEsts;		// used for estimating number of sequences, lengths
CSAMfile SAMFile;		// used if SAM or BAM and not fasta or fastq
int QSSchema;			// quality scoring schema - guestimated scoring schema - 0: no scoring, 1: Solexa, 2: Illumina 1.3+, 3: Illumina 1.5+, 4: Illumina 1.8+ or could be Sanger 

// make class instance copy of parameters
m_PMode = PMode;
m_bStrand = bStrand;
m_Trim5 = Trim5;
m_Trim3 = Trim3;
m_MaxKMerLen = MaxKMerLen;
m_KMerCCC = KMerCCC;
m_MaxContamSubRate = MaxContamSubRate;
m_MinContamLen = MinContamLen;
m_NumThreads = NumThreads;
m_bAffinity = bAffinity;
m_NumPE1InputFiles = NumPE1InputFiles;
m_ppszInPE1files = pszInPE1files;	
m_NumPE2InputFiles =NumPE2InputFiles;
m_ppszInPE2files = pszInPE2files;
m_pszContaminantFile = pszContaminantFile;
m_pszOutDistFile =pszOutDistFile;	
m_pszOutHTMLFile =pszOutHTMLFile;

// create/truncate any existing qscore/kmerdist/pearson report files
strcpy(m_szKMerDistRptFile,pszOutDistFile);
strcat(m_szKMerDistRptFile,".kmerdist.csv");
#ifdef _WIN32
if((m_hKMerDistRptFile = open(m_szKMerDistRptFile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
if((m_hKMerDistRptFile = open(m_szKMerDistRptFile,O_RDWR | O_CREAT |O_TRUNC, S_IREAD | S_IWRITE))==-1)
#endif
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessReadsetDist: (Instance %d) Unable to create or truncate %s - %s",ProcessingID,m_szKMerDistRptFile,strerror(errno));
	Reset();
	return(eBSFerrCreateFile);
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"(Instance %d) Output KMer distributions report file created/truncated: '%s'",ProcessingID,m_szKMerDistRptFile);

strcpy(m_szPearsonDistRptFile,pszOutDistFile);
strcat(m_szPearsonDistRptFile,".pearsondist.csv");
#ifdef _WIN32
if((m_hPearsonDistRptFile = open(m_szPearsonDistRptFile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
if((m_hPearsonDistRptFile = open(m_szPearsonDistRptFile,O_RDWR | O_CREAT |O_TRUNC, S_IREAD | S_IWRITE))==-1)
#endif
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessReadsetDist: (Instance %d) Unable to create or truncate %s - %s",ProcessingID,m_szPearsonDistRptFile,strerror(errno));
	Reset();
	return(eBSFerrCreateFile);
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"(Instance %d) Output Pearson KMer concordance distributions report file created/truncated: '%s'",ProcessingID,m_szPearsonDistRptFile);

strcpy(m_szQScoreDistRptFile,pszOutDistFile);
strcat(m_szQScoreDistRptFile,".qscoredist.csv");
#ifdef _WIN32
if((m_hQScoreDistRptFile = open(m_szQScoreDistRptFile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
if((m_hQScoreDistRptFile = open(m_szQScoreDistRptFile,O_RDWR | O_CREAT |O_TRUNC, S_IREAD | S_IWRITE))==-1)
#endif
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessReadsetDist: (Instance %d) Unable to create or truncate %s - %s",ProcessingID,m_szQScoreDistRptFile,strerror(errno));
	Reset();
	return(eBSFerrCreateFile);
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"(Instance %d) Output Phred quality score distributions report file created/truncated: '%s'",ProcessingID,m_szQScoreDistRptFile);

strcpy(m_szErrFreeReadDistRptFile,pszOutDistFile);
strcat(m_szErrFreeReadDistRptFile,".errfreedist.csv");
#ifdef _WIN32
if((m_hErrFreeReadDistRptFile = open(m_szErrFreeReadDistRptFile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
if((m_hErrFreeReadDistRptFile = open(m_szErrFreeReadDistRptFile,O_RDWR | O_CREAT |O_TRUNC, S_IREAD | S_IWRITE))==-1)
#endif
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessReadsetDist: (Instance %d) Unable to create or truncate %s - %s",ProcessingID,m_szErrFreeReadDistRptFile,strerror(errno));
	Reset();
	return(eBSFerrCreateFile);
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"(Instance %d) Output Phred derived error free probability distribution report file created/truncated: '%s'",ProcessingID,m_szQScoreDistRptFile);


strcpy(m_szDuplicatesDistRptFile,pszOutDistFile);
strcat(m_szDuplicatesDistRptFile,".duplicatesdist.csv");
#ifdef _WIN32
if((m_hDuplicatesDistRptFile = open(m_szDuplicatesDistRptFile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
if((m_hDuplicatesDistRptFile = open(m_szDuplicatesDistRptFile,O_RDWR | O_CREAT |O_TRUNC, S_IREAD | S_IWRITE))==-1)
#endif
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessReadsetDist: (Instance %d) Unable to create or truncate %s - %s",ProcessingID,m_szDuplicatesDistRptFile,strerror(errno));
	Reset();
	return(eBSFerrCreateFile);
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"(Instance %d) Output duplicate read distributions report file created/truncated: '%s'",ProcessingID,m_szDuplicatesDistRptFile);

strcpy(m_szReadLenDistRptFile,pszOutDistFile);
strcat(m_szReadLenDistRptFile,".readlendist.csv");
#ifdef _WIN32
if((m_hReadLenDistRptFile = open(m_szReadLenDistRptFile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
if((m_hReadLenDistRptFile = open(m_szReadLenDistRptFile,O_RDWR | O_CREAT |O_TRUNC, S_IREAD | S_IWRITE))==-1)
#endif
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessReadsetDist: (Instance %d) Unable to create or truncate %s - %s",ProcessingID,m_szReadLenDistRptFile,strerror(errno));
	Reset();
	return(eBSFerrCreateFile);
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"(Instance %d) Output read length distributions report file created/truncated: '%s'",ProcessingID,m_szReadLenDistRptFile);

// if putative contaminant processing required then create report file and load contaminants
if (pszContaminantFile != NULL && pszContaminantFile[0] != '\0')
	{
	strcpy(m_szContamRptFile,pszOutDistFile);
	strcat(m_szContamRptFile,".contaminants.csv");
#ifdef _WIN32
	if((m_hContamRptFile = open(m_szContamRptFile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
	if((m_hContamRptFile = open(m_szContamRptFile,O_RDWR | O_CREAT |O_TRUNC, S_IREAD | S_IWRITE))==-1)
#endif
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessReadsetDist: (Instance %d) Unable to create or truncate %s - %s",ProcessingID,m_szContamRptFile,strerror(errno));
		Reset();
		return(eBSFerrCreateFile);
		}

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"(Instance %d) Output contaminants report file created/truncated: '%s'",ProcessingID,m_szContamRptFile);
	if(m_pContaminates == NULL)
		{
		if((m_pContaminates = new CContaminants)==NULL)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "ProcessReadsetDist: (Instance %d) Unable to instantiate CContaminants",ProcessingID);
			Reset();
			return(eBSFerrObj);
			}
		}
	if ((Rslt=m_pContaminates->LoadContaminantsFile(pszContaminantFile)) < eBSFSuccess)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "ProcessReadsetDist: (Instance %d) Unable to load putative contaminants from file '%s", ProcessingID, pszContaminantFile);
		Reset();
		return(Rslt);
		}

	}

m_bPEProc = NumPE2InputFiles > 0 ? true : false;
if (m_bPEProc)
	ReqMaxDupSeeds = ((ReqMaxDupSeeds - 1) & ~0x01) + 2;	// round up to be even

m_ReqMaxDupSeeds = ReqMaxDupSeeds;
m_MinMeanPhredScore = MinPhredScore;

// single pass over input readsets is required
CSimpleGlob glob(SG_GLOB_FULLSORT);

	// estimate total number of reads in each file
EstTotNumPE1Reads = 0;
EstTotNumPE2Reads = 0;
EstTotPE1SeqLen = 0;
EstTotPE2SeqLen = 0;
EstMaxSeqLen = 0;
m_NumInFiles = 0;
m_NumPE1InFiles = 0;
m_NumPE2InFiles = 0;
m_EstMaxSeqLen = 0;
for(Idx = 0; Idx < NumPE1InputFiles; Idx++)
	{
	pszInFile = NULL;
	glob.Init();
	if(glob.Add(pszInPE1files[Idx]) < SG_SUCCESS)
		{
		gDiagnostics.DiagOut(eDLWarn,gszProcName,"ProcessReadsetDist: (Instance %d) Unable to glob '%s",ProcessingID,pszInPE1files[Idx]);
		return(eBSFerrOpnFile);	// treat as though unable to open file
		}

	if(glob.FileCount() <= 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessReadsetDist: (Instance %d) Unable to locate any source reads file matching '%s",ProcessingID,pszInPE1files[Idx]);
		continue;
		}
	Rslt = eBSFSuccess;
	for (int FileID = 0; Rslt >= eBSFSuccess &&  FileID < glob.FileCount(); ++FileID)
		{
		pszInFile = glob.File(FileID);

			// get estimate of number of sequences, mean sequence length, and mean descriptor length
		// check if a SAM or BAM file
		if(CSAMfile::IsSAM(pszInFile))
			{
			bIsSAM = true;
			if(m_bPEProc)		// can handle PE with ends in separate files
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessReadsetDist: (Instance %d) Unable to process paired end reads in multiple separate SAM or BAM files: '%s'",ProcessingID,pszInFile);
				return(eBSFerrOpnFile);
				}
			}
		else
			bIsSAM = false;

		if(bIsSAM)
			{
			if ((EstNumReads = SAMFile.EstSizes(pszInFile, NULL, NULL, NULL, &EstMaxSeqLen, &EstSeqLen,&QSSchema)) == 0)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessReadsetDist: (Instance %d) Unable to estimate number of reads in file '%s'",ProcessingID,pszInFile);
				return(eBSFerrOpnFile);
				}	 
			}
		else
			if ((EstNumReads = FastaEsts.FastaEstSizes(pszInFile, NULL, NULL, NULL, &EstMaxSeqLen, &EstSeqLen,&QSSchema)) == 0)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessReadsetDist: (Instance %d) Unable to estimate number of reads in file '%s'",ProcessingID,pszInFile);
				return(eBSFerrOpnFile);
				}	 
		memset(&m_InReadsFiles[m_NumInFiles], 0, sizeof(m_InReadsFiles[m_NumInFiles]));
		EstTotNumPE1Reads += EstNumReads;
		if(EstMaxSeqLen > m_EstMaxSeqLen)
			m_EstMaxSeqLen = EstMaxSeqLen;
		EstTotPE1SeqLen += (UINT64)EstSeqLen * EstNumReads;
		m_InReadsFiles[m_NumInFiles].FileID = m_NumInFiles + 1;
		strncpy(m_InReadsFiles[m_NumInFiles].szFileName, pszInFile, _MAX_PATH);
		m_InReadsFiles[m_NumInFiles].QSSchema = QSSchema;
		m_InReadsFiles[m_NumInFiles].szFileName[_MAX_PATH - 1] = '\0';
		m_InReadsFiles[m_NumInFiles].EstMeanReadLen = EstSeqLen;
		m_InReadsFiles[m_NumInFiles++].EstNumReads = EstNumReads;
		}
	}
if(EstTotNumPE1Reads < 1)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessReadsetDist: (Instance %d) Unable to estimate number of reads in file '%s'",ProcessingID,pszInFile == NULL ? "Any" : pszInFile);
	return(eBSFerrOpnFile);
	}
m_EstPE1MeanReadLen = (int)((EstTotPE1SeqLen + EstTotNumPE1Reads - 1) / EstTotNumPE1Reads);
m_NumPE1InFiles = m_NumInFiles;


if (m_bPEProc)
	{
	for(Idx = 0; Idx < NumPE2InputFiles; Idx++)
		{
		glob.Init();
		if(glob.Add(pszInPE2files[Idx]) < SG_SUCCESS)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessReadsetDist: Unable to glob '%s",pszInPE2files[Idx]);
			return(eBSFerrOpnFile);	// treat as though unable to open file
			}

		if(glob.FileCount() <= 0)
			{
			gDiagnostics.DiagOut(eDLWarn,gszProcName,"ProcessReadsetDist: (Instance %d) Unable to locate any source reads file matching '%s",ProcessingID,pszInPE2files[Idx]);
			continue;
			}
		Rslt = eBSFSuccess;
		for (int FileID = 0; Rslt >= eBSFSuccess &&  FileID < glob.FileCount(); ++FileID)
			{
			pszInFile = glob.File(FileID);
			// get estimate of number of sequences and mean sequence length
			if ((EstNumReads = FastaEsts.FastaEstSizes(pszInFile, NULL, NULL, NULL, &EstMaxSeqLen, &EstSeqLen, &QSSchema)) == 0)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessReadsetDist: (Instance %d) Unable to estimate number of reads in file '%s'",ProcessingID,pszInFile);
				return(eBSFerrOpnFile);
				}
			memset(&m_InReadsFiles[m_NumInFiles], 0, sizeof(m_InReadsFiles[m_NumInFiles]));
			EstTotNumPE2Reads += EstNumReads;
			if(EstMaxSeqLen > m_EstMaxSeqLen)
				m_EstMaxSeqLen = EstMaxSeqLen;
			EstTotPE2SeqLen += (UINT64)EstSeqLen * EstNumReads;
			m_InReadsFiles[m_NumInFiles].FileID = m_NumInFiles + 1;
			strncpy(m_InReadsFiles[m_NumInFiles].szFileName, pszInFile, _MAX_PATH);
			m_InReadsFiles[m_NumInFiles].QSSchema = QSSchema;
			m_InReadsFiles[m_NumInFiles].szFileName[_MAX_PATH - 1] = '\0';
			m_InReadsFiles[m_NumInFiles].EstMeanReadLen = EstSeqLen;
			m_InReadsFiles[m_NumInFiles++].EstNumReads = EstNumReads;
			m_NumPE2InFiles += 1;
			}
		}
	m_EstPE2MeanReadLen = (int)((EstTotPE2SeqLen + EstTotNumPE2Reads - 1) / EstTotNumPE2Reads);
	}
else
	m_EstPE2MeanReadLen = 0;

// double check that if PE processing then the number of PE1 files must match the number of PE2 files
if (m_bPEProc && (m_NumPE1InFiles != m_NumPE2InFiles))
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "ProcessReadsetDist: (Instance %d) Expected number of PE1 files (%d) to be same as PE2 files (%d) when PE readset processing", ProcessingID,m_NumPE1InFiles, m_NumPE2InFiles);
	return(eBSFerrOpnFile);
	}

// have work to do!
// preallocate for distribution counts
if((m_pBaseNs = new UINT32 [cMaxRSSeqLen])==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName," (Instance %d) ProcessReadsetDist: Unable to allocate memory for indeterminate base counts",ProcessingID);
	Reset();
	return(eBSFerrMem);
	}
memset(m_pBaseNs,0,sizeof(UINT32) * cMaxRSSeqLen);

if((m_pScores = new UINT32 [cMaxRSSeqLen * 42])==NULL)			// Phred scores can range from 0 to 41 inclusive
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"(Instance %d) ProcessReadsetDist: Unable to allocate memory for Phred quality score distributions",ProcessingID);
	Reset();
	return(eBSFerrMem);
	}
memset(m_pScores,0,sizeof(UINT32) * cMaxRSSeqLen * 42);

if ((m_pSeqHashes = new UINT32[cMaxHashArrayEntries]) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "ProcessReadsetDist: (Instance %d) Memory allocation of %lld bytes for hashes - %s", ProcessingID,(INT64)cMaxHashArrayEntries*sizeof(UINT32), strerror(errno));
	Reset();
	return(eBSFerrMem);
	}
memset(m_pSeqHashes, 0, sizeof(UINT32)*cMaxHashArrayEntries);
m_NumSeqHashes = 0;

size_t memreq;

m_AllocdMaxReadLen = (m_EstMaxSeqLen * 120)/100;	// allocate for longer than the estimated maximum read length to reduce the chances of having to later realloc
memreq = 0;
int Pow;
Pow = 1;
for(int Kmer = 1; Kmer <= m_MaxKMerLen; Kmer++)
	{
	Pow <<= 2;
	memreq += Pow;
	}
m_KMerCntsEls = (int)memreq;
memreq *= m_AllocdMaxReadLen * sizeof(UINT32);
 
#ifdef _WIN32
m_pKMerCnts = (UINT32 *)malloc(memreq);	// initial and with any luck perhaps the only allocation
if (m_pKMerCnts == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "ProcessReadsetDist: (Instance %d) Memory allocation of %lld bytes for K-mer distributions - %s", ProcessingID, (INT64)memreq, strerror(errno));
	Reset();
	return(eBSFerrMem);
	}
#else
	// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
m_pKMerCnts = (UINT32 *)mmap(NULL, memreq, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
if (m_pKMerCnts == MAP_FAILED)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "ProcessReadsetDist: (Instance %d) Memory allocation of %lld bytes  for K-mer distributions through mmap()  failed - %s", ProcessingID, (INT64)memreq, strerror(errno));
	m_pKMerCnts = NULL;
	Reset();
	return(eBSFerrMem);
	}
#endif
m_AllocdKMerCntsMem = memreq;
memset(m_pKMerCnts,0,memreq);

m_AllocdSampledSeqWrds = ((sizeof(tsSampledSeq)+3) / 4) + (((m_EstPE1MeanReadLen + m_EstPE2MeanReadLen + 5) + 15) / 16);	// 2bits per base and 16 bases per 32bit word, also allowing additional 5 bases as lengths are estimated	
m_AllocdSampledSeqWrds *= ReqMaxDupSeeds;
memreq = m_AllocdSampledSeqWrds * 4;

#ifdef _WIN32
m_pSampledSeqs = (UINT32 *)malloc(memreq);	// initial and with any luck perhaps the only allocation
if (m_pSampledSeqs == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "ProcessReadsetDist: (Instance %d) Memory allocation of %lld bytes - %s", ProcessingID, (INT64)memreq, strerror(errno));
	Reset();
	return(eBSFerrMem);
	}
#else
	// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
m_pSampledSeqs = (UINT32 *)mmap(NULL, memreq, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
if (m_pSampledSeqs == MAP_FAILED)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "ProcessReadsetDist: (Instance %d) Memory allocation of %lld bytes through mmap()  failed - %s",ProcessingID, (INT64)memreq, strerror(errno));
	m_pSampledSeqs = NULL;
	Reset();
	return(eBSFerrMem);
	}
#endif
m_AllocdSampledSeqMem = memreq;
m_UsedSampledSeqWrds = 0;
m_bTerminate = false;

// initialise thread contexts
CreateMutexes();

if (NumThreads > m_NumPE1InFiles)
	NumThreads = m_NumPE1InFiles;
tsThreadNGSQCPars *pThreads;
tsThreadNGSQCPars *pThread;
if ((pThreads = new tsThreadNGSQCPars[NumThreads]) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "ProcessReadsetDist: (Instance %d) Memory allocation for thread context failed",ProcessingID);
	return(eBSFerrMem);
	}
memset(pThreads,0,sizeof(tsThreadNGSQCPars) * NumThreads);
int ThreadIdx;
pThread = pThreads;
for (ThreadIdx = 1; ThreadIdx <= NumThreads; ThreadIdx++, pThread++)
	{
	pThread->ThreadIdx = ThreadIdx;
	pThread->pThis = this;
	pThread->ProcessingID = ProcessingID;

#ifdef _WIN32
	pThread->threadHandle = (HANDLE)_beginthreadex(NULL, 0x0fffff, ThreadedNGSQC, pThread, 0, &pThread->threadID);
#else
	pThread->threadRslt = pthread_create(&pThread->threadID, NULL, ThreadedNGSQC, pThread);
#endif
	}

// allow threads a few seconds to startup
#ifdef _WIN32
Sleep(3000);
#else
sleep(3);
#endif

pThread = pThreads;
for (ThreadIdx = 0; ThreadIdx < NumThreads; ThreadIdx++, pThread++)
	{
#ifdef _WIN32
	while (WAIT_TIMEOUT == WaitForSingleObject(pThread->threadHandle, 60000)) {	};
	CloseHandle(pThread->threadHandle);
#else
	struct timespec ts;
	int JoinRlt;
	clock_gettime(CLOCK_REALTIME, &ts);
	ts.tv_sec += 60;
	while ((JoinRlt = pthread_timedjoin_np(pThread->threadID, NULL, &ts)) != 0)
		{
		ts.tv_sec += 60;
		}
#endif
	}

DeleteMutexes();

if (m_bTerminate)		// early termination because of some problem?
	{
	delete [] pThreads;
	Reset();
	return(-1);
	}


// report on the number of reads processed
pThread = pThreads;
INT64 TotNumSEReads;
INT64 TotNumPEReads;
INT64 NotProcNs;
INT64 NotProcQS;
INT64 NotProcUL;


TotNumSEReads = 0;
TotNumPEReads = 0;
NotProcNs = 0;
NotProcQS = 0;
NotProcUL = 0;

for (ThreadIdx = 0; ThreadIdx < NumThreads; ThreadIdx++, pThread++)
	{
	TotNumSEReads += pThread->TotNumSEReads;
	TotNumPEReads += pThread->TotNumPEReads;

	NotProcNs += pThread->SeqCharacteristics.NotProcNs;
	NotProcQS += pThread->SeqCharacteristics.NotProcQS;
	NotProcUL += pThread->SeqCharacteristics.NotProcUL;
	}

if (!m_bPEProc)
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "(Instance %d) Total of %llu SE reads processed", ProcessingID, TotNumSEReads);
else
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "(Instance %d) Total of %llu reads, %llu PE pairs processed", ProcessingID, TotNumSEReads, TotNumPEReads);

gDiagnostics.DiagOut(eDLInfo, gszProcName, "(Instance %d) Total of %llu reads not accepted for processing as they were underlength", ProcessingID, NotProcUL);

gDiagnostics.DiagOut(eDLInfo, gszProcName, "(Instance %d) Total of %lu reads used as seed duplicates", ProcessingID, m_ActMaxDupSeeds);
gDiagnostics.DiagOut(eDLInfo, gszProcName, "(Instance %d) Total of %llu reads containing 'N's and/or %llu below minimum Phred scores not accepted for duplicate processing", ProcessingID, NotProcNs,NotProcQS);

if(m_ActMaxDupSeeds == 0)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName,"(Instance %d) Unable to accept any reads for duplicate and K-mer analysis, you may need to check on 'N's or lower the Phred score acceptance threshold",ProcessingID);
	Reset();
	return(-1);
	}

// report on duplicate counts for sampled reads
UINT32 DupDists[2001];			// duplicate counts limited 1..1999 and 2000+
UINT32 DupDist10[11];			// duplicate counts limited 1..9 and 10+

UINT32 TotalDupReads;			// total number of reads in samples
	
memset(DupDists,0,sizeof(DupDists));
memset(DupDist10, 0, sizeof(DupDist10));
TotalDupReads = 0;
if (m_pSampledSeqs != NULL)
	{
	int HashIdx;
	UINT32 HashSeqsOfs;
	tsSampledSeq *pSampledSeq;
	for (HashIdx = 0; HashIdx < cMaxHashArrayEntries; HashIdx++)
		{
		if ((HashSeqsOfs = m_pSeqHashes[HashIdx]) != 0)		// seen at least one sequence with this hash?
			{
			pSampledSeq = (tsSampledSeq *)&m_pSampledSeqs[HashSeqsOfs - 1];
			do
				{
				if (pSampledSeq->NumInstances > 2000)
					DupDists[1999] += pSampledSeq->NumInstances;
				else
					DupDists[pSampledSeq->NumInstances-1] += pSampledSeq->NumInstances;

				if (pSampledSeq->NumInstances > 10)
					DupDist10[9] += pSampledSeq->NumInstances;
				else
					DupDist10[pSampledSeq->NumInstances - 1] += pSampledSeq->NumInstances;
				TotalDupReads += pSampledSeq->NumInstances;
				}
			while (pSampledSeq->NxtSeq != 0 && (pSampledSeq = (tsSampledSeq *)&m_pSampledSeqs[pSampledSeq->NxtSeq - 1]));
			}
		}
	}

// log 1..10 duplicate instance counts
for (Idx = 0; Idx < 10; Idx++)
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "(Instance %d)  %d: %2.2f", ProcessingID, Idx + 1, (100.0 * DupDist10[Idx]) / TotalDupReads);

gDiagnostics.DiagOut(eDLInfo, gszProcName, "(Instance %d) Contaminants PE1 Checked: %u Contaminated: %u Percentage: %1.4f", ProcessingID,m_NumChkdPE1ContamHits,m_NumPE1ContamHits,m_NumChkdPE1ContamHits ? (m_NumPE1ContamHits*100.0)/m_NumChkdPE1ContamHits : 0.0);
gDiagnostics.DiagOut(eDLInfo, gszProcName, "(Instance %d) Contaminants PE2 Checked: %u Contaminated: %u Percentage: %1.4f", ProcessingID,m_NumChkdPE2ContamHits,m_NumPE2ContamHits,m_NumChkdPE2ContamHits ? (m_NumPE2ContamHits*100.0)/m_NumChkdPE2ContamHits : 0.0);

CBKPLPlot *pPlots;
if((pPlots = new CBKPLPlot)==NULL)
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "(Instance %d) Unable to instantiate CBKPlot");


// report actual duplicate counts for downstream analytics

if(pPlots != NULL)
	{
	char szInsts[10];
	pPlots->Init(1,1,10,(char *)"Duplicate Reads Distribution",(char *)"Instances",(char *)"Proportion");
	pPlots->SetWorldCoords(0.0,(PLFLT)10,0.0,1.0);
	SeriesID = pPlots->AddSeries((char *)"Duplicates",true,BKPLwheat);
	pPlots->SetAxis(0,1.0,0,(char *)"abgit");
	pPlots->SetAxis(1,0.1,5,(char *)"abignstv");
	strcpy(m_szSVGFile,pszOutDistFile);
	strcat(m_szSVGFile,".duplicatesdist.svg");
	for(Idx = 1; Idx <= 10; Idx++)
		{
		if(Idx < 10)
			sprintf(szInsts,"%d",Idx);
		else
			sprintf(szInsts,"%d plus",Idx);
		pPlots->AddPoint(SeriesID,DupDists[Idx-1]/(PLFLT)TotalDupReads,szInsts);
		}
	pPlots->PlotBarChartGraph(m_szSVGFile);
	}

char szRptBuff[10000];
int BuffIdx;
BuffIdx = sprintf(szRptBuff,"\"Instances\",\"Count\"\n");
for (Idx = 0; Idx < 1999; Idx++)
	{
	BuffIdx += sprintf(&szRptBuff[BuffIdx],"\"%d\",%u\n",Idx+1,DupDists[Idx]);
	if(BuffIdx + 100 > sizeof(szRptBuff))
		{
		CUtility::SafeWrite(m_hDuplicatesDistRptFile,szRptBuff,BuffIdx);
		BuffIdx = 0;
		}
	}
BuffIdx += sprintf(&szRptBuff[BuffIdx],"\"%d+\",%u\n",Idx+1,DupDists[Idx]);
if(BuffIdx)
	CUtility::SafeWrite(m_hDuplicatesDistRptFile,szRptBuff,BuffIdx);
#ifdef _WIN32
_commit(m_hDuplicatesDistRptFile);
#else
fsync(m_hDuplicatesDistRptFile);
#endif
close(m_hDuplicatesDistRptFile);
m_hDuplicatesDistRptFile = -1;


// report read length distributions for downstream analytics
if(pPlots != NULL)
	{
	char szLen[10];
	UINT64 SumCnts; 
	UINT64 CummulativeSum;
	int LenRange = 1 + m_MaxReadLen - m_MinReadLen;
	pPlots->Init(1,1,LenRange+10,(char *)"Read Length Distribution",(char *)"Length (bp)",(char *)"Proportion of Reads");
	pPlots->SetWorldCoords((PLFLT)max(1,m_MinReadLen - 5),(PLFLT)m_MaxReadLen+5,0.0,1.05);
	SeriesID = pPlots->AddSeries((char *)"Lengths",false,BKPLwheat);
	pPlots->SetAxis(0,5.0,5,(char *)"abgist");
	pPlots->SetAxis(1,0.1,5,(char *)"abignstv");
	strcpy(m_szSVGFile,pszOutDistFile);
	strcat(m_szSVGFile,".readlendist.svg");
	SumCnts = 0;
	CummulativeSum = 0;
	for (int Idx = m_MinReadLen; Idx <= m_MaxReadLen; Idx++)
		CummulativeSum += m_ReadLenDist[Idx];

	for (int Idx = max(1,m_MinReadLen - 5); Idx <= m_MaxReadLen+5; Idx++)
			{
			sprintf(szLen,"%d",Idx);
			pPlots->AddPoint(SeriesID,m_ReadLenDist[Idx]/(PLFLT)CummulativeSum,szLen);
			}
	pPlots->PlotBarChartGraph(m_szSVGFile);
	}

BuffIdx = sprintf(szRptBuff,"\"ReadLength\",\"Count\"\n");
for (Idx = 1; Idx <= m_MaxReadLen; Idx++)
	{
	BuffIdx += sprintf(&szRptBuff[BuffIdx],"\"%d\",%u\n",Idx,m_ReadLenDist[Idx]);
	if(BuffIdx + 100 > sizeof(szRptBuff))
		{
		CUtility::SafeWrite(m_hReadLenDistRptFile,szRptBuff,BuffIdx);
		BuffIdx = 0;
		}
	}
BuffIdx += sprintf(&szRptBuff[BuffIdx],"\"%d+\",%u\n",Idx,m_ReadLenDist[Idx]);
if(BuffIdx)
	CUtility::SafeWrite(m_hReadLenDistRptFile,szRptBuff,BuffIdx);
#ifdef _WIN32
_commit(m_hReadLenDistRptFile);
#else
fsync(m_hReadLenDistRptFile);
#endif
close(m_hReadLenDistRptFile);
m_hReadLenDistRptFile = -1;

// report Phred quality scores for downstream analytics
if(pPlots != NULL)
	{
	char szPhred[10];
	UINT64 SumCnts; 
	UINT64 CummulativeSum;
	int PhredIDs[42];
	int PhredIdx;
	int Idy;
	pPlots->Init(1,42,m_MaxReadLen,(char *)"Phred Quality Score Distribution",(char *)"Read Base Offset",(char *)"Phred Score");
	pPlots->SetWorldCoords(0.0,(PLFLT)m_MaxReadLen,0.0,(PLFLT)42.05);
	for(PhredIdx = 0; PhredIdx < 42; PhredIdx++)
		{
		sprintf(szPhred,"Phred: %d",PhredIdx);
		PhredIDs[PhredIdx] = pPlots->AddSeries(szPhred,false,BKPLwheat);
		}
	pPlots->SetAxis(0,10.0,5,(char *)"abgistn");
	pPlots->SetAxis(1,5.0,5,(char *)"abignstv");
	strcpy(m_szSVGFile,pszOutDistFile);
	strcat(m_szSVGFile,".qscoredist.svg");
	SumCnts = 0;
	CummulativeSum = 0;

	for (Idx = 0; Idx < m_MaxReadLen; Idx++)
		{
		CummulativeSum = 0;
		for (Idy = 0; Idy < 42; Idy++)
			CummulativeSum += m_pScores[Idy + (Idx * 42)];
		for (Idy = 0; Idy < 42; Idy++)
			pPlots->AddPoint(PhredIDs[Idy],m_pScores[Idy + (Idx * 42)]/(double)CummulativeSum);
		}
	pPlots->PlotPhredScores(m_szSVGFile);
	}

BuffIdx = sprintf(szRptBuff,"\"Phred\"");
for (Idx = 1; Idx <= m_MaxReadLen; Idx++)
	{
	BuffIdx += sprintf(&szRptBuff[BuffIdx],",\"Base:%d\"",Idx);
	if(BuffIdx + 100 > sizeof(szRptBuff))
		{
		CUtility::SafeWrite(m_hQScoreDistRptFile,szRptBuff,BuffIdx);
		BuffIdx = 0;
		}
	}
if(BuffIdx > 0)
	{
	CUtility::SafeWrite(m_hQScoreDistRptFile,szRptBuff,BuffIdx);
	BuffIdx = 0;
	}
for (int Idy = 0; Idy < 42; Idy++)
	{
	BuffIdx += sprintf(&szRptBuff[BuffIdx],"\n%d",Idy);
	for (Idx = 0; Idx < m_MaxReadLen; Idx++)
		{
		BuffIdx += sprintf(&szRptBuff[BuffIdx],",%d",m_pScores[Idy + (Idx * 42)]);
		if(BuffIdx + 100 > sizeof(szRptBuff))
			{
			CUtility::SafeWrite(m_hQScoreDistRptFile,szRptBuff,BuffIdx);
			BuffIdx = 0;
			}
		}
	}
if(BuffIdx)
	CUtility::SafeWrite(m_hQScoreDistRptFile,szRptBuff,BuffIdx);
#ifdef _WIN32
_commit(m_hQScoreDistRptFile);
#else
fsync(m_hQScoreDistRptFile);
#endif
close(m_hQScoreDistRptFile);
m_hQScoreDistRptFile = -1;


if(pPlots != NULL)
	{
	UINT64 SumCnts; 
	UINT64 CummulativeSum;
	int CumulativeSeriesID;
	pPlots->Init(2,2,100,(char *)"Phred Derived Error Free Probability Distribution",(char *)"Error Free Probability",(char *)"Proportion of Reads");
	pPlots->SetWorldCoords(0.0,(PLFLT)1.01,0.0,1.05);
	SeriesID = pPlots->AddSeries((char *)"Discrete",false,BKPLwheat);
	CumulativeSeriesID = pPlots->AddSeries((char *)"Cumulative",false,BKPLblue);
	pPlots->SetAxis(0,0.1,5,(char *)"abgistn");
	pPlots->SetAxis(1,0.1,5,(char *)"abignstv");
	strcpy(m_szSVGFile,pszOutDistFile);
	strcat(m_szSVGFile,".errfreedist.svg");
	SumCnts = 0;
	CummulativeSum = 0;
	for (int Idx = 0; Idx < 100; Idx++)
		SumCnts += m_ProbNoReadErrDist[Idx];
	if(SumCnts == 0)
		SumCnts = 1;
	for (int Idx = 0; Idx < 100; Idx++)
		{
		pPlots->AddXYPoint(SeriesID,(double)Idx/100,m_ProbNoReadErrDist[Idx]/(PLFLT)SumCnts);
		CummulativeSum += m_ProbNoReadErrDist[Idx];
		pPlots->AddXYPoint(CumulativeSeriesID,(double)Idx/100,1.0 - (CummulativeSum/(PLFLT)SumCnts));
		}
	pPlots->PlotLineGraph(m_szSVGFile);
	}

// report error free read probabilities distributions for downstream analytics
BuffIdx = sprintf(szRptBuff,"\"ProbabilityBin\",\"Count\"");
CUtility::SafeWrite(m_hErrFreeReadDistRptFile,szRptBuff,BuffIdx);
BuffIdx = 0;
for (int Idx = 0; Idx < 100; Idx++)
	BuffIdx += sprintf(&szRptBuff[BuffIdx],"\n%1.2f,%lld",(double)Idx*0.01,m_ProbNoReadErrDist[Idx]);
if(BuffIdx)
	CUtility::SafeWrite(m_hErrFreeReadDistRptFile,szRptBuff,BuffIdx);
#ifdef _WIN32
_commit(m_hErrFreeReadDistRptFile);
#else
fsync(m_hErrFreeReadDistRptFile);
#endif
close(m_hErrFreeReadDistRptFile);
m_hErrFreeReadDistRptFile = -1;


// report on the Kmer distributions, from 1 upto m_MaxKMerLen
int KMerLen;
int KMerOfs;
int NumEls;
int BaseIdx;
UINT32 *pCtrlKMerCnts;
UINT32 *pExprKMerCnts;

if(pPlots != NULL)
	{
	UINT64 SumCnts; 
	UINT32 *pKMerCnts;
	pPlots->Init(2,4,m_MaxReadLen,(char *)"Compositional Distribution",(char *)"Read Base Offset",(char *)"Proportion");
	pPlots->SetWorldCoords(0.0,(PLFLT)(PLFLT)m_MaxReadLen+0.05,0.0,1.05);
	pPlots->AddSeries((char *)"A",false,BKPLaquamarine);
	pPlots->AddSeries((char *)"C",false,BKPLcyan);
	pPlots->AddSeries((char *)"G",false,BKPLgreen);
	pPlots->AddSeries((char *)"T",false,BKPLturquoise);
	pPlots->SetAxis(0,10.0,5,(char *)"abgistn");
	pPlots->SetAxis(1,0.1,5,(char *)"abignstv");
	strcpy(m_szSVGFile,pszOutDistFile);
	strcat(m_szSVGFile,".acgtdist.svg");
	pCtrlKMerCnts = m_pKMerCnts;
	for(Idx = 0; Idx < m_MaxReadLen; Idx++,pCtrlKMerCnts += m_KMerCntsEls)
		{
		pKMerCnts = pCtrlKMerCnts;
		SumCnts = 0;
		for(int Idy=0;Idy < 4; Idy++,pKMerCnts++)
			SumCnts += (UINT64)*pKMerCnts;
		if(SumCnts == 0)
			SumCnts = 1;
		pKMerCnts = pCtrlKMerCnts;
		for(int Idy=0;Idy < 4; Idy++,pKMerCnts++)
			pPlots->AddXYPoint(Idy+1,(double)Idx,*pKMerCnts/(double)SumCnts);
		}
	if(pPlots != NULL)
		pPlots->PlotLineGraph(m_szSVGFile);
	}

KMerOfs = 0;
NumEls = 4;
pCtrlKMerCnts = m_pKMerCnts;
BuffIdx = 0;
for(KMerLen = 1; KMerLen <= m_MaxKMerLen; KMerLen++)
	{
	for(int Idy=0;Idy < NumEls; Idy++,pCtrlKMerCnts++)
		{
		BuffIdx += sprintf(&szRptBuff[BuffIdx],"\"");
		for(int Shr = KMerLen-1; Shr >= 0; Shr--)
			{
			BaseIdx = Idy >> (Shr*2);
			BaseIdx &= 0x03;
			switch(BaseIdx) {
				case 0:
					szRptBuff[BuffIdx++] = 'a';
					break;
				case 1:
					szRptBuff[BuffIdx++] = 'c';
					break;
				case 2:
					szRptBuff[BuffIdx++] = 'g';
					break;
				case 3:
					szRptBuff[BuffIdx++] = 't';
					break;
				}			
			}
		BuffIdx += sprintf(&szRptBuff[BuffIdx],"\"");
		pExprKMerCnts = pCtrlKMerCnts;
		for (Idx = 1; Idx <= m_MaxReadLen; Idx++, pExprKMerCnts += m_KMerCntsEls)
			{
			if((m_MaxReadLen - Idx) < (KMerLen-1))
				BuffIdx += sprintf(&szRptBuff[BuffIdx],",%u",*pExprKMerCnts);
			else
				BuffIdx += sprintf(&szRptBuff[BuffIdx],",%u",*pExprKMerCnts);
			if(BuffIdx + 100 > sizeof(szRptBuff))
				{
				CUtility::SafeWrite(m_hKMerDistRptFile,szRptBuff,BuffIdx);
				BuffIdx = 0;
				}
			}
		BuffIdx += sprintf(&szRptBuff[BuffIdx],"\n");
		}
	KMerOfs += NumEls;
	NumEls <<= 2;
	}
if(BuffIdx)
	CUtility::SafeWrite(m_hKMerDistRptFile,szRptBuff,BuffIdx);
#ifdef _WIN32
_commit(m_hKMerDistRptFile);
#else
fsync(m_hKMerDistRptFile);
#endif
close(m_hKMerDistRptFile);
m_hKMerDistRptFile = -1;


// report on the Pearson concordance distributions
if(pPlots != NULL)
	{
	char szSeries[20];
    char szPearsonTitle[80];
	sprintf(szPearsonTitle,"Pearson K-mer (1..%d) Concordance Distribution",m_KMerCCC);
	pPlots->Init(2,m_KMerCCC,m_MaxReadLen,szPearsonTitle,(char *)"Read Base Offset",(char *)"Pearson");
	pPlots->SetWorldCoords(0.0,(PLFLT)m_MaxReadLen,-1.0,1.0);
	for(Idx = 0; Idx < m_KMerCCC; Idx++)
		{
		sprintf(szSeries,"K-mer %d",Idx+1);
		SeriesID = pPlots->AddSeries(szSeries,false,(etBKPlotColor)((int)BKPLbrown + Idx));
		}
	pPlots->SetAxis(0,10.0,5,(char *)"abgistn");
	pPlots->SetAxis(1,0.2,5,(char *)"abignstv");
	strcpy(m_szSVGFile,pszOutDistFile);
	strcat(m_szSVGFile,".pearsondist.svg");
	}
double Pearson;
int ControlBase;
BuffIdx = sprintf(szRptBuff,"\"Pearson K-mer\"");
for(Idx = 0; Idx < m_MaxReadLen; Idx++)
	{
	BuffIdx += sprintf(&szRptBuff[BuffIdx],",\"Base:%d\"",Idx);
	if(BuffIdx + 100 > sizeof(szRptBuff))
		{
		CUtility::SafeWrite(m_hPearsonDistRptFile,szRptBuff,BuffIdx);
		BuffIdx = 0;
		}
	}
ControlBase = m_MaxReadLen/3;
for(int Idk = 1; Idk <= m_KMerCCC; Idk++)
	{
	BuffIdx += sprintf(&szRptBuff[BuffIdx],"\nK-mer: %d",Idk);
	pExprKMerCnts = m_pKMerCnts;
	pCtrlKMerCnts = &m_pKMerCnts[ControlBase * m_KMerCntsEls];
	KMerOfs = 0;
	NumEls = 4;
	for(Idx = 1; Idx < Idk; Idx++)
		{
		KMerOfs += NumEls;
		NumEls <<= 2;
		}
	pExprKMerCnts += KMerOfs;
	pCtrlKMerCnts += KMerOfs;
	for (Idx = 1; Idx <= m_MaxReadLen; Idx++,pExprKMerCnts+=m_KMerCntsEls)
		{
		Pearson = Pearsons(Idk,pCtrlKMerCnts,pExprKMerCnts);
		if(pPlots != NULL)
			pPlots->AddXYPoint(Idk,(double)Idx,Pearson);	
		BuffIdx += sprintf(&szRptBuff[BuffIdx],",%f",Pearson);
		if(BuffIdx + 100 > sizeof(szRptBuff))
			{
			CUtility::SafeWrite(m_hPearsonDistRptFile,szRptBuff,BuffIdx);
			BuffIdx = 0;
			}
		}
	}
if(BuffIdx)
	CUtility::SafeWrite(m_hPearsonDistRptFile,szRptBuff,BuffIdx);
#ifdef _WIN32
_commit(m_hPearsonDistRptFile);
#else
fsync(m_hPearsonDistRptFile);
#endif
close(m_hPearsonDistRptFile);
m_hPearsonDistRptFile = -1;
if(pPlots != NULL)
	pPlots->PlotLineGraph(m_szSVGFile);

if(m_hContamRptFile != -1)		// reporting on contaminants?
	{
	int Dist[cMaxContaminantLen+1];
	int ContamHits;
	int DistIdx;
	int NumContaminants;
	int NumFlankContaminants;
	int NumVectContaminants;
	int MaxContamLen;
	int ContamID;
	bool bPE1Contam;

	teContamType ContamType;
	teContamClass ContamClass;

	UINT32 NumChecks;
	int ContamLen;
	char *pszContamName;
	NumContaminants = m_pContaminates->NumOfContaminants();
	NumFlankContaminants = m_pContaminates->NumOfContaminants(eCCFlankContam);
	NumVectContaminants = m_pContaminates->NumOfContaminants(eCCVectContam);

	if(NumFlankContaminants > 0)
		{
		MaxContamLen = m_pContaminates->MaxContaminantLen(eCCFlankContam);
		BuffIdx = sprintf(szRptBuff,"\"Name\",\"Type\",\"Length\",\"Checked\",\"TotOverlaps\"");
		for(DistIdx = 1; DistIdx <= MaxContamLen; DistIdx++)
			BuffIdx += sprintf(&szRptBuff[BuffIdx],",\"Overlap:%d\"",DistIdx);
		BuffIdx += sprintf(&szRptBuff[BuffIdx],"\n");
		CUtility::SafeWrite(m_hContamRptFile,szRptBuff,BuffIdx);
		BuffIdx = 0;

		for(ContamID = 1; ContamID <= NumContaminants; ContamID++)
			{
			ContamClass = m_pContaminates->ContaminantClass(ContamID);
			if(ContamClass != eCCFlankContam)
				continue;
			pszContamName = m_pContaminates->ContaminantName(ContamID);
			ContamType = m_pContaminates->ContaminantType(ContamID);
		
			if(ContamType == eAOF5PE1Targ || ContamType == eAOF5PE2Targ)
				bPE1Contam = true;
			else
				bPE1Contam = false;
			ContamLen =  m_pContaminates->ContaminantLen(ContamID);
			NumChecks = m_pContaminates->NumChecks(ContamType);

			ContamHits = m_pContaminates->ContaminantDist(ContamID,cMaxContaminantLen,Dist);
			for(DistIdx = ContamLen-1; DistIdx > 0; DistIdx--)
				Dist[DistIdx-1] += Dist[DistIdx];

			BuffIdx = sprintf(szRptBuff,"\"%s\",\"%s\",%d,%u,%d",pszContamName,m_pContaminates->ContaminateType2Txt(ContamType),ContamLen,NumChecks,ContamHits);
			for(DistIdx = 0; DistIdx < ContamLen; DistIdx++)
				BuffIdx += sprintf(&szRptBuff[BuffIdx],",%d",Dist[DistIdx]);
			for(; DistIdx < MaxContamLen; DistIdx++)
				BuffIdx += sprintf(&szRptBuff[BuffIdx],",0");
			BuffIdx += sprintf(&szRptBuff[BuffIdx],"\n");
			CUtility::SafeWrite(m_hContamRptFile,szRptBuff,BuffIdx);
			BuffIdx = 0;
			}
		}

	if(NumVectContaminants > 0)
		{
		MaxContamLen = m_pContaminates->MaxContaminantLen(eCCVectContam);
		BuffIdx = sprintf(szRptBuff,"\"Name\",\"Type\",\"Length\",\"Checked\",\"TotContained\"\n");
		CUtility::SafeWrite(m_hContamRptFile,szRptBuff,BuffIdx);
		BuffIdx = 0;
		for(ContamID = 1; ContamID <= NumContaminants; ContamID++)
			{
			ContamClass = m_pContaminates->ContaminantClass(ContamID);
			if(ContamClass != eCCVectContam)
				continue;
			pszContamName = m_pContaminates->ContaminantName(ContamID);
			ContamType = m_pContaminates->ContaminantType(ContamID);
		
			ContamLen =  m_pContaminates->ContaminantLen(ContamID);
			NumChecks = m_pContaminates->NumChecks(ContamType);
			ContamHits = m_pContaminates->ContaminantDist(ContamID,cMaxContaminantLen,Dist);

			BuffIdx = sprintf(szRptBuff,"\"%s\",\"%s\",%d,%u,%d\n",pszContamName,m_pContaminates->ContaminateType2Txt(ContamType),ContamLen,NumChecks,ContamHits);
			CUtility::SafeWrite(m_hContamRptFile,szRptBuff,BuffIdx);
			BuffIdx = 0;
			}
		}

#ifdef _WIN32
	_commit(m_hContamRptFile);
#else
	fsync(m_hContamRptFile);
#endif
	close(m_hContamRptFile);
	m_hContamRptFile = -1;
	}

if(pPlots != NULL)
	delete pPlots;
return(eBSFSuccess);
}

int
CReadStats::CreateMutexes(void)
{
if (m_bMutexesCreated)
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
if (!InitializeCriticalSectionAndSpinCount(&m_hSCritSect, 1000))
	{
#else
if(pthread_spin_init(&m_hSpinLock,PTHREAD_PROCESS_PRIVATE)!=0)
	{
	pthread_rwlock_destroy(&m_hRwLock);
#endif
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Fatal: unable to create mutex");
	return(eBSFerrInternal);
	}


#ifdef _WIN32
if (!InitializeCriticalSectionAndSpinCount(&m_hSCritSectScores, 1000))
	{
#else
if(pthread_spin_init(&m_hSpinLockScores,PTHREAD_PROCESS_PRIVATE)!=0)
	{
	pthread_rwlock_destroy(&m_hRwLock);
#endif
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Fatal: unable to create mutex");
	return(eBSFerrInternal);
	}

#ifdef _WIN32
if (!InitializeCriticalSectionAndSpinCount(&m_hSCritSectKMers, 1000))
	{
#else
if(pthread_spin_init(&m_hSpinLockKMers,PTHREAD_PROCESS_PRIVATE)!=0)
	{
	pthread_rwlock_destroy(&m_hRwLock);
#endif
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Fatal: unable to create mutex");
	return(eBSFerrInternal);
	}

#ifdef _WIN32
if ((m_hMtxMHReads = CreateMutex(NULL, false, NULL)) == NULL)
	{
#else
if(pthread_mutex_init (&m_hMtxMHReads,NULL)!=0)
	{
#endif
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: unable to create mutex");
#ifndef _WIN32
	pthread_rwlock_destroy(&m_hRwLock);
#endif
	return(eBSFerrInternal);
	}

m_bMutexesCreated = true;
return(eBSFSuccess);
}

void
CReadStats::DeleteMutexes(void)
{
	if (!m_bMutexesCreated)
		return;
#ifdef _WIN32
	CloseHandle(m_hMtxMHReads);

#else
	pthread_mutex_destroy(&m_hMtxMHReads);
	pthread_rwlock_destroy(&m_hRwLock);

#endif
	m_bMutexesCreated = false;
}

void
CReadStats::AcquireSerialise(void)
{
	int SpinCnt = 1000;
#ifdef _WIN32
	while(!TryEnterCriticalSection(&m_hSCritSect))
	{
		if (SpinCnt -= 1)
			continue;
		SwitchToThread();
		SpinCnt = 100;
	}
#else
	while (pthread_spin_trylock(&m_hSpinLock) == EBUSY)
	{
		if (SpinCnt -= 1)
			continue;
		pthread_yield();
		SpinCnt = 100;
	}
#endif
}

void
CReadStats::ReleaseSerialise(void)
{
#ifdef _WIN32
	LeaveCriticalSection(&m_hSCritSect);
#else
	pthread_spin_unlock(&m_hSpinLock);
#endif
}


void
CReadStats::AcquireSerialiseScores(void)
{
int SpinCnt = 1000;
#ifdef _WIN32
while(!TryEnterCriticalSection(&m_hSCritSectScores))
	{
	if (SpinCnt -= 1)
		continue;
	SwitchToThread();
	SpinCnt = 100;
	}
#else
while (pthread_spin_trylock(&m_hSpinLockScores) == EBUSY)
	{
	if (SpinCnt -= 1)
		continue;
	pthread_yield();
	SpinCnt = 100;
	}
#endif
}

void
CReadStats::ReleaseSerialiseScores(void)
{
#ifdef _WIN32
	LeaveCriticalSection(&m_hSCritSectScores);
#else
	pthread_spin_unlock(&m_hSpinLockScores);
#endif
}

void
CReadStats::AcquireSerialiseKMers(void)
{
int SpinCnt = 1000;
#ifdef _WIN32
while(!TryEnterCriticalSection(&m_hSCritSectKMers))
	{
	if (SpinCnt -= 1)
		continue;
	SwitchToThread();
	SpinCnt = 100;
	}
#else
while (pthread_spin_trylock(&m_hSpinLockKMers) == EBUSY)
	{
	if (SpinCnt -= 1)
		continue;
	pthread_yield();
	SpinCnt = 100;
	}
#endif
}

void
CReadStats::ReleaseSerialiseKMers(void)
{
#ifdef _WIN32
	LeaveCriticalSection(&m_hSCritSectKMers);
#else
	pthread_spin_unlock(&m_hSpinLockKMers);
#endif
}

void
inline CReadStats::AcquireLock(bool bExclusive)
{
#ifdef _WIN32
	if (bExclusive)
		AcquireSRWLockExclusive(&m_hRwLock);
	else
		AcquireSRWLockShared(&m_hRwLock);
#else
	if (bExclusive)
		pthread_rwlock_wrlock(&m_hRwLock);
	else
		pthread_rwlock_rdlock(&m_hRwLock);
#endif
}

void
inline CReadStats::ReleaseLock(bool bExclusive)
{

#ifdef _WIN32
	if (bExclusive)
		ReleaseSRWLockExclusive(&m_hRwLock);
	else
		ReleaseSRWLockShared(&m_hRwLock);
#else
	pthread_rwlock_unlock(&m_hRwLock);
#endif
}


int
CReadStats::ProcNGSQC(tsThreadNGSQCPars *pThread)
{
int Rslt;
int Idx;
char *pszInFile;
int NumReads;
UINT32 NumPE1Reads;
UINT32 NumPE2Reads;
tsInReadsFile *pInPE1File;
tsInReadsFile *pInPE2File;
NumReads = 0;
NumPE1Reads = 0;
NumPE2Reads = 0;

pThread->NumInputFilesProcessed = 0;
pThread->TotNumSEReads = 0;
pThread->TotNumPEReads = 0;


if (!m_bPEProc)								// if SE then process each input reads file separately
	{
	pInPE1File = &m_InReadsFiles[0];
	for (Idx = 0; Idx < m_NumPE1InFiles; Idx++, pInPE1File++)
		{
		AcquireLock(true);
		if (m_bTerminate)		// check if requested to early exit
			{
			ReleaseLock(true);
			gDiagnostics.DiagOut(eDLInfo, gszProcName, "(Instance %d) Thread %d: Requested to early terminate", pThread->ProcessingID, pThread->ThreadIdx);
			return(eBSFerrInternal);
			}
		if (pInPE1File->ProcByThreadIdx != 0)
			{
			ReleaseLock(true);
			continue;
			}
		pInPE1File->ProcByThreadIdx = pThread->ThreadIdx + 1;
		ReleaseLock(true);

		pszInFile = pInPE1File->szFileName;
		pThread->NumInputFilesProcessed += 1;

		gDiagnostics.DiagOut(eDLInfo, gszProcName, "(Instance %d) Thread %d: Processing single ended reads from input read file '%s'", pThread->ProcessingID, pThread->ThreadIdx, pszInFile);
		if ((Rslt = ProcessReads(pThread,pInPE1File, NULL)) < eBSFSuccess)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "(Instance %d) Thread %d: Load failed for input sequences file '%s'", pThread->ProcessingID, pThread->ThreadIdx,pszInFile);
			AcquireLock(true);
			m_bTerminate = true;
			ReleaseLock(true);
			return((teBSFrsltCodes)Rslt);
			}
		NumReads = (int)pInPE1File->SeqCharacteristics.NumReads;

		if (Rslt == 0)
			gDiagnostics.DiagOut(eDLWarn, gszProcName, "(Instance %d) Thread %d: No sequences processed from input single ended reads file '%s'", pThread->ProcessingID, pThread->ThreadIdx, pszInFile);
		else
			{
			gDiagnostics.DiagOut(eDLInfo, gszProcName, "(Instance %d) Thread %d: Processed %u sequences from single ended reads file '%s'", pThread->ProcessingID, pThread->ThreadIdx, NumReads, pszInFile);
			gDiagnostics.DiagOut(eDLInfo, gszProcName, "(Instance %d) Thread %d: Unable to accept - %d were underlength, %d contained 'N's, %d below minimum Phred - reads from input reads file", pThread->ProcessingID,
				pInPE1File->SeqCharacteristics.NotProcUL,pInPE1File->SeqCharacteristics.NotProcNs,pInPE1File->SeqCharacteristics.NotProcQS);

			pThread->SeqCharacteristics.NotProcUL += pInPE1File->SeqCharacteristics.NotProcUL;
			pThread->SeqCharacteristics.NotProcNs += pInPE1File->SeqCharacteristics.NotProcNs;
			pThread->SeqCharacteristics.NotProcQS += pInPE1File->SeqCharacteristics.NotProcQS;
			if(pThread->SeqCharacteristics.MaxReadLen <  pInPE1File->SeqCharacteristics.MaxReadLen)
				pThread->SeqCharacteristics.MaxReadLen =  pInPE1File->SeqCharacteristics.MaxReadLen;
			if(pThread->SeqCharacteristics.MinReadLen == 0 || pThread->SeqCharacteristics.MinReadLen >  pInPE1File->SeqCharacteristics.MinReadLen)
				pThread->SeqCharacteristics.MinReadLen =  pInPE1File->SeqCharacteristics.MinReadLen;
			}

		pThread->TotNumSEReads += NumReads;
		}

	if (pThread->TotNumSEReads < 1)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "(Instance %d) Thread %d: No single end reads processed from input reads files - no further processing", pThread->ProcessingID, pThread->ThreadIdx);
		AcquireLock(true);
		m_bTerminate = true;
		ReleaseLock(true);
		return(eBSFerrNoEntries);
		}
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "(Instance %d) Thread %d: Processed total of %u reads from input single ended sequences files", pThread->ProcessingID, pThread->ThreadIdx, pThread->TotNumSEReads);
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "(Instance %d) Thread %d: Unable to accept - %lld were underlength, %lld contained 'N's, %lld below minimum Phred - total reads from input reads files", pThread->ProcessingID,
					pThread->ThreadIdx, pThread->SeqCharacteristics.NotProcUL,pThread->SeqCharacteristics.NotProcNs,pThread->SeqCharacteristics.NotProcQS);
	}
else // else must be paired end processing so process PE1 and PE2 as paired files
	{
	pInPE1File = &m_InReadsFiles[0];
	pInPE2File = &m_InReadsFiles[m_NumPE1InFiles];
	for (Idx = 0; Idx < m_NumPE1InFiles; Idx++, pInPE1File+=1, pInPE2File+=1)
		{
		AcquireLock(true);
		if(m_bTerminate)		// check if requested to early exit
			{
			ReleaseLock(true);
			gDiagnostics.DiagOut(eDLInfo, gszProcName, "(Instance %d) Thread %d: Requested to early terminate", pThread->ProcessingID, pThread->ThreadIdx);
			return(eBSFerrInternal);
			}
		if (pInPE1File->ProcByThreadIdx != 0)
			{
			ReleaseLock(true);
			continue;
			}

		pInPE1File->ProcByThreadIdx = pThread->ThreadIdx;
		pInPE2File->ProcByThreadIdx = pThread->ThreadIdx;
		ReleaseLock(true);

		pThread->NumInputFilesProcessed += 2;
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "(Instance %d) Thread %d: Processing paired end reads from input reads files '%s' and '%s'", pThread->ProcessingID, pThread->ThreadIdx, pInPE1File->szFileName, pInPE2File->szFileName);

		if ((Rslt = ProcessReads(pThread,pInPE1File, pInPE2File)) < eBSFSuccess)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "(Instance %d) Thread %d: Load failed for paired end sequences files '%s' and '%s'", pThread->ProcessingID, pThread->ThreadIdx, pInPE1File->szFileName, pInPE2File->szFileName);
			AcquireLock(true);
			m_bTerminate = true;
			ReleaseLock(true);
			return((teBSFrsltCodes)Rslt);
			}
		NumReads = (UINT32)Rslt;
		pThread->TotNumPEReads += NumReads;
		if (NumReads == 0)
			gDiagnostics.DiagOut(eDLWarn, gszProcName, "(Instance %d) Thread %d: No reads processed from input paired end files '%s' and '%s'", pThread->ProcessingID, pThread->ThreadIdx, pInPE1File->szFileName, pInPE2File->szFileName);
		else
			{
			gDiagnostics.DiagOut(eDLInfo, gszProcName, "(Instance %d) Thread %d: Processed %u paired ends from input paired end files '%s' and '%s'", pThread->ProcessingID, pThread->ThreadIdx, NumReads / 2, pInPE1File->szFileName, pInPE2File->szFileName);
			gDiagnostics.DiagOut(eDLInfo, gszProcName, "(Instance %d) Thread %d: Unable to accept - %d were underlength, %d contained 'N's, %d below minimum Phred - reads from input reads files", pThread->ProcessingID,
					pThread->ThreadIdx, pInPE1File->SeqCharacteristics.NotProcUL + pInPE2File->SeqCharacteristics.NotProcUL,
										pInPE1File->SeqCharacteristics.NotProcNs + pInPE2File->SeqCharacteristics.NotProcNs,
										pInPE1File->SeqCharacteristics.NotProcQS + pInPE2File->SeqCharacteristics.NotProcQS);

			pThread->SeqCharacteristics.NotProcUL += pInPE1File->SeqCharacteristics.NotProcUL + pInPE2File->SeqCharacteristics.NotProcUL;
			pThread->SeqCharacteristics.NotProcNs += pInPE1File->SeqCharacteristics.NotProcNs + pInPE2File->SeqCharacteristics.NotProcNs;
			pThread->SeqCharacteristics.NotProcQS += pInPE1File->SeqCharacteristics.NotProcQS + pInPE2File->SeqCharacteristics.NotProcQS;
			if(pThread->SeqCharacteristics.MaxReadLen <  pInPE1File->SeqCharacteristics.MaxReadLen)
				pThread->SeqCharacteristics.MaxReadLen =  pInPE1File->SeqCharacteristics.MaxReadLen;
			if(pThread->SeqCharacteristics.MinReadLen == 0 || pThread->SeqCharacteristics.MinReadLen >  pInPE1File->SeqCharacteristics.MinReadLen)
				pThread->SeqCharacteristics.MinReadLen =  pInPE1File->SeqCharacteristics.MinReadLen;
			}
		}

	if (pThread->TotNumPEReads < 2)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "(Instance %d) Thread %d: No paired end reads processed from input reads files - no further processing", pThread->ProcessingID, pThread->ThreadIdx, pThread->ThreadIdx);
		AcquireLock(true);
		m_bTerminate = true;
		ReleaseLock(true);
		return(eBSFerrNoEntries);
		}
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "(Instance %d) Thread %d: Processed total of %u paired ends (%u sequences) from all input reads files", pThread->ProcessingID, pThread->ThreadIdx, pThread->TotNumPEReads / 2, pThread->TotNumPEReads);
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "(Instance %d) Thread %d: Unable to accept - %lld were underlength, %lld contained 'N's, %lld below minimum Phred - total reads from input reads files", pThread->ProcessingID,
					pThread->ThreadIdx, pThread->SeqCharacteristics.NotProcUL,pThread->SeqCharacteristics.NotProcNs,pThread->SeqCharacteristics.NotProcQS);
	}
return(eBSFSuccess);
}

// ProcessReads
// Load and process reads from fasta, fastq formated raw reads file
teBSFrsltCodes
CReadStats::ProcessReads(tsThreadNGSQCPars *pThread, // thread specific processing state and context
				tsInReadsFile *pPE1File,		// file containing PE1 or SE reads
				tsInReadsFile *pPE2File)		// file containing PE2 reads if PE processing
{
static int FileNamesOfs = 0;
teBSFrsltCodes Rslt;
UINT32 NumSampledReads;

bool bIsSAMfile;
bool bIsPE1Fastq;
int NumPE1DescrReads;
UINT64 PE1TotReadsLen;
int PE1MinReadLen;
int PE1MaxReadLen;
int NumPE1AcceptedReads;
int NumPE1ParsedReads;
int PE1ReadLen;
int PE1DescrLen;
UINT8 szPE1DescrBuff[1024];

bool bIsPE2Fastq;
int NumPE2DescrReads;

int NumPE2AcceptedReads;
int NumPE2ParsedReads;
UINT64 PE2TotReadsLen;
int PE2MinReadLen;
int PE2MaxReadLen;
int PE2ReadLen;
int PE2DescrLen;
UINT8 szPE2DescrBuff[1024];

CFasta FastaPE1;
CFasta FastaPE2;

CSAMfile SAMfile;

bIsSAMfile = false;
bIsPE1Fastq = false;

// if not paired ends then check if a SAM/BAM file 
if(m_bPEProc == false)
	bIsSAMfile = CSAMfile::IsSAM(pPE1File->szFileName);

if(bIsSAMfile)
	{
	if((Rslt = (teBSFrsltCodes)SAMfile.Open(pPE1File->szFileName)) != eBSFSuccess)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "(Instance %d) Thread %d: Unable to open '%s' [%s] %s", pThread->ProcessingID, pThread->ThreadIdx, pPE1File->szFileName, FastaPE1.ErrText((teBSFrsltCodes)Rslt), FastaPE1.GetErrMsg());
		return(Rslt);
		}
	}
else
	{
	if ((Rslt = (teBSFrsltCodes)FastaPE1.Open(pPE1File->szFileName, true, cDfltStageBuffSize)) != eBSFSuccess)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "(Instance %d) Thread %d: Unable to open '%s' [%s] %s", pThread->ProcessingID, pThread->ThreadIdx, pPE1File->szFileName, FastaPE1.ErrText((teBSFrsltCodes)Rslt), FastaPE1.GetErrMsg());
		return(Rslt);
		}

	if(FastaPE1.IsSOLiD())
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "(Instance %d) Thread %d: Unable to load '%s', SOLiD colorspace not supported", pThread->ProcessingID, pThread->ThreadIdx, pPE1File->szFileName);
		return(eBSFerrFileType);
		}
	bIsPE1Fastq = FastaPE1.IsFastq();

	if (m_bPEProc)
		{
		if ((Rslt = (teBSFrsltCodes)FastaPE2.Open(pPE2File->szFileName, true, cDfltStageBuffSize)) != eBSFSuccess)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "(Instance %d) Thread %d: Unable to open '%s' [%s] %s", pThread->ProcessingID, pThread->ThreadIdx, pPE2File->szFileName, FastaPE2.ErrText((teBSFrsltCodes)Rslt), FastaPE2.GetErrMsg());
			FastaPE1.Close();
			return(Rslt);
			}

		if(FastaPE2.IsSOLiD())
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "(Instance %d) Thread %d: Unable to load '%s',SOLiD colorspace not supported", pThread->ProcessingID, pThread->ThreadIdx,  pPE2File->szFileName);
			FastaPE1.Close();
			FastaPE2.Close();
			return(eBSFerrFileType);
			}

		bIsPE2Fastq = FastaPE2.IsFastq();
		if(bIsPE1Fastq != bIsPE2Fastq)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "(Instance %d) Thread %d: Paired end file formats not of same type: '%s' is %s, '%s' is %s" , pThread->ProcessingID, pThread->ThreadIdx,
								 bIsPE1Fastq ? "Fastq" : "Fasta", pPE1File->szFileName, bIsPE2Fastq ? "Fastq" : "Fasta", pPE2File->szFileName);
			FastaPE1.Close();
			FastaPE2.Close();
			return(eBSFerrFileType);		
			}
		}
	}

NumSampledReads = 0;
NumPE1ParsedReads = 0;
NumPE1AcceptedReads = 0;
NumPE1DescrReads = 0;
PE1ReadLen = 0;
PE1DescrLen = 0;
NumPE2ParsedReads = 0;
NumPE2AcceptedReads = 0;
NumPE2DescrReads = 0;
PE2ReadLen = 0;
PE2DescrLen = 0;

pThread->PE1RawReadLen = 0;
pThread->PE2RawReadLen = 0;
pThread->PE1QScoresLen = 0;
pThread->PE1QScoresLen = 0;

PE1MinReadLen = 0;
PE1MaxReadLen = 0;
PE1TotReadsLen = 0;
PE2MinReadLen = 0;
PE2MaxReadLen = 0;
PE2TotReadsLen = 0;

memset(&pPE1File->SeqCharacteristics,0,sizeof(tsSeqCharacteristics));

if (m_bPEProc)
	memset(&pPE2File->SeqCharacteristics,0,sizeof(tsSeqCharacteristics));

time_t Started = time(0);
do {
	if(bIsSAMfile)
		Rslt = (teBSFrsltCodes)(PE1ReadLen = SAMfile.ReadSequence(pThread->PE1RawReadsBuff, cMaxRSSeqLen, true, false));
	else
		Rslt = (teBSFrsltCodes)(PE1ReadLen = FastaPE1.ReadSequence(pThread->PE1RawReadsBuff, cMaxRSSeqLen, true, false));

	if(Rslt <= eBSFSuccess)
		break;

	NumPE1ParsedReads += 1;
	if(!(NumPE1ParsedReads % 20000))
		{
		time_t Now = time(0);
		unsigned long ElapsedSecs = (unsigned long) (Now - Started);
		if(ElapsedSecs >= 120)
			{
			gDiagnostics.DiagOut(eDLInfo, gszProcName, "(Instance %d) Thread %d: Parsed %d reads %s", pThread->ProcessingID, pThread->ThreadIdx, NumPE1ParsedReads, m_bPEProc ? "paired sequences" : "sequences");
			Started = Now;
			}
		AcquireLock(true);
		if (m_bTerminate)		// check if requested to early exit
			{
			ReleaseLock(true);
			if(bIsSAMfile)
				SAMfile.Close();
			else
				{
				FastaPE1.Close();
				if (m_bPEProc)
					FastaPE2.Close();
				}
			gDiagnostics.DiagOut(eDLInfo, gszProcName, "(Instance %d) Thread %d: Requested to terminate", pThread->ProcessingID, pThread->ThreadIdx);
			return(eBSFerrInternal);
			}
		ReleaseLock(true);
		}

	if(PE1ReadLen == eBSFFastaDescr)		// just read a descriptor line which would be as expected for multifasta or fastq
		{
		NumPE1DescrReads += 1;
		if(bIsSAMfile)
			PE1DescrLen = SAMfile.ReadDescriptor((char *)szPE1DescrBuff,sizeof(szPE1DescrBuff)-1);
		else
			PE1DescrLen = FastaPE1.ReadDescriptor((char *)szPE1DescrBuff,sizeof(szPE1DescrBuff)-1);
		szPE1DescrBuff[sizeof(szPE1DescrBuff)-1] = '\0';

		if(bIsPE1Fastq)
			{
			pThread->PE1QScoresLen = FastaPE1.ReadQValues((char *)pThread->PE1QScoresBuff, cMaxRSSeqLen);
			pThread->PE1QScoresBuff[pThread->PE1QScoresLen] = '\0';
			}
		else
			{
			pThread->PE1QScoresLen = 0;
			pThread->PE1QScoresBuff[0] = '\0';
			}

		if(bIsSAMfile)
			PE1ReadLen = SAMfile.ReadSequence(pThread->PE1RawReadsBuff, cMaxRSSeqLen,true,false);
		else
			PE1ReadLen = FastaPE1.ReadSequence(pThread->PE1RawReadsBuff, cMaxRSSeqLen);
		if(PE1ReadLen < 1 || PE1ReadLen > cMaxRSSeqLen)
			{
			szPE1DescrBuff[80] = '\0';
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "(Instance %d) Thread %d: Problem parsing sequence after %d sequences parsed from file '%s'", pThread->ProcessingID, pThread->ThreadIdx, NumPE1DescrReads, pPE1File->szFileName);
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "(Instance %d) Thread %d: Last descriptor parsed: %s", pThread->ProcessingID, pThread->ThreadIdx, szPE1DescrBuff);
			if(bIsSAMfile)
				SAMfile.Close();
			else
				{
				FastaPE1.Close();
				if (m_bPEProc)
					FastaPE2.Close();
				}
			return(eBSFerrParse);
			}
		pThread->PE1RawReadLen = PE1ReadLen;
		PE1TotReadsLen += PE1ReadLen;
		if (PE1MinReadLen == 0 || PE1MinReadLen < PE1ReadLen)
			PE1MinReadLen = PE1ReadLen;
		if (PE1MaxReadLen < PE1ReadLen)
			PE1MaxReadLen = PE1ReadLen;

		if (pThread->PE1QScoresLen > 0 && pThread->PE1RawReadLen != pThread->PE1QScoresLen)
			{
			szPE1DescrBuff[80] = '\0';
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "(Instance %d) Thread %d: Problem parsing sequence after %d sequences parsed from file '%s'", pThread->ProcessingID, pThread->ThreadIdx, NumPE1DescrReads, pPE1File->szFileName);
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "(Instance %d) Thread %d: Expected sequence length (%d) to be same as quality length (%d). Last descriptor parsed: %s", pThread->ProcessingID, pThread->ThreadIdx,pThread->PE1RawReadLen,pThread->PE1QScoresLen, szPE1DescrBuff);
			if(bIsSAMfile)
				SAMfile.Close();
			else
				{
				FastaPE1.Close();
				if (m_bPEProc)
					FastaPE2.Close();
				}
			return(eBSFerrParse);
			}

		if (m_bPEProc)
			{
			Rslt = (teBSFrsltCodes)(PE2ReadLen = FastaPE2.ReadSequence(pThread->PE2RawReadsBuff, cMaxRSSeqLen, true, false));
			if(Rslt <= eBSFSuccess)
				{
				if(Rslt < eBSFSuccess)
					{
					gDiagnostics.DiagOut(eDLFatal, gszProcName, "(Instance %d) Thread %d: Errors processing file: %s ", pThread->ProcessingID, pThread->ThreadIdx, pPE2File->szFileName);
					while(FastaPE2.NumErrMsgs())
						gDiagnostics.DiagOut(eDLFatal,gszProcName,FastaPE2.GetErrMsg());
					}
				else
					gDiagnostics.DiagOut(eDLFatal, gszProcName, "(Instance %d) Thread %d: Insufficent reads in file: %s ", pThread->ProcessingID, pThread->ThreadIdx, pPE2File->szFileName);
				FastaPE1.Close();
				FastaPE2.Close();
				return(Rslt);
				}
			NumPE2ParsedReads += 1;
			if(PE2ReadLen == eBSFFastaDescr)		// just read a descriptor line which would be as expected for multifasta or fastq
				{
				NumPE2DescrReads += 1;
				PE2DescrLen = FastaPE2.ReadDescriptor((char *)szPE2DescrBuff,sizeof(szPE2DescrBuff)-1);
				szPE2DescrBuff[sizeof(szPE2DescrBuff)-1] = '\0'; 
				if(bIsPE2Fastq)
					{
					pThread->PE2QScoresLen = FastaPE2.ReadQValues((char *)pThread->PE2QScoresBuff, cMaxRSSeqLen);
					pThread->PE2QScoresBuff[pThread->PE2QScoresLen] = '\0';
					}
				else
					{
					pThread->PE2QScoresLen = 0;
					pThread->PE2QScoresBuff[0] = '\0';
					}
				PE2ReadLen = FastaPE2.ReadSequence(pThread->PE2RawReadsBuff, cMaxRSSeqLen);
				if(PE2ReadLen < 1 || PE2ReadLen > cMaxRSSeqLen)
					{
					szPE2DescrBuff[80] = '\0';
					gDiagnostics.DiagOut(eDLFatal, gszProcName, "(Instance %d) Thread %d: Problem parsing sequence after %d sequences parsed from file '%s'", pThread->ProcessingID, pThread->ThreadIdx, NumPE2DescrReads, pPE2File->szFileName);
					gDiagnostics.DiagOut(eDLFatal, gszProcName, "(Instance %d) Thread %d: Last descriptor parsed: %s", pThread->ProcessingID, pThread->ThreadIdx, szPE2DescrBuff);
					FastaPE1.Close();
					FastaPE2.Close();
					return(eBSFerrParse);
					}
				pThread->PE2RawReadLen = PE2ReadLen;
				PE2TotReadsLen += PE2ReadLen;
				if (PE2MinReadLen == 0 || PE2MinReadLen < PE2ReadLen)
					PE2MinReadLen = PE2ReadLen;
				if (PE2MaxReadLen < PE2ReadLen)
					PE2MaxReadLen = PE2ReadLen;
				if (pThread->PE2QScoresLen > 0 && pThread->PE2RawReadLen != pThread->PE2QScoresLen)
					{
					szPE2DescrBuff[80] = '\0';
					gDiagnostics.DiagOut(eDLFatal, gszProcName, "(Instance %d) Thread %d: Problem parsing sequence after %d sequences parsed from file '%s'", pThread->ProcessingID, pThread->ThreadIdx, NumPE2DescrReads, pPE2File->szFileName);
					gDiagnostics.DiagOut(eDLFatal, gszProcName, "(Instance %d) Thread %d: Expected sequence length to be same as quality length. Last descriptor parsed: %s", pThread->ProcessingID, pThread->ThreadIdx, szPE2DescrBuff);
					FastaPE1.Close();
					if (m_bPEProc)
						FastaPE2.Close();
					return(eBSFerrParse);
					}
				}
			else
				{
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "(Instance %d) Thread %d: Raw sequence file '%s' processing error: %s ", pThread->ProcessingID, pThread->ThreadIdx, pPE2File->szFileName,
							Rslt > 0 ? "over length read" : "not a multifasta short reads or fastq file");
				FastaPE1.Close();
				FastaPE2.Close();
				return(eBSFerrParse);
				}
			}

		if (!m_bPEProc)
			{
			if ((Rslt = AnalyseReads(pThread,PE1ReadLen, pThread->PE1RawReadsBuff, pThread->PE1QScoresLen, pThread->PE1QScoresBuff, pPE1File)) != eBSFSuccess)
				break;
			NumPE1AcceptedReads += 1;
			}
		else
			{
			if ((Rslt = AnalyseReads(pThread,PE1ReadLen, pThread->PE1RawReadsBuff, pThread->PE1QScoresLen, pThread->PE1QScoresBuff, pPE1File, PE2ReadLen, pThread->PE2RawReadsBuff, pThread->PE2QScoresLen, pThread->PE2QScoresBuff, pPE2File)) != eBSFSuccess)
				break;
			NumPE1AcceptedReads += 1;
			NumPE2AcceptedReads += 1;
			}

		continue;
		}
	else
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "(Instance %d) Thread %d: Raw sequence file '%s' processing error: %s ", pThread->ProcessingID, pThread->ThreadIdx, pPE1File->szFileName,
				Rslt > 0 ? "over length read" : "not a multifasta short reads or fastq file");
		if(bIsSAMfile)
			SAMfile.Close();
		else
			{
			FastaPE1.Close();
			if (m_bPEProc)
				FastaPE2.Close();
			}
		return(eBSFerrParse);
		}
	}
while(1);

if(Rslt < eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "(Instance %d) Thread %d: Errors processing file: %s ", pThread->ProcessingID, pThread->ThreadIdx, pPE1File->szFileName);
	while(FastaPE1.NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,FastaPE1.GetErrMsg());
	if(bIsSAMfile)
		SAMfile.Close();
	else
		{
		FastaPE1.Close();
		if (m_bPEProc)
			FastaPE2.Close();
		}
	return(Rslt);
	}

if(bIsSAMfile)
	SAMfile.Close();
else
	{
	FastaPE1.Close();
	if (m_bPEProc)
		FastaPE2.Close();
	}
gDiagnostics.DiagOut(eDLInfo, gszProcName, "(Instance %d) Thread %d: Parsed %d %s", pThread->ProcessingID, pThread->ThreadIdx, NumPE1ParsedReads, m_bPEProc ? "paired sequences" : "sequences");

pPE1File->SeqCharacteristics.MaxReadLen = PE1MaxReadLen;
pPE1File->SeqCharacteristics.MinReadLen = PE1MinReadLen;
pPE1File->SeqCharacteristics.MeanReadLen = (int)(PE1TotReadsLen / NumPE1AcceptedReads);

if(m_bPEProc)
	{
	pPE2File->SeqCharacteristics.MaxReadLen = PE2MaxReadLen;
	pPE2File->SeqCharacteristics.MinReadLen = PE2MinReadLen;
	pPE2File->SeqCharacteristics.MeanReadLen = (int)(PE2TotReadsLen / NumPE2AcceptedReads);
	}

return((teBSFrsltCodes)(NumPE1AcceptedReads + NumPE2AcceptedReads));
}

// Sampled read instances are only accepted if the read contains no indeterminates and within the max number of unique instances allowed
int				// returns 0 if not accepted as read instance, 1 if this is the first instance of an accepted sampled read, 2..N if multiple instances exist
CReadStats::AddReadInst(tsThreadNGSQCPars *pThread, // thread specific processing state and context
				int PE1ReadLen,				// number of bases in PE1 read
				UINT8 *pPE1RawRead,			// PE1 read sequence
				int PE2ReadLen,				// number of bases in PE2 read
				UINT8 *pPE2RawRead)			// PE2 read sequence
{
UINT32 PackedSeqs[((cMaxRSSeqLen+15)/16)*2];	// to hold both PE1 and PE2 (if PE processing) packed and concatenated together
UINT32 PackedBases;								// bases are packed 16 per word
int PartialPacked;								// number of bases currently packed into PackedBases
UINT32 *pPackedSeqs;
UINT32 SeqHash;									// sequence hash
UINT32 RevCplSeqHash;
UINT32 HashSeqsOfs;
UINT32 TermRevCplNxtSeq;
UINT32 TermNxtSeq;
tsSampledSeq *pSampledSeq;
tsSampledSeq *pNewSampledSeq;
int Idx;
int NumPackedWrds;
etSeqBase SeqBase;

UINT32 RevCplPackedSeqs[((cMaxRSSeqLen + 15) / 16) * 2];	// to hold RevCpl of both PE1 and PE2 (if PE processing) packed and concatenated together
UINT8 *pPE1RevCplRawRead;
UINT8 *pPE2RevCplRawRead;
UINT32 *pRevCplPackedSeqs;
tsSampledSeq *pRevCplSampledSeq;
UINT32 RevCplHashSeqsOfs;
UINT32 RevCplPackedBases;
int RevCplNumPackedWrds;
int RevCplPartialPacked;
int NumInstances;

if (!m_bStrand)
	{
	pPE1RevCplRawRead = &pPE1RawRead[PE1ReadLen-1];
	if(m_bPEProc)
		pPE2RevCplRawRead = &pPE2RawRead[PE2ReadLen - 1];
	else
		pPE2RevCplRawRead = NULL;
	}


pPackedSeqs = PackedSeqs;
PackedBases = 0;
NumPackedWrds = 0;
PartialPacked = 0;
SeqHash = 0;
for(Idx=0;Idx<PE1ReadLen;Idx++,pPE1RawRead++)
	{
	if((SeqBase = (etSeqBase)(*pPE1RawRead & 0x07)) > eBaseT)		// can only accept cannonical bases
		return(0);
	PackedBases <<= 2;
	PackedBases |= SeqBase;
	PartialPacked += 1;
	if(Idx > 0 && !(Idx % 16))
		{
		*pPackedSeqs++ = PackedBases;
		SeqHash <<= 1;
		SeqHash ^= PackedBases;
		NumPackedWrds += 1; 
		PackedBases = 0;
		PartialPacked = 0;
		}
	}

if(m_bPEProc)
	{
	for (Idx = 0; Idx < PE2ReadLen; Idx++, pPE2RawRead++)
		{
		if ((SeqBase = (etSeqBase)(*pPE2RawRead & 0x07)) > eBaseT)		// can only accept cannonical bases
			return(0);
		PackedBases <<= 2;
		PackedBases |= SeqBase;
		PartialPacked += 1;
		if(Idx > 0 && !(Idx % 16))
			{
			*pPackedSeqs++ = PackedBases;
			SeqHash <<= 1;
			SeqHash ^= PackedBases;
			NumPackedWrds += 1;
			PackedBases = 0;
			PartialPacked = 0;
			}
		}
	}
if (PartialPacked)
	{
	*pPackedSeqs = PackedBases;
	SeqHash <<= 1;
	SeqHash ^= PackedBases;
	NumPackedWrds += 1;
	}

SeqHash &= cHashMask;

// if not strand dependent then worth the cost of reverse complementing the reads even if not needed because identical read sequences already sampled
if (!m_bStrand)
	{
	pRevCplPackedSeqs = RevCplPackedSeqs;
	RevCplPackedBases = 0;
	RevCplNumPackedWrds = 0;
	RevCplPartialPacked = 0;
	RevCplSeqHash = 0;

	if (m_bPEProc)
		{
		for (Idx = 0; Idx < PE2ReadLen; Idx++, pPE2RevCplRawRead--)
			{
			if ((SeqBase = (etSeqBase)(*pPE2RevCplRawRead & 0x07)) > eBaseT)		// can only accept cannonical bases
				return(0);
			switch (SeqBase) {
				case eBaseA:
					SeqBase = eBaseT;
					break;
				case eBaseC:
					SeqBase = eBaseG;
					break;
				case eBaseG:
					SeqBase = eBaseC;
					break;
				case eBaseT:
					SeqBase = eBaseA;
					break;
				}
			RevCplPackedBases <<= 2;
			RevCplPackedBases |= SeqBase;
			RevCplPartialPacked += 1;
			if (Idx > 0 && !(Idx % 16))
				{
				*pRevCplPackedSeqs++ = RevCplPackedBases;
				RevCplSeqHash <<= 1;
				RevCplSeqHash ^= RevCplPackedBases;
				RevCplNumPackedWrds += 1;
				RevCplPackedBases = 0;
				RevCplPartialPacked = 0;
				}
			}
		}
	for (Idx = 0; Idx < PE1ReadLen; Idx++, pPE1RevCplRawRead--)
		{
		if ((SeqBase = (etSeqBase)(*pPE1RevCplRawRead & 0x07)) > eBaseT)		// can only accept cannonical bases
			return(0);
		switch (SeqBase)
			{
			case eBaseA:
				SeqBase = eBaseT;
				break;
			case eBaseC:
				SeqBase = eBaseG;
				break;
			case eBaseG:
				SeqBase = eBaseC;
				break;
			case eBaseT:
				SeqBase = eBaseA;
				break;
			}
		RevCplPackedBases <<= 2;
		RevCplPackedBases |= SeqBase;
		RevCplPartialPacked += 1;
		if (Idx > 0 && !(Idx % 16))
			{
			*pRevCplPackedSeqs++ = RevCplPackedBases;
			RevCplSeqHash <<= 1;
			RevCplSeqHash ^= RevCplPackedBases;
			RevCplNumPackedWrds += 1;
			RevCplPackedBases = 0;
			RevCplPartialPacked = 0;
			}
		}
	if (RevCplPartialPacked)
		{
		*pRevCplPackedSeqs = RevCplPackedBases;
		RevCplSeqHash <<= 1;
		RevCplSeqHash ^= RevCplPackedBases;
		RevCplNumPackedWrds += 1;
		}

	RevCplSeqHash &= cHashMask;
	}

// check if this packed sequence has been previously accepted as a sample
pSampledSeq = NULL;
TermNxtSeq = 0;
TermRevCplNxtSeq = 0;
AcquireSerialise();
if (m_ActMaxDupSeeds && (HashSeqsOfs = m_pSeqHashes[SeqHash]) != 0)		// seen at least one sequence previously with same hash?
	{
	TermNxtSeq = HashSeqsOfs;
	pSampledSeq = (tsSampledSeq *)&m_pSampledSeqs[HashSeqsOfs-1];
	do {
		if(pSampledSeq->PE1ReadLen == PE1ReadLen && (!m_bPEProc || (pSampledSeq->PE2ReadLen == PE2ReadLen)))
			{
			if(pSampledSeq->PackedSeqs[0] == PackedSeqs[0])
				{
				if((pSampledSeq->NumPackedWrds == 1) || !memcmp(pSampledSeq->PackedSeqs,PackedSeqs,pSampledSeq->NumPackedWrds * sizeof(UINT32)))
					{
					pSampledSeq->NumInstances += 1;
					NumInstances = (int)pSampledSeq->NumInstances;
					ReleaseSerialise();
					return(NumInstances);
					}
				}
			}
		if(pSampledSeq->NxtSeq != 0)
			TermNxtSeq = pSampledSeq->NxtSeq;
		}
	while (pSampledSeq->NxtSeq != 0 && (pSampledSeq = (tsSampledSeq *)&m_pSampledSeqs[pSampledSeq->NxtSeq - 1]));
	}

// if no success with finding existing sense natch and not strand dependent then try with sequence reverse complemented
if (!m_bStrand && m_ActMaxDupSeeds)
	{
	// check if this packed sequence has been previously accepted as a sample
	pRevCplSampledSeq = NULL;
	if ((RevCplHashSeqsOfs = m_pSeqHashes[RevCplSeqHash]) != 0)		// seen at least one sequence previously with same hash?
		{
		TermRevCplNxtSeq = RevCplHashSeqsOfs;
		pRevCplSampledSeq = (tsSampledSeq *)&m_pSampledSeqs[RevCplHashSeqsOfs - 1];
		do
			{
			if (pRevCplSampledSeq->PE1ReadLen == PE1ReadLen && (!m_bPEProc || (pRevCplSampledSeq->PE2ReadLen == PE2ReadLen)))
				{
				if (pRevCplSampledSeq->PackedSeqs[0] == RevCplPackedSeqs[0])
					{
					if((pRevCplSampledSeq->NumPackedWrds == 1) || !memcmp(pRevCplSampledSeq->PackedSeqs, RevCplPackedSeqs, pRevCplSampledSeq->NumPackedWrds * sizeof(UINT32)))
						{
						pRevCplSampledSeq->NumInstances += 1;
						pRevCplSampledSeq->NumRevCplInstances += 1;
						NumInstances = (int)pRevCplSampledSeq->NumInstances;
						ReleaseSerialise();
						return(NumInstances);
						}
					}
				}
			if(pRevCplSampledSeq->NxtSeq != 0)
				TermRevCplNxtSeq = pRevCplSampledSeq->NxtSeq;
			}
		while (pRevCplSampledSeq->NxtSeq != 0 && (pRevCplSampledSeq = (tsSampledSeq *)&m_pSampledSeqs[pRevCplSampledSeq->NxtSeq - 1]));
		}
	}

if (m_ActMaxDupSeeds >= m_ReqMaxDupSeeds)
	{
	ReleaseSerialise();
	return(0);
	}

// realloc memory if required
if((m_AllocdSampledSeqWrds - m_UsedSampledSeqWrds) < ((cMaxRSSeqLen+3)/4)*10)
	{
	UINT32 *pAllocd;
	size_t memreq = m_AllocdSampledSeqMem + (cReallocPackedWrds * sizeof(UINT32));
#ifdef _WIN32
	pAllocd = (UINT32 *)realloc(m_pSampledSeqs, memreq);
#else
	pAllocd = (UINT32 *)mremap(m_pSampledSeqs,m_AllocdSampledSeqMem,memreq,MREMAP_MAYMOVE);
	if(pAllocd == MAP_FAILED)
		pAllocd = NULL;
#endif
	if(pAllocd == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Thread %d: Memory re-allocation to %d bytes - %s",pThread->ThreadIdx, memreq,strerror(errno));
		ReleaseSerialise();
		return(eBSFerrMem);
		}
	m_pSampledSeqs = pAllocd;
	m_AllocdSampledSeqMem = memreq;
	m_AllocdSampledSeqWrds += cReallocPackedWrds;
	}

pNewSampledSeq =(tsSampledSeq *)&m_pSampledSeqs[m_UsedSampledSeqWrds];
pNewSampledSeq->NumInstances = 1;
pNewSampledSeq->Flags = 0;
pNewSampledSeq->NumRevCplInstances = 0;
pNewSampledSeq->NxtSeq = 0;
pNewSampledSeq->PE1ReadLen = PE1ReadLen;
pNewSampledSeq->PE2ReadLen = m_bPEProc ? PE2ReadLen : 0;
pNewSampledSeq->NumPackedWrds = NumPackedWrds;
memcpy(pNewSampledSeq->PackedSeqs,PackedSeqs,NumPackedWrds * sizeof(UINT32));
if (TermNxtSeq == 0)
	{
	 m_pSeqHashes[SeqHash] = m_UsedSampledSeqWrds + 1;
	 m_NumSeqHashes += 1;
	}
else
	((tsSampledSeq *)&m_pSampledSeqs[TermNxtSeq - 1])->NxtSeq = m_UsedSampledSeqWrds + 1;
m_UsedSampledSeqWrds += NumPackedWrds + (((sizeof(tsSampledSeq)-sizeof(UINT32)) + 3)/4);
m_ActMaxDupSeeds += 1;
ReleaseSerialise();
return(1);
}

double									// returned prob of read being error free
GenProbErrFreeRead(int QSSchema,		// quality scoring schema - guestimated scoring schema - 0: no scoring, 1: Solexa, 2: Illumina 1.3+, 3: Illumina 1.5+, 4: Illumina 1.8+ or Sanger
				   int ReadLen,			// read length
				   UINT8 *pQScores)		// Phred quality scores
{
int SeqOfs;
int Score;
double ProbErr;
double ProbNoReadErr;
if(QSSchema == 0 || ReadLen < 1 || pQScores == NULL || pQScores[0] == '\0')		// if can't score then have to assume the best - read is sequencing base call error free
	return(1.0);

ProbNoReadErr = 1.0;
for (SeqOfs = 0; SeqOfs < ReadLen; SeqOfs++, pQScores++)
	{
	switch(QSSchema) { // quality scoring schema - guestimated scoring schema - 0: no scoring, 1: Solexa, 2: Illumina 1.3+, 3: Illumina 1.5+, 4: Illumina 1.8+ or Sanger 
		case 1:	// Solexa  ';' (-5) to 'h' (40)
			if(*pQScores <= '@')
				Score = 0;
			else
				Score = (int)(*pQScores - (UINT8)'@');
			break;

		case 2: // Illumina 1.3+ '@' (0) to 'h' (40)
			Score = (int)(*pQScores - (UINT8)'@');
			break;

		case 3: // Illumina 1.5+ 'B' (2) to 'h' (40)
			Score = (int)(*pQScores - (UINT8)'@');
			break;

		case 4: // Illumina 1.8+ or Sanger. Sanger is '!' (0) to 'J' (41) and Illumina is '#' (2) to 'J' (41)
			Score = (int)(*pQScores - (UINT8)'!');
			break;
			}
	if(Score < 0)			// force scores to always be in the range 0..41 --- shouldn't be outside 0..41 but some qscores have been observed to be outside the expected range!
		Score = 0;
	else
		if(Score > 41)
			Score = 41;

	ProbErr = 1.0 / pow(10.0,(double)Score/10.0);
	ProbNoReadErr *= (1.0 - ProbErr); 
	}
if(ProbNoReadErr < 0.0)
	ProbNoReadErr = 0.0;
else
	if(ProbNoReadErr > 1.0)
		ProbNoReadErr = 1.0;
return(ProbNoReadErr);
}

// accumulate quality scores at each base offset along the read length
int				// minimum score at any offset within the sequence
CReadStats::AccumQScores(tsThreadNGSQCPars *pThread, // thread specific processing state and context
				int QSSchema,					// quality scoring schema - guestimated scoring schema - 0: no scoring, 1: Solexa, 2: Illumina 1.3+, 3: Illumina 1.5+, 4: Illumina 1.8+ or Sanger 
				int ReadLen,					// number of bases in read
				UINT8 *pSeq,					// read sequence
				UINT8 *pQScores)				// read quality scores (NULL if no associated quality scores)
{
int MinScore;
int SeqOfs;
int Score;
int SumBaseScores;
int MeanReadScore;
double ProbErr;
double ProbNoReadErr;
UINT32 Scores[cMaxRSSeqLen];
UINT32 Ns[cMaxRSSeqLen];
UINT32 *pScores;
UINT32 *pNs;
UINT32 *pScore;
UINT32 *pBaseNs;
bool bTrunc;

if(pThread == NULL || QSSchema < 0 || QSSchema > 4 || ReadLen < cMinRSSeqLen || pSeq == NULL)
	return(eBSFerrParams);
if(ReadLen > cMaxRSSeqLen)					// silently truncate if overlength
	{
	ReadLen = cMaxRSSeqLen;
	bTrunc = true;
	}
else
	bTrunc = false;

// scores are parsed into a local copy then when sequence completely parsed then the local copy is updated into the global scores within a single serialisation lock
MinScore = 0;
ProbNoReadErr = 1.0;
SumBaseScores = 0;
pScore = Scores;
pNs = Ns;
for (SeqOfs = 0; SeqOfs < ReadLen; SeqOfs++, pSeq++, pNs++, pScore++,pQScores++)
	{
	if(*pSeq > eBaseT)
		*pNs = 1;
	else
		*pNs = 0;
	
	switch(QSSchema) { // quality scoring schema - guestimated scoring schema - 0: no scoring, 1: Solexa, 2: Illumina 1.3+, 3: Illumina 1.5+, 4: Illumina 1.8+ or Sanger 
		case 1:	// Solexa  ';' (-5) to 'h' (40)
			if(*pQScores <= '@')
				Score = 0;
			else
				Score = (int)(*pQScores - (UINT8)'@');
			break;

		case 2: // Illumina 1.3+ '@' (0) to 'h' (40)
			Score = (int)(*pQScores - (UINT8)'@');
			break;

		case 3: // Illumina 1.5+ 'B' (2) to 'h' (40)
			Score = (int)(*pQScores - (UINT8)'@');
			break;

		case 4: // Illumina 1.8+ or Sanger. Sanger is '!' (0) to 'J' (41) and Illumina is '#' (2) to 'J' (41)
			Score = (int)(*pQScores - (UINT8)'!');
			break;
			
		case 0:	// no scoring
			Score = cDfltPhredScore;
			break;
			}
	if(Score < 0)			// force scores to always be in the range 0..41 --- shouldn't be outside 0..41 but some qscores have been observed to be outside the expected range!
		Score = 0;
	else
		if(Score > 41)
			Score = 41;
	*pScore = (UINT32)Score;

	ProbErr = 1.0 / pow(10.0,(double)Score/10.0);
	ProbNoReadErr *= (1.0 - ProbErr); 
	SumBaseScores += Score;
	if(MinScore == 0 || MinScore > Score)
		MinScore = Score;
	}
MeanReadScore = SumBaseScores /  ReadLen;

pScores = m_pScores;
pScore = Scores;
pNs = Ns;
pBaseNs = m_pBaseNs;
AcquireSerialiseScores();
for (SeqOfs = 0; SeqOfs < ReadLen; SeqOfs++, pScore++, pBaseNs++, pNs++, pScores += 42)
	{
	pScores[*pScore] += 1;
	*pBaseNs += *pNs;
	}
if(ReadLen > m_MaxReadLen)
	m_MaxReadLen = ReadLen;
if(m_MinReadLen == 0 || m_MinReadLen > ReadLen)
	m_MinReadLen = ReadLen;
m_ReadLenDist[ReadLen] += 1;
if(bTrunc)
	m_ReadLenDist[cMaxRSSeqLen+1] += 1;

int ProbNoReadErrBinIdx;
ProbNoReadErrBinIdx = (int)(ProbNoReadErr * 100);
if(ProbNoReadErrBinIdx > 99)
	ProbNoReadErrBinIdx = 99;
m_ProbNoReadErrDist[ProbNoReadErrBinIdx] += 1;
ReleaseSerialiseScores();
return(MinScore);
}

// accumulate K-mer counts at each base offset along the read length
int				
CReadStats::AccumKMers(tsThreadNGSQCPars *pThread, // thread specific processing state and context
				int ReadLen,					// number of bases in read
				UINT8 *pSeq)					// read sequence
{
int SeqOfs;
UINT32 Base;
UINT32 PackedKMer;
UINT32 PackedKMerMsk;
UINT32 KMerMsk;
int KMerLen;
int KmerCntsOfs;
UINT32 *pKMerCntsSeqOfs;
UINT32 KMerOfs;
UINT32 KMerLenOfs;
int Pow;

if(pThread == NULL || ReadLen < cMinRSSeqLen || pSeq == NULL)
	return(eBSFerrParams);
if(ReadLen > cMaxRSSeqLen)					// silently truncate if overlength
	ReadLen = cMaxRSSeqLen;
PackedKMer = 0;
PackedKMerMsk = ~(0x03 << (m_MaxKMerLen * 2));
pKMerCntsSeqOfs = pThread->KMerCntOfs;
KMerOfs = 0;

for (SeqOfs = 0; SeqOfs < (ReadLen + m_MaxKMerLen); SeqOfs++, pSeq++)
	{
	if(SeqOfs < ReadLen)
		{
		if((Base = (UINT32)*pSeq) > eBaseT)			// if an indeterminate (shouldn't be any as only sequences accepted for duplicate processing should be k-mer processed)
			Base = (UINT32)(SeqOfs % 4);			// substitute a cannonical base for indeterminate - should be near random as expecting indeterminates to be uniformly distributed along read length
		PackedKMer <<= 2;
		PackedKMer &= PackedKMerMsk;
		PackedKMer |= Base;
		if(SeqOfs < (m_MaxKMerLen-1))
			continue;
		}
	else
		{
		PackedKMer <<= 2;							// effectively imputes 'a's but these cnts are ignored when accumulating counts 
		PackedKMer &= PackedKMerMsk;
		}
	KMerMsk = 0x03 << ((m_MaxKMerLen-1) * 2);
	KMerLenOfs = 0;
	Pow = 1;
	for(KMerLen = 0; KMerLen < m_MaxKMerLen; KMerLen++)
		{
		KmerCntsOfs = (PackedKMer & KMerMsk) >> ((m_MaxKMerLen - KMerLen - 1) * 2);
		*pKMerCntsSeqOfs++ = KmerCntsOfs + KMerLenOfs + KMerOfs;
		Pow <<= 2;
		KMerLenOfs += Pow;
		KMerMsk |= (KMerMsk >> 2);
		}
	KMerOfs += m_KMerCntsEls;	
	}

AcquireSerialiseKMers();
// ensure allocation for KMerCnts is sufficent for read lengths
if(ReadLen > m_AllocdMaxReadLen)
	{
	UINT32 *pAllocd;
	m_AllocdMaxReadLen = (ReadLen * 120) / 100;
	size_t memreq = m_KMerCntsEls * m_AllocdMaxReadLen * sizeof(UINT32);
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "AddReadInst: Memory re-allocation to %d bytes", memreq);
#ifdef _WIN32
	pAllocd = (UINT32 *)realloc(m_pKMerCnts, memreq);
#else
	pAllocd = (UINT32 *)mremap(m_pKMerCnts,m_AllocdKMerCntsMem,memreq,MREMAP_MAYMOVE);
	if(pAllocd == MAP_FAILED)
		pAllocd = NULL;
#endif
	if(pAllocd == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Thread %d: Memory re-allocation to %d bytes - %s",pThread->ThreadIdx, memreq,strerror(errno));
		ReleaseSerialise();
		return(eBSFerrMem);
		}
	m_pKMerCnts = pAllocd;
	m_AllocdKMerCntsMem = memreq;
	}

pKMerCntsSeqOfs = pThread->KMerCntOfs;
for(SeqOfs = 0; SeqOfs < ReadLen; SeqOfs++)
	{
	for(KMerLen = 0; KMerLen < m_MaxKMerLen; KMerLen++,pKMerCntsSeqOfs++)
		if(SeqOfs + KMerLen < ReadLen)					// sloughing cnts for k-mers which would have extended out past the end of reads
			m_pKMerCnts[*pKMerCntsSeqOfs] += 1;		
	}	
ReleaseSerialiseKMers();
return(0);
}

// Pearson sample correlation coefficient
// Defined as the covariance of the two variables divided by the product of their standard deviations
// A count of 1 is added to both control and experiment so as to provide for Laplaces smoothing and prevent divide by zero errors
double									// returned Pearson
CReadStats::Pearsons(int KMerLen,		// KMer length for which counts were accumulated 
					UINT32 *pCtrlKMerCnts,  // KMer counts for control
					UINT32 *pExprKMerCnts)	// KMer counts for experiment	
{
UINT32 *pCtrl;
UINT32 *pExpr;

int Idx;
int NumLociCnts;
double MeanC;
double MeanE;
double TmpC;
double TmpE;
double SumNum;
double SumDenC;
double SumDenE;
double Correl;
int NumEls;

MeanC = 0.0;
MeanE = 0.0;
NumLociCnts = 0;

// calc number of elements for specified KmerLen as pow(4,KMerLen)
// 1->4, 2->16, 3->64, 4->256
NumEls = 1;
for(Idx = 0; Idx < KMerLen; Idx++)
	NumEls <<= 2;

// calc the means, assuming same number of reads from which KMer counts derived then means should be identical..
pCtrl = pCtrlKMerCnts;
pExpr = pExprKMerCnts;
for(Idx = 0; Idx < NumEls; Idx++, pCtrl++,pExpr++)
	{
	MeanC += *pCtrl+1.0;						// adding a psuedocount of 1 is Laplace's smoothing
	MeanE += *pExpr+1.0;							// which also conveniently ensures can never have divide by zero errors!
	}

MeanC /= NumEls;								// can now determine the means
MeanE /= NumEls;

SumNum = 0.0;
SumDenC = 0.0;
SumDenE = 0.0;
pCtrl = pCtrlKMerCnts;
pExpr = pExprKMerCnts;
for(Idx = 0; Idx < NumEls; Idx++, pCtrl++,pExpr++)
	{
	TmpC = *pCtrl - MeanC;
	TmpE = *pExpr - MeanE;
	SumNum += (TmpC * TmpE);
	SumDenC += (TmpC * TmpC);
	SumDenE += (TmpE * TmpE);
	}

if(SumDenC < 0.000001)			// set a floor so as to prevent the chance of a subsequent divide by zero error
	SumDenC = 0.000001;
if(SumDenE < 0.000001)
	SumDenE = 0.000001;
Correl = SumNum/sqrt(SumDenC*SumDenE);
return(Correl);
}

teBSFrsltCodes
CReadStats::AnalyseReads(tsThreadNGSQCPars *pThread, // thread specific processing state and context
		int PE1ReadLen,				    // number of bases in PE1 read
		UINT8 *pPE1RawRead,				// PE1 read sequence
		int PE1QScoreLen,				// number of quality scores in PE1 read
		UINT8 *pPE1QScores,				// PE1 read quality scores (0 if no associated quality scores)
		tsInReadsFile *pPE1File,		// file containing PE1 or SE reads
		int PE2ReadLen,					// number of bases in PE2 read
		UINT8 *pPE2RawRead,				// PE2 read sequence
		int PE2QScoreLen,				// number of quality scores in PE2 read (0 if no associated quality scores)
		UINT8 *pPE2QScores,				// PE2 read quality scores
		tsInReadsFile *pPE2File)		// file containing PE2
{
int NumInsts;
int SeqOfs;
UINT8 *pBase;
UINT8 *pScore;
bool bPE1Contaminated;
bool bPE2Contaminated;

// remove any flags which may be present in the sequences
if(PE1ReadLen > 0 && pPE1RawRead != NULL)
	{
	pBase = pPE1RawRead;
	for(SeqOfs = 0; SeqOfs < PE1ReadLen; SeqOfs++, pBase++)
		*pBase &= 0x07;
	}
if(PE2ReadLen > 0 && pPE2RawRead != NULL)
	{
	pBase = pPE2RawRead;
	for(SeqOfs = 0; SeqOfs < PE2ReadLen; SeqOfs++, pBase++)
		*pBase &= 0x07;
	}

// check that read lengths, after any applied end trimming, would still be an acceptable length
if ((PE1ReadLen - (m_Trim5 + m_Trim3) < cMinRSSeqLen) ||
	(m_bPEProc && (PE2ReadLen - (m_Trim5 + m_Trim3) < cMinRSSeqLen)))
	{
	pPE1File->SeqCharacteristics.NumNoneCntd += 1;
	if ((PE1ReadLen - (m_Trim5 + m_Trim3) < cMinRSSeqLen))
		pPE1File->SeqCharacteristics.NotProcUL += 1;
	if (m_bPEProc)
		{
		pPE2File->SeqCharacteristics.NumNoneCntd += 1;
		if ((PE2ReadLen - (m_Trim5 + m_Trim3) < cMinRSSeqLen))
			pPE2File->SeqCharacteristics.NotProcUL += 1;
		}
	pPE1File->SeqCharacteristics.NumReads += 1;
	if (m_bPEProc)
		pPE2File->SeqCharacteristics.NumReads += 1;
	return(eBSFSuccess);
	}

// apply any end trimmming as may have been requested by user
if (m_Trim5)
	{
	pPE1RawRead += m_Trim5;
	pPE1QScores += m_Trim5;
	PE1ReadLen -= m_Trim5;
	if (m_bPEProc)
		{
		pPE2RawRead += m_Trim5;
		pPE2QScores += m_Trim5;
		PE2ReadLen -= m_Trim5;
		}
	}

if (m_Trim3)
	{
	PE1ReadLen -= m_Trim3;
	if(m_bPEProc)
		PE2ReadLen -= m_Trim3;
	}

// accumulate quality score counts
AccumQScores(pThread,pPE1File->QSSchema,PE1ReadLen,pPE1RawRead,pPE1QScores);
if(m_bPEProc)
	AccumQScores(pThread,pPE2File->QSSchema,PE2ReadLen,pPE2RawRead,pPE2QScores);

// use this read as a sample?
// will only sample reads which have no bases less than a Phred score of 10 and a mean thred score at least m_MinMeanPhredScore, and which also contain no 'N' indeterminates
int BaseScore;
int SumReadScores;
int MeanReadScores;
UINT8 QScore0;
switch(pPE1File->QSSchema) { // quality scoring schema - guestimated scoring schema - 0: no scoring, 1: Solexa, 2: Illumina 1.3+, 3: Illumina 1.5+, 4: Illumina 1.8+ or Sanger 
	case 1:	// Solexa  ';' (-5) to 'h' (40)
		QScore0 = (UINT8)'@';
		break;

	case 2: // Illumina 1.3+ '@' (0) to 'h' (40)
		QScore0 =  (UINT8)'@';
		break;

	case 3: // Illumina 1.5+ 'B' (1) to 'h' (40)
		QScore0 =  (UINT8)'@';
		break;

	case 4: // Illumina 1.8+ or Sanger. Sanger is '!' (0) to 'J' (41) and Illumina is '#' (2) to 'J' (41)
		QScore0 = (UINT8)'!';
		break;
			
	case 0:	// no scoring
		QScore0 = 0;
		break;
	}

SumReadScores = 0;
MeanReadScores = 0;
pBase = pPE1RawRead;
pScore = pPE1QScores;
for (SeqOfs = 0; SeqOfs < PE1ReadLen; SeqOfs++, pBase++, pScore++)
	{
	if(*pBase > eBaseT)
		{
		pPE1File->SeqCharacteristics.NotProcNs += 1;
		break;
		}
	if (pPE1File->QSSchema == 0 || m_MinMeanPhredScore == 0 || PE1QScoreLen == 0 || pPE1QScores == NULL)
		continue;

	BaseScore = (int)*pScore - (int)QScore0;
	if (BaseScore < cMinBasePhredScore)   // not interested if any base in the read has score below cMinBasePhredScore
		{
		pPE1File->SeqCharacteristics.NotProcQS += 1;
		break;
		}
	SumReadScores += BaseScore;
	}
if(SeqOfs == PE1ReadLen && m_MinMeanPhredScore > 0)
	MeanReadScores =  SumReadScores / PE1ReadLen;

if (SeqOfs != PE1ReadLen || MeanReadScores < m_MinMeanPhredScore)
	{
	pPE1File->SeqCharacteristics.NumNoneCntd += 1;
	if (m_bPEProc)
		pPE2File->SeqCharacteristics.NumNoneCntd += 1;
	pPE1File->SeqCharacteristics.NumReads += 1;
	if (m_bPEProc)
		pPE2File->SeqCharacteristics.NumReads += 1;
	return(eBSFSuccess);
	}

if (m_bPEProc)
	{
	SumReadScores = 0;
	MeanReadScores = 0;
	pBase = pPE2RawRead;
	pScore = pPE2QScores;
	for (SeqOfs = 0; SeqOfs < PE2ReadLen; SeqOfs++, pBase++, pScore++)
		{
		if (*pBase > eBaseT)
			{
			pPE2File->SeqCharacteristics.NotProcNs += 1;
			break;
			}
		if (pPE2File->QSSchema == 0 || m_MinMeanPhredScore == 0 || PE2QScoreLen == 0 || pPE2QScores == NULL)
			continue;
		BaseScore = (int)*pScore - (int)QScore0;
		if (BaseScore < cMinBasePhredScore)   // not interested if any base in the read has score below cMinBasePhredScore
			{
			pPE2File->SeqCharacteristics.NotProcQS += 1;
			break;
			}
		SumReadScores += BaseScore;
		}
	if(SeqOfs == PE1ReadLen && m_MinMeanPhredScore > 0)
		MeanReadScores =  SumReadScores / PE2ReadLen;
	if (SeqOfs != PE2ReadLen || MeanReadScores < m_MinMeanPhredScore)
		{
		pPE1File->SeqCharacteristics.NumNoneCntd += 1;
		pPE2File->SeqCharacteristics.NumNoneCntd += 1;
		pPE1File->SeqCharacteristics.NumReads += 1;
		pPE2File->SeqCharacteristics.NumReads += 1;
		return(eBSFSuccess);
		}
	}

NumInsts = AddReadInst(pThread,PE1ReadLen, pPE1RawRead, PE2ReadLen, pPE2RawRead);
if (NumInsts < 0)
	return((teBSFrsltCodes)NumInsts);

if (NumInsts >= 1)
	{
	AccumKMers(pThread,PE1ReadLen, pPE1RawRead);
	if (m_bPEProc)
		AccumKMers(pThread,PE2ReadLen, pPE2RawRead);
	if(LocateContaminentMatch(PE1ReadLen,pPE1RawRead, false))
		bPE1Contaminated = true;
	else
		bPE1Contaminated = false;
	if (m_bPEProc)
		{
		if(LocateContaminentMatch(PE2ReadLen,pPE2RawRead, true))
			bPE2Contaminated = true;
		else
			bPE2Contaminated = false;
		}
	AcquireSerialise();
	m_NumChkdPE1ContamHits += 1;
	if(bPE1Contaminated)
		m_NumPE1ContamHits += 1;
	if (m_bPEProc)
		{
		m_NumChkdPE2ContamHits += 1;
		if(bPE2Contaminated)
			m_NumPE2ContamHits += 1;
		}
	ReleaseSerialise();
	}

if(NumInsts == 1)
	{
	pPE1File->SeqCharacteristics.NumSampled += 1;
	if(m_bPEProc)
		pPE2File->SeqCharacteristics.NumSampled += 1;
	}
if(NumInsts > 0)
	{
	pPE1File->SeqCharacteristics.NumCntd += 1;
	if(m_bPEProc)
		pPE2File->SeqCharacteristics.NumCntd += 1;
	}
else
	{
	pPE1File->SeqCharacteristics.NumNoneCntd += 1;
	if(m_bPEProc)
		pPE2File->SeqCharacteristics.NumNoneCntd += 1;
	}

pPE1File->SeqCharacteristics.NumReads += 1;
if (m_bPEProc)
	pPE2File->SeqCharacteristics.NumReads += 1;
return(eBSFSuccess);
}



