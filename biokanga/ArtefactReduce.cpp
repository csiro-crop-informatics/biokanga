// Copyright 2013 CSIRO  ( http://www.csiro.au/ ) m_bRawDedupe
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

#include "./biokanga.h"
#include "./Kangadna.h"
#include "./ArtefactReduce.h"


int
ProcessArtefactReduce(etARPMode PMode,	// processing mode; currently eAR2Fasta,  eAR2Packed or eARPacked2fasta
	    char *pszCheckpointFile,		// if file of this name exists and is a checkpoint then resume processing from this checkpoint, otherwise create a checkpoint file 
		etSfxSparsity SfxSparsity,		// suffix sparsity
		int IterativePasses,			// iterative passes of overlap processing
		int MinPhredScore,				// only accept reads for filtering if their mean Phred score is at least this threshold
		bool bNoDedupe,					// if true then do not remove all duplicate reads as per bDedupeIndependent parameter
		bool bStrand,					// true if read strand specific filtering
		int MaxNs,						// filter out input sequences having higher than this number of indeterminate bases per 100bp (default is 0, range 0..5)
		int Trim5,						// trim this number of 5' bases from input sequences (default is 0, range 0..20)
		int Trim3,						// trim this number of 3' bases from input sequences (default is 0, range 0..20)
		int MinSeqLen,		            // filter out input sequences (after any trimming) which are less than this length (default is 50bp, range 20..10000)
		int TrimSeqLen,					// trim sequences to be no longer than this length (default is 0 for no length trimming, MinSeqLen...10000)
		int MinOverlap,					// minimum required overlap (in bp) if -1 then no overlap processing
		int MinFlankLen,				// non-overlapping flank must be at least this length (defults to 1, else range 1bp to 25bp, only applies if overlap processing)
		int SampleNth,					// process every Nth reads
		int Zreads,						// maximum number of reads to accept for processing from any file
		bool bDedupeIndependent,		// if paired end preprocessing then treat as if single ended when deuping
		int NumThreads,					// number of worker threads to use
		bool bAffinity,					// thread to core affinity
		int NumPE1InputFiles,			// number of PE1 input files
		char *pszInPE1files[],			// input PE1 5' read files
		int NumPE2InputFiles,			// number of PE2 input files
		char *pszInPE2files[],		    // input PE2 3' read files
		char *pszContaminantFile,		// contaminants fasta file
		char *pszOutFile,				// where to write filtered sequences
		char *pszDupDistFile);			// write duplicate sequence distributions to this file

#ifdef _WIN32
int ArtefactReduce(int argc, char* argv[])
{
// determine my process name
_splitpath(argv[0],NULL,NULL,gszProcName,NULL);
#else
int
ArtefactReduce(int argc, char** argv)
{
// determine my process name
CUtility::splitpath((char *)argv[0],NULL,gszProcName);
#endif
int iScreenLogLevel;		// level of screen diagnostics
int iFileLogLevel;			// level of file diagnostics
char szLogFile[_MAX_PATH];	// write diagnostics to this file
int Rslt = 0;   			// function result code >= 0 represents success, < 0 on failure
int Idx;

int NumberOfProcessors;		// number of installed CPUs
int NumThreads;				// number of threads (0 defaults to number of CPUs)
bool bAffinity;				// thread to core affinity

int PMode;					// processing mode
int MinPhredScore;			// only accept reads for filtering if their mean Phred score is at least this threshold
int SfxSparsity;			// suffix sparsity
bool bStrand;				// if true then strand specific filtering

int MaxNs;					// filter out input sequences having higher than this number of indeterminate bases per 100bp (default is 0, range 0..5)
int Trim5;					// trim this number of 5' bases from input sequences (default is 0, range 0..20)
int Trim3;					// trim this number of 3' bases from input sequences (default is 0, range 0..20)
int MinSeqLen;              // filter out input sequences (after any trimming) which are less than this length (default is 50bp, range 20..500)
int TrimSeqLen;				// trim sequences to be no longer than this length (default is 0 for no length trimming, MinSeqLen...5000)
int MinOverlap;				// minimum required overlap (in bp) ( if -1 then no overlap processing required)
int MinFlankLen;			// non-overlapping flank must be at least this length (defults to 1, range 1bp to 25bp, only applies if overlap processing)
int SampleNth;				// process every Nth reads
int Zreads;					// maximum number of reads to accept for processing from any file 
int IterativePasses;		// iterative passes repeating overlap processing (default 1, range 1..10)

char szContaminantFile[_MAX_PATH];		// contaminants fasta file

bool bNoDedupe;					// do not dedupe exactly matching input sequences - default is to only retain a single copy of duplicate sequences
bool bDedupeIndependent;		// if paired end preprocessing then treat as if single ended when deuping

char szCheckpointFile[_MAX_PATH];	// if file of this name exists and is a checkpoint then resume processing from this checkpoint, otherwise create a checkpoint file
char szOutFile[_MAX_PATH];	// packed and deduped sequences written to this file
char szDupDistFile[_MAX_PATH];	// write duplicate sequence distributions to this file

int NumPE1InputFiles;			// number of PE1 input files
char *pszInPE1files[cKDNAMaxInFileSpecs];  // input PE1 5' read files
int NumPE2InputFiles;			// number of PE2 input files
char *pszInPE2files[cKDNAMaxInFileSpecs];  // input PE2 3' read files

char szSQLiteDatabase[_MAX_PATH];	// results summaries to this SQLite file
char szExperimentName[cMaxDatasetSpeciesChrom+1];			// experiment name
char szExperimentDescr[1000];		// describes experiment

// command line args
struct arg_lit  *help    = arg_lit0("h","help",                 "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to screen and logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"diagnostics log file");

struct arg_int *pmode = arg_int0("m","mode","<int>",		    "processing mode: 0 - filter reads into multifasta, 1 - filter reads into packed sequence file for subsequent de Novo assembly, 2 - load packed sequence file and output as multifasta");
struct arg_int *minphredscore = arg_int0("p", "minphred", "<int>", "only accept reads for filtering if their mean Phred scores are at least this threshold (default 0 for no Phred filtering, range 15..40)");
struct arg_int *iterativepasses = arg_int0("P", "iterativepasses", "<int>", "iterative passes repeating overlap processing (default 1, range 1..5)");

struct arg_file *inpe1files = arg_filen("i","inpe1","<file>",1,cKDNAMaxInFileSpecs,"Load single ended, or 5' if paired end, reads from fasta or fastq file(s); if single ended then wildcards allowed");
struct arg_file *inpe2files = arg_filen("I","inpe2","<file>",0,cKDNAMaxInFileSpecs,"Load 3' if paired end reads from fasta or fastq file(s)");
struct arg_file *outfile = arg_file0("o","out","<file>",		"Output multifasta ('-m0') or packed sequences ('-m1') to this file");
struct arg_file *dupdistfile = arg_file0("O","dupdist","<file>","Output duplicate distributions to this file");

struct arg_file *contaminantfile = arg_file0("c", "contaminants", "<file>", "Putative contaminant sequences fasta file");

struct arg_int *maxns = arg_int0("n","indeterminates","<int>",  "filter out input sequences having higher than this number of indeterminate bases (default is 0, range 0 to 5)");
struct arg_int *trim5 = arg_int0("x","trim5","<int>",			"trim this number of 5' bases from each input sequence (default is 0, range 0..20)");
struct arg_int *trim3 = arg_int0("X","trim3","<int>",			"trim this number of 3' bases from each input sequence (default is 0, range 0..20)");
struct arg_int *trimseqlen = arg_int0("L","trimseqlen","<int>",	"trim sequences to be no longer than this length (default is 0 for no length trimming, MinSeqLen to 1000)");
struct arg_int *minoverlap = arg_int0("y","minoverlap","<int>",	"accept as overlapping if overlaps are at least this length (defaults to 80%% of mean read length, range 25bp to 150bp, if -1 then no overlap processing)");
struct arg_int *minflanklen = arg_int0("Y","minflanklen","<int>","non-overlapping flank must be at least this length (default 1, range 1bp to 25bp, only applies if overlap processing)");

struct arg_int *minseqlen= arg_int0("l","minlen","<int>",       "filter out input sequences (after any trimming) which are less than this length (default is 50bp, range 20..5000)");
struct arg_int *threads = arg_int0("T","threads","<int>",		"number of processing threads 0..64 (defaults to 0 which sets threads to number of cores)");
struct arg_lit  *strand   = arg_lit0("S","strand",              "strand specific filtering - filter reads with read orientation - default is for non-strand specific");
struct arg_lit  *nodedupe   = arg_lit0("D","nodedupe",          "do not dedupe exactly matching input sequences - default is to only retain a single copy of duplicate sequences");
struct arg_lit  *dedupepe   = arg_lit0("d","dedupepe",          "if paired end preprocessing then treat ends as independent when deduping");

struct arg_file *summrslts = arg_file0("q","sumrslts","<file>",		"Output results summary to this SQLite3 database file");
struct arg_str *experimentname = arg_str0("w","experimentname","<str>",		"experiment name SQLite3 database file");
struct arg_str *experimentdescr = arg_str0("W","experimentdescr","<str>",	"experiment description SQLite3 database file");


struct arg_end *end = arg_end(200);

void *argtable[] = {help,version,FileLogLevel,LogFile,
	                pmode,minphredscore,strand,maxns,iterativepasses,trim5,trim3,contaminantfile,minseqlen,trimseqlen,minoverlap,minflanklen,nodedupe,dedupepe,inpe1files,inpe2files,outfile,dupdistfile,
					summrslts,experimentname,experimentdescr,
					threads,
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
		return(1);
		}

	iScreenLogLevel = iFileLogLevel = FileLogLevel->count ? FileLogLevel->ival[0] : eDLInfo;
	if(iFileLogLevel < eDLNone || iFileLogLevel > eDLDebug)
		{
		printf("\nError: FileLogLevel '-l%d' specified outside of range %d..%d\n",iFileLogLevel,eDLNone,eDLDebug);
		return(1);
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
		return(1);
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


	PMode = pmode->count ? pmode->ival[0] : (int)eAR2Fasta;
	if(PMode < 0 || PMode > eARPacked2fasta)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Processing mode '-m%d' must be in range 0..%d",PMode,eARPacked2fasta);
		return(1);
		}

	szCheckpointFile[0] = '\0';
	SfxSparsity = (etSfxSparsity)cDfltSuffixSparsity;
	bStrand = false;
	MinPhredScore = 0;
	strcpy(szDupDistFile,"DupInstDist");
	NumPE2InputFiles = 0;
	NumPE1InputFiles = 0;
	pszInPE1files[0] = NULL;
	pszInPE2files[0] = NULL;
	szContaminantFile[0] = '\0';
	Trim5 = 0;
	Trim3 = 0;
	MaxNs = 0;
	SampleNth = 1;
	Zreads = 0;
	MinSeqLen = 50;
	TrimSeqLen = 0;
	MinOverlap = 0;
	MinFlankLen = 0;
	bDedupeIndependent = false;
	bNoDedupe = false;
	IterativePasses = 1;

	if(PMode != eARPacked2fasta)
		{
		IterativePasses = iterativepasses->count ? iterativepasses->ival[0] : 1;
		if (IterativePasses < 1 || IterativePasses > 5)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: Iterative passes '-P%d' specified outside of range 1..5\n", IterativePasses);
			exit(1);
			}

		MinPhredScore = minphredscore->count ? minphredscore->ival[0] : 0;
		if (MinPhredScore != 0 && (MinPhredScore < 15 || MinPhredScore > 40))
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: Minimum mean Phred score '-p%d' specified outside of range 15..40\n", MinPhredScore);
			exit(1);
			}

 		bStrand = strand->count ? true : false;

		if (contaminantfile->count)
			{
			strncpy(szContaminantFile, contaminantfile->filename[0], _MAX_PATH);
			szContaminantFile[_MAX_PATH - 1] = '\0';
			CUtility::TrimQuotedWhitespcExtd(szContaminantFile);
			}
		if(dupdistfile->count)
			{
			strcpy(szDupDistFile,dupdistfile->filename[0]);
			CUtility::TrimQuotedWhitespcExtd(szDupDistFile);
			}
		if(!inpe1files->count && szCheckpointFile[0]=='\0')
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: No input file(s) specified with with '-i<filespec>' option)");
			return(1);
			}
		
		if(inpe1files->count)
			{
			for(NumPE1InputFiles=Idx=0;NumPE1InputFiles < cKDNAMaxInFileSpecs && Idx < inpe1files->count; Idx++)
				{
				pszInPE1files[Idx] = NULL;
				if(pszInPE1files[NumPE1InputFiles] == NULL)
					pszInPE1files[NumPE1InputFiles] = new char [_MAX_PATH];
				strncpy(pszInPE1files[NumPE1InputFiles],inpe1files->filename[Idx],_MAX_PATH);
				pszInPE1files[NumPE1InputFiles][_MAX_PATH-1] = '\0';
				CUtility::TrimQuotedWhitespcExtd(pszInPE1files[NumPE1InputFiles]);
				if(pszInPE1files[NumPE1InputFiles][0] != '\0')
					NumPE1InputFiles++;
				}

			if(!NumPE1InputFiles)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: After removal of whitespace, no input file(s) specified with '-i<filespec>' option");
				return(1);
				}
			}
	
		if(NumPE1InputFiles > 0 && inpe2files->count)
			{
			for(NumPE2InputFiles=Idx=0;NumPE2InputFiles < cKDNAMaxInFileSpecs && Idx < inpe2files->count; Idx++)
				{
				pszInPE2files[Idx] = NULL;
				if(pszInPE2files[NumPE2InputFiles] == NULL)
					pszInPE2files[NumPE2InputFiles] = new char [_MAX_PATH];
				strncpy(pszInPE2files[NumPE2InputFiles],inpe2files->filename[Idx],_MAX_PATH);
				pszInPE2files[NumPE2InputFiles][_MAX_PATH-1] = '\0';
				CUtility::TrimQuotedWhitespcExtd(pszInPE2files[NumPE2InputFiles]);
				if(pszInPE2files[NumPE2InputFiles][0] != '\0')
					NumPE2InputFiles++;
				}
			if(!NumPE2InputFiles)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: After removal of whitespace, no PE2 3' input file(s) specified with '-I<filespec>' option");
				return(1);
				}

			if(NumPE1InputFiles != NumPE2InputFiles)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Paired end processing, expected number of PE1 files (%d) to equal number of PE2 files (%d)",NumPE1InputFiles,NumPE2InputFiles);
				return(1);
				}
			}

		if(NumPE1InputFiles > 0)
			{
			Trim5 = trim5->count ? trim5->ival[0] : 0;
			if(Trim5 < 0 || Trim5 > 20)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: 5' trim '-x%d' must be in range 0..20 bases",Trim5);
				return(1);
				}

			Trim3 = trim3->count ? trim3->ival[0] : 0;
			if(Trim3 < 0 || Trim3 > 20)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: 3' trim '-x%d' must be in range 0..20 bases",Trim3);
				return(1);
				}


			MaxNs = maxns->count ? maxns->ival[0] : 0;
			if(MaxNs < 0 || MaxNs > 5)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: max percentage of indeterminate bases '-n%d' must be in range 0..5",MaxNs);
				return(1);
				}

			SampleNth = 1;
			Zreads = 0;

			MinSeqLen = minseqlen->count ? minseqlen->ival[0] : 50;
			if(MinSeqLen < 20 || MinSeqLen > 10000)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Minimum sequence length '-l%d' must be in range 30..10000 bases",MinSeqLen);
				return(1);
				}

			TrimSeqLen = trimseqlen->count ? trimseqlen->ival[0] : 0;
			if(TrimSeqLen != 0 &&  (TrimSeqLen < MinSeqLen || TrimSeqLen > 1000))
				{
				if(TrimSeqLen < MinSeqLen && MinSeqLen > 20)
					gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Trimming sequence length '-L%d' must be in range %d..1000 bases",TrimSeqLen,MinSeqLen);
				else
					gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Trimming sequence length '-L%d' must be in range 20..1000 bases",TrimSeqLen);
				return(1);
				}
			}
	
		MinOverlap = minoverlap->count ? minoverlap->ival[0] : 0;
		if(!(MinOverlap == 0 || MinOverlap == -1) && (MinOverlap < 25 || MinOverlap > 150))
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Minimum overlap '-y%d' must be in range 25..150 or 0 for auto or -1 to disable overlap processing",MinOverlap);
			return(1);
			}

		if(MinOverlap != -1) // -1 if no overlap processing
			{
			MinFlankLen = minflanklen->count ? minflanklen->ival[0] : 1;
			if(MinFlankLen < 1 || MinFlankLen > 25)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Minimum flanking length '-Y%d' must be in range 1..25",MinFlankLen);
				return(1);
				}
			}
		else
			MinFlankLen = 0;

		bNoDedupe = nodedupe->count ? true : false;
		if(bNoDedupe)
			bDedupeIndependent = false;
		else
			bDedupeIndependent = dedupepe->count ? true : false;
		}
	else
		{
		if(!inpe1files->count)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: No input packed file specified with with '-i<filespec>' option)");
			return(1);
			}
		pszInPE1files[0] = new char [_MAX_PATH];
		strncpy(pszInPE1files[0],inpe1files->filename[0],_MAX_PATH);
		pszInPE1files[0][_MAX_PATH-1] = '\0';
		CUtility::TrimQuotedWhitespcExtd(pszInPE1files[0]);
		if(pszInPE1files[0][0] == '\0')
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: After removal of whitespace, no input packed sequence file(s) specified with '-i<filespec>' option\n");
			return(1);
			}
		}

	if(!outfile->count)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Expected output file '-o<outfile>' to be specified");
		return(1);
		}
	strncpy(szOutFile,outfile->filename[0],_MAX_PATH);
	szOutFile[_MAX_PATH-1] = '\0';


#ifdef _WIN32
	SYSTEM_INFO SystemInfo;
	GetSystemInfo(&SystemInfo);
	NumberOfProcessors = SystemInfo.dwNumberOfProcessors;
	bAffinity = false;
#else
	bAffinity = false;
	NumberOfProcessors = sysconf(_SC_NPROCESSORS_ONLN);
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

	if(bAffinity)
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
		for(AffinityIdx=0; AffinityIdx < CPU_SETSIZE; AffinityIdx++)	// determine number of CPUs to which affinity is enabled
			{
			if(CPU_ISSET(AffinityIdx,&CurCpuSet))
				CoresPerCPU += 1;
			}
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Info: %d CPUs (%d cores) available",CoresPerCPU, NumberOfProcessors);
		// number of CPUs with affinity known, estimate cores per CPU

		CoresPerCPU = NumberOfProcessors/CoresPerCPU;
	
		// now cores per CPU has been guestimated then try and set affinity for CPU's
		ReqAffinity = 0;
		for(AffinityIdx=0; AffinityIdx < CPU_SETSIZE; AffinityIdx++)
			{
			if(CPU_ISSET(AffinityIdx,&CurCpuSet))
				{
				CPU_SET(AffinityIdx,&NewCpuSet);
				ReqCores -= CoresPerCPU;
				ReqAffinity += 1;
				if(ReqCores < 1)
					break;
				}
			}
		pthread_setaffinity_np(Me,sizeof(NewCpuSet), &NewCpuSet);
#endif
		}

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing parameters:");
	const char *pszDescr;
	switch((etARPMode)PMode) {
		case eAR2Fasta:
			pszDescr = "Load, parse, filter raw reads files into multifasta files";
			break;
		case eAR2Packed:
			pszDescr = "Load, parse and filter raw reads files into packed sequence file for assembly";
			break;
		case eARPacked2fasta:
			pszDescr = "Load a previously generated packed and output as multifasta file(s)";
		}

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Processing mode is : '%s'",pszDescr);

	if(PMode != eARPacked2fasta)
		{
		gDiagnostics.DiagOutMsgOnly(eDLInfo, "Min mean Phred threshold score for processing acceptance : %d", MinPhredScore);
	
		gDiagnostics.DiagOutMsgOnly(eDLInfo, "Iterative overlap processing passes : %d", IterativePasses);
	

		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Filter out duplicate sequences : '%s'",bNoDedupe ? "No" : "Yes");

		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Strand assembly : '%s'",bStrand ? "dependent" : "independent");

		if (szContaminantFile[0] != '\0')
			gDiagnostics.DiagOutMsgOnly(eDLInfo, "Contaminant sequences file: %s", szContaminantFile);

		if(NumPE1InputFiles > 0)
			{
			if(SampleNth > 1)
				gDiagnostics.DiagOutMsgOnly(eDLInfo,"Process every : %d reads",SampleNth);

			if(Zreads > 0)
				gDiagnostics.DiagOutMsgOnly(eDLInfo,"Processing up to: %d reads per file",Zreads);
	
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Accept input sequences if percentage of indeterminate bases (Ns) no more than : %d",MaxNs);
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Trim 5' input sequences by : %d",Trim5);
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Trim 3' input sequences by : %d",Trim3);
			if(TrimSeqLen)
				gDiagnostics.DiagOutMsgOnly(eDLInfo,"Trim input sequence length to be no more than : %d",TrimSeqLen);
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Accept input sequences (after any trim) if at least this length : %d",MinSeqLen);
			}

		if(MinOverlap == 0)
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Minimum sequence overlap: Use 80%% of mean read length");
		else
			{
			if(MinOverlap == -1)
				gDiagnostics.DiagOutMsgOnly(eDLInfo,"Minimum sequence overlaps: No overlap processing");
			else
				gDiagnostics.DiagOutMsgOnly(eDLInfo,"Minimum sequence flank overlap: %d bp",MinOverlap);
			}
		if(MinOverlap != -1)
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Non-overlapping flank must be at least this length: %dbp",MinFlankLen);

		if(NumPE2InputFiles)
			{
			if(!bNoDedupe)
				gDiagnostics.DiagOutMsgOnly(eDLInfo,"Independently dedupe paired ends: %s",bDedupeIndependent ? "Yes" : "No");
			for(Idx=0; Idx < NumPE1InputFiles; Idx++)
				gDiagnostics.DiagOutMsgOnly(eDLInfo,"Input PE 5' reads file (%d) : '%s'",Idx+1,pszInPE1files[Idx]);

			for(Idx=0; Idx < NumPE2InputFiles; Idx++)
				gDiagnostics.DiagOutMsgOnly(eDLInfo,"Input PE 3' reads file (%d) : '%s'",Idx+1,pszInPE2files[Idx]);
			if(PMode == eAR2Fasta)
				{
				gDiagnostics.DiagOutMsgOnly(eDLInfo,"Output PE1 multifasta file : '%s%s'",szOutFile,".R1.fasta");
				gDiagnostics.DiagOutMsgOnly(eDLInfo,"Output PE2 multifasta file : '%s%s'",szOutFile,".R2.fasta");
				}
			else
				gDiagnostics.DiagOutMsgOnly(eDLInfo,"Output packed sequence file : '%s'",szOutFile);
			}
		else
			{
			for(Idx=0; Idx < NumPE1InputFiles; Idx++)
				gDiagnostics.DiagOutMsgOnly(eDLInfo,"Input single ended reads file (%d) : '%s'",Idx+1,pszInPE1files[Idx]);
			if(PMode == eAR2Fasta)
				gDiagnostics.DiagOutMsgOnly(eDLInfo,"Output SE multifasta file : '%s'",szOutFile);
			else
				gDiagnostics.DiagOutMsgOnly(eDLInfo,"Output packed sequence file : '%s'",szOutFile);
			}
	
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Duplicate sequence instances distribution file name prefix : '%s'",szDupDistFile);
		}
	else
		{
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Input packed reads file : '%s'",pszInPE1files[0]);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Output multifasta file prefix : '%s'",szOutFile);
		}

	if(szExperimentName[0] != '\0')
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"This processing reference: %s",szExperimentName);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Number of threads : %d",NumThreads);

	if(gExperimentID > 0)
		{
		int ParamID;
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szLogFile),"log",szLogFile);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,sizeof(PMode),"mode",&PMode);
		if(PMode != eARPacked2fasta)
			{
			ParamID = gSQLiteSummaries.AddParameter(gProcessingID, ePTInt32, (int)sizeof(MinPhredScore), "minphredscore", &MinPhredScore);
			ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,sizeof(MaxNs),"indeterminates",&MaxNs);
			ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,sizeof(IterativePasses),"passes",&IterativePasses);
			ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,sizeof(Trim5),"trim5",&Trim5);
			ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,sizeof(Trim3),"trim3",&Trim3);

			if (szContaminantFile[0] != '\0')
				ParamID = gSQLiteSummaries.AddParameter(gProcessingID, ePTText, (int)strlen(szContaminantFile), "contaminantfile", szContaminantFile);
			ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,sizeof(TrimSeqLen),"trimseqlen",&TrimSeqLen);
			ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,sizeof(MinOverlap),"minoverlap",&MinOverlap);
			ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,sizeof(MinFlankLen),"minflanklen",&MinFlankLen);
			ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,sizeof(MinSeqLen),"minlen",&MinSeqLen);
			ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTBool,sizeof(bStrand),"strand",&bStrand);
			ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTBool,sizeof(bNoDedupe),"nodedupe",&bNoDedupe);
			ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTBool,sizeof(bDedupeIndependent),"dedupepe",&bDedupeIndependent);

			ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,sizeof(NumPE1InputFiles),"NumPE1InputFiles",&NumPE1InputFiles);
			for(Idx=0; Idx < NumPE1InputFiles; Idx++)
				ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(pszInPE1files[Idx]),"inpe1",pszInPE1files[Idx]);

			ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,sizeof(NumPE2InputFiles),"NumPE2InputFiles",&NumPE2InputFiles);
			for(Idx=0; Idx < NumPE2InputFiles; Idx++)
				ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(pszInPE2files[Idx]),"inpe2",pszInPE2files[Idx]);

			ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szOutFile),"out",szOutFile);
			ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szDupDistFile),"dupdist",szDupDistFile);
			}
		else
			{
			ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(pszInPE1files[0]),"inpe1",pszInPE1files[0]);
			ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szOutFile),"out",szOutFile);
			}

		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,sizeof(NumThreads),"threads",&NumThreads);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,sizeof(NumberOfProcessors),"cpus",&NumberOfProcessors);

		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szSQLiteDatabase),"sumrslts",szSQLiteDatabase);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szExperimentName),"experimentname",szExperimentName);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szExperimentDescr),"experimentdescr",szExperimentDescr);
		}

#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	gStopWatch.Start();
	Rslt = ProcessArtefactReduce((etARPMode)PMode,szCheckpointFile,(etSfxSparsity)SfxSparsity,IterativePasses,MinPhredScore,bNoDedupe,bStrand,MaxNs,Trim5,Trim3, MinSeqLen,TrimSeqLen,MinOverlap,MinFlankLen,SampleNth,Zreads,bDedupeIndependent,NumThreads,bAffinity,
							NumPE1InputFiles,pszInPE1files,NumPE2InputFiles,pszInPE2files,szContaminantFile, szOutFile, szDupDistFile);
	Rslt = Rslt >=0 ? 0 : 1;
	if(gExperimentID > 0)
		{
		if(gProcessingID)
			gSQLiteSummaries.EndProcessing(gProcessingID,Rslt);
		gSQLiteSummaries.EndExperiment(gExperimentID);
		}
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Exit code: %d Total processing time: %s",Rslt,gStopWatch.Read());
	return(Rslt);
	}
else
	{
	printf("\n%s %s %s, Version %s\n", gszProcName,gpszSubProcess->pszName,gpszSubProcess->pszFullDescr,cpszProgVer);
	arg_print_errors(stdout,end,gszProcName);
	arg_print_syntax(stdout,argtable,"\nUse '-h' to view option and parameter usage\n");
	return(1);
	}
return 0;
}

int
ProcessArtefactReduce(etARPMode PMode,	// processing mode, currently eAR2Fasta,  eAR2Packed
	    char *pszCheckpointFile,		// if file of this name exists and is a checkpoint then resume processing from this checkpoint, otherwise create a checkpoint file 
		etSfxSparsity SfxSparsity,		// suffix sparsity
		int IterativePasses,			// iterative passes of overlap processing
		int MinPhredScore,				// only accept reads for filtering if their mean Phred score is at least this threshold
		bool bNoDedupe,					// if true then do not remove all duplicate reads as per bDedupeIndependent parameter
		bool bStrand,					// true if read strand specific filtering
		int MaxNs,						// filter out input sequences having higher than this number of indeterminate bases per 100bp (default is 0, range 0..5)
		int Trim5,						// trim this number of 5' bases from input sequences (default is 0, range 0..20)
		int Trim3,						// trim this number of 3' bases from input sequences (default is 0, range 0..20)
		int MinSeqLen,		            // filter out input sequences (after any trimming) which are less than this length (default is 50bp, range 20..10000)
		int TrimSeqLen,					// trim sequences to be no longer than this length (default is 0 for no length trimming, MinSeqLen...10000)
		int MinOverlap,					// minimum required overlap (in bp) or -1 if no overlap processing
		int MinFlankLen,				// non-overlapping flank must be at least this length (defults 1, else range 1bp to 25bp, only applies if overlap processing)
		int SampleNth,					// process every Nth reads
		int Zreads,						// maximum number of reads to accept for processing from any file
		bool bDedupeIndependent,		// if paired end preprocessing then treat as if single ended when deuping
		int NumThreads,					// number of worker threads to use
		bool bAffinity,					// thread to core affinity
		int NumPE1InputFiles,			// number of PE1 input files
		char *pszInPE1files[],			// input PE1 5' read files
		int NumPE2InputFiles,			// number of PE2 input files
		char *pszInPE2files[],		    // input PE2 3' read files
		char *pszContaminantFile,		// contaminants fasta file
		char *pszOutFile,				// where to write filtered sequences
		char *pszDupDistFile)			// write duplicate sequence distributions to this file
{
int Rslt;
CArtefactReduce *pArtefactReduce;

if((pArtefactReduce = new CArtefactReduce) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate instance of CArtefactReduce");
	return(-1);
	}

Rslt = pArtefactReduce->Process(PMode,pszCheckpointFile,SfxSparsity,IterativePasses,MinPhredScore,bNoDedupe,bStrand,MaxNs,Trim5,Trim3,MinSeqLen,TrimSeqLen,MinOverlap,MinFlankLen,
								SampleNth,Zreads,bDedupeIndependent,NumThreads,	bAffinity, NumPE1InputFiles,pszInPE1files,NumPE2InputFiles,pszInPE2files,pszContaminantFile,pszOutFile,pszDupDistFile);

delete pArtefactReduce;
return(Rslt);
}


#ifdef _WIN32
unsigned __stdcall ThreadedIdentDups(void * pThreadPars)
#else
void * ThreadedIdentDups(void * pThreadPars)
#endif
{
int Rslt = 0;
tsThreadIdentDuplicatePars *pPars = (tsThreadIdentDuplicatePars *)pThreadPars; // makes it easier not having to deal with casts!
CArtefactReduce *pThis = (CArtefactReduce *)pPars->pThis;
Rslt = pThis->ProcIdentDuplicates(pPars);
pPars->Rslt = Rslt;
#ifdef _WIN32
_endthreadex(0);
return(eBSFSuccess);
#else
pthread_exit(NULL);
#endif
}

#ifdef _WIN32
unsigned __stdcall ThreadedIdentOverlaps(void * pThreadPars)
#else
void * ThreadedIdentOverlaps(void * pThreadPars)
#endif
{
int Rslt = 0;
tsThreadIdentOverlapPars *pPars = (tsThreadIdentOverlapPars *)pThreadPars; // makes it easier not having to deal with casts!
CArtefactReduce *pThis = (CArtefactReduce *)pPars->pThis;
Rslt = pThis->ProcIdentOverlaps(pPars);
pPars->Rslt = Rslt;
#ifdef _WIN32
_endthreadex(0);
return(eBSFSuccess);
#else
pthread_exit(NULL);
#endif
}

// generally relies on base classes constructors
CArtefactReduce::CArtefactReduce(void)
{
ARInit();
}

// generally relies on base classes destructors
CArtefactReduce::~CArtefactReduce(void)
{
ARReset();
}

void
CArtefactReduce::ARInit(void)
{
m_pKMerSeqHashes = NULL; 
m_pKMerSeqs = NULL;
ARReset();
}

void
CArtefactReduce::ARReset(void) 
{
if(m_pKMerSeqHashes != NULL)
	{
	delete m_pKMerSeqHashes;
	m_pKMerSeqHashes = NULL;
	}
if(m_pKMerSeqs != NULL)
	{
#ifdef _WIN32
	free(m_pKMerSeqs);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pKMerSeqs != MAP_FAILED)
		munmap(m_pKMerSeqs,m_AllocdKMerSeqInstsMem);
#endif
	m_pKMerSeqs = NULL;
	}

m_KMerSeqLen = 0;
m_KMerSeqInstSize = 0;
m_NumKMerSeqHashes = 0;
m_UsedKMerSeqInsts = 0; 
m_AllocdKMerSeqInsts = 0;
m_AllocdKMerSeqInstsMem = 0;

}

int
CArtefactReduce::Process(etARPMode PMode,	// processing mode, currently eAR2Fasta,  eAR2Packed or eARPacked2fasta
		char *pszCheckpointFile,		// if file of this name exists and is a checkpoint then resume processing from this checkpoint, otherwise create a checkpoint file 
		etSfxSparsity SfxSparsity,		// suffix sparsity
		int IterativePasses,			// iterative passes of overlap processing
		int MinPhredScore,				// only accept reads for filtering if their mean Phred score is at least this threshold
		bool bNoDedupe,					// if true then do not remove all duplicate reads as per bDedupeIndependent parameter
		bool bStrand,					// true if read strand specific assembly
		int MaxNs,						// filter out input sequences having higher than this number of indeterminate bases per 100bp (default is 0, range 0..5)
		int Trim5,						// trim this number of 5' bases from input sequences (default is 0, range 0..20)
		int Trim3,						// trim this number of 3' bases from input sequences (default is 0, range 0..20)
		int MinSeqLen,		            // filter out input sequences (after any trimming) which are less than this length (default is 50bp, range 20..10000)
		int TrimSeqLen,					// trim sequences to be no longer than this length (default is 0 for no length trimming, MinSeqLen...10000)
		int MinOverlap,					// minimum required overlap (in bp) or -1 if no overlap processing
		int MinFlankLen,				// non-overlapping flank must be at least this length (defults to 1, else range 1bp to 25bp, only applies if overlap processing)
		int SampleNth,					// process every Nth reads
		int Zreads,						// maximum number of reads to accept for processing from any file
		bool bDedupeIndependent,		// if paired end preprocessing then treat as if single ended when deuping
		int NumThreads,					// number of worker threads to use
		bool bAffinity,					// thread to core affinity
		int NumPE1InputFiles,			// number of PE1 input files
		char *pszInPE1files[],			// input PE1 5' read files
		int NumPE2InputFiles,			// number of PE2 input files
		char *pszInPE2files[],		    // input PE2 3' read files
		char *pszContaminantFile,		// contaminants fasta file
		char *pszOutFile,				// where to write packed sequences
		char *pszDupDistFile)			// write duplicate sequence distributions to this file
{
int Rslt;
int Idx;
int NumInputFilesProcessed;
UINT32 TotNumSeqsAccepted;
UINT32 NumReads;
UINT32 TotNumSEReads;
UINT32 TotNumPEReads;
UINT32 NumPE1Reads;
UINT32 NumPE2Reads;

int SeqWrdBytes;
UINT64 CumulativeMemory;
UINT32 CumulativeSequences;
char *pszInFile;

char szPEDupDistFile[_MAX_PATH];

if(pszDupDistFile == NULL || pszDupDistFile[0] == '\0')
	szPEDupDistFile[0]  = '\0';
else
	{
	strcpy(szPEDupDistFile,pszDupDistFile);
	strcat(szPEDupDistFile,(char *)".dist.csv");
	}

Reset(false);
SetPMode(PMode);
SetNumThreads(NumThreads,bAffinity);
SetDedupePE(bDedupeIndependent);
SetSfxSparsity(SfxSparsity);

if(PMode == eARPacked2fasta)
	{
	if((Rslt = LoadPackedSeqsFromFile(pszInPE1files[0])) < eBSFSuccess)
		return((teBSFrsltCodes)Rslt);
	Rslt = SaveAssembSeqs(pszOutFile);
	return((teBSFrsltCodes)Rslt);
	}

NumInputFilesProcessed = 0;
TotNumSeqsAccepted = 0;
TotNumSEReads = 0;
TotNumPEReads = 0;
CumulativeMemory = 0;

SeqWrdBytes = GetSeqWrdBytes();

Rslt = eBSFerrOpnFile;		// assumes no check point resumption
if(pszCheckpointFile != NULL && pszCheckpointFile[0] != '\0')
	{
	int hFile = open(pszCheckpointFile,O_READSEQ);		// if file exists and can be opened for reading then worth trying to load...
	if(hFile != -1)
		{
		close(hFile);
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Checkpoint file '%s' exists, loading...",pszCheckpointFile);
		Rslt = LoadPackedSeqsFromFile(pszCheckpointFile);
		if(Rslt >= eBSFSuccess)
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"Checkpoint file '%s' load completed",pszCheckpointFile);
		else
			{
			if(NumPE1InputFiles)
				gDiagnostics.DiagOut(eDLInfo,gszProcName,"Unable to load from checkpoint file '%s', will load and process raw sequence datasets",pszCheckpointFile);
			else
				{
				gDiagnostics.DiagOut(eDLInfo,gszProcName,"Unable to load from checkpoint file '%s', and no raw sequence files specified, nothing to do...",pszCheckpointFile);
				return(eBSFerrOpnFile);	// treat as though unable to open file
				}
			}
		}	

	}	

if(Rslt < eBSFSuccess)
	{
	// if requested then load and prepare for contaminant sequence processing
	if(pszContaminantFile != NULL && pszContaminantFile[0] != '\0')
		{
		if((m_pContaminants = new CContaminants)==NULL)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Unable to instantiate CContaminants");
			Reset();
			return(eBSFerrObj);
			}
		if ((Rslt=m_pContaminants->LoadContaminantsFile(pszContaminantFile)) < eBSFSuccess)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Unable to load putative contaminants from file '%s", pszContaminantFile);
			Reset();
			return(Rslt);
			}
		}

	CSimpleGlob glob(SG_GLOB_FULLSORT);

	// estimate mean read lengths plus total number of reads in each file so can more efficently allocate memory as a single block...
	for(Idx = 0; Idx < NumPE1InputFiles; Idx++)
		{
		glob.Init();
		if(glob.Add(pszInPE1files[Idx]) < SG_SUCCESS)
			{
			gDiagnostics.DiagOut(eDLWarn,gszProcName,"Unable to glob '%s",pszInPE1files[Idx]);
			return(eBSFerrOpnFile);	// treat as though unable to open file
			}

		if(glob.FileCount() <= 0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to locate any source reads file matching '%s",pszInPE1files[Idx]);
			continue;
			}
		Rslt = eBSFSuccess;
		for (int FileID = 0; Rslt >= eBSFSuccess &&  FileID < glob.FileCount(); ++FileID)
			{
			pszInFile = glob.File(FileID);
			if((Rslt = EstMemReq(pszInFile)) != eBSFSuccess)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to estimate memory required for '%s",pszInFile);
				return(eBSFerrOpnFile);	// treat as though unable to open file
				}
			}
		}


	if(NumPE2InputFiles)
		{
		for(Idx = 0; Idx < NumPE2InputFiles; Idx++)
			{
			glob.Init();
			if(glob.Add(pszInPE2files[Idx]) < SG_SUCCESS)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to glob '%s",pszInPE2files[Idx]);
				return(eBSFerrOpnFile);	// treat as though unable to open file
				}

			if(glob.FileCount() <= 0)
				{
				gDiagnostics.DiagOut(eDLWarn,gszProcName,"Unable to locate any source reads file matching '%s",pszInPE2files[Idx]);
				continue;
				}
			Rslt = eBSFSuccess;
			for (int FileID = 0; Rslt >= eBSFSuccess &&  FileID < glob.FileCount(); ++FileID)
				{
				pszInFile = glob.File(FileID);
				if((Rslt = EstMemReq(pszInFile)) != eBSFSuccess)
					{
					gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to estimate memory required for '%s",pszInFile);
					return(eBSFerrOpnFile);	// treat as though unable to open file
					}
				}
			}
		}

	CumulativeMemory = EstMemReq(Trim5,Trim3,TrimSeqLen,SampleNth,Zreads);
	CumulativeSequences = GetEstSeqsToProc();
	if(CumulativeSequences < cMinSeqs2Assemb)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to estimate number of sequences to filter; require at least 100 sequences to process");
		return(eBSFerrOpnFile);	// treat as though unable to open file
		}

	CumulativeMemory += CumulativeSequences * 12; // very rough estimate allowing for sparse suffix and flags requirements
	if(CumulativeMemory < 1000000000)	// 100M
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Estimated total cumulative minimum required memory: %1.1f MB",(double)CumulativeMemory /1000000);
	else
		{
		if(CumulativeMemory < 1000000000) // 1000M
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"Estimated total cumulative minimum required memory: %1.0f MB",(double)CumulativeMemory /1000000);
		else
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"Estimated total cumulative minimum required memory: %1.1f GB",(double)CumulativeMemory /1000000000);
		}

	// allocate to hold est cumulative memory upfront - may as well know now rather than later if there is insufficent memory...
	if((Rslt = AllocSeqs2AssembMem((CumulativeMemory * 110)/100))!= eBSFSuccess)	// add 10% , reduces chances of having to later realloc...
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to continue");
		Reset(false);
		return(Rslt);
		}

	// now load the PE1 raw read sequences applying any trimming and Phred score filtering 
	// if no PE2 files to load then can allow wildcards, if PE2 to load then can't allow wildcards as
	// wouldn't be able to reliably associate PE1 reads with the PE2 reads
	if(!NumPE2InputFiles)
		{
		for(Idx = 0; Idx < NumPE1InputFiles; Idx++)
			{
			glob.Init();
			if(glob.Add(pszInPE1files[Idx]) < SG_SUCCESS)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to glob '%s",pszInPE1files[Idx]);
				return(eBSFerrOpnFile);	// treat as though unable to open file
				}

			if(glob.FileCount() <= 0)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to locate any source reads file matching '%s",pszInPE1files[Idx]);
				continue;
				}
			Rslt = eBSFSuccess;
			for (int FileID = 0; Rslt >= eBSFSuccess &&  FileID < glob.FileCount(); ++FileID)
				{
				pszInFile = glob.File(FileID);
				NumInputFilesProcessed += 1;
				gDiagnostics.DiagOut(eDLInfo,gszProcName,"Process: Parsing and filtering single ended reads from input read file '%s'",pszInFile);
				if((Rslt = LoadReads(MaxNs,MinPhredScore,Trim5,Trim3,MinSeqLen,TrimSeqLen,SampleNth,Zreads,pszInFile,NULL)) < eBSFSuccess)
					{
					gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: Load failed for input sequences file '%s'",pszInFile);
					Reset(false);
					return((teBSFrsltCodes)Rslt);
					}
				NumReads = (UINT32)Rslt;
				TotNumSEReads += NumReads;
				if(NumReads == 0)
					gDiagnostics.DiagOut(eDLWarn,gszProcName,"Process: No sequences accepted from input single ended reads file '%s'",pszInFile);
				else
					gDiagnostics.DiagOut(eDLInfo,gszProcName,"Process: Accepted %u sequences from single ended reads file '%s'",NumReads,pszInFile);
				}
			}
		if(TotNumSEReads < 1)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: No single end reads accepted from input reads files - nothing to assemble");
			Reset(true);
			return(eBSFerrNoEntries);
			}
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Process: Accepted total of %u reads from input single ended sequences files",TotNumSEReads);
		}
	else // else must be paired end processing
		{
		for(Idx = 0; Idx < NumPE1InputFiles; Idx++)
			{
			NumInputFilesProcessed += 2;
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"Process: Parsing and filtering paired end reads from input reads files '%s' and '%s'",pszInPE1files[Idx], pszInPE2files[Idx]);
			if((Rslt = LoadReads(MaxNs,MinPhredScore,Trim5,Trim3,MinSeqLen,TrimSeqLen,SampleNth,Zreads,pszInPE1files[Idx],pszInPE2files[Idx])) < eBSFSuccess)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: Load failed for paired end sequences files '%s' and '%s'",pszInPE1files[Idx], pszInPE2files[Idx]);
				Reset(false);
				return((teBSFrsltCodes)Rslt);
				}
			NumReads = (UINT32)Rslt;
			TotNumPEReads += NumReads;
			if(NumReads == 0)
				gDiagnostics.DiagOut(eDLWarn,gszProcName,"Process: No reads accepted from input paired end files '%s' and '%s'",pszInPE1files[Idx], pszInPE2files[Idx]);
			else
				gDiagnostics.DiagOut(eDLInfo,gszProcName,"Process: Accepted %u paired ends from input paired end files '%s' and '%s'",NumReads/2,pszInPE1files[Idx], pszInPE2files[Idx]);
			}

		if(TotNumPEReads < 2)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: No paired end reads accepted from input reads files - nothing to assemble");
			Reset(true);
			return(eBSFerrNoEntries);
			}
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Process: Accepted total of %u paired ends (%u sequences) from input reads files",TotNumPEReads/2,TotNumPEReads);
		}


	// if user has requested it then save checkpoint to file
	if(pszCheckpointFile != NULL && pszCheckpointFile[0] != '\0')
		{
		Rslt = SavePackedSeqsToFile(pszCheckpointFile);
		if(Rslt < eBSFSuccess)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: Unable to create checkpoint point file '%s'",pszCheckpointFile);
			Reset(true);
			return(Rslt);
			}
		}
	}

// now identify and remove all exact duplicates retaining a single copy of each
UINT32 TotSeqsParsed;
UINT32 TotSeqsUnderlength;
UINT32 TotSeqsExcessNs;
UINT32 MeanSeqLen;
GetNumReads(&NumPE1Reads,&NumPE2Reads,NULL,NULL,&TotSeqsParsed,&TotSeqsUnderlength,&TotSeqsExcessNs,&MeanSeqLen);

if(gProcessingID > 0)
	{
	gSQLiteSummaries.AddResult(gProcessingID,(char *)"SEReads",ePTUint32,sizeof(TotSeqsParsed),"Parsed",&TotSeqsParsed);
	gSQLiteSummaries.AddResult(gProcessingID,(char *)"SEReads",ePTUint32,sizeof(TotSeqsUnderlength),"Underlength",&TotSeqsUnderlength);
	gSQLiteSummaries.AddResult(gProcessingID,(char *)"SEReads",ePTUint32,sizeof(TotSeqsExcessNs),"ExcessNs",&TotSeqsExcessNs);
	if(!NumPE2Reads)
		gSQLiteSummaries.AddResult(gProcessingID,(char *)"SEReadsLoaded",ePTUint32,sizeof(NumPE1Reads),"Cnt",&NumPE1Reads);
	else
		{
		gSQLiteSummaries.AddResult(gProcessingID,(char *)"PE1ReadsLoaded",ePTUint32,sizeof(NumPE1Reads),"Cnt",&NumPE1Reads);
		gSQLiteSummaries.AddResult(gProcessingID,(char *)"PE2ReadsLoaded",ePTUint32,sizeof(NumPE2Reads),"Cnt",&NumPE2Reads);
		}
	}

if(NumPE1Reads < cMinSeqs2Assemb)				// arbitary lower limit on number of reads required to continue processing
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: Unable to continue, only %d sequences loaded, require at least %d",NumPE1Reads,cMinSeqs2Assemb);
	if(gProcessingID > 0)
		gSQLiteSummaries.AddLog(gProcessingID,"Unable to continue, only %d sequences loaded, require at least %d",NumPE1Reads,cMinSeqs2Assemb);

	Reset(true);
	return(eBSFerrFastqSeq);
	}

if(MinOverlap == 0)
	MinOverlap = (79+(MeanSeqLen * 80))/100;
if(MinOverlap != -1)
	{
	if(MinOverlap < 25)
		MinOverlap = 25;
	else
		if(MinOverlap > 150)
			MinOverlap = 150;
	}
if(MinOverlap != -1)
	{
	if(MinFlankLen < 1)
		MinFlankLen = 1;
	else
		if(MinFlankLen > 25)
			MinFlankLen = 25;
	}
else
	MinFlankLen = 0;

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Process: Mean sequence length was %d, using minimum required overlap length of %d and non-overlap flank length of %d",MeanSeqLen,MinOverlap,MinFlankLen);

GenSeqStarts(true);
GenRdsSfx(1);			// first SeqWrd only requires indexing

if(!bNoDedupe)
	{
	RemoveDuplicates(NumPE2Reads > 0 ? true : false,bStrand,szPEDupDistFile);
	FreeSfx();
	FreeSeqStarts();
	GetNumReads(NULL,NULL,&NumPE1Reads,&NumPE2Reads);

	if(gProcessingID > 0)
		{
		if(!NumPE2Reads)
			gSQLiteSummaries.AddResult(gProcessingID,(char *)"SEReadsDeduped",ePTUint32,sizeof(NumPE1Reads),"Cnt",&NumPE1Reads);
		else
			{
			gSQLiteSummaries.AddResult(gProcessingID,(char *)"PE1ReadsDeduped",ePTUint32,sizeof(NumPE1Reads),"Cnt",&NumPE1Reads);
			gSQLiteSummaries.AddResult(gProcessingID,(char *)"PE2ReadsDeduped",ePTUint32,sizeof(NumPE2Reads),"Cnt",&NumPE2Reads);
			}
		}

	if(NumPE1Reads < cMinSeqs2Assemb)				// arbitary lower limit on number of reads required to continue processing
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: Unable to continue, after deduping insufficent remaining sequences %d, require at least %d",NumPE1Reads,cMinSeqs2Assemb);
		if(gProcessingID > 0)
			gSQLiteSummaries.AddLog(gProcessingID,"Unable to continue, after deduping insufficent remaining sequences %d, require at least %d",NumPE1Reads,cMinSeqs2Assemb);
		Reset(true);
		}
	GenSeqStarts(true);
	GenRdsSfx(1);	// first SeqWrd only requires indexing
	}

#ifdef USEKMERCNTS     // K-mer count filtering is still being developed so currently is not fully implemented

int KMerLenFilter =  MinSeqLen - 10; // arbitary, should be a configurable parameter

if(KMerLenFilter)
	{
	// allocate for KMer sequence instance processing
	m_KMerSeqLen = KMerLenFilter;			
	m_KMerSeqInstSize = (sizeof(tsKMerSeqInst) - 1) +  ((m_KMerSeqLen + 3) / 4);
	m_AllocdKMerSeqInsts = cInitAllocKMerSeqInsts;
	m_AllocdKMerSeqInstsMem = m_AllocdKMerSeqInsts * m_KMerSeqInstSize;
#ifdef _WIN32
	m_pKMerSeqs = (tsKMerSeqInst *)malloc(m_AllocdKMerSeqInstsMem);	// initial and with any luck perhaps the only allocation
	if (m_pKMerSeqs == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "ProcessReadsetDist: Memory allocation of %lld bytes for K-mer distributions - %s", (INT64)m_AllocdKMerSeqInstsMem, strerror(errno));
		m_AllocdKMerSeqInstsMem = 0;
		ARReset();
		Reset();
		return(eBSFerrMem);
		}
#else
		// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
	m_pKMerSeqs = (tsKMerSeqInst *)mmap(NULL, m_AllocdKMerSeqInstsMem, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
	if (m_pKMerSeqs == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "ProcessReadsetDist: Memory allocation of %lld bytes  for K-mer distributions through mmap()  failed - %s", (INT64)m_AllocdKMerSeqInstsMem, strerror(errno));
		m_pKMerSeqs = NULL;
		m_AllocdKMerSeqInstsMem = 0;
		ARReset();
		Reset();
		return(eBSFerrMem);
		}
#endif
	memset(m_pKMerSeqs,0,m_AllocdKMerSeqInstsMem);
	m_NumKMerSeqHashes = 0;
	if((m_pKMerSeqHashes = new UINT32 [cMaxKMerSeqHashArrayEntries]) == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "ProcessReadsetDist: Memory allocation of %d entries for K-mer hashes failed", cMaxKMerSeqHashArrayEntries);
		m_AllocdKMerSeqInstsMem = 0;
		ARReset();
		Reset();
		return(eBSFerrMem);
		}
	}

#else
m_pKMerSeqs = NULL;
m_KMerSeqLen = 0;
m_KMerSeqInstSize = 0;
m_AllocdKMerSeqInsts = 0;
m_AllocdKMerSeqInstsMem = 0;
#endif

	// now identify those reads which are not overlapped on both 5' and 3' by some other read
	// if not overlapped then remove as these are likely to contain sequencer errors
if(MinOverlap != -1)
	RemoveNonOverlaps(MinOverlap,MinFlankLen,IterativePasses);

gSQLiteSummaries.AddResult(gProcessingID,(char *)"Retained",ePTUint32,sizeof(UINT32),"Cnt",&m_Sequences.NumSeqs2Assemb);

FreeSfx();
FreeSeqStarts();

Rslt = SavePackedSeqsToFile(pszOutFile);
ARReset();
Reset(Rslt == eBSFSuccess ? true : false);
return(Rslt);
}

int
CArtefactReduce::RemoveDuplicates(bool bPEdups,				// request that duplicates are for both PE1 and PE2 being duplicates 
									bool bStrand,			// if true then strand specific duplicates
									char *pszDupDist)		// optionally specify a file to which duplicate distributions are to be written
{
int Rslt;

m_bRawDedupe = true;				// removing all duplicate reads

// identify and mark all duplicates
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Identifying duplicate sequences...");
if((Rslt = IdentifyDuplicates(bPEdups,bStrand,pszDupDist))!= eBSFSuccess)	
	return(Rslt);

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Removing duplicate sequences...");
if((Rslt = RemoveMarkedSeqs(cFlgSeqNthDup | cFlgSeqRemove, 0 ,true))!= eBSFSuccess)	
	return(Rslt);

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Removal of duplicate sequences completed");
return(eBSFSuccess);
}

int
CArtefactReduce::RemoveNonOverlaps(int MinOverlap,			// minimum required overlap (in bp)
						int MinFlankLen,            // minimum required non-overlap flank (in bp)
						int NumIterations) 	// because of artefact errors tending to be at end of reads (both 5' and 3') then by default 2 iterations of passes are utilised 
{
int Rslt;
int Iteration;

if(MinOverlap < 25 || MinOverlap > 150 || NumIterations < 1 || NumIterations > 10)
	return(eBSFerrParams);

for(Iteration = 1; Iteration <= NumIterations; Iteration++)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Identification of non-overlapping sequences processing starting, iteration %d ...",Iteration);
	if(Iteration > 1)
		{
		UpdateAllSeqHeaderFlags(0,~(cFlgSeqPE | cFlgSeqPE2),false);
		GenSeqStarts(true);
		GenRdsSfx(1);	// first SeqWrd only requires indexing
		}
	else
		GenSeqStarts(true);

	// start with phase eOvlpSenseToSense
	if((Rslt = IdentifyOverlaps(eOvlpSenseToSense,MinOverlap,MinFlankLen,false))!=0)	
		return(Rslt);

	// followup with  phase eOvlpAntiSenseToSense
	if((Rslt = IdentifyOverlaps(eOvlpAntiSenseToSense,MinOverlap,MinFlankLen,false))!=0)	
		return(Rslt);

	// lastly with phase eOvlpSenseToAntiSense
	PackedRevCplAll(false);
	GenRdsSfx(1);
	if((Rslt = IdentifyOverlaps(eOvlpSenseToAntiSense,MinOverlap,MinFlankLen,true))!=0)	
		return(Rslt);
	PackedRevCplAll();

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Non-overlapping sequences identified, now removing...");
	if(m_Sequences.bPESeqs && !m_bDedupeIndependent)	
		{
		// check that PE's both ends are overlapping; if not then mark both for removal
		UINT16 *pPE1SeqFlags;
		UINT16 *pPE2SeqFlags;
		UINT16 PE1Flags;
		UINT16 PE2Flags;
		tSeqID SeqID;

		pPE1SeqFlags = m_Sequences.pSeqFlags;
		pPE2SeqFlags = pPE1SeqFlags+1;
		for(SeqID = 1; SeqID < m_Sequences.NumSeqs2Assemb; SeqID+=2,pPE1SeqFlags+=2,pPE2SeqFlags+=2)
			{
			// if either not marked as having both 5' and 3' overlaps then both to be deleted
			// if either marked for deletion then both are to be deleted
			PE1Flags = *pPE1SeqFlags;
			PE2Flags = *pPE2SeqFlags;
			if((PE1Flags & (cFlg5Prime | cFlg3Prime)) != (cFlg5Prime | cFlg3Prime))
				PE1Flags |= cFlgSeqRemove;
			if((PE2Flags & (cFlg5Prime | cFlg3Prime)) != (cFlg5Prime | cFlg3Prime))
				PE2Flags |= cFlgSeqRemove;
			if(PE1Flags & (cFlgSeqRemove  | cFlgSeqNthDup))
				PE2Flags |= cFlgSeqRemove;
			if(PE2Flags & (cFlgSeqRemove | cFlgSeqNthDup))
				PE1Flags |= cFlgSeqRemove;
			*pPE1SeqFlags = PE1Flags;
			*pPE2SeqFlags = PE2Flags;
			}
		}

	if((Rslt = RemoveMarkedSeqs(cFlgSeqRemove | cFlgSeqNthDup, cFlg5Prime | cFlg3Prime ,true))!= eBSFSuccess)	
		return(Rslt);
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Removal of non-overlapping sequences completed, iteration %d",Iteration);
	}
return(eBSFSuccess);
}

// IdentifyDuplicates
int
CArtefactReduce::IdentifyDuplicates(bool bPEdups,			// optionally request that duplicates are for both PE1 and PE2 being duplicates
								bool bStrand,			// if true then strand specific duplicates
								char *pszDupDist)		// optionally specify a file to which duplicate distributions are to be written
{
tsThreadIdentDuplicatePars *pThreadParams;
tsThreadIdentDuplicatePars *pCurThread;
UINT32 *pNumDupInstances;
int ThreadIdx;
int NumThreads;
tSeqID CurStartSeqID;
int hDupDistFile;
UINT32 Idx;
UINT32 CurNumProcessed;
UINT32 PrevNumProcessed = 0;
UINT32 CurNumDuplicates = 0;
UINT32 PrevNumDuplicates = 0;
UINT32 MaxDuplicates = 0;
UINT32 MaxDuplicateInsts = 0;

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Starting duplicate %s sequences identification...",bPEdups ? "paired end" : "single ended");

NumThreads = m_NumThreads;

// balance number threads vs the number of sequences so as to minimise the thread startup costs
if(m_Sequences.NumSeqs2Assemb < 10000)
	NumThreads = 1;
else
	{
	NumThreads = (m_Sequences.NumSeqs2Assemb + 9999) / 10000;
	if(NumThreads > m_NumThreads)
		NumThreads = m_NumThreads;
	}

if((pThreadParams = new tsThreadIdentDuplicatePars[NumThreads])==NULL)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Unable to allocate memory for threads...");
	Reset(false);
	return(eBSFerrMem);
	}
memset(pThreadParams,0,sizeof(tsThreadIdentDuplicatePars) * NumThreads);

if((pNumDupInstances = new UINT32[cMaxDupInstances+1])==NULL)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Unable to allocate memory for threads...");
	Reset(false);
	return(eBSFerrMem);
	}
memset(pNumDupInstances,0,sizeof(UINT32) * (cMaxDupInstances+1));

CurStartSeqID = 1;

m_Sequences.NumProcessed = 0;
m_Sequences.NumDuplicates = 0;
m_Sequences.NumOverlapping = 0;
m_Sequences.NumOverlapped = 0;

m_FinalProcSeqID = m_Sequences.NumSeqs2Assemb;
m_NextProcSeqID = 1;
m_StartProcSeqID = 1;
m_NumProcSeqIDs = cMaxMultiSeqFlags;
m_ThreadsProcessing = NumThreads;
pCurThread = pThreadParams;
#ifndef _WIN32
	// increase the default stack of just 2MB
	size_t defaultStackSize;
	struct timespec ts;
	int JoinRslt;
	pthread_attr_t threadattr; 
	pthread_attr_init(&threadattr);
	pthread_attr_getstacksize(&threadattr, &defaultStackSize);
	if(defaultStackSize != cWorkThreadStackSize)
		pthread_attr_setstacksize(&threadattr, cWorkThreadStackSize);
#endif
for(ThreadIdx = 1; ThreadIdx <= NumThreads; ThreadIdx++,pCurThread++)
	{
	pCurThread->ThreadIdx = ThreadIdx;
	pCurThread->bPEdups = bPEdups;
	pCurThread->bDedupeIndependent = m_bDedupeIndependent;
	pCurThread->bStrand = bStrand;
	pCurThread->pThis = this;
	pCurThread->pProbeSubSeq = NULL;
	pCurThread->pPE1SeqWrds = NULL;
	pCurThread->pPE2SeqWrds = NULL;
	if(!bStrand)
		{
		pCurThread->AllocMemPE1SeqWrds = cMaxOvrlapSeqWrds * sizeof(tSeqWrd4);
		pCurThread->pPE1SeqWrds = new UINT8 [pCurThread->AllocMemPE1SeqWrds];
		if(bPEdups)
			{
			pCurThread->AllocMemPE2SeqWrds = cMaxOvrlapSeqWrds * sizeof(tSeqWrd4);
			pCurThread->pPE2SeqWrds = new UINT8 [pCurThread->AllocMemPE2SeqWrds];
			}
		}

#ifdef _WIN32
	pCurThread->threadHandle = (HANDLE)_beginthreadex(NULL,cWorkThreadStackSize,ThreadedIdentDups,pCurThread,0,&pCurThread->threadID);
#else
	pCurThread->threadRslt = pthread_create (&pCurThread->threadID , &threadattr , ThreadedIdentDups , pCurThread);
#endif
	}
#ifndef _WIN32
pthread_attr_destroy(&threadattr);		// no longer required
#endif

// allow threads a few seconds to startup
#ifdef _WIN32
	Sleep(3000);
#else
	sleep(3);
#endif
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Progress: 0 sequences processed");
CurNumProcessed = 0;
PrevNumProcessed = 0;
CurNumDuplicates = 0;
PrevNumDuplicates = 0;
MaxDuplicates = 0;
pCurThread = pThreadParams;
for(ThreadIdx = 0; ThreadIdx < NumThreads; ThreadIdx++,pCurThread++)
	{
#ifdef _WIN32
	while(WAIT_TIMEOUT == WaitForSingleObject(pCurThread->threadHandle, 60000 * 10))
		{
		AcquireLock(false);
		CurNumProcessed = m_Sequences.NumProcessed;
		CurNumDuplicates = m_Sequences.NumDuplicates;
		ReleaseLock(false);
		if(CurNumProcessed > PrevNumProcessed)
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"Progress: %u sequences processed",CurNumProcessed);
		PrevNumProcessed = CurNumProcessed;
		PrevNumDuplicates = CurNumDuplicates;
		}
	CloseHandle( pCurThread->threadHandle);
#else
	struct timespec ts;
	int JoinRlt;
	clock_gettime(CLOCK_REALTIME, &ts);
	ts.tv_sec += 60 * 10;
	while((JoinRlt = pthread_timedjoin_np(pCurThread->threadID, NULL, &ts)) != 0)
		{
		AcquireLock(false);
		CurNumProcessed = m_Sequences.NumProcessed;
		CurNumDuplicates = m_Sequences.NumDuplicates;
		ReleaseLock(false);
		if(CurNumProcessed > PrevNumProcessed)
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"Progress: %u sequences processed",CurNumProcessed);
		PrevNumProcessed = CurNumProcessed;
		PrevNumDuplicates = CurNumDuplicates;
		ts.tv_sec += 60;
		}
#endif
	if(pCurThread->pProbeSubSeq != NULL)
		delete (UINT8 *)pCurThread->pProbeSubSeq;
	if(pCurThread->pPE1SeqWrds != NULL)
		delete (UINT8 *)pCurThread->pPE1SeqWrds;
	if(pCurThread->pPE2SeqWrds != NULL)
		delete (UINT8 *)pCurThread->pPE2SeqWrds;
	if(pCurThread->MaxDuplicates > MaxDuplicateInsts)
		MaxDuplicateInsts = pCurThread->MaxDuplicates;
	for(Idx = 0; Idx <= cMaxDupInstances; Idx++)
		{
		if(pCurThread->NumDupInstances[Idx] && Idx > MaxDuplicates) 
			MaxDuplicates = Idx;

		pNumDupInstances[Idx] += pCurThread->NumDupInstances[Idx];
		}
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Completed: %u sequences processed, provisionally %d are duplicates",m_Sequences.NumProcessed,m_Sequences.NumDuplicates);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Maximum number of provisional duplicates for any sequence was %u",MaxDuplicateInsts);

if(gProcessingID > 0)
	{
	gSQLiteSummaries.AddResult(gProcessingID,(char *)"ReadsDuplicated",ePTUint32,sizeof(m_Sequences.NumDuplicates),"Cnt",&m_Sequences.NumDuplicates);
	gSQLiteSummaries.AddResult(gProcessingID,(char *)"ReadsMaxDuplicate",ePTUint32,sizeof(MaxDuplicateInsts),"Cnt",&MaxDuplicateInsts);
	}

if(pszDupDist != NULL && pszDupDist[0] != '\0')
	{
#ifdef _WIN32
	if((hDupDistFile = open(pszDupDist, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
	if((hDupDistFile = open(pszDupDist, O_RDWR | O_CREAT |O_TRUNC, S_IREAD | S_IWRITE))==-1)
#endif
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create or truncate %s - %s",pszDupDist,strerror(errno));
		if(pThreadParams)
			delete pThreadParams;
		if(pNumDupInstances)
			delete pNumDupInstances;
		Reset();
		return(eBSFerrCreateFile);
		}

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Output duplicate sequence distribution file created/truncated: '%s'",pszDupDist);
	}
else
	hDupDistFile = -1;


if(hDupDistFile != -1)
	{
	char szBuff[4096];
	int BuffIdx = 0;
	BuffIdx += sprintf(&szBuff[BuffIdx],"\"Copies\",\"NumInstances\"\n");
	for(Idx = 0; Idx <= MaxDuplicates; Idx++)
		{
		if(pNumDupInstances[Idx] == 0)
			continue;
		BuffIdx += sprintf(&szBuff[BuffIdx],"%d,%d\n",Idx+1,pNumDupInstances[Idx]);
		if(BuffIdx + 100 > sizeof(szBuff))
			{
			BuffIdx = write(hDupDistFile,szBuff,BuffIdx);
			BuffIdx = 0;
			}
		}
	if(BuffIdx > 0)
		BuffIdx = write(hDupDistFile,szBuff,BuffIdx);
#ifdef _WIN32
	_commit(hDupDistFile);
#else
	fsync(hDupDistFile);
#endif
	close(hDupDistFile);
	hDupDistFile = -1;
	}


if(gProcessingID > 0)
	{
	for(Idx = 1; (Idx <= MaxDuplicates+1); Idx++)
		{
		if(pNumDupInstances[Idx] == 0)
			continue;
		gSQLiteSummaries.AddResultXY(gProcessingID,(char *)"DupReadsDist",ePTUint32,sizeof(UINT32),"Copies",&Idx,ePTUint32,sizeof(UINT32),"NumInstances",&pNumDupInstances[Idx-1]);
		}
	}

if(pThreadParams)
	delete pThreadParams;
if(pNumDupInstances)
	delete pNumDupInstances;

if(!m_bDedupeIndependent && bPEdups)
	{
	// iterate PE1 and PE2 flags and if either has cFlgSeqNthDup then set the other end to also have cFlgSeqNthDup
	UINT16 *pPE1SeqFlags;
	UINT16 *pPE2SeqFlags;

	tSeqID SeqID;
	pPE1SeqFlags = m_Sequences.pSeqFlags;
	pPE2SeqFlags = pPE1SeqFlags + 1;
	for(SeqID = 1; SeqID < m_Sequences.NumSeqs2Assemb; SeqID+=2,pPE1SeqFlags+=2,pPE2SeqFlags+=2)
		{
		if(*pPE1SeqFlags & cFlgSeqNthDup)
			*pPE2SeqFlags |= cFlgSeqNthDup;
		else
			if(*pPE2SeqFlags & cFlgSeqNthDup)
				*pPE1SeqFlags |= cFlgSeqNthDup;
		}
	}

return(eBSFSuccess);
}


// IdentifyOverlaps
int
CArtefactReduce::IdentifyOverlaps(etOvlFlankPhase OvlFlankPhase,	// overlap flank processing phase
							int MinOverlap,			// sequences must overlap by at least this number of bases
							int MinFlankLen,            // minimum required non-overlap flank (in bp)
							bool bRevCpl)			// if true then all sequences have been reverse complemented and sfx index is over these sequences
							
{
tsThreadIdentOverlapPars *pThreadParams;
tsThreadIdentOverlapPars *pCurThread;
int ThreadIdx;
int NumThreads;
char *pszPhaseDescr;
UINT32 ProvNum2Retain;
UINT32 Idx;
UINT16 *pSeqFlags;
tSeqID CurStartSeqID;
UINT32 CurNumProcessed;
UINT32 PrevNumProcessed = 0;

switch(OvlFlankPhase) {
	case eOvlpSenseToSense:					// it's sense overlap sense flank processing (probe cFlg5Prime and target cFlg3Prime set if overlap) 
		pszPhaseDescr = (char *)"Sense overlap Sense (3' - 5')";
		break;
	case eOvlpAntiSenseToSense:				// it's antisense overlap sense flank processing (probe cFlg5Prime and target cFlg5Prime set if overlap)
		pszPhaseDescr = (char *)"Antisense overlap Sense (5' - 5')";
		break;
	case eOvlpSenseToAntiSense:				// it's sense overlap antisense flank processing (probe cFlg3Prime and target cFlg3Prime set if overlap)
		pszPhaseDescr = (char *)"Sense overlap Antisense (3' - 3')";
		break;
	};

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Starting flank %s sequence overlap identification",pszPhaseDescr);

// ensure that cFlgNoProc has been globally reset
pSeqFlags = m_Sequences.pSeqFlags;
for(Idx = 0; Idx < m_Sequences.NumSeqs2Assemb; Idx++,pSeqFlags++)
	*pSeqFlags &= ~cFlgNoProc;

NumThreads = m_NumThreads;

// balance the number threads vs the number of sequences to minimise the thread startup costs
if(m_Sequences.NumSeqs2Assemb < 10000)
	NumThreads = 1;
else
	{
	NumThreads = (m_Sequences.NumSeqs2Assemb + 9999) / 10000;
	if(NumThreads > m_NumThreads)
		NumThreads = m_NumThreads;
	}

if((pThreadParams = new tsThreadIdentOverlapPars[NumThreads])==NULL)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Unable to allocate memory for threads...");
	Reset(false);
	return(eBSFerrMem);
	}
memset(pThreadParams,0,sizeof(tsThreadIdentOverlapPars) * NumThreads);

CurStartSeqID = 1;
m_Sequences.NumProcessed = 0;
m_Sequences.NumDuplicates = 0;
m_Sequences.NumOverlapping = 0;
m_Sequences.NumOverlapped = 0;
m_FinalProcSeqID = m_Sequences.NumSeqs2Assemb;
m_NextProcSeqID = 1;
m_StartProcSeqID = 1;
m_NumProcSeqIDs = cMaxMultiSeqFlags;
pCurThread = pThreadParams;

#ifndef _WIN32
	// increase the default stack of just 2MB
	size_t defaultStackSize;
	struct timespec ts;
	int JoinRslt;
	pthread_attr_t threadattr; 
	pthread_attr_init(&threadattr);
	pthread_attr_getstacksize(&threadattr, &defaultStackSize);
	if(defaultStackSize < cWorkThreadStackSize)
		pthread_attr_setstacksize(&threadattr, cWorkThreadStackSize);
#endif
for(ThreadIdx = 1; ThreadIdx <= NumThreads; ThreadIdx++,pCurThread++)
	{
	pCurThread->ThreadIdx = ThreadIdx;
	pCurThread->pThis = this;
	pCurThread->bRevCpl = bRevCpl;
	pCurThread->MinOverlap = MinOverlap;
	pCurThread->MinFlankLen = MinFlankLen;
	pCurThread->OvlFlankPhase = OvlFlankPhase;
	pCurThread->AllocMemOverlapSeq = cMaxOvrlapSeqWrds * sizeof(tSeqWrd4);	
	pCurThread->pOverlapSeq = new UINT8 [pCurThread->AllocMemOverlapSeq];
	pCurThread->pOverlapFlankSeq = new UINT8 [pCurThread->AllocMemOverlapSeq];
#ifdef _WIN32
	pCurThread->threadHandle = (HANDLE)_beginthreadex(NULL,cWorkThreadStackSize,ThreadedIdentOverlaps,pCurThread,0,&pCurThread->threadID);
#else
	pCurThread->threadRslt = pthread_create (&pCurThread->threadID , &threadattr , ThreadedIdentOverlaps , pCurThread);
#endif
	}
#ifndef _WIN32
pthread_attr_destroy(&threadattr);		// no longer required
#endif

// allow threads a few seconds to startup
#ifdef _WIN32
	Sleep(5000);
#else
	sleep(5);
#endif

CurNumProcessed = 0;
PrevNumProcessed = 0;
pCurThread = pThreadParams;
for(ThreadIdx = 0; ThreadIdx < NumThreads; ThreadIdx++,pCurThread++)
	{
#ifdef _WIN32
	while(WAIT_TIMEOUT == WaitForSingleObject(pCurThread->threadHandle, 60000))
		{
		AcquireLock(false);
		CurNumProcessed = m_Sequences.NumProcessed;
		ReleaseLock(false);
		if(CurNumProcessed > PrevNumProcessed)
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"Progress: %u sequences processed",CurNumProcessed);
		PrevNumProcessed = CurNumProcessed;
		}
	CloseHandle( pCurThread->threadHandle);
#else
	struct timespec ts;
	int JoinRlt;
	clock_gettime(CLOCK_REALTIME, &ts);
	ts.tv_sec += 60 * 10;
	while((JoinRlt = pthread_timedjoin_np(pCurThread->threadID, NULL, &ts)) != 0)
		{
		AcquireLock(false);
		CurNumProcessed = m_Sequences.NumProcessed;
		ReleaseLock(false);
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Progress: %u sequences processed",CurNumProcessed);
		PrevNumProcessed = CurNumProcessed;
		ts.tv_sec += 60;
		}
#endif
	delete (UINT8 *)pCurThread->pOverlapSeq;
	delete (UINT8 *)pCurThread->pOverlapFlankSeq;
	}

// report provisional number of reads which will be retained
ProvNum2Retain = 0;
pSeqFlags = m_Sequences.pSeqFlags;
for(Idx = 0; Idx < m_Sequences.NumSeqs2Assemb; Idx++,pSeqFlags++)
	{
	*pSeqFlags &= ~cFlgNoProc;
	if(*pSeqFlags & cFlg5Prime && *pSeqFlags & cFlg3Prime)
		ProvNum2Retain += 1;
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Completed flank %s sequence overlap identification, %u sequences processed, provisionally %u will be retained",pszPhaseDescr,m_Sequences.NumProcessed,ProvNum2Retain);

if(pThreadParams)
	delete pThreadParams;
return(eBSFSuccess);
}

// ProcIdentDuplicates
// Identify and mark identical sequence duplicates
// If paired end processing then both the 5' and 3' ends normally be identical to another paired end's 5' and 3' to be marked as being identical
// although there is an option for paired ends to be treated as if single ended, i.e. the 5' and 3' ends are independent of each other but still
// require both ends to be marked as duplicates.
int
CArtefactReduce::ProcIdentDuplicates(tsThreadIdentDuplicatePars *pPars)
{
int CmpRslt;
tSeqID StartingSeqID;
tSeqID EndingSeqID;
tsMultiSeqFlags MultiSeqFlags[cMaxMultiSeqFlags+1];			// always allow for 1 extra!
int TooMAnyWarnings;
tSeqID SeqID;
tSeqID MatchTargID;
UINT32 NumProcessed;
UINT32 NumUniques;
UINT32 NumDuplicates;
UINT32 NumRevCplDups;
UINT32 CurDuplicates;
UINT32 MaxDuplicates;

UINT64 SfxWrdIdx;
UINT64 SfxWrdLowHit;
UINT32 ProbeLen;
UINT32 ProbeFlags;
UINT32 TargFlags;
UINT32 TargLen;
tSeqWrd4 *pStartSeqWrd;
tSeqWrd4 *pTarg;

tSeqWrd4 *pPE2StartSeqWrd;
UINT32 PE2ProbeFlags;
UINT32 PE2ProbeLen;
UINT32 PE2TargLen;
UINT32 PE2TargFlags;
tSeqWrd4 *pPE2Targ;

int DiagLevel = gDiagnostics.GetFileDiagLevel();

gDiagnostics.DiagOut(eDLDebug,gszProcName,"Thread %d startup for duplicate %s identification...",pPars->ThreadIdx, 
										pPars->bPEdups ? "paired end sequences" : "sequences");


NumProcessed = 0;
NumUniques = 0;
NumDuplicates = 0;
CurDuplicates = 0;
NumRevCplDups = 0;
MaxDuplicates = 0;
TooMAnyWarnings = 0;

time_t Started = time(0);

while(GetSeqProcRange(&StartingSeqID,&EndingSeqID,cMaxMultiSeqFlags) > 0)
	{
	// get thread local copy of flags for range of sequences to be processed
	GetSeqFlags(StartingSeqID,EndingSeqID,MultiSeqFlags,true);
	for(SeqID = StartingSeqID; SeqID <= EndingSeqID; SeqID++)
		{
		NumProcessed+=1;
		if(!(NumProcessed % 2000))
			{
			time_t Now = time(0);
			unsigned long ElapsedSecs = (unsigned long) (Now - Started);
			if(ElapsedSecs >= 30)
				{
				pPars->NumProcessed += NumProcessed;
				pPars->NumDuplicates += NumDuplicates;
				AcquireLock(true);
				m_Sequences.NumProcessed += NumProcessed;
				m_Sequences.NumDuplicates += NumDuplicates;
				ReleaseLock(true);
				NumDuplicates = 0;
				NumProcessed = 0;
				Started = Now;
				}
			}

		if(!pPars->bDedupeIndependent && (MultiSeqFlags[SeqID - StartingSeqID].CurFlags & cFlgSeqPE2)) 
			continue;

		if((ProbeFlags = MultiSeqFlags[SeqID - StartingSeqID].CurFlags) & cFlgSeqNthDup)		// if already marked as a duplicate then skip
			{
			if(!pPars->bDedupeIndependent && pPars->bPEdups)									// if paired end duplicate detection then also mark 3' end as a duplicate
				MultiSeqFlags[1 + SeqID - StartingSeqID].SetFlags |= cFlgSeqNthDup;          
//			NumDuplicates += 1;
			continue;
			}
		if(!pPars->bDedupeIndependent && pPars->bPEdups)										// if paired end duplicate detection and 3' end marked then mark the 5' end as a duplicate
			{
			if((PE2ProbeFlags = MultiSeqFlags[1 + SeqID - StartingSeqID].CurFlags) & cFlgSeqNthDup)		
				{
				MultiSeqFlags[SeqID - StartingSeqID].SetFlags |= cFlgSeqNthDup;
				continue;
				}
			}

		pStartSeqWrd = GetSeqHeader(SeqID,NULL,NULL,&ProbeLen,false);
		if(pStartSeqWrd == NULL)				// serious problem if can't get header, currently just writes to log and continues unless there are simply too many...
			{
			if((TooMAnyWarnings+=1) < 10)
				gDiagnostics.DiagOut(eDLWarn,gszProcName,"Thread %d Couldn't find sequence header for known sequence %d...",pPars->ThreadIdx,SeqID);
			continue;
			}
		
		if(!pPars->bDedupeIndependent && pPars->bPEdups)					// if paired end processing then also need header for 3' end
			{
			if((pPE2StartSeqWrd = GetSeqHeader(SeqID+1,NULL,NULL,&PE2ProbeLen,false))==NULL)
				{
				if((TooMAnyWarnings+=1) < 10)
					gDiagnostics.DiagOut(eDLWarn,gszProcName,"Thread %d Couldn't find PE2 sequence header for SeqID %u...",pPars->ThreadIdx,SeqID+1);
				continue;
				}
			}

		int StrandPhase;
		bool bRevCpl;
		bool bSeenSelf;
		if(!pPars->bStrand)	// if strand independent then need to try for duplicates on both strands...
			StrandPhase = 2;
		else
			StrandPhase = 1;

		bRevCpl = false;
		bSeenSelf = false;
		while(StrandPhase > 0)		// 1: 1st if strand dependent, 2nd pass if stand independent. 2: 1st if strand independent
			{
			SfxWrdIdx =	LocateFirstExact(m_Sequences.SfxElSize,		// sizeof elements in pSfxArray - currently will be either 4 or 5 bytes
							  pStartSeqWrd,							// pts to probe sequence
							  ProbeLen,								// probe length (in bases, not tSeqWrd4's) to exactly match over
							  m_Sequences.pSeqs2Assemb,				// target sequences
							  (UINT8 *)m_Sequences.pSuffixArray,	// target sequence suffix array
							  0,									// low index in pSfxArray
							  m_Sequences.NumSuffixEls-1);			// high index in pSfxArray

			SfxWrdLowHit = SfxWrdIdx;
			if(SfxWrdIdx == 0)
				{
				if(bRevCpl)		// after reverse complementing then it's to be expected that there may be no exact matches
					{
					StrandPhase -= 1;
					continue;
					}
				
				if((TooMAnyWarnings+=1) < 10)	// if not reverse complementing then surely must have at least one exact match on self
					{
					gDiagnostics.DiagOut(eDLWarn,gszProcName,"Thread %d (W%d) Couldn't locate any instance of matching sequence for SeqID %u...",pPars->ThreadIdx,TooMAnyWarnings,SeqID);
					gDiagnostics.DiagOut(eDLWarn,gszProcName,"SeqID(%d) %s",SeqID,AsciifySequence(SeqID,200),false);
					}
				break;
				}

			CmpRslt = 0;
			CurDuplicates = 0;
			SfxWrdIdx -= 1;				// LocateFirstExact returned located index + 1
			do  {
				SfxWrdIdx += 1;
				if((pTarg = SfxIdxToFirstSeqWrd(SfxWrdIdx))==NULL)
					break;
				pTarg = GetSeqHeader(pTarg,&MatchTargID,NULL,NULL,&TargLen,false);
				if(MatchTargID <= SeqID)
					{
					if(!bRevCpl && MatchTargID == SeqID)		// expect to have seen self!
						bSeenSelf = true;
					continue;
					}

				if(!pPars->bDedupeIndependent && pPars->bPEdups && !(MatchTargID & 0x00000001))   // if paired ends then sequences with even numbered SeqIDs are the 3' ends so skip these
					continue;

				if((CmpRslt = CmpPackedSeqs(pStartSeqWrd,pTarg,ProbeLen))==0)
					{
					// As sequences are processed in sequence identifier order then can safely slough sequences with lower order sequence identifiers
					if(MatchTargID <= SeqID || TargLen != ProbeLen)
						continue;
					TargFlags = m_Sequences.pSeqFlags[MatchTargID-1];	// not too concerned with serialising flag reads as flags of interest are contained within a single byte so atomically updated 
					// if already processed by some other thread then look for another matching sequence
					if(TargFlags & (cFlgSeqNthDup | cFlgSeq1stDup | cFlgSeqUnique)) 
						continue;

					// if paired end processing then check the 3' PE2 to see if that is also a duplicate
					// only mark paired ends as duplicates if both are identical to the probe PE1
					if(!pPars->bDedupeIndependent && pPars->bPEdups)
						{
						pPE2Targ = GetSeqHeader(MatchTargID+1,NULL,NULL,&PE2TargLen,false);  // 3' end is immediately after the 5' end so simply add 1 on to the 5' seqid
						if(PE2TargLen != PE2ProbeLen)
							continue;
						if((CmpRslt = CmpPackedSeqs(pPE2StartSeqWrd,pPE2Targ,PE2ProbeLen))!=0)
							{
							CmpRslt = 0;	// try next PE1
							continue;
							}

						// both PE1 and PE2 are exact matches so mark both as being duplicates to be removed
						PE2TargFlags = m_Sequences.pSeqFlags[MatchTargID+1]; // not too concerned with serialising flag reads as flags of interest are contained within a single byte
						if(!(PE2TargFlags & cFlgSeqNthDup))				// serialisation is required so only update flags if not already marked as a duplicate
							{
							NumDuplicates += 1;
							UpdateSeqFlags(MatchTargID + 1,cFlgSeqNthDup,0,true);		
							}
						}
			
					// mark PE1 as being a duplicate
					if(!(TargFlags & cFlgSeqNthDup))						// serialisation is required so only update flags if not already marked as a duplicate
						{
						NumDuplicates += 1;
						UpdateSeqFlags(MatchTargID,cFlgSeqNthDup,0,true);	
						}
					if(bRevCpl)
						NumRevCplDups += 1;
					CurDuplicates += 1;
					}
				}
			while(!CmpRslt);

			if(CurDuplicates > MaxDuplicates)				// interested in stats of max number duplicates for any sequence instance
				MaxDuplicates = CurDuplicates;
			if(CurDuplicates >  cMaxDupInstances)
				CurDuplicates = cMaxDupInstances;
			pPars->NumDupInstances[CurDuplicates] += 1;
			int ProcFlags;
			if(CurDuplicates == 0)
				{
				NumUniques += 1;
				ProcFlags = cFlgSeqUnique;
				}
			else
				ProcFlags = cFlgSeq1stDup;

			MultiSeqFlags[SeqID - StartingSeqID].SetFlags = ProcFlags;
			if(!pPars->bDedupeIndependent && pPars->bPEdups)
				MultiSeqFlags[1 + SeqID - StartingSeqID].SetFlags = ProcFlags;

			if(StrandPhase == 2)		// if just finished phase 2 (must have been strand insensitive duplicate processign) then that means phase 1 is with reverse complemented sequences
				{
				bRevCpl = true;
				int Len;
				// make copy and reverse complement
				GetPackedSeq(0,pStartSeqWrd,(tSeqWrd4 *)pPars->pPE1SeqWrds);
				pStartSeqWrd = (tSeqWrd4 *)pPars->pPE1SeqWrds;
				Len=PackedRevCpl(pStartSeqWrd);
				if(!pPars->bDedupeIndependent && pPars->bPEdups)
					{
					GetPackedSeq(0,pPE2StartSeqWrd,(tSeqWrd4 *)pPars->pPE2SeqWrds);
					pPE2StartSeqWrd = (tSeqWrd4 *)pPars->pPE2SeqWrds;
					Len=PackedRevCpl(pPE2StartSeqWrd);
					pTarg = pStartSeqWrd;
					pStartSeqWrd = pPE2StartSeqWrd;
					pPE2StartSeqWrd = pTarg;
					ProbeLen = Len;
					}
				}
			StrandPhase -= 1;
			}
		if(!bSeenSelf)
			{
			if((TooMAnyWarnings+=1) < 10)
				{
				gDiagnostics.DiagOut(eDLWarn,gszProcName,"Thread %d Couldn't locate self instance of matching sequence for SeqID %u...",pPars->ThreadIdx,SeqID);
				gDiagnostics.DiagOut(eDLWarn,gszProcName,"SeqID(%d) %s",SeqID,AsciifySequence(SeqID,200),false);
				}
			}
		}
	UpdateSeqFlags(1 + EndingSeqID - StartingSeqID,MultiSeqFlags,true);			// serialisation required
	}
pPars->MaxDuplicates = MaxDuplicates;
pPars->NumProcessed += NumProcessed;
pPars->NumDuplicates += NumDuplicates;
AcquireLock(true);
m_Sequences.NumProcessed += NumProcessed;
m_Sequences.NumDuplicates += NumDuplicates;
ReleaseLock(true);

gDiagnostics.DiagOut(eDLDebug,gszProcName,"Thread %d completed duplicate %s identification",pPars->ThreadIdx, pPars->bPEdups ? "paired end sequences" : "sequences");
gDiagnostics.DiagOut(eDLDebug,gszProcName,"Thread %d discovered NumRevCplDups %d",pPars->ThreadIdx,NumRevCplDups);

return(1);		// success
}

// ProcIdentOverlaps
// Iterates all sequences which are not marked as duplicates and identifies those which overlap at least one other sequence
// Phase eOvlpSenseToSense:
// Assumes cFlgNoProc is reset on all sequences
// Probe:			aaaacgtagatgatagat					Probe will be marked cFlg3Prime + cFlgNoProc 
// Target:				 gtagatgatagataaaatttttt		Target will be marked cFlg5Prime
// If a given Probe is is identified as overlapping any other sequence then, as an optimisation, all other sequences which are identical to the given Probe
// are given the same overlap flag and thus when iterating Probe sequences if marked with cFlgNoProc then these sequences need not be checked
// for any otherlaps in this phase

//
// Phase eOvlpAntiSenseToSense:
// Assumes cFlgNoProc is reset on all sequences
// When Probe is antisense to target then this is detected by looking for probe antisense overlaps onto reverse complemented targets
// RevCpl Probe:    aaaacgtagatgatagat				    will be marked  cFlg5Prime + cFlgNoProc
// Target:                  gtagatgatagataaaatttttt		will be marked  cFlg5Prime
// If a given Probe is is identified as overlapping any other sequence then, as an optimisation, all other sequences which are identical to the given Probe
// are given the same overlap flags and thus when iterating Probe sequences if marked with cFlgNoProc then these sequences need not be checked
// for any otherlaps in this phase
//
//
// eOvlpSenseToAntiSense:
// Assumes cFlgNoProc is reset on all sequences
// When Probe overlays target coming from the antisense strand then this is detected by looking for probe sense overlaps onto reverse complemented targets
// Probe:      		aaaacgtagatgatagat					Probe will be marked cFlg3Prime + cFlgNoProc
// RevCpl Target:        gtagatgatagataaaatttttt		Target will be marked cFlg3Prime
// If a given Probe is is identified as overlapping any other sequence then, as an optimisation, all other sequences which are identical to the given Probe
// are given the same overlap flags and thus when iterating Probe sequences if marked with cFlgNoProc then that sequences need not be checked
// for any otherlaps in this phase
// 
// When all overlaps have been identified then any those sequence without both cFlg5Prime and cFlg3Prime will be filtered out from subsequent processing
//
int
CArtefactReduce::ProcIdentOverlaps(tsThreadIdentOverlapPars *pPars)
{
int CmpRslt;
int SubOfs;
int MinOverlap;
int MinFlankLen;
int ProcFlags;
int FlgOverlapping;
int FlgOverlapped;
tSeqID StartingSeqID;
tSeqID EndingSeqID;
int TooMAnyWarnings;
tSeqID SeqID;
tSeqID MatchTargID;
UINT32 NumProcessed;

UINT32 NumOverlapping;
UINT32 NumOverlapped;

tsMultiSeqFlags MultiSeqFlags[cMaxMultiSeqFlags+1];		// allow for 1 extra!

UINT64 SfxWrdIdx;
UINT32 ProbeLen;
UINT32 PrevProbeLen;
UINT32 ProbeFlags;
UINT32 TargFlags;
UINT32 TargLen;
tSeqWrd4 *pStartSeqWrd;

tSeqWrd4 *pTarg;

gDiagnostics.DiagOut(eDLDebug,gszProcName,"Thread %d startup for overlap identification...",pPars->ThreadIdx);

NumProcessed = 0;
TooMAnyWarnings = 0;
NumOverlapping = 0;
NumOverlapped = 0;
switch(pPars->OvlFlankPhase) {
	case eOvlpSenseToSense:					// it's sense overlap sense flank processing (probe cFlg3Prime and target cFlg5Prime set if overlap) 
		FlgOverlapping = cFlg3Prime;
		FlgOverlapped = cFlg5Prime;
		break;
	case eOvlpAntiSenseToSense:				// it's antisense overlap sense flank processing (probe cFlg5Prime and target cFlg5Prime set if overlap)
		FlgOverlapping = cFlg5Prime;
		FlgOverlapped = cFlg5Prime;
		break;
	case eOvlpSenseToAntiSense:				// it's sense overlap antisense flank processing (probe cFlg3Prime and target cFlg3Prime set if overlap)
		FlgOverlapping = cFlg3Prime;
		FlgOverlapped = cFlg3Prime;
		break;
	}

time_t Started = time(0);
while(GetSeqProcRange(&StartingSeqID,&EndingSeqID,cMaxMultiSeqFlags) > 0)
	{
	PrevProbeLen = 0;

		// get thread local copy of flags for range of sequences to be processed
	GetSeqFlags(StartingSeqID,EndingSeqID,MultiSeqFlags,true);
	for(SeqID = StartingSeqID; SeqID <= EndingSeqID; SeqID++)
		{
		NumProcessed+=1;
		if(!(NumProcessed % 1000))
			{
			time_t Now = time(0);
			unsigned long ElapsedSecs = (unsigned long) (Now - Started);
			if(ElapsedSecs >= 30)
				{
				pPars->NumProcessed += NumProcessed;
				pPars->NumOverlapping += NumOverlapping;
				pPars->NumOverlapped += NumOverlapped;
				AcquireLock(true);
				m_Sequences.NumProcessed += NumProcessed;
				m_Sequences.NumOverlapping += NumOverlapping;
				m_Sequences.NumOverlapped += NumOverlapped;
				ReleaseLock(true);
				NumProcessed = 0;
				NumOverlapping = 0;
				NumOverlapped = 0;
				Started = Now;
				}
			}

		// no further processing required on this sequence if cFlgNoProc set
		if((ProbeFlags = MultiSeqFlags[SeqID - StartingSeqID].CurFlags) & cFlgNoProc)		
			continue;
		
		pStartSeqWrd = GetSeqHeader(SeqID,NULL,NULL,&ProbeLen,false);
		if(pStartSeqWrd == NULL)				// serious problem if can't get header...
			{
			if((TooMAnyWarnings+=1) < 10)
				gDiagnostics.DiagOut(eDLWarn,gszProcName,"Thread %d Couldn't find sequence header for known sequence %d...",pPars->ThreadIdx,SeqID);
			continue;
			}
	
		if(pPars->MinOverlap == 0)			// if not specified then default to be 80% of the probe length
			MinOverlap = (79+(ProbeLen * 80)) / 100; // but normally expect the default to have been set to be 80% of read length mean
		else
			MinOverlap = pPars->MinOverlap;
		if(MinOverlap < 25)
			MinOverlap = 25;
		else
			if(MinOverlap > 150)
				MinOverlap = 150;

		MinFlankLen = pPars->MinFlankLen;
		if(MinFlankLen < 1)
			MinFlankLen = 1;
		else
			if(MinFlankLen > 25)
				MinFlankLen = 25;

		if((MinOverlap + MinFlankLen) > (int)ProbeLen)			// no point in processing this sequence further if too short to overlap any other sequence with at least 1 base overhang
			continue;

		// make a copy of probe sequence as will be slicing and dicing when subsequencing...
		GetSeqWrdSubSeq(0,ProbeLen, pStartSeqWrd, (tSeqWrd4 *)pPars->pOverlapSeq); 

		// need to revcpl the probe sequence if not processing sense overlapping sense
		if(pPars->OvlFlankPhase != eOvlpSenseToSense)
			PackedRevCpl((tSeqWrd4 *)pPars->pOverlapSeq);

		// determine if this read overlaps any other reads by at least MinOverlap 
		CmpRslt = 0;
		ProcFlags = 0;
		for(SubOfs = MinFlankLen; SubOfs <= ((int)ProbeLen - MinOverlap); SubOfs++)
			{
			GetSeqWrdSubSeq(SubOfs,ProbeLen - SubOfs, (tSeqWrd4 *)pPars->pOverlapSeq, (tSeqWrd4 *)pPars->pOverlapFlankSeq); 
	
			SfxWrdIdx =	LocateFirstExact(m_Sequences.SfxElSize,		// sizeof elements in pSfxArray - currently will be either 4 or 5 bytes
						pPars->pOverlapFlankSeq,			// pts to probes flank subsequence
						ProbeLen - SubOfs,					// probe length (in bases, not tSeqWrd4's) to exactly match over
						m_Sequences.pSeqs2Assemb,			// target sequence
						(UINT8 *)m_Sequences.pSuffixArray,	// target sequence suffix array
						0,									// low index in pSfxArray
						m_Sequences.NumSuffixEls-1);			// high index in pSfxArray
	
			CmpRslt = 0;
			do  {
				if(!SfxWrdIdx || SfxWrdIdx > m_Sequences.NumSuffixEls)
					break;
				if((pTarg = SfxIdxToFirstSeqWrd(SfxWrdIdx++))==NULL)
					break;

				MatchTargID = 0;
				pTarg = GetSeqHeader(pTarg,&MatchTargID,NULL,NULL,&TargLen,false);

				if(pTarg == NULL || MatchTargID == 0 || MatchTargID > m_Sequences.NumSeqs2Assemb)
					{
					gDiagnostics.DiagOut(eDLWarn,gszProcName,"Thread %d Couldn't find header word for %u ...",pPars->ThreadIdx,MatchTargID);
					break;
					}

				if((CmpRslt = CmpPackedSeqs((tSeqWrd4 *)pPars->pOverlapFlankSeq,pTarg,ProbeLen - SubOfs))!=0)
					break;
				
				if(MatchTargID == SeqID)			// if self then try any next target sequence
					continue;
				// if already known probe is overlapping and target also overlapped then no need to check target again if overlapped
				TargFlags = m_Sequences.pSeqFlags[MatchTargID-1];
				if(ProcFlags & FlgOverlapping && TargFlags & FlgOverlapped)		
					continue;

				if(!(ProcFlags & FlgOverlapping))
					{
					ProcFlags |= FlgOverlapping;
					NumOverlapping += 1;
					}

				if(!(TargFlags & FlgOverlapped))
					{
					UpdateSeqFlags(MatchTargID,FlgOverlapped,0,true);
					NumOverlapped += 1;
					}
				}
			while(!CmpRslt);
			}

		if(ProcFlags)		// if any flags set for this sequence then can propagate these to all other identical sequences
			{
			MultiSeqFlags[SeqID - StartingSeqID].SetFlags = ProcFlags | cFlgNoProc;

			SfxWrdIdx =	LocateFirstExact(m_Sequences.SfxElSize,	// sizeof elements in pSfxArray - currently will be either 4 or 5 bytes
						  pStartSeqWrd,							// pts to probe sequence
						  ProbeLen,								// probe length (in bases, not tSeqWrd4's) to exactly match over
						  m_Sequences.pSeqs2Assemb,				// target sequence
						  (UINT8 *)m_Sequences.pSuffixArray,	// target sequence suffix array
						  0,									// low index in pSfxArray
						  m_Sequences.NumSuffixEls-1);			// high index in pSfxArray
			CmpRslt = 0;
			do {
				if(!SfxWrdIdx || SfxWrdIdx > m_Sequences.NumSuffixEls)
					break;
				if((pTarg = SfxIdxToFirstSeqWrd(SfxWrdIdx++)) == NULL)
					break;

				MatchTargID = 0;
				pTarg = GetSeqHeader(pTarg,&MatchTargID,NULL,NULL,&TargLen,false);
				if(pTarg == NULL || MatchTargID == 0 || MatchTargID > m_Sequences.NumSeqs2Assemb)
					{
					gDiagnostics.DiagOut(eDLWarn,gszProcName,"Thread %d Couldn't find header word for %u ...",pPars->ThreadIdx,MatchTargID);
					break;
					}

				if(MatchTargID == SeqID || TargLen != ProbeLen)	// if self, or not same length, then try next matching sequence
					continue;

				if((CmpRslt = CmpPackedSeqs(pStartSeqWrd,pTarg,ProbeLen))!=0)
					break;

				TargFlags = m_Sequences.pSeqFlags[MatchTargID-1];
				if(TargFlags & cFlgNoProc)						// if already marked then try next matching sequence
					continue;
				UpdateSeqFlags(MatchTargID,ProcFlags | cFlgNoProc,0,true);
				NumOverlapping += 1;
				}
			while(!CmpRslt);
			}
		}
	UpdateSeqFlags(1 + EndingSeqID - StartingSeqID,MultiSeqFlags,true);
	}

pPars->NumProcessed += NumProcessed;
pPars->NumOverlapping += NumOverlapping;
pPars->NumOverlapped += NumOverlapped;
AcquireLock(true);
m_Sequences.NumProcessed += NumProcessed;
m_Sequences.NumOverlapping += NumOverlapping;
m_Sequences.NumOverlapped += NumOverlapped;
ReleaseLock(true);
gDiagnostics.DiagOut(eDLDebug,gszProcName,"Thread %d overlapping sequence identification completed",pPars->ThreadIdx);
return(1);		// success
}

// AddReadKMers
// Adds all unique KMer instances of length m_KMerSeqLen from pRead to m_pKMerSeqs
int									// returns number of KMers of length m_KMerSeqLen accepted from pRead, 0 if none, < 0 if errors
CArtefactReduce::AddReadKMers(int ReadLen,		// number of bases in read
				etSeqBase *pRead)				// read sequence

{
int Rslt;
int Idx;
etSeqBase *pStart;
int CurKMerLen;
int NumKMers;

if(ReadLen < m_KMerSeqLen || ReadLen < cMinKMerDistLen)
	return(0);

pStart = pRead;
CurKMerLen = 0;
NumKMers = 0;
for(Idx = 0; Idx < ReadLen; Idx++, pRead++)
	{
	if(*pRead > eBaseT)		// all accepted KMers must only contain A..T bases
		{
		CurKMerLen = 0;
		pStart = pRead+1;
		continue;
		}
	if(++CurKMerLen >= m_KMerSeqLen)
		{
		if((Rslt = AddReadKMer(pStart++)) < 1)
			return(Rslt);
		NumKMers += 1;
		}
	}
return(NumKMers);
}

int									 // returns < 0 if errors, 1 if this is the first instance of the KMer sequence, 2..N if multiple instances previously added
CArtefactReduce::AddReadKMer(etSeqBase *pKMerSeq)	// KMer bases, expected to be of at least length m_KMerSeqLen
{

UINT8 PackedSeqs[((cMaxKMerDistLen+3)/4)];	// to hold packed read bases of upto cMaxKMerDistLen in length
UINT32 NxtSeqInst;
size_t NxtSeqInstOfs;
etSeqBase *pKMerBase;
UINT8 *pPacked;
UINT8 Packed;
UINT32 KMerHash;
UINT32 *pKMerWrd;
int Idx;
tsKMerSeqInst *pKMerSeqInst;
tsKMerSeqInst *pNewKMerSeq;
int NumInstances;

pKMerBase = pKMerSeq;
pPacked = PackedSeqs;
KMerHash = 0;
Packed = 0;
for(Idx = 0; Idx < m_KMerSeqLen; Idx++,pKMerBase++)
	{
	if(!(Idx % 4))
		{
		if(Idx)
			{
			*pPacked = Packed;
			pPacked += 1;
			}
		Packed = 0;
		}
	else
		Packed <<= 2;
	Packed |= *pKMerBase;
	if(Idx == m_KMerSeqLen-1)
		*pPacked = Packed;
	KMerHash <<= 1;
	KMerHash += *pKMerBase;
	}
KMerHash &= cKMerSeqHashMask;

// check if this packed Kmer sequence has been previously processed
pKMerWrd = NULL;
AcquireSerialise();
if (m_NumKMerSeqHashes && (NxtSeqInst = m_pKMerSeqHashes[KMerHash]) != 0)		// seen at least one sequence previously with same hash?
	{
	do {
		NxtSeqInstOfs = (size_t)((NxtSeqInst - 1) * (UINT64)m_KMerSeqInstSize);
		pKMerSeqInst = (tsKMerSeqInst *)((UINT8 *)m_pKMerSeqs + NxtSeqInstOfs);
		if((UINT32 *)pKMerSeqInst->PackedSeqs == (UINT32 *)PackedSeqs)	
			{
			if(!memcmp(pKMerSeqInst->PackedSeqs,PackedSeqs,(m_KMerSeqLen + 3) / 4))
				{
				pKMerSeqInst->NumInstances += 1;
				NumInstances = (int)pKMerSeqInst->NumInstances;
				ReleaseSerialise();
				return(NumInstances);
				}
			}
		}
	while ((NxtSeqInst = pKMerSeqInst->NxtSeq) != 0);
	}

// must be the first instance of this Kmer
if(m_UsedKMerSeqInsts >= cMaxKMerSeqInsts)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Too many unique KMer instances: %u",cMaxKMerSeqInsts);
	return(eBSFerrMaxEntries);
	}

// realloc memory as required to hold this new instance
if(m_AllocdKMerSeqInsts >= m_UsedKMerSeqInsts)
	{
	tsKMerSeqInst *pAllocd;
	UINT32 ReallocInstsTo;
	ReallocInstsTo = min((size_t)(UINT64)m_AllocdKMerSeqInsts+cReallocKMerSeqInsts,cMaxKMerSeqInsts);

	size_t memreq = ReallocInstsTo * (size_t)m_KMerSeqInstSize;

#ifdef _WIN32
	pAllocd = (tsKMerSeqInst *)realloc(m_pKMerSeqs, memreq);
#else
	pAllocd = (tsKMerSeqInst *)mremap(m_pKMerSeqs,m_AllocdKMerSeqInstsMem,memreq,MREMAP_MAYMOVE);
	if(pAllocd == MAP_FAILED)
		pAllocd = NULL;
#endif
	if(pAllocd == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Memory re-allocation to %d bytes - %s",memreq,strerror(errno));
		ReleaseSerialise();
		return(eBSFerrMem);
		}
	m_pKMerSeqs = pAllocd;
	m_AllocdKMerSeqInstsMem = memreq;
	m_AllocdKMerSeqInsts = ReallocInstsTo;
	}

pNewKMerSeq = &m_pKMerSeqs[m_UsedKMerSeqInsts];
pNewKMerSeq->NumInstances = 1;
pNewKMerSeq->Flags = 0;
pNewKMerSeq->NxtSeq = 0;
memcpy(pNewKMerSeq->PackedSeqs,PackedSeqs,(m_KMerSeqLen + 3) / 4);
if (pKMerSeqInst == NULL)
	{
	 m_pKMerSeqHashes[KMerHash] = m_UsedKMerSeqInsts + 1;
	 m_NumKMerSeqHashes += 1;
	}
else
	pKMerSeqInst->NxtSeq = m_UsedKMerSeqInsts + 1;
m_UsedKMerSeqInsts += 1;
ReleaseSerialise();
return(1);
}




