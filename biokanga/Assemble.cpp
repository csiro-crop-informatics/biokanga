/*
 * CSIRO Open Source Software License Agreement (GPLv3)
 * Copyright (c) 2017, Commonwealth Scientific and Industrial Research Organisation (CSIRO) ABN 41 687 119 230.
 * See LICENSE for the complete license information (https://github.com/csiro-crop-informatics/biokanga/LICENSE)
 * Contact: Alex Whan <alex.whan@csiro.au>
 */

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
#include "./deNovoAssemb.h"


int
deNovoAssemble(etdeNovoPMode PMode,		// processing mode, currently either eAMEAssemble (default), eAMESAssemble (stringent) or eAMQAssemble (quick)
		int TrimEnds,					// trim all input sequences, both 5' and 3' ends by this many bases
		int MinSeqLen,					// only accept input sequences which are of at least this length after any trimming
		int TrimPE2SE,					// trim PEs both 5' and 3' ends by this many bases before treating as SE
		bool bAllowSE2PE,				// if true then if SE overlaps PE1 or PE2 ends singularly but not both then merge the SE with the relevant overlapped PE end and retain as a PE		
		bool bSenseStrandOnly,			// sequences from sense strand specific
		bool bSingleEnded,				// treat all sequences as being single ended even if loaded as paired ends
		int MaxPasses,					// limit number of de Novo assembly passes to this maximum (quick mode defaults to 10, exhaustive defaults to 50) set to 0 for no limit
		int PassThres,					// pass threshold at which to output intermediate assembled sequences (0 if to only write final assemblies)
		int NReduceThresSteps,			// reduce overlap thresholds over this many steps (defaults: 3 quick, 5 standard, 8 stringent assemble)");
		int Subs100bp,					// allow this many induced substitutions per 100bp overlapping sequence fragments (defaults to 1, range 0..5)
		int MaxEnd12Subs,				// allow at the initial 12bp (2 x hexamer primer artefacts) of the 5' or 3' of overlap to have this many base mismatches in addition to the overall allowed MaxSubs 
		int InitSEOvlp,					// initial minimal SE overlap required to merge SEs
		int FinSEOvlp,					// final minimal SE overlap required to merge SEs
		int InitPEOvlp,					// initial minimal PE overlap required to merge PEs
		int FinPEOvlp,					// final minimal PE overlap required to merge PEs
		int MinPE2SEOvlp,				// minimal overlap of PE1 onto PE2 required to merge as SE
		int PE2SESteps,					// when less than this many steps remaining then treat PE1 and PE2 as individual SE sequences if excessive lengths (defaults to 2, set 0 to disable)");
		int OrientatePE,				// PE end orientations 0: sense/antisense, 1: sense/sense, 2: antisense/sense, 3: antisense/antisense 
		int NumThreads,					// number of worker threads to use
		bool bAffinity,					// thread to core affinity
		char *pszPE1File,				// optional input high confidence seed PE1 sequences file
		char *pszPE2File,				// optional input high confidence seed PE2 sequences file
		char *pszSeedContigsFile,		// optional input high confidence seed SE contigs file
		char *pszInArtReducfile,		// optional input preprocessed artefact reduced packed reads from this file
		char *pszOutFile);				// where to write assembled sequence fragments ("SE","PE1","PE2" appended)


#ifdef _WIN32
int Assemble(int argc, char* argv[])
{
// determine my process name
_splitpath(argv[0],NULL,NULL,gszProcName,NULL);
#else
int
Assemble(int argc, char** argv)
{
// determine my process name
CUtility::splitpath((char *)argv[0],NULL,gszProcName);
#endif
int iScreenLogLevel;		// level of screen diagnostics
int iFileLogLevel;			// level of file diagnostics
char szLogFile[_MAX_PATH];	// write diagnostics to this file
int Rslt = 0;   			// function result code >= 0 represents success, < 0 on failure

int NumberOfProcessors;		// number of installed CPUs
int NumThreads;				// number of threads (0 defaults to number of CPUs)
bool bAffinity;				// thread to core affinity

int PMode;					// processing mode, currently either eAMEAssemble (default), eAMESAssemble (stringent) or eAMQAssemble (quick)

bool bNoAutoOvlp;           // if false then using overlap threshold defaults derived from estimated sequence lengths; true to use hard coded defaults

int TrimEnds;				// trim input sequences, both 5' and 3' ends by this many bases
int MinSeqLen;				// only accept input sequences which are of at least this length after any trimming
int TrimPE2SE;				// trim PEs both 5' and 3' ends by this many bases before treating as SE

int AllowSE2PE;				// if true then if SE overlaps PE1 or PE2 ends singularly but not both then merge the SE with the relevant overlapped PE end and retain as a PE

int SenseStrandOnly;		// sequences from sense strand specific
int SingleEnded;			// treat all sequences as being single ended even if loaded as paired ends

int MaxPasses;				// limit number of de Novo assembly passes to this maximum (quick mode defaults to 30, standard defaults to 50) set to 0 for no limit
int PassThres;				// pass threshold at which to output intermediate assembled sequences (0 if to only write final assemblies)

int Subs100bp;				// allow this many induced substitutions per 100bp overlapping sequence fragments (defaults to 1, range 0..5)
int End12Subs;				// allow at the initial 12bp (2 x hexamer primer artefacts) of the 5' or 3' of overlap to have this many base mismatches in addition to the overall allowed MaxSubs 

int InitSEOvlp;				// initial minimal SE overlap required to merge SEs
int FinSEOvlp;				// final minimal SE overlap required to merge SEs
int InitPEOvlp;				// initial minimal PE overlap required to merge PEs
int FinPEOvlp;				// final minimal PE overlap required to merge PEs
int MinPE2SEOvlp;			// minimal overlap of PE1 onto PE2 required to merge as SE
int NReduceThresSteps;		// reduce overlap thresholds over this many steps (defaults: 3 quick, 5 standard, 8 stringent assemble, range 2..10)");
int PE2SESteps;				// when less than this many steps remaining then treat PE1 and PE2 as individual SE sequences if excessive lengths (defaults to 2, set 0 to disable)");

int OrientatePE;		    // PE end orientations 0: sense/antisense, 1: sense/sense, 2: antisense/sense, 3: antisense/antisense 

char szInArtReducfile[_MAX_PATH];	// preprocessed artefact reduced packed reads from this file
char szOutFile[_MAX_PATH];			// assemblies to this prefix file

char szSeedContigsFile[_MAX_PATH];	// optional high confidence seed contig sequences file
char szPE1File[_MAX_PATH];			// optional high confidence seed contig sequences file
char szPE2File[_MAX_PATH];			// optional high confidence seed contig sequences file

char szSQLiteDatabase[_MAX_PATH];	// results summaries to this SQLite file
char szExperimentName[cMaxDatasetSpeciesChrom+1];			// experiment name
char szExperimentDescr[1000];		// describes experiment

// command line args
struct arg_lit  *help    = arg_lit0("h","help",                 "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to screen and logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"diagnostics log file");

struct arg_int *pmode = arg_int0("m","mode","<int>",		    "processing mode: 0 - standard, 1 - higher stringency, 2 - quick de Novo assembly");

struct arg_lit  *noautoovlp = arg_lit0("l", "noautoovlp",		"if not specified then using overlap threshold defaults derived from estimated sequence lengths; otherwise use hard coded defaults");

struct arg_int *trimends = arg_int0("t","trimends","<int>",     "when loading reads or high confidence seed contigs then trim 5' and 3' ends by this many bases (default 0, range 1..50)");
struct arg_int *minseqlen = arg_int0("X","minseqlen","<int>",   "only accept reads or high confidence seed contigs, after any trimming, of at least this length (default 70, range 50..500)");

struct arg_int *trimpe2se = arg_int0("x","trimpe2se","<int>",	"trim PEs both 5' and 3' ends by this many bases when treating as SE (default 10, range 0..50)");

struct arg_lit  *sensestrandonly = arg_lit0("E","senseonly",    "process sequences as strand specific");
struct arg_lit  *singleended    = arg_lit0("e","singleend",     "process all paired ends as if single ended");

struct arg_int *subs100bp = arg_int0("s","maxsubs100bp","<int>",  "allow max induced substitutions per 100bp overlapping sequence fragments (defaults to 1, range 0..5)");
struct arg_int *end12subs = arg_int0("S","maxendsubs","<int>",    "allow max induced substitutions in overlap 12bp ends (defaults to 0, range 0..6)");

struct arg_int *initseovlp = arg_int0("j","initseovlp","<int>",   "initial minimal SE overlap required to merge SEs (defaults to 150, range 20..500)");
struct arg_int *finseovlp = arg_int0("J","finseovlp","<int>",     "final minimal SE overlap required to merge SEs (defaults to 25, range 20..initseovlp)");
struct arg_int *initpeovlp = arg_int0("k","initpeovlp","<int>",  "initial minimal PE total sum of end overlaps required to merge PEs (defaults to 150, range 35..200)");
struct arg_int *finpeovlp = arg_int0("K","finpeovlp","<int>",    "final minimal PE total sum of end overlaps required to merge PEs (defaults to 35, range 35..initpeovlp)");
struct arg_int *minpe2seovlp = arg_int0("g","minpe2seovlp","<int>", "minimal overlap of PE1 onto PE2 required to merge as SE (defaults to 20, range 16..100)");
struct arg_int *reducethressteps = arg_int0("r","reducethressteps","<int>", "reduce overlap thresholds over this many steps (defaults: 3 quick, 5 standard, 8 stringent assemble, range 2..10)");
struct arg_int  *pe2sesteps = arg_int0("R","pe2sesteps","<int>",   "when less or equal to this remaining threshold steps then treat PE1 and PE2 as individual SE sequences if excessive lengths (defaults to 2, set 0 to disable)");

struct arg_int *passthres = arg_int0("P","passthres","<int>",   "de Novo assembly process pass threshold at which to start writing intermediate checkpoint assemblies (defaults to 0, only write completed assemblies)");
struct arg_int *maxpasses = arg_int0("p","maxpasses","<int>",   "limit number of de Novo assembly processing passes to this maximum (defaults: standard 50, stringent 75, quick 30) range 20..10000");

struct arg_lit  *allowse2pe = arg_lit0("Z","allowse2pe",         "if SE overlaps PE1 or PE2 ends singularly but not both then merge the SE with the relevant overlapped PE end and retain as a PE (default no)");

struct arg_int *orientatepe = arg_int0("M","orientatepe","<int>",    "read pair end orientations, 0: sense/antisense (PE short insert), 1: sense/sense (MP Roche 454), 2: antisense/sense (MP Illumina circularized), 3: antisense/antisense (MP SOLiD)");

struct arg_file *inpe1file = arg_file0("a","inpe1","<file>","optionally Load 5' paired end previously assembled fragments or filtered PE1 reads from fasta file");
struct arg_file *inpe2file = arg_file0("A","inpe2","<file>","optionally Load 3' paired end previously assembled fragments or filtered PE2 reads from fasta file");

struct arg_file *inartreducfile = arg_file0("i","inartreducfile","<file>",	"load from previously generated artefact reduced packed reads file");

struct arg_file *seedcontigsfile = arg_file0("c","seedcontigsfile","<file>", "optionally load high confidence fasta seed contigs or previously assembled fasta SE sequences file");

struct arg_file *outfile = arg_file0("o","out","<file>",					"Output assembled contigs to this file");

struct arg_file *summrslts = arg_file0("q","sumrslts","<file>",				"Output results summary to this SQLite3 database file");
struct arg_str *experimentname = arg_str0("w","experimentname","<str>",		"experiment name SQLite3 database file");
struct arg_str *experimentdescr = arg_str0("W","experimentdescr","<str>",	"experiment description SQLite3 database file");

struct arg_int *threads = arg_int0("T","threads","<int>",		"number of processing threads 0..128 (defaults to 0 which sets threads to number of CPU cores)");

struct arg_end *end = arg_end(200);

void *argtable[] = {help,version,FileLogLevel,LogFile,
	                pmode,noautoovlp,trimends,minseqlen,trimpe2se,allowse2pe,sensestrandonly,singleended,maxpasses,reducethressteps,passthres,subs100bp,end12subs,
					initseovlp,finseovlp,initpeovlp,finpeovlp,minpe2seovlp,pe2sesteps,
					orientatepe,inpe1file,inpe2file,seedcontigsfile,inartreducfile,outfile,
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
		printf("\n      To invoke this parameter file then precede its name with '@'");
		printf("\n      e.g. %s %s @myparams.txt\n",gszProcName,gpszSubProcess->pszName);
		printf("\nPlease report any issues regarding usage of %s at https://github.com/csiro-crop-informatics/biokanga/issues\n\n",gszProcName);
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


	PMode = pmode->count ? pmode->ival[0] : (int)eAMEAssemble;
	if(PMode < 0 || PMode > eAMQAssemble)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Processing mode '-m%d' must be in range 0..%d",PMode,eAMQAssemble);
		return(1);
		}

	bNoAutoOvlp = noautoovlp->count ? true : false;

	// can have various combinations of input files
	// if pe inputs then both 5' and 3' files must be specified
	// se file is optional
	// packed reads file is optional
	// BUT at least one of PE or SE or Packed must be specified
	if (inartreducfile->count)
		{
		strcpy(szInArtReducfile, inartreducfile->filename[0]);
		CUtility::TrimQuotedWhitespcExtd(szInArtReducfile);
		}
	else
		szInArtReducfile[0] = '\0';

	if (inpe1file->count)
		{
		strcpy(szPE1File, inpe1file->filename[0]);
		CUtility::TrimQuotedWhitespcExtd(szPE1File);
		}
	else
		szPE1File[0] = '\0';

	if (inpe2file->count)
		{
		strcpy(szPE2File, inpe2file->filename[0]);
		CUtility::TrimQuotedWhitespcExtd(szPE2File);
		}
	else
		szPE2File[0] = '\0';
	if ((szPE1File[0] == '\0' && szPE2File[0] != '\0') ||
		(szPE1File[0] != '\0' && szPE2File[0] == '\0'))
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: If either PE previously assembled fragments file specified then both PE files must be specified");
		return(1);
		}

	if (seedcontigsfile->count)
		{
		strcpy(szSeedContigsFile, seedcontigsfile->filename[0]);
		CUtility::TrimQuotedWhitespcExtd(szSeedContigsFile);
		}
	else
		szSeedContigsFile[0] = '\0';

	// check that if artifact reduced is specified then the PE files are not allowed to also be specified
	if (szInArtReducfile[0] != '\0')
		{
		if (szPE1File[0] != '\0')
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: If artefact reduced file '-i%s' specified then PE1/PE2 can't also be specified", szInArtReducfile);
			return(1);
			}
		}

	// check that PE and/or SE and/or Packed have been specified
	if (szPE1File[0] == '\0' && szInArtReducfile[0] == '\0' && szSeedContigsFile[0] == '\0')
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: At least one of PE '-a and -A' or SE '-c' or artefact reduced packed '-i' files must be specified");
		return(1);
		}

	if (!outfile->count)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: Expected output file '-o<outfile>' to be specified");
		return(1);
		}
	strncpy(szOutFile, outfile->filename[0], _MAX_PATH);
	szOutFile[_MAX_PATH - 1] = '\0';



	TrimEnds = 0;
	MinSeqLen = 70;
	TrimPE2SE = 10;
	InitPEOvlp = cDfltInitPEOvlp;
	FinPEOvlp = cDfltFinPEOvlp;
	MinPE2SEOvlp = cDfltMinPE1PE2ToSEOvlp;
	InitSEOvlp = cDfltInitSEOvlp;
	FinSEOvlp = cDfltFinSEOvlp;

	if(!bNoAutoOvlp)
		{
		CFasta Fasta;
		UINT32 PE1EstNumSeqs;
		INT32 PE1EstMaxSeqLen;
		INT32 PE1EstMeanSeqLen;
		UINT32 PE2EstNumSeqs;
		INT32 PE2EstMaxSeqLen;
		INT32 PE2EstMeanSeqLen;
		UINT32 SEEstNumSeqs;
		INT32 SEEstMaxSeqLen;
		INT32 SEEstMeanSeqLen;

		PE1EstNumSeqs = 0;
		PE1EstMaxSeqLen = 0;
		PE1EstMeanSeqLen = 0;
		PE2EstNumSeqs = 0;
		PE2EstMaxSeqLen = 0;
		PE2EstMeanSeqLen = 0;
		SEEstNumSeqs = 0;
		SEEstMaxSeqLen = 0;
		SEEstMeanSeqLen = 0;
		if(szPE1File[0] != '\0')
			{
			if((PE1EstNumSeqs = Fasta.FastaEstSizes(szPE1File,NULL,NULL,NULL,&PE1EstMaxSeqLen,&PE1EstMeanSeqLen)) == 0)
				{
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: Unable to estimate sequence sizes for PE1 file: '%s'", szPE1File);
				return(1);
                }

			if((PE2EstNumSeqs = Fasta.FastaEstSizes(szPE1File, NULL, NULL, NULL, &PE2EstMaxSeqLen, &PE2EstMeanSeqLen)) == 0)
				{
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: Unable to estimate sequence sizes for PE2 file: '%s'", szPE2File);
				return(1);
				}

			MinSeqLen = (80 * (PE1EstMeanSeqLen + PE2EstMeanSeqLen)) / 200;  // use 80% of mean sequence lengths as the default minimum
			if(MinSeqLen < 50)
				MinSeqLen = 50;
			else
				if(MinSeqLen > 500)
					MinSeqLen = 500;

			TrimPE2SE = min(MinSeqLen / 20, 50);							// trimming by 5% of sequence length9

			InitPEOvlp = max(cMinPEOvlp,(80 * (PE1EstMeanSeqLen + PE2EstMeanSeqLen)) / 100);       // requiring initially 80% overlap summed over both ends
			if(InitPEOvlp > cMaxPEOvlp)
				InitPEOvlp = cMaxPEOvlp;


			FinPEOvlp = max(cMinPEOvlp, (30 * (PE1EstMeanSeqLen + PE2EstMeanSeqLen)) / 100);		// requiring finally 30% overlap summed over both ends
			if (FinPEOvlp > InitPEOvlp)
				FinPEOvlp = InitPEOvlp;
		
			MinPE2SEOvlp = max(16,MinSeqLen / 4);
			if(MinPE2SEOvlp > 100)
				MinPE2SEOvlp = 100;
		
			InitSEOvlp = (1 + PE1EstMeanSeqLen + PE2EstMeanSeqLen) / 2;

			}
		else
			{
			if (szSeedContigsFile[0] != '\0')
				{
				if((SEEstNumSeqs = Fasta.FastaEstSizes(szSeedContigsFile, NULL, NULL, NULL, &SEEstMaxSeqLen, &SEEstMeanSeqLen)) == 0)
					{
					gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: Unable to estimate sequence sizes for seed contigs file: '%s'", szSeedContigsFile);
					return(1);
					}
				}
			InitSEOvlp = (SEEstMeanSeqLen + 2) / 3;
			}

 
		if(InitSEOvlp > cMaxSEOvlp)
			InitSEOvlp = cMaxSEOvlp;
		else
			if (InitSEOvlp < cMinSEOvlp)
				InitSEOvlp = cMinSEOvlp;

		FinSEOvlp = (InitSEOvlp + 3) / 4;        
		if (FinSEOvlp > InitSEOvlp)
			FinSEOvlp = InitSEOvlp;
		else
			if (FinSEOvlp < cMinSEOvlp)
				FinSEOvlp = cMinSEOvlp;
		}

	AllowSE2PE = allowse2pe->count ? 1 : 0; 
	TrimEnds = trimends->count ? trimends->ival[0] : TrimEnds;
	if(TrimEnds < 0 || TrimEnds > 50)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Expected end trimming '-t%d' to be in range 0..50",TrimEnds);
		return(1);
		} 
	MinSeqLen = minseqlen->count ? minseqlen->ival[0] : MinSeqLen;
	if(MinSeqLen < 50 || MinSeqLen > 500)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Expected minimum sequence length '-X%d' to be in range 50..500",MinSeqLen);
		return(1);
		} 

	TrimPE2SE = trimpe2se->count ? trimpe2se->ival[0] : TrimPE2SE;
	if(TrimPE2SE < 0 || TrimPE2SE > 50)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Expected PE to SE end trimming '-x%d' to be in range 0..50",TrimPE2SE);
		return(1);
		} 

	SenseStrandOnly = sensestrandonly->count ? (int)true : int(false);
	SingleEnded = singleended->count ? (int)true : int(false);

	OrientatePE = orientatepe->count ? orientatepe->ival[0] : 0;
	if(OrientatePE < 0 || OrientatePE > 3)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Expected paired end orientation '-M%d' to be in range 0..3",OrientatePE);
		return(1);
		}

	MaxPasses = maxpasses->count ? maxpasses->ival[0] : 0;
	if((MaxPasses != 0) && (MaxPasses < 20 || MaxPasses > 1000))
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Max allowed passes '-p%d' must be in range 20..1000",MaxPasses);
		return(1);
		}
	if(MaxPasses == 0)
		{
		switch(PMode) {
			case eAMEAssemble:		// standard de Novo assemble
				MaxPasses = 50;
				break;

			case eAMESAssemble:		// more stringent de Novo assemble
				MaxPasses = 75;
				break;

			case eAMQAssemble:		// quick assemble packed reads with low stringency
				MaxPasses = 30;
				break;
			}
		}

	NReduceThresSteps = reducethressteps->count ? reducethressteps->ival[0] : 0;
	if((NReduceThresSteps != 0) && (NReduceThresSteps < 2 || NReduceThresSteps > 10))
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Threshold reduction steps '-r%d' must be in range 2..10",NReduceThresSteps);
		return(1);
		}

	if(NReduceThresSteps == 0)		// if not specified then use defaults
		{
		switch(PMode) {
			case eAMEAssemble:		// standard de Novo assemble
				NReduceThresSteps = 5;
				break;

			case eAMESAssemble:		// more stringent de Novo assemble
				NReduceThresSteps = 8;
				break;

			case eAMQAssemble:		// quick assemble packed reads with low stringency
				NReduceThresSteps = 3;
				break;
			}
		}

	PassThres = passthres->count ? passthres->ival[0] : 0;
	if((PassThres != 0) && (PassThres < 1 || PassThres > MaxPasses))
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Intermediate assemblies '-P%d' must be in range 1..%d",PassThres,MaxPasses);
		return(1);
		}

	Subs100bp = subs100bp->count ? subs100bp->ival[0] : 1;
	if(Subs100bp < 0 || Subs100bp > 5)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Max induced substitutions '-s%d' per 100bp overlapping must be in range 0..5",Subs100bp);
		return(1);
		}

	End12Subs = end12subs->count ? end12subs->ival[0] : 0;
	if(End12Subs < 0 || End12Subs > 6)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Max induced substitutions in overlap 12bp ends '-E%d' must be in range 0..6",End12Subs);
		return(1);
		}


	if(szPE1File[0] != '\0' && szPE2File[0] != '\0') // if processing PEs then PE parameterisation needs to be parsed
		{
		InitPEOvlp = initpeovlp->count ? initpeovlp->ival[0] : InitPEOvlp;
		if(InitPEOvlp < cMinPEOvlp || InitPEOvlp > cMaxPEOvlp)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Initial min required PE merge overlap '-k%d' must be in range %d..%d",InitPEOvlp,cMinPEOvlp,cMaxPEOvlp);
			return(1);
			}

		FinPEOvlp = finpeovlp->count ? finpeovlp->ival[0] : FinPEOvlp;
		if(FinPEOvlp < cMinPEOvlp || FinPEOvlp > InitPEOvlp)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Final min required PE merge overlap '-K%d' must be in range %d..%d",FinPEOvlp,cMinPEOvlp,InitPEOvlp);
			return(1);
			}

		MinPE2SEOvlp = minpe2seovlp->count ? minpe2seovlp->ival[0] : MinPE2SEOvlp;
		if(MinPE2SEOvlp < cMinPE1PE2ToSEOvlp || MinPE2SEOvlp > cMaxPE1PE2ToSEOvlp)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: minimal overlap of PE1 end onto PE2 end required to merge as SE '-g%d' must be in range %d..%d",MinPE2SEOvlp,cMinPE1PE2ToSEOvlp,cMaxPE1PE2ToSEOvlp);
			return(1);
			}

		PE2SESteps = pe2sesteps->count ? pe2sesteps->ival[0] : min(2,NReduceThresSteps);
		if(PE2SESteps < 0 || PE2SESteps > NReduceThresSteps)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Remaining steps before excessive PE end length '-z%d' checking starts must be in range 1..%d",PE2SESteps,NReduceThresSteps);
			return(1);
			}
		}
	else
		{
		InitPEOvlp = 0;
		FinPEOvlp = 0;
		MinPE2SEOvlp = 0;
		PE2SESteps = 0;
		}

	InitSEOvlp = initseovlp->count ? initseovlp->ival[0] : InitSEOvlp;
	if(InitSEOvlp < cMinSEOvlp || InitSEOvlp > cMaxSEOvlp)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: initial minimal SE overlap required to merge SEs '-j%d' must be in the range %d..%d",InitSEOvlp,cMinSEOvlp,cMaxSEOvlp);
			return(1);
			}
	FinSEOvlp = finseovlp->count ? finseovlp->ival[0] : min(cDfltFinSEOvlp,InitSEOvlp);
	if(FinSEOvlp < cMinSEOvlp || FinSEOvlp > InitSEOvlp)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: initial minimal SE overlap required to merge SEs '-j%d' must be in the range %d..%d",FinSEOvlp,cMinSEOvlp,InitSEOvlp);
			return(1);
			}

// show user current resource limits
#ifndef _WIN32
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Resources: %s",CUtility::ReportResourceLimits());
#endif

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
	switch((etdeNovoPMode)PMode) {
		case eAMEAssemble:
			pszDescr = "standard de Novo assemble";
			break;
		case eAMESAssemble:
			pszDescr = "high stringency de Novo assemble";
			break;
		case eAMQAssemble:
			pszDescr = "quick de Novo assemble";
			break;
		}

	gDiagnostics.DiagOutMsgOnly(eDLInfo, "Using overlap threshold defaults derived from estimated sequence lengths : '%s'",bNoAutoOvlp ? "No" : "Yes");

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Processing mode is : '%s'",pszDescr);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"End trimming by: %dbp",TrimEnds);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Accept input sequences, after any trimming, which are at least: %dbp",MinSeqLen);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"PE to SE end trimming by: %dbp",TrimPE2SE);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Allow SE conversion into PE: '%s'", AllowSE2PE == 0 ? "No" : "Yes");

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Process sequences as strand specific: %s",SenseStrandOnly ? "Yes" : "No");
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Process sequences as always single end: %s",SingleEnded ? "Yes" : "No");
	
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Initial minimal SE overlap required to merge SEs: %d",InitSEOvlp);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Final minimal SE overlap required to merge SEs: %d",FinSEOvlp);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Initial minimal sum of PE end overlaps required to merge PEs: %d",InitPEOvlp);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Final minimal sum of PE end overlaps required to merge PEs: %d",FinPEOvlp);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Minimal overlap of PE1 onto PE2 required to merge as SE: %d",MinPE2SEOvlp);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Limit number of de Novo assembly processing passes to: %d",MaxPasses);

	if(PassThres)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Output intermediate assemblies to file starting from pass: %d",PassThres);
	else
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"No intermediate assemblies output to file");

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Allow max induced substitutions per 100bp overlapping sequence fragments: %d",Subs100bp);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Allow max induced substitutions end 12bp of overlaps: %d",End12Subs);


	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Threhold reduction steps: %d",NReduceThresSteps);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Remaining steps before excessive PE end length checking: %d",PE2SESteps);

	if(szInArtReducfile[0] != '\0')
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Input artefact reduced packed reads file : '%s'",szInArtReducfile);

	if(szPE1File[0] != '\0')
		{
		char *pszPE1Orientate;
		char *pszPE2Orientate;
		switch(OrientatePE) {
			case 0:			// the default is for standard short inserts with PE1 sense and PE2 antisense
				pszPE1Orientate = (char *)"Sense";			// ---->   <-----
				pszPE2Orientate = (char *)"Antisense";
				break;
			case 1:			// mate pair Roche 454  // ---->   ----->
				pszPE1Orientate = (char *)"Sense";
				pszPE2Orientate = (char *)"Sense";
				break;
			case 2:			// mate pair Illuminia circularised library 
				pszPE1Orientate = (char *)"Antisense";      // <----   ------>
				pszPE2Orientate = (char *)"Sense";
				break;
			case 3:			// mate pair SOLiD      // <-----  <------
				pszPE1Orientate = (char *)"Antisense";
				pszPE2Orientate = (char *)"Antisense";
				break;
			}
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"5' PE1 is : '%s' and 3' PE2 is : '%s'",pszPE1Orientate,pszPE2Orientate);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Input 5' paired end or previously assembled fragments file : '%s'",szPE1File);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Input 3' paired end or previously assembled fragments file : '%s'",szPE2File);
		}

	if(szSeedContigsFile[0] != '\0')
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"input high confidence seed contig or previously assembled SE sequences file: '%s'",szSeedContigsFile);
	
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Output assembly multifasta file : '%s'",szOutFile);

	if(szExperimentName[0] != '\0')
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"This processing reference: %s",szExperimentName);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Number of threads : %d",NumThreads);

	if(gExperimentID > 0)
		{
		int ParamID;
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTText,(int)strlen(szLogFile),"log",szLogFile);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTInt32,sizeof(PMode),"mode",&PMode);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTInt32,sizeof(TrimEnds),"trimends",&TrimEnds);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTInt32,sizeof(MinSeqLen),"minseqlen",&MinSeqLen);

		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTInt32,sizeof(TrimPE2SE),"trimpe2se",&TrimPE2SE);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTInt32,sizeof(SenseStrandOnly),"sensestrandonly",&SenseStrandOnly);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTInt32,sizeof(SingleEnded),"singleended",&SingleEnded);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTInt32,sizeof(MaxPasses),"maxpasses",&MaxPasses);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTInt32,sizeof(PassThres),"passthres",&PassThres);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTInt32,sizeof(Subs100bp),"subs100bp",&Subs100bp);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTInt32,sizeof(End12Subs),"maxendsubs",&End12Subs);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTInt32,sizeof(NReduceThresSteps),"reducethressteps",&NReduceThresSteps);

		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTInt32,sizeof(InitSEOvlp),"initseovlp",&InitSEOvlp);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTInt32,sizeof(FinSEOvlp),"finseovlp",&FinSEOvlp);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTInt32,sizeof(InitPEOvlp),"initpeovlp",&InitPEOvlp);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTInt32,sizeof(FinPEOvlp),"finpeovlp",&FinPEOvlp);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTInt32,sizeof(MinPE2SEOvlp),"minpe2seovlp",&MinPE2SEOvlp);

		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTInt32,sizeof(PE2SESteps),"pe2sesteps",&PE2SESteps);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTInt32,sizeof(AllowSE2PE),"allowse2pe",&AllowSE2PE);

		if(szSeedContigsFile[0] != '\0')
			ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTText,(int)strlen(szSeedContigsFile),"seedcontigsfile",szSeedContigsFile);
		if(szInArtReducfile[0] != '\0')
			ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTText,(int)strlen(szInArtReducfile),"in",szInArtReducfile);
		if(szPE1File[0] != '\0')
			{
			ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTInt32,sizeof(int),"orientatepe",&OrientatePE);
			ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTText,(int)strlen(szPE1File),"inpe1",szPE1File);
			ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTText,(int)strlen(szPE2File),"inpe2",szPE2File);
			}

		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTText,(int)strlen(szOutFile),"out",szOutFile);

		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTInt32,sizeof(NumThreads),"threads",&NumThreads);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTInt32,sizeof(NumberOfProcessors),"cpus",&NumberOfProcessors);

		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTText,(int)strlen(szSQLiteDatabase),"sumrslts",szSQLiteDatabase);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTText,(int)strlen(szExperimentName),"experimentname",szExperimentName);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTText,(int)strlen(szExperimentDescr),"experimentdescr",szExperimentDescr);
		}

#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	gStopWatch.Start();
	Rslt = deNovoAssemble((etdeNovoPMode)PMode,TrimEnds,MinSeqLen,TrimPE2SE,AllowSE2PE == 0 ? false : true,SenseStrandOnly ? true : false,SingleEnded ? true : false,MaxPasses,PassThres,NReduceThresSteps,Subs100bp,End12Subs,
							InitSEOvlp,FinSEOvlp,InitPEOvlp,FinPEOvlp,MinPE2SEOvlp,PE2SESteps,OrientatePE,NumThreads,bAffinity,szPE1File,szPE2File,szSeedContigsFile,szInArtReducfile,szOutFile);
	Rslt = Rslt >=0 ? 0 : 1;
	if(gExperimentID > 0)
		{
		if(gProcessingID)
			gSQLiteSummaries.EndProcessing(gExperimentID, gProcessingID,Rslt);
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
deNovoAssemble(etdeNovoPMode PMode,			// processing mode, currently either eAMEAssemble (default), eAMESAssemble (stringent) or eAMQAssemble (quick)
		int TrimEnds,						// trim input sequences, both 5' and 3' ends by this many bases
		int MinSeqLen,						// only accept input sequences which are of at least this length after any trimming
		int TrimPE2SE,						// trim PEs both 5' and 3' ends by this many bases before treating as SE
		bool bAllowSE2PE,					// if true then if SE overlaps PE1 or PE2 ends singularly but not both then merge the SE with the relevant overlapped PE end and retain as a PE		
		bool bSenseStrandOnly,				// sequences from sense strand specific
		bool bSingleEnded,					// treat all sequences as being single ended even if loaded as paired ends
		int MaxPasses,						// limit number of de Novo assembly passes to this maximum (quick mode defaults to 10, exhaustive defaults to 100) set to 0 for no limit
		int PassThres,						// pass threshold at which to output intermediate assembled sequences (0 if to only write final assemblies)
		int NReduceThresSteps,				// reduce overlap thresholds over this many steps (defaults: 3 quick, 5 standard, 8 stringent assemble)");
		int Subs100bp,						// allow this many induced substitutions per 100bp overlapping sequence fragments (defaults to 0, range 1..5)
		int MaxEnd12Subs,					// allow at the initial 12bp (2 x hexamer primer artefacts) of the 5' or 3' of overlap to have this many base mismatches in addition to the overall allowed MaxSubs 
		int InitSEOvlp,						// initial minimal SE overlap required to merge SEs
		int FinSEOvlp,						// final minimal SE overlap required to merge SEs
		int InitPEOvlp,						// initial minimal PE overlap required to merge PEs
		int FinPEOvlp,						// final minimal PE overlap required to merge PEs
		int MinPE2SEOvlp,					// minimal overlap of PE1 onto PE2 required to merge as SE
		int PE2SESteps,						// when less than this many steps remaining then treat PE1 and PE2 as individual SE sequences if excessive lengths (defaults to 2, set 0 to disable)");
		int OrientatePE,					// PE end orientations 0: sense/antisense, 1: sense/sense, 2: antisense/sense, 3: antisense/antisense 
		int NumThreads,						// number of worker threads to use
		bool bAffinity,						// thread to core affinity
		char *pszPE1File,					// optional input high confidence seed PE1 sequences file
		char *pszPE2File,					// optional input high confidence seed PE2 sequences file
		char *pszSeedContigsFile,			// optional input high confidence seed SE contigs file
		char *pszInArtReducfile,			// optional input preprocessed artefact reduced packed reads from this file
		char *pszOutFile)					// where to write assembled sequence fragments ("SE","PE1","PE2" appended)
{
int Rslt;
int SeqWrdBytes;
INT64 CumulativeMemory;
UINT32 CumulativeSequences;

CdeNovoAssemb *pAssemble;
pAssemble = new CdeNovoAssemb();
pAssemble->Reset(false);
pAssemble->SetPMode(PMode);
pAssemble->SetNumThreads(NumThreads,bAffinity);
pAssemble->SetSfxSparsity(eSSparsity15);
SeqWrdBytes = pAssemble->GetSeqWrdBytes();

// if not loading artefact reduced reads then need to preallocate memory 
if(pszInArtReducfile == NULL || pszInArtReducfile[0] == '\0')
	{
	CumulativeMemory = 0;
	CumulativeSequences = 0;

	if(pszPE1File != NULL && pszPE1File[0] != '\0')
		{
		if((Rslt = pAssemble->EstMemReq(pszPE1File)) != eBSFSuccess)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to estimate memory required for '%s",pszPE1File);
			return(eBSFerrOpnFile);	// treat as though unable to open file
			}
		if((Rslt = pAssemble->EstMemReq(pszPE2File)) != eBSFSuccess)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to estimate memory required for '%s",pszPE1File);
			return(eBSFerrOpnFile);	// treat as though unable to open file
			}
		}

	if(pszSeedContigsFile != NULL && pszSeedContigsFile[0] != '\0')
		{
		if((Rslt = pAssemble->EstMemReq(pszSeedContigsFile)) != eBSFSuccess)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to estimate memory required for '%s",pszSeedContigsFile);
			return(eBSFerrOpnFile);	// treat as though unable to open file
			}
		}

	CumulativeMemory = pAssemble->EstMemReq(TrimEnds,TrimEnds,0,0,0);
	CumulativeSequences = pAssemble->GetEstSeqsToProc();
	CumulativeMemory += CumulativeSequences * 12; // very rough estimate allowing for sparse suffix and flags requirements
	// NOTE: seems there may be a problem with memory allocations if estimating memory of readsets containing large variations in read lengths 
	// Following hack is to increase allocation size by 50% to see if that fixes the problem
	CumulativeMemory *= 3;
	CumulativeMemory /= 2;
	// end of temp hack

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
	if((Rslt = pAssemble->AllocSeqs2AssembMem((CumulativeMemory * 105)/100))!= eBSFSuccess)	// add 5% for safety...
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to continue");
		pAssemble->Reset(false);
		delete pAssemble;
		return(Rslt);
		}
	}

if((Rslt = pAssemble->AssembReads(PMode,TrimEnds,MinSeqLen,TrimPE2SE,bSenseStrandOnly,bAllowSE2PE,bSingleEnded,MaxPasses,PassThres,NReduceThresSteps,Subs100bp * 10,MaxEnd12Subs,
									InitSEOvlp,FinSEOvlp,InitPEOvlp,FinPEOvlp,MinPE2SEOvlp,PE2SESteps,
									OrientatePE,pszPE1File,pszPE2File,pszSeedContigsFile,pszInArtReducfile,pszOutFile)) < eBSFSuccess)
	return(Rslt);

// write out sequences here
if(Rslt >= eBSFSuccess)
	Rslt = pAssemble->SaveAssembSeqs(pszOutFile,0);

delete pAssemble;
return(0);
}
