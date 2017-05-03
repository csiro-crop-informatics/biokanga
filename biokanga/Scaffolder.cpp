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
#include "./AssembGraph.h"

#include "./Scaffolder.h"


int
ProcessScaffoldedContigs(int PMode,			// processing mode: 0 - output scaffold multifasta with edge report
		int NumThreads,						// number of worker threads to use
		bool bAffinity,						// thread to core affinity
		int Subs100bp,						// allow this many induced substitutions per 100bp overlapping sequence fragments (defaults to 0, range 0..5)
		int MaxEnd12Subs,					// allow at the initial 12bp (2 x hexamer primer artefacts) of the 5' or 3' of overlap to have this many base mismatches in addition to the overall allowed MaxSubs 
		int MinPEReadlen,					// only accept PE read lengths of at least this many bp
		int MinPEInsertSize,				// PE sequences are of this minimum insert size
		int MaxPEInsertSize,				// PE sequences are of this maximum insert size
		int MinScaffoldedSeqLen,			// reported scaffolded sequences are at least this length
		int OrientatePE,					// PE end orientations 0: sense/antisense, 1: sense/sense, 2: antisense/sense, 3: antisense/antisense 
		char *pszPE1File,					// input PE1 sequences file
		char *pszPE2File,					// input PE2 sequences file
		char *pszContigsFile,				// input SE contigs file
		char *pszScaffoldedFile);			// where to write scaffolded contigs


#ifdef _WIN32
int ScaffoldContigs(int argc, char* argv[])
{
// determine my process name
_splitpath(argv[0],NULL,NULL,gszProcName,NULL);
#else
int
ScaffoldContigs(int argc, char** argv)
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

int PMode;					// processing mode
int Subs100bp;				// allow this many induced substitutions per 100bp overlapping sequence fragments (defaults to 0, range 0..5)
int End12Subs;				// allow at the initial 12bp (2 x hexamer primer artefacts) of the 5' or 3' of overlap to have this many base mismatches in addition to the overall allowed MaxSubs 
int MinPEReadlen;			// only accept PE read lengths of at least this many bp
int MinPEInsertSize;		// PE sequences are of this minimum insert size
int MaxPEInsertSize;		// PE sequences are of this maximum insert size

int MinScaffoldedSeqLen;	// reported scaffolded sequences are to be at least this length

int OrientatePE;			// PE end orientations 0: sense/antisense, 1: sense/sense, 2: antisense/sense, 3: antisense/antisense 

char szOutFile[_MAX_PATH];			// scaffolded contigs to this file

char szContigsFile[_MAX_PATH];		// input SE contigs file
char szPE1File[_MAX_PATH];			// input PE1 sequences file
char szPE2File[_MAX_PATH];			// input PE2 sequences file

char szSQLiteDatabase[_MAX_PATH];	// results summaries to this SQLite file
char szExperimentName[cMaxDatasetSpeciesChrom+1];			// experiment name
char szExperimentDescr[1000];		// describes experiment

// command line args
struct arg_lit  *help    = arg_lit0("h","help",                 "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to screen and logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"diagnostics log file");

struct arg_int *pmode = arg_int0("m","mode","<int>",		    "processing mode: 0 - default");

struct arg_int *subs100bp = arg_int0("s","maxsubs100bp","<int>", "allow max induced substitutions per 100bp overlapping sequence fragments (defaults to 0, range 1..5)");
struct arg_int *end12subs = arg_int0("E","maxendsubs","<int>",   "allow max induced substitutions in overlap 12bp ends (defaults to 0, range 0..6)");
struct arg_int *minpereadlen = arg_int0("L", "minpereadlen", "<int>", "min length of any PE read (defaults to 80, range 30..5000)");

struct arg_int *minpeinsertsize = arg_int0("p","minpeinsert","<int>","minimum PE insert size (default 110, range 100..50000)");
struct arg_int *maxpeinsertsize = arg_int0("P","maxpeinsert","<int>","maximum PE insert size (default 1000, range minpeinsert..50000)");
struct arg_int *minscafflen = arg_int0("l","minscafflen","<int>",    "report scaffolded sequences of at least this minimum length (default 300, range 100 .. 5000)");

struct arg_int *orientatepe = arg_int0("M","orientatepe","<int>",    "read pair end orientations, 0: sense/antisense (PE short insert), 1: sense/sense (MP Roche 454), 2: antisense/sense (MP Illumina circularized), 3: antisense/antisense (MP SOLiD)");

struct arg_file *inpe1file = arg_file0("a","inpe1","<file>",		 "Load 5' paired end PE1 reads or partially assembled fragments from fasta file");
struct arg_file *inpe2file = arg_file0("A","inpe2","<file>",		 "Load 3' paired end PE1 reads or partially assembled fragments from fasta file");

struct arg_file *contigsfile = arg_file1("c","contigsfile","<file>", "load SE assembled contigs fasta file");

struct arg_file *outfile = arg_file1("o","out","<file>",			 "Output scaffolded contigs to this file");

struct arg_file *summrslts = arg_file0("q","sumrslts","<file>",		 "Output results summary to this SQLite3 database file");
struct arg_str *experimentname = arg_str0("w","experimentname","<str>","experiment name SQLite3 database file");
struct arg_str *experimentdescr = arg_str0("W","experimentdescr","<str>","experiment description SQLite3 database file");

struct arg_int *threads = arg_int0("T","threads","<int>",		"number of processing threads 0..128 (defaults to 0 which sets threads to number of CPU cores)");

struct arg_end *end = arg_end(200);

void *argtable[] = {help,version,FileLogLevel,LogFile,
	                pmode,subs100bp,end12subs,minpereadlen,minpeinsertsize,maxpeinsertsize,minscafflen,orientatepe,
					inpe1file,inpe2file,contigsfile,outfile,
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


	PMode = pmode->count ? pmode->ival[0] : (int)0;
	if(PMode < 0 || PMode > 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Processing mode '-m%d' must be 0",PMode);
		return(1);
		}

	OrientatePE = orientatepe->count ? orientatepe->ival[0] : 0;
	if(OrientatePE < 0 || OrientatePE > 3)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Expected paired end orientation '-M%d' to be in range 0..3",OrientatePE);
		return(1);
		}

	Subs100bp = subs100bp->count ? subs100bp->ival[0] : 0;
	if(Subs100bp < 0 || Subs100bp > cMaxSubs100bp)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Max allowed induced substitutions '-s%d' per 100bp overlapping must be in range 0..%d",Subs100bp,cMaxSubs100bp);
		return(1);
		}

	End12Subs = end12subs->count ? end12subs->ival[0] : 0;
	if(End12Subs < 0 || End12Subs > 6)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Max allowed induced substitutions in overlap 12bp ends '-E%d' must be in range 0..6",End12Subs);
		return(1);
		}

	MinPEReadlen = minpereadlen->count ? minpereadlen->ival[0] : cMinDfltSeqLenToAssemb;
	if (MinPEReadlen < 30 || MinPEReadlen > 5000)
	{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: Minimum accpted PE read length '-L%d' must be in range 30..5000", MinPEReadlen);
		return(1);
	}


	MinPEInsertSize = minpeinsertsize->count ? minpeinsertsize->ival[0] : 110;
	if(MinPEInsertSize < 100 || MinPEInsertSize > 50000)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: minimum allowed PE insert size '-p%d' must be in range 100..50000",MinPEInsertSize);
		return(1);
		}

	MaxPEInsertSize = maxpeinsertsize->count ? maxpeinsertsize->ival[0] : 1000;
	if(MaxPEInsertSize < MinPEInsertSize || MinPEInsertSize > 50000)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: maximum allowed PE insert size '-P%d' must be in range %d..50000",MaxPEInsertSize,MinPEInsertSize);
		return(1);
		}


	MinScaffoldedSeqLen = minscafflen->count ? minscafflen->ival[0] : cDfltMinScaffLen;
	if(MinScaffoldedSeqLen < cMinScaffLen || MinScaffoldedSeqLen > cMaxScaffLen)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Minimum reported scaffolded sequence length '-P%d' must be in range %d..%d",MinScaffoldedSeqLen,cMinScaffLen,cMaxScaffLen);
		return(1);
		}

	if(inpe1file->count)
		{
		strcpy(szPE1File,inpe1file->filename[0]);
		CUtility::TrimQuotedWhitespcExtd(szPE1File);
		}
	else
		szPE1File[0] = '\0';

	if(inpe2file->count)
		{
		strcpy(szPE2File,inpe2file->filename[0]);
		CUtility::TrimQuotedWhitespcExtd(szPE2File);
		}
	else
		szPE2File[0] = '\0';

	if(contigsfile->count)
		{
		strcpy(szContigsFile,contigsfile->filename[0]);
		CUtility::TrimQuotedWhitespcExtd(szContigsFile);
		}
	else
		szContigsFile[0] = '\0';

	// check that if one PE specified then the other PE is also specified
	if((szPE1File[0] == '\0' && szPE2File[0] != '\0') ||
		(szPE1File[0] != '\0' && szPE2File[0] == '\0'))
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Both 5' and 3' paired ends must be specified with '-a' and '-A'");
		return(1);
		}
	// both or either PE SE must be specified
	if(szPE1File[0] == '\0' && szContigsFile[0] == '\0')
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: No PE or SE/Contig input files specified");
		return(1);
		}

	strncpy(szOutFile,outfile->filename[0],_MAX_PATH);
	szOutFile[_MAX_PATH-1] = '\0';

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
	switch(PMode) {
		case 0:
			pszDescr = "Output scaffold multifasta with edge report";
			break;
		}

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Processing mode is : '%s'",pszDescr);

	gDiagnostics.DiagOutMsgOnly(eDLInfo, "Minimum PE read length: %d", MinPEReadlen);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Allow max induced substitutions per 100bp overlapping sequence fragments: %d",Subs100bp);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Allow max induced substitutions end 12bp of overlaps: %d",End12Subs);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Minimum PE insert size: %d",MinPEInsertSize);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Maximum PE insert size: %d",MaxPEInsertSize);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Minimum reported scaffolded sequence length: %d",MinScaffoldedSeqLen);

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
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Input 5' paired end reads or fragments file : '%s'",szPE1File);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Input 3' paired end reads or fragments file : '%s'",szPE2File);
		}
	if(szContigsFile[0] != '\0')
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Input contigs file: '%s'",szContigsFile);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Output scaffolded contigs multifasta file : '%s'",szOutFile);
	if(szExperimentName[0] != '\0')
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"This processing reference: %s",szExperimentName);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Number of threads : %d",NumThreads);

	if(gExperimentID > 0)
		{
		int ParamID;
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szLogFile),"log",szLogFile);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,sizeof(PMode),"mode",&PMode);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,sizeof(Subs100bp),"subskbp",&Subs100bp);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,sizeof(End12Subs),"maxendsubs",&End12Subs);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID, ePTInt32, sizeof(MinPEReadlen), "minpereadlen", &MinPEReadlen);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,sizeof(MinPEInsertSize),"minpeinsert",&MinPEInsertSize);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,sizeof(MaxPEInsertSize),"maxpeinsert",&MaxPEInsertSize);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,sizeof(MinScaffoldedSeqLen),"minscafflen",&MinScaffoldedSeqLen);


		if(szPE1File[0] != '\0')
			{
			ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,sizeof(int),"orientatepe",&OrientatePE);
			ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szPE1File),"inpe1",szPE1File);
			ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szPE2File),"inpe2",szPE2File);
			}
		if(szContigsFile[0] != '\0')
			ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szContigsFile),"contigsfile",szContigsFile);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szOutFile),"out",szOutFile);

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
	Rslt = ProcessScaffoldedContigs(PMode,NumThreads,bAffinity,Subs100bp,End12Subs, MinPEReadlen,MinPEInsertSize,MaxPEInsertSize,MinScaffoldedSeqLen,OrientatePE,szPE1File,szPE2File,szContigsFile,szOutFile);
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
GenScaffoldedContigs(int PMode,				// processing mode: 0 - output scaffold multifasta with edge report
		int NumThreads,						// number of worker threads to use
		bool bAffinity,						// thread to core affinity
		int Subs100bp,						// allow this many induced substitutions per 100bp overlapping sequence fragments (defaults to 1, range 0..5)
		int MaxEnd12Subs,					// allow at the initial 12bp (2 x hexamer primer artefacts) of the 5' or 3' of overlap to have this many base mismatches in addition to the overall allowed MaxSubs 
        int MinPEReadlen,					// only accept PE read lengths of at least this many bp
		int MinPEInsertSize,				// PE sequences are of this minimum insert size
		int MaxPEInsertSize,				// PE sequences are of this maximum insert size
		int MinScaffoldedSeqLen,			// reported scaffolded sequences are at least this length
		int OrientatePE,					// PE end orientations 0: sense/antisense, 1: sense/sense, 2: antisense/sense, 3: antisense/antisense 
		char *pszPE1File,					// input PE1 sequences file
		char *pszPE2File,					// input PE2 sequences file
		char *pszContigsFile,				// input SE contigs file
		char *pszScaffoldedFile)			// where to write scaffolded contigs
{
int Rslt;

CScaffolder *pScaffolder;
pScaffolder = new CScaffolder();
pScaffolder->Reset(false);
pScaffolder->SetPMode(PMode);
pScaffolder->SetNumThreads(NumThreads);
pScaffolder->SetSfxSparsity(eSSparsity15);
Rslt = pScaffolder->GenScaffoldedContigs(PMode,Subs100bp,MaxEnd12Subs, MinPEReadlen,MinPEInsertSize,MaxPEInsertSize,MinScaffoldedSeqLen,OrientatePE,pszPE1File,pszPE2File,pszContigsFile,pszScaffoldedFile);
delete pScaffolder;
return(Rslt < eBSFSuccess ? Rslt : 0);
}

// main processing thread startup - has increased stack size because of deep recursion when graph processing
#ifdef _WIN32
unsigned __stdcall ThreadedScaffolds(void * pThreadPars)
#else
void *ThreadedScaffolds(void * pThreadPars)
#endif
{
int Rslt = 0;
tsScaffoldPars *pPars = (tsScaffoldPars *)pThreadPars; // makes it easier not having to deal with casts!
Rslt = GenScaffoldedContigs(pPars->PMode,pPars->NumThreads,pPars->bAffinity,pPars->Subs100bp,pPars->MaxEnd12Subs, pPars->MinPEReadlen,
						pPars->MinPEInsertSize,pPars->MaxPEInsertSize,pPars->MinScaffoldedSeqLen,pPars->OrientatePE,
						pPars->pszPE1File,pPars->pszPE2File,pPars->pszContigsFile,pPars->pszScaffoldedFile);
pPars->Rslt = Rslt;
#ifdef _WIN32
_endthreadex(0);
return(eBSFSuccess);
#else
pthread_exit(NULL);
#endif
}

int
ProcessScaffoldedContigs(int PMode,				// processing mode: 0 - output scaffold multifasta with edge report
		int NumThreads,						// number of worker threads to use
		bool bAffinity,						// thread to core affinity
		int Subs100bp,						// allow this many induced substitutions per 100bp overlapping sequence fragments (defaults to 1, range 0..5)
		int MaxEnd12Subs,					// allow at the initial 12bp (2 x hexamer primer artefacts) of the 5' or 3' of overlap to have this many base mismatches in addition to the overall allowed MaxSubs 
        int MinPEReadlen,					// only accept PE read lengths of at least this many bp
		int MinPEInsertSize,				// PE sequences are of this minimum insert size
		int MaxPEInsertSize,				// PE sequences are of this maximum insert size
		int MinScaffoldedSeqLen,			// reported scaffolded sequences are at least this length
		int OrientatePE,					// PE end orientations 0: sense/antisense, 1: sense/sense, 2: antisense/sense, 3: antisense/antisense 
		char *pszPE1File,					// input PE1 sequences file
		char *pszPE2File,					// input PE2 sequences file
		char *pszContigsFile,				// input SE contigs file
		char *pszScaffoldedFile)			// where to write scaffolded contigs
{
tsScaffoldPars ScaffoldPars;

ScaffoldPars.PMode = PMode;				
ScaffoldPars.NumThreads = NumThreads;					
ScaffoldPars.bAffinity = bAffinity;						
ScaffoldPars.Subs100bp = Subs100bp;						
ScaffoldPars.MaxEnd12Subs = MaxEnd12Subs;	
ScaffoldPars.MinPEReadlen = MinPEReadlen;
ScaffoldPars.MinPEInsertSize = MinPEInsertSize;	
ScaffoldPars.MaxPEInsertSize = MaxPEInsertSize;
ScaffoldPars.MinScaffoldedSeqLen = MinScaffoldedSeqLen;	
ScaffoldPars.OrientatePE = OrientatePE;			
ScaffoldPars.pszPE1File = pszPE1File;					
ScaffoldPars.pszPE2File = pszPE2File;					
ScaffoldPars.pszContigsFile = pszContigsFile;			
ScaffoldPars.pszScaffoldedFile = pszScaffoldedFile;		

#ifdef _WIN32
ScaffoldPars.threadHandle = (HANDLE)_beginthreadex(NULL,cMainThreadStackSize, ThreadedScaffolds, &ScaffoldPars, 0, &ScaffoldPars.threadID);
Sleep(3000);
while (WAIT_TIMEOUT == WaitForSingleObject(ScaffoldPars.threadHandle, 60000)) 
	{  // one day could put a reporting progress function in here ....
	};
CloseHandle(ScaffoldPars.threadHandle);
#else
	// could recurse very deeply, need to increase the default stack of just 2MB
size_t defaultStackSize;
struct timespec ts;
int JoinRslt;
pthread_attr_t threadattr; 
pthread_attr_init(&threadattr);
pthread_attr_getstacksize(&threadattr, &defaultStackSize);
if(defaultStackSize < cMainThreadStackSize)
	pthread_attr_setstacksize(&threadattr, cMainThreadStackSize);
ScaffoldPars.threadRslt = pthread_create(&ScaffoldPars.threadID, &threadattr, ThreadedScaffolds, &ScaffoldPars);
pthread_attr_destroy(&threadattr);
sleep(3);
clock_gettime(CLOCK_REALTIME, &ts);
ts.tv_sec += 60;
while ((JoinRslt = pthread_timedjoin_np(ScaffoldPars.threadID, NULL, &ts)) != 0)
	{ // one day could put a reporting progress function in here ....
	ts.tv_sec += 60;
	}
#endif

return(ScaffoldPars.Rslt);
}


// relies on base classes constructors
CScaffolder::CScaffolder(void)
{
ScaffolderInit();
}

// relies on base classes destructors
CScaffolder::~CScaffolder(void)
{
ScaffolderReset(false);
}


int
CScaffolder::ScaffolderInit(void)
{
m_pSeqEdges = NULL;
m_ppToSeqEdges = NULL;
m_pSeqVertices = NULL;
m_pszScaffoldBuff = NULL;
m_hScaffoldFasta = -1;
CKangadna::Reset(false);
CdeNovoAssemb::Reset(false);
ScaffolderReset(false);
return(eBSFSuccess);
}

int 
CScaffolder::ScaffolderReset(bool bSync)
{
if(m_hScaffoldFasta != -1)
	{
	if(bSync)
		{
#ifdef _WIN32
		_commit(m_hScaffoldFasta);
#else
		fsync(m_hScaffoldFasta);
#endif
		}
	close(m_hScaffoldFasta);
	m_hScaffoldFasta = -1;
	}

if(m_pSeqEdges != NULL)
	{
#ifdef _WIN32
	free(m_pSeqEdges);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pSeqEdges != MAP_FAILED)
		munmap(m_pSeqEdges,m_AllocMemSeqEdges);
#endif	
	m_pSeqEdges = NULL;
	}

if(m_ppToSeqEdges != NULL)
	{
#ifdef _WIN32
	free(m_ppToSeqEdges);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_ppToSeqEdges != MAP_FAILED)
		munmap(m_ppToSeqEdges,m_AllocMemToSeqEdges);
#endif	
	m_ppToSeqEdges = NULL;
	}

if(m_pSeqVertices != NULL)
	{
#ifdef _WIN32
	free(m_pSeqVertices);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pSeqVertices != MAP_FAILED)
		munmap(m_pSeqVertices,m_AllocMemSeqVertices);
#endif	
	m_pSeqVertices = NULL;
	}

if(m_pszScaffoldBuff != NULL)
	{
	delete m_pszScaffoldBuff;
	m_pszScaffoldBuff = NULL;
	}

CKangadna::Reset(bSync);
CdeNovoAssemb::Reset(bSync);
m_AllocMemSeqEdges = 0;
m_NumSeqEdges = 0;
m_AllocMemToSeqEdges = 0;
m_AllocMemSeqVertices = 0;
m_NumSeqVertices = 0;
m_LowestEdgeFromSeqID = 0;			 
m_HighestEdgeFromSeqID = 0;				
m_LowestEdgeToSeqID = 0;					 
m_HighestEdgeToSeqID = 0;	
m_NumPEOverlapping = 0;
m_NumSEOverlapping = 0;
m_NumPEOverlapped = 0;
m_NumSEOverlapped = 0;

m_szScaffoldSetsFile[0] = '\0';
m_ScaffoldBuffLen = 0;
m_AllocScaffoldBuff = 0;
m_NumScaffoldSets = 0;
m_ScaffoldSetScore = 0;
m_ScaffoldSetGaps = 0;
return(eBSFSuccess);
}

teBSFrsltCodes
CScaffolder::SetNumThreads(int maxThreads)
{
if(maxThreads < 0 || maxThreads > cMaxWorkerThreads)
		return(eBSFerrParams);
CdeNovoAssemb::SetNumThreads(maxThreads);
m_NumThreads = maxThreads;
m_MTqsort.SetMaxThreads(maxThreads);
return(eBSFSuccess);
}

int
CScaffolder::ScaffoldAssemble(bool bSenseStrandOnly,			// sequences from sense strand specific
		bool bSingleEnded,				// treat all sequences as being single ended even if loaded as paired ends
		int OrientatePE,					// PE end orientations 0: sense/antisense, 1: sense/sense, 2: antisense/sense, 3: antisense/antisense 
		int NumThreads,					// number of worker threads to use
		bool bAffinity,					// thread to core affinity
		int MinPEReadLen,				// PE reads must be at least this many bp long
		char *pszPE1File,				// input PE1 sequences file
		char *pszPE2File,				// input PE2 sequences file
		char *pszSeedContigsFile)		// input SE contigs to be scaffolded file
{
int Rslt;
int EstNumCtgs;
int SeqWrdBytes;
INT64 CumulativeMemory;
UINT32 CumulativeSequences;

Reset(false);
SetPMode(0);
SetNumThreads(m_NumThreads);
SetSfxSparsity(eSSparsity15);
SeqWrdBytes = GetSeqWrdBytes();

CumulativeMemory = 0;
CumulativeSequences = 0;

if(pszSeedContigsFile != NULL && pszSeedContigsFile[0] != '\0')
	{
	if((Rslt = EstMemReq(pszSeedContigsFile)) != eBSFSuccess)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to estimate memory required for '%s",pszSeedContigsFile);
		return(eBSFerrOpnFile);	// treat as though unable to open file
		}
	EstMemReq(0,0,0,0,0);
	EstNumCtgs = GetEstSeqsToProc();
	}
else
	EstNumCtgs = 0;


if(pszPE1File != NULL && pszPE1File[0] != '\0')
	{
	if((Rslt =EstMemReq(pszPE1File)) != eBSFSuccess)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to estimate memory required for '%s",pszPE1File);
		return(eBSFerrOpnFile);	// treat as though unable to open file
		}
	if((Rslt = EstMemReq(pszPE2File)) != eBSFSuccess)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to estimate memory required for '%s",pszPE1File);
		return(eBSFerrOpnFile);	// treat as though unable to open file
		}
	}

CumulativeMemory = EstMemReq(0,0,0,0,0);
CumulativeSequences = GetEstSeqsToProc();
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
if((Rslt = AllocSeqs2AssembMem((CumulativeMemory * 110)/100))!= eBSFSuccess)	// add 10% to reduce risk of later having to realloc....
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to continue");
	Reset(false);
	return(Rslt);
	}

	// contigs may actually be scaffolds so need to alloc memory for holding loci of indeterminate blocks - 'N's were used to separate contigs in same scaffold in 
    // previous iterations. Estimate number of blocks likely to be required, if fewer contigs then guestimate these are fewer because they are the result of earlier
    // scaffolding iterations so likely to have more contained 'N's
if(EstNumCtgs < 100000)
	EstNumCtgs *= 25;
else
	{
	if(EstNumCtgs < 1000000)
		EstNumCtgs *= 10;
	else
		{
		if(EstNumCtgs < 10000000)
			EstNumCtgs *= 5;
		else
			{
			if(EstNumCtgs < 50000000)
				EstNumCtgs *= 2;
			}
		}
	}	
 
if((Rslt = AllocBlockNsLoci(EstNumCtgs))!= eBSFSuccess)	
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to continue");
	Reset(false);
	return(Rslt);
	}

if((Rslt = LoadSeqsOnly(bSenseStrandOnly,bSingleEnded,OrientatePE, MinPEReadLen,pszPE1File, pszPE2File, pszSeedContigsFile)) < eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"AssembReads: Unable to continue");
	Reset(false);
	return(Rslt);
	}

return(eBSFSuccess);
}

teBSFrsltCodes
CScaffolder::GenScaffoldedContigs(int PMode,	//  processing mode: 0 - output single multifasta with edge report
					int Subs100bp,				// allow this many induced substitutions per 100bp overlapping sequence fragments (defaults to 0, range 0..5)
					int MaxEnd12Subs,			// allow at the initial 12bp (2 x hexamer primer artefacts) of the 5' or 3' of overlap to have this many base mismatches in addition to the overall allowed MaxSubs 
					int MinPEReadlen,			// only accept PE read lengths of at least this many bp
					int MinPEInsertSize,		// PE sequences are of this minimum insert size
					int MaxPEInsertSize,		// PE sequences are of this maximum insert size
					int MinScaffoldedSeqLen,	// reported scaffolded sequences are at least this length
					int OrientatePE,			// PE end orientations 0: sense/antisense, 1: sense/sense, 2: antisense/sense, 3: antisense/antisense 
					char *pszPE1File,			// input PE1 sequences file
					char *pszPE2File,			// input PE2 sequences file
					char *pszContigsFile,		// input SE contigs file
					char *pszScaffoldedFile)	// where to write scaffolded contigs
{
int Rslt;
int Subs1Kbp;
UINT32 CurTotNumPEs;	// total number of paired ends
int CurPE1MinLen;		// returned PE1 min length
int CurPE1MeanLen;		// returned PE1 mean length
int CurPE1MaxLen;		// returned PE1 max length
int CurPE2MinLen;		// returned PE2 min length
int CurPE2MeanLen;		// returned PE2 mean length
int CurPE2MaxLen;		// returned PE2 max length
UINT32 CurTotNumSEs;	// total number of single ends
int CurSEMinLen;		// returned SE min length
int CurSEMeanLen;		// returned SE mean length
int CurSEMaxLen;		// returned SE max length

UINT32 RemovedPEs;

m_MinReqPESepDist = MinPEInsertSize;
m_MaxReqPESepDist = MaxPEInsertSize;
m_MinScaffoldedSeqLen = MinScaffoldedSeqLen;
if(Subs100bp < 0 || Subs100bp > cMaxSubs100bp ||
   pszScaffoldedFile == NULL || pszScaffoldedFile[0] == '\0')
   return(eBSFerrParams);
Subs1Kbp = Subs100bp * 10;

	// load scaffolding PEs plus contigs to be scaffolded
if((Rslt = ScaffoldAssemble(false,false,OrientatePE,m_NumThreads,m_bAffinity, MinPEReadlen,pszPE1File,pszPE2File,pszContigsFile)) < eBSFSuccess)
	return((teBSFrsltCodes)Rslt);

	// get min, mean, max lengths for PE1, PE2 and SEs and report on these
GetSeqLenDist(&CurTotNumPEs,
				&CurPE1MinLen,&CurPE1MeanLen,&CurPE1MaxLen,
				&CurPE2MinLen,&CurPE2MeanLen,&CurPE2MaxLen,
				&CurTotNumSEs,
				&CurSEMinLen,&CurSEMeanLen,&CurSEMaxLen);

if(CurTotNumSEs < 2)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"GenScaffoldedContigs: Need at least 2 SE sequences to scaffold, loaded %u",CurTotNumSEs);
	return(eBSFerrNoEntries);
	}

if(CurTotNumPEs < 1)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"GenScaffoldedContigs: Need at least 1 PE to scaffold, loaded 0");
	return(eBSFerrNoEntries);
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"GenScaffoldedContigs: Allowing %d substitutions per overlapping Kbp, %d subs in overlap 12bp ends, processing total of %u sequences...",Subs1Kbp,MaxEnd12Subs,(2 * CurTotNumPEs) + CurTotNumSEs);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"GenScaffoldedContigs: Num PEs: %u PE1 lens - min: %u mean: %u max: %u  ",CurTotNumPEs,CurPE1MinLen,CurPE1MeanLen,CurPE1MaxLen);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"GenScaffoldedContigs: Num PEs: %u PE2 lens - min: %u mean: %u max: %u  ",CurTotNumPEs,CurPE2MinLen,CurPE2MeanLen,CurPE2MaxLen);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"GenScaffoldedContigs: Num SEs: %u  SE lens - min: %u mean: %u max: %u  ",CurTotNumSEs,CurSEMinLen,CurSEMeanLen,CurSEMaxLen);

RemovedPEs = RemoveContainedSeqs(Subs1Kbp,MaxEnd12Subs,MinPEInsertSize,MaxPEInsertSize,false);
if(RemovedPEs == CurTotNumPEs)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"GenScaffoldedContigs: Nothing to do, need at least 1 PE to scaffold, all PE sequences were contained within SE sequences");
	return(eBSFerrNoEntries);
	}

		// only flags of interest to retain are those for paired end or single end sequences
m_NumSeqEdges = 0;
UpdateAllSeqHeaderFlags(0,~(cFlgSeqPE2 | cFlgSeqPE | cFlgNonOverlap),false);

		// generate array of sequence starts plus array of flags from sequence headers
if((Rslt=GenSeqStarts(true,false)) < eBSFSuccess)
	return((teBSFrsltCodes)Rslt);

		// index all sequences but not the PEs
if((Rslt=GenRdsSfx(0,2,true)) < eBSFSuccess)
	return((teBSFrsltCodes)Rslt);

		// generate edges for sense PE probe sequences against sense SE targets
gDiagnostics.DiagOut(eDLInfo,gszProcName,"GenScaffoldedContigs: Processing for sense PE individual ends onto sense SE (vertices) as edges ...");
GenSeqEdges(0,Subs1Kbp,MaxEnd12Subs,MinPEInsertSize,MaxPEInsertSize);

		// generate graph for antisense PE probe sequences against sense SE targets
gDiagnostics.DiagOut(eDLInfo,gszProcName,"GenScaffoldedContigs: Processing for antisense PE individual end onto sense SE (vertices) as edges ...");
GenSeqEdges(1,Subs1Kbp,MaxEnd12Subs,MinPEInsertSize,MaxPEInsertSize);

	// now generate a graph with SE and PEs as the vertices and overlaps between same as the edges
gDiagnostics.DiagOut(eDLInfo,gszProcName,"GenScaffoldedContigs: Generating SE vertices with connecting PE edges scaffold graph");

if((Rslt=GenerateScaffoldGraph())!=eBSFSuccess)
	{
	Reset(false);
	return((teBSFrsltCodes)Rslt);
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"GenScaffoldedContigs: Generated scaffold graph has %u vertices with %u edges",m_NumSeqVertices,m_NumSeqEdges);

Rslt = ReportScaffoldSets(pszScaffoldedFile);

return(Rslt < 0 ? (teBSFrsltCodes)Rslt: eBSFSuccess);
}

static int m_RecurseDepth = 0;
static int m_DeepestRecurseDepth = 0;


int
CScaffolder::RemoveContainedSeqs(int Subs1Kbp, // allow this many induced substitutions per Kbp overlapping sequence fragments (defaults to 10, range 0..50)
					int MaxEnd12Subs,			 // allow at the initial 12bp (2 x hexamer primer artefacts) of the 5' or 3' of overlap to have this many base mismatches in addition to the overall allowed MaxSubs 
					int MinPEInsertSize,		// PE sequences are of this minimum insert size
					int MaxPEInsertSize,		// PE sequences are of this maximum insert size
					bool bPEOnly)				 // if false then identify both SE and PE (both ends) sequences which are fully contained, if true then identify contained PEs only 
{
int Rslt;
int TotContained;

		// only flags of interest to retain are those for paired end or single end sequences
m_NumSeqEdges = 0;

UpdateAllSeqHeaderFlags(0,~(cFlgSeqPE2 | cFlgSeqPE | cFlgNonOverlap),false);

				// generate array of sequence starts plus array of flags from sequence headers
if((Rslt=GenSeqStarts(true,false)) < eBSFSuccess)
	return((teBSFrsltCodes)Rslt);

	// index all sequences
if((Rslt=GenRdsSfx(0,2,true)) < eBSFSuccess)
	return((teBSFrsltCodes)Rslt);

			// identify sense probe sequences against sense targets
gDiagnostics.DiagOut(eDLInfo,gszProcName,"RemoveContainedSeqs: Identifying for removal sense PE paired ends contained within same SE ...");
if((Rslt = MarkContainedSeqs(0,bPEOnly,Subs1Kbp,MaxEnd12Subs,MinPEInsertSize,MaxPEInsertSize)) < 0)
	return(Rslt);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"RemoveContainedSeqs: Identified %u sense PE paired ends contained within same SE ...", Rslt);
TotContained = Rslt;

	// identify antisense probe sequences against sense targets
gDiagnostics.DiagOut(eDLInfo,gszProcName,"RemoveContainedSeqs: Identifying for removal antisense PE paired ends contained within same SE ...");
if((Rslt = MarkContainedSeqs(1,bPEOnly,Subs1Kbp,MaxEnd12Subs,MinPEInsertSize,MaxPEInsertSize)) < 0)
	return(Rslt);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"RemoveContainedSeqs: Identified %u antisense PE paired ends contained within same SE ...", Rslt);
TotContained += Rslt;
gDiagnostics.DiagOut(eDLInfo,gszProcName,"RemoveContainedSeqs: Removing %u PE paired ends identified as contained within same SE ...", TotContained);

if(TotContained)
	RemoveMarkedSeqs(cFlgContainRemove);


UINT32 CurTotNumPEs;
GetSeqLenDist(&CurTotNumPEs);	// also updates min, mean, max lengths for PE1, PE2 and SEs counts


gDiagnostics.DiagOut(eDLInfo,gszProcName,"RemoveContainedSeqs: Removed %u PE paired ends", TotContained);

return(TotContained);
}





int							// all vertices marked as requiring revcpl to make sequences sense consistent are marked as not to be processed for overlaps
CScaffolder::IdentifyNonSenseSeqs(void)
{
int TotNumNonSense;
UINT32 CurSeqVertexIdx;
tsSeqVertex *pCurSeqVertex;

if(m_pSeqVertices == NULL || m_NumSeqVertices < 1)	
	return(0);	


TotNumNonSense = 0;
pCurSeqVertex = m_pSeqVertices;
for(CurSeqVertexIdx = 0; CurSeqVertexIdx < m_NumSeqVertices; CurSeqVertexIdx++,pCurSeqVertex++)
	{
	if(!pCurSeqVertex->FlgReqRevCpl)		
		continue;
	UpdateSeqHeaderFlags(pCurSeqVertex->SeqID,cFlgNonOverlap,0);
	TotNumNonSense += 1;
	}
return(TotNumNonSense);
}


teBSFrsltCodes
CScaffolder::GenerateScaffoldGraph(void)		// generates scaffold graph over all edges (PE end overlaps) + vertices (SE)
{
UINT64 ToIdx;
UINT64 FromIdx;
tsSeqEdge *pToEdge;
tsSeqEdge *pFromEdge;

tSeqID CurToSeqID;
tSeqID SeqID;

UINT32 SeqFlgs;
UINT32 SeqLen;

tsSeqVertex SEVertex;
tVertID SeqVertID;


int NumRemovePEs;
int NumRemoveEdges;
bool bRemovePE;
tSeqID PE1SeqID;
UINT32 PE1Flags;
tsSeqEdge *pPE1Edge;
tEdgeID PE1EdgeID;
int NumPE1Edges;

tSeqID PE2SeqID;
tsSeqEdge *pPE2Edge;
tEdgeID PE2EdgeID;
int NumPE2Edges;

m_LowestEdgeFromSeqID = 0;			 
m_HighestEdgeFromSeqID = 0;				
m_LowestEdgeToSeqID = 0;					 
m_HighestEdgeToSeqID = 0;
m_NumPEOverlapping = 0;
m_NumSEOverlapping = 0;
m_NumPEOverlapped = 0;
m_NumSEOverlapped = 0;
m_NumSeqVertices = 0;

	// sort edges by FromSeqID ascending, OverlapLen descending
if(m_NumSeqEdges > 1)
	m_MTqsort.qsort(m_pSeqEdges,m_NumSeqEdges,sizeof(tsSeqEdge),SortSeqEdgesFromSeqID);

pToEdge = m_pSeqEdges;
m_LowestEdgeFromSeqID = 0;
m_HighestEdgeFromSeqID = 0;

for(ToIdx = 1; ToIdx <= m_NumSeqEdges; ToIdx++,pToEdge++)
	{
	if(m_LowestEdgeFromSeqID == 0 || pToEdge->ProbeSeqID < m_LowestEdgeFromSeqID)
		m_LowestEdgeFromSeqID = pToEdge->ProbeSeqID;
	if(pToEdge->ProbeSeqID  > m_HighestEdgeFromSeqID)
		m_HighestEdgeFromSeqID = pToEdge->ProbeSeqID;
	pToEdge->SeqEdgeID = ToIdx;
	}

// remove edges which are unsupported by the partner PE; i.e both ends of a PE must be overlapping onto a SE and they must overlap onto different SE
// at same time mark all PE sequences to be removed as these are no longer required
NumRemovePEs = 0;
NumRemoveEdges = 0;
bRemovePE = false;
for(SeqID = 1; SeqID <= m_Sequences.NumSeqs2Assemb; SeqID++)
	{
	GetSeqHeader(SeqID,					// 32 bit sequence identifier
						NULL,			// returned 5 bit source file identifier in bits 4..0
						&PE1Flags,		// returned 16 bit sequence flags in bits 20..0
						&SeqLen,		// returned 30 bit sequence length
						false);			// set true if access to headers are required to be serialised
	if(!((PE1Flags & (cFlgSeqPE | cFlgSeqPE2)) == cFlgSeqPE)) // looking for PE1 sequences
		continue;

	PE1SeqID = SeqID++;					// located a PE1
	PE2SeqID = SeqID;					// PE2 will immediately follow the PE1

	// iterate over all PE1 and PE2 edges looking for PE1 and PE2 overlapping the same target SE sequence - both PE1 and PE2 will be marked for deletion
    // at same time look for PE1 or PE2 not overlapping any SE - if either does not overlap then both PE1 and PE2 will be marked for deletion
	bRemovePE = false;
	pPE1Edge = NULL;
	PE1EdgeID = 0;
	NumPE1Edges = 0;
	NumPE2Edges = 0;
	while(!bRemovePE && ((pPE1Edge = IterateEdgeFromSeqID(PE1SeqID,&PE1EdgeID)) != NULL))
		{
		NumPE1Edges += 1;			// at least one PE1 overlap onto a SE
		PE2EdgeID = 0;
		while((pPE2Edge = IterateEdgeFromSeqID(PE2SeqID,&PE2EdgeID)) != NULL)
			{
			NumPE2Edges += 1;		// at least one PE2 overlap onto a SE
			if(pPE1Edge->TargetSeqID == pPE2Edge->TargetSeqID)  // are both PE1 and PE2 overlaying same SE?
				{
				bRemovePE = true;
				break;
				} 
			}
		}

	if(NumPE1Edges == 0 || NumPE2Edges == 0 || bRemovePE) // if either end had no overlaps or both overlapping onto same SE (bRemove set) then mark all effected PE edges with FlgSlough and delete the PEs 
		{
		PE1EdgeID = 0;
		pPE1Edge = NULL;
		while((pPE1Edge = IterateEdgeFromSeqID(PE1SeqID,&PE1EdgeID)) != NULL)
			{
			pPE1Edge->FlgSlough = 1;
			NumRemoveEdges += 1;
			}
		PE2EdgeID = 0;
		pPE2Edge = NULL;
		while((pPE2Edge = IterateEdgeFromSeqID(PE2SeqID,&PE2EdgeID)) != NULL)
			{
			pPE2Edge->FlgSlough = 1;
			NumRemoveEdges += 1;
			}
		}

	NumRemovePEs += 1;
	UpdateSeqHeaderFlags(PE1SeqID,cFlgContainRemove,0,false);
	UpdateSeqHeaderFlags(PE2SeqID,cFlgContainRemove,0,false);
	}

// now actually remove the PEs and also remove all contained or orphan edges
 if(NumRemovePEs || NumRemoveEdges)
	{
	if(NumRemovePEs)		// count of any contained or PE ends not overlapping onto a SE 
		RemoveMarkedSeqs(cFlgContainRemove);
	if(NumRemoveEdges)      // count of all edges which are to deleted, these will have had FlgSlough set
		{
		m_LowestEdgeFromSeqID = 0;
		m_HighestEdgeFromSeqID = 0;
		pToEdge = m_pSeqEdges;
		pFromEdge = pToEdge;
		ToIdx = 0;
		for(FromIdx = 1; FromIdx <= m_NumSeqEdges; FromIdx++,pFromEdge++)
			{
			if(pFromEdge->FlgSlough)
				continue;
			if(FromIdx != ToIdx+1)
				*pToEdge = *pFromEdge;
			if(m_LowestEdgeFromSeqID == 0 || pToEdge->ProbeSeqID < m_LowestEdgeFromSeqID)
				m_LowestEdgeFromSeqID = pToEdge->ProbeSeqID;
			if(pToEdge->ProbeSeqID  > m_HighestEdgeFromSeqID)
				m_HighestEdgeFromSeqID = pToEdge->ProbeSeqID;
			pToEdge->SeqEdgeID = ++ToIdx;
			pToEdge += 1;
			}
		m_NumSeqEdges = ToIdx;
		}
	}

// non-relevant SEs and PEs have been removed 
// now create index over remaining edges sorted by ToSeqID
if(m_ppToSeqEdges != NULL)
	{
#ifdef _WIN32
	free(m_ppToSeqEdges);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_ppToSeqEdges != MAP_FAILED)
		munmap(m_ppToSeqEdges,m_AllocMemToSeqEdges);
#endif	
	m_ppToSeqEdges = NULL;
	m_AllocMemToSeqEdges = 0;
	}

if(m_NumSeqEdges)
	{
	m_AllocMemToSeqEdges = (size_t)sizeof(tsSeqEdge *) * m_NumSeqEdges;
#ifdef _WIN32
	m_ppToSeqEdges = (tsSeqEdge **)malloc((size_t)m_AllocMemToSeqEdges);	
	if(m_ppToSeqEdges == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"GenerateGraph: Unable to allocate memory of %llu bytes - %s",m_AllocMemToSeqEdges,strerror(errno));
		Reset(false);
		return(eBSFerrMem);
		}
#else
	if((m_ppToSeqEdges = (tsSeqEdge **)mmap(NULL,m_AllocMemToSeqEdges, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0)) == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"GenerateGraph: Unable to allocate memory of %llu bytes - %s",m_AllocMemToSeqEdges,strerror(errno));
		m_ppToSeqEdges = NULL;
		Reset(false);
		return(eBSFerrMem);
		}
#endif

	pToEdge = m_pSeqEdges;
	for(ToIdx = 0; ToIdx < m_NumSeqEdges; ToIdx++,pToEdge++)
		m_ppToSeqEdges[ToIdx] = pToEdge;
	m_MTqsort.qsort(m_ppToSeqEdges,m_NumSeqEdges,sizeof(tsSeqEdge *),SortSeqEdgesToSeqID);

	CurToSeqID = 0;
	for(ToIdx = 0; ToIdx < m_NumSeqEdges; ToIdx++)
		{
		pToEdge = m_ppToSeqEdges[ToIdx];
		if(CurToSeqID != pToEdge->TargetSeqID)
			{
			CurToSeqID = pToEdge->TargetSeqID;
			if(m_LowestEdgeToSeqID == 0 || m_LowestEdgeToSeqID > CurToSeqID)
				m_LowestEdgeToSeqID = CurToSeqID;
			if(m_HighestEdgeToSeqID < CurToSeqID)
				m_HighestEdgeToSeqID = CurToSeqID;
			}
		}
	}

// now create and initialse the vertices, one for each SE sequence
tSeqWrd4 *pCurSeq;
pCurSeq = NULL;
while((pCurSeq = IterSeqHeaders(pCurSeq,	// iterate to next sequence following this
			&SeqID,				// returned sequence identifier
			NULL,				// returned 8 bit source file identifier
			&SeqFlgs,			// returned 16 bit sequence flags
			&SeqLen,false))!=NULL)			// returned 30 bit sequence length
	{
	if(SeqFlgs & cFlgSeqPE)
		continue;
	memset(&SEVertex,0,sizeof(SEVertex));
	SEVertex.SeqID = SeqID;
	SEVertex.SeqLen = SeqLen;
	GetPredSuccSeq(&SEVertex);
	SeqVertID = AddSeqVertex(&SEVertex);
	}

// create index over vertices on the SeqID ascending
m_MTqsort.qsort(m_pSeqVertices,m_NumSeqVertices,sizeof(tsSeqVertex),SortVerticesSeqID);

// now check for inconsistent vertices predecessor and successor linkages
UINT32 VertIdx;

tSeqID TmpSeqID;
UINT16 TmpCnts;
UINT16 TmpScore;
tsSeqVertex *pVertex;
tsSeqVertex *pPredVertex;
tsSeqVertex *pCurVertex;
tsSeqVertex *pSuccVertex;

int SuccProblems = 0;
int PredProblems = 0;
int SuccOk = 0;
int PredOk = 0;


// attempt to make all scaffold set vertices the same relative sense
// 
pVertex = m_pSeqVertices;
for(VertIdx = 0; VertIdx < m_NumSeqVertices; VertIdx++,pVertex++)
	{
	if(pVertex->FlgSenseChkd || (pVertex->PredSeqID == 0 && pVertex->SuccSeqID == 0))	// check next vertex if this is an isolated vertex or has already been sense checked
		continue;

	// this vertex is treated as being the sense reference vertex with all other vertices having their linkages sense relative to this refererence
	// first follow all predecessor linkages adjusting linkages as may be required for those which are antisense to this reference
	if(pVertex->PredSeqID != 0)
		{
		pCurVertex = pVertex;
		do {
			if(pCurVertex->FlgSenseChkd)							// would only be set if there is a circular linkage
				{
				if(pCurVertex->SuccSeqID == pSuccVertex->SeqID)
					pCurVertex->SuccSeqID = 0;
				if(pCurVertex->PredSeqID == pSuccVertex->SeqID)
					pCurVertex->PredSeqID = 0;
				pSuccVertex->PredSeqID = 0;							// breaks the linkage
				break;
				}
			else
				pCurVertex->FlgSenseChkd = 1;						// mark as sense checked so can determine when iterating if a circular linkage 

			if(pCurVertex->PredSeqID)
				{
				pPredVertex = LocateVerticesSeqID(pCurVertex->PredSeqID);

				if(pPredVertex->PredSeqID == pCurVertex->SeqID)
					{
					TmpSeqID = pPredVertex->PredSeqID;					// exchange pred/succ linkage and mark to be RevCpl when writing sequence to multifasta file
					pPredVertex->PredSeqID = pPredVertex->SuccSeqID;
					pPredVertex->SuccSeqID = TmpSeqID;
					TmpCnts = pPredVertex->PredNumPEs;
					pPredVertex->PredNumPEs = pPredVertex->SuccNumPEs;
					pPredVertex->SuccNumPEs = TmpSeqID;
					TmpScore = pPredVertex->PredScore;
					pPredVertex->PredScore = pPredVertex->SuccScore;
					pPredVertex->SuccScore = TmpScore;
					pPredVertex->FlgReqRevCpl = 1;
					}

				if(pPredVertex->SuccSeqID != pCurVertex->SeqID ||  pCurVertex->PredSeqID != pPredVertex->SeqID)
					{
					pPredVertex->SuccSeqID = 0;
					pCurVertex->PredSeqID = 0;
					break;
					}
				

				}
			else
				pPredVertex = NULL;
			pSuccVertex = pCurVertex;									// copy the current vertex in case the predecessor shows up as circular
			}
		while((pCurVertex = pPredVertex)!=NULL);
		}

	if(pVertex->SuccSeqID != 0)								
		{
		pCurVertex = pVertex;
		pCurVertex->FlgSenseChkd = 0;								// don't want the check for circular linkage reporting yes because marked during predecessor processing!
		do {
			if(pCurVertex->FlgSenseChkd)							// would only be set if there is a circular linkage
				{
				if(pCurVertex->PredSeqID == pSuccVertex->SeqID)
					pCurVertex->PredSeqID = 0;
				if(pCurVertex->SuccSeqID == pSuccVertex->SeqID)
					pCurVertex->SuccSeqID = 0;
				pPredVertex->SuccSeqID = 0;							// breaks the linkage
				break;
				}
			else
				pCurVertex->FlgSenseChkd = 1;						// mark as sense checked so can determine when iterating if a circular linkage 

			if(pCurVertex->SuccSeqID)
				{
				pSuccVertex = LocateVerticesSeqID(pCurVertex->SuccSeqID);
				if(pSuccVertex->SuccSeqID == pCurVertex->SeqID)
					{
					TmpSeqID = pSuccVertex->PredSeqID;						// exchange pred/succ linkage and mark to be RevCpl when writing sequence to multifasta file
					pSuccVertex->PredSeqID = pSuccVertex->SuccSeqID;
					pSuccVertex->SuccSeqID = TmpSeqID;
					TmpCnts = pSuccVertex->PredNumPEs;
					pSuccVertex->PredNumPEs = pSuccVertex->SuccNumPEs;
					pSuccVertex->SuccNumPEs = TmpSeqID;
					TmpScore = pSuccVertex->PredScore;
					pSuccVertex->PredScore = pSuccVertex->SuccScore;
					pSuccVertex->SuccScore = TmpScore;
					pSuccVertex->FlgReqRevCpl = 1;
					}
				if(pSuccVertex->PredSeqID != pCurVertex->SeqID || pCurVertex->SuccSeqID != pSuccVertex->SeqID)
					{
					pSuccVertex->PredSeqID = 0;
					pCurVertex->SuccSeqID = 0;
					break;
					}

				}
			else
				pSuccVertex = NULL;
			pPredVertex = pCurVertex;	// copy the current vertex in case the successor shows up as circular	
			}
		while((pCurVertex = pSuccVertex)!=NULL);
		}
	}

// expecting that if a vertex has successor linkage then that successor will have predecessor linka
pVertex = m_pSeqVertices;
for(VertIdx = 0; VertIdx < m_NumSeqVertices; VertIdx++,pVertex++)
	{
	if(pVertex->PredSeqID == 0 && pVertex->SuccSeqID == 0)	// check next vertex if no linkage
		continue;
	
	if(pVertex->PredSeqID != 0)
		{
		pPredVertex = LocateVerticesSeqID(pVertex->PredSeqID);
		if(pPredVertex->SuccSeqID != pVertex->SeqID)		// if the predecessor does not have the current vertex SeqID as it's successor then there is a problem!
			{
			pVertex->PredSeqID = 0;
			pVertex->PredNumPEs = 0;
			pVertex->PredScore = 0;
			SuccProblems++;
			}
		else
			SuccOk++;
		}

	if(pVertex->SuccSeqID != 0)
		{
		pSuccVertex = LocateVerticesSeqID(pVertex->SuccSeqID);
		if(pSuccVertex->PredSeqID != pVertex->SeqID)		// if the successor does not have the current vertex SeqID as it's predecessor then there is a problem!
			{
			pVertex->SuccSeqID = 0;
			pVertex->SuccNumPEs = 0;
			pVertex->SuccScore = 0;
			PredProblems++;
			}
		else
			PredOk++;
		}
	}

// identify any circular references and break any detected
tSeqID CurSuccID;
pVertex = m_pSeqVertices;
for(VertIdx = 0; VertIdx < m_NumSeqVertices; VertIdx++,pVertex++)
	{
	if(((CurSuccID = pVertex->SuccSeqID) == 0) || pVertex->FlgCircChkd)	// check next vertex if no successor linkage and already checked for circular references
		continue;
	while((pSuccVertex = LocateVerticesSeqID(CurSuccID))!=NULL)
		{
		pSuccVertex->FlgCircChkd = 1;
		if(pSuccVertex->SuccSeqID == pVertex->SeqID)
			{
			pSuccVertex->SuccSeqID = 0;
			pSuccVertex->SuccNumPEs = 0;
			pSuccVertex->SuccScore = 0;
			pVertex->PredSeqID = 0;
			pVertex->PredNumPEs = 0;
			pVertex->PredScore = 0;
			break;
			}
		CurSuccID = pSuccVertex->SuccSeqID;
		}
	}

return(eBSFSuccess);
}



int																	// total number of reported scaffolds 
CScaffolder::ReportScaffoldSets(char *pszScaffoldFile)				// report to this multifasta pszScaffoldFile
{
int NumReported;
int LineLen;
int BaseIdx;
int NIdx;
char *pszAsciiBases;
UINT32 VertIdx;
UINT32 NumScaffCtgs;
tSeqID CurPredID;
tSeqID CurSuccID;
tsSeqVertex *pVertex;
tsSeqVertex *pPredVertex;
tsSeqVertex *pSuccVertex;
int ScaffLen;

strcpy(m_szScaffoldSetsFile,pszScaffoldFile);

// allocate buffering for buffering sequences as ascii multifasta
if(m_pszScaffoldBuff == NULL)
	{
	if((m_pszScaffoldBuff = new char [cAllocScaffoldBuffSize]) == NULL)
		return(eBSFerrMem);
	m_AllocScaffoldBuff = cAllocScaffoldBuffSize;
	}
m_ScaffoldBuffLen = 0;

// create/truncate scaffolded sequences file
#ifdef _WIN32
if((m_hScaffoldFasta = open(m_szScaffoldSetsFile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
if((m_hScaffoldFasta = open(m_szScaffoldSetsFile, O_RDWR | O_CREAT |O_TRUNC, S_IREAD | S_IWRITE))==-1)
#endif
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"ReportScaffoldGraph: Unable to create or truncate %s - %s",m_szScaffoldSetsFile,strerror(errno));
	Reset();
	return(eBSFerrCreateFile);
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"ReportScaffoldGraph: Output scaffolded sequence set multifasta file created/truncated: '%s'",m_szScaffoldSetsFile);


NumReported = 0;
m_AcceptedScaffolds = 0;
m_UnderlenScaffolds = 0;
pVertex = m_pSeqVertices;
for(VertIdx = 0; VertIdx < m_NumSeqVertices; VertIdx++,pVertex++)
	{
	if(pVertex->FlgPathed)
		continue;

	// locate starting contig for this scaffold
	pPredVertex = pVertex;
	while((CurPredID = pPredVertex->PredSeqID) != 0)
		pPredVertex = LocateVerticesSeqID(CurPredID);

	// pSuccVertex now pts to starting vertex
	// calc total length of this scaffold allowing for N-mer contig sequence separators
	ScaffLen = 0;
	NumScaffCtgs = 0;
	pSuccVertex = pPredVertex;
	do {
		ScaffLen += pSuccVertex->SeqLen;
		NumScaffCtgs += 1;
		CurSuccID = pSuccVertex->SuccSeqID;
		if(CurSuccID != 0)
			{
			ScaffLen += cDfltScaffoldGapNs;
			pSuccVertex = LocateVerticesSeqID(CurSuccID);
			}
		}
	while(CurSuccID != 0);

	if(m_MinScaffoldedSeqLen && ScaffLen < m_MinScaffoldedSeqLen)
		{
		m_UnderlenScaffolds += 1;
		continue;
		}

	// length is known, can now output the descriptor line
	m_ScaffoldBuffLen += sprintf(&m_pszScaffoldBuff[m_ScaffoldBuffLen],">Scaff%d %d|%d\n",NumReported+1,NumScaffCtgs,ScaffLen);

	// can now output the scaffold sequence
	LineLen = 0;
	pSuccVertex = pPredVertex;
	do {
		pSuccVertex->FlgPathed = 1;

		// will write sequence out here ....
		pszAsciiBases = AsciifySequence(pSuccVertex->SeqID,			// sequence identifier
				pSuccVertex->SeqLen,			// asciify at most MaxSeqLen bases (limited to at most cMaxDiagAsciiSeqLen)
				NumScaffCtgs > 1 && pSuccVertex->FlgReqRevCpl == 1 ? true : false);			// true if sequence to be reverse complemented before asciifying

		for(BaseIdx = 0; BaseIdx < (int)pSuccVertex->SeqLen; BaseIdx++,pszAsciiBases++)
			{
			m_pszScaffoldBuff[m_ScaffoldBuffLen++] = *pszAsciiBases;
			LineLen += 1;
			if(LineLen > 75)
				{
				m_ScaffoldBuffLen += sprintf(&m_pszScaffoldBuff[m_ScaffoldBuffLen],"\n");
				LineLen = 0;
				}
			if((m_ScaffoldBuffLen + 2000) > m_AllocScaffoldBuff)
				{
				CUtility::SafeWrite(m_hScaffoldFasta,m_pszScaffoldBuff,m_ScaffoldBuffLen);
				m_ScaffoldBuffLen = 0;
				}
			}

		CurSuccID = pSuccVertex->SuccSeqID;
		if(CurSuccID != 0)
			{
			// write out separator N-mer here
			for(NIdx = 0; NIdx < cDfltScaffoldGapNs; NIdx++)
				{
				m_pszScaffoldBuff[m_ScaffoldBuffLen++] = 'N';
				LineLen += 1;
				if(LineLen > 75)
					{
					m_ScaffoldBuffLen += sprintf(&m_pszScaffoldBuff[m_ScaffoldBuffLen],"\n");
					LineLen = 0;
					}
				}

			pSuccVertex = LocateVerticesSeqID(CurSuccID);
			}
		}
	while(CurSuccID != 0);

	if(LineLen > 0)
		m_ScaffoldBuffLen += sprintf(&m_pszScaffoldBuff[m_ScaffoldBuffLen],"\n");
	
	if((m_ScaffoldBuffLen + 2000) > m_AllocScaffoldBuff)
		{
		CUtility::SafeWrite(m_hScaffoldFasta,m_pszScaffoldBuff,m_ScaffoldBuffLen);
		m_ScaffoldBuffLen = 0;
		}
	m_AcceptedScaffolds += 1;
	NumReported += 1;
	}

if(m_ScaffoldBuffLen)
	CUtility::SafeWrite(m_hScaffoldFasta,m_pszScaffoldBuff,m_ScaffoldBuffLen);

#ifdef _WIN32
_commit(m_hScaffoldFasta);
#else
fsync(m_hScaffoldFasta);
#endif

close(m_hScaffoldFasta);
m_hScaffoldFasta = -1;

delete m_pszScaffoldBuff;
m_pszScaffoldBuff = NULL;
m_ScaffoldBuffLen = 0;
m_AllocScaffoldBuff = 0;

gDiagnostics.DiagOut(eDLInfo,gszProcName,"ReportScaffoldGraph: Wrote %u scaffold sets, %u scaffold sets rejected because underlength %d",m_AcceptedScaffolds,m_UnderlenScaffolds,m_MinScaffoldedSeqLen);
return(NumReported);
}

#ifdef _WIN32
unsigned __stdcall ThreadedProcScaffolds(void * pThreadPars)
#else
void * ThreadedProcScaffolds(void * pThreadPars)
#endif
{
int Rslt = 0;
tsThreadProcScaffoldsPars *pPars = (tsThreadProcScaffoldsPars *)pThreadPars; // makes it easier not having to deal with casts!
CScaffolder *pThis = (CScaffolder *)pPars->pThis;
Rslt = pThis->ProcScaffolds(pPars);
pPars->Rslt = Rslt;
#ifdef _WIN32
_endthreadex(0);
return(eBSFSuccess);
#else
pthread_exit(NULL);
#endif
}


#ifdef _WIN32
unsigned __stdcall ThreadedProcContained(void * pThreadPars)
#else
void * ThreadedProcContained(void * pThreadPars)
#endif
{
int Rslt = 0;
tsThreadProcScaffoldsPars *pPars = (tsThreadProcScaffoldsPars *)pThreadPars; // makes it easier not having to deal with casts!
CScaffolder *pThis = (CScaffolder *)pPars->pThis;
Rslt = pThis->ProcContained(pPars);
pPars->Rslt = Rslt;
#ifdef _WIN32
_endthreadex(0);
return(eBSFSuccess);
#else
pthread_exit(NULL);
#endif
}

// GenSeqEdges
// generate sequence edges - these are overlaps of PE sequences onto SE sequences
int
CScaffolder::GenSeqEdges(int CurPhase,				// 0 if lookinging for sense PE overlaps onto sense targets; 1 if looking for antisense (RevCpl'd) PE overlaps onto sense targets
						int SubsKbp,				// allow this many induced substitutions per 1Kbp overlapping sequence fragments (defaults to 0, range 1..10)
						int MaxEnd12Subs,			// allow at the initial 12bp (2 x hexamer primer artefacts) of the 5' or 3' of overlap to have this many base mismatches in addition to the overall allowed MaxSubs 
						int MinReqPESepDist,		// PE start sites expected to be separated by at least this many bases
						int MaxReqPESepDist)		// PE start sites expected be separated by no more than this many bases
{
tsThreadProcScaffoldsPars *pThreadParams;
tsThreadProcScaffoldsPars *pCurThread;
int ThreadIdx;
int NumThreads;
tSeqID CurStartSeqID;

UINT32 CurNumProcessed;
UINT32 PrevNumProcessed = 0;
UINT32 NumOverlapped = 0;

gDiagnostics.DiagOut(eDLDiag,gszProcName,"Starting sequence edge build ...");

NumThreads = m_NumThreads;

// balance the number threads vs the number of sequences to minimise the thread startup costs
if(m_Sequences.NumSeqs2Assemb < 5)
	NumThreads = 1;
else
	{
	NumThreads = (m_Sequences.NumSeqs2Assemb + 4) / 5;
	if(NumThreads > m_NumThreads)
		NumThreads = m_NumThreads;
	}

if((pThreadParams = new tsThreadProcScaffoldsPars[NumThreads])==NULL)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Unable to allocate memory for threads...");
	Reset(false);
	return(eBSFerrMem);
	}
memset(pThreadParams,0,sizeof(tsThreadProcScaffoldsPars) * NumThreads);

CurStartSeqID = 1;

m_Sequences.NumProcessed = 0;
m_Sequences.NumDuplicates = 0;
m_Sequences.NumOverlapping = 0;
m_Sequences.NumOverlapped = 0;

m_FinalProcSeqID = m_Sequences.NumSeqs2Assemb;
m_NextProcSeqID = 1;
m_StartProcSeqID = 1;
if(NumThreads == 1)
	m_NumProcSeqIDs = cMaxMultiSeqFlags;
else
	{
	m_NumProcSeqIDs = min(cMaxMultiSeqFlags,(m_Sequences.NumSeqs2Assemb + NumThreads - 1)/(NumThreads * 10));
	if(m_NumProcSeqIDs < 16)
		m_NumProcSeqIDs = 16;
	}


m_ThreadsProcessing = NumThreads;
pCurThread = pThreadParams;
#ifndef _WIN32
	// set the default stack 
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
	pCurThread->bPESeqs = m_Sequences.bPESeqs ? true : false;
	pCurThread->pThis = this;

	pCurThread->CurPhase = CurPhase;						// 0 if targets indexed on sense; 1 if PE probes are to be RevCpl
	pCurThread->MaxSubs1K = SubsKbp;						// allowing for sub max per 1K overlapping bases
	pCurThread->MaxEnd12Subs = 	MaxEnd12Subs;				// allowing for end hexamer primer artefacts in end 12 bases

	pCurThread->MinReqPESepDist = MinReqPESepDist;			// PE start sites expected to be separated by at least this many bases
	pCurThread->MaxReqPESepDist = MaxReqPESepDist;			// PE start sites expected be separated by no more than this many bases

	pCurThread->MaxProbeSeqWrds = cMaxOvrlapSeqWrds;
	pCurThread->AllocProbeSeqLen = pCurThread->MaxProbeSeqWrds * sizeof(tSeqWrd4);
	pCurThread->MaxMateSeqWrds = 0;
	pCurThread->AllocMateSeqLen = 0;

	pCurThread->pProbeSeqWrds = new UINT8 [pCurThread->AllocProbeSeqLen];
	pCurThread->pMateSeqWrds = NULL;
	if(CurPhase > 0)
		{
		pCurThread->pRevCplProbeSeqWrds = new UINT8 [pCurThread->AllocProbeSeqLen];
		pCurThread->pRevCplMateSeqWrds = NULL;
		}
	else
		{
		pCurThread->pRevCplProbeSeqWrds = NULL;
		pCurThread->pRevCplMateSeqWrds = NULL;
		}

#ifdef _WIN32
	pCurThread->threadHandle = (HANDLE)_beginthreadex(NULL,cWorkThreadStackSize,ThreadedProcScaffolds,pCurThread,0,&pCurThread->threadID);
#else
	pCurThread->threadRslt = pthread_create (&pCurThread->threadID , &threadattr , ThreadedProcScaffolds , pCurThread);
#endif
	}
#ifndef _WIN32
pthread_attr_destroy(&threadattr);		// no longer required
#endif

// allow threads a few seconds to startup
#ifdef _WIN32
	Sleep(2000);
#else
	sleep(2);
#endif
gDiagnostics.DiagOut(eDLDiag,gszProcName,"Progress: 0 sequences processed");
CurNumProcessed = 0;
PrevNumProcessed = 0;
NumOverlapped = 0;
pCurThread = pThreadParams;
for(ThreadIdx = 0; ThreadIdx < NumThreads; ThreadIdx++,pCurThread++)
	{
#ifdef _WIN32
	while(WAIT_TIMEOUT == WaitForSingleObject(pCurThread->threadHandle, 60000 * 10))
		{
		AcquireLock(false);
		CurNumProcessed = m_NextProcSeqID;
		NumOverlapped = (UINT32)m_Sequences.NumOverlapping;
		ReleaseLock(false);
		if(CurNumProcessed > PrevNumProcessed)
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"Progress: %u sequences processed with %u edges",CurNumProcessed,NumOverlapped);
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
		CurNumProcessed = m_NextProcSeqID;
		NumOverlapped = m_Sequences.NumOverlapping;
		ReleaseLock(false);
		if(CurNumProcessed > PrevNumProcessed)
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"Progress: %u sequences processed with %u edges",CurNumProcessed,NumOverlapped);
		PrevNumProcessed = CurNumProcessed;
		ts.tv_sec += 60;
		}
#endif
	if(pCurThread->pProbeSeqWrds != NULL)
		delete (UINT8 *)pCurThread->pProbeSeqWrds;
	if(pCurThread->pRevCplProbeSeqWrds != NULL)
		delete (UINT8 *)pCurThread->pRevCplProbeSeqWrds;
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Completed: %u sequences processed with %u edges",m_Sequences.NumProcessed,m_Sequences.NumOverlapping);

if(pThreadParams)
	delete pThreadParams;

return(eBSFSuccess);
}


// MarkContainedSeqs
// Mark sequences fully contained by other SE sequences
int
CScaffolder::MarkContainedSeqs(int CurPhase,		// 0 if targets indexed on sense; 1 if targets have been RevCpl and then indexed
						bool bPEOnly,				// if false then identify both SE and PE (both ends) sequences which are fully contained, if true then identify contained PEs only 
						int SubsKbp,				// allow this many induced substitutions per 1Kbp overlapping sequence fragments (defaults to 0, range 1..10)
						int MaxEnd12Subs,			// allow at the initial 12bp (2 x hexamer primer artefacts) of the 5' or 3' of overlap to have this many base mismatches in addition to the overall allowed MaxSubs 
						int MinReqPESepDist,		// PE start sites expected to be separated by at least this many bases
						int MaxReqPESepDist)		// PE start sites expected be separated by no more than this many bases
{
tsThreadProcScaffoldsPars *pThreadParams;
tsThreadProcScaffoldsPars *pCurThread;
int ThreadIdx;
int NumThreads;
tSeqID CurStartSeqID;

UINT32 CurNumProcessed;
UINT32 PrevNumProcessed = 0;
UINT32 NumOverlapped = 0;

NumThreads = m_NumThreads;

// balance the number threads vs the number of sequences to minimise the thread startup costs
if(m_Sequences.NumSeqs2Assemb < 5)
	NumThreads = 1;
else
	{
	NumThreads = (m_Sequences.NumSeqs2Assemb + 4) / 5;
	if(NumThreads > m_NumThreads)
		NumThreads = m_NumThreads;
	}

if((pThreadParams = new tsThreadProcScaffoldsPars[NumThreads])==NULL)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Unable to allocate memory for threads...");
	Reset(false);
	return(eBSFerrMem);
	}
memset(pThreadParams,0,sizeof(tsThreadProcScaffoldsPars) * NumThreads);

CurStartSeqID = 1;

m_Sequences.NumProcessed = 0;
m_Sequences.NumDuplicates = 0;
m_Sequences.NumOverlapping = 0;
m_Sequences.NumOverlapped = 0;

m_FinalProcSeqID = m_Sequences.NumSeqs2Assemb;
m_NextProcSeqID = 1;
m_StartProcSeqID = 1;
if(NumThreads == 1)
	m_NumProcSeqIDs = cMaxMultiSeqFlags;
else
	{
	m_NumProcSeqIDs = min(cMaxMultiSeqFlags,(m_Sequences.NumSeqs2Assemb + NumThreads - 1)/(NumThreads * 10));
	if(m_NumProcSeqIDs < 16)
		m_NumProcSeqIDs = 16;
	}


m_ThreadsProcessing = NumThreads;
pCurThread = pThreadParams;
#ifndef _WIN32
	// set the default stack 
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
	pCurThread->bPESeqs = m_Sequences.bPESeqs ? true : false;
	pCurThread->pThis = this;

	pCurThread->CurPhase = CurPhase;						// 0 if targets indexed on sense; 1 if targets have been RevCpl and then indexed
	pCurThread->bPEOnly = bPEOnly;							// if false then identify both SE and PE (both ends) sequences which are fully contained, if true then identify contained PEs only

	pCurThread->MaxSubs1K = SubsKbp;						// allowing for sub max per 1K overlapping bases
	pCurThread->MaxEnd12Subs = 	MaxEnd12Subs;				// allowing for end hexamer primer artefacts in end 12 bases

	pCurThread->MinReqPESepDist = MinReqPESepDist;			// PE start sites expected to be separated by at least this many bases
	pCurThread->MaxReqPESepDist = MaxReqPESepDist;			// PE start sites expected be separated by no more than this many bases

	pCurThread->MaxProbeSeqWrds = cMaxOvrlapSeqWrds;
	pCurThread->AllocProbeSeqLen = pCurThread->MaxProbeSeqWrds * sizeof(tSeqWrd4);
	pCurThread->MaxMateSeqWrds = cMaxOvrlapSeqWrds;
	pCurThread->AllocMateSeqLen = pCurThread->MaxMateSeqWrds * sizeof(tSeqWrd4);

	pCurThread->pProbeSeqWrds = new UINT8 [pCurThread->AllocProbeSeqLen];
	pCurThread->pMateSeqWrds = new UINT8 [pCurThread->AllocMateSeqLen];
	if(CurPhase > 0)
		{
		pCurThread->pRevCplProbeSeqWrds = new UINT8 [pCurThread->AllocProbeSeqLen];
		pCurThread->pRevCplMateSeqWrds = new UINT8 [pCurThread->AllocMateSeqLen];
		}
	else
		{
		pCurThread->pRevCplProbeSeqWrds = NULL;
		pCurThread->pRevCplMateSeqWrds = NULL;
		}


#ifdef _WIN32
	pCurThread->threadHandle = (HANDLE)_beginthreadex(NULL,cWorkThreadStackSize,ThreadedProcContained,pCurThread,0,&pCurThread->threadID);
#else
	pCurThread->threadRslt = pthread_create (&pCurThread->threadID , &threadattr , ThreadedProcContained , pCurThread);
#endif
	}
#ifndef _WIN32
pthread_attr_destroy(&threadattr);		// no longer required
#endif

// allow threads a few seconds to startup
#ifdef _WIN32
	Sleep(2000);
#else
	sleep(2);
#endif
gDiagnostics.DiagOut(eDLDiag,gszProcName,"Progress: 0 sequences processed");
CurNumProcessed = 0;
PrevNumProcessed = 0;
NumOverlapped = 0;
pCurThread = pThreadParams;
for(ThreadIdx = 0; ThreadIdx < NumThreads; ThreadIdx++,pCurThread++)
	{
#ifdef _WIN32
	while(WAIT_TIMEOUT == WaitForSingleObject(pCurThread->threadHandle, 60000 * 10))
		{
		AcquireLock(false);
		CurNumProcessed = m_NextProcSeqID;
		NumOverlapped = (UINT32)m_Sequences.NumOverlapping;
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
		CurNumProcessed = m_NextProcSeqID;
		NumOverlapped = m_Sequences.NumOverlapping;
		ReleaseLock(false);
		if(CurNumProcessed > PrevNumProcessed)
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"Progress: %u sequences processed",CurNumProcessed);
		PrevNumProcessed = CurNumProcessed;
		ts.tv_sec += 60;
		}
#endif
	if(pCurThread->pProbeSeqWrds != NULL)
		delete (UINT8 *)pCurThread->pProbeSeqWrds;
	if(pCurThread->pRevCplProbeSeqWrds != NULL)
		delete (UINT8 *)pCurThread->pRevCplProbeSeqWrds;
	if(pCurThread->pMateSeqWrds != NULL)
		delete (UINT8 *)pCurThread->pMateSeqWrds;
	if(pCurThread->pRevCplMateSeqWrds != NULL)
		delete (UINT8 *)pCurThread->pRevCplMateSeqWrds;
	}

if(pThreadParams)
	delete pThreadParams;

return(m_Sequences.NumOverlapping);
}


int
CScaffolder::ProcContained(tsThreadProcScaffoldsPars *pPars)
{
int CmpRslt;
int SubOfs;

int MinOvrlp;
int RelSfxWrdOfs;
int TargOverlapLen;

tSeqID StartingSeqID;
tSeqID EndingSeqID;
int TooMAnyWarnings;
tSeqID SeqID;

tSeqID ProbeSeqID;
UINT32 ProbeLen;
UINT32 ProbeFlags;
tSeqWrd4 *pProbeHdr;
tSeqWrd4 *pProbeSeqWrds;

tSeqID MateSeqID;
UINT32 MateLen;
UINT32 MateFlags;
tSeqWrd4 *pMateHdr;
tSeqWrd4 *pMateSeqWrds;

int PutMatchOfs;
int PutMatchLen;

int MaxSubOfs;

bool bTargetContainsProbe;
bool bProbeContainsTarget;

int MateTargOfs;
bool bPE2Match;
bool bContained;
UINT32 NumProcessed;
UINT32 NumContained;

UINT64 SfxWrdIdx;


UINT32 PrevProbeLen;

tSeqID CurTargSeqID;
UINT32 CurTargSeqLen;
UINT32 TargFlags;
tSeqWrd4 *pHit;
tSeqWrd4 *pCurTargStartSeqWrd;
int TargSubs;

gDiagnostics.DiagOut(eDLDebug,gszProcName,"Thread %d startup for contained sequence identification...",pPars->ThreadIdx);

NumProcessed = 0;
TooMAnyWarnings = 0;
NumContained = 0;

time_t Started = time(0);
while(GetSeqProcRange(&StartingSeqID,&EndingSeqID,cMaxMultiSeqFlags) > 0)
	{
	PrevProbeLen = 0;

	for(SeqID = StartingSeqID; SeqID <= EndingSeqID; SeqID++)
		{
		NumProcessed+=1;				// processing PE1, PE2, or SE, may be a little presumptuous but would only not be fully processed if a fatal error whilst processing ..
		if(!(NumProcessed % 100))
			{
			time_t Now = time(0);
			unsigned long ElapsedSecs = (unsigned long) (Now - Started);
			if(ElapsedSecs >= 30)
				{
				pPars->NumProcessed += NumProcessed;
				pPars->NumOverlapped += NumContained;
				AcquireLock(true);
				m_Sequences.NumProcessed += NumProcessed;
				m_Sequences.NumOverlapping += NumContained;
				ReleaseLock(true);
				NumProcessed = 0;
				NumContained = 0;
				Started = Now;
				}
			}

		pProbeHdr = (tSeqWrd4 *)GetSeqHeader(SeqID,NULL,&ProbeFlags,&ProbeLen,false);
		if(pProbeHdr == NULL)				// serious problem if can't get header...
			{
			if((TooMAnyWarnings+=1) < 10)
				gDiagnostics.DiagOut(eDLWarn,gszProcName,"Thread %d Couldn't find sequence header for known sequence %u...",pPars->ThreadIdx,SeqID);
			continue;
			}

		if(ProbeFlags & cFlgNoProc)			// no further processing required on this sequence if cFlgNoProc
			continue;

		if(ProbeFlags & cFlgContainRemove)			// no further processing required on this sequence cFlgContainRemove already set
			{	
			if((ProbeFlags & (cFlgSeqPE | cFlgSeqPE2)) == cFlgSeqPE)	// initially only interested in PE1, containments on PE1 by a SE are then explored to see if PE2 also contained
				{
				SeqID++;								// skip PE2
				NumProcessed += 1;
				}
			continue;
			}

		if((ProbeFlags & (cFlgSeqPE | cFlgSeqPE2)) != cFlgSeqPE)	// initially only interested in PE1, containments on PE1 by a SE are then explored to see if PE2 also contained
			continue;

		if(pPars->CurPhase == 0)			// exchange sequence identifiers if not processing sense overlapping sense, once these sequences have been loaded then they will be RevCpl'd
			{
			ProbeSeqID = SeqID++;
			MateSeqID = SeqID;
			pMateHdr = (tSeqWrd4 *)GetSeqHeader(MateSeqID,NULL,&MateFlags,&MateLen,false);
			if(pMateHdr == NULL)				// serious problem if can't get header...
				{
				if((TooMAnyWarnings+=1) < 10)
					gDiagnostics.DiagOut(eDLWarn,gszProcName,"Thread %d Couldn't find sequence header for known sequence %u...",pPars->ThreadIdx,MateSeqID);
				continue;
				}
			}
		else
			{
			pMateHdr = pProbeHdr;
			MateFlags = ProbeFlags;
			MateLen = ProbeLen;
			MateSeqID = SeqID++;
			ProbeSeqID = SeqID;
			pProbeHdr = (tSeqWrd4 *)GetSeqHeader(ProbeSeqID,NULL,&ProbeFlags,&ProbeLen,false);
			if(pMateHdr == NULL)				// serious problem if can't get header...
				{
				if((TooMAnyWarnings+=1) < 10)
					gDiagnostics.DiagOut(eDLWarn,gszProcName,"Thread %d Couldn't find mate sequence header for known sequence %u...",pPars->ThreadIdx,ProbeSeqID);
				continue;
				}
			}
		NumProcessed+=1;							// processing both PE1 and Pe2, may be a little presumptuous but would only not be fully processed if a fatal error whilst processing ..
		pProbeSeqWrds = pProbeHdr;
		pProbeHdr -= 3;
		pMateSeqWrds = pMateHdr;
		pMateHdr -= 3;

		// make a copy of probe  sequence as will be slicing and dicing when subsequencing...
		GetSeqWrdSubSeq(0,ProbeLen, pProbeSeqWrds, (tSeqWrd4 *)pPars->pProbeSeqWrds);

		// sequences loaded, RevCpl if not processing sense overlapping sense
		if(pPars->CurPhase > 0)
			{
			PackedRevCpl((tSeqWrd4 *)pPars->pProbeSeqWrds);
			GetSeqWrdSubSeq(0,ProbeLen, (tSeqWrd4 *)pPars->pProbeSeqWrds, (tSeqWrd4 *)pPars->pRevCplProbeSeqWrds);
			pProbeSeqWrds = (tSeqWrd4 *)pPars->pRevCplProbeSeqWrds;

			GetSeqWrdSubSeq(0,MateLen, pMateSeqWrds, (tSeqWrd4 *)pPars->pMateSeqWrds);
			PackedRevCpl((tSeqWrd4 *)pPars->pMateSeqWrds);
			pMateSeqWrds = (tSeqWrd4 *)pPars->pMateSeqWrds;
			}

		if(pPars->MaxSubs1K || pPars->MaxEnd12Subs)
			{
			MinOvrlp = min(max(ProbeLen,30),60);			// looking for minimum overlaps of at least 60bp as seed exact matches which will be then extended out to full length allowing any mismatches
			MaxSubOfs =  min(90,(int)ProbeLen - MinOvrlp);
			}
		else
			{
			MaxSubOfs =  15;
			MinOvrlp = ProbeLen - MaxSubOfs;	
			}
		bContained = false;
		for(SubOfs = 0; SubOfs < MaxSubOfs && !bContained; SubOfs++) 
			{
			if(SubOfs)
				ShfLeftPackedSeq(ProbeLen + 1 - SubOfs,(tSeqWrd4 *)pPars->pProbeSeqWrds);
			
			SfxWrdIdx =	LocateFirstExact(m_Sequences.SfxElSize,	// sizeof elements in pSfxArray - currently will be either 4 or 5 bytes
							pPars->pProbeSeqWrds,				// pts to probes flank subsequence
							MinOvrlp,							// probe length (in bases, not tSeqWrd4's) required to minimally exactly match over
							m_Sequences.pSeqs2Assemb,			// target sequence
							(UINT8 *)m_Sequences.pSuffixArray,	// target sequence suffix array
							0,									// low index in pSfxArray
							m_Sequences.NumSuffixEls-1);		// high index in pSfxArray

			if(!SfxWrdIdx)		// will be 0 if unable to find any targets prefix sequences exactly matching that of the probes subsequence (pPars->pOverlapSeq)
				continue;

			// have at least one exact seed matching, iterate over these seed matches extending and allowing any mismatches
			while(!bContained) {
				CmpRslt = 0;
				if((pHit = (tSeqWrd4 *)SfxIdxToFirstSeqWrd(SfxWrdIdx++,&RelSfxWrdOfs))==NULL)  // get ptr to targets starting (immediately following header) seqword containing 5' bases 
					break;									   // should only be a NULL if exhusted all potential overlaps

				// if overlap is less than required minimum overlap then must have exhusted all potential overlaps for current probe
				if((TargOverlapLen = GetExactMatchLen((tSeqWrd4 *)pPars->pProbeSeqWrds,&pHit[RelSfxWrdOfs],MinOvrlp)) < MinOvrlp)
					break;

				// which target was hit and what length is it?
				pCurTargStartSeqWrd = (tSeqWrd4 *)GetSeqHeader(pHit,&CurTargSeqID,NULL,&TargFlags,(UINT32 *)&CurTargSeqLen,false);

				// can't accept self hits
				if(CurTargSeqID == ProbeSeqID)				// if self hit then not interested, try for another target sequence
					continue;

				if(TargFlags & cFlgSeqPE)					// if hit was to another PE then not interested
					continue;

				bTargetContainsProbe = false;
				bProbeContainsTarget = false;

				PutMatchOfs =  (RelSfxWrdOfs * 15) - SubOfs;	// if PutMatchOfs < 0 then probe was overlapping onto target, we need the target to be overlapping the probe
				if(PutMatchOfs < 0 || (PutMatchOfs + ProbeLen) > CurTargSeqLen)
					continue;

				PutMatchLen = min(ProbeLen,CurTargSeqLen);

				// check if putative match really is a full match
				if(!IsMatched(PutMatchLen,0,ProbeLen,pProbeSeqWrds,PutMatchOfs,CurTargSeqLen,pCurTargStartSeqWrd,pPars->MaxSubs1K,pPars->MaxEnd12Subs,&TargSubs))
					continue;

				// check if the other mate PE end would also fit within target
				int MinPutMateOfs;
				int MaxPutMateOfs;

				
				if(pPars->MinReqPESepDist > (int)MateLen)
					MinPutMateOfs = pPars->MinReqPESepDist - (int)MateLen;
				else
					MinPutMateOfs = 0;
				MinPutMateOfs += PutMatchOfs;
					
				if(((UINT32)MinPutMateOfs +  MateLen) > CurTargSeqLen)
					continue;

				MaxPutMateOfs = PutMatchOfs + pPars->MaxReqPESepDist - (int)MateLen;

				for(MateTargOfs = MinPutMateOfs; MateTargOfs <= MaxPutMateOfs; MateTargOfs++) 
					if(bPE2Match = IsMatched(MateLen,0,MateLen,pMateSeqWrds,MateTargOfs,CurTargSeqLen,pCurTargStartSeqWrd,pPars->MaxSubs1K,pPars->MaxEnd12Subs,&TargSubs))
						break;
				if(!bPE2Match)
					continue;

				UpdateSeqHeaderFlags(pProbeHdr,cFlgContainRemove,0,true);
				UpdateSeqHeaderFlags(pMateHdr,cFlgContainRemove,0,true);
				NumContained += 1;
				bContained = true;
				}
			}
		}
	}

pPars->NumProcessed += NumProcessed;
pPars->NumOverlapped += NumContained;
AcquireLock(true);
m_Sequences.NumProcessed += NumProcessed;
m_Sequences.NumOverlapping += NumContained;
ReleaseLock(true);
gDiagnostics.DiagOut(eDLDebug,gszProcName,"Thread %d contained sequence identification completed",pPars->ThreadIdx);
return(1);		// success
}

// ProcScaffolds
// Identified overlaps used to construct scaffolding graph
int
CScaffolder::ProcScaffolds(tsThreadProcScaffoldsPars *pPars)
{
int CmpRslt;
int SubOfs;

int MinOvrlp;

int RelSfxWrdOfs;
int TargOverlapLen;

tSeqID StartingSeqID;
tSeqID EndingSeqID;
int TooMAnyWarnings;
tSeqID SeqID;
tSeqID ProbeSeqID;
UINT32 ProbeLen;
UINT32 ProbeFlags;
tSeqWrd4 *pProbeHdr;
tSeqWrd4 *pProbeSeqWrds;

UINT32 MateSeqID;

UINT32 NumProcessed;

UINT32 NumOverlapping;
UINT32 NumOverlapped;

UINT64 SfxWrdIdx;
UINT32 PrevProbeLen;

int PutMatchOfs;
int PutMatchLen;

int MaxSubOfs;

int ProvTargSeqID;
int	ProvTargSubs;

bool bTargetContainsProbe;
bool bProbeContainsTarget;


tSeqID CurTargSeqID;
UINT32 CurTargSeqLen;
UINT32 TargFlags;
tSeqWrd4 *pHit;
tSeqWrd4 *pCurTargStartSeqWrd;
int TargSubs;

tsSeqEdge *pProvAcceptedEdge;
tsSeqEdge ProvAcceptedEdges[cMaxOverlapsProbeTarg];
int ProvIdx;
int NumProvAcceptedEges;
UINT32 NumAdd2Graph;

gDiagnostics.DiagOut(eDLDebug,gszProcName,"Thread %d startup for overlap identification...",pPars->ThreadIdx);

NumProcessed = 0;
TooMAnyWarnings = 0;
NumOverlapping = 0;
NumOverlapped = 0;
NumAdd2Graph = 0;

time_t Started = time(0);
while(GetSeqProcRange(&StartingSeqID,&EndingSeqID,cMaxMultiSeqFlags) > 0)
	{
	PrevProbeLen = 0;
	for(SeqID = StartingSeqID; SeqID <= EndingSeqID; SeqID++)
		{
		NumProcessed+=1;
		if(!(NumProcessed % 100))
			{
			time_t Now = time(0);
			unsigned long ElapsedSecs = (unsigned long) (Now - Started);
			if(ElapsedSecs >= 30)
				{
				pPars->NumProcessed += NumProcessed;
				pPars->NumOverlapped += NumOverlapped;
				AcquireLock(true);
				m_Sequences.NumProcessed += NumProcessed;
				m_Sequences.NumOverlapping += NumAdd2Graph;
				ReleaseLock(true);
				NumProcessed = 0;
				NumAdd2Graph = 0;
				Started = Now;
				}
			}

		pProbeHdr = (tSeqWrd4 *)GetSeqHeader(SeqID,NULL,&ProbeFlags,&ProbeLen,false);
		if(pProbeHdr == NULL)						// serious problem if can't get header...
			{
			if((TooMAnyWarnings+=1) < 10)
				gDiagnostics.DiagOut(eDLWarn,gszProcName,"Thread %d Couldn't find sequence header for known sequence %u...",pPars->ThreadIdx,SeqID);
			continue;
			}

		if(ProbeFlags & cFlgNoProc)					// no further processing required on this sequence if cFlgNoProc
			continue;

		if(!(ProbeFlags & (cFlgSeqPE | cFlgSeqPE2))) // looking for PE overlaps onto SE, so if SE as a probe then not interested
			continue;

		ProbeSeqID = SeqID;
		MateSeqID = ProbeFlags & cFlgSeqPE2 ? ProbeSeqID - 1 : ProbeSeqID + 1;
		pProbeSeqWrds = pProbeHdr;
		pProbeHdr -= 3;

		// make a copy of probe  sequence as will be slicing and dicing when subsequencing...
		GetSeqWrdSubSeq(0,ProbeLen, pProbeSeqWrds, (tSeqWrd4 *)pPars->pProbeSeqWrds);

		// sequence loaded, RevCpl if not processing sense overlapping sense
		if(pPars->CurPhase > 0)
			{
			PackedRevCpl((tSeqWrd4 *)pPars->pProbeSeqWrds);
			GetSeqWrdSubSeq(0,ProbeLen, (tSeqWrd4 *)pPars->pProbeSeqWrds, (tSeqWrd4 *)pPars->pRevCplProbeSeqWrds);
			pProbeSeqWrds = (tSeqWrd4 *)pPars->pRevCplProbeSeqWrds;
			}


		if(pPars->MaxSubs1K || pPars->MaxEnd12Subs)
			{
			MinOvrlp = min(max(ProbeLen,30),60);			// looking for minimum overlaps of at least 60bp as seed exact matches which will be then extended out to full length allowing any mismatches
			MaxSubOfs =  min(90,(int)ProbeLen - MinOvrlp);
			}
		else
			{
			MaxSubOfs =  15;
			MinOvrlp = ProbeLen - MaxSubOfs;	
			}

		ProvTargSeqID = 0;
		ProvTargSubs = 0;

		NumProvAcceptedEges = 0;
		for(SubOfs = 0; SubOfs < MaxSubOfs && NumProvAcceptedEges != -1; SubOfs++) 
			{
			if(SubOfs)
				ShfLeftPackedSeq(ProbeLen + 1 - SubOfs,(tSeqWrd4 *)pPars->pProbeSeqWrds);

			SfxWrdIdx =	LocateFirstExact(m_Sequences.SfxElSize,	// sizeof elements in pSfxArray - currently will be either 4 or 5 bytes
							pPars->pProbeSeqWrds,				// pts to probes flank subsequence
							MinOvrlp,							// probe length (in bases, not tSeqWrd4's) required to minimally exactly match over
							m_Sequences.pSeqs2Assemb,			// target sequence
							(UINT8 *)m_Sequences.pSuffixArray,	// target sequence suffix array
							0,									// low index in pSfxArray
							m_Sequences.NumSuffixEls-1);		// high index in pSfxArray

			if(!SfxWrdIdx)		// will be 0 if unable to find any targets prefix sequences exactly matching that of the probes subsequence (pPars->pOverlapSeq)
				continue;

			// iterating over all overlapped target sequences
			do {
				CmpRslt = 0;
				if((pHit = (tSeqWrd4 *)SfxIdxToFirstSeqWrd(SfxWrdIdx++,&RelSfxWrdOfs))==NULL)  // get ptr to targets starting (immediately following header) seqword containing 5' bases 
					break;									   // should only be a NULL if exhusted all potential overlaps

				// if overlap is less than required minimum overlap then must have exhusted all potential overlaps for current probe
				if((TargOverlapLen = GetExactMatchLen((tSeqWrd4 *)pPars->pProbeSeqWrds,&pHit[RelSfxWrdOfs],MinOvrlp)) < MinOvrlp)
					break;

				// which target was hit and what length is it?
				pCurTargStartSeqWrd = (tSeqWrd4 *)GetSeqHeader(pHit,&CurTargSeqID,NULL,&TargFlags,(UINT32 *)&CurTargSeqLen,false);

				// can't accept self hits
				if(CurTargSeqID == ProbeSeqID)				// if self hit then not interested, try for another target sequence
					continue;

				if(TargFlags & cFlgSeqPE)					// if hit was to another PE then not interested
					continue;

				bTargetContainsProbe = false;
				bProbeContainsTarget = false;

				PutMatchOfs =  (RelSfxWrdOfs * 15) - SubOfs;	// if PutMatchOfs < 0 then probe was overlapping onto target, we need the traget to be overlapping the probe
				if(PutMatchOfs < 0 || (PutMatchOfs + ProbeLen) > CurTargSeqLen)
					continue;

				PutMatchLen = min(ProbeLen,CurTargSeqLen);

				// check if putative match really is a full match
				if(!IsMatched(PutMatchLen,0,ProbeLen,pProbeSeqWrds,PutMatchOfs,CurTargSeqLen,pCurTargStartSeqWrd,pPars->MaxSubs1K,pPars->MaxEnd12Subs,&TargSubs))
					continue;

				if(pPars->CurPhase == 0)	// sense onto sense, PE1 needs to be within MaxReqPESepDist of 3' SE and PE2 within MaxReqPESepDist of 5' SE
					{
					if(ProbeFlags & cFlgSeqPE2)
						{
						if(((UINT32)PutMatchOfs + ProbeLen) > (UINT32)pPars->MaxReqPESepDist)		// needs to have overlapped within the maximum PE insert size
							continue;
						}
					else
						if((CurTargSeqLen - (UINT32)PutMatchOfs) > (UINT32)pPars->MaxReqPESepDist)		// needs to have overlapped within the maximum PE insert size
							continue;
					}
				else						// antisense onto sense, PE1 needs to be within MaxReqPESepDist of 5' SE and PE2 within MaxReqPESepDist of 3' SE                          
					{
					if(ProbeFlags & cFlgSeqPE2)
						{
						if((CurTargSeqLen - (UINT32)PutMatchOfs) > (UINT32)pPars->MaxReqPESepDist)		// needs to have overlapped within the maximum PE insert size
							continue;
						}
					else
						if(((UINT32)PutMatchOfs + ProbeLen) > (UINT32)pPars->MaxReqPESepDist)		// needs to have overlapped within the maximum PE insert size
							continue;
					}


				if(ProvTargSeqID != 0 &&  TargSubs > ProvTargSubs)
					continue;

				// provisionally accepted this hit of a PE end onto SE
				// can only finally accept if hit was to same SE as other hits from same PE or if fewer subs would be required
				if(ProvTargSeqID != 0 &&  CurTargSeqID != ProvTargSeqID)
					{
					if(TargSubs == 0 && ProvTargSubs == 0)
						{
						NumProvAcceptedEges = -1;
						break;
						}

					NumProvAcceptedEges = 0;
					if(TargSubs == ProvTargSubs) // TargSubs same as ProvTargSubs and both more than 0 so keep exploring for a target requiring fewer subs
						{
						ProvTargSeqID = 0;
						ProvTargSubs = 0;
						continue;
						}
					}

				if(NumProvAcceptedEges)		// check if already provisionally accepted this overlap
					{
					pProvAcceptedEdge = ProvAcceptedEdges;
					for(ProvIdx = 0; ProvIdx < NumProvAcceptedEges; ProvIdx++, pProvAcceptedEdge++)
						if(pProvAcceptedEdge->zFlankLen == PutMatchOfs)
							break;
					if(ProvIdx < NumProvAcceptedEges)		// same sequence + loci being hit, or hits to more than one target SE
						continue;
					}
				else
					{
					ProvTargSeqID = CurTargSeqID;
					ProvTargSubs = TargSubs;
					}
				if(TargSubs < ProvTargSubs)
					{			
					ProvTargSubs = TargSubs;
					NumProvAcceptedEges = 0;
					}

					// if already at limit of allowed overlaps of a PE sequence onto a SE sequences then continue looking in case the same PE overlaps onto a different SE 
				if(NumProvAcceptedEges >= cMaxOverlapsProbeTarg)
					continue;

				pProvAcceptedEdge = &ProvAcceptedEdges[NumProvAcceptedEges++];
				memset(pProvAcceptedEdge,0,sizeof(tsSeqEdge));
				pProvAcceptedEdge->TargetSeqID = CurTargSeqID;
				pProvAcceptedEdge->TargetLen = CurTargSeqLen;
				pProvAcceptedEdge->zFlankLen = PutMatchOfs;
				pProvAcceptedEdge->ProbeMateSeqID = MateSeqID;
				pProvAcceptedEdge->ProcPhase = pPars->CurPhase;
				pProvAcceptedEdge->ProbeSeqID = SeqID;
				pProvAcceptedEdge->ProbeLen = ProbeLen;
				pProvAcceptedEdge->FlgRevCpltd = pPars->CurPhase == 0 ? 0 : 1;
				if(ProbeFlags & cFlgSeqPE2)
					pProvAcceptedEdge->FlgIsPE2 = 1;
				pProvAcceptedEdge->zOverlapLen = ProbeLen;
				pProvAcceptedEdge->NumSubs = TargSubs;
				}
			while(1);
			}

		if(NumProvAcceptedEges > 0)
			{
			pProvAcceptedEdge = ProvAcceptedEdges;
			for(ProvIdx = 0; ProvIdx < NumProvAcceptedEges; ProvIdx++, pProvAcceptedEdge++)
				{
				if(pProvAcceptedEdge->NumSubs <= ProvTargSubs)
					{
					AddOverlapEdge(pProvAcceptedEdge);
					NumAdd2Graph += 1;
					NumOverlapped += 1;
					}
				}
			}
		NumProvAcceptedEges = 0;
		}
	}

pPars->NumProcessed += NumProcessed;
pPars->NumOverlapped += NumOverlapped;
AcquireLock(true);
m_Sequences.NumProcessed += NumProcessed;
m_Sequences.NumOverlapping += NumAdd2Graph;
ReleaseLock(true);
gDiagnostics.DiagOut(eDLDebug,gszProcName,"Thread %d overlapping sequence identification completed",pPars->ThreadIdx);
return(1);		// success
}

teBSFrsltCodes
CScaffolder::AddOverlapEdge(tsSeqEdge *pEdge)  // pts to pre-initialised edge, except for SeqEdgeID which is set within this function
{
tsSeqEdge *pSeqEdge;

AcquireSerialiseSeqFlags();
if(m_pSeqEdges == NULL)
	{
	m_AllocMemSeqEdges = (size_t)sizeof(tsSeqEdge) * cSeqEdges2Alloc;
#ifdef _WIN32
	m_pSeqEdges = (tsSeqEdge *)malloc((size_t)m_AllocMemSeqEdges);	
	if(m_pSeqEdges == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddOverlap: Unable to allocate memory of %llu bytes - %s",m_AllocMemSeqEdges,strerror(errno));
		Reset(false);
		return(eBSFerrMem);
		}
#else
	if((m_pSeqEdges = (tsSeqEdge *)mmap(NULL,m_AllocMemSeqEdges, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0)) == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddOverlap: Unable to allocate memory of %llu bytes - %s",m_AllocMemSeqEdges,strerror(errno));
		m_pSeqEdges = NULL;
		ReleaseSerialiseSeqFlags();
		Reset(false);
		return(eBSFerrMem);
		}
#endif
	m_NumSeqEdges = 0;
	}
else
	{
	if(((m_NumSeqEdges + 10) * sizeof(tsSeqEdge)) >= m_AllocMemSeqEdges)
		{
		size_t	memreq = (size_t)((m_AllocMemSeqEdges * 150) / (UINT64)100);	// increase current allocation by 50%
#ifdef _WIN32
		pSeqEdge = (tsSeqEdge *)realloc(m_pSeqEdges,memreq);
#else
		pSeqEdge = (tsSeqEdge *)mremap(m_pSeqEdges,m_AllocMemSeqEdges,memreq,MREMAP_MAYMOVE);
		if(pSeqEdge == MAP_FAILED)
			pSeqEdge = NULL;
#endif
		if(pSeqEdge == NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddOverlap: Memory re-allocation to %d bytes - %s",memreq,strerror(errno));
			ReleaseSerialiseSeqFlags();
			return(eBSFerrMem);
			}

		m_pSeqEdges = pSeqEdge;
		m_AllocMemSeqEdges = (UINT64)memreq;
		}
	}

pSeqEdge = &m_pSeqEdges[m_NumSeqEdges++];
pEdge->SeqEdgeID = m_NumSeqEdges;
*pSeqEdge = *pEdge;
ReleaseSerialiseSeqFlags();
return(eBSFSuccess);
} 


tVertID							// returned unique Vertex identifier
CScaffolder::AddSeqVertex(tsSeqVertex *pSE)	// SE
{
tsSeqVertex *pSeqVertex;

AcquireSerialiseSeqFlags();
if(m_pSeqVertices == NULL)
	{
	m_AllocMemSeqVertices = (size_t)sizeof(tsSeqVertex) * (100 + m_Sequences.NumPE1Seqs2Assemb);
#ifdef _WIN32
	m_pSeqVertices = (tsSeqVertex *)malloc((size_t)m_AllocMemSeqVertices);	
	if(m_pSeqVertices == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddSeqVertex: Unable to allocate memory of %llu bytes - %s",m_AllocMemSeqVertices,strerror(errno));
		Reset(false);
		return(eBSFerrMem);
		}
#else
	if((m_pSeqVertices = (tsSeqVertex *)mmap(NULL,m_AllocMemSeqVertices, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0)) == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddSeqVertex: Unable to allocate memory of %llu bytes - %s",m_AllocMemSeqVertices,strerror(errno));
		m_pSeqVertices = NULL;
		ReleaseSerialiseSeqFlags();
		Reset(false);
		return(eBSFerrMem);
		}
#endif
	m_NumSeqVertices = 0;
	}
else
	{
	if((m_NumSeqVertices * sizeof(tsSeqVertex)) >= m_AllocMemSeqVertices)
		{
		size_t	memreq = (size_t)((m_AllocMemSeqVertices * 125) / (UINT64)100);	// increase current allocation by 25%
#ifdef _WIN32
		pSeqVertex = (tsSeqVertex *)realloc(m_pSeqVertices,memreq);
#else
		pSeqVertex = (tsSeqVertex *)mremap(m_pSeqVertices,m_AllocMemSeqVertices,memreq,MREMAP_MAYMOVE);
		if(pSeqVertex == MAP_FAILED)
			pSeqVertex = NULL;
#endif
		if(pSeqVertex == NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddSeqVertex: Memory re-allocation to %d bytes - %s",memreq,strerror(errno));
			ReleaseSerialiseSeqFlags();
			return(eBSFerrMem);
			}

		m_pSeqVertices = pSeqVertex;
		m_AllocMemSeqVertices = (UINT64)memreq;
		}
	}

pSeqVertex = &m_pSeqVertices[m_NumSeqVertices++];
memcpy(pSeqVertex,pSE,sizeof(tsSeqVertex));
pSeqVertex->SeqVertID = m_NumSeqVertices;
ReleaseSerialiseSeqFlags();
return(pSeqVertex->SeqVertID);
}


tSeqID			// 0 if no sequence identifier after FromSeqID, otherwise the next ordered edge FromSeqID 
CScaffolder::IterateNextFromSeqID(tSeqID FromSeqID)						// 0 to start from 1st
{
UINT64 EdgeIdx;
tsSeqEdge *pEdge;
if(m_NumSeqEdges == 0 || m_pSeqEdges == NULL)
	return(0);
if(FromSeqID == 0)
	return(m_pSeqEdges[0].ProbeSeqID);

if(FromSeqID < m_LowestEdgeFromSeqID || FromSeqID >= m_HighestEdgeFromSeqID)
	return(0);

if((EdgeIdx = LocateFirstEdgeFromSeqID(FromSeqID)) == 0)
	return(0);
pEdge = &m_pSeqEdges[EdgeIdx-1];
for(;EdgeIdx < m_NumSeqEdges; EdgeIdx += 1, pEdge += 1)
	if(pEdge->ProbeSeqID != FromSeqID)
		return(pEdge->ProbeSeqID);
return(0);
}

tSeqID			// 0 if no sequence identifier after ToSeqID, otherwise the next ordered edge ToSeqID 
CScaffolder::IterateNextToSeqID(tSeqID ToSeqID)	// 0 to start from 1st
{
UINT64 SeqIdx;
tsSeqEdge *pEdge;
if(m_NumSeqEdges == 0 || m_pSeqEdges == NULL || m_ppToSeqEdges == NULL)
	return(0);
if(ToSeqID == 0)
	return(m_ppToSeqEdges[0]->TargetSeqID);

if(ToSeqID < m_LowestEdgeToSeqID || ToSeqID >= m_HighestEdgeToSeqID)
	return(0);

if((SeqIdx = LocateFirstEdgeToSeqID(ToSeqID)) == 0)
	return(0);

for(;SeqIdx < m_NumSeqEdges; SeqIdx += 1)
	{
	pEdge = m_ppToSeqEdges[SeqIdx-1];
	if(pEdge->TargetSeqID != ToSeqID)
		return(pEdge->TargetSeqID);
	}
return(0);
}

tsSeqEdge *			// NULL if unable to locate any with matching FromSeqID				 
CScaffolder::IterateEdgeFromSeqID(tSeqID FromSeqID,	// iterate over sequence edges with matching FromSeqID
								  tEdgeID *pEdgeId)	// set *pSeqIdx to 0 to return 1st sequence edge, returns index to use on next iteration of IterateEdgeFromSeqID()
{
UINT64 EdgeIdx;
tsSeqEdge *pEdge;
if(m_NumSeqEdges == 0 || m_pSeqEdges == NULL || FromSeqID < m_LowestEdgeFromSeqID || FromSeqID > m_HighestEdgeFromSeqID)
	return(NULL);
if(pEdgeId != NULL && (EdgeIdx = *pEdgeId) > 0)
	{
	*pEdgeId = 0;
	if(EdgeIdx >= m_NumSeqEdges)
		return(NULL);
	pEdge = &m_pSeqEdges[EdgeIdx];
	if(pEdge->ProbeSeqID == FromSeqID)
		{
		*pEdgeId = EdgeIdx+1;
		return(pEdge);
		}
	return(NULL);
	}
EdgeIdx = LocateFirstEdgeFromSeqID(FromSeqID);
if(pEdgeId != NULL)
	*pEdgeId = EdgeIdx;
return(EdgeIdx > 0 ? &m_pSeqEdges[EdgeIdx-1] : NULL);
}


tsSeqEdge *			// NULL if unable to locate any with matching ToSeqID				 
CScaffolder::IterateEdgeToSeqID(tSeqID ToSeqID,	// iterate over sequence edges with matching ToSeqID
									tEdgeID *pEdgeID)	// set *pEdgeID to 0 to return 1st sequence edge, returns edge identifier to use on next iteration of IterateEdgeToSeqID()
{
tEdgeID EdgeID;
tsSeqEdge *pEdge;
if(m_NumSeqEdges == 0 || m_pSeqEdges == NULL || m_ppToSeqEdges == NULL || ToSeqID < m_LowestEdgeToSeqID || ToSeqID > m_HighestEdgeToSeqID)
	return(NULL);
if(pEdgeID != NULL && (EdgeID = *pEdgeID) > 0)
	{
	*pEdgeID = 0;
	if(EdgeID >= m_NumSeqEdges)
		return(NULL);
	pEdge = m_ppToSeqEdges[EdgeID];
	if(pEdge->TargetSeqID == ToSeqID)
		{
		*pEdgeID = EdgeID+1;
		return(pEdge);
		}
	return(NULL);
	}
EdgeID = LocateFirstEdgeToSeqID(ToSeqID);
if(pEdgeID != NULL)
	*pEdgeID = EdgeID;
return(EdgeID > 0 ? m_ppToSeqEdges[EdgeID-1] : NULL);
}


const int cMaxLinkedSeqs = 5;		// if more than this number of up/dnstream PE linked different sequences then deem as unresolvable 
typedef struct TAG_sLinkedSeq {
	int SeqID;			// this sequence is linked
	int NumSubs;		// number of subs
	int MateNumSubs;    // number of subs in linked sequence 
	int SenseCnt;       // number of links which are sense relative
	int AntisenseCnt;   // number of links which are antisense relative
} tsLinkedSeq;

// A: antisense
// S: sense
//			
//	          Ref  Mate PredSense SuccSense		
// R1)  PE1	   S	S		         S
// R2)  PE1	   S	A		         A  XXXX
// R3)  PE1	   A	S	   A	        XXXX
// R4)  PE1	   A	A	   S	
// R5)  PE2	   S	S	   S	
// R6)  PE2	   S	A	   A	         XXXX
// R7)  PE2	   A	S		        A    XXXX
// R8)  PE2	   A	A		        S

int											//  0 if none, 1 if predecessor only, 2 if successor only, 3 if both predecessor and successor identified
CScaffolder::GetPredSuccSeq(tsSeqVertex *pSeqVert)		// sequence vertex to update with predecessor and successor sequence identifiers and scores
{
int Rslt;
int LinkedSeqsIdx;
int NumLinkedSeqs;
tsLinkedSeq LinkedSeqs[cMaxLinkedSeqs];
tsLinkedSeq *pLinkedSeq;
tsLinkedSeq *pRetSeq;

int NumEdges;
tsSeqEdge *pEdge;
tsSeqEdge *pMateEdge;
tSeqID LinkedSeqID;
tEdgeID EdgeID;
tEdgeID MateEdgeID;

int ReqNumSubs;

double Diff;
double HiDiff;

if(pSeqVert == NULL)
	return(0);

pSeqVert->PredSeqID = 0;
pSeqVert->PredScore = 0;
pSeqVert->PredNumPEs = 0;
pSeqVert->SuccSeqID = 0;
pSeqVert->SuccScore = 0;
pSeqVert->SuccNumPEs = 0;


// firstly determine most probable predecessor to the current vertex
// iterate over all incomming edges to this sequence and determine most likely predecessor using rules R3, R4, R5, R6
//	          Ref  Mate PredSense		
// R3)  PE1    A    S	   A	 
// R4)  PE1	   A	A	   S	
// R5)  PE2	   S	S	   S	
// R6)  PE2	   S	A	   A	

ReqNumSubs = -1;
if((EdgeID = LocateFirstEdgeToSeqID(pSeqVert->SeqID)) != 0)			// if no edges onto this vertex sequence then can't be any linkage from any putative predecessor sequence
	{
	LinkedSeqID = 0;
	NumLinkedSeqs = 0;

	NumEdges = 0;
	do {
		LinkedSeqsIdx = 0;
		pEdge = m_ppToSeqEdges[EdgeID-1];
		if(pEdge->TargetSeqID != pSeqVert->SeqID)					// if no longer matching then must have exhusted all PEs linking to this sequence
			break;

		if((!pEdge->FlgIsPE2 && pEdge->FlgRevCpltd) ||				// only interested in matches (R3, R4, R5, R6) for PE1 being antisense or PE2 being sense in matches onto the reference sequence SeqID 
			 (pEdge->FlgIsPE2 && !pEdge->FlgRevCpltd)) 
			{
			// iterate over the PE mate edges as these will be overlapping onto putative predecessor SEs 
			MateEdgeID = LocateFirstEdgeFromSeqID(pEdge->ProbeMateSeqID);
			do
				{
				pMateEdge = &m_pSeqEdges[MateEdgeID-1];				// PE mate has the predecessor sequence identifier
				if(pEdge->ProbeMateSeqID != pMateEdge->ProbeSeqID)	
					break;

				if(ReqNumSubs == -1 || (pEdge->NumSubs + pMateEdge->NumSubs) < ReqNumSubs) // selecting for PE links with minimal number of subs required for a match
					{
					ReqNumSubs = pEdge->NumSubs + pMateEdge->NumSubs;
					NumLinkedSeqs = 0;
					}
				if((pEdge->NumSubs + pMateEdge->NumSubs) > ReqNumSubs)
					continue;

				pLinkedSeq = LinkedSeqs;								// check if this mate already known 
				for(LinkedSeqsIdx = 0; LinkedSeqsIdx < NumLinkedSeqs; LinkedSeqsIdx += 1, pLinkedSeq += 1)
					{
					if(pLinkedSeq->SeqID == pMateEdge->TargetSeqID)
						break;
					}
				if(LinkedSeqsIdx == cMaxLinkedSeqs)						// if unknown predecessor sequence but too many putative upstream already known then treat as if none upstream
					{
					NumLinkedSeqs = 0;
					break;
					}
				if(LinkedSeqsIdx == NumLinkedSeqs)						// will be equal if this predecessor sequence not already known
					{
					pLinkedSeq->SenseCnt = 0;
					pLinkedSeq->AntisenseCnt = 0;
					pLinkedSeq->SeqID = pMateEdge->TargetSeqID;
					NumLinkedSeqs += 1;
					}

				pLinkedSeq->MateNumSubs = pMateEdge->NumSubs;
				pLinkedSeq->NumSubs = pEdge->NumSubs;               
				if((!pEdge->FlgIsPE2 && pMateEdge->FlgRevCpltd) ||	 
					(pEdge->FlgIsPE2 && !pMateEdge->FlgRevCpltd))   
					pLinkedSeq->SenseCnt += 1;                        // R4 or R5, predecessor is relative sense
				else                                                
					pLinkedSeq->AntisenseCnt += 1;                    // R3 or R6, predecessor is relative antisense
				}
			while(++MateEdgeID <= m_NumSeqEdges);
			}
		}
	while(LinkedSeqsIdx < cMaxLinkedSeqs && ++EdgeID <= m_NumSeqEdges);
	if(NumLinkedSeqs != 0)			// if linking to at least one predecessor (5' upstream) sequence linking, choose the predecessor with the most PE linkage support
		{
		pLinkedSeq = LinkedSeqs;
		HiDiff = 0.0;
		pRetSeq = NULL;
		for(LinkedSeqsIdx = 0; LinkedSeqsIdx < NumLinkedSeqs; LinkedSeqsIdx += 1, pLinkedSeq += 1)
			{
			Diff = (double)abs(pLinkedSeq->AntisenseCnt - pLinkedSeq->SenseCnt)/(double)(pLinkedSeq->AntisenseCnt + pLinkedSeq->SenseCnt);
			if(Diff >= HiDiff && Diff <= (HiDiff + 0.1)) // if very little difference from previously accepted HiDiff then treat as if same and thus need to look for a clear winner
				pRetSeq = NULL;
			else
				if(Diff > 0.75 && Diff > HiDiff)	// needs to be a substantial differential between sense and antisense counts before putative acceptance
					{
					HiDiff = Diff;
					pRetSeq = pLinkedSeq;
					}
			}

		if(pRetSeq != NULL)
			{
			pSeqVert->PredSeqID = pRetSeq->SeqID;
			pSeqVert->PredNumPEs = max(pRetSeq->SenseCnt,pRetSeq->AntisenseCnt);
			pSeqVert->PredScore = max(1,10 - ReqNumSubs) * pSeqVert->PredNumPEs;
			if(pRetSeq->SenseCnt < pRetSeq->AntisenseCnt)
				pSeqVert->PredScore *= -1;
			}
		}
	}

// now determine most probable successor to the current vertex
// iterate over all incomming edges to this sequence and determine most likely successor using rules R1, R2, R7, R8
//	          Ref  Mate  SuccSense		
// R1)  PE1	   S	S		S
// R2)  PE1	   S	A		A
// R7)  PE2	   A	S		A
// R8)  PE2	   A	A		S

ReqNumSubs = -1;
if((EdgeID = LocateFirstEdgeToSeqID(pSeqVert->SeqID)) != 0)	// if none then no successor sequence
	{
	LinkedSeqID = 0;
	NumLinkedSeqs = 0;
	NumEdges = 0;
	do {
		LinkedSeqsIdx = 0;
		pEdge = m_ppToSeqEdges[EdgeID-1];
		if(pEdge->TargetSeqID != pSeqVert->SeqID)					// if no longer matching then exhusted PEs linking to this sequence
			break;

		if((!pEdge->FlgIsPE2 && !pEdge->FlgRevCpltd) ||	// interested in matches (R1, R2, R7, R8) for PE1 being sense or PE2 being antisense in matches onto the reference sequence SeqID 
			 (pEdge->FlgIsPE2 && pEdge->FlgRevCpltd)) 
			{
			// iterate over the PE mate edges
			MateEdgeID = LocateFirstEdgeFromSeqID(pEdge->ProbeMateSeqID);
			do
				{
				pMateEdge = &m_pSeqEdges[MateEdgeID-1];			// PE mate has the successor sequence identifier
				if(pEdge->ProbeMateSeqID != pMateEdge->ProbeSeqID)	
					break;

				if(ReqNumSubs == -1 || (pEdge->NumSubs + pMateEdge->NumSubs) < ReqNumSubs) // selecting for PE links with minimal number of subs required for a match
					{
					ReqNumSubs = pEdge->NumSubs + pMateEdge->NumSubs;
					NumLinkedSeqs = 0;
					}
				if((pEdge->NumSubs + pMateEdge->NumSubs) > ReqNumSubs)
					continue;

				pLinkedSeq = LinkedSeqs;						// check if this mate already known 
				for(LinkedSeqsIdx = 0; LinkedSeqsIdx < NumLinkedSeqs; LinkedSeqsIdx += 1, pLinkedSeq += 1)
					{
					if(pLinkedSeq->SeqID == pMateEdge->TargetSeqID)
						break;
					}
				if(LinkedSeqsIdx == cMaxLinkedSeqs)						// if unknown successor sequence but too many putative upstream already known then treat as if none upstream
					{
					NumLinkedSeqs = 0;
					break;
					}
				if(LinkedSeqsIdx == NumLinkedSeqs)						// will be equal if this successor sequence not already known
					{
					pLinkedSeq->SenseCnt = 0;
					pLinkedSeq->AntisenseCnt = 0;
					pLinkedSeq->SeqID = pMateEdge->TargetSeqID;
					NumLinkedSeqs += 1;
					}

				pLinkedSeq->MateNumSubs = pMateEdge->NumSubs;
				pLinkedSeq->NumSubs = pEdge->NumSubs;

				if((!pEdge->FlgIsPE2 && !pMateEdge->FlgRevCpltd) ||     // R1
					(pEdge->FlgIsPE2 && pMateEdge->FlgRevCpltd))        // R8
					pLinkedSeq->SenseCnt += 1;							// R1 or R8, successor is relative sense
				else
					pLinkedSeq->AntisenseCnt += 1;                      // R2 or R7, successor is relative antisense
				}
			while(++MateEdgeID <= m_NumSeqEdges);
			}
		}
	while(LinkedSeqsIdx < cMaxLinkedSeqs && ++EdgeID <= m_NumSeqEdges);
	if(NumLinkedSeqs > 0)  // if linking to at least one successor (3' downstream) sequence linking, choose the successor with the most PE linkage support
		{
		pLinkedSeq = LinkedSeqs;
		HiDiff = 0.0; 
		pRetSeq = NULL;
		for(LinkedSeqsIdx = 0; LinkedSeqsIdx < NumLinkedSeqs; LinkedSeqsIdx += 1, pLinkedSeq += 1)
			{
			Diff = (double)abs(pLinkedSeq->AntisenseCnt - pLinkedSeq->SenseCnt)/(double)(pLinkedSeq->AntisenseCnt + pLinkedSeq->SenseCnt);
			if(Diff >= HiDiff && Diff <= (HiDiff + 0.1)) // if very little difference from previously accepted HiDiff then treat as if same and thus need to look for a clear winner
				pRetSeq = NULL;
			else
				if(Diff > 0.75 && Diff > HiDiff)	// needs to be a substantial differential between sense and antisense counts before putative acceptance
					{
					HiDiff = Diff;
					pRetSeq = pLinkedSeq;
					}
			}
		if(pRetSeq != NULL)
			{
			pSeqVert->SuccSeqID = pRetSeq->SeqID;
			pSeqVert->SuccNumPEs = max(pRetSeq->SenseCnt,pRetSeq->AntisenseCnt);
			pSeqVert->SuccScore = max(1,10 - ReqNumSubs) * pSeqVert->SuccNumPEs;
			if(pRetSeq->SenseCnt < pRetSeq->AntisenseCnt)
				pSeqVert->SuccScore *= -1;
			}
		}
	}

Rslt = pSeqVert->PredSeqID ? 0x01 : 0;
Rslt |= pSeqVert->SuccSeqID ? 0x02 : 0;

return(Rslt);
}


int				// assuming that SeqID is correctly orientated as sense then the number of edges (in) from other sequences onto SeqID
CScaffolder::NumEdgesIn(tSeqID SeqID)		// edges into this sequence
{
int NumEdges;
tsSeqEdge *pEdge;
tEdgeID EdgeID;

if(m_NumSeqEdges == 0 || m_pSeqEdges == NULL || m_ppToSeqEdges == NULL || SeqID < m_LowestEdgeToSeqID || SeqID > m_HighestEdgeToSeqID)
	return(0);

if((EdgeID = LocateFirstEdgeToSeqID(SeqID)) == 0)
	return(0);

NumEdges = 0;
do {
	pEdge = m_ppToSeqEdges[EdgeID-1];
	if(pEdge->TargetSeqID != SeqID)
		break;
	if((pEdge->FlgIsPE2 && !pEdge->FlgRevCpltd) || (!pEdge->FlgIsPE2 && pEdge->FlgRevCpltd))
		NumEdges += 1;
	}
while(++EdgeID <= m_NumSeqEdges);
return(NumEdges);
}

int				// number of edges (out) from SeqID onto other sequences
CScaffolder::NumEdgesOut(tSeqID SeqID)      // edges out from this sequence
{
int NumEdges;
tsSeqEdge *pEdge;
tEdgeID EdgeID;

if(m_NumSeqEdges == 0 || m_pSeqEdges == NULL || m_ppToSeqEdges == NULL || SeqID < m_LowestEdgeToSeqID || SeqID > m_HighestEdgeToSeqID)
	return(0);

if((EdgeID = LocateFirstEdgeToSeqID(SeqID)) == 0)
	return(0);

NumEdges = 0;
do {
	pEdge = m_ppToSeqEdges[EdgeID-1];
	if(pEdge->TargetSeqID != SeqID)
		break;
	if((pEdge->FlgIsPE2 && pEdge->FlgRevCpltd) || (!pEdge->FlgIsPE2 && !pEdge->FlgRevCpltd))
		NumEdges += 1;
	}
while(++EdgeID <= m_NumSeqEdges);
return(NumEdges);
}

tEdgeID			// 0 if unable to locate any with matching FromSeqID, otherwise m_pSeqEdges[index-1] of lowest ordered sequence edge				 
CScaffolder::LocateFirstEdgeFromSeqID(tSeqID FromSeqID)	// find lowest ordered matching sequence edge with matching FromSeqID 
{
tsSeqEdge *pEl2;
int CmpRslt;
INT64 Mark;
INT64 TargPsn;
INT64 NodeLo = 0;
INT64 NodeHi = m_NumSeqEdges-1;
if(m_NumSeqEdges == 0 || m_pSeqEdges == NULL)
	return(0);

do {
	TargPsn = (NodeLo + NodeHi) / 2L;
	pEl2 = &m_pSeqEdges[TargPsn];

	if(FromSeqID > pEl2->ProbeSeqID)
		CmpRslt = 1;
	else
		if(FromSeqID < pEl2->ProbeSeqID)
			CmpRslt = -1;
		else
			CmpRslt = 0;

	if(!CmpRslt)	// if a match then may not be the lowest indexed match
		{
		if(TargPsn == 0 || NodeLo == TargPsn) // check if already lowest
			return(TargPsn+1);
		// iterate until lowest located
		while(1) {
			if(CmpRslt == 0)
				{
				Mark = TargPsn;
				if(Mark == 0)
					return(Mark+1);
				NodeHi = TargPsn - 1;
				}	
			TargPsn = ((UINT64)NodeLo + NodeHi) / 2L;

			pEl2 = &m_pSeqEdges[TargPsn];
			if(FromSeqID > pEl2->ProbeSeqID)
				CmpRslt = 1;
			else
				if(FromSeqID < pEl2->ProbeSeqID)
					CmpRslt = -1;
				else
					{
					CmpRslt = 0;
					continue;
					}
			NodeLo = TargPsn + 1;
			if(NodeLo == Mark)
				return(Mark+1);
			}
		}

	if(CmpRslt < 0)
		{
		if(TargPsn == 0)
			break;
		NodeHi = TargPsn - 1;
		}
	else
		NodeLo = TargPsn+1;
	}
while(NodeHi >= NodeLo);
return(0);	// unable to locate any instance of FromSeqID
}

tEdgeID			// 0 if unable to locate any with matching ToSeqID, otherwise m_ppToSeqEdges[index-1] of lowest ordered sequence edge				 
CScaffolder::LocateFirstEdgeToSeqID(tSeqID ToSeqID)	// find lowest ordered matching sequence edge with matching ToSeqID 
{
tsSeqEdge *pEl2;
int CmpRslt;
INT64 Mark;
INT64 TargPsn;
INT64 NodeLo = 0;
INT64 NodeHi = m_NumSeqEdges-1;
if(m_NumSeqEdges == 0 || m_pSeqEdges == NULL || m_ppToSeqEdges == NULL)
	return(0);

do {
	TargPsn = (NodeLo + NodeHi) / 2L;
	pEl2 = m_ppToSeqEdges[TargPsn];

	if(ToSeqID > pEl2->TargetSeqID)
		CmpRslt = 1;
	else
		if(ToSeqID < pEl2->TargetSeqID)
			CmpRslt = -1;
		else
			CmpRslt = 0;


	if(!CmpRslt)	// if a match then may not be the lowest indexed match
		{
		if(TargPsn == 0 || NodeLo == TargPsn) // check if already lowest
			return(TargPsn+1);
		// iterate until lowest located
		while(1) {
			if(CmpRslt == 0)
				{
				Mark = TargPsn;
				if(Mark == 0)
					return(Mark+1);
				NodeHi = TargPsn - 1;
				}	
			TargPsn = ((UINT64)NodeLo + NodeHi) / 2L;

			pEl2 = m_ppToSeqEdges[TargPsn];
			if(ToSeqID > pEl2->TargetSeqID)
				CmpRslt = 1;
			else
				if(ToSeqID < pEl2->TargetSeqID)
					CmpRslt = -1;
				else
					{
					CmpRslt = 0;
					continue;
					}
			NodeLo = TargPsn + 1;
			if(NodeLo == Mark)
				return(Mark+1);
			}
		}

	if(CmpRslt < 0)
		{
		if(TargPsn == 0)
			break;
		NodeHi = TargPsn - 1;
		}
	else
		NodeLo = TargPsn+1;
	}
while(NodeHi >= NodeLo);
return(0);	// unable to locate any instance of ToSeqID
}


tsSeqVertex *									// NULL if unable to locate any with matching SeqID				 
CScaffolder::LocateVerticesSeqID(tSeqID SeqID)	// find vertex with this SeqID 
{
tsSeqVertex *pEl2;
int CmpRslt;
INT64 TargPsn;
INT64 NodeLo = 0;
INT64 NodeHi = m_NumSeqVertices-1;
if(m_NumSeqVertices == 0 || m_pSeqVertices == NULL)
	return(0);

do {
	TargPsn = (NodeLo + NodeHi) / 2L;
	pEl2 = &m_pSeqVertices[TargPsn];

	if(SeqID > pEl2->SeqID)
		CmpRslt = 1;
	else
		if(SeqID < pEl2->SeqID)
			CmpRslt = -1;
		else
			CmpRslt = 0;

	if(!CmpRslt)	
		return(pEl2);

	if(CmpRslt < 0)
		{
		if(TargPsn == 0)
			break;
		NodeHi = TargPsn - 1;
		}
	else
		NodeLo = TargPsn+1;
	}
while(NodeHi >= NodeLo);
return(NULL);	
}






int // function for sorting vertices on their SeqID ascending
CScaffolder::SortVerticesSeqID(const void *arg1, const void *arg2)
{
tsSeqVertex *pV1 = (tsSeqVertex *)arg1;
tsSeqVertex *pV2 = (tsSeqVertex *)arg2;
if(pV1->SeqID < pV2->SeqID)
	return(-1);
else
	if(pV1->SeqID > pV2->SeqID)
		return(1);
return(0);
}

int  // function for sorting sequence edges on FromSeqID ascending, OverlapLen descending 
CScaffolder::SortSeqEdgesFromSeqID(const void *arg1, const void *arg2)
{
tsSeqEdge *pSeqEdge1 = (tsSeqEdge *)arg1;
tsSeqEdge *pSeqEdge2 = (tsSeqEdge *)arg2;

if(pSeqEdge1->ProbeSeqID < pSeqEdge2->ProbeSeqID)
	return(-1);
else
	if(pSeqEdge1->ProbeSeqID > pSeqEdge2->ProbeSeqID)
		return(1);
if(pSeqEdge1->zOverlapLen > pSeqEdge2->zOverlapLen)
	return(-1);
else
	if(pSeqEdge1->zOverlapLen < pSeqEdge2->zOverlapLen)
		return(1);
return(0);
}

int  // function for sorting sequence edges on FromSeqID ascending, TargetSeqID ascending 
CScaffolder::SortSeqEdgesFromToSeqID(const void *arg1, const void *arg2)
{
tsSeqEdge *pSeqEdge1 = (tsSeqEdge *)arg1;
tsSeqEdge *pSeqEdge2 = (tsSeqEdge *)arg2;

if(pSeqEdge1->ProbeSeqID < pSeqEdge2->ProbeSeqID)
	return(-1);
else
	if(pSeqEdge1->ProbeSeqID > pSeqEdge2->ProbeSeqID)
		return(1);

if(pSeqEdge1->TargetSeqID < pSeqEdge2->TargetSeqID)
	return(-1);
else
	if(pSeqEdge1->TargetSeqID > pSeqEdge2->TargetSeqID)
		return(1);

return(0);
}

int  // function for sorting m_ppToSeqEdges[] sequence edges on ToSeqID ascending, OverlapLen descending 
CScaffolder::SortSeqEdgesToSeqID(const void *arg1, const void *arg2)
{
tsSeqEdge *pSeqEdge1 = *(tsSeqEdge **)arg1;
tsSeqEdge *pSeqEdge2 = *(tsSeqEdge **)arg2;

if(pSeqEdge1->TargetSeqID < pSeqEdge2->TargetSeqID)
	return(-1);
else
	if(pSeqEdge1->TargetSeqID > pSeqEdge2->TargetSeqID)
		return(1);
if(pSeqEdge1->zOverlapLen > pSeqEdge2->zOverlapLen)
	return(-1);
else
	if(pSeqEdge1->zOverlapLen < pSeqEdge2->zOverlapLen)
		return(1);
return(0);
}