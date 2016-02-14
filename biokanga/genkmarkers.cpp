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

// locmarkers.cpp : Defines the entry point for the console application.
// Uses a pregenerated suffix array built over all concatenated reads, psueudo-chromosome for each species/cultivar, to identify 
// potential marker sequences of user specified length in a specified cultivar. These are k-mers which are unique to a species/cultivar and are 
// at least H Hammings away from any other species k-mer.
// These identified sequences are potential markers for the specified species.
//
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

#include "biokanga.h"
#include "./LocKMers.h"

int
LocMarkers(etPMode PMode,				// processing mode - defaults to 0
		  int KMerLen,					// this length K-mers
	  	  int PrefixLen,				// inter-cultivar shared prefix length
		  int SuffixLen,				// cultivar specific suffix length
		  int MinWithPrefix,			// minimum number of cultivars required to have the shared prefix
		  int MinHamming,				// must be at least this Hamming away from any other K-mer in other cultivars
		  char *pszCultivarName,		// targeted cultivar name for which K-mer markers are required
		  int NumPartialCultivars,		// there are this many pseudo chromosomes for targeted cultivar for which K-mer markers are required
		  char *ppszPartialCultivars[],	// pseudo chromosome names which identify targeted cultivar
		  char *pszSfxPseudoGenome,		// contains pregenerated suffix over psuedochromosomes for each cultivar
		  char *pszMarkerFile,			// output potential markers to this file
		  char *pszMarkerReadsFile,		// optionally output reads containing potential markers to this file
		  int NumThreads);				// max number of threads allowed

#ifdef _WIN32
int genkmarkers(int argc, char* argv[])
{
// determine my process name
_splitpath(argv[0],NULL,NULL,gszProcName,NULL);
#else
int
genkmarkers(int argc, char** argv)
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

int PMode;					// processing mode
int KMerLen;				// this length K-mers
int PrefixLen;				// inter-cultivar shared prefix length
int SuffixLen;				// cultivar specific suffix length
int MinWithPrefix;			// minimum number of cultivars required to have the shared prefix
int MinHamming;				// must be at least this Hamming away from any other K-mer in other species
char szCultivarName[cMaxDatasetSpeciesChrom+1];		// cultivar name
char szPartialCultivarsList[(cMaxDatasetSpeciesChrom + 10) * cMaxTargCultivarChroms];		// individual species parsed from this comma/tab/space separated list
int NumPartialCultivars;	// there are this many pseudo chromosomes for targeted cultivar for which K-mer markers are required
char *pszPartialCultivars[cMaxTargCultivarChroms+1];	// pseudo chromosome names which identify targeted cultivar
char szSfxPseudoGenome[_MAX_PATH];		// contains assembly + suffix array over all psuedo-chromosomes for all cultivars
char szMarkerFile[_MAX_PATH];			// output potential markers to this file
char szMarkerReadsFile[_MAX_PATH];		// output reads containing potential markers to this file


char szSQLiteDatabase[_MAX_PATH];	// results summaries to this SQLite file
char szExperimentName[cMaxDatasetSpeciesChrom+1];			// experiment name
char szExperimentDescr[1000];		// describes experiment

// command line args
struct arg_lit  *help    = arg_lit0("h","help",                 "Print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"Print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to screen and logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"Diagnostics log file");

struct arg_int *pmode = arg_int0("m","mode","<int>",			"Processing mode : 0 - default with K-mer extension, 1 - no K-mer extension, 2 - inter-cultivar shared prefix sequences ");
struct arg_int *kmerlen = arg_int0("k","kmer","<int>",			"Cultivar specific K-mers of this length (default 50, range 25..100)");

struct arg_int *prefixlen = arg_int0("p","prefixlen","<int>",	"Cultivar specific K-mers to contain inter-cultivar shared prefix sequences of this length (Mode 2 only");
struct arg_int *minwithprefix = arg_int0("s","minshared","<int>","Inter-cultivar shared prefix sequences must be present in this many cultivars (Mode 2 only, default all)");

struct arg_int *minhamming = arg_int0("K","minhamming","<int>", "Minimum Hamming separation distance in other non-target cultivars (default 2, range 1..5)");
struct arg_str *cultivar = arg_str1("c","cultivar","<str>",		"Cultivar name to associate with identified marker K-mers");
struct arg_str *chromnames = arg_str1("C","chromnames","<str>",	"Comma/space separated list of pseudo-chrom names specific to cultivar for which markers are required");
struct arg_file *infile = arg_file1("i","in","<file>",		    "Use this suffix indexed pseudo-chromosomes file");
struct arg_file *outfile = arg_file1("o","markers","<file>",		"Output accepted marker K-mer sequences to this multifasta file");
struct arg_file *outreadsfile = arg_file0("O","markerreads","<file>",	"Output reads containing accepted marker K-mers to this multifasta file");
struct arg_int *numthreads = arg_int0("T","threads","<int>",		"number of processing threads 0..128 (defaults to 0 which sets threads to number of CPU cores)");

struct arg_file *summrslts = arg_file0("q","sumrslts","<file>",		"Output results summary to this SQLite3 database file");
struct arg_str *experimentname = arg_str0("w","experimentname","<str>",		"experiment name SQLite3 database file");
struct arg_str *experimentdescr = arg_str0("W","experimentdescr","<str>",	"experiment description SQLite3 database file");
struct arg_end *end = arg_end(200);

void *argtable[] = {help,version,FileLogLevel,LogFile,
					summrslts,experimentname,experimentdescr,
	                pmode,kmerlen,prefixlen,minwithprefix,minhamming,cultivar,chromnames,infile,outfile,outreadsfile,
					numthreads,
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

// show user current resource limits
#ifndef _WIN32
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Resources: %s",CUtility::ReportResourceLimits());
#endif

#ifdef _WIN32
	SYSTEM_INFO SystemInfo;
	GetSystemInfo(&SystemInfo);
	NumberOfProcessors = SystemInfo.dwNumberOfProcessors;
#else
	NumberOfProcessors = sysconf(_SC_NPROCESSORS_CONF);
#endif
	int MaxAllowedThreads = min(cMaxWorkerThreads,NumberOfProcessors);	// limit to be at most cMaxWorkerThreads
	if((NumThreads = numthreads->count ? numthreads->ival[0] : MaxAllowedThreads)==0)
		NumThreads = MaxAllowedThreads;
	if(NumThreads < 0 || NumThreads > MaxAllowedThreads)
		{
		gDiagnostics.DiagOut(eDLWarn,gszProcName,"Warning: Number of threads '-T%d' specified was outside of range %d..%d",NumThreads,1,MaxAllowedThreads);
		gDiagnostics.DiagOut(eDLWarn,gszProcName,"Warning: Defaulting number of threads to %d",MaxAllowedThreads);
		NumThreads = MaxAllowedThreads;
		}


	PMode = (etPMode)(pmode->count ? pmode->ival[0] : ePMExtdKMers);
	if(PMode < ePMExtdKMers || PMode >= ePMplaceholder)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Processing mode '-m%d' specified outside of range %d..%d\n",PMode,ePMExtdKMers,(int)ePMplaceholder-1);
		exit(1);
		}

	KMerLen = kmerlen->count ? kmerlen->ival[0] : cDfltKMerLen;
	if(KMerLen < cMinKMerLen || KMerLen > cMaxKMerLen)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: K-mer core length '-k%d' must be in range %d..%d",KMerLen,cMinKMerLen,cMaxKMerLen);
		return(1);
		}

	if(PMode == 2)
		{
		PrefixLen = prefixlen->count ? prefixlen->ival[0] : KMerLen/2;
		if(PrefixLen < cMinKMerLen/2 || PrefixLen > KMerLen-(cMinKMerLen/2))
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Prefix length '-p%d' must be in range %d..%d",KMerLen,cMinKMerLen/2,KMerLen-(cMinKMerLen/2));
			return(1);
			}
		SuffixLen = KMerLen - PrefixLen;
		MinWithPrefix = minwithprefix->count ? minwithprefix->ival[0] : 0;
		if(MinWithPrefix != 0 && MinWithPrefix < 1)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Minimum cultivars sharing prefix sequence '-s%d' must be either 0 (all) or at least 1",MinWithPrefix);
			return(1);
			}
		}

	MinHamming = minhamming->count ? minhamming->ival[0] : cDfltHamming;
	if(MinHamming < cMinHamming || MinHamming > cMaxHamming)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Minimum Hamming separation '-K%d' must be in range %d..%d",MinHamming,cMinHamming,cMaxHamming);
		return(1);
		}

	strncpy(szCultivarName,cultivar->sval[0],cMaxDatasetSpeciesChrom);
	szCultivarName[cMaxDatasetSpeciesChrom]= '\0';
	CUtility::TrimQuotedWhitespcExtd(szCultivarName);
	if(strlen(szCultivarName) < 1)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Expected cultivar name '-c<name>' is empty");
		return(1);
		}

	strncpy(szPartialCultivarsList,chromnames->sval[0],sizeof(szPartialCultivarsList));
	szPartialCultivarsList[sizeof(szPartialCultivarsList)-1] = '\0';
	CUtility::TrimQuotedWhitespcExtd(szPartialCultivarsList);
	CUtility::ReduceWhitespace(szPartialCultivarsList);
	char *pChr = szPartialCultivarsList;
	char *pStartChr;
	char Chr;
	int CurSpeciesLen;
	NumPartialCultivars=0;
	CurSpeciesLen = 0;
	pStartChr = pChr;
	while((Chr = *pChr++) != '\0')
		{
		if(Chr == ' ' || Chr == '\t' || Chr == ',')	// treat any of these as delimiters
			{
			pChr[-1] = '\0';
			if(CurSpeciesLen != 0)
				{
				pszPartialCultivars[NumPartialCultivars++] = pStartChr;
				CurSpeciesLen = 0;
				}
			pStartChr = pChr;
			continue;
			}
		CurSpeciesLen += 1;
		}
	if(CurSpeciesLen)
		pszPartialCultivars[NumPartialCultivars++] = pStartChr;

	if(!NumPartialCultivars)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Expected at least one ('-C<names>') targeted cultivar chromosme name");
		return(1);
		}

	strcpy(szSfxPseudoGenome,infile->filename[0]);
	CUtility::TrimQuotedWhitespcExtd(szSfxPseudoGenome);
	if(strlen(szSfxPseudoGenome) < 1)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Expected input pseudo-genome suffix array filename '-i<name>' is empty");
		return(1);
		}

	strcpy(szMarkerFile,outfile->filename[0]);
	CUtility::TrimQuotedWhitespcExtd(szMarkerFile);
	if(strlen(szMarkerFile) < 1)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Expected marker file to generate filename '-o<name>' is empty");
		return(1);
		}

	if(outreadsfile->count)
		{
		strcpy(szMarkerReadsFile,outreadsfile->filename[0]);
		CUtility::TrimQuotedWhitespcExtd(szMarkerReadsFile);
		if(strlen(szMarkerReadsFile) < 1)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Reads containing markers filename '-O<name>' is empty");
			return(1);
			}
		}
	else
		szMarkerReadsFile[0] = '\0';

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing parameters:");
	const char *pszDescr;
	switch(PMode) {
		case ePMExtdKMers:
			pszDescr = "Extended K-mer markers";
			break;
		case ePMNoExtdKMers:
			pszDescr = "Non-extended K-mer markers";
			break;
		case ePMPrefixKMers:
			pszDescr = "K-mers to share prefix sequence with other cultivars";
			break;
		}

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Processing mode is : '%s'",pszDescr);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Targeted cultivar name : '%s'",szCultivarName);

	for(int Idx = 0; Idx < NumPartialCultivars; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Targeted cultivar chromosome name (%d) : '%s'", Idx + 1,pszPartialCultivars[Idx]);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Core K-mer length : %d",KMerLen);
	if(PMode == ePMPrefixKMers)
		{
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Inter-cultivar shared prefix sequence length : %d",PrefixLen);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Cultivar specific suffix sequence length : %d",SuffixLen);
		if(MinWithPrefix)
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Min number cultivars sharing prefix : %d",MinWithPrefix);
		else
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Min number cultivars sharing prefix : 'All'");
		}
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Minimum Hamming separation : %d",MinHamming);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Input indexed pseudo-genome file: '%s'",szSfxPseudoGenome);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Write marker K-mers to file: '%s'",szMarkerFile);
	if(szMarkerReadsFile[0] != '\0')
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Write marker containing reads to file: '%s'",szMarkerReadsFile);

	if(szExperimentName[0] != '\0')
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"This processing reference: %s",szExperimentName);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Number of processing threads: %d",NumThreads);

	if(gExperimentID > 0)
		{
		int ParamID;

		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szLogFile),"log",szLogFile);

		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,sizeof(PMode),"mode",&PMode);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,sizeof(KMerLen),"kmer",&KMerLen);

		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,sizeof(PrefixLen),"prefixlen",&PrefixLen);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,sizeof(SuffixLen),"suffixlen",&SuffixLen);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,sizeof(MinWithPrefix),"minwithprefix",&MinWithPrefix);

		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,sizeof(MinHamming),"minhamming",&MinHamming);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,sizeof(NumThreads),"threads",&NumThreads);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,sizeof(NumberOfProcessors),"cpus",&NumberOfProcessors);

		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szCultivarName),"cultivar",szCultivarName);

		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,sizeof(NumPartialCultivars),"NumPartialCultivars",&NumPartialCultivars);
		for(int Idx = 0; Idx < NumPartialCultivars; Idx++)
			ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(pszPartialCultivars[Idx]),"chromnames",pszPartialCultivars[Idx]);

		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szSfxPseudoGenome),"in",szSfxPseudoGenome);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szMarkerFile),"markers",szMarkerFile);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szMarkerReadsFile),"markerreads",szMarkerReadsFile);

		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szSQLiteDatabase),"sumrslts",szSQLiteDatabase);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szExperimentName),"experimentname",szExperimentName);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szExperimentDescr),"experimentdescr",szExperimentDescr);
		}


#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	gStopWatch.Start();
	Rslt = LocMarkers((etPMode)PMode,KMerLen,PrefixLen,SuffixLen,MinWithPrefix,MinHamming,szCultivarName,NumPartialCultivars,pszPartialCultivars,szSfxPseudoGenome,szMarkerFile,szMarkerReadsFile,NumThreads);
	Rslt = Rslt >=0 ? 0 : 1;
	if(gExperimentID > 0)
		{
		if(gProcessingID)
			gSQLiteSummaries.EndProcessing(gProcessingID,Rslt);
		gSQLiteSummaries.EndExperiment(gExperimentID);
		}
	gStopWatch.Stop();
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
LocMarkers(etPMode PMode,				// processing mode - defaults to 0
		  int KMerLen,					// this length K-mers
	  	  int PrefixLen,				// inter-cultivar shared prefix length
		  int SuffixLen,				// cultivar specific suffix length
		  int MinWithPrefix,			// minimum number of cultivars required to have the shared prefix
		  int MinHamming,				// must be at least this Hamming away from any other K-mer in other cultivars
		  char *pszCultivarName,		// targeted cultivar name for which K-mer markers are required
		  int NumPartialCultivars,		// there are this many pseudo chromosomes for targeted cultivar for which K-mer markers are required
		  char *ppszPartialCultivars[],	// pseudo chromosome names which identify targeted cultivar
		  char *pszSfxPseudoGenome,		// contains pregenerated suffix over psuedochromosomes for each cultivar
		  char *pszMarkerFile,			// output potential markers to this file
		  char *pszMarkerReadsFile,		// optionally output reads containing potential markers to this file
		  int NumThreads)				// max number of threads allowed
{
int Rslt;
CLocKMers Markers;

Rslt = Markers.LocKMers(PMode,KMerLen,PrefixLen,SuffixLen,MinHamming,MinWithPrefix,pszCultivarName,NumPartialCultivars,ppszPartialCultivars,pszSfxPseudoGenome,pszMarkerFile,pszMarkerReadsFile,NumThreads);

Markers.Reset();
return(Rslt);
}