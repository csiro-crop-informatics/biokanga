// locmarkers.cpp : Defines the entry point for the console application.
// Uses a pregenerated suffix array built over all concatenated reads, psueudo-chromosome for each species/cultivar, to identify 
// potential marker sequences of user specified length in a specified cultivar. These are k-mers which are unique to a species/cultivar and are 
// at least H Hammings away from any other species k-mer.
// These identified sequences are potential markers for the specified species.
// 0.9.4		Option to report reads containing markers

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

#include "./LocKMers.h"

const char *cpszProgVer = "0.9.5";		// increment with each release
const char *cpszProcOverview = "Locate alignment-less markers";

CStopWatch gStopWatch;
CDiagnostics gDiagnostics;				// for writing diagnostics messages to log file

char gszProcName[_MAX_FNAME];			// process name

int
LocMarkers(etPMode Mode,				// processing mode - defaults to 0
		  int KMerLen,					// this length K-mers
		  int MinHamming,				// must be at least this Hamming away from any other K-mer in other cultivars
		  char *pszCultivarName,		// targeted cultivar name for which K-mer markers are required
		  int NumPartialCultivars,		// there are this many pseudo chromosomes for targeted cultivar for which K-mer markers are required
		  char *pszPartialCultivars[],	// pseudo chromosome names which identify targeted cultivar
		  char *pszSfxPseudoGenome,		// contains pregenerated suffix over psuedochromosomes for each cultivar
		  char *pszMarkerFile,			// output potential markers to this file
		  char *pszMarkerReadsFile,		// optionally output reads containing potential markers to this file
		  int NumThreads);				// max number of threads allowed

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
	return _T("locmarkers");
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

int NumberOfProcessors;		// number of installed CPUs
int NumThreads;				// number of threads (0 defaults to number of CPUs)

etPMode PMode;				// processing mode
int KMerLen;				// this length K-mers
int MinHamming;				// must be at least this Hamming away from any other K-mer in other species
char szCultivarName[cMaxDatasetSpeciesChrom+1];		// cultivar name
char szPartialCultivarsList[(cMaxDatasetSpeciesChrom + 10) * cMaxTargCultivarChroms];		// individual species parsed from this comma/tab/space separated list
int NumPartialCultivars;	// there are this many pseudo chromosomes for targeted cultivar for which K-mer markers are required
char *pszPartialCultivars[cMaxTargCultivarChroms+1];	// pseudo chromosome names which identify targeted cultivar
char szSfxPseudoGenome[_MAX_PATH];		// contains assembly + suffix array over all psuedo-chromosomes for all cultivars
char szMarkerFile[_MAX_PATH];			// output potential markers to this file
char szMarkerReadsFile[_MAX_PATH];		// output reads containing potential markers to this file

// command line args
struct arg_lit  *help    = arg_lit0("h","help",                 "Print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"Print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to screen and logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"Diagnostics log file");

struct arg_int *pmode = arg_int0("m","mode","<int>",			"Processing mode : 0 - default with K-mer extension, 1 - no K-mer extension ");
struct arg_int *kmerlen = arg_int0("k","kmer","<int>",			"Marker K-mers of this core length (default 50, range 25..100)");
struct arg_int *minhamming = arg_int0("K","minhamming","<int>", "Minimum Hamming separation distance in other non-target cultivars (default 2, range 1..5)");
struct arg_str *cultivar = arg_str1("c","cultivar","<str>",		"Cultivar name to associate with identified marker K-mers");
struct arg_str *chromnames = arg_str1("C","chromnames","<str>",	"Comma/space separated list of pseudo-chrom names specific to cultivar for which markers are required");
struct arg_file *infile = arg_file1("i","in","<file>",		    "Use this suffix indexed pseudo-chromosomes file");
struct arg_file *outfile = arg_file1("o","markers","<file>",		"Output accepted marker K-mer sequences to this multifasta file");
struct arg_file *outreadsfile = arg_file0("O","markerreads","<file>",	"Output reads containing accepted marker K-mers to this multifasta file");
struct arg_int *numthreads = arg_int0("T","threads","<int>",		"number of processing threads 0..128 (defaults to 0 which sets threads to number of CPU cores)");

struct arg_end *end = arg_end(80);

void *argtable[] = {help,version,FileLogLevel,LogFile,
	                pmode,kmerlen,minhamming,cultivar,chromnames,infile,outfile,numthreads,outreadsfile,
					end};

char **pAllArgs;
int argerrors;
argerrors = CUtility::arg_parsefromfile(argc,(char **)argv,&pAllArgs);
if(argerrors >= 0)
	argerrors = arg_parse(argerrors,pAllArgs,argtable);

/* special case: '--help' takes precedence over error reporting */
if (help->count > 0)
        {
		printf("\n%s Locate putative K-mer marker sequences, Version %s\nOptions ---\n", gszProcName,cpszProgVer);
        arg_print_syntax(stdout,argtable,"\n");
        arg_print_glossary(stdout,argtable,"  %-25s %s\n");
		printf("\nNote: Parameters can be entered into a parameter file, one parameter per line.");
		printf("\n      To invoke this parameter file then precede it's name with '@'");
		printf("\n      e.g. %s @myparams.txt\n",gszProcName);
		printf("\nPlease report any issues regarding usage of %s at https://github.com/csiro-crop-informatics/biokanga/issues\n\n",gszProcName);
		return(1);
        }

    /* special case: '--version' takes precedence error reporting */
if (version->count > 0)
        {
		printf("\n%s Version %s\n",gszProcName,cpszProgVer);
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
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Version: %s",cpszProgVer);

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
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Expected input pseudo-genome suffix array file '-i<name>' is empty");
		return(1);
		}

	strcpy(szMarkerFile,outfile->filename[0]);
	CUtility::TrimQuotedWhitespcExtd(szMarkerFile);
	if(strlen(szMarkerFile) < 1)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Expected marker file to generate '-o<name>' is empty");
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

// show user current resource limits
#ifndef _WIN32
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Resources: %s",CUtility::ReportResourceLimits());
#endif

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing parameters:");
	const char *pszDescr;
	switch(PMode) {
		case ePMExtdKMers:
			pszDescr = "Extended K-mer markers";
			break;
		case ePMNoExtdKMers:
			pszDescr = "Non-extended K-mer markers";
			break;
		}

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Processing mode is : '%s'",pszDescr);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Targeted cultivar name : '%s'",szCultivarName);

	for(int Idx = 0; Idx < NumPartialCultivars; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Targeted cultivar chromosome name (%d) : '%s'", Idx + 1,pszPartialCultivars[Idx]);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Core K-mer length : %d",KMerLen);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Minimum Hamming separation : %d",MinHamming);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Input indexed pseudo-genome file: '%s'",szSfxPseudoGenome);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Write marker K-mers to file: '%s'",szMarkerFile);
	if(szMarkerReadsFile[0] != '\0')
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Write marker containing reads to file: '%s'",szMarkerReadsFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Number of processing threads: %d",NumThreads);

#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	gStopWatch.Start();
	Rslt = LocMarkers(PMode,KMerLen,MinHamming,szCultivarName,NumPartialCultivars,pszPartialCultivars,szSfxPseudoGenome,szMarkerFile,szMarkerReadsFile,NumThreads);
	gStopWatch.Stop();
	Rslt = Rslt >=0 ? 0 : 1;
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Exit code: %d Total processing time: %s",Rslt,gStopWatch.Read());
	return(Rslt);
	}
else
	{
	printf("\n%s SNPs %s, Version %s\n",gszProcName,cpszProcOverview,cpszProgVer);
	arg_print_errors(stdout,end,gszProcName);
	arg_print_syntax(stdout,argtable,"\nUse '-h' to view option and parameter usage\n");
	return(1);
	}
return 0;
}

int
LocMarkers(etPMode PMode,				// processing mode - currently unused, defaults to 0
		  int KMerLen,					// this length K-mers
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

Rslt = Markers.LocKMers(PMode,KMerLen,MinHamming,pszCultivarName,NumPartialCultivars,ppszPartialCultivars,pszSfxPseudoGenome,pszMarkerFile,pszMarkerReadsFile,NumThreads);

Markers.Reset();
return(Rslt);
}