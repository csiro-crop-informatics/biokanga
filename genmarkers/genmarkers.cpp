// genmarkers.cpp : Defines the entry point for the console application.
// 
// 1.0.1     only report markers if no other species has more than a user specified max number of counts at the putative SNP loci 

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

#include "./genmarkers.h"
#include "./Markers.h"

const char *cpszProgVer = "1.0.3";		// increment with each release
const char *cpszProcOverview = "Generate Markers";

CStopWatch gStopWatch;
CDiagnostics gDiagnostics;				// for writing diagnostics messages to log file
char gszProcName[_MAX_FNAME];			// process name


int Process(etPMode PMode,					// processing mode
			int MinCovBases,				// accept SNPs with at least this number covering bases
			double MaxPValue,				// accept SNPs with at most this P-value
			int MinSpeciesTotCntThres,		// individual species must have at least this number of total bases at SNP loci to count as SNP - 0 if no threshold
			int MinSpeciesWithCnts,			// only report markers where at least this number of species has SNP at the SNP loci
			int AltSpeciesMaxCnt,			// only report markers if no other species has more than this number of counts at the putative SNP loci
			char *pszRefGenome,				// reference genome assembly against which other species were aligned
			int NumRelGenomes,				// number of relative genome names
			char *pszRelGenomes[],			// relative genome names
			int NumThreads,					// number of worker threads to use
    		int NumSNPFiles,				// number of input SNP files
			char *pszSNPFiles[],			// names of input files,
			int NumAlignFiles,				// number of input alignment files
			char *pszAlignFiles[],			// names of alignment files
			char *pszMarkerFile);			// output markers to this file				

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
	return _T("genmarkers");
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

int iFileLogLevel;			// level of file diagnostics
int iScreenLogLevel;		// level of file diagnostics
char szLogFile[_MAX_PATH];	// write diagnostics to this file

int Rslt = 0;   			// function result code >= 0 represents success, < 0 on failure
int Idx;

etPMode PMode;				// processing mode

int MinCovBases;				// accept SNPs with at least this number covering bases
double MaxPValue;				// accept SNPs with at most this P-value
int MinSpeciesTotCntThres;		// individual species must have at least this number of total bases at SNP loci to count as SNP - 0 if no threshold
int AltSpeciesMaxCnt;			// only report markers if no other species has more than this number of counts at the putative SNP loci
int MinSpeciesWithCnts;			// only report markers where at least this number of species has SNP at the SNP loci


char szRefGenome[cMaxLenName+1];	// reference genome against which other relative genomes were aligned

int NumRelGenomes;			// number of relative genome names
char *pszRelGenomes[cRRMaxInFileSpecs];  // names of relative genome names

int NumSNPFiles;			// number of input SNP files
char *pszSNPFiles[cRRMaxInFileSpecs];  // input SNP files

int NumAlignFiles;			// number of input alignment files
char *pszAlignFiles[cRRMaxInFileSpecs];  // input alignment files

char szMarkerFile[_MAX_PATH];		// write markers to this file

int NumberOfProcessors;		// number of installed CPUs
int NumThreads;				// number of threads (0 defaults to number of CPUs)

// command line args
struct arg_lit  *help    = arg_lit0("hH","help",                "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"diagnostics log file");

struct arg_int *pmode = arg_int0("m","mode","<int>",		    "Marker processing mode: 0 - default");
struct arg_int *mincovbases=arg_int0("b", "MinCovBases","<int>","Filter out SNPs with less than this number of covering bases (default 5)");
struct arg_dbl *maxpvalue=arg_dbl0("p", "MaxPValue","<dbl>",	"Filter out SNPs with P-Value higher (default 0.05)");

struct arg_int *mintotcntthres = arg_int0("z","mintotcntthres","<int>",	"Species must have at least this number of total bases covering marker loci");
struct arg_int *mincovspecies=arg_int0("Z", "mincovspecies","<int>","Do not report markers unless this minimum number of species have SNP at same loci");
struct arg_int *altspeciesmaxcnt = arg_int0("a","altspeciesmaxcnt","<int>",	"Only report markers if no other species has more than this number of counts at the putative SNP loci (defaults to 1)");

struct arg_str *refgenome=arg_str1("r", "refgenome","<str>",	 "alignments and SNPs of relative genomes were against this genome assembly (default 'RefGenome')");
struct arg_str *relgenomes = arg_strn("R","relgenomes","<relgenomes>",1,cRRMaxInFileSpecs,"alignments and SNPs from these genomes");
struct arg_file *snpfiles = arg_filen("i","insnps","<file>",1,cRRMaxInFileSpecs,"Load SNPs from file(s)");
struct arg_file *alignfiles = arg_filen("I","inaligns","<file>",1,cRRMaxInFileSpecs,"Load alignments from file(s)");

struct arg_file *markerfile = arg_file1("o","out","<file>",		"Output markers to this file");

struct arg_int *threads = arg_int0("T","threads","<int>",		"number of processing threads 0..128 (defaults to 0 which sets threads to number of CPU cores)");
struct arg_end *end = arg_end(100);

void *argtable[] = {help,version,FileLogLevel,LogFile,
	pmode,mincovbases,maxpvalue,mintotcntthres,altspeciesmaxcnt,mincovspecies,refgenome,relgenomes,snpfiles,alignfiles,markerfile,threads,
	end};
char **pAllArgs;
int argerrors;
argerrors = CUtility::arg_parsefromfile(argc,(char **)argv,&pAllArgs);
if(argerrors >= 0)
	argerrors = arg_parse(argerrors,pAllArgs,argtable);

/* special case: '--help' takes precedence over error reporting */
if (help->count > 0)
        {
		printf("\n%s Generate Markers, Version %s\nOptions ---\n", gszProcName,cpszProgVer);
        arg_print_syntax(stdout,argtable,"\n");
        arg_print_glossary(stdout,argtable,"  %-25s %s\n");
		printf("\nNote: Parameters can be entered into a parameter file, one parameter per line.");
		printf("\n      To invoke this parameter file then precede its name with '@'");
		printf("\n      e.g. %s @myparams.txt\n",gszProcName);
		printf("\nPlease report any issues regarding usage of %s at https://github.com/csiro-crop-informatics/biokanga/issues\n\n",gszProcName);
		exit(1);
        }

    /* special case: '--version' takes precedence error reporting */
if (version->count > 0)
        {
		printf("\n%s Version %s\n",gszProcName,cpszProgVer);
		exit(1);
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

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Version: %s",cpszProgVer);

	PMode = (etPMode)(pmode->count ? pmode->ival[0] : ePMdefault);
	if(PMode < ePMdefault || PMode >= ePMplaceholder)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Processing mode '-m%d' specified outside of range %d..%d\n",PMode,ePMdefault,(int)ePMplaceholder-1);
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

	MinCovBases = mincovbases->count ? mincovbases->ival[0] : 5;
	if(MinCovBases < 1 || MinCovBases > 10000)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Minimum covering bases '-b%d' must be in range 1..1000",MinCovBases);
		exit(1);
		}

	AltSpeciesMaxCnt = altspeciesmaxcnt->count ? altspeciesmaxcnt->ival[0] : 1;
	if(AltSpeciesMaxCnt < 1 || AltSpeciesMaxCnt > 10000)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Max alternative species coverage '-a%d' at putative marker loci must be in range 1..1000",AltSpeciesMaxCnt);
		exit(1);
		}


	MaxPValue = maxpvalue->count ? maxpvalue->dval[0] : 0.05;
	if(MaxPValue < 0.0 || MaxPValue > 0.25)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Maximum P-Value '-p%1.4f' must be in range 0.0..0.25",MaxPValue);
		exit(1);
		}

	if(refgenome->count)
		{
		strncpy(szRefGenome,refgenome->sval[0],cMaxLenName);
		szRefGenome[cMaxLenName-1]= '\0';
		}
	else
		strcpy(szRefGenome,"RefGenome");

	if(!relgenomes->count)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: No alignment from genome name(s) specified with with '-R<relgenomes>' option)");
		exit(1);
		}
	for(NumRelGenomes=Idx=0;NumRelGenomes < cRRMaxInFileSpecs && Idx < relgenomes->count; Idx++)
		{
		pszRelGenomes[Idx] = NULL;
		if(pszRelGenomes[NumRelGenomes] == NULL)
			pszRelGenomes[NumRelGenomes] = new char [_MAX_PATH];
		strncpy(pszRelGenomes[NumRelGenomes],relgenomes->sval[Idx],_MAX_PATH);
		pszRelGenomes[NumRelGenomes][_MAX_PATH-1] = '\0';
		CUtility::TrimQuotedWhitespcExtd(pszRelGenomes[NumRelGenomes]);
		if(pszRelGenomes[NumRelGenomes][0] != '\0')
			NumRelGenomes++;
		}


	strncpy(szMarkerFile,markerfile->filename[0],_MAX_PATH);
	szMarkerFile[_MAX_PATH-1] = '\0';

	if(!snpfiles->count)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: No input SNP file(s) specified with with '-i<filespec>' option)");
		exit(1);
		}

	for(NumSNPFiles=Idx=0;NumSNPFiles < cRRMaxInFileSpecs && Idx < snpfiles->count; Idx++)
		{
		pszSNPFiles[Idx] = NULL;
		if(pszSNPFiles[NumSNPFiles] == NULL)
			pszSNPFiles[NumSNPFiles] = new char [_MAX_PATH];
		strncpy(pszSNPFiles[NumSNPFiles],snpfiles->filename[Idx],_MAX_PATH);
		pszSNPFiles[NumSNPFiles][_MAX_PATH-1] = '\0';
		CUtility::TrimQuotedWhitespcExtd(pszSNPFiles[NumSNPFiles]);
		if(pszSNPFiles[NumSNPFiles][0] != '\0')
			NumSNPFiles++;
		}

	if(!NumSNPFiles)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: After removal of whitespace, no input SNP file(s) specified with '-i<filespec>' option");
		exit(1);
		}

	if(!alignfiles->count)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: No input alignment file(s) specified with with '-I<filespec>' option)");
		exit(1);
		}

	for(NumAlignFiles=Idx=0;NumAlignFiles < cRRMaxInFileSpecs && Idx < alignfiles->count; Idx++)
		{
		pszAlignFiles[Idx] = NULL;
		if(pszAlignFiles[NumAlignFiles] == NULL)
			pszAlignFiles[NumAlignFiles] = new char [_MAX_PATH];
		strncpy(pszAlignFiles[NumAlignFiles],alignfiles->filename[Idx],_MAX_PATH);
		pszAlignFiles[NumAlignFiles][_MAX_PATH-1] = '\0';
		CUtility::TrimQuotedWhitespcExtd(pszAlignFiles[NumAlignFiles]);
		if(pszAlignFiles[NumAlignFiles][0] != '\0')
			NumAlignFiles++;
		}

	if(!NumAlignFiles)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: After removal of whitespace, no input alignment file(s) specified with '-I<filespec>' option");
		exit(1);
		}

	// number of alignment files must be same as the number of SNP files and genome names!
	if(NumAlignFiles != NumSNPFiles && NumAlignFiles != NumRelGenomes)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Expected same number of genome names, alignment files and SNP files, %d genome names, %d alignment files, %d SNP files",NumRelGenomes,NumAlignFiles,NumSNPFiles);
		exit(1);
		}
	
	MinSpeciesTotCntThres = mincovspecies->count ? mincovspecies->ival[0] : 1;
	if(MinSpeciesTotCntThres < 1 || MinSpeciesTotCntThres > 10000)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Minimum species total bases '-z%d' must be in range 1..10000",MinSpeciesTotCntThres);
		exit(1);
		}

	MinSpeciesWithCnts = mincovspecies->count ? mincovspecies->ival[0] : 1;
	if(MinSpeciesWithCnts < 1 || MinSpeciesWithCnts > NumAlignFiles)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Minimum species to call marker '-Z%d' must be in range 1..%d",NumAlignFiles);
		exit(1);
		}

// show user current resource limits
#ifndef _WIN32
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Resources: %s",CUtility::ReportResourceLimits());
#endif

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing parameters:");
	const char *pszDescr;

	switch(PMode) {
		case ePMdefault:
			pszDescr = "Default marker processing";
			break;
		}

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Processing mode is : '%s'",pszDescr);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Minimum coverage : %d",MinCovBases);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Maximum P-Value : %1.4f'",MaxPValue);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Reference genome assembly name : '%s'",szRefGenome);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Minimum total bases for species at SNP call loci : %d",MinSpeciesTotCntThres);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Maximum alternative species coverage at SNP call loci : %d",AltSpeciesMaxCnt);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Minimum number of species with SNP call at same loci : %d",MinSpeciesWithCnts);

	for(Idx=0; Idx < NumRelGenomes; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Alignments and SNPs from this genome (%d) : '%s'",Idx+1,pszRelGenomes[Idx]);

	for(Idx=0; Idx < NumSNPFiles; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Input SNP file (%d) : '%s'",Idx+1,pszSNPFiles[Idx]);

	for(Idx=0; Idx < NumAlignFiles; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Input alignment file (%d) : '%s'",Idx+1,pszAlignFiles[Idx]);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Output markers to file : '%s'",szMarkerFile);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"number of threads : %d",NumThreads);

	#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	gStopWatch.Start();
	Rslt = Process(PMode,MinCovBases,MaxPValue,MinSpeciesTotCntThres,MinSpeciesWithCnts,AltSpeciesMaxCnt,szRefGenome,NumRelGenomes,pszRelGenomes,NumThreads,NumSNPFiles,pszSNPFiles,NumAlignFiles,pszAlignFiles,szMarkerFile);
	gStopWatch.Stop();
	Rslt = Rslt >=0 ? 0 : 1;
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Exit code: %d Total processing time: %s",Rslt,gStopWatch.Read());
	exit(Rslt);
	}
else
	{
	printf("\n%s Generate Markers, Version %s\n",gszProcName,cpszProgVer);
	arg_print_errors(stdout,end,gszProcName);
	arg_print_syntax(stdout,argtable,"\nUse '-h' to view option and parameter usage\n");
	exit(1);
	}
return 0;
}


int Process(etPMode PMode,					// processing mode
			int MinCovBases,				// accept SNPs with at least this number covering bases
			double MaxPValue,				// accept SNPs with at most this P-value
			int MinSpeciesTotCntThres,		// individual species must have at least this number of total bases at SNP loci to count as SNP - 0 if no threshold
			int MinSpeciesWithCnts,			// only report markers where at least this number of species has SNP at the SNP loci
			int AltSpeciesMaxCnt,			// only report markers if no other species has more than this number of counts at the putative SNP loci
			char *pszRefGenome,				// reference genome assembly against which other species were aligned
			int NumRelGenomes,				// number of relative genome names
			char *pszRelGenomes[],			// relative genome names
			int NumThreads,					// number of worker threads to use
    		int NumSNPFiles,				// number of input SNP files
			char *pszSNPFiles[],			// names of input files,
			int NumAlignFiles,				// number of input alignment files
			char *pszAlignFiles[],			// names of alignment files
			char *pszMarkerFile)			// output markers to this file	
{
char *pszSNPFile;
char *pszAlignFile;
char szProbeSpecies[cMaxLenName];
int FileIdx;
int Rslt;
CMarkers *pMarkers = NULL;

if((pMarkers = new CMarkers)==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CMarkers");
	return(eBSFerrObj);
	}

for(FileIdx = 0; FileIdx < NumSNPFiles; FileIdx++)
	{
	pszSNPFile = pszSNPFiles[FileIdx];
	sprintf(szProbeSpecies,"ProbeSpecies%d",FileIdx+1);

	Rslt = pMarkers->LoadSNPFile(MinCovBases,	// accept SNPs with at least this number covering bases
					MaxPValue,					// accept SNPs with at most this P-value
					pszRefGenome,				// this is the reference species 
					szProbeSpecies,				// this species reads were aligned to the reference species from which SNPs were called 
					pszSNPFile);				// SNP file to parse and load
	if(Rslt < 0)
		break;
	}

// if no SNPs are being called was it because there was no coverage?
// sort the SNPs by refseq.loci.probespecies ascending
// iterate looking for missing probespecies at each refseq.loci
// load the probespecies alignments and check for coverage at the missing refseq.loci

for(FileIdx = 0; FileIdx < NumAlignFiles; FileIdx++)
	{
	pszAlignFile = pszAlignFiles[FileIdx];
	sprintf(szProbeSpecies,"ProbeSpecies%d",FileIdx+1);

	Rslt = pMarkers->AddImputedAlignments(MinCovBases,	// must be at least this number of reads covering the SNP loci
					pszRefGenome,				// this is the reference species 
					szProbeSpecies,				// this species reads were aligned to the reference species from which SNPs were called 
					pszAlignFile);				// alignment file to parse and load
	if(Rslt < 0)
		break;
	}

pMarkers->SortTargSeqLociSpecies();



pMarkers->IdentSpeciesSpec(AltSpeciesMaxCnt,	// max count allowed for base being processed in any other species
						MinCovBases,			// min count required for base being processed in species
						cMinBaseThres);			// to be processed a base must be at least this proportion of total

pMarkers->Report(pszRefGenome,NumRelGenomes,pszRelGenomes,pszMarkerFile,MinSpeciesWithCnts,MinSpeciesTotCntThres);

delete pMarkers;
return(0);
}

