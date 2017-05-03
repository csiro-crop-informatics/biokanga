/*
 * CSIRO Open Source Software License Agreement (GPLv3)
 * Copyright (c) 2017, Commonwealth Scientific and Industrial Research Organisation (CSIRO) ABN 41 687 119 230.
 * See LICENSE for the complete license information (https://github.com/csiro-crop-informatics/biokanga/LICENSE)
 * Contact: Alex Whan <alex.whan@csiro.au>
 */

// gensnpmarkers.cpp : Defines the entry point for the console application.
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
#include "./Markers.h"

// report processing modes
typedef enum TAG_eRPMode {
	eRPMdefault,				// default processing mode
	eRPMInterCultOnly,			// must be inter-cultivar SNP markers
	eRPMplaceholder				// used to set the enumeration range
	} etRPMode;


int GSMProcess(etRPMode PMode,				// report processing mode
			int MinCovBases,				// accept SNPs with at least this number covering bases
			double MaxPValue,				// accept SNPs with at most this P-value
		    double SNPMmajorPC,				// only accept for processing if more/equal than this percentage number of reads are major SNP at putative SNP loci (defaults to 50.0) 
			int MinSpeciesTotCntThres,		// individual species must have at least this number of total bases at SNP loci to count as SNP - 0 if no threshold
			int MinSpeciesWithCnts,			// only report markers where at least this number of species has SNP at the SNP loci
			int AltSpeciesMaxCnt,			// only report markers if no other species has more than this number of counts at the putative SNP loci, 0 if no limit
			char *pszRefGenome,				// reference genome assembly against which other species were aligned
			int NumRelGenomes,				// number of relative genome names
			char *pszRelGenomes[],			// relative genome names
    		int NumSNPFiles,				// number of input SNP files
			char *pszSNPFiles[],			// names of input files,
			int NumAlignFiles,				// number of input alignment files
			char *pszAlignFiles[],			// names of alignment files
			char *pszMarkerFile);			// output markers to this file				

#ifdef _WIN32
int gensnpmarkers(int argc, char* argv[])
{
// determine my process name
_splitpath(argv[0],NULL,NULL,gszProcName,NULL);
#else
int gensnpmarkers(int argc, char** argv)
{
// determine my process name
CUtility::splitpath((char *)argv[0],NULL,gszProcName);
#endif

int iFileLogLevel;			// level of file diagnostics
int iScreenLogLevel;		// level of file diagnostics
char szLogFile[_MAX_PATH];	// write diagnostics to this file

int Rslt = 0;   			// function result code >= 0 represents success, < 0 on failure
int Idx;

int PMode;				// processing mode

int MinCovBases;				// accept SNPs with at least this number covering bases
double MaxPValue;				// accept SNPs with at most this P-value
double SNPMmajorPC;				// only accept for processing if more/equal than this percentage number of reads are major SNP at putative SNP loci (defaults to 50.0) 
int MinSpeciesTotCntThres;		// individual species must have at least this number of total bases at SNP loci to count as SNP - 0 if no threshold
int AltSpeciesMaxCnt;			// only report markers if no other species has more than this number of counts at the putative SNP loci - 0 if no limit
int MinSpeciesWithCnts;			// only report markers where at least this number of species has SNP at the SNP loci


char szRefGenome[cMaxLenName+1];	// reference genome against which other relative genomes were aligned

int NumRelGenomes;			// number of relative genome names
char *pszRelGenomes[cMaxMarkerSpecies+1];  // names of relative genome names

int NumSNPFiles;			// number of input SNP files
char *pszSNPFiles[cMaxMarkerSpecies+1];  // input SNP files

int NumAlignFiles;			// number of input alignment files
char *pszAlignFiles[cMaxMarkerSpecies+1];  // input alignment files

char szMarkerFile[_MAX_PATH];		// write markers to this file


char szSQLiteDatabase[_MAX_PATH];	// results summaries to this SQLite file
char szExperimentName[cMaxDatasetSpeciesChrom+1];			// experiment name
char szExperimentDescr[1000];		// describes experiment


int NumberOfProcessors;		// number of installed CPUs

// command line args
struct arg_lit  *help    = arg_lit0("hH","help",                "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"diagnostics log file");

struct arg_int *pmode = arg_int0("m","mode","<int>",		    "Marker reporting mode: 0 - SNP markers if either inter-species/cultivar or relative to reference, 1 - report SNP markers only if inter-species/cultivar differences (default 0)");
struct arg_int *mincovbases=arg_int0("b", "mincovbases","<int>","Filter out SNPs with less than this number of covering bases (default 5)");
struct arg_dbl *maxpvalue=arg_dbl0("p", "maxpvalue","<dbl>",	"Filter out SNPs with P-Value higher (default 0.05)");
struct arg_dbl *snpmajorpc = arg_dbl0("P", "snpmajorpc", "<dbl>", "Min percentage Major SNP at putative loci (defaults to 50.0, range 15.0 to 90.0)");

struct arg_int *mintotcntthres = arg_int0("z","mintotcntthres","<int>",	"Species must have at least this number of total bases covering marker loci (defaults to 0 for no limit)");
struct arg_int *mincovspecies=arg_int0("Z", "mincovspecies","<int>","Do not report marker SNPs unless this minimum number of species have SNP at same loci");
struct arg_int *altspeciesmaxcnt = arg_int0("a","altspeciesmaxcnt","<int>",	"Only report markers if no other species has more than this number of counts at the putative SNP loci (defaults to 0 for no limit)");

struct arg_str *refgenome=arg_str1("r", "refgenome","<str>",	 "alignments and SNPs of relative genomes were against this reference genome assembly (default 'RefGenome')");
struct arg_str *relgenomes = arg_strn("R","relgenomes","<relgenomes>",1,cMaxMarkerSpecies,"alignments and SNPs from these species or cultivars");
struct arg_file *snpfiles = arg_filen("i","insnps","<file>",1,cMaxMarkerSpecies,"Load SNPs from file(s)");
struct arg_file *alignfiles = arg_filen("I","inaligns","<file>",1,cMaxMarkerSpecies,"Load alignments from file(s)");

struct arg_file *markerfile = arg_file1("o","out","<file>",		"Output marker SNP loci to this file");
struct arg_file *summrslts = arg_file0("q","sumrslts","<file>",		"Output results summary to this SQLite3 database file");
struct arg_str *experimentname = arg_str0("w","experimentname","<str>",		"experiment name SQLite3 database file");
struct arg_str *experimentdescr = arg_str0("W","experimentdescr","<str>",	"experiment description SQLite3 database file");
struct arg_end *end = arg_end(200);

void *argtable[] = {help,version,FileLogLevel,LogFile,
	summrslts,experimentname,experimentdescr,
	pmode,mincovbases,maxpvalue,snpmajorpc,mintotcntthres,altspeciesmaxcnt,mincovspecies,refgenome,relgenomes,snpfiles,alignfiles,markerfile,
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
	PMode = (etRPMode)(pmode->count ? pmode->ival[0] : eRPMdefault);
	if(PMode < eRPMdefault || PMode >= eRPMplaceholder)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Report processing mode '-m%d' specified outside of range %d..%d\n",PMode,eRPMdefault,(int)eRPMplaceholder-1);
		exit(1);
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
	MinCovBases = mincovbases->count ? mincovbases->ival[0] : 5;
	if(MinCovBases < 1 || MinCovBases > 10000)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Minimum covering bases '-b%d' must be in range 1..10000",MinCovBases);
		exit(1);
		}

	AltSpeciesMaxCnt = altspeciesmaxcnt->count ? altspeciesmaxcnt->ival[0] : 0;
	if(AltSpeciesMaxCnt < 0 || AltSpeciesMaxCnt > 10000)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Max alternative species coverage '-a%d' at putative marker loci must be in range 0..10000",AltSpeciesMaxCnt);
		exit(1);
		}


	MaxPValue = maxpvalue->count ? maxpvalue->dval[0] : 0.05;
	if(MaxPValue < 0.0 || MaxPValue > 0.25)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Maximum P-Value '-p%1.4f' must be in range 0.0..0.25",MaxPValue);
		exit(1);
		}

	SNPMmajorPC = snpmajorpc->count ? snpmajorpc->dval[0] : 50.0;
	if (SNPMmajorPC < 15.0 || SNPMmajorPC > 90.0)
	{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: Major SNP minimum non-ref '-P%f' for controlling SNP FDR must be in range 15.0 to 90.0\n", SNPMmajorPC);
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
	for(NumRelGenomes=Idx=0;NumRelGenomes < cMaxMarkerSpecies && Idx < relgenomes->count; Idx++)
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

	for(NumSNPFiles=Idx=0;NumSNPFiles < cMaxMarkerSpecies && Idx < snpfiles->count; Idx++)
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

	for(NumAlignFiles=Idx=0;NumAlignFiles < cMaxMarkerSpecies && Idx < alignfiles->count; Idx++)
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
	
	MinSpeciesTotCntThres = mintotcntthres->count ? mintotcntthres->ival[0] : 0;
	if(MinSpeciesTotCntThres < 0 || MinSpeciesTotCntThres > 10000)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Minimum species total bases '-z%d' must be in range 0..10000",MinSpeciesTotCntThres);
		exit(1);
		}

	MinSpeciesWithCnts = mincovspecies->count ? mincovspecies->ival[0] : 1;
	if(MinSpeciesWithCnts < 1 || MinSpeciesWithCnts > NumAlignFiles)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Minimum species to call marker '-Z%d' must be in range 1..%d",NumAlignFiles);
		exit(1);
		}

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing parameters:");
	const char *pszDescr;

	switch(PMode) {
		case eRPMdefault:
			pszDescr = "Reference relative or inter-cultivar SNP markers";
			break;
		case eRPMInterCultOnly:
			pszDescr = "Inter-cultivar SNP markers only";
			break;
		}

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Processing mode is : '%s'",pszDescr);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Minimum coverage : %d",MinCovBases);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Maximum P-Value : %1.4f'",MaxPValue);
	gDiagnostics.DiagOutMsgOnly(eDLInfo, "Minimum Major SNP percentage : %1.4f'", SNPMmajorPC);


	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Reference genome assembly name : '%s'",szRefGenome);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Minimum total bases for species at SNP call loci : %d",MinSpeciesTotCntThres);
	if(AltSpeciesMaxCnt)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Maximum alternative species coverage at SNP call loci : %d",AltSpeciesMaxCnt);
	else
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Maximum alternative species coverage at SNP call loci : unlimited");

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Minimum number of species with SNP call at same loci : %d",MinSpeciesWithCnts);

	for(Idx=0; Idx < NumRelGenomes; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Alignments and SNPs from this species or cultivar (%d) : '%s'",Idx+1,pszRelGenomes[Idx]);

	for(Idx=0; Idx < NumSNPFiles; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Input SNP file (%d) : '%s'",Idx+1,pszSNPFiles[Idx]);

	for(Idx=0; Idx < NumAlignFiles; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Input alignment file (%d) : '%s'",Idx+1,pszAlignFiles[Idx]);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Output markers to file : '%s'",szMarkerFile);

	if(szExperimentName[0] != '\0')
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"This processing reference: %s",szExperimentName);

	if(gExperimentID > 0)
		{
		int ParamID;
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szLogFile),"log",szLogFile);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,sizeof(PMode),"mode",&PMode);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,sizeof(MinCovBases),"mincovbases",&MinCovBases);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTDouble,sizeof(MaxPValue),"maxpvalue",&MaxPValue);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID, ePTDouble, sizeof(SNPMmajorPC), "snpmajorpc", &SNPMmajorPC);

		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,sizeof(AltSpeciesMaxCnt),"altspeciesmaxcnt",&AltSpeciesMaxCnt);

		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,sizeof(MinSpeciesWithCnts),"mincovspecies",&MinSpeciesWithCnts);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,sizeof(MinSpeciesTotCntThres),"mintotcntthres",&MinSpeciesTotCntThres);

		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szRefGenome),"cultivar",szRefGenome);

		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,sizeof(NumRelGenomes),"NumRelGenomes",&NumRelGenomes);
		for(Idx=0; Idx < NumRelGenomes; Idx++)
			ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(pszRelGenomes[Idx]),"relgenomes",pszRelGenomes[Idx]);

		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,sizeof(NumSNPFiles),"NumSNPFiles",&NumSNPFiles);
		for(Idx=0; Idx < NumSNPFiles; Idx++)
			ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(pszSNPFiles[Idx]),"insnps",pszSNPFiles[Idx]);

		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,sizeof(NumAlignFiles),"NumAlignFiles",&NumAlignFiles);
		for(Idx=0; Idx < NumAlignFiles; Idx++)
			ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(pszAlignFiles[Idx]),"inaligns",pszAlignFiles[Idx]);

		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szMarkerFile),"out",szMarkerFile);

		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,sizeof(NumberOfProcessors),"cpus",&NumberOfProcessors);

		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szSQLiteDatabase),"sumrslts",szSQLiteDatabase);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szExperimentName),"experimentname",szExperimentName);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szExperimentDescr),"experimentdescr",szExperimentDescr);
		}


#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	gStopWatch.Start();
	Rslt = GSMProcess((etRPMode)PMode,MinCovBases,MaxPValue, SNPMmajorPC,MinSpeciesTotCntThres,MinSpeciesWithCnts,AltSpeciesMaxCnt,szRefGenome,NumRelGenomes,pszRelGenomes,NumSNPFiles,pszSNPFiles,NumAlignFiles,pszAlignFiles,szMarkerFile);
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

typedef struct TAG_sSAMFileEls {
	UINT32 SeqLen;				// estimated mean sequence length
	UINT32 NumAlignments;		// estimated number of alignments		
	} tsSAMFileEls;

int GSMProcess(etRPMode PMode,				// processing mode
			int MinCovBases,				// accept SNPs with at least this number covering bases
			double MaxPValue,				// accept SNPs with at most this P-value
		    double SNPMmajorPC,				// only accept for processing if more/equal than this percentage number of reads are major SNP at putative SNP loci (defaults to 50.0) 
			int MinSpeciesTotCntThres,		// individual species must have at least this number of total bases at SNP loci to count as SNP - 0 if no threshold
			int MinSpeciesWithCnts,			// only report markers where at least this number of species has SNP at the SNP loci
			int AltSpeciesMaxCnt,			// only report markers if no other species has more than this number of counts at the putative SNP loci, 0 if no limit
			char *pszRefGenome,				// reference genome assembly against which other species were aligned
			int NumRelGenomes,				// number of relative genome names
			char *pszRelGenomes[],			// relative genome names
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
INT64 Rslt64;
tsSAMFileEls EstSAMFileEls[cMaxMarkerSpecies];

INT64 PrevAlignLoci;
INT64 CurAlignLoci;
INT64 InitalAlignLoci;
INT64 TotSNPRows;
INT64 TotAlignments;
INT64 SumMeanSeqLens;
INT32 MeanSeqLen;
CMarkers *pMarkers = NULL;
CCSVFile *pCSV = NULL;
CSAMfile *pSAM = NULL;
size_t TotMemToAlloc;

if((pMarkers = new CMarkers)==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CMarkers");
	return(eBSFerrObj);
	}

if((pCSV = new CCSVFile)==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CCSVFile");
	return(eBSFerrObj);
	}

// try to guestimate minimum memory requirements
TotSNPRows = 0;
for(FileIdx = 0; FileIdx < NumSNPFiles; FileIdx++)
	{
	UINT32 NumRows;
	INT64 FileSize;
	int MaxChrsPerRow;
	int MaxFields;
	int MeanFields;
	int MeanChrsRow;
	
	pszSNPFile = pszSNPFiles[FileIdx];
	if(((NumRows = pCSV->CSVEstSizes(pszSNPFile,&FileSize,&MaxFields,&MeanFields,&MaxChrsPerRow,&MeanChrsRow))==0) || MeanFields < 23)
		{
		delete pCSV;
		delete pMarkers;
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to estimate memory required for SNPs in file: '%s'",pszSNPFile);
		if(MeanFields < 23)
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expected at least 23 fields per row if a SNP file");
		return(eBSFerrParse);
		}
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Estimating %d SNPs in file: '%s'",NumRows,pszSNPFile);		
	TotSNPRows += NumRows; 
	}
delete pCSV;

if((pSAM = new CSAMfile)==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CSAMfile");
	return(eBSFerrObj);
	}
TotAlignments = 0;
SumMeanSeqLens = 0;
int MaxEstSAMFileIdx = 0;
size_t EstSAMFileMem;
size_t MaxEstSAMFileMem;

MaxEstSAMFileMem = 0;
for(FileIdx = 0; FileIdx < NumAlignFiles; FileIdx++)
	{
	UINT32 NumAlignments;
	EstSAMFileEls[FileIdx].NumAlignments = 0;
	EstSAMFileEls[FileIdx].SeqLen = 0;
	pszAlignFile = pszAlignFiles[FileIdx];
	
	if((NumAlignments = pSAM->EstSizes(pszAlignFile,NULL,NULL,NULL,NULL,&MeanSeqLen,NULL))==0)
		{
		delete pSAM;
		delete pMarkers;
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to estimate memory required for alignments in file: '%s'",pszAlignFile);
		return(eBSFerrParse);
		}
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Estimating %u alignments with mean sequence length %u in file: '%s'",NumAlignments,MeanSeqLen,pszAlignFile);

	NumAlignments = (UINT32)(((UINT64)NumAlignments * 105)/100); // allow 5% extra when allocating memory in case underestimating
	MeanSeqLen = (UINT32)(((UINT64)MeanSeqLen * 105)/100);


	EstSAMFileEls[FileIdx].NumAlignments = NumAlignments;
	EstSAMFileEls[FileIdx].SeqLen = MeanSeqLen;
	EstSAMFileMem = NumAlignments * (sizeof(tsHyperElement) + MeanSeqLen);
	if(EstSAMFileMem > MaxEstSAMFileMem)
		{
		MaxEstSAMFileIdx = FileIdx;
		MaxEstSAMFileMem = EstSAMFileMem;
		}
	TotAlignments += NumAlignments;
	}
delete pSAM;
// estimate a minimum total memory required
TotSNPRows *= 11;					// assume a 10x overhead for additional SNPs from imputed alignments
TotMemToAlloc = MaxEstSAMFileMem;
TotMemToAlloc += TotSNPRows * sizeof(tsAlignLoci);
// allow 15% overhead for indexes, reallocs, other allocations etc
TotMemToAlloc = (TotMemToAlloc * 115)/100;
gDiagnostics.DiagOut(eDLFatal,gszProcName,"Estimating minimum total memory requirements to be: %dGB",(int)((TotMemToAlloc + 0x040000000 - 1)/0x040000000));

// try to prealloc for SNPs 
gDiagnostics.DiagOut(eDLFatal,gszProcName,"Pre-allocating memory for %lld SNP loci allowing 10x additional for impuned SNPS",TotSNPRows);

if((Rslt = pMarkers->PreAllocSNPs(TotSNPRows)) != eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to pre-allocate %dGB memory for SNP loci",TotSNPRows * sizeof(tsAlignLoci));
	delete pMarkers;
	return(eBSFerrMem);
	}

// then try an allocation of memory for the maximal sized SAM file sequences to check if allocations whilst imputing SNPs are likely to be successful
CHyperEls *pHyperEls;
if((pHyperEls = new CHyperEls) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CHyperEls");
	delete pMarkers;
	return(eBSFerrObj);
	}
gDiagnostics.DiagOut(eDLFatal,gszProcName,"Pre-allocating memory for maximal %u SAM alignments with mean sequence length of %u in any SAM file allowing 5%% additional for underestimation",EstSAMFileEls[MaxEstSAMFileIdx].NumAlignments, EstSAMFileEls[MaxEstSAMFileIdx].SeqLen);

if((Rslt = pHyperEls->PreAllocMem(EstSAMFileEls[MaxEstSAMFileIdx].NumAlignments,EstSAMFileEls[MaxEstSAMFileIdx].SeqLen)) != eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to pre-allocate memory for SAM alignments");
	delete pMarkers;
	return(eBSFerrMem);
	}
delete pHyperEls;

// load all aligner identified  SNPS
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
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadSNPFile(%s) returned error",pszSNPFile);
		delete pMarkers;
		return(Rslt);
		}
	}

if(Rslt == 0)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Nothing to do - no SNPs to process for markers!");
	delete pMarkers;
	return(Rslt);
	}

CurAlignLoci = pMarkers->NumAlignLoci();			// report on total number of accepted SNP alignments
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Accepted a total of %lld SNP alignments to further process for markers",CurAlignLoci);

// if no SNPs are being called at a loci for one cultivar then was it because there was no coverage?
// sort the SNPs by refseq.loci.probespecies ascending
// iterate looking for missing probespecies at each refseq.loci
// load the probespecies alignments and check for coverage at the missing refseq.loci
gDiagnostics.DiagOut(eDLFatal,gszProcName,"Now checking for imputed alignments where no SNP called in one or more cultivars...");

InitalAlignLoci = PrevAlignLoci = CurAlignLoci;
for(FileIdx = 0; FileIdx < NumAlignFiles; FileIdx++)
	{
	pszAlignFile = pszAlignFiles[FileIdx];
	sprintf(szProbeSpecies,"ProbeSpecies%d",FileIdx+1);
	Rslt64 = pMarkers->AddImputedAlignments(MinCovBases,	// must be at least this number of reads covering the SNP loci
					pszRefGenome,				// this is the reference species 
					szProbeSpecies,				// this species reads were aligned to the reference species from which SNPs were called 
					pszAlignFile,0,true,EstSAMFileEls[FileIdx].NumAlignments, EstSAMFileEls[FileIdx].SeqLen);		// alignment file to parse and load
	if(Rslt64 < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddImputedAlignments('%s') returned error %d",pszAlignFile,(int)Rslt64);
		delete pMarkers;
		return((int)Rslt64);
		}
	CurAlignLoci =  pMarkers->NumAlignLoci();
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"AddImputedAlignments('%s') imputed %lld alignments, subtotal %lld",pszAlignFile,CurAlignLoci - PrevAlignLoci, CurAlignLoci - InitalAlignLoci);
	PrevAlignLoci = CurAlignLoci;
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Total alignments %lld or which %lld are imputed alignments",CurAlignLoci, CurAlignLoci - InitalAlignLoci);

Rslt64 = pMarkers->SortTargSeqLociSpecies();
pMarkers->IdentSpeciesSpec(AltSpeciesMaxCnt,	// max count allowed for base being processed in any other species, 0 if no limit
						MinCovBases,			// min count required for base being processed in species
						   SNPMmajorPC);		// to be processed major putative SNP base must be at least this proportion of total

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Reporting markers to '%s'...",pszMarkerFile);
Rslt64 = pMarkers->Report(pszRefGenome,NumRelGenomes,pszRelGenomes,pszMarkerFile,MinSpeciesWithCnts,MinSpeciesTotCntThres,PMode == eRPMInterCultOnly ? true : false);
if(Rslt64 < 0)
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Reporting of markers to '%s' error %d",pszMarkerFile,(int)Rslt64);
else
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Reporting of %lld markers to '%s' completed",Rslt64,pszMarkerFile);
delete pMarkers;
return(Rslt64 < 0 ? (int)Rslt64 : eBSFSuccess);
}

