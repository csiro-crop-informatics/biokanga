/*
 * CSIRO Open Source Software License Agreement (variation of the BSD / MIT License)
 * Copyright (c) 2017, Commonwealth Scientific and Industrial Research Organisation (CSIRO) ABN 41 687 119 230.
 * See LICENCE for the complete licence information (https://github.com/csiro-crop-informatics/biokanga/LICENSE)
 * Contact: Alex Whan <alex.whan@csiro.au>
 */

// maploci2features.cpp : Defines the entry point for the console application.

#include "stdafx.h"
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "../libbiokanga/commhdrs.h"

#include "biokanga.h"
#include "MapLoci2Feat.h"

int MLFProcess(etMLFPMode PMode,		// processing mode
			bool bDedupe,				// true if input elements are to be deduped
			bool bFeatinsts,			// true if input elements are to be associated to individual features, false if to all features at that locus
			bool bOneCntRead,			// true if one count per read rule to be applied (functional regions are prioritised with CDS as the highest) 
			etISOFProc IsoformRprt,		// feature isoform reporting mode
			etStrandProc StrandProc,	// how to process read + element strand
			int Ftype,					// input element file format: 0 - auto, 1 - CSV, 2 - BED, 3 - SAM
			teCSVFormat CSVFormat,		// expected input CSV loci file format
			char *pszInLociFile,		// input CSV loci file
			char *pszInBEDFile,			// input BED file
			char *pszRsltsFile,			// optional output loci mapping file
			char *pszFeatRsltsFile,		// optional output feature mapping results file
			char *pszSummRsltsFile,		// optional output chrom summary results file
			int RegRegionLen,			// regulatory region length
			int MinLength,				// minimum element length
			int MaxLength,				// maximum element length
			int JoinOverlap);			// deduping join overlap

char *CSVFormat2Text(teCSVFormat Format);


#ifdef _WIN32
int maploci2features(int argc, char* argv[])
{
// determine my process name
_splitpath(argv[0],NULL,NULL,gszProcName,NULL);
#else
int
maploci2features(int argc, char** argv)
{
// determine my process name
CUtility::splitpath((char *)argv[0],NULL,gszProcName);
#endif
int iScreenLogLevel;		// level of screen diagnostics
int iFileLogLevel;			// level of file diagnostics
char szLogFile[_MAX_PATH];	// write diagnostics to this file

int Rslt;

int PMode;				// processing mode
int iCSVFormat;
int FType;					// expected input element file type - auto, CSV, BED or SAM
int iRegLen;				// up/down stream regulatory region length
int iMinLength;
int iMaxLength;
int iJoinOverlap;
bool bDedupe;				// if true then dedupe input elements
bool bFeatinsts;			// if true then associate to individual features
bool bOneCntRead;			// true if one count per read rule to be applied (functional regions are prioritised with CDS as the highest) 
etISOFProc IsoformRprt;			// feature isoform report processing
etStrandProc StrandProc;	// how to process read + element strand

char szInLociFile[_MAX_PATH];	// input element loci from this file
char szInBEDFile[_MAX_PATH];	// input bed file containing gene features
char szRsltsFile[_MAX_PATH];	// output loci + features to this file
char szFeatRsltsFile[_MAX_PATH]; // optional features output file
char szSummRsltsFile[_MAX_PATH]; // optional summary output file


char szSQLiteDatabase[_MAX_PATH];	// results summaries to this SQLite file
char szExperimentName[cMaxDatasetSpeciesChrom+1];			// experiment name
char szExperimentDescr[1000];		// describes experiment

// command line args
struct arg_lit  *help    = arg_lit0("hH","help",                "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"diagnostics log file");

struct arg_int *pmode = arg_int0("m","mode","<int>",		    "processing mode: 0 - full length attributes, 1 - 5' attributes associations, 2 - 5' dyad attribute associations (default = 0)");

struct arg_int *ftype = arg_int0("t","filetype","<int>",		"input element file format: 0 - auto, 1 - CSV, 2 - BED, 3 - SAM (default = 0)");
struct arg_lit  *dedupe = arg_lit0("d","dedupe",                "dedupe the input elements (default is to process all input elements)");
struct arg_int  *isoformrprt = arg_int0("x","isoform","<int>",  "feature isoform processing: 0 - report highest RPKM, 1 - report highest total reads, 2 - report all)");
struct arg_int  *strandproc = arg_int0("s","strandproc","<int>","strand processing: 0 - independent, 1 - sense, 2 - antisense (default is independent)");
struct arg_int  *CSVFormat = arg_int0("c","informat","<int>",	"if input CSV file processing 0:Loci 1:locsfx probe 2:locsfx target");
struct arg_file *InLociFile = arg_file1("i","inloci","<file>",	"input element loci (CSV, BED, SAM) file");
struct arg_file *InBEDFile = arg_file1("I","inbed","<file>",	"input gene or feature biobed BED file");
struct arg_file *RsltsFile = arg_file0("o","output","<file>",	 "optional mapped element output file");
struct arg_file *featrsltsfile = arg_file0("O","feats","<file>", "optional mapped features output file");
struct arg_file *summrsltsfile = arg_file0("C","chromsumm","<file>","optional chrom mapping summary output file");
struct arg_int  *RegLen = arg_int0("r","updnstream","<int>",	"length of 5'up or 3'down  stream regulatory region length (default = 2000) 0..1000000");
struct arg_int  *MinLength = arg_int0("l","minlength","<int>",   "minimum element length (default 0)");
struct arg_int  *MaxLength = arg_int0("L","maxlength","<int>",   "maximum element length (default 1000000)");
struct arg_int  *JoinOverlap = arg_int0("j","joinoverlap","<int>","joins cores which only differ by this many bases in start loci (default 0)");
struct arg_lit  *featinsts    = arg_lit0("z","featinsts",        "associate to individual features");
struct arg_lit  *onecntread   = arg_lit0("y","onecntread",       "prioritise functional regions so that there is only one count per read");
struct arg_file *summrslts = arg_file0("q","sumrslts","<file>",		"Output results summary to this SQLite3 database file");
struct arg_str *experimentname = arg_str0("w","experimentname","<str>",		"experiment name SQLite3 database file");
struct arg_str *experimentdescr = arg_str0("W","experimentdescr","<str>",	"experiment description SQLite3 database file");

struct arg_end *end = arg_end(200);

void *argtable[] = {help,version,FileLogLevel,LogFile,
					summrslts,experimentname,experimentdescr,
					pmode,dedupe,ftype,featinsts,isoformrprt,onecntread,strandproc,CSVFormat,InLociFile,InBEDFile,RsltsFile,featrsltsfile,summrsltsfile,RegLen,
					MinLength,MaxLength,JoinOverlap,
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
		printf("\nError: FileLogLevel '-f%d' specified but no logfile '-F<logfile>'",FileLogLevel->ival[0]);
		exit(1);
		}

	iScreenLogLevel = iFileLogLevel = FileLogLevel->count ? FileLogLevel->ival[0] : eDLInfo;
	if(iFileLogLevel < eDLNone || iFileLogLevel > eDLDebug)
		{
		printf("\nError: FileLogLevel '-l%d' specified outside of range %d..%d",iFileLogLevel,eDLNone,eDLDebug);
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

	if(!gDiagnostics.Open(szLogFile,(etDiagLevel)iScreenLogLevel,(etDiagLevel)iFileLogLevel,true))
		{
		printf("\nError: Unable to start diagnostics subsystem.");
		if(szLogFile[0] != '\0')
			printf(" Most likely cause is that logfile '%s' can't be opened/created",szLogFile);
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
	PMode = (etMLFPMode)(pmode->count ? pmode->ival[0] : ePMdefault);
	if(PMode < ePMdefault || PMode >= ePMplaceholder)
		{
		printf("\nError: Processing mode '-m%d' specified outside of range %d..%d",PMode,0,(int)ePMplaceholder-1);
		exit(1);
		}
	bDedupe = dedupe->count ? true : false;
	bFeatinsts = featinsts->count ? true : false;
	bOneCntRead = onecntread->count ? true : false;

	FType = ftype->count ? ftype->ival[0] : 0;
	if(FType < 0 || FType >= 4)
		{
		printf("\nError: Expected input element file format '-t%d' specified outside of range %d..%d",FType,0,3);
		exit(1);
		}

	IsoformRprt = (etISOFProc)(isoformrprt->count ? isoformrprt->ival[0] : eISOFPRPKM);
	if(IsoformRprt < eISOFPRPKM || IsoformRprt >= eISOFPlaceholder)
		{
		printf("\nError: Isoform report mode '-x%d' must be in range %d..%d",IsoformRprt,eISOFPRPKM,eISOFPlaceholder-1);
		exit(1);
		}

	StrandProc = (etStrandProc)(strandproc->count ? strandproc->ival[0] : eStrandDflt);
	if(StrandProc < eStrandDflt || StrandProc >= eStrandPlaceholder)
		{
		printf("\nError: Strand processing mode '-s%d' must be in range %d..%d",StrandProc,eStrandDflt,eStrandAnti);
		exit(1);
		}
	
	if(FType <= 1)
		{
		iCSVFormat = CSVFormat->count ? CSVFormat->ival[0] : 0;
		if(iCSVFormat < eCSVFdefault || iCSVFormat > eCSVFtarget)
			{
			printf("\nError: expected input CSV format specified '-c%d' must be in range 0..2",iCSVFormat);
			exit(1);
			}
		}
	else
		iCSVFormat = 0;

	iMinLength = MinLength->count ? MinLength->ival[0] : 0;
	if(iMinLength < 0 || iMinLength > cMaxLengthRange)
		{
		printf("Error: Minimum element length '-m%d' is not in range 0..%d",iMinLength,cMaxLengthRange);
		exit(1);
		}

	iMaxLength = MaxLength->count ? MaxLength->ival[0] : cMaxLengthRange;
	if(iMaxLength < iMinLength || iMaxLength > cMaxLengthRange)
		{
		printf("Error: Maximum element length '-M%d' is not in range %d..%d",iMaxLength,iMinLength,cMaxLengthRange);
		exit(1);
		}

	iJoinOverlap = JoinOverlap->count ? JoinOverlap->ival[0] : cDfltJoinOverlap;
	if(iJoinOverlap < 0 || iJoinOverlap > cMaxJoinOverlap)
		{
		printf("Error: join overlap length '-j%d' is not in range 0..%d",iJoinOverlap,cMaxJoinOverlap);
		exit(1);
		}

	if(InBEDFile->count)
		{
		strncpy(szInBEDFile,InBEDFile->filename[0],_MAX_PATH);
		szInBEDFile[_MAX_PATH-1] = '\0';
		iRegLen = RegLen->count ? RegLen->ival[0] : cDfltRegLen;
		if(iRegLen < cMinRegLen)
			{
			printf("\nRegulatory region length '-L%d' less than minimum %d, assuming you meant to use '-L%d'",iRegLen,cMinRegLen,cMinRegLen);
			iRegLen = cMinRegLen;
			}
		else
			{
			if(iRegLen > cMaxRegLen)
				{
				printf("\nRegulatory region length '-L%d' more than maximum %d, assuming you meant to use '-L%d'",iRegLen,cMaxRegLen,cMaxRegLen);
				iRegLen = cMaxRegLen;
				}
			}
		}
	else
		{
		printf("Error: No BED gene/feature file specified");
		exit(1);
		}

	strncpy(szInLociFile,InLociFile->filename[0],_MAX_PATH);
	szInLociFile[_MAX_PATH-1] = '\0';

	if(RsltsFile->count)
		{
		strncpy(szRsltsFile,RsltsFile->filename[0],_MAX_PATH);
		szRsltsFile[_MAX_PATH-1] = '\0';
		}
	else
		szRsltsFile[0] = '\0';

	if(featrsltsfile->count)
		{
		strncpy(szFeatRsltsFile,featrsltsfile->filename[0],_MAX_PATH);
		szFeatRsltsFile[_MAX_PATH-1] = '\0';
		}
	else
		szFeatRsltsFile[0] = '\0';

	if(szFeatRsltsFile[0] == '\0' && szRsltsFile[0] == '\0')
		{
		printf("Error: at least one of '-o<rsltsfile>' or '-O<featrsltsfile>' must be specified");
		exit(1);
		}

	if(summrsltsfile->count)
		{
		strncpy(szSummRsltsFile,summrsltsfile->filename[0],_MAX_PATH);
		szSummRsltsFile[_MAX_PATH-1] = '\0';
		}
	else
		szSummRsltsFile[0] = '\0';

// show user current resource limits
#ifndef _WIN32
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Resources: %s",CUtility::ReportResourceLimits());
#endif

	const char *pszDescr;
	switch(PMode) {
		case ePMdefault:
			pszDescr = "associate features along full length of cores";
			break;
		case ePMstarts:
			pszDescr = "associate features with 5' core starts only";
			break;
		case ePMdyad:
			pszDescr = "associate features with assumed nucleosome dyad positioned 73nt downstream of 5' core start";
			break;

		}
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Processing mode is : '%s'",pszDescr);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Dedupe input elements : '%s'",bDedupe ? "yes" : "no");

	switch(IsoformRprt) {
		case eISOFPRPKM:
			pszDescr = "report the feature isoform with maximal RPKM";
			break;
		case eISOFReads:
			pszDescr = "report the feature isoform with maximal total reads";
			break;
		case eISOFPall:
			pszDescr = "report all feature isoforms";
			break;
		}
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Feature isoform report mode is : '%s'",pszDescr);

	switch(StrandProc) {
		case eStrandDflt:
			pszDescr = "ignore strand";
			break;
		case eStrandSense:
			pszDescr = "process for sense";
			break;
		case eStrandAnti:			
			pszDescr = "process for antisense";
			break;
		}

	
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Feature association is : '%s'",bFeatinsts ? "individual" : "all");
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"One count per read : '%s'",bOneCntRead ? "yes" : "no");
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Strand processing : '%s'",pszDescr);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Input element loci file: '%s'",szInLociFile);
	switch(FType) {
		case 0:
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Auto-classify input element file as either CSV, BED or SAM");
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"If CSV then as: %s",CSVFormat2Text((teCSVFormat)iCSVFormat));
			break;

		case 1:
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Expecting input element file to be CSV format as: %s",CSVFormat2Text((teCSVFormat)iCSVFormat));
			break;

		case 2:
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Expecting input element file to be BED format");
			break;

		case 3:
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Expecting input element file to be SAM format");
			break;
		}
	
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Element mapping results to file: '%s'",szRsltsFile[0]=='\0'?"None Specified":szRsltsFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Feature mapping results to file: '%s'",szFeatRsltsFile[0]=='\0'?"None Specified":szFeatRsltsFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Chrom summary results to file: '%s'",szSummRsltsFile[0]=='\0'?"None Specified":szSummRsltsFile);


	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Input biobed BED gene feature file: '%s'",szInBEDFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Feature mappings with regulatory region length: %d",iRegLen);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Minimum element length: %d",iMinLength);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Maximum element length: %d",iMaxLength);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Join overlap: %d",iJoinOverlap);

	if(szExperimentName[0] != '\0')
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"This processing reference: %s",szExperimentName);

#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	gStopWatch.Start();
	Rslt = MLFProcess((etMLFPMode)PMode,bDedupe,bFeatinsts,bOneCntRead,IsoformRprt,StrandProc,FType,(teCSVFormat)iCSVFormat,szInLociFile,szInBEDFile,szRsltsFile,szFeatRsltsFile,szSummRsltsFile,iRegLen,iMinLength,iMaxLength,iJoinOverlap);
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
}

char *
CSVFormat2Text(teCSVFormat Format)
{
switch(Format) {
	case eCSVFdefault:
		return((char *)"Default 9+ field loci");
	case eCSVFprobe:
		return((char *)"locsfx probe loci");
	case eCSVFtarget:
		return((char *)"locsfx target loci");
	default:
		break;
	}
return((char *)"Unsupported");
}

int MLFProcess(etMLFPMode PMode,				// processing mode
			   bool bDedupe,				// true if input elements are to be deduped
			   bool bFeatinsts,			// true if input elements are to be associated to individual features, false if to all features at that locus
			   bool bOneCntRead,			// true if one count per read rule to be applied (functional regions are prioritised with CDS as the highest) 
			   etISOFProc IsoformRprt,		// feature isoform reporting mode
			   etStrandProc StrandProc,	// how to process read + element strand
			   int Ftype,					// input element file format: 0 - auto, 1 - CSV, 2 - BED, 3 - SAM
			   teCSVFormat CSVFormat,		// expected input CSV loci file format
			   char *pszInLociFile,		// input CSV loci file
			   char *pszInBEDFile,			// input BED file
			   char *pszRsltsFile,			// optional output loci mapping file
			   char *pszFeatRsltsFile,		// optional output feature mapping results file
			   char *pszSummRsltsFile,		// optional output chrom summary results file
			   int RegRegionLen,			// regulatory region length
			   int MinLength,				// minimum element length
			   int MaxLength,				// maximum element length
			   int JoinOverlap)			// deduping join overlap
{
int Rslt;
CMapLoci2Feat *pMapLoci2Feat;
if((pMapLoci2Feat = new CMapLoci2Feat) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "MLFProcess - error Unable to create instance of CMapLoci2Feat");
	return(eBSFerrObj);
	}
Rslt = pMapLoci2Feat->MLFProcess(PMode, bDedupe, bFeatinsts, bOneCntRead, IsoformRprt, StrandProc, Ftype, CSVFormat, pszInLociFile, pszInBEDFile, pszRsltsFile,
						  pszFeatRsltsFile, pszSummRsltsFile, RegRegionLen, MinLength, MaxLength, JoinOverlap);
delete pMapLoci2Feat;
return(Rslt);
}


