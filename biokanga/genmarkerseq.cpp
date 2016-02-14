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

#include "./biokanga.h"
#include "./MarkerSeq.h"

int
 ProcessMarkerSeqs(int PMode,					// currently default processing only is supported
					int Ftype,					// input element file format: 0 - auto, 1 - CSV, 2 - BED, 3 - SAM
					teCSVFormat CSVFormat,		// expected input CSV loci file format
					int MaxNs,					// filter out marker sequences having higher than this number of indeterminate bases (default is 0, range 0..5)
					int Extend5,				// extend markers 5' bases from centroid of loci (default is 25, range 0..250)
					int Extend3,				// extend markers 3' bases from centroid of loci (default is 25, range 0..250)
					int MinSeqLen,				// filter out marker sequences which are less than this length (default is 10bp, range 10..1000)
					int MaxSeqLen,				// filter out marker sequences which are longer than this length (default is 1000bp, range minseqlen..1000)
					bool bNoOverlaps,			// filter out marker sequences which overlap with other marker sequences
					char *pszInLociFile,		// Loci file specifying the SNP or region loci in assembly from which marker sequences are to be generated (CSV, BED or SAM)
					char *pszInFastaFile,		// multifasta assembly file from which marker sequences are to be generated
					char *pszOutFile);			// marker sequences written to this file
 
#ifdef _WIN32
int GenMarkerSeqs(int argc, char* argv[])
{
// determine my process name
_splitpath(argv[0],NULL,NULL,gszProcName,NULL);
#else
int
GenMarkerSeqs(int argc, char** argv)
{
// determine my process name
CUtility::splitpath((char *)argv[0],NULL,gszProcName);
#endif
int iScreenLogLevel;		// level of screen diagnostics
int iFileLogLevel;			// level of file diagnostics
char szLogFile[_MAX_PATH];	// write diagnostics to this file
int Rslt = 0;   			// function result code >= 0 represents success, < 0 on failure

int PMode;					// currently default processing only is supported
int CSVFormat;				// input CSV format
int FType;					// expected input element file type - auto, CSV, BED or SAM
int MaxNs;					// filter out input sequences having higher than this number of indeterminate bases per 100bp (default is 0, range 0..5)
int Extend5;				// extend marker 5' bases from centre of loci (default is 50, range 5..250)
int Extend3;				// extend marker 3' bases from centre of loci (default is 50, range 5..250)
int MinSeqLen;				// filter out marker sequences which are less than this length (default is 50bp, range 10..1000)
int MaxSeqLen;				// trim marker sequences to be no longer than this length (default is 0 for no length trimming, 10..1000)
bool bNoOverlaps;			// filter out marker sequences which overlap with other marker sequences

char szInLociFile[_MAX_PATH];  // Loci file specifying the SNP or region loci in assembly from which marker sequences are to be generated (CSV, BED or SAM)
char szInFastaFile[_MAX_PATH];		// multifasta assembly file from which marker sequences are to be generated
char szOutFile[_MAX_PATH];	// marker sequences written to this file

char szSQLiteDatabase[_MAX_PATH];	// results summaries to this SQLite file
char szExperimentName[cMaxDatasetSpeciesChrom+1];			// experiment name
char szExperimentDescr[1000];		// describes experiment

// command line args
struct arg_lit  *help    = arg_lit0("h","help",                 "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to screen and logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"diagnostics log file");

struct arg_int *pmode = arg_int0("m","mode","<int>",		    "processing mode: 0 default processing");
struct arg_int *ftype = arg_int0("t","filetype","<int>",		"input element file format: 0 - auto, 1 - CSV, 2 - BED, 3 - SAM (default = 0)");
struct arg_int  *csvformat = arg_int0("c","informat","<int>",	"if input CSV file processing 0:Loci 1:locsfx probe 2:locsfx target");
struct arg_file *infastafile = arg_file1("I","infasta","<file>",	"Marker sequences from this multifasta assembly or sequences file");
struct arg_file *inlocifile = arg_file1("i","inloci","<file>",	"Loci file specifying the SNP or region loci in assembly from which marker sequences are to be generated (CSV, BED)");
struct arg_file *outfile = arg_file0("o","out","<file>",		"Output multifasta marker sequences to this file");
struct arg_int *extend5 = arg_int0("x","extend5","<int>",		"extend markers 5' bases from centroid of loci (default is 25, range 0..250)");
struct arg_int *extend3 = arg_int0("X","extend3","<int>",		"extend markers 3' bases from centroid of loci (default is 25, range 0..250)");
struct arg_int *maxns = arg_int0("n","indeterminates","<int>",  "filter out marker sequences having higher than this number of indeterminate bases (default is 0, range 0..5)");
struct arg_int *minseqlen= arg_int0("l","minlen","<int>",       "filter out marker sequences which are less than this length (default is 10bp, range 10..1000)");
struct arg_int *maxseqlen = arg_int0("L","maxlen","<int>",		"filter out marker sequences which are longer than this length (default is 1000bp, range minlen..1000)");
struct arg_lit  *nooverlaps    = arg_lit0("O","nooverlaps",     "filter out marker sequences with 5' flank extensions which overlap other marker 3' flank extensions");

struct arg_file *summrslts = arg_file0("q","sumrslts","<file>",		"Output results summary to this SQLite3 database file");
struct arg_str *experimentname = arg_str0("w","experimentname","<str>",		"experiment name SQLite3 database file");
struct arg_str *experimentdescr = arg_str0("W","experimentdescr","<str>",	"experiment description SQLite3 database file");


struct arg_end *end = arg_end(200);

void *argtable[] = {help,version,FileLogLevel,LogFile,
	                pmode,ftype,csvformat,maxns,extend5,extend3,minseqlen,maxseqlen,nooverlaps,infastafile,inlocifile,outfile,
					summrslts,experimentname,experimentdescr,
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


	PMode = 1;

	if(!inlocifile->count)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Expected input loci file '-i<locifile>' to be specified");
		return(1);
		}
	strncpy(szInLociFile,inlocifile->filename[0],_MAX_PATH);
	szInLociFile[_MAX_PATH-1] = '\0';

	if(!infastafile->count)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Expected input multifasta assembly or sequences file '-I<infasta>' to be specified");
		return(1);
		}
	strncpy(szInFastaFile,infastafile->filename[0],_MAX_PATH);
	szInFastaFile[_MAX_PATH-1] = '\0';

	if(!outfile->count)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Expected output multifasta markers file '-o<outfile>' to be specified");
		return(1);
		}
	strncpy(szOutFile,outfile->filename[0],_MAX_PATH);
	szOutFile[_MAX_PATH-1] = '\0';

	FType = ftype->count ? ftype->ival[0] : 0;
	if(FType < 0 || FType >= 4)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Expected input element file format '-t%d' specified outside of range %d..%d",FType,0,3);
		exit(1);
		}

	bNoOverlaps =	nooverlaps->count ? true : false;

	if(FType <= 1)
		{
		CSVFormat = csvformat->count ? csvformat->ival[0] : 0;
		if(CSVFormat < eCSVFdefault || CSVFormat > eCSVFtarget)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: expected input CSV format specified '-c%d' must be in range 0..2",CSVFormat);
			exit(1);
			}
		}
	else
		CSVFormat = 0;

	Extend5 = extend5->count ? extend5->ival[0] : 25;
	if(Extend5 < 0 || Extend5 > 250)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: 5' extend '-x%d' must be in range 0..250 bases",Extend5);
		return(1);
		}

	Extend3 = extend3->count ? extend3->ival[0] : 25;
	if(Extend3 < 0 || Extend3 > 250)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: 3' extend '-X%d' must be in range 0..250 bases",Extend3);
		return(1);
		}

	MinSeqLen = minseqlen->count ? minseqlen->ival[0] : 10;
	if(MinSeqLen < 10 || MinSeqLen > 1000)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Minimum sequence length '-l%d' must be in range 10..1000 bases",MinSeqLen);
		return(1);
		}


	MaxSeqLen = maxseqlen->count ? maxseqlen->ival[0] : 1000;
	if(MaxSeqLen < MinSeqLen || MaxSeqLen > 1000)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: trim sequence length '-L%d' must be in range %d..1000 bases",MaxSeqLen,MinSeqLen);
		return(1);
		}

	MaxNs = maxns->count ? maxns->ival[0] : 0;
	if(MaxNs < 0 || MaxNs > 5)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: max percentage of indeterminate bases '-n%d' must be in range 0..5",MaxNs);
		return(1);
		}

// show user current resource limits
#ifndef _WIN32
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Resources: %s",CUtility::ReportResourceLimits());
#endif

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing parameters:");
	const char *pszDescr = "Generate marker sequences";

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Processing mode is : '%s'",pszDescr);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"CSV format : %d",CSVFormat);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Input loci file type : %d",FType);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Extend 5' : %d",Extend5);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Extend 3' : %d",Extend3);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Maximum indeterminant 'N's : %d",MaxNs);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Minimum marker sequence length : %d",MinSeqLen);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Maximum marker sequence length length : %d",MaxSeqLen);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Filter out marker sequences with 5' flank extensions which overlap other marker 3' flank extensions: %s",bNoOverlaps == true ? "Yes" : "No");


	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Input loci file : '%s'",szInLociFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Input assembly file : '%s'",szInFastaFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Output multifasta file : '%s'",szOutFile);
	
	if(szExperimentName[0] != '\0')
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"This processing reference: %s",szExperimentName);

	if(gExperimentID > 0)
		{
		int ParamID;
		int NoOverlaps;
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szLogFile),"log",szLogFile);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,sizeof(PMode),"mode",&PMode);

		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,sizeof(CSVFormat),"informat",&CSVFormat);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,sizeof(FType),"filetype",&FType);

		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,sizeof(MaxNs),"indeterminates",&MaxNs);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,sizeof(Extend5),"extend5",&Extend5);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,sizeof(Extend3),"extend3",&Extend3);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,sizeof(MaxSeqLen),"maxlen",&MaxSeqLen);
		NoOverlaps = bNoOverlaps ? 1 : 0;
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,sizeof(NoOverlaps),"nooverlaps",&bNoOverlaps);

		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,sizeof(MaxNs),"indeterminates",&MaxNs);

		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szInFastaFile),"infasta",szInFastaFile);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szInLociFile),"inloci",szInLociFile);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szOutFile),"out",szOutFile);

		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szSQLiteDatabase),"sumrslts",szSQLiteDatabase);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szExperimentName),"experimentname",szExperimentName);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szExperimentDescr),"experimentdescr",szExperimentDescr);
		}


#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	gStopWatch.Start();
	Rslt = ProcessMarkerSeqs(PMode,FType,(teCSVFormat)CSVFormat,MaxNs,Extend5,Extend3,MinSeqLen,MaxSeqLen,bNoOverlaps,szInLociFile,szInFastaFile,szOutFile);
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
 ProcessMarkerSeqs(int PMode,		// currently default processing only is supported
					int Ftype,					// input element file format: 0 - auto, 1 - CSV, 2 - BED, 3 - SAM
					teCSVFormat CSVFormat,		// expected input CSV loci file format
					int MaxNs,					// filter out marker sequences having higher than this number of indeterminate bases (default is 0, range 0..5)
					int Extend5,				// extend markers 5' bases from centroid of loci (default is 25, range 0..250)
					int Extend3,				// extend markers 3' bases from centroid of loci (default is 25, range 0..250)
					int MinSeqLen,				// filter out marker sequences which are less than this length (default is 10bp, range 10..1000)
					int MaxSeqLen,				// filter out marker sequences which are longer than this length (default is 1000bp, range minseqlen..1000)
					bool bNoOverlaps,			// filter out marker sequences which overlap with other marker sequences
					char *pszInLociFile,		// Loci file specifying the SNP or region loci in assembly from which marker sequences are to be generated (CSV, BED or SAM)
					char *pszInFastaFile,		// multifasta assembly file from which marker sequences are to be generated
					char *pszOutFile)			// marker sequences written to this file
{
int Rslt;
CMarkerSeq MarkerSeq;

Rslt = MarkerSeq.ProcessMarkerSeqs(PMode,Ftype,CSVFormat,MaxNs,Extend5,Extend3,MinSeqLen,MaxSeqLen,bNoOverlaps,pszInLociFile,pszInFastaFile,pszOutFile);

MarkerSeq.Reset();
return(Rslt);
}