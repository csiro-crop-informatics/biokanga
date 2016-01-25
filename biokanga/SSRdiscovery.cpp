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

// SSRdiscovery.cpp : Defines the entry point for the console application.
// Purpose -
// Loads multifasta files and discovers/reports on all SSRs contained within the fasta sequences

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

#include "SSRdiscovery.h"
#include "SQLiteSummaries.h"

#include "biokanga.h"


const int cMaxInFileSpecs = 100;			// allow at most this number of input files

int
Process(int PMode,					// processing mode
		teRptSSRsFromat RptSSRsFormat,	// report SSRs in this file format 
		int MinRepElLen,			// identify repeating elements of this minimum length
		int MaxRepElLen,			// ranging upto this maximum length
		int MinTandemRpts,			// minimum number of tandem repeats
		int MaxTandemRpts,			// maximum number of repeats
		int SSRFlankLen,			// SSR flanking sequence length
		int NumInputFiles,			// number of input filespecs
		char *pszInputFiles[],		// input multifasta files
		char *pszKMerFreqFile,		// optional, output element KMer freq to this file
		char *pszOutFile);			// output SSRs to this file


#ifdef _WIN32
int SSRdiscovery(int argc, char* argv[])
{
// determine my process name
_splitpath(argv[0],NULL,NULL,gszProcName,NULL);
#else
int SSRdiscovery(int argc, char** argv)
{
// determine my process name
CUtility::splitpath((char *)argv[0],NULL,gszProcName);
#endif

int iScreenLogLevel;		// level of screen diagnostics
int iFileLogLevel;			// level of file diagnostics
char szLogFile[_MAX_PATH];	// write diagnostics to this file
int Rslt = 0;   			// function result code >= 0 represents success, < 0 on failure

int PMode;					// processing mode

int MinRepElLen;			// identify repeating elements of this minimum length
int MaxRepElLen;			// ranging upto this maximum length
int MinTandemRpts;			// minimum number of tandem repeats
int MaxTandemRpts;			// maximum number of repeats
int SSRFlankLen;			// SSR flanking sequence length

int NumInputFiles;			// number of intput filespecs
char *pszInputFiles[cMaxInFileSpecs];	// names of input files (wildcards allowed) containing sequences to be processed for SSRs

char szKMerFreqFile[_MAX_PATH];	// optional, output element KMer freq to this file
char szOutFile[_MAX_PATH];		// write identified SSRs to this file

char szSQLiteDatabase[_MAX_PATH];	// results summaries to this SQLite file
char szExperimentName[cMaxDatasetSpeciesChrom+1];			// experiment name
char szExperimentDescr[1000];		// describes experiment

// command line args
struct arg_lit  *help    = arg_lit0("h","help",                 "Print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"Print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to screen and logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"Diagnostics log file");

struct arg_int *mode = arg_int0("m","mode","<int>",								"Processing mode: 0 default");

struct arg_int *minrepellen = arg_int0("k","minrepellen","<int>",				"identify repeating K-mer elements of this minimum length (2 default)");
struct arg_int *maxrepellen = arg_int0("K","maxrepellen","<int>",				"identify repeating K-mer elements of this maximum length (5 default)");
struct arg_int *mintandemrpts = arg_int0("r","mintandemrpts","<int>",			"minimum number of tandem element repeats (5 default)");
struct arg_int *maxtandemrpts = arg_int0("R","maxtandemrpts","<int>",			"maximum number of tandem element repeats (10 default)");

struct arg_int *flanklen = arg_int0("l","flanklen","<int>",						"report SSR flanking lengths (100bp default)");

struct arg_file *inputfiles = arg_filen("i","in","<file>",0,cMaxInFileSpecs,	"Input file(s) containing sequences to process for SSRs");
struct arg_file *kmerfreq = arg_file0("O","outkmerfreq","<file>",				"Output K-mer element freq to this file");
struct arg_file *outfile = arg_file1("o","out","<file>",						"Output SSRs to this file");

struct arg_file *summrslts = arg_file0("q","sumrslts","<file>",					"Output results summary to this SQLite3 database file");
struct arg_str *experimentname = arg_str0("w","experimentname","<str>",			"experiment name SQLite3 database file");
struct arg_str *experimentdescr = arg_str0("W","experimentdescr","<str>",		"experiment description SQLite3 database file");

struct arg_end *end = arg_end(200);

void *argtable[] = {help,version,FileLogLevel,LogFile,
					summrslts,experimentname,experimentdescr,
	                mode,minrepellen, maxrepellen, mintandemrpts, maxtandemrpts, flanklen, inputfiles, kmerfreq, outfile,
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
		gProcessID = gSQLiteSummaries.AddProcess((char *)gszProcName,(char *)gszProcName,(char *)szExperimentDescr);
		if(gProcessID < 1)
			return(1);
		gProcessingID = gSQLiteSummaries.StartProcessing(gExperimentID,gProcessID,(char *)cpszProgVer);
		if(gProcessingID < 1)
			return(1);
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Initialised SQLite database '%s' for SSR results summary",szSQLiteDatabase);
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"SQLite database experiment identifier for '%s' is %d",szExperimentName,gExperimentID);
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"SQLite database process identifier for '%s' is %d",(char *)gszProcName,gProcessID);
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"SQLite database processing instance identifier is %d",gProcessingID);
		}
	else
		{
		szSQLiteDatabase[0] = '\0';
		szExperimentDescr[0] = '\0';
		}

	PMode = mode->count ? mode->ival[0] : 0;
	if(PMode < 0 || PMode > 1)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Processing mode '-m%d' must be in range 0..0",PMode);
		return(1);
		}


	MinRepElLen = minrepellen->count ? minrepellen->ival[0] : cDfltMinRepElLen;
	if(MinRepElLen < cMinRepElLen || MinRepElLen > cMaxRepElLen)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Repeating K-mer element minimum length '-k%d' must be in range %d..%d",MinRepElLen,cMinRepElLen,cMaxRepElLen);
		return(1);
		}

	MaxRepElLen = maxrepellen->count ? maxrepellen->ival[0] : max(cDfltMaxRepElLen,MinRepElLen);
	if(MaxRepElLen < MinRepElLen || MaxRepElLen > cMaxRepElLen)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Repeating K-mer element maximum length '-K%d' must be in range %d..%d",MaxRepElLen,MinRepElLen,cMaxRepElLen);
		return(1);
		}

	MinTandemRpts = mintandemrpts->count ? mintandemrpts->ival[0] : cDfltMinTandemRpts;
	if(MinTandemRpts < cMinTandemRpts || MinTandemRpts > cMaxTandemRpts)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Repeating K-mer element tandem min repeats '-r%d' must be in range %d..%d",MinTandemRpts,cMinTandemRpts,cMaxTandemRpts);
		return(1);
		}

	MaxTandemRpts = maxtandemrpts->count ? maxtandemrpts->ival[0] : max(cDfltMaxTandemRpts,MinTandemRpts);
	if(MaxTandemRpts < MinTandemRpts || MaxTandemRpts > cMaxTandemRpts)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Repeating K-mer element tandem max repeats '-R%d' must be in range %d..%d",MaxTandemRpts,MinTandemRpts,cMaxTandemRpts);
		return(1);
		}

	SSRFlankLen = flanklen->count ? flanklen->ival[0] : cDfltSSRFlankLen;
	if(SSRFlankLen < cMinSSRFlankLen || SSRFlankLen > cMaxSSRFlankLen)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error:SSR flanking length '-l%d' must be in range %d..%d",SSRFlankLen,cMinSSRFlankLen,cMaxSSRFlankLen);
		return(1);
		}

	int Idx;
	for(NumInputFiles=Idx=0;NumInputFiles < cMaxInFileSpecs && Idx < inputfiles->count; Idx++)
		{
		pszInputFiles[Idx] = NULL;
		if(pszInputFiles[NumInputFiles] == NULL)
			pszInputFiles[NumInputFiles] = new char [_MAX_PATH];
		strncpy(pszInputFiles[NumInputFiles],inputfiles->filename[Idx],_MAX_PATH);
		pszInputFiles[NumInputFiles][_MAX_PATH-1] = '\0';
		CUtility::TrimQuotedWhitespcExtd(pszInputFiles[NumInputFiles]);
		if(pszInputFiles[NumInputFiles][0] != '\0')
			NumInputFiles++;
		}

	if(!NumInputFiles)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: After removal of whitespace, no input file(s) specified with '-i<filespec>' option)\n");
		exit(1);
		}

	if(kmerfreq->count)
		{
		if(!(MinRepElLen == MaxRepElLen && MinRepElLen <= 10))
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: SSR element frequency counting not supported for range of element lengths or if element length > 10\n");
			exit(1);
			}
		strcpy(szKMerFreqFile,kmerfreq->filename[0]);
		CUtility::TrimQuotedWhitespcExtd(szKMerFreqFile);
		}
	else
		szKMerFreqFile[0] = '\0';

	strcpy(szOutFile,outfile->filename[0]);
	CUtility::TrimQuotedWhitespcExtd(szOutFile);


	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing parameters:");
	const char *pszDescr;
	switch(PMode) {
		case 0:
			pszDescr = "Identify SSRs";
			break;
		}

// show user current resource limits
#ifndef _WIN32
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Resources: %s",CUtility::ReportResourceLimits());
#endif

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Processing mode is : '%s'",pszDescr);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"identify repeating K-mer elements of this minimum length: %d",MinRepElLen);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"identify repeating K-mer elements of this maximum length: %d",MaxRepElLen);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"minimum number of tandem element repeats: %d",MinTandemRpts);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"maximum number of tandem element repeats: %d",MaxTandemRpts);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"SSR flanking length: %dbp",SSRFlankLen);


	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Experiment name : '%s'",szExperimentName);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Experiment description : '%s'",szExperimentDescr);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"SQLite database file: '%s'",szSQLiteDatabase);

	for(Idx = 0; Idx < NumInputFiles; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Input file (%d): '%s'",Idx+1,pszInputFiles[Idx]);
	if(szKMerFreqFile[0] != '\0')
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"SSR element K-Mer freq file: '%s'",szKMerFreqFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"SSRs to file: '%s'",szOutFile);

	if(gExperimentID > 0)
		{
		int ParamID;
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szLogFile),"log",szLogFile);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,sizeof(PMode),"mode",&PMode);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,sizeof(MinRepElLen),"minrepellen",&MinRepElLen);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,sizeof(MaxRepElLen),"maxrepellen",&MaxRepElLen);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,sizeof(MinTandemRpts),"mintandemrpts",&MinTandemRpts);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,sizeof(MaxTandemRpts),"maxtandemrpts",&MaxTandemRpts);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,sizeof(SSRFlankLen),"flanklen",&SSRFlankLen);

	    ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,sizeof(NumInputFiles),"NumInputFiles",&NumInputFiles);
		for(Idx=0; Idx < NumInputFiles; Idx++)
			ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(pszInputFiles[Idx]),"in",pszInputFiles[Idx]);

		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szOutFile),"out",szOutFile);
		if(szKMerFreqFile[0] != '\0')
			ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szKMerFreqFile),"outkmerfreq",szKMerFreqFile);

		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szSQLiteDatabase),"sumrslts",szSQLiteDatabase);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szExperimentName),"experimentname",szExperimentName);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szExperimentDescr),"experimentdescr",szExperimentDescr);
		}


#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	gStopWatch.Start();

	Rslt = Process(PMode,eRFCsv,MinRepElLen,MaxRepElLen,MinTandemRpts,MaxTandemRpts,SSRFlankLen,NumInputFiles,pszInputFiles,szKMerFreqFile,szOutFile);
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
Process(int PMode,					// processing mode
		teRptSSRsFromat RptSSRsFormat,	// report SSRs in this file format 
		int MinRepElLen,			// identify repeating elements of this minimum length
		int MaxRepElLen,			// ranging upto this maximum length
		int MinTandemRpts,			// minimum number of tandem repeats
		int MaxTandemRpts,			// maximum number of repeats
		int SSRFlankLen,			// SSR flanking lengths
		int NumInputFiles,			// number of input filespecs
		char *pszInputFiles[],		// input multifasta files
		char *pszKMerFreqFile,		// optional, output element KMer freq to this file
		char *pszOutFile)			// output SSRs to this file
{
int Rslt;
CSSRDiscovery CSSRDiscovery;
Rslt = CSSRDiscovery.Process(PMode,RptSSRsFormat,MinRepElLen,MaxRepElLen,MinTandemRpts,MaxTandemRpts,SSRFlankLen,NumInputFiles,pszInputFiles,pszKMerFreqFile,pszOutFile);
return(Rslt);
}




CSSRDiscovery::CSSRDiscovery(void)
{
m_pSeqBuff = NULL;
m_pKMerDist = NULL;
m_pszRptSSRsBuff = NULL;
Init();
}

CSSRDiscovery::~CSSRDiscovery(void)
{
if(m_pSeqBuff != NULL)
	{
#ifdef _WIN32
	free(m_pSeqBuff);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pSeqBuff != MAP_FAILED)
		munmap(m_pSeqBuff,m_AllocdSeqBuffMem);
#endif
	m_pSeqBuff = NULL;
	}
if(m_pKMerDist != NULL)
	{
#ifdef _WIN32
	free(m_pKMerDist);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pKMerDist != MAP_FAILED)
		munmap(m_pKMerDist,m_AllocdKMerFreqMem);
#endif
	m_pKMerDist = NULL;
	}
if(m_pszRptSSRsBuff != NULL)
	delete m_pszRptSSRsBuff;
}


void
CSSRDiscovery::Init(void)
{
m_pSeqBuff = NULL;
m_pKMerDist = NULL;
m_pszRptSSRsBuff = NULL;
m_hOutFile = -1;
m_hOutKMerFreqFile = -1;
Reset();
m_CurTime.Start();
}

void
CSSRDiscovery::Reset(void)
{
if(m_hOutFile != -1)
	{
#ifdef _WIN32
	_commit(m_hOutFile);
#else
	fsync(m_hOutFile);
#endif
	close(m_hOutFile);
	m_hOutFile = -1;
	}

if(m_hOutKMerFreqFile != -1)
	{
#ifdef _WIN32
	_commit(m_hOutKMerFreqFile);
#else
	fsync(m_hOutKMerFreqFile);
#endif
	close(m_hOutKMerFreqFile);
	m_hOutKMerFreqFile = -1;
	}

if(m_pszRptSSRsBuff != NULL)
	{
	delete m_pszRptSSRsBuff;
	m_pszRptSSRsBuff = NULL;
	}

if(m_pSeqBuff != NULL)
	{
#ifdef _WIN32
	free(m_pSeqBuff);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pSeqBuff != MAP_FAILED)
		munmap(m_pSeqBuff,m_AllocdSeqBuffMem);
#endif
	m_pSeqBuff = NULL;
	}
if(m_pKMerDist != NULL)
	{
#ifdef _WIN32
	free(m_pKMerDist);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pKMerDist != MAP_FAILED)
		munmap(m_pKMerDist,m_AllocdKMerFreqMem);
#endif
	m_pKMerDist = NULL;
	}

m_AllocdSeqBuffMem = 0;
m_SeqBuffLen = 0;
m_AllocdKMerFreqMem = 0;
m_KMerFreqLen = 0;
m_IdxRptSSRs = 0;
m_CurTime.Stop();
}


etSeqBase *
CSSRDiscovery::AllocSeqBuff(size_t SeqLen)				// allocate for at least this sequence length
{
size_t memreq;
etSeqBase *pTmp;

if(m_pSeqBuff != NULL && m_AllocdSeqBuffMem >= SeqLen)
	return(m_pSeqBuff);

if(m_pSeqBuff == NULL)
	{
	memreq = max(SeqLen,(size_t)cMaxAllocBuffChunk);
#ifdef _WIN32
	m_pSeqBuff = (etSeqBase *) malloc(SeqLen);	// initial and perhaps the only allocation
	if(m_pSeqBuff == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AllocSeqBuff: Memory allocation of %lld bytes - %s",(INT64)SeqLen,strerror(errno));
		return(NULL);
		}
#else
	// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
	m_pSeqBuff = (etSeqBase *)mmap(NULL,memreq, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
	if(m_pSeqBuff == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AllocSeqBuff: Memory allocation of %lld bytes through mmap()  failed - %s",(INT64)memreq,strerror(errno));
		m_pSeqBuff = NULL;
		return(NULL);
		}
#endif
	}
else
	{
	memreq = SeqLen + cMaxAllocBuffChunk;
#ifdef _WIN32
	pTmp = (etSeqBase *) realloc(m_pSeqBuff,memreq);
#else
	pTmp = (etSeqBase *)mremap(m_pSeqBuff,m_AllocdSeqBuffMem,memreq,MREMAP_MAYMOVE);
	if(pTmp == MAP_FAILED)
		pTmp = NULL;
#endif
	if(pTmp == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AllocSeqBuff: Memory re-allocation to %lld bytes - %s",(INT64)memreq,strerror(errno));
		return(NULL);
		}
	m_pSeqBuff = pTmp;
	}
m_AllocdSeqBuffMem = memreq;
return(m_pSeqBuff);
}

// ProcessBioseqFile
// Process input biosequence file
int
CSSRDiscovery::ProcessBioseqFile(char *pszFile)
{
CBioSeqFile BioSeqFile;
char szSource[cBSFSourceSize];
char szDescription[cBSFDescriptionSize];
UINT32 SeqLen;
int Rslt;
tBSFEntryID CurEntryID;

if((Rslt=BioSeqFile.Open(pszFile,cBSFTypeSeq,false))!=eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessBioseqFile - error %d %s",Rslt,BioSeqFile.GetErrMsg());
	return(Rslt);
	}

CurEntryID = 0;
while((Rslt = CurEntryID = BioSeqFile.Next(CurEntryID)) > eBSFSuccess)
	{
	BioSeqFile.GetNameDescription(CurEntryID,cBSFSourceSize-1,(char *)&szSource,
											cBSFDescriptionSize-1,(char *)&szDescription);
	SeqLen = BioSeqFile.GetDataLen(CurEntryID);
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing %s|%s",szSource,szDescription);

	if(!SeqLen)
		continue;

	if(AllocSeqBuff(SeqLen) == NULL)
		{
		Rslt = eBSFerrMem;
		break;
		}

	if((Rslt = BioSeqFile.GetData(CurEntryID,eSeqBaseType,0,m_pSeqBuff,SeqLen)) != SeqLen)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessBioseqFile - error %d %s",Rslt,BioSeqFile.GetErrMsg());
		break;
		}

	if((Rslt=IdentifySSRs(szSource,pszFile,SeqLen,m_pSeqBuff)) < eBSFSuccess)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessBioseqFile - error %d",Rslt);
		break;
		}

	ReportProgress();
	}
if(Rslt == eBSFerrEntry)
	Rslt = eBSFSuccess;

BioSeqFile.Close();
return(Rslt);
}

// ProcessFastaFile
// Parse input fasta format file
int
CSSRDiscovery::ProcessFastaFile(char *pszFile)
{
CFasta Fasta;

size_t AvailBuffSize;
char szName[cBSFSourceSize];
char szDescription[cBSFDescriptionSize];
UINT32 SeqLen;
int Descrlen;
bool bFirstEntry;
bool bEntryCreated;
int Rslt;
int SeqID;

if((Rslt=Fasta.Open(pszFile,true))!=eBSFSuccess)
	{
	if(Rslt != eBSFerrNotFasta)
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessFastaFile: Unable to open '%s' [%s] %s",pszFile,Fasta.ErrText((teBSFrsltCodes)Rslt),Fasta.GetErrMsg());
	return(Rslt);
	}

if(m_pSeqBuff == NULL)				// if not already allocated then allocate to hold cMaxAllocBuffChunk bases 
	{
	SeqLen = cMaxAllocBuffChunk;
	if(AllocSeqBuff(SeqLen) == NULL)
		{
		Rslt = eBSFerrMem;
		Fasta.Close();
		return(Rslt);
		}
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"ProcessFastaFile:- Adding %s..",pszFile);

bFirstEntry = true;
bEntryCreated = false;
SeqID = 0;
m_SeqBuffLen = 0;
AvailBuffSize = m_AllocdSeqBuffMem;
while((Rslt = SeqLen = Fasta.ReadSequence(&m_pSeqBuff[m_SeqBuffLen],(int)min(AvailBuffSize,(size_t)cMaxAllocBuffChunk),true,false)) > eBSFSuccess)
	{
	if(SeqLen == eBSFFastaDescr)		// just read a descriptor line
		{
		SeqID++;
		if(bEntryCreated)				// add any previous entry
			{
			if((Rslt=IdentifySSRs(szName,pszFile,m_SeqBuffLen,m_pSeqBuff)) < eBSFSuccess)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessFastaFile - error %d",Rslt);
				break;
				}
			}
		Descrlen = Fasta.ReadDescriptor(szDescription,cBSFDescriptionSize);
		// An assumption - will one day bite real hard - is that the
		// fasta descriptor line starts with some form of unique identifier.
		// Use this identifier as the entry name.
		if(sscanf(szDescription," %s[ ,]",szName)!=1)
			sprintf(szName,"%s.%d",pszFile,++SeqID);

		bFirstEntry = false;
		bEntryCreated = true;
		m_SeqBuffLen = 0;
		continue;
		}
	else
		if(bFirstEntry)	// if there was no descriptor then dummy up one...
			{
			SeqID++;
			sprintf(szName,"%s.%d",pszFile,SeqID);
			strcpy(szDescription,"No Description provided");
			bFirstEntry = false;
			bEntryCreated = true;
			}

	m_SeqBuffLen += SeqLen;
	AvailBuffSize -= SeqLen;
	if(AvailBuffSize < (size_t)(cMaxAllocBuffChunk / 8))
		{
		if(AllocSeqBuff(m_AllocdSeqBuffMem + SeqLen) == NULL)
			{
			Rslt = eBSFerrMem;
			Fasta.Close();
			return(Rslt);
			}
		AvailBuffSize = m_AllocdSeqBuffMem - m_SeqBuffLen;
		}
	}

if(Rslt >= eBSFSuccess && bEntryCreated && m_SeqBuffLen > 0)			// close entry
	if((Rslt=IdentifySSRs(szName,pszFile,m_SeqBuffLen,m_pSeqBuff)) < eBSFSuccess)
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"IdentifySSRs - error %d",Rslt);
	else
		Rslt = eBSFSuccess;

return(Rslt);
}


int
CSSRDiscovery::ReportCSV(int RepElLen,	// identified SSR contains repeat elements of this length
	   int NumTandemEls,			// each repeat element is repeated this many times
	   INT64 SSRStartOfs,			// repeat element starts at this offset within the targeted sequence
	    char *pszDescr,				// descriptor for the targeted sequence
		char *pszInFile,			// sequence parsed from this file
	   INT64 TargSeqLen,			// targeted sequence is this length
	   etSeqBase *pTargSeq)			// targeted sequence within which the SSR has been located
{
static bool bFirst = true;
int Flank5Len;
INT64 Flank5Ofs;
int Flank3Len;
INT64 Flank3Ofs;

if(SSRStartOfs < m_SSRFlankLen)
	{
	Flank5Len = (int)SSRStartOfs;
	Flank5Ofs = 0;
	}
else
	{
	Flank5Len = m_SSRFlankLen;
	Flank5Ofs = SSRStartOfs - m_SSRFlankLen;
	}

Flank3Ofs = SSRStartOfs + (NumTandemEls * RepElLen);

if((Flank3Ofs +  m_SSRFlankLen) > TargSeqLen)
	Flank3Len = (int)(TargSeqLen - Flank3Ofs);
else
	Flank3Len = m_SSRFlankLen;

if(bFirst)
	{
	m_IdxRptSSRs = sprintf(m_pszRptSSRsBuff,"\"ID\",\"Proc\",\"Species\",\"Chrom\",\"SeqStart\",\"SSRStart\",\"SSREnd\",\"SSRLen\",\"SeqEnd\",\"Strand\",\"KMer\",\"Rpts\",\"Seq\"\n");
	bFirst = false;
	}
if(m_IdxRptSSRs > (cMaxAllocRptSSRs - 2000))
	{
	CUtility::SafeWrite(m_hOutFile,m_pszRptSSRsBuff,m_IdxRptSSRs);
	m_IdxRptSSRs = 0;
	}

m_IdxRptSSRs += sprintf(&m_pszRptSSRsBuff[m_IdxRptSSRs],"%d,\"SSRs\",\"N/A\",\"%s\",%lld,%lld,%lld,%d,%lld,\"+\",%d,%d,\"",m_TotNumAcceptedSSRs,
						pszDescr,Flank5Ofs,SSRStartOfs,SSRStartOfs+(RepElLen*NumTandemEls)-1,RepElLen*NumTandemEls,Flank3Ofs + Flank3Len - 1,RepElLen,NumTandemEls);

if(Flank5Len > 0)
	{
	CSeqTrans::MapSeq2LCAscii(&pTargSeq[Flank5Ofs],Flank5Len,&m_pszRptSSRsBuff[m_IdxRptSSRs]);
	m_IdxRptSSRs += Flank5Len;
	}
CSeqTrans::MapSeq2UCAscii(&pTargSeq[SSRStartOfs],RepElLen * NumTandemEls,&m_pszRptSSRsBuff[m_IdxRptSSRs]);
m_IdxRptSSRs += RepElLen * NumTandemEls;

if(Flank3Len > 0)
	{
	CSeqTrans::MapSeq2LCAscii(&pTargSeq[Flank3Ofs],Flank3Len,&m_pszRptSSRsBuff[m_IdxRptSSRs]);
	m_IdxRptSSRs += Flank3Len;
	}
m_IdxRptSSRs += sprintf(&m_pszRptSSRsBuff[m_IdxRptSSRs],"\"\n");
return(m_IdxRptSSRs);
}

int
CSSRDiscovery::ReportBED(int RepElLen,	// identified SSR contains repeat elements of this length
	   int NumTandemEls,			// each repeat element is repeated this many times
	   INT64 SSRStartOfs,			// repeat element starts at this offset within the targeted sequence
	    char *pszDescr,				// descriptor for the targeted sequence
		char *pszInFile,			// sequence parsed from this file
	   INT64 TargSeqLen,			// targeted sequence is this length
	   etSeqBase *pTargSeq)			// targeted sequence within which the SSR has been located
{
return(0);
}

int
CSSRDiscovery::ReportSAM(int RepElLen,	// identified SSR contains repeat elements of this length
	   int NumTandemEls,			// each repeat element is repeated this many times
	   INT64 SSRStartOfs,			// repeat element starts at this offset within the targeted sequence
	    char *pszDescr,				// descriptor for the targeted sequence
		char *pszInFile,			// sequence parsed from this file
	   INT64 TargSeqLen,			// targeted sequence is this length
	   etSeqBase *pTargSeq)			// targeted sequence within which the SSR has been located
{
return(0);
}

// Reporting as CSV, BED or SAM
int
CSSRDiscovery::Report(int RepElLen,	// identified SSR contains repeat elements of this length
	   int NumTandemEls,			// each repeat element is repeated this many times
	   INT64 SSRStartOfs,			// repeat element starts at this offset within the targeted sequence
	    char *pszDescr,				// descriptor for the targeted sequence
		char *pszInFile,			// sequence parsed from this file
	   INT64 TargSeqLen,			// targeted sequence is this length
	   etSeqBase *pTargSeq)			// targeted sequence within which the SSR has been located
{
int Rslt;
switch(m_RptSSRsFormat) {
	case eRFCsv:			// CSV
		Rslt = ReportCSV(RepElLen,NumTandemEls,SSRStartOfs,pszDescr,pszInFile,TargSeqLen,pTargSeq);
		break;

	case eRFBed:			// BED
		Rslt = ReportBED(RepElLen,NumTandemEls,SSRStartOfs,pszDescr,pszInFile,TargSeqLen,pTargSeq);
		break;

	case eRFSam:			// SAM
		Rslt = ReportSAM(RepElLen,NumTandemEls,SSRStartOfs,pszDescr,pszInFile,TargSeqLen,pTargSeq);
		break;
		break;
	}
return(Rslt);
}


UINT16
GenSeqHash16(int SeqLen,	// hash this length sequence
			etSeqBase *pSeq) // sequence to hash
{
UINT32 Hash;
int Idx;

Hash = 19937;			// circular prime as hash seed

// hash over the element
for(Idx = 0; Idx < SeqLen; Idx++)
	{
	Hash = (Hash ^ (*pSeq & 0x07)) * 3119;	// a circular prime
	Hash ^= (Hash >> 13);
	Hash &= 0x0ffff;
	}
if(Hash == 0)			// 0 reserved as an error indicator so force hash to be 1
	Hash += 1;
return(Hash);
}


// if kmer is less than 10bp then count number of instances of KMer

int
CSSRDiscovery::CntKMer(int KMerLen,	 // count this KMer
					int Rpts,		 // tandem repeat counts
					etSeqBase *pSeq) // sequence 
{
int Idx;
UINT32 FreqIdx;
tsKMerDist *pKMerDist;

if(m_pKMerDist == NULL || KMerLen != m_KMerFreqLen)
	return(0);

FreqIdx = 0;
for(Idx = 0; Idx < KMerLen; Idx++)
	{
	FreqIdx <<= 2;
	FreqIdx |= 0x03 & *pSeq++;
	}

pKMerDist = (tsKMerDist *)((UINT8 *)m_pKMerDist + (m_SizeOfKMerDist * FreqIdx));

pKMerDist->Cnt += 1;
pKMerDist->TandemRpts[Rpts - m_MinTandemRpts] += 1;
return(pKMerDist->Cnt);
}


int
CSSRDiscovery::IdentifySSRs(char *pszDescr,	// descriptor for the targeted sequence
			 char *pszInFile,			// sequence parsed from this file
			 INT64 TargSeqLen,			// sequence length of targeted sequence within which to search for SSRs
			 etSeqBase *pTargSeq)		// identify SSRs in this targeted sequence
{
bool bSlough;
etSeqBase BaseA;
etSeqBase BaseB;
INT64 Ofs;
INT64 SSRStartOfs;
int RepElLen;
int RepElOfs;
int NumTandemEls;
int NumSSRs;
etSeqBase *pBase;
etSeqBase *pRepElBase;
etSeqBase *pRepMark;
int MarkLen;

// ensure no flag bits set in most significan nibble
pBase = pTargSeq;
for(Ofs = 0; Ofs < TargSeqLen; Ofs++,pBase++)
	*pBase &= 0x07;

NumSSRs = 0;
for(RepElLen = m_MinRepElLen; RepElLen <= m_MaxRepElLen; RepElLen++)
	{
	pBase = pTargSeq;
	SSRStartOfs = 0;
	NumTandemEls = 0;
	bSlough = false;
	for(Ofs = 0; Ofs < (TargSeqLen - RepElLen); Ofs++,pBase++)
		{
		pRepElBase = pBase;
		for(RepElOfs = 0; RepElOfs < RepElLen; RepElOfs++,pRepElBase++)
			{
			BaseA = *pRepElBase;
			BaseB = pRepElBase[RepElLen];
			if(BaseA > eBaseT || BaseB > eBaseT)	// not interested in potential SSR if it contains any non-cannonical base
				{
				bSlough = true;
				continue;
				}
			if(BaseA != BaseB)
				break;
			}
		if(RepElOfs == RepElLen)
			{
			if(NumTandemEls == 0)
				{
				SSRStartOfs = Ofs;
				NumTandemEls = 1;
				}
			NumTandemEls += 1;
			Ofs += RepElLen - 1;
			pBase = &pTargSeq[Ofs];
			}
		else
			{
			if(!bSlough && m_MinRepElLen > 1 && NumTandemEls)
				{
				int Zdx;
				int SubK;
				etSeqBase *pZBase;
				etSeqBase *pYBase;
				etSeqBase ZBase;
				pZBase = &pTargSeq[SSRStartOfs];
				ZBase = *pZBase++;
				for(Zdx = 1; Zdx < RepElLen; Zdx++,pZBase++)
					{
					if(*pZBase != ZBase)
						break;
					}
				if(Zdx == RepElLen)
					bSlough = true;
				if(!bSlough && RepElLen >= 4)
					{
					for(SubK = 2; SubK < RepElLen; SubK++)
						{
						if(RepElLen % SubK)
							continue;
						pZBase = &pTargSeq[SSRStartOfs];
						pYBase = pZBase + SubK;
						for(Zdx = 0;Zdx < SubK;Zdx++)
							{
							if(*pZBase++ != *pYBase++)
								break;
							}
						if(Zdx == SubK)
							{
							bSlough = true;
							break;
							}
						}
					}

				}
			if(!bSlough && NumTandemEls >= m_MinTandemRpts)
				{
				if(NumTandemEls > m_MaxTandemRpts)
					m_TotNumExcessiveTandemSSRs += 1;
				else
					{
					m_TotNumAcceptedSSRs += 1;
					Report(RepElLen,NumTandemEls,SSRStartOfs,pszDescr,pszInFile,TargSeqLen,pTargSeq);
					NumSSRs += 1;
					m_TotNumAcceptedKmerSSRs[RepElLen] += 1;
					CntKMer(RepElLen,NumTandemEls,&pTargSeq[SSRStartOfs]);

					// mark this sequence so it doesn't get accepted as being a longer K-mer SSR
					pRepMark = &pTargSeq[SSRStartOfs];
					MarkLen = RepElLen * NumTandemEls;
					do {
						*pRepMark++ |= 0x08;
						}
					while(MarkLen--);
					}
				}
			NumTandemEls = 0;
			bSlough = false;
			}
		if(!(Ofs % 100))
			ReportProgress();
		}
	}
return(NumSSRs);
}	

int
CSSRDiscovery::ReportKMers(char *pszKMerFreqFile)	// report SSR repeating element K-mer frequencies to this file
{
int Idx;
tsKMerDist *pFreqDist;
UINT32 SeqIdx;
char Base;
int Idy;
int Rpts;
char szKmer[100];
int BuffIdx;
char szBuff[0x03fff];

if(m_hOutKMerFreqFile != -1 && m_KMerFreqLen != 0)
	{
	BuffIdx = 0;
	pFreqDist = m_pKMerDist;
	for(Idx = 0; Idx < m_NumKMers; Idx+=1)
		{
		SeqIdx = Idx;
		for(Idy = 0; Idy < m_KMerFreqLen; Idy++)
			{
			switch(SeqIdx & 0x03) {
				case 0:
					Base = 'A';
					break;
				case 1:
					Base = 'C';
					break;
				case 2:
					Base = 'G';
					break;
				case 3:
					Base = 'T';
					break;
				}
			szKmer[m_KMerFreqLen - (1+Idy)] = Base;
			SeqIdx >>= 2;
			}
		szKmer[m_KMerFreqLen] = '\0';
		BuffIdx += sprintf(&szBuff[BuffIdx],"\"%s\",%d",szKmer,pFreqDist->Cnt);
		for(Rpts = cMinTandemRpts; Rpts <= m_MaxTandemRpts; Rpts++)
			if(Rpts < m_MinTandemRpts)
				BuffIdx += sprintf(&szBuff[BuffIdx],",0");
			else
				BuffIdx += sprintf(&szBuff[BuffIdx],",%d",pFreqDist->TandemRpts[Rpts - m_MinTandemRpts]);
		BuffIdx += sprintf(&szBuff[BuffIdx],"\n");
		if((BuffIdx + 1000) > (int)sizeof(szBuff))
			{
			CUtility::SafeWrite(m_hOutKMerFreqFile,szBuff,BuffIdx);
			BuffIdx = 0;
			}
		pFreqDist = (tsKMerDist *)((UINT8 *)pFreqDist + m_SizeOfKMerDist);
		}
	if(BuffIdx > 0)
		CUtility::SafeWrite(m_hOutKMerFreqFile,szBuff,BuffIdx);
#ifdef _WIN32
	_commit(m_hOutKMerFreqFile);
#else
	fsync(m_hOutKMerFreqFile);
#endif
	close(m_hOutKMerFreqFile);
	m_hOutKMerFreqFile = -1;
	}
return(eBSFSuccess);
}


int
CSSRDiscovery::ReportProgress(bool bForce)	// let user know that there is processing activity, normally progress reported every 30 sec unless bForce set true
{
static long PrevSecs = 0;
long CurSecs = m_CurTime.ReadUSecs();
if(bForce || CurSecs < PrevSecs || ((CurSecs - PrevSecs) > 60))
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Progress: Identified %u SSRs",m_TotNumAcceptedSSRs);
	PrevSecs = CurSecs;
	}
return(0);
}

int
CSSRDiscovery::Process(int PMode,			// procesisng mode - currently unused..
		teRptSSRsFromat RptSSRsFormat,		// report SSRs in this file format
		int MinRepElLen,					// identify repeating elements of this minimum length
		int MaxRepElLen,					// ranging upto this maximum length
		int MinTandemRpts,					// minimum number of tandem repeats
		int MaxTandemRpts,					// maximum number of repeats
		int SSRFlankLen,					// SSR flanking sequences lengths to report
		int NumInFileSpecs,					// number of input, could be wildcarded, file specs
		char *pszInFiles[],					// files to be processed
		char *pszKMerFreqFile,				// optional, output element KMer freq to this file
		char *pszOutFile)					// SSRs to this file
{
int Rslt;
int Idx;
UINT32 TotNumAcceptedSSRs;
UINT32 TotNumExcessiveTandemSSRs;

char *pszInFile;
CSimpleGlob glob(SG_GLOB_FULLSORT);

m_RptSSRsFormat = RptSSRsFormat;
m_MinRepElLen = MinRepElLen;
m_MaxRepElLen = MaxRepElLen;
m_MinTandemRpts = MinTandemRpts;
m_MaxTandemRpts = MaxTandemRpts;
m_SSRFlankLen = SSRFlankLen;

m_TotNumExcessiveTandemSSRs = 0;
m_TotNumAcceptedSSRs = 0;

TotNumAcceptedSSRs = 0;
TotNumExcessiveTandemSSRs = 0;

memset(m_TotNumAcceptedKmerSSRs,0,sizeof(m_TotNumAcceptedKmerSSRs));

// if single KMer element length being processed and that KMer is <= 10 then allocate to hold instance counts
if(pszKMerFreqFile != NULL && pszKMerFreqFile[0] != '\0' &&
   MinRepElLen == MaxRepElLen && MinRepElLen <= 10)
	{
	size_t memreq;
	int Power2 = MinRepElLen;
	m_SizeOfKMerDist = (int)(sizeof(tsKMerDist) + (sizeof(UINT32) * (MaxTandemRpts - MinTandemRpts)));
	m_NumKMers = 1;
	while(Power2--)
		m_NumKMers <<= 2;
	memreq = (size_t)m_NumKMers * (size_t)m_SizeOfKMerDist;
#ifdef _WIN32
	m_pKMerDist = (tsKMerDist *) malloc(memreq);	// initial and only allocation
	if(m_pKMerDist == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: Memory allocation of %lld bytes - %s",(INT64)memreq,strerror(errno));
		Reset();
		return(eBSFerrMem);
		}
#else
	// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
	m_pKMerDist = (tsKMerDist *)mmap(NULL,memreq, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
	if(m_pKMerDist == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: Memory allocation of %lld bytes through mmap()  failed - %s",(INT64)memreq,strerror(errno));
		m_pKMerDist = NULL;
		Reset();
		return(eBSFerrMem);
		}
#endif
	m_AllocdKMerFreqMem = memreq;
	memset(m_pKMerDist,0,memreq);
	m_KMerFreqLen = MaxRepElLen;
#ifdef _WIN32
	if((m_hOutKMerFreqFile = open(pszKMerFreqFile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
	if((m_hOutKMerFreqFile = open(pszKMerFreqFile,O_RDWR | O_CREAT |O_TRUNC, S_IREAD | S_IWRITE))==-1)
#endif
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create or truncate %s - %s",pszKMerFreqFile,strerror(errno));
		Reset();
		return(eBSFerrCreateFile);
		}
	char szBuff[2000];
	int BuffIdx;
	BuffIdx = sprintf(szBuff,"\"K-Mer\",\"Cnt\"");
	for(Idx = cMinTandemRpts; Idx <= m_MaxTandemRpts; Idx++)
		BuffIdx += sprintf(&szBuff[BuffIdx],",\"Rpts:%d\"",Idx);
	BuffIdx += sprintf(&szBuff[BuffIdx],"\n");
	CUtility::SafeWrite(m_hOutKMerFreqFile,szBuff,BuffIdx);
	}
else
	{
	m_AllocdKMerFreqMem = 0;
	m_SizeOfKMerDist = 0;
	m_KMerFreqLen = 0;
	m_NumKMers = 0;
	}

if((m_pszRptSSRsBuff = new char [cMaxAllocRptSSRs]) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate %d bytes buffering for SSRs output",cMaxAllocRptSSRs);
	Reset();
	return(eBSFerrMem);
	}
m_IdxRptSSRs = 0;

#ifdef _WIN32
if((m_hOutFile = open(pszOutFile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
if((m_hOutFile = open(pszOutFile,O_RDWR | O_CREAT |O_TRUNC, S_IREAD | S_IWRITE))==-1)
#endif
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create or truncate %s - %s",pszOutFile,strerror(errno));
	Reset();
	return(eBSFerrCreateFile);
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Output results file created/truncated: '%s'",pszOutFile);

for(Idx = 0; Idx < NumInFileSpecs; Idx++)
	{
	glob.Init();
	if(glob.Add(pszInFiles[Idx]) < SG_SUCCESS)
		{
		gDiagnostics.DiagOut(eDLWarn,gszProcName,"Unable to glob '%s",pszInFiles[Idx]);
		return(eBSFerrOpnFile);	// treat as though unable to open file
		}

	if(glob.FileCount() <= 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to locate any input file matching '%s",pszInFiles[Idx]);
		continue;
		}
	Rslt = eBSFSuccess;
	for (int FileID = 0; Rslt >= eBSFSuccess &&  FileID < glob.FileCount(); ++FileID)
		{
		pszInFile = glob.File(FileID);
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing file '%s'",pszInFile);
		if((Rslt = ProcessFastaFile(pszInFile)) < eBSFSuccess)
			{
			gDiagnostics.DiagOut(eDLWarn,gszProcName,"Failed procesisng file '%s",pszInFile);
			return(Rslt);
			}
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Accepted %u SSRs, rejected %u excessively long",
							m_TotNumAcceptedSSRs - TotNumAcceptedSSRs,
							m_TotNumExcessiveTandemSSRs - TotNumExcessiveTandemSSRs);
		TotNumAcceptedSSRs = m_TotNumAcceptedSSRs;
		TotNumExcessiveTandemSSRs = m_TotNumExcessiveTandemSSRs;
		}
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Total SSRs accepted %u, rejected %u tandem repeats as excessively long",m_TotNumAcceptedSSRs,m_TotNumExcessiveTandemSSRs);
for(Idx=MinRepElLen; Idx <= MaxRepElLen; Idx++)
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"K-mer %d, Number %u",Idx,m_TotNumAcceptedKmerSSRs[Idx]);



if(m_hOutFile != -1)
	{
	if(m_IdxRptSSRs)
		{
		CUtility::SafeWrite(m_hOutFile,m_pszRptSSRsBuff,m_IdxRptSSRs);
		m_IdxRptSSRs = 0;
		}
#ifdef _WIN32
	_commit(m_hOutFile);
#else
	fsync(m_hOutFile);
#endif
	close(m_hOutFile);
	m_hOutFile = -1;
	}


Rslt = ReportKMers(pszKMerFreqFile);
Reset();
return(Rslt);
}