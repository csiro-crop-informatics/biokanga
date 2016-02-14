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
// PEScaffold.cpp : Defines the entry point for the console application.
// Purpose -
// Loads SAM alignments, whereby PE1 and PE2 were independently aligned to a targeted de Novo assembly, and attempts to
// report scaffolding potential for the assembly contigs in a CSV file
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
#include "PEScaffold.h"

int Process(int PMode,					// processing mode
		char *pszSeqIDTerm,				// pair sequence identifiers until this terminating character(s) - defaults to none terminating
		char *pszInPE1File,				// input PE1 file
		char *pszInPE2File,				// input PE2 file
		char *pszOutFile);				// output corelations file

#ifdef _WIN32
int pescaffold(int argc, char* argv[])
{
// determine my process name
_splitpath(argv[0],NULL,NULL,gszProcName,NULL);
#else
int pescaffold(int argc, char** argv)
{
// determine my process name
CUtility::splitpath((char *)argv[0],NULL,gszProcName);
#endif

int iScreenLogLevel;		// level of screen diagnostics
int iFileLogLevel;			// level of file diagnostics
char szLogFile[_MAX_PATH];	// write diagnostics to this file
int Rslt = 0;   			// function result code >= 0 represents success, < 0 on failure

int PMode;					// processing mode

char szSeqIDTerm[50];				// pair sequence identifiers until this terminating character(s) - defaults to none terminating
char szInPE1File[_MAX_PATH];	// parse PE1 alignments from this file
char szInPE2File[_MAX_PATH];	// parse PE1 alignments from this file

char szOutFile[_MAX_PATH];		// write corelations to this file

char szSQLiteDatabase[_MAX_PATH];	// results summaries to this SQLite file
char szExperimentName[cMaxDatasetSpeciesChrom+1];			// experiment name
char szExperimentDescr[1000];		// describes experiment

// command line args
struct arg_lit  *help    = arg_lit0("h","help",                 "Print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"Print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to screen and logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"Diagnostics log file");

struct arg_int *mode = arg_int0("m","mode","<int>",				"Processing mode: 0 default");


struct arg_str *seqidterm = arg_str0("p","seqidterm","<str>",	"Truncate sequence identifiers prefixes at rightmost instance of any one of these chars (default is for no truncation)");
struct arg_file *inpe1file = arg_file1("i","in","<file>",		"Input SAM file containing PE1 alignments");
struct arg_file *inpe2file = arg_file1("I","in","<file>",		"Input SAM file containing PE2 alignments");

struct arg_file *outfile = arg_file1("o","out","<file>",		"Output corelations to this file");

struct arg_file *summrslts = arg_file0("q","sumrslts","<file>",		"Output results summary to this SQLite3 database file");
struct arg_str *experimentname = arg_str0("w","experimentname","<str>",		"experiment name SQLite3 database file");
struct arg_str *experimentdescr = arg_str0("W","experimentdescr","<str>",	"experiment description SQLite3 database file");

struct arg_end *end = arg_end(200);

void *argtable[] = {help,version,FileLogLevel,LogFile,
					summrslts,experimentname,experimentdescr,
	                mode,seqidterm,inpe1file,inpe2file, outfile,
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

	PMode = mode->count ? mode->ival[0] : 0;
	if(PMode < 0 || PMode > 1)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Processing mode '-m%d' must be in range 0..0",PMode);
		return(1);
		}

	szSeqIDTerm[0] = '\0';
	if(seqidterm->count > 0)
		{
		strncpy(szSeqIDTerm,seqidterm->sval[0],sizeof(szSeqIDTerm));	
		szSeqIDTerm[sizeof(szSeqIDTerm)-1] = '\0';
		CUtility::TrimQuotedWhitespcExtd(szSeqIDTerm);
		}
	
	strcpy(szInPE1File,inpe1file->filename[0]);
	CUtility::TrimQuotedWhitespcExtd(szInPE1File);

	strcpy(szInPE2File,inpe2file->filename[0]);
	CUtility::TrimQuotedWhitespcExtd(szInPE2File);

	strcpy(szOutFile,outfile->filename[0]);
	CUtility::TrimQuotedWhitespcExtd(szOutFile);

// show user current resource limits
#ifndef _WIN32
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Resources: %s",CUtility::ReportResourceLimits());
#endif

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing parameters:");
	const char *pszDescr;
	switch(PMode) {
		case 0:
			pszDescr = "Corelate PE1 and PE2 alignments";
			break;
		}

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Processing mode is : '%s'",pszDescr);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Experiment name : '%s'",szExperimentName);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Experiment description : '%s'",szExperimentDescr);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"SQLite database file: '%s'",szOutFile);
	if(szSeqIDTerm[0] != '\0')
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Truncate sequence identifiers at last : '%s'",szSeqIDTerm);
	else
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Truncate sequence identifiers at last : 'No truncation'");
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"PE1 SAM file: '%s'",szInPE1File);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"PE2 SAM file: '%s'",szInPE2File);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Corelations to file: '%s'",szOutFile);

	if(gExperimentID > 0)
		{
		int ParamID;
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szLogFile),"log",szLogFile);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,sizeof(PMode),"pmode",&PMode);
		if(szSeqIDTerm[0] != '\0')
			ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szSeqIDTerm),"seqidtrim",szSeqIDTerm);
		else
			ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,3,"seqidtrim",(void *)"N/A");
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szInPE1File),"inpe1file",szInPE1File);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szInPE2File),"inpe2file",szInPE2File);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szOutFile),"outfile",szOutFile);

		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szSQLiteDatabase),"sumrslts",szSQLiteDatabase);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szExperimentName),"experimentname",szExperimentName);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szExperimentDescr),"experimentdescr",szExperimentDescr);
		}

#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	gStopWatch.Start();

	Rslt = Process(PMode,szSeqIDTerm,szInPE1File,szInPE2File,szOutFile);
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


int Process(int PMode,					// processing mode
		char *pszSeqIDTerm,				// pair sequence identifiers until this terminating character(s) - defaults to none terminating
		char *pszInPE1File,				// input PE1 file
		char *pszInPE2File,				// input PE2 file
		char *pszOutFile)				// output corelations file
{
int Rslt;
CPEScaffold *pPEScaffold;
if((pPEScaffold = new CPEScaffold)==NULL)
	return(-1);
Rslt = pPEScaffold->Process(PMode,pszSeqIDTerm,pszInPE1File,pszInPE2File,pszOutFile);
delete pPEScaffold;
return(Rslt);
}


CPEScaffold::CPEScaffold(void)
{
Init();
}


CPEScaffold::~CPEScaffold(void)
{
Reset();
}

void 
CPEScaffold::Reset(void)
{
if(m_hOutFile != -1)
	{
	close(m_hOutFile);
	m_hOutFile = -1;
	}
if(m_pPEIdents != NULL)
	{
#ifdef _WIN32
	free(m_pPEIdents);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pPEIdents != MAP_FAILED)
		munmap(m_pPEIdents,m_AllocdPEIdentsMem);
#endif
	m_pPEIdents = NULL;
	}

if(m_pHashPEIdents)
	{
	delete m_pHashPEIdents;
	m_pHashPEIdents = NULL;
	}

if(m_pScaffoldContigs != NULL)
	{
#ifdef _WIN32
	free(m_pScaffoldContigs);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pScaffoldContigs != MAP_FAILED)
		munmap(m_pScaffoldContigs,m_AllocdScaffoldContigsMem);
#endif
	m_pScaffoldContigs = NULL;
	}

if(m_pHashContigs)
	{
	delete m_pHashContigs;
	m_pHashContigs = NULL;
	}

if(m_pScaffolds != NULL)
	{
#ifdef _WIN32
	free(m_pScaffolds);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pScaffolds != MAP_FAILED)
		munmap(m_pScaffolds,m_AllocdScaffoldsMem);
#endif
	m_pScaffolds = NULL;
	}


if(m_ppPE2Scaffolds != NULL)
	{
#ifdef _WIN32
	free(m_ppPE2Scaffolds);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_ppPE2Scaffolds != MAP_FAILED)
		munmap(m_ppPE2Scaffolds,m_AllocdPE2ScaffoldsMem);
#endif
	m_ppPE2Scaffolds = NULL;
	}

m_szSeqIDTermChrs[0] = '\0';
m_AllocdNumPEIdents = 0;
m_NumPEIdents = 0;
m_AllocdPEIdentsMem = 0;

m_AllocdNumScaffoldContigs = 0;
m_AllocdScaffoldContigsMem = 0;
m_NumScaffoldContigs = 0;

m_AllocdScaffoldsMem = 0;
m_AllocdNumScaffolds = 0;
m_NumScaffolds = 0;

m_AllocdPE2ScaffoldsMem = 0;
}

void 
CPEScaffold::Init(void)
{
m_pPEIdents = NULL;
m_pHashPEIdents = NULL;
m_pScaffoldContigs = NULL;
m_pHashContigs = NULL;
m_pScaffolds = NULL;
m_ppPE2Scaffolds = NULL;
m_hOutFile = -1;
Reset();
}

char *
CPEScaffold::TrimWhitespace(char *pTxt)
{
char *pStart;
char Chr;
	// strip leading whitespace
while(Chr = *pTxt++)
	if(!isspace(Chr))
			break;
if(Chr == '\0')					// empty line?
	return(pTxt-1);
pStart = pTxt-1;
while(Chr = *pTxt)			// fast forward to line terminator
	pTxt++;
pTxt-=1;
while(Chr = *pTxt--)
	if(!isspace(Chr))
		break;
pTxt[2] = '\0';
return(pStart);
}

// AddContigName
// If chrom already known then return existing chrom identifier otherwise add to m_pScaffoldContigs
int
CPEScaffold::AddContigName(char *pszContigName)
{
static int PrevContigID = 0;
int Hash;
int HashIdx;
tsPEScaffoldContig *pScaffoldContig;

// check to see if chromosome name already known
if(PrevContigID != 0)
	{
	pScaffoldContig = &m_pScaffoldContigs[PrevContigID-1];
	if(!stricmp(pszContigName,pScaffoldContig->szContig))
		return(pScaffoldContig->ContigID);
	}

Hash = CUtility::GenHash24(pszContigName);
if((HashIdx = m_pHashContigs[Hash]) != 0)
	{
	do {
		pScaffoldContig = &m_pScaffoldContigs[HashIdx-1];
		if(!stricmp(pszContigName,pScaffoldContig->szContig))
			{
			PrevContigID = pScaffoldContig->ContigID;
			return(pScaffoldContig->ContigID);
			}
		HashIdx = pScaffoldContig->HashNext;
		}
	while(HashIdx > 0);
	}

// its a new chrom not previously seen
// realloc as may be required to hold this new chrom
if(m_NumScaffoldContigs == m_AllocdNumScaffoldContigs)
	{
	size_t memreq = m_AllocdScaffoldContigsMem + (cAllocContigNames * sizeof(tsPEScaffoldContig));
#ifdef _WIN32
	pScaffoldContig = (tsPEScaffoldContig *) realloc(m_pScaffoldContigs,memreq);
#else
	pScaffoldContig = (tsPEScaffoldContig *)mremap(m_pScaffoldContigs,m_AllocdScaffoldContigsMem,memreq,MREMAP_MAYMOVE);
	if(pScaffoldContig == MAP_FAILED)
		pScaffoldContig = NULL;
#endif
	if(pScaffoldContig == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddContigName: Memory re-allocation to %d bytes - %s",memreq,strerror(errno));
		return(eBSFerrMem);
		}
	m_pScaffoldContigs = pScaffoldContig;
	m_AllocdScaffoldContigsMem = memreq;
	m_AllocdNumScaffoldContigs += cAllocContigNames;
	}

pScaffoldContig = &m_pScaffoldContigs[m_NumScaffoldContigs++];
pScaffoldContig->ContigID = m_NumScaffoldContigs;
pScaffoldContig->HashNext = m_pHashContigs[Hash];
pScaffoldContig->NumClustered = 0;
pScaffoldContig->ClusterID = 0;
m_pHashContigs[Hash] = m_NumScaffoldContigs;
strncpy(pScaffoldContig->szContig,pszContigName,cMaxContigNameLen);
pScaffoldContig->szContig[cMaxContigNameLen] = '\0';
PrevContigID = m_NumScaffoldContigs;
return(m_NumScaffoldContigs);
}

char *
CPEScaffold::GetContigName(int ContigID)
{
tsPEScaffoldContig *pScaffoldContig;
pScaffoldContig = &m_pScaffoldContigs[ContigID-1];
return(pScaffoldContig->szContig);
}

// AddPEIdent
// If PE name already known then return existing identifier otherwise add to m_pScaffoldContigs
int
CPEScaffold::AddPEIdent(char *pszIdentName)
{
static INT32 PrevPEIdentID = 0;
tsPEIdent *pPEIdent;
int Hash;
int HashIdx;
int IdentLen;
char Chr;
char TermChr;
char *pIdentChr;
char *pTermChr;

// identifiers are only significant up until user specified set of terminating chars
if(m_szSeqIDTermChrs[0] != '\0')
	{
	if((IdentLen = (int)strlen(pszIdentName)) >= 4)		// don't bother triming identifiers which are too short..
		{
		pIdentChr = &pszIdentName[IdentLen-1];			// search for terminating char starts from last char ...
		IdentLen -= 3;									// ensures identifiers after trimming are at least 3 chrs long
		while(IdentLen--)
			{
			Chr = *pIdentChr--;
			pTermChr = m_szSeqIDTermChrs;
			while((TermChr = *pTermChr++) != '\0')
				{
				if(Chr == TermChr)
					{
					pIdentChr[1] = '\0';
					break;
					}
				}
			if(pIdentChr[1] == '\0')
				break;
			}
		}
	}

// check to see if PEIdent name already known
if(PrevPEIdentID != 0)
	{
	pPEIdent = &m_pPEIdents[PrevPEIdentID-1];
	if(!stricmp(pszIdentName,pPEIdent->szIdent))
		return(pPEIdent->IdentID);
	}

Hash = CUtility::GenHash24(pszIdentName);
if((HashIdx = m_pHashPEIdents[Hash]) != 0)
	{
	do {
		pPEIdent = &m_pPEIdents[HashIdx-1];
		if(!stricmp(pszIdentName,pPEIdent->szIdent))
			{
			PrevPEIdentID = pPEIdent->IdentID;
			return( pPEIdent->IdentID);
			}
		HashIdx = pPEIdent->HashNext;
		}
	while(HashIdx > 0);
	}

// its a new PE identifier not previously seen
// realloc as may be required to hold this new chrom
if(m_NumPEIdents == m_AllocdNumPEIdents)
	{
	size_t memreq = m_AllocdPEIdentsMem + (cAllocPENames * sizeof(tsPEIdent));
#ifdef _WIN32
	pPEIdent = (tsPEIdent *) realloc(m_pPEIdents,memreq);
#else
	pPEIdent = (tsPEIdent *)mremap(m_pPEIdents,m_AllocdPEIdentsMem,memreq,MREMAP_MAYMOVE);
	if(pPEIdent == MAP_FAILED)
		pPEIdent = NULL;
#endif
	if(pPEIdent == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddPEIdent: Memory re-allocation to %d bytes - %s",memreq,strerror(errno));
		return(eBSFerrMem);
		}
	m_pPEIdents = pPEIdent;
	m_AllocdPEIdentsMem = memreq;
	m_AllocdNumPEIdents += cAllocPENames;
	}

pPEIdent = &m_pPEIdents[m_NumPEIdents++];
pPEIdent->IdentID = m_NumPEIdents;
pPEIdent->PEScafoldID = 0;
pPEIdent->HashNext = m_pHashPEIdents[Hash];
m_pHashPEIdents[Hash] = m_NumPEIdents;
strncpy(pPEIdent->szIdent,pszIdentName,cMaxContigNameLen);
pPEIdent->szIdent[cMaxContigNameLen] = '\0';
PrevPEIdentID = m_NumPEIdents;
return(m_NumPEIdents);
}

char *
CPEScaffold::GetSeqName(int SeqID)
{
tsPEIdent *pIdent;
pIdent = &m_pPEIdents[SeqID-1];
return(pIdent->szIdent);
}

int
CPEScaffold::AddScaffold(bool bPE2,						// if false then PE1, if true then PE2
			char *pszPEIdent,				// paired end indentifier used to corelate paired ends
			char *pszContig,					// PE aligns onto this chromosome
			char Strand)					// '+' or '-'
{
static int PrevScafoldID = 0;
tsPEScaffold *pPEScaffold;
tsPEIdent *pPEIdent;
int ContigID;
int PEIdentID;

ContigID = AddContigName(pszContig);
PEIdentID = AddPEIdent(pszPEIdent);

// if processing PE2 then check to see if already have scaffold with same PEIdent
if(bPE2)
	{
	if(PrevScafoldID != 0)
		{
		pPEScaffold = &m_pScaffolds[PrevScafoldID-1];
		if(pPEScaffold->PE12SeqID == PEIdentID)
			{
			pPEScaffold->PE2ContigID = ContigID;
			pPEScaffold->PE2Sense = Strand == '+' ? 1 : 0;
			return(PrevScafoldID);
			}
		}

	pPEIdent = &m_pPEIdents[PEIdentID-1];
	if((PrevScafoldID = pPEIdent->PEScafoldID) != 0)
		{
		pPEScaffold = &m_pScaffolds[PrevScafoldID-1];
		pPEScaffold->PE2ContigID = ContigID;
		pPEScaffold->PE2Sense = Strand == '+' ? 1 : 0;
		return(PrevScafoldID);
		}
	}
else	// PE1 processing, normally expect that the PE1 identifiers are unique, here we are just confirming that they are unique
	{   // when satisfied PE1 identifiers are always unique then could skip this check...
	pPEIdent = &m_pPEIdents[PEIdentID-1];
	if((PrevScafoldID = pPEIdent->PEScafoldID) != 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddScaffold: duplicate PE1 identifer - %s onto %s",pszPEIdent,pszContig);
		return(0);
		}
	}

// new scaffold required
// realloc as may be required to hold this new scaffold
if(m_NumScaffolds == m_AllocdNumScaffolds)
	{
	size_t memreq = m_AllocdScaffoldsMem + (cAllocScafolds * sizeof(tsPEScaffold));
#ifdef _WIN32
	pPEScaffold = (tsPEScaffold *) realloc(m_pScaffolds,memreq);
#else
	pPEScaffold = (tsPEScaffold *)mremap(m_pScaffolds,m_AllocdScaffoldsMem,memreq,MREMAP_MAYMOVE);
	if(pPEScaffold == MAP_FAILED)
		pPEScaffold = NULL;
#endif
	if(pPEScaffold == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddScaffold: Memory re-allocation to %d bytes - %s",memreq,strerror(errno));
		return(eBSFerrMem);
		}
	m_pScaffolds = pPEScaffold;
	m_AllocdScaffoldsMem = memreq;
	m_AllocdNumScaffolds += cAllocScafolds;
	}
pPEScaffold = &m_pScaffolds[m_NumScaffolds++];
memset(pPEScaffold,0,sizeof(tsPEScaffold));
pPEScaffold->PEScafoldID = m_NumScaffolds;
pPEIdent->PEScafoldID = m_NumScaffolds;
pPEScaffold->PE12SeqID = PEIdentID;
if(!bPE2)
	{
	pPEScaffold->PE1ContigID = ContigID;
	pPEScaffold->PE1Sense = Strand == '+' ? 1 : 0;
	}
else
	{
	pPEScaffold->PE2ContigID = ContigID;
	pPEScaffold->PE2Sense = Strand == '+' ? 1 : 0;
	}
return(m_NumScaffolds);
}

int
CPEScaffold::LoadSAM(bool bPE2,			// false if loading PE1, true if loading PE2
		char *pszSAMFile)	// load alignments from this SAM file
{
int Rslt;
etClassifyFileType FileType;
int NumParsedElLines;
int NumAcceptedEls;
char szLine[16000];				// buffer input lines
char *pTxt;
char szDescriptor[250];			// parsed out descriptor
int Flags;						// parsed out flags
char szContig[128];				// parsed out chrom
int StartLoci;					// start loci
int NumUnmappedEls;
int ScaffoldID;
CSAMfile BAMfile;
int LineLen;

	// open SAM for reading
if(pszSAMFile == NULL || *pszSAMFile == '\0')
	return(eBSFerrParams);

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Loading alignments for %s from: '%s'",bPE2 ? "PE2" : "PE1", pszSAMFile);

FileType = CUtility::ClassifyFileType(pszSAMFile);
if(FileType != eCFTSAM)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to classify file as SAM formated: '%s'",pszSAMFile);
	return(eBSFerrOpnFile);
	}

if((Rslt = (teBSFrsltCodes)BAMfile.Open(pszSAMFile)) != eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"ParseSAMFileElements: Unable to open SAM format file %s",pszSAMFile);
	return((teBSFrsltCodes)Rslt);
	}

NumParsedElLines = 0;
NumAcceptedEls = 0;
NumUnmappedEls = 0;

while((LineLen = BAMfile.GetNxtSAMline(szLine)) > 0)
	{
	NumParsedElLines += 1;
	if(!(NumParsedElLines % 1000000) || NumParsedElLines == 1)
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Parsed %d SAM lines",NumParsedElLines);

	szLine[sizeof(szLine)-1] = '\0';
	pTxt = TrimWhitespace(szLine);
	if(*pTxt=='\0' || *pTxt=='@')	// simply slough lines which were just whitespace or start with '@'
		continue;
	
	// expecting to parse as "%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t", szDescriptor, Flags, m_szSAMTargContigName, StartLoci+1,MAPQ,szCigar,pszRNext,PNext,TLen);
	// interest is in the descriptor,chromname,flags
	sscanf(szLine,"%s\t%d\t%s\t%d\t",szDescriptor, &Flags, szContig, &StartLoci);
		// check if element has been mapped, if not then slough ...
	if(StartLoci == 0 || Flags & 0x04)	// will be set if unmapped
		{
		NumUnmappedEls += 1;
	    continue;
		}

	if((ScaffoldID = AddScaffold(bPE2,szDescriptor,szContig,Flags & 0x10 ? '+' : '-')) < 1)
		{
		BAMfile.Close();
		return(ScaffoldID);
		}
	NumAcceptedEls += 1;
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Loading alignments (%d) for %s from: '%s' completed",NumAcceptedEls,bPE2 ? "PE2" : "PE1", pszSAMFile);
BAMfile.Close();
return(NumAcceptedEls);
}

// locate scaffold having matching PE1ContigID and PE2ContigID
// assumes scaffolds have been sorted PE1ContigID.PE2ContigID ascending
tsPEScaffold *
CPEScaffold::LocateMateScaffold(int PE1ContigID,int PE2ContigID)
{
int Mid;
int Hi;
int Lo;
tsPEScaffold *pEl1;
tsPEScaffold *pEl2;
Lo = 0;
Hi = m_NumScaffolds-1;
do {
	Mid = (Hi + Lo) / 2;
	pEl1 = &m_pScaffolds[Mid];
	if(pEl1->PE1ContigID == PE1ContigID)
		{
		if(pEl1->PE2ContigID == PE2ContigID)
			{
			// drill down to locate 1st instance - linear drill down as not expecting too many instances
			pEl2 = pEl1--;
			while(Mid-- && (pEl1->PE1ContigID == PE1ContigID && pEl1->PE2ContigID == PE2ContigID))
					pEl2 = pEl1--;
			return(pEl2);
			}
		if(pEl1->PE2ContigID > PE2ContigID)
			Hi = Mid - 1;
		else
			Lo = Mid + 1;
		}
	else
		{
		if(pEl1->PE1ContigID > PE1ContigID)
			Hi = Mid - 1;
		else
			Lo = Mid + 1;
		}
	}
while(Hi >= Lo);
return(NULL);		// unable to match any
}

// locate scaffold having matching ContigID as PE1
// assumes scaffolds have been sorted PE1ContigID.PE2ContigID ascending
int
CPEScaffold::LocatePE1Scaffold(int ContigID) // locate scaffold having matching ContigID for PE1
{
int Mark;
int Mid;
int Hi;
int Lo;
tsPEScaffold *pEl1;
if(m_ppPE2Scaffolds == NULL || m_NumScaffolds == 0)		// better safe than sorry...
	return(0);
Lo = 0;
Hi = m_NumScaffolds-1;
do {
	Mid = (Hi + Lo) / 2;
	pEl1 = &m_pScaffolds[Mid];
	if(pEl1->PE1ContigID == ContigID)
		{
		// drill down to locate 1st instance - linear drill down as not expecting too many instances
		Mark = Mid;
		while(Mid--)
			{
			pEl1 -= 1;
			if(pEl1->PE1ContigID != ContigID)
				break;
			Mark = Mid;
			}
		return(Mark+1);
		}
	else
		{
		if(pEl1->PE1ContigID > ContigID)
			Hi = Mid - 1;
		else
			Lo = Mid + 1;
		}
	}
while(Hi >= Lo);
return(0);		// unable to match any
}

// locate scaffold having matching ContigID as PE2
// assumes scaffolds have been sorted PE1ContigID.PE2ContigID ascending
int											// returns 0 if unable to locate
CPEScaffold::LocatePE2Scaffold(int ContigID) // locates scaffold having matching ContigID as PE2
{
int Mark;
int Mid;
int Hi;
int Lo;
tsPEScaffold *pEl1;
if(m_ppPE2Scaffolds == NULL || m_NumScaffolds == 0)		// better safe than sorry...
	return(0);
Lo = 0;
Hi = m_NumScaffolds-1;
do {
	Mid = (Hi + Lo) / 2;
	pEl1 = m_ppPE2Scaffolds[Mid];
	if(pEl1->PE2ContigID == ContigID)
		{
		// drill down to locate 1st instance - linear drill down as not expecting too many instances
		Mark = Mid;
		while(Mid--)
			{
			pEl1 = m_ppPE2Scaffolds[Mid];
			if(pEl1->PE2ContigID != ContigID)
				break;
			Mark = Mid;
			}
		return(Mark+1);
		}
	else
		{
		if(pEl1->PE2ContigID > ContigID)
			Hi = Mid - 1;
		else
			Lo = Mid + 1;
		}
	}
while(Hi >= Lo);
return(0);		// unable to match any
}


int						// returned next CurIdx to use on a subsequent call to IterPE1s (0 if all scaffords have been iterated)
CPEScaffold::IterPE1s(int CurIdx,	// iterates, in PE1ContigID.PE2ContigID ascending order, all tsPEScaffolds, to start iterations
		int ContigID,				// iterate for this PE1 chrom/contig
		tsPEScaffold **ppPEScaffold)	// returned scafford or NULL if all scaffords have been iterated
{
tsPEScaffold *pScaffold;

if(ppPEScaffold != NULL)
	*ppPEScaffold = NULL;
if(CurIdx < 0 || CurIdx >= m_NumScaffolds)
	return(0);
if(CurIdx == 0)
	{
	CurIdx = LocatePE1Scaffold(ContigID);
	if(CurIdx == 0)
		return(0);
	CurIdx -= 1;
	pScaffold = &m_pScaffolds[CurIdx];
	}
else
	{
	pScaffold = &m_pScaffolds[CurIdx];
	if(pScaffold->PE1ContigID != ContigID)
		return(0);
	}

if(ppPEScaffold != NULL)
	*ppPEScaffold = pScaffold;
return(CurIdx + 1);
}

int						// returned next CurIdx to use on a subsequent call to IterPE2s (0 if all scaffords have been iterated)
CPEScaffold::IterPE2s(int CurIdx,	// iterates, in PE2ContigID.PE1ContigID ascending order, all tsPEScaffolds, 0 to start iterations
		int ContigID,				// iterate for this PE2 chrom/contig
		tsPEScaffold **ppPEScaffold)	// returned scafford or NULL if all scaffords have been iterated
{
tsPEScaffold *pScaffold;

if(ppPEScaffold != NULL)
	*ppPEScaffold = NULL;
if(CurIdx < 0 || CurIdx >= m_NumScaffolds)
	return(0);

if(CurIdx == 0)
	{
	CurIdx = LocatePE2Scaffold(ContigID);
	if(CurIdx == 0)
		return(0);
	CurIdx -= 1;
	}
pScaffold = m_ppPE2Scaffolds[CurIdx];
if(pScaffold->PE2ContigID != ContigID)
	return(0);

if(ppPEScaffold != NULL)
	*ppPEScaffold = pScaffold;
return(CurIdx + 1);
}

int
CPEScaffold::RecurseCluster(int ClusterID,		// identifies the cluster (1..N) to be associated with all scaffolded chrom/contigs
			   int ContigID)						// starting from this chrom/contig (1..N)
{
int CurIdx;
tsPEScaffoldContig *pScaffoldContig;
tsPEScaffold *pScaffold;
if(ContigID == 0)
	return(0);

pScaffoldContig = &m_pScaffoldContigs[ContigID-1];
if(pScaffoldContig->ClusterID > 0)		
	return(0);
pScaffoldContig->ClusterID = ClusterID;

CurIdx = 0;
while((CurIdx = IterPE1s(CurIdx,ContigID,&pScaffold)) > 0)
	if(pScaffold->PE1ContigID != pScaffold->PE2ContigID)
		RecurseCluster(ClusterID,pScaffold->PE2ContigID);	
CurIdx = 0;
while((CurIdx = IterPE2s(CurIdx,ContigID,&pScaffold)) > 0)
	if(pScaffold->PE1ContigID != pScaffold->PE2ContigID)
		RecurseCluster(ClusterID,pScaffold->PE1ContigID);	
return(0);
}


int												// returns largest cluster size
CPEScaffold::IdentifyClusters(void)				// indentify cluster sizes
{
int ClusterID;
int MaxNumClustered;
int Idx;
int *pNumClusterCnts;
tsPEScaffoldContig *pContig;

ClusterID = 0;
pContig = m_pScaffoldContigs;
for(Idx = 0; Idx < m_NumScaffoldContigs; Idx++,pContig++)
	if(pContig->ClusterID == 0)
		{
		ClusterID += 1;
		RecurseCluster(ClusterID,pContig->ContigID);
		}

if(ClusterID)
	{
	// clusters identified, for each cluster count the number of members
	MaxNumClustered = 0;
	pNumClusterCnts = new int [ClusterID];
	memset(pNumClusterCnts,0,sizeof(int) * ClusterID);
	pContig = m_pScaffoldContigs;
	for(Idx = 0; Idx < m_NumScaffoldContigs; Idx++,pContig++)
		if(pContig->ClusterID > 0)
			pNumClusterCnts[pContig->ClusterID-1] += 1;
	pContig = m_pScaffoldContigs;
	for(Idx = 0; Idx < m_NumScaffoldContigs; Idx++,pContig++)
		if(pContig->ClusterID > 0)
			{
			pContig->NumClustered = pNumClusterCnts[pContig->ClusterID-1];
			if(MaxNumClustered < pContig->NumClustered)
				MaxNumClustered = pContig->NumClustered;
			}
	delete pNumClusterCnts;
	}
m_NumClusters = ClusterID;
m_MaxNumClustered = MaxNumClustered;

return(MaxNumClustered);
}

int
CPEScaffold::ReportCorelationships(char *pszOutFile)
{
char *pszPE1Contig;
char *pszPE2Contig;
char *pszSeqID;
int PrevPE1ContigID;
int PrevPE2ContigID;
int NumPEAligned;
int NumSEAligned;
int NumSenseSense;
int NumSenseAnti;
int	RevNumSenseSense;
int	RevNumSenseAnti;
int NumUnpaired;
int NumPaired;
int NumIntraPaired;
int NumInterPaired;
int ClusterID;
int ClusterSize;
int MaxClusterSize;
int NumClusters;

tsPEScaffold *pPEScaffold;
bool bRevMate;
tsPEScaffold *pMateScaffold;
int Idx;
int BuffIdx;
char szBuff[16000];

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Identifying and writing corelations to file: '%s'",pszOutFile);

PrevPE1ContigID = -1;
PrevPE2ContigID = -1;
NumPEAligned = 0;
NumSEAligned = 0;
NumSenseSense = 0;
NumSenseAnti = 0;
RevNumSenseSense = 0;
RevNumSenseAnti = 0;
NumUnpaired = 0;
NumPaired = 0;
BuffIdx = 0;
NumIntraPaired = 0;
NumInterPaired = 0;
MaxClusterSize = 0;
NumClusters = 0;
BuffIdx = sprintf(&szBuff[0],"\"PE1\",\"PE2\",\"NumAligned\",\"NumSenseSense\",\"NumSenseAnti\",\"Paired\",\"Self\",\"RevMate\",\"RevNumSenseSense\",\"RevNumSenseAnti\",\"ClusterID\",\"ClusterSize\"\n");
pPEScaffold = m_pScaffolds;
for(Idx = 0; Idx < m_NumScaffolds; Idx++,pPEScaffold++)
	{
	if(pPEScaffold->PE1ContigID != PrevPE1ContigID || pPEScaffold->PE2ContigID != PrevPE2ContigID) // starting a different scaffold?
		{
		if(PrevPE1ContigID > 0 || PrevPE2ContigID > 0)
			{
			// write out any prev scaffold info here
			if(BuffIdx > (sizeof(szBuff)-500))
				{
				CUtility::SafeWrite(m_hOutFile,szBuff,BuffIdx);
				BuffIdx = 0;
				}
			if(PrevPE1ContigID > 0 && PrevPE2ContigID > 0)
				{
				BuffIdx += sprintf(&szBuff[BuffIdx],"\"%s\",\"%s\",%d,%d,%d,\"Y\",\"%s\",\"%s\",%d,%d,%d,%d\n",
							pszPE1Contig,pszPE2Contig,NumPEAligned,NumSenseSense,NumSenseAnti,PrevPE1ContigID==PrevPE2ContigID ? "Y" : "N",bRevMate ? "Y" : "N",
							RevNumSenseSense,RevNumSenseAnti,ClusterID,ClusterSize);
				if(PrevPE1ContigID==PrevPE2ContigID)
					NumIntraPaired += 1;
				else
					NumInterPaired += 1;
				}
			else
				{
				if(PrevPE1ContigID > 0)
					BuffIdx += sprintf(&szBuff[BuffIdx],"\"%s\",\"%s\",%d,%d,%d,\"N\",\"N\",\"N\",0,0,%d,%d\n",pszPE1Contig,"N/A",NumSEAligned,0,0,ClusterID,ClusterSize);
				else
					BuffIdx += sprintf(&szBuff[BuffIdx],"\"%s\",\"%s\",%d,%d,%d,\"N\",\"N\",\"N\",0,0,%d,%d\n","N/A",pszPE2Contig,NumSEAligned,0,0,ClusterID,ClusterSize);
				}
			}

		NumPEAligned = 0;					// number aligning as paired ends
		NumSEAligned = 0;					// number aligning as single ended only
		NumSenseSense = 0;					// number aligning sense to sense
		NumSenseAnti = 0;					// number aligning sense to antisense (or antisense to sense)
		if(pPEScaffold->PE1ContigID > 0 && pPEScaffold->PE2ContigID > 0)
			{
			NumPEAligned = 1;
			if(pPEScaffold->PE1Sense == pPEScaffold->PE2Sense)
				NumSenseSense = 1;
			else
				NumSenseAnti = 1;
			NumPaired += 1;
			}
		else
			{
			NumSEAligned = 1;
			NumUnpaired += 1;
			}
		pszSeqID = GetSeqName(pPEScaffold->PE12SeqID);
		if(pPEScaffold->PE1ContigID)
			{
			pszPE1Contig = GetContigName(pPEScaffold->PE1ContigID);
			ClusterID = m_pScaffoldContigs[pPEScaffold->PE1ContigID-1].ClusterID;
			ClusterSize = m_pScaffoldContigs[pPEScaffold->PE1ContigID-1].NumClustered;
			}
		else
			pszPE1Contig = (char *)"N/A";

		if(pPEScaffold->PE2ContigID)
			{
			pszPE2Contig = GetContigName(pPEScaffold->PE2ContigID);
			ClusterID = m_pScaffoldContigs[pPEScaffold->PE2ContigID-1].ClusterID;
			ClusterSize = m_pScaffoldContigs[pPEScaffold->PE2ContigID-1].NumClustered;
			}
		else
			pszPE2Contig = (char *)"N/A";

		if(ClusterSize > MaxClusterSize)
			MaxClusterSize = ClusterSize;
		if(ClusterSize >= 2)
			NumClusters += 1;

		PrevPE1ContigID = pPEScaffold->PE1ContigID;
		PrevPE2ContigID = pPEScaffold->PE2ContigID;
	
		// if not to self then check if the PE2 chrom has any PE2 linking back to this PE1 and count linking sense/antisense
		bRevMate = false;
		RevNumSenseSense = 0;
		RevNumSenseAnti = 0;
		if(PrevPE1ContigID != PrevPE2ContigID && (PrevPE1ContigID > 0 && PrevPE2ContigID > 0))
			{
			pMateScaffold = LocateMateScaffold(PrevPE2ContigID,PrevPE1ContigID);
			if(pMateScaffold != NULL)
				{
				do {
					if(pMateScaffold->PE1Sense == pMateScaffold->PE2Sense)
						RevNumSenseSense += 1;
					else
						RevNumSenseAnti += 1;
					if(pMateScaffold->PEScafoldID == m_NumScaffolds)
						break;
					pMateScaffold += 1;
					}
				while(pMateScaffold->PE1ContigID == PrevPE2ContigID && pMateScaffold->PE2ContigID == PrevPE1ContigID);

				bRevMate = true;
				}
			}
		continue;
		}

	// same pair of chromosomes so accumulate counts
	if(pPEScaffold->PE1ContigID > 0 && pPEScaffold->PE2ContigID > 0)
		{
		NumPEAligned += 1;
		if(pPEScaffold->PE1Sense == pPEScaffold->PE2Sense)
			NumSenseSense += 1;
		else
			NumSenseAnti += 1;
		NumPaired += 1;
		}
	else
		{
		NumSEAligned += 1;
		NumUnpaired += 1;
		}
	}

if(BuffIdx > 0 || NumPEAligned > 0 || NumSEAligned > 0)
	{
	if(PrevPE1ContigID > 0 && PrevPE2ContigID > 0)		// if PE1 and PE2 aligned then paired
		{
		BuffIdx += sprintf(&szBuff[BuffIdx],"\"%s\",\"%s\",%d,%d,%d,\"Y\",\"%s\",\"%s\",%d,%d,%d,%d\n",
							pszPE1Contig,pszPE2Contig,NumPEAligned,NumSenseSense,NumSenseAnti,PrevPE1ContigID==PrevPE2ContigID ? "Y" : "N",bRevMate ? "Y" : "N",
							RevNumSenseSense,RevNumSenseAnti,ClusterID,ClusterSize);
		if(PrevPE1ContigID==PrevPE2ContigID)
			NumIntraPaired += 1;
		else
			NumInterPaired += 1;
		}
	else
		{
		if(PrevPE1ContigID > 0)
			BuffIdx += sprintf(&szBuff[BuffIdx],"\"%s\",\"%s\",%d,%d,%d,\"N\",\"N\",\"N\",0,0,%d,%d\n",pszPE1Contig,"N/A",NumSEAligned,0,0,ClusterID,ClusterSize);
		else
			BuffIdx += sprintf(&szBuff[BuffIdx],"\"%s\",\"%s\",%d,%d,%d,\"N\",\"N\",\"N\",0,0,%d,%d\n","N/A",pszPE2Contig,NumSEAligned,0,0,ClusterID,ClusterSize);
		}

	CUtility::SafeWrite(m_hOutFile,szBuff,BuffIdx);
	}

close(m_hOutFile);
m_hOutFile = -1;

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Completed writing corelations (%d paired - intra %d inter %d - with %d orphaned) to file",NumPaired,NumIntraPaired,NumInterPaired,NumUnpaired);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Max contigs in any cluster: %d, number of clusters with 2 or more contigs: %d",MaxClusterSize,NumClusters);
return(0);
}

int
CPEScaffold::Process(int PMode,		// processing mode
		char *pszSeqIDTerm,			// pair sequence identifiers until this terminating character(s) - defaults to none terminating
		char *pszInPE1File,			// input PE1 file
		char *pszInPE2File,			// input PE2 file
		char *pszOutFile)			// output corelations file
{
int Rslt;
size_t memreq;
Init();

if(pszSeqIDTerm != NULL && pszSeqIDTerm[0] != '\0')
	{
	strncpy(m_szSeqIDTermChrs,pszSeqIDTerm,sizeof(m_szSeqIDTermChrs)-1);
	m_szSeqIDTermChrs[sizeof(m_szSeqIDTermChrs)-1] = '\0';
	}
else
	m_szSeqIDTermChrs[0] = '\0';

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

memreq = (size_t)(sizeof(tsPEScaffoldContig) * cAllocContigNames);	
#ifdef _WIN32
m_pScaffoldContigs = (tsPEScaffoldContig *) malloc(memreq);	// initial and perhaps the only allocation

if(m_pScaffoldContigs == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadReads: Memory allocation of %lld bytes for m_pScaffoldContigs - %s",(INT64)memreq,strerror(errno));
	Reset();
	return(eBSFerrMem);
	}
#else
	// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
m_pScaffoldContigs = (tsPEScaffoldContig *)mmap(NULL,memreq, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
if(m_pScaffoldContigs == MAP_FAILED)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadReads: Memory allocation of %lld bytes for m_pScaffoldContigs through mmap()  failed - %s",(INT64)memreq,strerror(errno));
	m_pScaffoldContigs = NULL;
	Reset();
	return(eBSFerrMem);
	}
#endif
m_AllocdScaffoldContigsMem = memreq;
m_AllocdNumScaffoldContigs = cAllocContigNames;
m_NumScaffoldContigs = 0;

m_pHashContigs = new int [cHashSize];
if(m_pHashContigs == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadReads: Memory allocation of %d bytes for m_pHashContigs - %s",cHashSize * sizeof(int),strerror(errno));
	Reset();
	return(eBSFerrMem);
	}
memset(m_pHashContigs,0,sizeof(int) * cHashSize);

memreq = (size_t)(sizeof(tsPEScaffold) * cAllocScafolds);	
#ifdef _WIN32
m_pScaffolds = (tsPEScaffold *) malloc(memreq);	// initial and perhaps the only allocation

if(m_pScaffolds == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadReads: Memory allocation of %lld bytes for m_pScaffolds - %s",(INT64)memreq,strerror(errno));
	Reset();
	return(eBSFerrMem);
	}
#else
	// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
m_pScaffolds = (tsPEScaffold *)mmap(NULL,memreq, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
if(m_pScaffolds == MAP_FAILED)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadReads: Memory allocation of %lld bytes for m_pScaffolds through mmap()  failed - %s",(INT64)memreq,strerror(errno));
	m_pScaffolds = NULL;
	Reset();
	return(eBSFerrMem);
	}
#endif
m_AllocdScaffoldsMem = memreq;
m_AllocdNumScaffolds = cAllocScafolds;
m_NumScaffolds = 0;

memreq = (size_t)(sizeof(tsPEIdent) * cAllocPENames);	
#ifdef _WIN32
m_pPEIdents = (tsPEIdent *) malloc(memreq);	// initial and perhaps the only allocation

if(m_pPEIdents == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadReads: Memory allocation of %lld bytes for m_pPEIdents - %s",(INT64)memreq,strerror(errno));
	Reset();
	return(eBSFerrMem);
	}
#else
	// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
m_pPEIdents = (tsPEIdent *)mmap(NULL,memreq, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
if(m_pPEIdents == MAP_FAILED)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadReads: Memory allocation of %lld bytes for m_pPEIdents through mmap()  failed - %s",(INT64)memreq,strerror(errno));
	m_pPEIdents = NULL;
	Reset();
	return(eBSFerrMem);
	}
#endif
m_AllocdPEIdentsMem = memreq;
m_AllocdNumPEIdents = cAllocPENames;
m_NumPEIdents = 0;

m_pHashPEIdents = new int [cHashSize];
if(m_pHashPEIdents == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadReads: Memory allocation of %d bytes for m_pHashPEIdents - %s",cHashSize * sizeof(int),strerror(errno));
	Reset();
	return(eBSFerrMem);
	}
memset(m_pHashPEIdents,0,sizeof(int) * cHashSize);

Rslt = LoadSAM(false,pszInPE1File);
if(Rslt < 1)
	return(Rslt);

Rslt = LoadSAM(true,pszInPE2File);
if(Rslt < 1)
	return(Rslt);

if(m_NumScaffolds > 1)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Sorting %d scaffolds by PE1ContigID.PE2ContigID ascending",m_NumScaffolds);
	m_qsort.qsort(m_pScaffolds,m_NumScaffolds,sizeof(tsPEScaffold),SortScaffolds);
	tsPEScaffold *pScaffold;
	int Idx;
	pScaffold = m_pScaffolds;
	for(Idx=1;Idx <= m_NumScaffolds; Idx++,pScaffold++)
		pScaffold->PEScafoldID = Idx;
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Sorting scaffolds by PE1ContigID.PE2ContigID ascending completed");
	}

// number of scaffolds known, allocate to hold PE2 indexes and then sort
if(m_NumScaffolds > 0)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Create index and sorting %d scaffolds by PE2ContigID.PE1ContigID ascending",m_NumScaffolds);
	memreq = (size_t)(sizeof(tsPEScaffold *) * m_NumScaffolds);	
#ifdef _WIN32
	m_ppPE2Scaffolds = (tsPEScaffold **) malloc(memreq);	// initial and perhaps the only allocation

	if(m_ppPE2Scaffolds == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Memory allocation of %lld bytes for m_ppPE2Scaffolds - %s",(INT64)memreq,strerror(errno));
		Reset();
		return(eBSFerrMem);
		}
#else
		// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
	m_ppPE2Scaffolds = (tsPEScaffold **)mmap(NULL,memreq, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
	if(m_ppPE2Scaffolds == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Memory allocation of %lld bytes for m_ppPE2Scaffolds through mmap()  failed - %s",(INT64)memreq,strerror(errno));
		m_ppPE2Scaffolds = NULL;
		Reset();
		return(eBSFerrMem);
		}
#endif
	m_AllocdPE2ScaffoldsMem = memreq;
	// initialise index and then sort
	for(int Idx = 0; Idx < m_NumScaffolds; Idx++)
		m_ppPE2Scaffolds[Idx] = &m_pScaffolds[Idx];
	if(m_NumScaffolds > 1)
		m_qsort.qsort(m_ppPE2Scaffolds,m_NumScaffolds,sizeof(tsPEScaffold **),SortPE2Scaffolds);
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Sorting scaffolds by PE2ContigID.PE1ContigID ascending completed");
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Identifying clusters of scaffolded chrom/contigs");
IdentifyClusters();
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Cluster identification of scaffolded chrom/contigs completed");
ReportCorelationships(pszOutFile);
Reset();
return(0);
}

// Sort scaffolds by PE1ContigID.PE2ContigID ascending
int
CPEScaffold::SortScaffolds(const void *arg1, const void *arg2)
{
tsPEScaffold *pEl1 = (tsPEScaffold *)arg1;
tsPEScaffold *pEl2 = (tsPEScaffold *)arg2;

if(pEl1->PE1ContigID > pEl2->PE1ContigID)
	return(1);
if(pEl1->PE1ContigID < pEl2->PE1ContigID)
	return(-1);
if(pEl1->PE2ContigID > pEl2->PE2ContigID)
	return(1);
if(pEl1->PE2ContigID < pEl2->PE2ContigID)
	return(-1);
return(0);
}

// Sort scaffolds by PE2ContigID.PE1ContigID ascending
int
CPEScaffold::SortPE2Scaffolds(const void *arg1, const void *arg2)
{
tsPEScaffold *pEl1 = *(tsPEScaffold **)arg1;
tsPEScaffold *pEl2 = *(tsPEScaffold **)arg2;

if(pEl1->PE2ContigID > pEl2->PE2ContigID)
	return(1);
if(pEl1->PE2ContigID < pEl2->PE2ContigID)
	return(-1);
if(pEl1->PE1ContigID > pEl2->PE1ContigID)
	return(1);
if(pEl1->PE1ContigID < pEl2->PE1ContigID)
	return(-1);
return(0);
}