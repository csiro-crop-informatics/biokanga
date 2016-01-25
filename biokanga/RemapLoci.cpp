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

// RemapLoci.cpp : Defines the entry point for the console application.

#include "stdafx.h"
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "../libbiokanga/commhdrs.h"

#include "biokanga.h"

#include "RemapLoci.h"


int
RemapLociProcess(int PMode,			// processing mode
				 int FType,			// alignment file type
				char *pszInAlignFile,	// alignment file with loci to be remapped
				char *pszInBEDFile,     // BED file containing loci remapping
				char *pszRemappedFile);	// write remapped alignments to this file

#ifdef _WIN32
int RemapLoci(int argc, char* argv[])
{
// determine my process name
_splitpath(argv[0],NULL,NULL,gszProcName,NULL);
#else
int
RemapLoci(int argc, char** argv)
{
// determine my process name
CUtility::splitpath((char *)argv[0],NULL,gszProcName);
#endif
int iScreenLogLevel;		// level of screen diagnostics
int iFileLogLevel;			// level of file diagnostics
char szLogFile[_MAX_PATH];	// write diagnostics to this file

int Rslt;

int PMode;				// processing mode
int FType;					// expected input element file type - auto, CSV, BED or SAM/BAM

char szInLociFile[_MAX_PATH];	// input element loci from this file
char szInBEDFile[_MAX_PATH];	// input bed file containing gene features
char szRemappedFile[_MAX_PATH];	// write remapped alignments to this file

char szSQLiteDatabase[_MAX_PATH];	// results summaries to this SQLite file
char szExperimentName[cMaxDatasetSpeciesChrom+1];			// experiment name
char szExperimentDescr[1000];		// describes experiment

// command line args
struct arg_lit  *help    = arg_lit0("hH","help",                "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"diagnostics log file");

struct arg_int *pmode = arg_int0("m","mode","<int>",		    "processing mode: 0 - full length attributes, 1 - 5' attributes associations, 2 - 5' dyad attribute associations (default = 0)");

struct arg_int *ftype = arg_int0("t","filetype","<int>",		"input element file format: 0 - auto, 1 - CSV (not currently supported), 2 - BED, 3 - SAM/BAM (default = 0)");
struct arg_file *InLociFile = arg_file1("i","inloci","<file>",	"input alignments file with loci to be remapped (BED, SAM/BAM) file");
struct arg_file *InBEDFile = arg_file1("I","inbed","<file>",	"input BED file containing remapping loci");
struct arg_file *RemappedFile = arg_file1("o","output","<file>", "write remapped alignments to this file, same format as input alignment file");
struct arg_file *summrslts = arg_file0("q","sumrslts","<file>",		"Output results summary to this SQLite3 database file");
struct arg_str *experimentname = arg_str0("w","experimentname","<str>",		"experiment name SQLite3 database file");
struct arg_str *experimentdescr = arg_str0("W","experimentdescr","<str>",	"experiment description SQLite3 database file");

struct arg_end *end = arg_end(200);

void *argtable[] = {help,version,FileLogLevel,LogFile,
					summrslts,experimentname,experimentdescr,
					pmode,ftype,InLociFile,InBEDFile,RemappedFile,
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
	PMode = pmode->count ? pmode->ival[0] : 0;
	if(PMode < 0 || PMode > 0)
		{
		printf("\nError: Processing mode '-m%d' specified outside of range %d..%d",PMode,0,(int)0);
		exit(1);
		}

	FType = ftype->count ? ftype->ival[0] : 0;
	if(FType < 0 || FType >= 4)
		{
		printf("\nError: Expected input element file format '-t%d' specified outside of range %d..%d",FType,0,3);
		exit(1);
		}

	if(FType == 1)
		{
		printf("\nError: Input CSV format alignments not supported");
		exit(1);
		}

	if(InBEDFile->count)
		{
		strncpy(szInBEDFile,InBEDFile->filename[0],_MAX_PATH);
		szInBEDFile[_MAX_PATH-1] = '\0';
		}
	else
		{
		printf("Error: No BED gene/feature file specified");
		exit(1);
		}

	strncpy(szInLociFile,InLociFile->filename[0],_MAX_PATH);
	szInLociFile[_MAX_PATH-1] = '\0';

	strncpy(szRemappedFile,RemappedFile->filename[0],_MAX_PATH);
	szRemappedFile[_MAX_PATH-1] = '\0';

// show user current resource limits
#ifndef _WIN32
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Resources: %s",CUtility::ReportResourceLimits());
#endif

	const char *pszDescr;
	switch(PMode) {
		case 0:
			pszDescr = "Remapping alignment loci";
			break;
		}

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Processing mode is : '%s'",pszDescr);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Input element loci file: '%s'",szInLociFile);
	switch(FType) {
		case 0:
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Auto-classify input element file as either BED or SAM/BAM");
			break;

		case 2:
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Expecting input element file to be BED format");
			break;

		case 3:
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Expecting input element file to be SAM format");
			break;
		}
	
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Remapped alignment locii to file: '%s'",szRemappedFile);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Input BED remapping locii file: '%s'",szInBEDFile);

	if(szExperimentName[0] != '\0')
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"This processing reference: %s",szExperimentName);

#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	gStopWatch.Start();
	Rslt = RemapLociProcess(PMode,FType,szInLociFile,szInBEDFile,szRemappedFile);
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


int
RemapLociProcess(int PMode,			// processing mode
				 int FType,			// alignment file type
				char *pszInAlignFile,	// alignment file with loci to be remapped
				char *pszInBEDFile,     // BED file containing loci remapping
				char *pszRemappedFile)	// write remapped alignments to this file
{
CRemapLoci Remapper;
return(Remapper.RemapLocii(PMode,FType,pszInAlignFile,pszInBEDFile,pszRemappedFile));
}

CRemapLoci::CRemapLoci()
{
m_pMappingBED = NULL;
m_pInBAMfile = NULL;
m_pOutBAMfile = NULL;
m_pszOutBuff = NULL;
m_hOutFile = -1;
}


CRemapLoci::~CRemapLoci()
{
if(m_pMappingBED != NULL)
	delete m_pMappingBED;
if(m_pInBAMfile != NULL)
	delete m_pInBAMfile;
if(m_pOutBAMfile != NULL)
	delete m_pOutBAMfile;
if(m_pszOutBuff != NULL)
	delete m_pszOutBuff;
}

void
CRemapLoci::Reset(void)
{
if(m_pMappingBED != NULL)
	{
	delete m_pMappingBED;
	m_pMappingBED = NULL;
	}
if(m_pInBAMfile != NULL)
	{
	delete m_pInBAMfile;
	m_pInBAMfile = NULL;
	}
if(m_pOutBAMfile != NULL)
	{
	delete m_pOutBAMfile;
	m_pOutBAMfile = NULL;
	}

if(m_hOutFile != -1)
	{
	close(m_hOutFile);
	m_hOutFile = -1;
	}

if(m_pszOutBuff != NULL)
	{
	delete m_pszOutBuff;
	m_pszOutBuff = NULL;
	}

m_hOutFile = -1;
m_OutBuffIdx = 0;
m_AllocOutBuff = 0;	
}

int
CRemapLoci::RemapLocii(int PMode,			// processing mode
				 int FType,			// alignment file type
				char *pszInAlignFile,	// alignment file with loci to be remapped
				char *pszInBEDFile,     // BED file containing loci remapping
				char *pszRemappedFile)	// write remapped alignments to this file
{
int Rslt;
etClassifyFileType FileType;

Reset();

// load remapping locii BED file
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Loading BED file containing loci remapping '%s'",pszInBEDFile);
if((m_pMappingBED = new CBEDfile) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CBEDfile");
	Reset();
	return(eBSFerrObj);
	}

if((Rslt=m_pMappingBED->Open(pszInBEDFile))!=eBSFSuccess)
	{
	while(m_pMappingBED->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,m_pMappingBED->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open BED file '%s'",pszInBEDFile);
	Reset();
	return(eBSFerrObj);
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Opened BED file '%s'",pszInBEDFile);


if(FType == 0)
	FileType = CUtility::ClassifyFileType(pszInAlignFile);
else
	FileType = (etClassifyFileType)(FType - 1);

switch(FileType) {
	case eCFTopenerr:		// unable to open file for reading
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open file: '%s'",pszInAlignFile);
		return(eBSFerrOpnFile);

	case eCFTlenerr:		// file length is insufficent to classify type
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to classify file type (insufficent data points): '%s'",pszInAlignFile);
		return(eBSFerrFileAccess);

	case eCFTunknown:		// unable to reliably classify
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to reliably classify file type: '%s'",pszInAlignFile);
		return(eBSFerrFileType);

	case eCFTCSV:			// file has been classified as being CSV but format is unsupported
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"File type is CSV which is not supported: '%s'",pszInAlignFile);
		return(eBSFerrFileType);

	case eCFTBED:			// file has been classified as being BED
		if((m_pszOutBuff = new char [cAllocOutBuff]) == NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate memory for output buffering");
			return(eBSFerrMem);
			}
		m_AllocOutBuff = cAllocOutBuff;
		m_OutBuffIdx = 0;

		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Parsing BED file: '%s'",pszInAlignFile);
		if((Rslt = RemapBEDLocii(pszInAlignFile,pszRemappedFile)) < 0)
			{	
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Processing errors");
			Reset();
			return(Rslt);
			}
		break;

	case eCFTSAM:			// file has been classified as being SAM
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing SAM/BAM file: '%s'",pszInAlignFile);
		if((Rslt = RemapSAMLocii(pszInAlignFile,pszRemappedFile)) < 0)
			{	
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Processing errors");
			Reset();
			return(Rslt);
			}
		break;
	}
Reset();
return(eBSFSuccess);
}



int
CRemapLoci::RemapBEDLocii(char *pszInAlignFile,		// BED alignment file with loci to be remapped
				char *pszRemappedFile)				// write remapped alignments to this file
{
teBSFrsltCodes Rslt;
int hOutFile;
CBEDfile *pBedFile;
int RelStartLoci;
char szTitle[128];
char szChrom[128];
char szFeatName[128];
int CurFeatureID;
int NumEls;
int NumFeatures;
int ContigID;

if((pBedFile = new CBEDfile)==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CBEDfile");
	Reset();
	return(eBSFerrObj);
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Starting to load from %s",pszInAlignFile);
if((Rslt=pBedFile->Open(pszInAlignFile,eBTAnyBed)) !=eBSFSuccess)
	{
	while(pBedFile->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,pBedFile->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open '%s' for processing",pszInAlignFile);
	delete pBedFile;
	Reset();
	return(eBSFerrOpnFile);
	}
NumFeatures = pBedFile->GetNumFeatures();
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Completed load from '%s' of %d features",pszInAlignFile,NumFeatures);
if(NumFeatures < 1)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Nothing to do, no features to map read loci against!");
	delete pBedFile;
	return(eBSFerrNoFeatures);
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Completed load from %s",pszInAlignFile);

#ifdef _WIN32
if((hOutFile = open(pszRemappedFile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
	if((hOutFile = open(pszRemappedFile, O_RDWR | O_CREAT | O_TRUNC, S_IREAD | S_IWRITE ))==-1)
#endif
	{
	gDiagnostics.DiagOutMsgOnly(eDLFatal,"Unable to create or truncate  output file %s error: %s",pszRemappedFile,strerror(errno));
	return(-1);
	}

NumEls = 0;
CurFeatureID = 0;
Rslt = eBSFSuccess;
pBedFile->GetTitle(sizeof(szTitle),szTitle);

m_OutBuffIdx = sprintf(m_pszOutBuff,"track type=bed name=\"Remaploci\" description=\"Alignment loci remapped to different coordinate system\"\n");
while(Rslt >= eBSFSuccess && (CurFeatureID = pBedFile->GetNextFeatureID(CurFeatureID)) > 0)
	{
	if(!(NumEls % 5000000))
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Parsed %d elements",NumEls);
	pBedFile->GetFeature(CurFeatureID,	szFeatName);

	if((ContigID = m_pMappingBED->LocateFeatureIDbyName(szFeatName)) < 1)
		continue;
	m_pMappingBED->GetFeature(ContigID,NULL,szChrom,&RelStartLoci);
	Rslt = (teBSFrsltCodes)pBedFile->GetRemappedBEDFormatedFeature(CurFeatureID,m_AllocOutBuff - m_OutBuffIdx,&m_pszOutBuff[m_OutBuffIdx],szChrom,szFeatName,RelStartLoci);
	if(Rslt < 1)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"GetRemappedBEDFormatedFeature failed");
		delete pBedFile;
		Reset();
		return(Rslt);
		}
	m_OutBuffIdx += (int)Rslt;
	m_OutBuffIdx += sprintf(&m_pszOutBuff[m_OutBuffIdx],"\n");
	if((m_OutBuffIdx + 10000) > m_AllocOutBuff)
		{
		CUtility::SafeWrite(m_hOutFile,m_pszOutBuff,m_OutBuffIdx);
		m_OutBuffIdx = 0;
		}
	NumEls++;
	}
delete pBedFile;
if(m_OutBuffIdx)
	{
	CUtility::SafeWrite(m_hOutFile,m_pszOutBuff,m_OutBuffIdx);
	m_OutBuffIdx = 0;
	}
#ifdef _WIN32
	_commit(m_hOutFile);
#else
	fsync(m_hOutFile);
#endif
close(m_hOutFile);
m_hOutFile = -1;
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Parsed %d elements",NumEls);
return(Rslt >= 0 ? NumEls : Rslt);
}

char *
CRemapLoci::TrimWhitespace(char *pTxt)
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

bool								// true if file to be generated compressed with gzopen/gzwrite/gzclose
CRemapLoci::FileReqWriteCompr(char *pszFile) // If last 3 chars of file name is ".gz" then this file is assumed to require compression
{
int Len;
if(pszFile == NULL || pszFile[0] == '\0')
	return(false);
if((Len = (int)strlen(pszFile)) < 4)
	return(false);
return(stricmp(".gz",&pszFile[Len-3]) == 0 ? true : false);
}


int
CRemapLoci::RemapSAMLocii(char *pszInAlignFile,		// SAM or BAM alignment file with loci to be remapped
				char *pszRemappedFile)		// write remapped alignments to this file
{
teBSFrsltCodes Rslt;
int NumParsedElLines;
int NumAcceptedEls;
int NumUnmappedEls;
int NumMappedChroms;
int LineLen;
char szLine[cMaxReadLen  * 3];				// buffer input lines
bool bFirstAlignment;
tsBAMalign ProvBAMalignment;
tsBAMalign AcceptedBAMalignment;

char szGenome[100];
char szChrom[100];
char szContig[100];
int ChromID;
int PrevChromID;
int ChromLen;
int ContigID;
int PrevContigID;
int ContigLen;
int LastChromContigID;

int Tmp = 0;

char *pTxt;

	// if SAM output format then could be BGZF compressed BAM; use the extension to determine which...
	// if extension is '.bam' then BGZF compressed BAM, any other extension is for SAM
int Len;
int SAMFormat = 0;
Len = (int)strlen(pszRemappedFile);
if(Len > 5)
	{
	if(!stricmp(".bam",&pszRemappedFile[Len-4]))
		SAMFormat = 1;
	}

// open SAM for reading
if(pszInAlignFile == NULL || *pszInAlignFile == '\0')
	{
	Reset();
	return(eBSFerrParams);
	}

if((m_pInBAMfile = new CSAMfile) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"RemapSAMLocii: Unable to instantiate class CSAMfile");
	Reset();
	return(eBSFerrInternal);
	}

if((Rslt = (teBSFrsltCodes)m_pInBAMfile->Open(pszInAlignFile)) != eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"RemapSAMLocii: Unable to open SAM/BAM format file %s",pszInAlignFile);
	delete m_pInBAMfile;
	m_pInBAMfile = NULL;
	Reset();
	return((teBSFrsltCodes)Rslt);
	}

if((m_pOutBAMfile = new CSAMfile) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"RemapSAMLocii: Unable to instantiate class CSAMfile");
	delete m_pInBAMfile;
	m_pInBAMfile = NULL;
	Reset();
	return(eBSFerrInternal);
	}

eSAMFileType FileType;
bool bgzOutFile = FileReqWriteCompr(pszRemappedFile);

switch(SAMFormat) {
	case 0:			// output SAM
		if(bgzOutFile)
			FileType = eSFTSAMgz;
		else
			FileType = eSFTSAM;
		break;
	case 1:				// output as BAM compressed with bgzf
		FileType = eSFTBAM_BAI;
		break;
	}

if((Rslt = (teBSFrsltCodes)m_pOutBAMfile->Create(FileType,pszRemappedFile,6,(char *)cpszProgVer)) < eBSFSuccess) // defaulting to compression level 6 if compressed BAM 
	{
	delete m_pInBAMfile;
	m_pInBAMfile = NULL;
	delete m_pOutBAMfile;
	m_pInBAMfile = NULL;
	Reset();
	return(Rslt);
	}

NumParsedElLines = 0;
NumAcceptedEls = 0;
NumUnmappedEls = 0;
bFirstAlignment = true;
Rslt = eBSFSuccess;
PrevContigID = 0;
PrevChromID = 0;
NumMappedChroms = 0;
while(Rslt >= eBSFSuccess && (LineLen = m_pInBAMfile->GetNxtSAMline(szLine)) > 0)
	{
	NumParsedElLines += 1;
	if(!(NumParsedElLines % 5000000) || NumParsedElLines == 1)
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Parsed %d element lines",NumParsedElLines);

	szLine[sizeof(szLine)-1] = '\0';
	pTxt = TrimWhitespace(szLine);
	if(*pTxt=='\0')			// simply slough lines which are just whitespace
		continue;

	if(*pTxt== '@')			// if a reference sequence dictionary entry then check if an entry in the mapping loci BED
		{
		if(pTxt[1] != 'S' && pTxt[2] != 'Q')
			continue;
		if(3 != sscanf(&pTxt[3]," AS:%s SN:%s LN:%d",szGenome,szContig,&ContigLen))
			continue;
		if((ContigID = m_pMappingBED->LocateFeatureIDbyName(szContig)) < 1)
			continue;
		if(ContigID == PrevContigID)
			continue;
		PrevContigID = ContigID;
		ChromID = m_pMappingBED->GetFeatureChromID(ContigID);
		if(ChromID == PrevChromID)
			continue;
		PrevChromID = ChromID;
		m_pMappingBED->GetChromosome(ChromID,szChrom,NULL,NULL,&LastChromContigID);
		m_pMappingBED->GetFeature(LastChromContigID,NULL,NULL,NULL,&ChromLen);
		m_pOutBAMfile->AddRefSeq(szGenome,szChrom,ChromLen);
		NumMappedChroms += 1;
		continue;
		}	
	if(bFirstAlignment)
		{
		m_pOutBAMfile->StartAlignments();	
		bFirstAlignment = false;
		}
	if(!NumMappedChroms)
		break;

	// primary interest is in the reference chromname, startloci, length
	if((Rslt = (teBSFrsltCodes)m_pOutBAMfile->ParseSAM2BAMalign(pTxt,&ProvBAMalignment,m_pMappingBED)) < eBSFSuccess)
		{
		NumUnmappedEls += 1;
		if(Rslt == eBSFerrFeature)
			Rslt = eBSFSuccess;
	    continue;
		}

		// check if read has been mapped, if not then slough ...
	if(ProvBAMalignment.refID == -1 || (ProvBAMalignment.flag_nc >> 16) & 0x04 || ProvBAMalignment.cigar[0] == '*')	// set if unmapped or Cigar is unknown
		{
		NumUnmappedEls += 1;
	    continue;
		}

	if(NumAcceptedEls > 0)
		if((Rslt = (teBSFrsltCodes)m_pOutBAMfile->AddAlignment(&AcceptedBAMalignment,false)) < eBSFSuccess)
			continue;
	AcceptedBAMalignment = ProvBAMalignment;
	NumAcceptedEls += 1;
	}

if(Rslt >= eBSFSuccess && NumAcceptedEls > 0)
	Rslt = (teBSFrsltCodes)m_pOutBAMfile->AddAlignment(&AcceptedBAMalignment,true);
if(Rslt >= eBSFSuccess)
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Parsed %d element lines, unmapped %d, accepted %d, ",NumParsedElLines,NumUnmappedEls,NumAcceptedEls);
if(m_pInBAMfile != NULL)
	{
	m_pInBAMfile->Close();
	delete m_pInBAMfile;
	m_pInBAMfile = NULL;
	}
if(m_pOutBAMfile != NULL)
	{
	m_pOutBAMfile->Close();
	delete m_pOutBAMfile;
	m_pOutBAMfile = NULL;
	}
return(Rslt >= 0 ? NumAcceptedEls : Rslt);
}