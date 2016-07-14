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
#include "FastaNxx.h"

int 
Process(etNxxPMode Mode,			// processing mode - 0  N50 distributions, 1 n-mer count distributions
        int MinLength,				// core elements must be of at least this length
		int MaxLength,				// will be truncated to this length
		int NumBins,				// when generating length distributions then use this many bins - 0 defaults to using 1000
		int BinDelta,				// when generating length distributions then each bin holds this length delta - 0 defaults to auto determine from NunBins and longest sequence length
		int NumInputFiles,			// number of input sequence files
		char **pszInFastaFile,		// names of input sequence files (wildcards allowed)
		char *pszRsltsFile);		// file to write results out into


#ifdef _WIN32
int fasta2nxx(int argc, char* argv[])
{
// determine my process name
_splitpath(argv[0],NULL,NULL,gszProcName,NULL);
#else
int
fasta2nxx(int argc, char** argv)
{
// determine my process name
CUtility::splitpath((char *)argv[0],NULL,gszProcName);
#endif

int iFileLogLevel;			// level of file diagnostics
int iScreenLogLevel;		// level of file diagnostics
char szLogFile[_MAX_PATH];	// write diagnostics to this file
int Rslt = 0;   			// function result code >= 0 represents success, < 0 on failure

etNxxPMode PMode;				// processing mode
int Idx;

int iMinLength;				// sequences must be of at least this length
int iMaxLength;				// and no longer than this length

int NumBins;				// when generating length distributions then use this many bins - 0 defaults to using 1000
int BinDelta;				// when generating length distributions then each bin holds this length delta - 0 defaults to auto determine from NunBins and longest sequence length

int NumInputFiles;							// number of input sequence files
char *pszInFastaFile[cMaxInFileSpecs];		// names of input sequence files (wildcards allowed)

char szRsltsFile[_MAX_PATH];	// output stats to this file

char szSQLiteDatabase[_MAX_PATH];	// results summaries to this SQLite file
char szExperimentName[cMaxDatasetSpeciesChrom+1];			// experiment name
char szExperimentDescr[1000];		// describes experiment

// command line args
struct arg_lit  *help    = arg_lit0("h","help",                 "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to screen and logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"diagnostics log file");

struct arg_int *pmode = arg_int0("m","mode","<int>",		    "Processing mode:  0 - N50 stats, 1 - mono/di/tri-mer distributions");

struct arg_file *pinputfiles = arg_filen("i","infasta","<file>",1,cMaxInFileSpecs,"input fasta file(s)");


struct arg_file *RsltsFile = arg_file0("o","output","<file>",	"output file");
struct arg_int  *MinLength = arg_int0("l","minlength","<int>",	"minimum fasta sequence length (default 10)");
struct arg_int  *MaxLength = arg_int0("L","maxlength","<int>",	"truncate fasta sequence length (default 200)");

struct arg_int  *numbins = arg_int0("b","numbins","<int>",	"when generating length distributions then use this many bins (defaults to 1000, range 10..10000)");
struct arg_int  *bindelta = arg_int0("B","bindelta","<int>","when generating length distributions then each bin holds this length delta (default 0 for auto-determination, range 1,2,5,10,25,50,100,250,500 or 1000)");

struct arg_file *summrslts = arg_file0("q","sumrslts","<file>",		"Output results summary to this SQLite3 database file");
struct arg_str *experimentname = arg_str0("w","experimentname","<str>",		"experiment name SQLite3 database file");
struct arg_str *experimentdescr = arg_str0("W","experimentdescr","<str>",	"experiment description SQLite3 database file");

struct arg_end *end = arg_end(100);

void *argtable[] = {help,version,FileLogLevel,LogFile,
					summrslts,experimentname,experimentdescr,
					pmode,pinputfiles,RsltsFile,MinLength,MaxLength,numbins,bindelta,
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

	PMode = (etNxxPMode)(pmode->count ? pmode->ival[0] : ePMdefault);
	if(PMode < ePMdefault || PMode > ePMKMerDist)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Processing mode '-m%d' specified outside of range %d..%d\n",PMode,ePMdefault,(int)ePMKMerDist);
		exit(1);
		}

	NumBins = 1000;
	BinDelta = 0;
	if(PMode == ePMdefault)
		{
		iMinLength = MinLength->count ? MinLength->ival[0] : cDfltMinLengthRange;
		if(iMinLength < 1 || iMinLength > cMaxN50Length)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Mininum contig length '-l%d' is not in range 1..%d",iMinLength,cMaxN50Length);
			exit(1);
			}

		iMaxLength = MaxLength->count ? MaxLength->ival[0] : cDfltN50Length;
		if(iMaxLength < iMinLength || iMaxLength > cMaxN50Length)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Maximum contig length '-L%d' is not in range %d..%d",iMaxLength,iMinLength,cMaxN50Length);
			exit(1);
			}
		}

	if(PMode == ePMKMerDist)
		{
		iMinLength = MinLength->count ? MinLength->ival[0] : cDfltMinLengthRange;
		if(iMinLength < 1 || iMinLength > cMaxTruncLength)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Mininum element length '-l%d' is not in range 1..%d",iMinLength,cMaxTruncLength);
			exit(1);
			}

		iMaxLength = MaxLength->count ? MaxLength->ival[0] : cDfltTruncLength;
		if(iMaxLength < iMinLength || iMaxLength > cMaxTruncLength)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Maximum element length '-L%d' is not in range %d..%d",iMaxLength,iMinLength,cMaxTruncLength);
			exit(1);
			}
		}


	NumInputFiles = 0;

	for(NumInputFiles=Idx=0;NumInputFiles < cMaxInFileSpecs && Idx < pinputfiles->count; Idx++)
		{
		pszInFastaFile[Idx] = NULL;
		if(pszInFastaFile[NumInputFiles] == NULL)
			pszInFastaFile[NumInputFiles] = new char [_MAX_PATH];
		strncpy(pszInFastaFile[NumInputFiles],pinputfiles->filename[Idx],_MAX_PATH);
		pszInFastaFile[NumInputFiles][_MAX_PATH-1] = '\0';
		CUtility::TrimQuotedWhitespcExtd(pszInFastaFile[NumInputFiles]);
		if(pszInFastaFile[NumInputFiles][0] != '\0')
			NumInputFiles++;
		}

	if(!NumInputFiles)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: After removal of whitespace, no input file(s) specified with '-i<filespec>' option)\n");
		exit(1);
		}

	if(PMode == ePMKMerDist)
		{
		if(!RsltsFile->count)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: No results file specified");
			exit(1);
			}
		strncpy(szRsltsFile,RsltsFile->filename[0],_MAX_PATH);
		szRsltsFile[_MAX_PATH-1] = '\0';
		}
	else
		{
		NumBins = numbins->count ? numbins->ival[0] : 1000;
		if(NumBins == 0)
			NumBins = 1000;
		if(NumBins < 10 || NumBins > 10000)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Number of bins '-b%d' is not in range 10..10000");
			exit(1);
			}
		BinDelta = bindelta->count ? bindelta->ival[0] : 0;
		switch(BinDelta) {
			case 0: case 1: case 2: case 5: case 10: case 25: case 50: case 100: case 250: case 500: case 1000:
				break;
			default:
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Bin length delta '-B%d' must be either 0 (auto),1,2,5,10,25,50,100,250,500 or 1000");
				exit(1);
			}

		szRsltsFile[0] = '\0';
		if(RsltsFile->count)
			{
			strncpy(szRsltsFile,RsltsFile->filename[0],_MAX_PATH);
			szRsltsFile[_MAX_PATH-1] = '\0';
			}
		}

// show user current resource limits
#ifndef _WIN32
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Resources: %s",CUtility::ReportResourceLimits());
#endif

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Version: %s Processing parameters:",cpszProgVer);

	for(Idx=0; Idx < NumInputFiles; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"input fasta sequence files (%d): '%s'",Idx+1,pszInFastaFile[Idx]);

	if(szRsltsFile[0] != '\0')
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Output to file: '%s'",szRsltsFile);

	if(PMode == ePMdefault)
		{
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Minimum contig sequence length: %d",iMinLength);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Maximum contig sequence length: %d",iMaxLength);
		if(szRsltsFile[0] != '\0')
			{
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Number of bins: %d",NumBins);
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Bin length delta: %d",BinDelta);
			}
		}
	if(PMode == ePMKMerDist)
		{
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Minimum fasta sequence length: %d",iMinLength);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"fasta sequences truncated length: %d",iMaxLength);
		}

	if(szExperimentName[0] != '\0')
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"This processing reference: %s",szExperimentName);

	if(gExperimentID > 0)
		{
		int ParamID;
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szLogFile),"log",szLogFile);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,sizeof(INT32),"mode",&PMode);
		
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,sizeof(INT32),"minlen",&iMinLength);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,sizeof(INT32),"maxlen",&iMaxLength);

		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,sizeof(INT32),"numbins",&numbins);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,sizeof(INT32),"bindelta",&bindelta);

		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,sizeof(INT32),"NumInputFiles",&NumInputFiles);
		for(Idx = 0; Idx < NumInputFiles; Idx++)
			ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(pszInFastaFile[Idx]),"in",pszInFastaFile[Idx]);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szRsltsFile),"out",szRsltsFile);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szSQLiteDatabase),"sumrslts",szSQLiteDatabase);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szExperimentName),"experimentname",szExperimentName);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szExperimentDescr),"experimentdescr",szExperimentDescr);
		}

#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	gStopWatch.Start();
	Rslt = Process(PMode,iMinLength,iMaxLength,NumBins,BinDelta,NumInputFiles,pszInFastaFile,szRsltsFile);
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

int 
Process(etNxxPMode Mode,			// processing mode - 0  N50 distributions, 1 n-mer count distributions
        int MinLength,				// core elements must be of at least this length
		int MaxLength,				// will be truncated to this length
		int NumBins,				// when generating length distributions then use this many bins - 0 defaults to using 1000
		int BinDelta,				// when generating length distributions then each bin holds this length delta - 0 defaults to auto determine from NunBins and longest sequence length
		int NumInputFiles,			// number of input sequence files
		char **pszInFastaFile,		// names of input sequence files (wildcards allowed)
		char *pszRsltsFile)			// file to write fasta into
{
int Rslt;
CFastaNxx *pFastaNxx;
if((pFastaNxx = new CFastaNxx) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create CFastaNxx object");
	return(eBSFerrObj);
	}
Rslt = pFastaNxx->Process(Mode,MinLength,MaxLength,NumBins,BinDelta,NumInputFiles,pszInFastaFile,pszRsltsFile);

delete pFastaNxx;
return(Rslt);
}

CFastaNxx::CFastaNxx()
{
m_pSeq = NULL;
m_pDimerCnts = NULL;
m_pTrimerCnts = NULL;
m_pTetramerCnts = NULL;
m_pBins = NULL;
m_pContigLengths = NULL;
Reset();
}


CFastaNxx::~CFastaNxx()
{
Reset();
}

void
CFastaNxx::Reset(void)
{
if(m_pDimerCnts != NULL)
	{
	delete(m_pDimerCnts);
	m_pDimerCnts = NULL;
	}
if(m_pTrimerCnts != NULL)
	{
	delete(m_pTrimerCnts);
	m_pTrimerCnts = NULL;
	}
if(m_pTetramerCnts != NULL)
	{
	delete(m_pTetramerCnts);
	m_pTetramerCnts = NULL;
	}
if(m_pSeq != NULL)
	{
	delete(m_pSeq);
	m_pSeq = NULL;
	}
if(m_pBins != NULL)
	{
	delete(m_pBins);
	m_pBins = NULL;
	}
if(m_pContigLengths != NULL)
	{
#ifdef _WIN32
	free(m_pContigLengths);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pContigLengths != MAP_FAILED)
		munmap(m_pContigLengths,m_AllocdCtgLens);
#endif	
	m_pContigLengths = NULL;
	}

m_NumContigLens = 0;
m_AllocdCtgLens = 0; 
m_TotLenCovered = 0;
m_MaxLengthRead = 0;

memset(m_BaseCnts,0,sizeof(m_BaseCnts));
memset(m_DistBaseCnts,0,sizeof(m_DistBaseCnts));
memset(m_SeqNsCnts,0,sizeof(m_SeqNsCnts));
}


int
CFastaNxx::ProcessFile(etNxxPMode Mode,	// processing mode - 0  N50 distributions
			int MinLength,				// core elements must be of at least this length
			int MaxLength,			    // process sequences of upto this length
			char *pszFastaFile)			// load sequences from this input file
{
int Rslt;
int CntIdx;
int SeqIdx;
int SeqLen;
bool bInSeq;
int NumAccepted;
int NumProcessed;
int PrevNumProcessed;
int NumUnderLen;
int NumOverLen;
int NumContigs;
int MaxLengthRead;
int SeqNsIdx;

etSeqBase *pSeqFwd;

CFasta *pFasta = NULL;

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Starting to process source fasta file '%s'",pszFastaFile);

if((pFasta = new CFasta())==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create CFasta object");
	return(eBSFerrObj);
	}

if((Rslt=pFasta->Open(pszFastaFile,true)) < eBSFSuccess)
	{
	while(pFasta->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,pFasta->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open fasta file '%s'",pszFastaFile);
	delete pFasta;
	return(Rslt);
	}


NumAccepted = 0;
NumProcessed = 0;
NumUnderLen = 0;
NumOverLen = 0;
NumContigs = 0;
MaxLengthRead = 0;
m_TotLenCovered = 0;
PrevNumProcessed = 0;
time_t Started = time(0);
bInSeq = false;
while((Rslt = SeqLen = pFasta->ReadSequence(m_pSeq,MaxLength + 10)) > eBSFSuccess)
	{
	if(SeqLen == eBSFFastaDescr)		// just read a descriptor line, slough, only interested in sequences
		{
		bInSeq = false;	
		continue;
		}
	if(bInSeq)							// only process initial 0..MaxLength of sequence
		continue;

	// every minute let the user know as to now many sequences have been processed

	// is it time to let user know progress
	if(NumProcessed > (PrevNumProcessed + 5000))
		{
		time_t Now = time(0);
		unsigned long ElapsedSecs = (unsigned long) (Now - Started);
		if(ElapsedSecs >= (60 * 10))
			{
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"Sequences loaded %d",NumProcessed);
			Started = Now;
			}
		PrevNumProcessed = NumProcessed;
		}
	

	NumProcessed += 1;
	if(Mode == ePMdefault)
		{
		bInSeq = true;

		// slough contig sequences not within requested length range 
		if(SeqLen < MinLength)
			{
			NumUnderLen += 1;
			continue;
			}
		if(SeqLen > MaxLength)			
			{
			NumOverLen += 1;
			continue;
			}

		if(SeqLen > MaxLengthRead)
			MaxLengthRead = SeqLen;

		NumContigs += 1;
		m_TotLenCovered += SeqLen;

		// is more memory required to hold this new contig length?
		if(((m_NumContigLens + 10) * sizeof(int)) > m_AllocdCtgLens)	// allow a little safety margin
			{
			int *pTmp;
			size_t AllocNeeded = m_AllocdCtgLens + (cAllocCtgLenDist * sizeof(int));
#ifdef _WIN32
			pTmp = (int *) realloc(m_pContigLengths,AllocNeeded);
#else
			pTmp = (int *)mremap(m_pContigLengths,m_AllocdCtgLens,AllocNeeded,MREMAP_MAYMOVE);
			if(pTmp == MAP_FAILED)
				pTmp = NULL;
#endif
			if(pTmp == NULL)
				{
				gDiagnostics.DiagOut(eDLInfo,gszProcName,"GenContig memory reallocation for contig lengths distribution of %d bytes failed..",AllocNeeded);
				return(eBSFerrMem);
				}
			m_pContigLengths = pTmp;
			m_AllocdCtgLens = AllocNeeded;
			}
		m_pContigLengths[m_NumContigLens++] = SeqLen;
		NumAccepted += 1;
		continue;
		}
	else
		{
		if(SeqLen < MinLength)
			{
			NumUnderLen += 1;
			continue;
			}

		if(SeqLen > MaxLength)			// overlength sequences are truncated
			{
			bInSeq = true;
			NumOverLen += 1;
			SeqLen = MaxLength;
			}
		m_TotLenCovered += SeqLen;
		}

	if(SeqLen > MaxLengthRead)
		MaxLengthRead = SeqLen;

	CSeqTrans::RemoveMasking(m_pSeq,SeqLen); // ensure nucleotides are not masked

	//GenStats(SeqLen,pSeq,&CompStats);
	SeqNsIdx = 0;
	CntIdx  = 0;
	pSeqFwd = m_pSeq; 
	for(SeqIdx = 0; SeqIdx < SeqLen; SeqIdx++,pSeqFwd++)
		{
		if(*pSeqFwd == eBaseN)
			SeqNsIdx += 1;

		m_BaseCnts[*pSeqFwd] += 1;
		m_DistBaseCnts[*pSeqFwd][SeqIdx] += 1;

		if(SeqIdx < SeqLen-1)
			{
			CntIdx = (((*pSeqFwd * 5) + pSeqFwd[1]) * cMaxFastQSeqLen)+SeqIdx;
			m_pDimerCnts[CntIdx] += 1;
			}

		if(SeqIdx < SeqLen-2)
			{
			CntIdx = (((((*pSeqFwd * 5) + pSeqFwd[1])*5) + pSeqFwd[2]) * cMaxFastQSeqLen)+SeqIdx;
			m_pTrimerCnts[CntIdx] += 1;
			}

		if(SeqIdx < SeqLen-3)
			{
			CntIdx = (((((((*pSeqFwd * 5) + pSeqFwd[1])*5) + pSeqFwd[2]) * 5) +  pSeqFwd[3]) * cMaxFastQSeqLen)+SeqIdx;
			m_pTetramerCnts[CntIdx] += 1;
			}
		}
	m_SeqNsCnts[SeqNsIdx] += 1;
	NumAccepted += 1;
	}
if(Rslt < eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Encountered errors after %d sequences loaded of which %d were processed. There were %d underlength and %d overlength",
					NumProcessed,NumAccepted,NumUnderLen,NumOverLen);
	while(pFasta->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,pFasta->GetErrMsg());
	delete pFasta;
	return(Rslt);
	}
delete pFasta;

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Sequences loaded %d of which %d were processed. There were %d underlength and %d overlength",
					NumProcessed,NumAccepted,NumUnderLen,NumOverLen);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Maximum length sequence loaded and accepted for processing was %d bp",MaxLengthRead);
if(NumAccepted > 0)
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Total length of all accepted sequences was %lld bp with mean length %d",m_TotLenCovered,m_TotLenCovered/NumAccepted);
if(MaxLengthRead > m_MaxLengthRead)
	m_MaxLengthRead = MaxLengthRead;
return(Rslt);
}

int 
CFastaNxx::Process(etNxxPMode Mode,		// processing mode - 0  N50 distributions, 1 n-mer count distributions
        int MinLength,				// core elements must be of at least this length
		int MaxLength,				// will be truncated to this length
		int NumBins,				// when generating length distributions then use this many bins - 0 defaults to using 1000
		int BinDelta,				// when generating length distributions then each bin holds this length delta - 0 defaults to auto determine from NunBins and longest sequence length
		int NumInputFiles,			// number of input sequence files
		char **pszInFastaFile,		// names of input sequence files (wildcards allowed)
		char *pszRsltsFile)			// file to write fasta into
{
int Rslt;
int WrtOfs;
char szWrtBuff[cMaxReadLen+1];
int hRslts;
int SeqIdx;

int CntIdx;

Reset();

if(Mode == ePMKMerDist)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Allocating memory for N-Mer count distributions...");
	if((m_pDimerCnts = new UINT32 [5*5*cMaxFastQSeqLen])==NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate memory for dimer counts",strerror(errno));
		return(eBSFerrMem);
		}
	
	if((m_pTrimerCnts = new UINT32 [5*5*5*cMaxFastQSeqLen])==NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate memory for trimer counts",strerror(errno));
		Reset();
		return(eBSFerrMem);
		}

	if((m_pTetramerCnts = new UINT32 [5*5*5*5*cMaxFastQSeqLen])==NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate memory for tetramer counts",strerror(errno));
		Reset();
		return(eBSFerrMem);
		}

	memset(m_BaseCnts,0,sizeof(m_BaseCnts));
	memset(m_pDimerCnts,0,sizeof(UINT32) * 5*5*cMaxFastQSeqLen);
	memset(m_pTrimerCnts,0,sizeof(UINT32) * 5*5*5*cMaxFastQSeqLen);
	memset(m_pTetramerCnts,0,sizeof(UINT32) * 5*5*5*5*cMaxFastQSeqLen);
	memset(m_SeqNsCnts,0,sizeof(m_SeqNsCnts));
	memset(m_DistBaseCnts,0,sizeof(m_DistBaseCnts));
	}

if((m_pSeq = new UINT8 [MaxLength+10])==NULL)	// allows to check if sequence has been loaded which is longer than MaxLength
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate memory to hold sequence of length %d",MaxLength);
	Reset();
	return(eBSFerrMem);
	}

if(Mode == ePMdefault)
	{
	size_t MemReq = cAllocCtgLenDist * sizeof(int);
#ifdef _WIN32
	m_pContigLengths = (int *) malloc(MemReq);	// initial and perhaps the only allocation
	if(m_pContigLengths == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"GenContigs: Memory allocation of %d bytes for contig length distributions - %s",MemReq,strerror(errno));
		Reset();
		return(eBSFerrMem);
		}
#else
	// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
	m_pContigLengths = (int *)mmap(NULL,MemReq, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
	if(m_pContigLengths == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadReads: Memory allocation of %d bytes through mmap()  failed - %s",MemReq,strerror(errno));
		m_pContigLengths = NULL;
		Reset();
		return(eBSFerrMem);
		}
#endif
	m_AllocdCtgLens = MemReq;
	m_NumContigLens = 0;
	m_TotLenCovered = 0;
	}

CSimpleGlob glob(SG_GLOB_FULLSORT);
int Idx;
char *pszInFile;
for(Idx = 0; Idx < NumInputFiles; Idx++)
	{
	glob.Init();
	if(glob.Add(pszInFastaFile[Idx]) < SG_SUCCESS)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to glob '%s",pszInFastaFile[Idx]);
		Reset();
		return(eBSFerrOpnFile);	// treat as though unable to open file
		}	
	if(glob.FileCount() <= 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to locate any source raw reads file matching '%s",pszInFastaFile[Idx]);
		continue;
		}
	Rslt = eBSFSuccess;
	for (int FileID = 0; Rslt >= eBSFSuccess &&  FileID < glob.FileCount(); ++FileID)
		{
		pszInFile = glob.File(FileID);
		if((Rslt = ProcessFile(Mode,MinLength,MaxLength,pszInFile)) < 0)
			{
			delete m_pSeq;
			m_pSeq = NULL;
			return(Rslt);
			}
		}

	}


if(Mode == ePMdefault)
	{
	INT64 NxLens[11];	// to hold contig lengths for all Nx where Nx varies from 10 to 100 in increments of N10
	INT64 SumCtgLens;
	INT64 NxSum;
	INT64 *pNxLens;
	int CtgLenIdx;
	int CtgNxIdx;
	int NumCtgs;
	int *pCnts;

	if(pszRsltsFile != NULL && pszRsltsFile[0] != '\0')
		{
		#ifdef _WIN32
		if((hRslts = open(pszRsltsFile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
		#else
		if((hRslts = open(pszRsltsFile, O_RDWR | O_CREAT | O_TRUNC, S_IREAD | S_IWRITE ))==-1)
		#endif
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create or truncate results file %s error: %s",pszRsltsFile,strerror(errno));
			Reset();
			return(eBSFerrCreateFile);
			}

		if(NumBins == 0)
			NumBins = 1000;			// default is 1000 bins for length distributions
		if(BinDelta == 0)			// needing to auto determine bin length delta's of 1,2,5,10,20,25,50,100,200,500 and 1000 dependent on the MaxLengthRead accepted
			{
			if(m_MaxLengthRead <= NumBins)
				BinDelta = 1;
			else
				{
				if(m_MaxLengthRead < 2 * NumBins)
					BinDelta = 2;
				else
					if(m_MaxLengthRead < 5 * NumBins)
						BinDelta = 5;
					else
						if(m_MaxLengthRead < 10 * NumBins)
							BinDelta = 10;
						else
							if(m_MaxLengthRead < 25 * NumBins)
								BinDelta = 25;
							else
								if(m_MaxLengthRead < 50 * NumBins)
									BinDelta = 50;
								else
									if(m_MaxLengthRead < 100 * NumBins)
										BinDelta = 100;
									else
										if(m_MaxLengthRead < 250 * NumBins)
											BinDelta = 250;
										else
											if(m_MaxLengthRead < 500 * NumBins)
												BinDelta = 500;
											else
												BinDelta = 1000;
				}
			}

		if((m_pBins = new UINT32 [NumBins+1])==NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"GenContigs: Memory allocation for %d bins failed - %s",NumBins,strerror(errno));
			Reset();
			return(eBSFerrMem);
			}
		memset(m_pBins,0,sizeof(UINT32) * (NumBins+1));
		}
	else
		m_pBins = NULL;

	SumCtgLens = 0;
	pCnts = m_pContigLengths;
	for(CtgLenIdx = 0; CtgLenIdx < m_NumContigLens; CtgLenIdx++,pCnts++)
		{
		SumCtgLens += *pCnts;
		if(m_pBins != NULL && *pCnts != 0)
			m_pBins[min(*pCnts / BinDelta,NumBins)] += 1;
		}	

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"  Total sequence length over accepted contigs: %lld",SumCtgLens);

	if(m_pBins != NULL)
		{
		WrtOfs = 0;
		pCnts = (int *)m_pBins;
		for(CtgLenIdx = 0; CtgLenIdx <= NumBins; CtgLenIdx++, pCnts++)
			{
			WrtOfs += sprintf(&szWrtBuff[WrtOfs],"%d,%d\n",CtgLenIdx * BinDelta,*pCnts);
			if(WrtOfs + 50 > sizeof(szWrtBuff))
				{
				CUtility::SafeWrite(hRslts,szWrtBuff,WrtOfs);
				WrtOfs = 0;
				}
			}
		if(WrtOfs > 0)
			{
			CUtility::SafeWrite(hRslts,szWrtBuff,WrtOfs);
			WrtOfs = 0;
			}
#ifdef _WIN32
		_commit(hRslts);
#else
		fsync(hRslts);
#endif
		hRslts = -1;
		delete m_pBins;
		m_pBins = NULL;
		}

	m_MTqsort.SetMaxThreads(4);
	m_MTqsort.qsort(m_pContigLengths,(UINT64)m_NumContigLens,sizeof(int),SortByContigLen);


	pNxLens = NxLens;
	for(CtgNxIdx = 1; CtgNxIdx <= 10; CtgNxIdx++,pNxLens++)
		*pNxLens = (SumCtgLens * CtgNxIdx)/10;

	NxSum = 0;
	CtgNxIdx = 0;
	NumCtgs = 0;
	pCnts = m_pContigLengths;
	for(int Idx = 0; Idx < m_NumContigLens; Idx++,pCnts++)
		{
		if(*pCnts)
			{
			NxSum += *pCnts;
			NumCtgs += 1;
			}

		if(NxSum >= NxLens[CtgNxIdx])
			{
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"  N%d contig length is %d, Contigs: %d (%1.2f%%)",(CtgNxIdx+1)*10,*pCnts,NumCtgs, ((double)NumCtgs * 100)/(double)m_NumContigLens);
			NxLens[CtgNxIdx++] = *pCnts;
			}
		}

#ifdef _WIN32
	free(m_pContigLengths);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pContigLengths != MAP_FAILED)
		munmap(m_pContigLengths,m_AllocdCtgLens);
#endif	
	m_pContigLengths = NULL;

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing completed");
	delete m_pSeq;
	m_pSeq = NULL;
	return(eBSFSuccess);
	}



// output results
gDiagnostics.DiagOut(eDLInfo,gszProcName,"About to write results into file '%s'",pszRsltsFile);

#ifdef _WIN32
if((hRslts = open(pszRsltsFile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
if((hRslts = open(pszRsltsFile, O_RDWR | O_CREAT | O_TRUNC, S_IREAD | S_IWRITE ))==-1)
#endif
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create or truncate results file %s error: %s",pszRsltsFile,strerror(errno));
	Reset();
	return(eBSFerrCreateFile);
	}

int ReadIdx;
WrtOfs = 0;

// momomer distribution as counts
// row 1 is the base position (1..MaxLengthRead) heading
for(SeqIdx = 1; SeqIdx <= m_MaxLengthRead; SeqIdx++)
	{
	WrtOfs += sprintf(&szWrtBuff[WrtOfs],",%d",SeqIdx);
	if(WrtOfs + 100 >= sizeof(szWrtBuff))
		{
		CUtility::SafeWrite(hRslts,szWrtBuff,WrtOfs);
		WrtOfs = 0;
		}
	}
for(ReadIdx = 0; ReadIdx < 5; ReadIdx++)
	{
	WrtOfs += sprintf(&szWrtBuff[WrtOfs],"\n%c",CSeqTrans::MapBase2Ascii(ReadIdx));
	if(WrtOfs + 100 >= sizeof(szWrtBuff))
		{
		CUtility::SafeWrite(hRslts,szWrtBuff,WrtOfs);
		WrtOfs = 0;
		}
	for(SeqIdx = 0; SeqIdx < m_MaxLengthRead; SeqIdx++)
		{
		WrtOfs += sprintf(&szWrtBuff[WrtOfs],",%d",m_DistBaseCnts[ReadIdx][SeqIdx]);
		if(WrtOfs + 100 >= sizeof(szWrtBuff))
			{
			CUtility::SafeWrite(hRslts,szWrtBuff,WrtOfs);
			WrtOfs = 0;
			}
		}
	}
WrtOfs += sprintf(&szWrtBuff[WrtOfs],"\n\n");
CUtility::SafeWrite(hRslts,szWrtBuff,WrtOfs);
WrtOfs = 0;

// repeat monomer distribution but this time as proportions
for(SeqIdx = 1; SeqIdx <= m_MaxLengthRead; SeqIdx++)
	{
	WrtOfs += sprintf(&szWrtBuff[WrtOfs],",%d",SeqIdx);
	if(WrtOfs + 100 >= sizeof(szWrtBuff))
		{
		CUtility::SafeWrite(hRslts,szWrtBuff,WrtOfs);
		WrtOfs = 0;
		}
	}
for(ReadIdx = 0; ReadIdx < 4; ReadIdx++)
	{
	WrtOfs += sprintf(&szWrtBuff[WrtOfs],"\n%c",CSeqTrans::MapBase2Ascii(ReadIdx));
	if(WrtOfs + 100 >= sizeof(szWrtBuff))
		{
		CUtility::SafeWrite(hRslts,szWrtBuff,WrtOfs);
		WrtOfs = 0;
		}
	for(SeqIdx = 0; SeqIdx < m_MaxLengthRead; SeqIdx++)
		{
		int ColTot = 0;
		for(int RowIdx = 0; RowIdx < 4; RowIdx++)
			ColTot += m_DistBaseCnts[RowIdx][SeqIdx];
		if(ColTot > 0)
			WrtOfs += sprintf(&szWrtBuff[WrtOfs],",%1.4f",(double)m_DistBaseCnts[ReadIdx][SeqIdx]/(double)ColTot);
		else
			WrtOfs += sprintf(&szWrtBuff[WrtOfs],"0.0000");
		if(WrtOfs + 100 >= sizeof(szWrtBuff))
			{
			CUtility::SafeWrite(hRslts,szWrtBuff,WrtOfs);
			WrtOfs = 0;
			}
		}
	}
WrtOfs += sprintf(&szWrtBuff[WrtOfs],"\n\n");
CUtility::SafeWrite(hRslts,szWrtBuff,WrtOfs);
WrtOfs = 0;


// dimer distributions
// row 1 is the base position (1..MaxLengthRead) heading
for(SeqIdx = 1; SeqIdx <= m_MaxLengthRead-1; SeqIdx++)
	{
	WrtOfs += sprintf(&szWrtBuff[WrtOfs],",%d",SeqIdx);
	if(WrtOfs + 100 >= sizeof(szWrtBuff))
		{
		CUtility::SafeWrite(hRslts,szWrtBuff,WrtOfs);
		WrtOfs = 0;
		}
	}

int Base1Idx;
for(Base1Idx = 0; Base1Idx < 25; Base1Idx++)
	{
	WrtOfs += sprintf(&szWrtBuff[WrtOfs],"\n%c",CSeqTrans::MapBase2Ascii(Base1Idx / 5));
	WrtOfs += sprintf(&szWrtBuff[WrtOfs],"%c",CSeqTrans::MapBase2Ascii(Base1Idx % 5));
	CntIdx = Base1Idx * cMaxFastQSeqLen;
	for(SeqIdx = 0; SeqIdx < m_MaxLengthRead-1; SeqIdx++)
		{
		WrtOfs += sprintf(&szWrtBuff[WrtOfs],",%d",m_pDimerCnts[CntIdx + SeqIdx]);
		if(WrtOfs + 100 >= sizeof(szWrtBuff))
			{
			CUtility::SafeWrite(hRslts,szWrtBuff,WrtOfs);
			WrtOfs = 0;
			}
		}
	}
WrtOfs += sprintf(&szWrtBuff[WrtOfs],"\n\n");
CUtility::SafeWrite(hRslts,szWrtBuff,WrtOfs);
WrtOfs = 0;

// trimer distributions
// row 1 is the base position (1..MaxLengthRead) heading
for(SeqIdx = 1; SeqIdx <= m_MaxLengthRead-2; SeqIdx++)
	{
	WrtOfs += sprintf(&szWrtBuff[WrtOfs],",%d",SeqIdx);
	if(WrtOfs + 100 >= sizeof(szWrtBuff))
		{
		CUtility::SafeWrite(hRslts,szWrtBuff,WrtOfs);
		WrtOfs = 0;
		}
	}


for(Base1Idx = 0; Base1Idx < (5*5*5); Base1Idx++)
	{
	WrtOfs += sprintf(&szWrtBuff[WrtOfs],"\n%c",CSeqTrans::MapBase2Ascii((Base1Idx/5)/5));
	WrtOfs += sprintf(&szWrtBuff[WrtOfs],"%c",CSeqTrans::MapBase2Ascii((Base1Idx/5)%5));
	WrtOfs += sprintf(&szWrtBuff[WrtOfs],"%c",CSeqTrans::MapBase2Ascii(Base1Idx%5));
	CntIdx = Base1Idx * cMaxFastQSeqLen;
	if(WrtOfs + 100 >= sizeof(szWrtBuff))
		{
		CUtility::SafeWrite(hRslts,szWrtBuff,WrtOfs);
		WrtOfs = 0;
		}
	for(SeqIdx = 0; SeqIdx < m_MaxLengthRead-2; SeqIdx++)
		{
		WrtOfs += sprintf(&szWrtBuff[WrtOfs],",%d",m_pTrimerCnts[CntIdx+SeqIdx]);
		if(WrtOfs + 100 >= sizeof(szWrtBuff))
			{
			CUtility::SafeWrite(hRslts,szWrtBuff,WrtOfs);
			WrtOfs = 0;
			}
		}
	}
WrtOfs += sprintf(&szWrtBuff[WrtOfs],"\n\n");
CUtility::SafeWrite(hRslts,szWrtBuff,WrtOfs);
WrtOfs = 0;

// tetramer distributions
// row 1 is the base position (1..MaxLengthRead) heading
for(SeqIdx = 1; SeqIdx <= m_MaxLengthRead-3; SeqIdx++)
	{
	WrtOfs += sprintf(&szWrtBuff[WrtOfs],",%d",SeqIdx);
	if(WrtOfs + 100 >= sizeof(szWrtBuff))
		{
		CUtility::SafeWrite(hRslts,szWrtBuff,WrtOfs);
		WrtOfs = 0;
		}
	}


for(Base1Idx = 0; Base1Idx < (5*5*5*5); Base1Idx++)
	{
	WrtOfs += sprintf(&szWrtBuff[WrtOfs],"\n%c",CSeqTrans::MapBase2Ascii(((Base1Idx/5)/5)/5));
	WrtOfs += sprintf(&szWrtBuff[WrtOfs],"%c",CSeqTrans::MapBase2Ascii(((Base1Idx/5)/5)%5));
	WrtOfs += sprintf(&szWrtBuff[WrtOfs],"%c",CSeqTrans::MapBase2Ascii((Base1Idx/5)%5));
	WrtOfs += sprintf(&szWrtBuff[WrtOfs],"%c",CSeqTrans::MapBase2Ascii(Base1Idx%5));
	CntIdx = Base1Idx * cMaxFastQSeqLen;

	if(WrtOfs + 100 >= sizeof(szWrtBuff))
		{
		CUtility::SafeWrite(hRslts,szWrtBuff,WrtOfs);
		WrtOfs = 0;
		}
	for(SeqIdx = 0; SeqIdx < m_MaxLengthRead-3; SeqIdx++)
		{
		WrtOfs += sprintf(&szWrtBuff[WrtOfs],",%d",m_pTetramerCnts[CntIdx+SeqIdx]);
		if(WrtOfs + 100 >= sizeof(szWrtBuff))
			{
			CUtility::SafeWrite(hRslts,szWrtBuff,WrtOfs);
			WrtOfs = 0;
			}
		}
	}
WrtOfs += sprintf(&szWrtBuff[WrtOfs],"\n");
CUtility::SafeWrite(hRslts,szWrtBuff,WrtOfs);
close(hRslts);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing completed");
delete m_pSeq;
m_pSeq = NULL;
return(eBSFSuccess);
}

// SortByContigLen
// Sorts contig lengths descending
int
CFastaNxx::SortByContigLen(const void *arg1, const void *arg2)
{
int *pEl1 = (int *)arg1;
int *pEl2 = (int *)arg2;

if(*pEl1 < *pEl2)
	return(1);
if(*pEl1 > *pEl2)
	return(-1);
return(0);
}



