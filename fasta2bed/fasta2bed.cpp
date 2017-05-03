// fasta2bed.cpp : Defines the entry point for the console application.
//
// 0.0.2 Cosmetic changes only

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

const char *cpszProgVer = "1.0.0";		// increment with each release
const int cMaxInFileSpecs = 50;			// allow at most this many input file specs

// processing modes
typedef enum TAG_ePMode {
	ePMdefault,					// default BED , no exons
	ePMplaceholder				// used to set the enumeration range
	} etPMode;

int 
Process(etPMode Mode,				// processing mode - 0 default BED , single exon
		int NumInputFiles,			// number of input sequence files
		char **pszInFastaFile,		// names of input sequence files (wildcards allowed)
		char *pszRsltsFile);		// file to write BED format into

CStopWatch gStopWatch;
CDiagnostics gDiagnostics;			// for writing diagnostics messages to log file
char gszProcName[_MAX_FNAME];		// process name


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

int Rslt;
etPMode PMode;				// processing mode
int Idx;

int NumInputFiles;							// number of input sequence files
char *pszInFastaFile[cMaxInFileSpecs];		// names of input sequence files (wildcards allowed)

char szRsltsFile[_MAX_PATH];	// output stats to this file


// command line args
struct arg_lit  *help    = arg_lit0("h","help",                 "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to screen and logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"diagnostics log file");

struct arg_int *pmode = arg_int0("m","mode","<int>",		    "Processing mode:  0 - BED output, single exon");
struct arg_file *pinputfiles = arg_filen("i","infasta","<file>",1,cMaxInFileSpecs,"input fasta file(s)");
struct arg_file *RsltsFile = arg_file0("o","output","<file>",	"output BED file");
struct arg_end *end = arg_end(20);

void *argtable[] = {help,version,FileLogLevel,LogFile,
					pmode,pinputfiles,RsltsFile,
					end};
char **pAllArgs;
int argerrors;
argerrors = CUtility::arg_parsefromfile(argc,(char **)argv,&pAllArgs);
if(argerrors >= 0)
	argerrors = arg_parse(argerrors,pAllArgs,argtable);

/* special case: '--help' takes precedence over error reporting */
if (help->count > 0)
        {
		printf("\n%s Generate BED file from multifasta file, Version %s\nOptions ---\n", gszProcName,cpszProgVer);
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

	if(!RsltsFile->count)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: No output BED file specified");
		exit(1);
		}
	strncpy(szRsltsFile,RsltsFile->filename[0],_MAX_PATH);
	szRsltsFile[_MAX_PATH-1] = '\0';

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing parameters:");
	for(Idx=0; Idx < NumInputFiles; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"input fasta sequence files (%d): '%s'",Idx+1,pszInFastaFile[Idx]);

	if(szRsltsFile[0] != '\0')
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Output to BED file: '%s'",szRsltsFile);

#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	// processing here...
	gStopWatch.Start();
	Rslt = Process(PMode,NumInputFiles,pszInFastaFile,szRsltsFile);

	gStopWatch.Stop();
	Rslt = Rslt >=0 ? 0 : 1;
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Exit code: %d Total processing time: %s",Rslt,gStopWatch.Read());
	exit(Rslt);
	}
else
	{
	printf("\n%s Generate BED file from multifasta file, Version %s\n",gszProcName,cpszProgVer);
	arg_print_errors(stdout,end,gszProcName);
	arg_print_syntax(stdout,argtable,"\nUse '-h' to view option and parameter usage\n");
	exit(1);
	}
return 0;
}

// globals
etSeqBase *m_pSeq;

int
ProcessFile(etPMode Mode,				// processing mode - 0  BED output, no exon features
			int hOutputBED,				// write to this BED file
			char *pszFastaFile)			// input file containing fasta sequences
{
int Rslt;
int SeqLen;
bool bInSeq;
int NumAccepted;
int NumProcessed;

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
char szFeature[128];
char szBEDFeature[256];
int LineLen;
int DescrLen;
bInSeq = false;
while((Rslt = SeqLen = pFasta->ReadSequence(NULL,0)) > eBSFSuccess)
	{
	if(SeqLen == eBSFFastaDescr)		// just read a descriptor line, parse out the feature identifier
		{
		DescrLen = pFasta->ReadDescriptor(szBEDFeature, sizeof(szBEDFeature));
		sscanf(szBEDFeature," %s[ ,]",szFeature);
		szFeature[35] = '\0';
		bInSeq = false;	
		continue;
		}
	LineLen = sprintf(szBEDFeature,"%s\t0\t%d\t%s\t0\t+\t0\t%d\t0\t1\t%d,\t0\n",szFeature,SeqLen,szFeature,SeqLen,SeqLen);
	Rslt=write(hOutputBED,szBEDFeature,LineLen);
	NumProcessed += 1;

	bInSeq = true;

	NumAccepted += 1;
	continue;
	}

delete pFasta;
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Sequences processed %d ",	NumProcessed);
return(Rslt);
}

void
Reset(void)
{
}

void
Init(void)
{

}

int 
Process(etPMode Mode,				// processing mode - 0  BED file output
		int NumInputFiles,			// number of input sequence files
		char **pszInFastaFile,		// names of input sequence files (wildcards allowed)
		char *pszRsltsFile)			// file to write BED format into
{
int Rslt;
int hRslts;

Init();

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
		if((Rslt = ProcessFile(Mode,hRslts,pszInFile)) < 0)
			{
			close(hRslts);
			return(Rslt);
			}
		}

	}

close(hRslts);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing completed");
return(eBSFSuccess);
}




