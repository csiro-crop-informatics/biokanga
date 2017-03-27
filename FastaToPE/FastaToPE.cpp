// FastaToPE.cpp : Defines the entry point for the console application.
// Loads multifasta file in which each entry is actually paired end reads concatenated together
// These concatenated entries are split into 2 separate multifasta files P1 and P2
// Option added whereby the P2 read can be RevCpl, complemented, or reversed
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

const char *cpszProgVer = "1.0.0";			// increment with each release

typedef enum TAG_ePMode {
	ePMdefault = 0,	// retain split reads in original orientation
	ePMReverseR2,	// reverse R2
	ePMComplementR2, // complement R2
	ePMRevCplR2,	// reverse complement R2
	ePMplaceholder	// used to determine the number of processing modes supported
	} etPMode;

const int cSeqBuffLen = 10000;		// buffer sequences in blocks of this length

int Process(int Mode,				    // processing mode
			char *pszInFile,			// input multifasta to split
				char *pszOutputFile);	// write splits to files with this name prefix

CStopWatch gStopWatch;
CDiagnostics gDiagnostics;				// for writing diagnostics messages to log file
char gszProcName[_MAX_FNAME];			// process name

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
	return _T("kangade");
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

int Rslt;
int iMode = 0;			// processing mode

char szInFile[_MAX_PATH];			// input fasta containing reads to split

char szOutputFile[_MAX_PATH];

// command line args
struct arg_lit  *help    = arg_lit0("hH","help",                "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to screen and logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"diagnostics log file");

struct arg_int  *Mode = arg_int0("m","mode","<int>",			"processing mode - 0: retain split reads in original orientation, 1: reverse R2, 2: complement R2, 3: RevCpl R2 (default 0)");
struct arg_file *InFile = arg_file1("i",NULL,"<file>",			"input file containing concatenated paired end reads to split");
struct arg_file *OutFile= arg_file1("o",NULL,"<file>",			"output split reads to two files using this as the base file name");
struct arg_end *end = arg_end(20);

void *argtable[] = {help,version,FileLogLevel,LogFile,Mode,InFile,OutFile,end};

char **pAllArgs;
int argerrors;
argerrors = CUtility::arg_parsefromfile(argc,(char **)argv,&pAllArgs);
if(argerrors >= 0)
	argerrors = arg_parse(argerrors,pAllArgs,argtable);

/* special case: '--help' takes precedence over error reporting */
if (help->count > 0)
        {
		printf("\n%s Split concatenated reads, Version %s\nOptions ---\n", gszProcName,cpszProgVer);
        arg_print_syntax(stdout,argtable,"\n");
        arg_print_glossary(stdout,argtable,"  %-25s %s\n");
		printf("\nNote: Parameters can be entered into a parameter file, one parameter per line.");
		printf("\n      To invoke this parameter file then precede it's name with '@'");
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

	iMode = Mode->count ? Mode->ival[0] : 0;
	if(iMode < 0 || iMode >= ePMplaceholder)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Unsupported Mode '-m%d' requested, must be in range 0..%d",iMode,ePMplaceholder-1);
		exit(1);
		}

	strncpy(szInFile,InFile->filename[0],_MAX_PATH);
	szInFile[_MAX_PATH] = '\0';

	strncpy(szOutputFile,OutFile->filename[0],_MAX_PATH);
	szOutputFile[_MAX_PATH] = '\0';

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing parameters:");

	char *pszMode;
	switch(iMode) {
		case ePMdefault:		// retain split reads in original orientation
			pszMode = (char *)"retain split reads in original orientation";
			break;
		case ePMReverseR2:		// reverse R2
			pszMode = (char *)"reverse R2";
			break;
		case ePMComplementR2:	// complement R2
			pszMode = (char *)"complement R2";
			break;
		case ePMRevCplR2:		// reverse complement R2
			pszMode = (char *)"reverse complement R2";
			break;
		}

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Mode: %s",pszMode);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Read concatenated paired reads in this file: '%s'",szInFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Output split paired reads to two files with prefixed file names: '%s'",szOutputFile);
	
	gStopWatch.Start();
#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	Rslt = Process(iMode,				    // processing mode
					szInFile,			// input multifasta to split
				szOutputFile);	// write splits to files with this name prefix
	Rslt = Rslt >=0 ? 0 : 1;
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Exit code: %d Total processing time: %s",Rslt,gStopWatch.Read());
	exit(Rslt);
	}
else
	{
	printf("\n%s Split concatenated reads, Version %s\n",gszProcName,cpszProgVer);
	arg_print_errors(stdout,end,gszProcName);
	arg_print_syntax(stdout,argtable,"\nUse '-h' to view option and parameter usage\n");
	exit(1);
	}
return 0;
}

CFasta *m_pInFasta;
CFasta *m_pOutP1;
CFasta *m_pOutP2;


void
Reset(void)
{
if(m_pInFasta != NULL)
	{
	delete m_pInFasta;
	m_pInFasta = NULL;
	}
if(m_pOutP1 != NULL)
	{
	m_pOutP1->Close();
	delete m_pOutP1;
	m_pOutP1 = NULL;
	}
if(m_pOutP2 != NULL)
	{
	m_pOutP2->Close();
	delete m_pOutP2;
	m_pOutP2 = NULL;
	}
}

void
Init(void)
{
m_pInFasta = NULL;
m_pOutP1 = NULL;
m_pOutP2 = NULL;
Reset();
}


const int cMaxReadSeqLen = cMaxReadLen;		// allow for paired reads of this length

int Process(int Mode,				    // processing mode
			char *pszInFile,			// input multifasta to split
				char *pszOutputFile)	// write splits to files with this name prefix
{
int Rslt;
int NumConcatReads;
int ConcatLen;
int MaxConcatLen;
int ReadLen;
char szDescr[120];
char szP1File[_MAX_PATH];
char szP2File[_MAX_PATH];
UINT8 ConcatSeq[(cMaxReadSeqLen + 100) * 2];
UINT8 ReadSeq[cMaxReadSeqLen + 1];

Init();

MaxConcatLen = cMaxReadSeqLen * 2;
// open input file

if((m_pInFasta = new CFasta()) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CFasta object");
	Reset();
	return(eBSFerrObj);
	}

if((Rslt=m_pInFasta->Open(pszInFile,true)) < eBSFSuccess)
	{
	while(m_pInFasta->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,m_pInFasta->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open fasta file '%s'",pszInFile);
	Reset();
	return(Rslt);
	}

// create output files

strcpy(szP1File,pszOutputFile);
strcat(szP1File,".R1.fasta");

strcpy(szP2File,pszOutputFile);
strcat(szP2File,".R2.fasta");

m_pOutP1 = new CFasta();
if(m_pOutP1 == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CFasta object");
		Reset();
		return(eBSFerrObj);
		}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Creating output file to contain PE 1 reads: '%s'",szP1File);
if((Rslt = m_pOutP1->Open(szP1File,false)) != eBSFSuccess)
	{
	while(m_pOutP1->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,m_pOutP1->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create fasta file '%s'",szP1File);
	Reset();
	return(Rslt);
	}

m_pOutP2 = new CFasta();
if(m_pOutP2 == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CFasta object");
		Reset();
		return(eBSFerrObj);
		}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Creating output file to contain PE 2 reads: '%s'",szP2File);
if((Rslt = m_pOutP2->Open(szP2File,false)) != eBSFSuccess)
	{
	while(m_pOutP2->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,m_pOutP2->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create fasta file '%s'",szP2File);
	Reset();
	return(Rslt);
	}

NumConcatReads = 0;
while((Rslt = ConcatLen = m_pInFasta->ReadSequence(ConcatSeq,sizeof(ConcatSeq)-1)) > eBSFSuccess)	
	{
	if(Rslt == eBSFFastaDescr)		// just read a descriptor line
		{
		m_pInFasta->ReadDescriptor(szDescr,sizeof(szDescr)-1);
		m_pOutP1->WriteDescriptor(szDescr);		
		m_pOutP2->WriteDescriptor(szDescr);
		NumConcatReads += 1;
		if(!(NumConcatReads % 5000000) || NumConcatReads == 1)
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing concatenated read %d",NumConcatReads);
		continue;
		}

	if(ConcatLen > MaxConcatLen)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Concatenated sequence read %d, of length %d, is longer than the max accepted length of %d",NumConcatReads,ConcatLen,MaxConcatLen);
		Reset();
		return(-1);
		}

	if(ConcatLen & 0x001)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Concatenated sequence read %d, length %d, can not be split into two equally sized reads",NumConcatReads,ConcatLen);
		Reset();
		return(-1);
		}

	ReadLen = ConcatLen/2;
	UINT8 *pInBase = ConcatSeq;
	UINT8 *pOutBase = ReadSeq;
	int BaseIdx;
	for(BaseIdx = 0; BaseIdx < ReadLen; BaseIdx++,pInBase++,pOutBase++)
		{
		switch(*pInBase & 0x07) {
			case 0:
				*pOutBase = 'a';
				break;
			case 1:
				*pOutBase = 'c';
				break;
			case 2:
				*pOutBase = 'g';
				break;
			case 3:
				*pOutBase = 't';
				break;
			default:
				*pOutBase = 'n';
				break;
			}
		}

	m_pOutP1->Write((char *)ReadSeq,ReadLen);

	pInBase = &ConcatSeq[ReadLen];

	// need as reverse/complement or RevCpl???
	switch(Mode) {
		case ePMdefault:			// retain split reads in original orientation
			break;
		case ePMReverseR2:			// reverse R2
			CSeqTrans::ReverseSeq(ReadLen,&ConcatSeq[ReadLen]);
			break;

		case ePMComplementR2:		// complement R2
			CSeqTrans::ComplementStrand(ReadLen,&ConcatSeq[ReadLen]);
			break;

		case ePMRevCplR2:			// reverse complement R2
			CSeqTrans::ReverseComplement(ReadLen,&ConcatSeq[ReadLen]);
			break;
		}

	pOutBase = ReadSeq;
	for(BaseIdx = 0; BaseIdx < ReadLen; BaseIdx++,pInBase++,pOutBase++)
		{
		switch(*pInBase & 0x07) {
			case 0:
				*pOutBase = 'a';
				break;
			case 1:
				*pOutBase = 'c';
				break;
			case 2:
				*pOutBase = 'g';
				break;
			case 3:
				*pOutBase = 't';
				break;
			default:
				*pOutBase = 'n';
				break;
			}
		}
	m_pOutP2->Write((char *)ReadSeq,ReadLen);
	}
m_pOutP1->Close();
m_pOutP2->Close();
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing completed, there were %d concatenated reads",NumConcatReads);
Reset();
return(0);
}
