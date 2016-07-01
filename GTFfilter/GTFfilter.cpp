// GTFfilter.cpp : Defines the entry point for the console application.
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

#include "./ChromMap.h"

const char *cpszProgVer = "1.2.0";		// increment with each release

// processing modes
typedef enum TAG_ePMode {
	ePMdefault,					// default processing mode
	ePMplaceholder				// used to set the enumeration range
	} etPMode;

int Process(etPMode PMode,		// processing mode
			char *pszGTFfile,	// input from this GTF file
			char *pszMapFile,	// input from this mapping file
			char *pszOutFile);	// output to this GTF or GFF file

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
	return _T("GTFfilter");
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
int iScreenLogLevel;			// level of file diagnostics
char szLogFile[_MAX_PATH];	// write diagnostics to this file
int Rslt = 0;   			// function result code >= 0 represents success, < 0 on failure

etPMode PMode;				// processing mode

char szInMapFile[_MAX_PATH];	// input from this contig to chrom mapping file
char szInGtfFile[_MAX_PATH];	// input from this GTF file
char szOutFile[_MAX_PATH];		// output to this GTF or GFF file


// command line args
struct arg_lit  *help    = arg_lit0("h","help",                 "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to screen and logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"diagnostics log file");

struct arg_int *pmode = arg_int0("m","mode","<int>",		    "processing mode: 0 - standard");
struct arg_file *ingtffile = arg_file1("i","in","<file>",		"input from this GTF file");
struct arg_file *inmapfile = arg_file0("I","in","<file>",		"optionally input from this contig to chrom mapping file");
struct arg_file *outfile = arg_file1("o","out","<file>",		"output to this GTF or GFF file");

struct arg_end *end = arg_end(20);

void *argtable[] = {help,version,FileLogLevel,LogFile,
					pmode,ingtffile,inmapfile,outfile,
					end};

char **pAllArgs;
int argerrors;
argerrors = CUtility::arg_parsefromfile(argc,(char **)argv,&pAllArgs);
if(argerrors >= 0)
	argerrors = arg_parse(argerrors,pAllArgs,argtable);

/* special case: '--help' takes precedence over error reporting */
if (help->count > 0)
        {
		printf("\n%s the the GTF file processor, Version %s\nOptions ---\n", gszProcName,cpszProgVer);
        arg_print_syntax(stdout,argtable,"\n");
        arg_print_glossary(stdout,argtable,"  %-25s %s\n");
		printf("\nNote: Parameters can be entered into a parameter file, one parameter per line.");
		printf("\n      To invoke this parameter file then precede it's name with '@'");
		printf("\n      e.g. %s @myparams.txt\n",gszProcName);
		printf("\nPlease report any issues regarding usage of %s to stuart.stephen@csiro.au\n\n",gszProcName);
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
			// now that command parameters have been parsed then initialise diagnostics log system
	if(!gDiagnostics.Open(szLogFile,(etDiagLevel)iScreenLogLevel,(etDiagLevel)iFileLogLevel,true))
		{
		printf("\nError: Unable to start diagnostics subsystem.");
		if(szLogFile[0] != '\0')
			printf(" Most likely cause is that logfile '%s' can't be opened/created",szLogFile);
		exit(1);
		}

	PMode = (etPMode)(pmode->count ? pmode->ival[0] : ePMdefault);
	if(PMode < ePMdefault || PMode >= ePMplaceholder)
		{
		printf("\nError: Processing mode '-m%d' specified outside of range %d..%d",PMode,ePMdefault,(int)ePMplaceholder-1);
		exit(1);
		}

	if(inmapfile->count != 0)
	strcpy(szInMapFile,inmapfile->filename[0]);
	else
		szInMapFile[0] = '\0';
	strcpy(szInGtfFile,ingtffile->filename[0]);
	strcpy(szOutFile,outfile->filename[0]);

// show user current resource limits
#ifndef _WIN32
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Resources: %s",CUtility::ReportResourceLimits());
#endif

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Version: %s Processing parameters:",cpszProgVer);

	const char *pszDescr;
	switch(PMode) {
		case ePMdefault:
			pszDescr = "standard processing";
			break;
		}
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Processing mode is : '%s'",pszDescr);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Process from this GTF file : '%s'",szInGtfFile);
	if(szInMapFile[0] != '\0')
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Use chromosome mapping in this file : '%s'",szInMapFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Write mapped loci to this file : '%s'",szOutFile);

#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	gStopWatch.Start();
	Rslt = Process(PMode,szInGtfFile,szInMapFile,szOutFile);
	gStopWatch.Stop();
	Rslt = Rslt >=0 ? 0 : 1;
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Exit code: %d Total processing time: %s",Rslt,gStopWatch.Read());
	exit(Rslt);
	}
else
	{
	printf("\n%s the GTF file processor, Version %s\n",gszProcName,cpszProgVer);
	arg_print_errors(stdout,end,gszProcName);
	arg_print_syntax(stdout,argtable,"\nUse '-h' to view option and parameter usage\n");
	exit(1);
	}
return 0;
}

CGTFFile *m_pGTFFile;
CChromMap *m_pChromMap;
int m_hRsltsFile;

void
Reset(void)
{
if(m_hRsltsFile != -1)
	{
	close(m_hRsltsFile);
	m_hRsltsFile = -1;
	}
if(m_pGTFFile != NULL)
	{
	delete m_pGTFFile;
	m_pGTFFile = NULL;
	}
if(m_pChromMap != NULL)
	{
	delete m_pChromMap;
	m_pChromMap = NULL;
	}
}

void
Init(void)
{
m_pGTFFile = NULL;
m_pChromMap = NULL;
m_hRsltsFile = -1;
Reset();
}

int Process(etPMode PMode,		// processing mode
			char *pszGTFfile,	// input from this GTF file
			char *pszMapFile,	// input from this mapping file
			char *pszOutFile)	// output to this GTF or GFF file
{
int Rslt;
int NumEntries;
int NumUnmapped;
tsGTFFields *pCurFields;
char szLineBuff[16000];
int BufOfs;
Init();
if(pszMapFile != NULL && pszMapFile[0] != '\0')
	{
m_pChromMap = new CChromMap;
if((Rslt = m_pChromMap->LoadMap(pszMapFile)) <= eBSFSuccess)
	{
	Reset();
	return(Rslt);
		}
	}
else
	m_pChromMap = NULL;
m_pGTFFile = new CGTFFile;
if((Rslt = m_pGTFFile->Open(pszGTFfile)) != eBSFSuccess)
	{
	Reset();
	return(Rslt);
	}


#ifdef _WIN32
if((m_hRsltsFile = open(pszOutFile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
if((m_hRsltsFile = open(pszOutFile, O_RDWR | O_CREAT | O_TRUNC, S_IREAD | S_IWRITE ))==-1)
#endif
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create output file: %s - %s",pszOutFile,strerror(errno));
	Reset();
	return(eBSFerrCreateFile);
	}

NumEntries = 0;
NumUnmapped = 0;
BufOfs = 0;
while((Rslt=m_pGTFFile->NextRecordOfType(eGGTFany)) > 0)
	{
	NumEntries += 1;
	pCurFields = m_pGTFFile->GetFields();

	if(m_pChromMap != NULL)
		{
	if((Rslt = m_pChromMap->Map(pCurFields->szSeqName,(int *)&pCurFields->Start,(int *)&pCurFields->End))!=eBSFSuccess)
		{
		if(NumUnmapped++ < 10)
			gDiagnostics.DiagOut(eDLWarn,gszProcName,"Unable to map contig '%s' onto any chromosome",pCurFields->szSeqName);
			}
		}
	if((BufOfs + 1000) > sizeof(szLineBuff))
		{
		CUtility::SafeWrite(m_hRsltsFile,szLineBuff,BufOfs);
		BufOfs = 0;
		}
	BufOfs += sprintf(&szLineBuff[BufOfs],"%s\t%s\t%s\t%d\t%d\t",pCurFields->szSeqName,pCurFields->szSeqSource,pCurFields->szFeature,pCurFields->Start,pCurFields->End);
	BufOfs += sprintf(&szLineBuff[BufOfs],"%s\n",&pCurFields->szRawLine[pCurFields->ScoreOfs]);
	}

if(Rslt >= eBSFSuccess && (BufOfs + 1000) > sizeof(szLineBuff))
	CUtility::SafeWrite(m_hRsltsFile,szLineBuff,BufOfs);
close(m_hRsltsFile);
m_hRsltsFile = 0;
if(Rslt >= eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Total entries processed %d of which %d were unmappable to any chromosome",NumEntries,NumUnmapped);
	}
Reset();
return(Rslt);
}