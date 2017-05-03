// GFFfilter.cpp : Defines the entry point for the console application.
// Filters input GFF3 files (usually as downloaded from TAIR) for a gene type of interest and outputs a GFF file containing
// only those records relevant to that gene type.
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

const unsigned int cProgVer = 100;		// increment with each release

// processing modes
typedef enum TAG_ePMode {
	ePMdefault,					// default processing mode is for GFF format output
	ePMbed,						// output as BED format
	ePMplaceholder				// used to set the enumeration range
	} etPMode;


int Process(etPMode PMode,				// processing mode
		etGGFFGeneType GeneType,		// output this GFF type
		double ScaleFact,				// (BED only) scale score by this factor
		char *pszEntityIdent,			// (BED only) this attribute contains identifying names
		char *pszInFile,				// input GFF file
		char *pszOutFile);				// output file

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
	return _T("GFF2Bed");
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
int Rslt = 0;   			// function result code >= 0 represents success, < 0 on failure

etPMode PMode;				// processing mode
etGGFFGeneType GeneType;
char szName[cMaxGeneNameLen];
double ScaleFact;
char szInFile[_MAX_PATH];
char szOutFile[_MAX_PATH];

// command line args
struct arg_lit  *help    = arg_lit0("h","help",                 "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_int *ScreenLogLevel=arg_int0("S", "ScreenLogLevel",	"<int>","Level of diagnostics written to screen 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"diagnostics log file");

struct arg_int *pmode = arg_int0("m","mode","<int>",		    "processing mode: 0 - output GFF format, 1 - output BED format (default 0)");
struct arg_int *genes = arg_int0("g","genes","<int>",		    "filter records: 0 - all, 1 - protein genes, 2 - transposons, 3 - miRNAs, 4 - snoRNAs, 5 - tRNAs, 6 - pseudogenes (default: 1)");
struct arg_str *name = arg_str0("n","name","<string>",          "if BED output, attribute containing name identifier (default 'Name')");
struct arg_dbl *scalefact=arg_dbl0("s", "scale","<dbl>",			"if BED output, score scaling factor (default 1.0)");
struct arg_file *infile = arg_file1("i","in","<file>",			"input from this GFF file");
struct arg_file *outfile = arg_file1("o","out","<file>",		"output to this file");

struct arg_end *end = arg_end(20);

void *argtable[] = {help,version,FileLogLevel,ScreenLogLevel,LogFile,
					pmode,genes,name,scalefact,infile,outfile,end};

char **pAllArgs;
int argerrors;
argerrors = CUtility::arg_parsefromfile(argc,(char **)argv,&pAllArgs);
if(argerrors >= 0)
	argerrors = arg_parse(argerrors,pAllArgs,argtable);

/* special case: '--help' takes precedence over error reporting */
if (help->count > 0)
        {
		printf("\n%s ", gszProcName);
        arg_print_syntax(stdout,argtable,"\n");
        arg_print_glossary(stdout,argtable,"  %-25s %s\n");
		printf("\nNote: Parameters can be entered into a parameter file, one parameter per line.");
		printf("\n      To invoke this parameter file then precede its name with '@'");
		printf("\n      e.g. %s @myparams.txt\n\n",gszProcName);
		exit(1);
        }

    /* special case: '--version' takes precedence error reporting */
if (version->count > 0)
        {
		printf("\n%s Version %d.%2.2d",gszProcName,cProgVer/100,cProgVer%100);
		exit(1);
        }

if (!argerrors)
	{
	iScreenLogLevel = ScreenLogLevel->count ? ScreenLogLevel->ival[0] : eDLInfo;
	if(iScreenLogLevel < eDLNone || iScreenLogLevel > eDLDebug)
		{
		printf("\nError: ScreenLogLevel '-S%d' specified outside of range %d..%d",iScreenLogLevel,eDLNone,eDLDebug);
		exit(1);
		}
	if(FileLogLevel->count && !LogFile->count)
		{
		printf("\nError: FileLogLevel '-f%d' specified but no logfile '-F<logfile>'",FileLogLevel->ival[0]);
		exit(1);
		}

	iFileLogLevel = FileLogLevel->count ? FileLogLevel->ival[0] : eDLInfo;
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

	PMode = (etPMode)(pmode->count ? pmode->ival[0] : 0);
	if(PMode < 0 || PMode >= ePMplaceholder)
		{
		printf("\nError: Processing mode '-m%d' specified outside of range %d..%d",PMode,0,(int)ePMplaceholder-1);
		exit(1);
		}

	if(PMode == ePMdefault)
		{
		GeneType = (etGGFFGeneType)(genes->count ? genes->ival[0] : 1);
		if(GeneType < 0 || GeneType >= eGGFFplaceholder)
			{
			printf("\nError: Gene type '-g%d' specified outside of range %d..%d",GeneType,0,(int)eGGFFplaceholder-1);
			exit(1);
			}
		szName[0] = '\0';
		ScaleFact = 1.0f;
		}
	else
		{
		GeneType = eGGFFany;
		ScaleFact = (double)(scalefact->count ? scalefact->dval[0] : 1.0f);
		if(name->count)
			{
			strncpy(szName,name->sval[0],cMaxGeneNameLen);
			szName[cMaxGeneNameLen-1]='\0';
			}
		else
			strcpy(szName,"Name");
		}

	strcpy(szOutFile,outfile->filename[0]);
	strcpy(szInFile,infile->filename[0]);

				// now that command parameters have been parsed then initialise diagnostics log system
	if(!gDiagnostics.Open(szLogFile,(etDiagLevel)iScreenLogLevel,(etDiagLevel)iFileLogLevel,true))
		{
		printf("\nError: Unable to start diagnostics subsystem.");
		if(szLogFile[0] != '\0')
			printf(" Most likely cause is that logfile '%s' can't be opened/created",szLogFile);
		exit(1);
		}

// show user current resource limits
#ifndef _WIN32
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Resources: %s",CUtility::ReportResourceLimits());
#endif

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Version: %d.%2.2d Processing parameters:",cProgVer/100,cProgVer%100);

	const char *pszDescr;
	switch(PMode) {
		case ePMdefault:
			pszDescr = "Output GFF format";
			break;
		case ePMbed:
			pszDescr = "Output BED format";
			break;
		}
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Processing mode is : '%s'",pszDescr);

	switch(GeneType) {
		case eGGFFany:	
			pszDescr = "all GFF records";
			break;
		case eGGFFgene:			// only protein gene related records
			pszDescr = "protein genes";
			break;
		case eGGFFtransposon:	// only transposon related records
			pszDescr = "transposons";
			break;
		case eGGFFmiRNA:		// only miRNA related records
			pszDescr = "miRNAs";
			break;
		case eGGFFsnoRNA:		// only snoRNA related records
			pszDescr = "snoRNAs";
			break;
		case eGGFFtRNA:			// only tRNA related records
			pszDescr = "tRNAs";
			break;
		case eGGFFpseudogene:	// only pseudogene related records
			pszDescr = "pseudogenes";
			break;
		case eGGFFncRNA:		// only ncRNA (also referenced as other_RNA)
			pszDescr = "ncRNA/other_RNA";
			break;
		}
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Outputing '%s' only",pszDescr);

	if(PMode == ePMbed)
		{
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Score scaling factor: %1.4f",ScaleFact);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Use name from this attribute: '%s'",szName);
		}

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Input GFF file: '%s'",szInFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"output to file: '%s'",szOutFile);

#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	gStopWatch.Start();
	Rslt = Process(PMode,GeneType,ScaleFact,szName,szInFile,szOutFile);
	gStopWatch.Stop();
	Rslt = Rslt >=0 ? 0 : 1;
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Exit code: %d Total processing time: %s",Rslt,gStopWatch.Read());
	exit(Rslt);
	}
else
	{
	arg_print_errors(stdout,end,gszProcName);
	arg_print_syntax(stdout,argtable,"\nend of help\n");
	exit(1);
	}
return 0;
}

void
Reset(void)
{
}

int
Process(etPMode PMode,					// processing mode
		etGGFFGeneType GeneType,		// output this GFF type
		double ScaleFact,				// (BED only) scale score by this factor
		char *pszEntityIdent,			// (BED only) this attribute contains identifying names
		char *pszInFile,				// input GFF file
		char *pszOutFile)				// output file
{
int Rslt;
int hOutFile;
char szBuff[cGFFMaxLineLen*2];
int BuffIdx;
int NumRecs;
int Start;
int End;
CGFFFile *pGFF = NULL;

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Opening GFF file '%s'", pszInFile);
if((pGFF = new CGFFFile()) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CGFFFile");
	Reset();
	return(eBSFerrObj);
	}

if((Rslt=pGFF->Open(pszInFile))!=eBSFSuccess)
	{
	while(pGFF->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,pGFF->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open input GFF file '%s'",pszInFile);
	Reset();
	return(Rslt);
	}


gDiagnostics.DiagOut(eDLInfo,gszProcName,"Creating/truncating output %s file '%s'", PMode == ePMbed ? "BED" : "GFF", pszOutFile);

#ifdef _WIN32
hOutFile = open(pszOutFile, O_CREATETRUNC);
#else
if((hOutFile = open64(pszOutFile,O_RDWR | O_CREAT,S_IREAD | S_IWRITE))!=-1)
     if(ftruncate(hOutFile,0)!=0)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate output file '%s' - %s",pszOutFile,strerror(errno));
	delete pGFF;
	close(hOutFile);
	return(eBSFerrCreateFile);
	}
#endif

if(hOutFile == -1)					// check if file open succeeded
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open/create output file '%s' - %s",pszOutFile,strerror(errno));
	delete pGFF;
	return(eBSFerrCreateFile);
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing...");

NumRecs = 0;
BuffIdx = 0;
while((Rslt = pGFF->NextRecordOfType(GeneType)) > 0)
	{
	if(PMode == ePMdefault)
		{
		if(NumRecs == 0)
			{
			BuffIdx = sprintf(szBuff,"##gff-version %d\n",pGFF->GetGFFversion());
			gDiagnostics.DiagOut(eDLInfo,gszProcName,szBuff);
			}
		NumRecs += 1;
		strcpy(&szBuff[BuffIdx],pGFF->GetRecord());
		BuffIdx += pGFF->GetCurLineLen();
		szBuff[BuffIdx++] = '\n';
		}
	else	// output in BED format
		{
		char szName[cMaxGeneNameLen];
		char *pszIdent;
		char *pszName;
		int NameLen;
		if(!stricmp(pszEntityIdent,"Source"))
			pszIdent = pGFF->GetSource();
		else
			if((pszIdent = pGFF->GetNoteValue(pszEntityIdent))==NULL)
				continue;
		// cleanup the identifier, remove any quotes, whitespace etc.
		NameLen = 0;
		pszName = &szName[0];
		while(NameLen < (sizeof(szName)-1) && (*pszName = *pszIdent++))
			{
			if(*pszName == '\'' || *pszName == '"')	// simply slough any quotes
				continue;
		
			if(*pszName != '.' && !isalnum(*pszName))
				{
				if(!NameLen || pszName[-1] == '_')
					continue;
				*pszName = '_';
				}
			NameLen += 1;
			pszName += 1;
			}
		*pszName = '\0';
	
		// do not know why, but have observed GFF records with starts > ends!
		Start = pGFF->GetStart();
		End = pGFF->GetEnd();
		if(Start > End)
			{
			End = Start;
			Start = pGFF->GetEnd();
			}
		Start -= 1;
		BuffIdx += sprintf(&szBuff[BuffIdx],"%s\t%d\t%d\t%s\t%d\t%c\n",
					    pGFF->GetSeqName(),
						Start,
						End,
						szName,
						pGFF->GetScore() < 0.0f ? 0 : (int)(pGFF->GetScore() * ScaleFact),
						pGFF->GetStrand() == '?' ? '+' : pGFF->GetStrand());
		NumRecs += 1;

		}
	if((BuffIdx + cGFFMaxLineLen) >= sizeof(szBuff))
		{
		CUtility::SafeWrite(hOutFile,szBuff,BuffIdx);
		BuffIdx = 0;
		}
	}
if(Rslt != eBSFSuccess)
	{
	while(pGFF->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,pGFF->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Stopped processing due to errors");
	}
else
	if(BuffIdx)
		CUtility::SafeWrite(hOutFile,szBuff,BuffIdx);
pGFF->Close();
delete pGFF;
if(hOutFile != -1)
	close(hOutFile);
if(Rslt == eBSFSuccess)
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing completed, wrote a total of %d records",NumRecs);
return(Rslt);
}

