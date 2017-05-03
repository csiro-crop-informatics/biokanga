// genbiobed.cpp : Defines the entry point for the console application.
//
//

#include "stdafx.h"
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#if _WIN32
#include "../libbiokanga/commhdrs.h"
#else
#include "../libbiokanga/commhdrs.h"
#endif

const char *cpszProgVer = "1.9.1";		// increment with each release

const int cMARawBuffSize = 0x0ffffff;
const int cMALineSize    = 0x0fffff;
const int cMaxExperSeqLen     = 0x0ffffff;

CStopWatch gStopWatch;
CDiagnostics gDiagnostics;				// for writing diagnostics messages to log file
char gszProcName[_MAX_FNAME];			// process name


int CreateBioBed(int NumInputFiles,char *pszInfileSpecs[],char *pszBioBedFile,char *pszDescr,char *pszTitle,int iBEDtype,int MinLen,int MaxLen);

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
int Idx;

char szOutputFileSpec[_MAX_PATH];
char szDescription[cMBSFFileDescrLen];
char szTitle[cMBSFShortFileDescrLen];
int iBEDtype;
int iMinLen;
int iMaxLen;
char *pszBedType;

int NumInputFiles;			// number of input files
char *pszInfileSpecs[cRRMaxInFileSpecs];  // input (wildcards allowed if generating .rsx output) read files


// command line args
struct arg_lit  *help    = arg_lit0("hH","help",                "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_int *ScreenLogLevel=arg_int0("S", "ScreenLogLevel",	"<int>","Level of diagnostics written to logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",	"diagnostics log file");

struct arg_int  *BEDType = arg_int1("b","bedtype","<int>",		"type of BED file, 1 == contains exon detail, 2 == all other BED, 3 == 'genultras.exe' ultra csv file, 4 == 'genultracores.exe' ultracores");
struct arg_file *infiles =  arg_filen("i","in","<file>",0,cRRMaxInFileSpecs,	"input from bed file(s)");
struct arg_file *OutFile = arg_file1("o",NULL,"<file>",			"output to biobed file");
struct arg_str *Descr = arg_str1("d","descr","<string>",		"full description");
struct arg_str *Title = arg_str1("t","title","<string>",		"short title");
struct arg_int *MinLen = arg_int0("m","minlen","<int>",			"minimum length accepted (b=4 only)");
struct arg_int *MaxLen = arg_int0("M","maxlen","<int>",			"maximum length accepted (b=4 only)");

struct arg_end *end = arg_end(20);

void *argtable[] = {help,version,FileLogLevel,ScreenLogLevel,LogFile,infiles,OutFile,Descr,Title,BEDType,MinLen,MaxLen,end};

char **pAllArgs;
int argerrors;
argerrors = CUtility::arg_parsefromfile(argc,(char **)argv,&pAllArgs);
if(argerrors >= 0)
	argerrors = arg_parse(argerrors,pAllArgs,argtable);

/* special case: '--help' takes precedence over error reporting */
if (help->count > 0)
        {
		printf("\n%s the K-mer Adaptive Next Generation Aligner BED File preprocessor, Version %s\nOptions ---\n", gszProcName,cpszProgVer);
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

				// now that log parameters have been parsed then initialise diagnostics log system
	if(!gDiagnostics.Open(szLogFile,(etDiagLevel)iScreenLogLevel,(etDiagLevel)iFileLogLevel,true))
		{
		printf("\nError: Unable to start diagnostics subsystem\n");
		if(szLogFile[0] != '\0')
			printf(" Most likely cause is that logfile '%s' can't be opened/created\n",szLogFile);
		exit(1);
		}
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Version: %s",cpszProgVer);


	if(!infiles->count)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"No input files specified");
		exit(1);
		}

	for(NumInputFiles=Idx=0;NumInputFiles < cRRMaxInFileSpecs && Idx < infiles->count; Idx++)
		{
		pszInfileSpecs[Idx] = NULL;
		if(pszInfileSpecs[NumInputFiles] == NULL)
			pszInfileSpecs[NumInputFiles] = new char [_MAX_PATH];
		strncpy(pszInfileSpecs[NumInputFiles],infiles->filename[Idx],_MAX_PATH);
		pszInfileSpecs[NumInputFiles][_MAX_PATH-1] = '\0';
		CUtility::TrimQuotedWhitespcExtd(pszInfileSpecs[NumInputFiles]);
		if(pszInfileSpecs[NumInputFiles][0] != '\0')
			NumInputFiles++;
		}

	if(!NumInputFiles)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"After removal of whitespace, no input file(s) specified with '-i<filespec>' option)");
		exit(1);
		}

	strcpy(szOutputFileSpec,OutFile->filename[0]);
	strncpy(szTitle,Title->sval[0],cMBSFShortFileDescrLen);
	szTitle[cMBSFShortFileDescrLen-1] = '\0';
	strncpy(szDescription,Descr->sval[0],cMBSFFileDescrLen);
	szDescription[cMBSFFileDescrLen-1] = '\0';

	iBEDtype = BEDType->ival[0];
	if(iBEDtype < eBTGeneExons)
		iBEDtype = eBTGeneExons;
	else
		if(iBEDtype > 4)
			iBEDtype = 4;
	if(iBEDtype == 4)
		{
		iMinLen = MinLen->count > 0 ? MinLen->ival[0] : 0;
		if(iMinLen < 0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Min length: %d must not be negative",iMinLen);
			exit(1);
			}
		iMaxLen = MaxLen->count > 0 ? MaxLen->ival[0] : INT_MAX;
		if(iMaxLen < iMinLen)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Max length: %d must not be less than minimum length: %d",iMaxLen,iMinLen);
			exit(1);
			}
		}
	else
		{
		iMaxLen = INT_MAX;
		iMinLen = 0;
		}


	switch(iBEDtype) {
		case eBTGeneExons:
			pszBedType = (char *)"gene loci + exon detail";
			break;
		case eBTGeneral:
			pszBedType = (char *)"feature loci only";
			break;
		case 3:
			pszBedType = (char *)"two species ultras";
			break;
		case 4:
			pszBedType = (char *)"multispecies hyperconserved cores";
			break;
	}


	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing parameters:");
	for(Idx=0; Idx < NumInputFiles; Idx++)
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"input bed file (%d): '%s'",Idx+1,pszInfileSpecs[Idx]);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Output to bioseq file: '%s'",szOutputFileSpec);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"UCSC .BED (type %d)contains: '%s'",iBEDtype,pszBedType);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Min Length accepted: %d",iMinLen);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Max Length accepted: %d",iMaxLen);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Title text: '%s'",szTitle);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Descriptive text: '%s'",szDescription);
#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif

	gStopWatch.Start();
	Rslt = CreateBioBed(NumInputFiles,pszInfileSpecs,szOutputFileSpec,szDescription,szTitle,iBEDtype,iMinLen,iMaxLen);
	gStopWatch.Stop();
	Rslt = Rslt >=0 ? 0 : 1;
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Exit code: %d Total processing time: %s",Rslt,gStopWatch.Read());
	exit(Rslt);
	}
else
	{
	printf("\n%s the K-mer Adaptive Next Generation Aligner BED File preprocessor, Version %s\n", gszProcName,cpszProgVer);
	arg_print_errors(stdout,end,gszProcName);
	arg_print_syntax(stdout,argtable,"\nUse '-h' to view option and parameter usage\n");
	exit(1);
	}
}


// CreateBioBed
// Create biobed file from specified source file
int
CreateBioBed(int NumInputFiles,char *pszInfileSpecs[],char *pszBioBedFile,char *pszDescr,char *pszTitle,int iBEDtype,int MinLen,int MaxLen)
{
int Rslt;
int Idx;
int NumInputFilesProcessed;
char *pszInFile;

CBEDfile *pBed = new CBEDfile;

Rslt = pBed->Open(pszBioBedFile,(teBEDFeatureType)(iBEDtype > eBTGeneral ? eBTGeneral : iBEDtype),true);
if(Rslt != eBSFSuccess)
	{
	while(pBed->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,pBed->GetErrMsg());

	gDiagnostics.DiagOut(eDLFatal,gszProcName,"\nUnable to open %s",pszBioBedFile);
	delete pBed;
	return(Rslt);
	}

pBed->SetDescription(pszDescr);
pBed->SetTitle(pszTitle);

Rslt = eBSFSuccess;
NumInputFilesProcessed = 0;
CSimpleGlob glob(SG_GLOB_FULLSORT);
for(Idx = 0; Rslt >= eBSFSuccess && Idx < NumInputFiles; Idx++)
	{
	glob.Init();
	if(glob.Add(pszInfileSpecs[Idx]) < SG_SUCCESS)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to glob '%s",pszInfileSpecs[Idx]);
		return(eBSFerrOpnFile);	// treat as though unable to open file
		}

	if(glob.FileCount() <= 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to locate any source reads file matching '%s",pszInfileSpecs[Idx]);
		continue;
		}
	Rslt = eBSFSuccess;
	for (int FileID = 0; Rslt >= eBSFSuccess &&  FileID < glob.FileCount(); ++FileID)
		{
		pszInFile = glob.File(FileID);
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing input file (%d) : %s",NumInputFilesProcessed,pszInFile);
		switch(iBEDtype) {
			case 4:
				Rslt = pBed->ProcessGroupCoreFile(pszInFile,MinLen,MaxLen);
				break;

			case 3:
				Rslt = pBed->ProcessUltraCoreFile(pszInFile);
				break;

			default:
				Rslt = pBed->ProcessBedFile(pszInFile);
				break;
			}


		if(Rslt < eBSFSuccess)
			break;
		NumInputFilesProcessed += 1;
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Cumulative chromosomes: %d, features: %d",pBed->GetProvNumChromosomes(),pBed->GetProvNumFeatures());
		}
	}
	
if(Rslt >= eBSFSuccess)
	{
    gDiagnostics.DiagOut(eDLInfo,gszProcName,"Generating indexes and saving to %s...",pszBioBedFile);
	Rslt = pBed->Close(true);
	if(Rslt >= eBSFSuccess)
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Save completed");
	}
if(Rslt < eBSFSuccess)
	{
	while(pBed->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,pBed->GetErrMsg());
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing errors in %s",pszBioBedFile);
	pBed->Close(false);
	delete pBed;
	return(Rslt);
	}
return(Rslt);
}

