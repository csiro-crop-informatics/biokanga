// genstructstats.cpp : Defines the entry point for the console application.
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

CStopWatch gStopWatch;

const char *cpszProgVer = "2.1.0";		// increment with each release


int CreateStructValues(char *pszConfParams,char *pszOutFile,bool bSort);
CDiagnostics gDiagnostics;				// for writing diagnostics messages to log file
char gszProcName[_MAX_FNAME];			// process name


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
bool bSort;
int Rslt;

char szOutputFileSpec[_MAX_PATH];
char szConfParams[_MAX_PATH];

// command line args
struct arg_lit  *help    = arg_lit0("hH","help",                "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_lit  *Sort   = arg_lit0("s","sort",					"sort order is by flanking inwards - default is sort by octamer");
struct arg_file *InFile = arg_file0("i",NULL,"<file>",			"input structure parameters file");
struct arg_file *OutFile = arg_file1("o",NULL,"<file>",			"output to results file as CSV");
struct arg_end *end = arg_end(20);

void *argtable[] = {help,version,Sort,InFile,OutFile,end};

char **pAllArgs;
int argerrors;
argerrors = CUtility::arg_parsefromfile(argc,(char **)argv,&pAllArgs);
if(argerrors >= 0)
	argerrors = arg_parse(argerrors,pAllArgs,argtable);

    /* special case: '--help' takes precedence over error reporting */
if (help->count > 0)
        {
		printf("\n%s Generate CSV file with conformational structure values with imputed major groove from input structure file, Version %s\nOptions ---\n", gszProcName,cpszProgVer);
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
		printf("\n%s: Version: %s\n",gszProcName,cpszProgVer);
		exit(1);
        }

if (!argerrors)
	{
	bSort = Sort->count ? true : false;
	strcpy(szOutputFileSpec,OutFile->filename[0]);

	if(InFile->count)
		{
		strncpy(szConfParams,InFile->filename[0],_MAX_PATH);
		szConfParams[_MAX_PATH-1] = '\0';
		}
	else
		szConfParams[0] = '\0';

	// creating conformation parameter values
	printf("\n%s: creating central step creating conformation parameter values into CSV file '%s'\n",
				gszProcName,szOutputFileSpec);

	gStopWatch.Start();
#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	Rslt = CreateStructValues(szConfParams,szOutputFileSpec,bSort);

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
}

bool
iSeq2Octamer(int Octamer,int Len,etSeqBase *pOctamer)
{
int Idx;
for(Idx = 0; Idx < Len; Idx++)
	{
	*pOctamer++ = Octamer % 4;
	Octamer /= 4;
	}
return(true);
}

char *Idx2Bases(int StepIdx,int NumBases,char *pszBuffer)
{
char *pDst = pszBuffer;
int CurBase;
int BasMsk = 0x03 << ((NumBases - 1) * 2);
while(NumBases--) {
	CurBase = (StepIdx & BasMsk) >> (NumBases * 2);
	BasMsk >>= 2;
	switch(CurBase) {
		case 0:
			*pDst++ = 'a';
			break;
		case 1:
			*pDst++ = 'c';
			break;
		case 2:
			*pDst++ = 'g';
			break;
		case 3:
			*pDst++ = 't';
			break;
		}
	}
*pDst = '\0';
return(pszBuffer);
}


int 
CreateStructValues(char *pszConfParams,char *pszOutFile,bool bSort)
{
bool bDefaultConfParams;

int Idx;
CConformation Conf;
int ParamValues[eSSNumStatParams];
int Param;
etSeqBase Octamer[8];
char szOctamer[9];
char szLine[512];
int hRsltsFile;
int Len;
int Rslt;
bDefaultConfParams = pszConfParams == NULL || pszConfParams[0] == '\0';

if((Rslt=Conf.LoadStructOctamersParams(pszConfParams))!=eBSFSuccess)
	{
	printf("Unable to open %s - %s",bDefaultConfParams ? "Defaults" : pszConfParams,strerror(errno));
	return(eBSFerrCreateFile);
	}

#ifdef _WIN32
if((hRsltsFile = open(pszOutFile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
if((hRsltsFile = open(pszOutFile, O_RDWR | O_CREAT | O_TRUNC, S_IREAD | S_IWRITE ))==-1)
#endif
	{
	printf("Unable to create %s - %s",pszOutFile,strerror(errno));
	return(eBSFerrCreateFile);
	}


Len=sprintf(szLine,"\"ID\",\"Octamer\",\"energy\",\"minor groove\",\"inferenced major groove\",\"twist\",\"roll\",\"tilt\",\"rise\",\"slide\",\"shift\",\"rmsd\",\"ORChID hydroxyl radical cleavage\"\n");
CUtility::SafeWrite(hRsltsFile,szLine,Len);


for(Idx = 0; Idx < cNumParamOctamers; Idx++)
	{
	iSeq2Octamer(Idx,8,Octamer);
	for(Param = 0; Param < eSSNumStatParams; Param++)
		ParamValues[Param] = Conf.StructValue((teOctStructStats)Param,		// which structural parameter value to return
			4,									// which step in sequence to return structural value for
			8,									// total length of sequence
			Octamer,							// sequence to be processed
			-1);								// value to return for undefined or indeterminate ('N') bases

	Idx2Bases(Idx,8,szOctamer);
	Len = sprintf(szLine,"%d,\"%s\"",Idx,szOctamer);
	for(Param = 0; Param < eSSNumStatParams; Param++)
		Len+=sprintf(&szLine[Len],",%d",ParamValues[Param]);
	Len += sprintf(&szLine[Len],"\n");
	CUtility::SafeWrite(hRsltsFile,szLine,Len);
	}
close(hRsltsFile);
return(0);
}

