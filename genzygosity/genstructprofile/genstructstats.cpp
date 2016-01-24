// genstructstats.cpp : Defines the entry point for the console application.
//
#include "stdafx.h"
#include <io.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
#include <shlwapi.h>

#include "../conservlib/commhdrs.h"

CStopWatch gStopWatch;

const unsigned int cProgVer = 202;		// increment with each release


int CreateStructValues(char *pszConfParams,char *pszOutFile,bool bSort);

int _tmain(int argc, _TCHAR* argv[])
{
// determine my process name
char szProcName[_MAX_FNAME];
_splitpath(argv[0],NULL,NULL,szProcName,NULL);
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
		printf("Usage: %s ", szProcName);
        arg_print_syntax(stdout,argtable,"\n");
        arg_print_glossary(stdout,argtable,"  %-25s %s\n");
		exit(1);
        }

    /* special case: '--version' takes precedence error reporting */
if (version->count > 0)
        {
		printf("\n%s: Version: %d.%2.2d\n",szProcName,cProgVer/100,cProgVer%100);
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
				szProcName,szOutputFileSpec);

	gStopWatch.Start();

	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
	Rslt = CreateStructValues(szConfParams,szOutputFileSpec,bSort);

	gStopWatch.Stop();
	printf("\n%s: Total processing time: %s\n",szProcName,gStopWatch.Read());
	return(Rslt);
	}
else
	{
	arg_print_errors(stdout,end,szProcName);
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

if((Rslt=Conf.LoadStructParams(pszConfParams))!=eBSFSuccess)
	{
	printf("Unable to open %s - %s",bDefaultConfParams ? "Defaults" : pszConfParams,strerror(errno));
	return(eBSFerrCreateFile);

	}
if((hRsltsFile = open(pszOutFile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
	{
	printf("Unable to create %s - %s",pszOutFile,strerror(errno));
	return(eBSFerrCreateFile);
	}


Len=sprintf(szLine,"\"ID\",\"Octamer\",\"energy\",\"minorgroove\",\"majorgroove\",\"twist\",\"roll\",\"tilt\",\"rise\",\"slide\",\"shift\",\"rmsd\",\"sum3qt\",\"sumqt\",\"isbistable\"\n");
write(hRsltsFile,szLine,Len);


for(Idx = 0; Idx < 0x0ffff; Idx++)
	{
	iSeq2Octamer(Idx,8,Octamer);
	for(Param = 0; Param < eSSNumStatParams; Param++)
		ParamValues[Param] = Conf.StructValue((teStructStats)Param,		// which structural parameter value to return
			4,									// which step in sequence to return structural value for
			8,									// total length of sequence
			Octamer,							// sequence to be processed
			-1);								// value to return for undefined or indeterminate ('N') bases

	Idx2Bases(Idx,8,szOctamer);
	Len = sprintf(szLine,"%d,\"%s\"",Idx,szOctamer);
	for(Param = 0; Param < eSSNumStatParams; Param++)
		Len+=sprintf(&szLine[Len],",%d",ParamValues[Param]);
	Len += sprintf(&szLine[Len],"\n");
	write(hRsltsFile,szLine,Len);
	}
close(hRsltsFile);
return(0);
}

