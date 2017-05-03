/*
 * CSIRO Open Source Software License Agreement (GPLv3)
 * Copyright (c) 2017, Commonwealth Scientific and Industrial Research Organisation (CSIRO) ABN 41 687 119 230.
 * See LICENSE for the complete license information (https://github.com/csiro-crop-informatics/biokanga/LICENSE)
 * Contact: Alex Whan <alex.whan@csiro.au>
 */

// pacbiokanga.cpp : Defines the entry point for the console application.

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

#include "pacbiokanga.h"

const char *cpszProgVer = "4.3.6";		// increment with each release
const char *cpszProcOverview = "BioKanga PacBio Processing Toolset";

// Subprocesses 
#ifdef USEPBfilter
extern int ProcFilter(int argc, char* argv[]);
#endif
extern int ProcErrCorrect(int argc, char* argv[]);
extern int ProcAssemb(int argc, char* argv[]);
extern int ProcECContigs(int argc, char* argv[]);
extern int ProcSWService(int argc, char* argv[]);
extern int ProcKMerDist(int argc, char* argv[]);

// inplace text cleaning; any leading/trailing or internal quote characters are removed; excessive whitespace is reduced to single
char *
RemoveQuotes(char *pszRawText)
{
char *pSrcChr;
char *pDstChr;
bool bInSpace;
char Chr;
CUtility::TrimQuotedWhitespcExtd(pszRawText);
pSrcChr = pszRawText;
pDstChr = pSrcChr;
bInSpace = false;
while((Chr = *pSrcChr++)!= '\0')
	{
	if(Chr == '\'' || Chr == '"')
		continue;
	if(Chr == ' ' || Chr == '\t')
		{
		if(bInSpace)
			continue;
		bInSpace = true;
		}
	else
		bInSpace = false;
	*pDstChr++ = Chr;
	}
*pDstChr = '\0';
return(pszRawText);
}

tsSubProcess SubProcesses[] = {
#ifdef USEPBfilter
	{"filter","Filter Reads", "Filter PacBio reads for retained hairpins", ProcFilter },
#endif
	{"ecreads","Error Correct Reads", "Error correct PacBio reads", ProcErrCorrect },
	{"contigs","Assemb Contigs","Assemble error corrected PacBio reads into contigs",ProcAssemb},
	{"eccontigs","Error Correct Contigs","Error correct assembled PacBio contigs",ProcECContigs},
	{"swservice","SW Service","Distributed computing SW service provider",ProcSWService},
	{"kmerdist","MAF K-mers","Generate exact matching K-mer distributions from MAF",ProcKMerDist}
	};
const int cNumSubProcesses = (sizeof(SubProcesses) / sizeof(tsSubProcess));

CStopWatch gStopWatch;
CDiagnostics gDiagnostics;				// for writing diagnostics messages to log file
CSQLiteSummaries gSQLiteSummaries;		// for writing processing result summaries to SQLite database
int	gExperimentID = 0;					// SQLite experiment identifier
int gProcessID = 0;						// SQLite process identifier
int	gProcessingID = 0;					// SQLite processing identifier

char gszProcName[_MAX_FNAME];			// this processes name
tsSubProcess *gpszSubProcess;			// selected subprocess

#ifdef UsEtHis
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
	return _T("biokanga");
}
// end of str library required code
#endif
#endif

void
GiveHelpSubProcesses(char *pszProcOverview)
{
int Idx;
printf("Help for %s Version %s\n",gszProcName,cpszProgVer);
printf("  %s\n",pszProcOverview);
printf("  Please specify one of the following subprocesses:\n");
for(Idx = 0; Idx < cNumSubProcesses; Idx++)
	printf("  %s\t\t%s\n",SubProcesses[Idx].pszName,SubProcesses[Idx].pszFullDescr);
printf("To obtain parameter help on any subprocess then enter that subprocess name e.g:\n%s %s -h\n",gszProcName,SubProcesses[0].pszName);
}

int
IsValidSubprocess(char *pszSubProcess)
{
int Idx;
tsSubProcess *pSubProcess;
pSubProcess = SubProcesses;
for(Idx = 0; Idx < cNumSubProcesses; Idx++,pSubProcess++)
	if(!stricmp(pSubProcess->pszName,pszSubProcess))
		return(Idx+1);
return(-1);
}

int
ExecSubProcess(int SubProcID,	// which subprocess as returned by IsValidSubprocess()
		int argc,char *argv[])	// parameters for this subprocess
{
int Idx;
char *pArg;
char *pszSubProc;
char szParam[1024];
if(SubProcID < 0 || SubProcID > cNumSubProcesses)
	{
	printf("\nInvalid SubProcID %d",SubProcID);
	return(-1);
	}

// process parameters and remove subprocess specifier
pArg = argv[1];
for(Idx = 1; Idx < argc; Idx++,pArg++)
	{
	if(pArg[0] == '-')
		{
		if(pArg[1] != 'p')
			continue;
		pszSubProc = &pArg[2];
		}
	else
		pszSubProc = pArg;
	strncpy(szParam,pszSubProc,sizeof(szParam) - 1);
	szParam[sizeof(szParam) - 1] = '\0';
	CUtility::TrimQuotedWhitespcExtd(szParam);
	if(!stricmp(szParam,SubProcesses[SubProcID-1].pszName))
		break;
	}
if(Idx != argc) // if located the subprocess specifier then remove from parameters
	{
	while(Idx < argc)
		{
		pArg = argv[Idx];
		argv[Idx] = argv[Idx+1];
		argv[++Idx] = pArg;
		}
	argc -= 1;
	}
gpszSubProcess = &SubProcesses[SubProcID-1];
srand((unsigned int)time(NULL));
return((*SubProcesses[SubProcID-1].SubFunct)(argc,argv));
}

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

// check if user has specified the subprocess required, each subprocess has it's own set of parameters
// the subprocess requested is either the first word immediately following the process name or the '-p<subprocess' option
// iterate all command line parameters checking for subprocess
int Rslt;
int Idx;
int SubProcID;
char szSubProc[1024];
char *pszSubProc;
char *pArg;

pArg  = (char *)argv[1];
SubProcID = 0;
for(Idx = 1; Idx < argc; Idx++)
	{
	if(pArg[0] == '-')
		{
		if(pArg[1] != 'p')
			continue;
		pszSubProc = &pArg[2];
		}
	else
		pszSubProc = pArg;
	strncpy(szSubProc,pszSubProc,sizeof(szSubProc) - 1);
	szSubProc[sizeof(szSubProc) - 1] = '\0';
	CUtility::TrimQuotedWhitespcExtd(szSubProc);
	if((SubProcID = IsValidSubprocess(szSubProc)) > 0)
		break;
	}

if(SubProcID > 0)
	Rslt = ExecSubProcess(SubProcID,argc,(char **)argv);
else
	{
	GiveHelpSubProcesses((char *)cpszProcOverview);
	Rslt = -1;
	}

return(Rslt);
}



