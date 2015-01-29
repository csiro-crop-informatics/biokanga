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

// csv2sqlite.cpp : Defines the entry point for the console application.
// generates a sqlite database from input csv file
// the schema of the generated database utilises, where appropriate, denomormised column entities and
// is dependent on the csv file type
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

#include "biokanga.h"

#include "./SQLiteMarkers.h"
#include "./SQLiteDE.h"

// Subprocesses 
int Markers2SQLite(int argc, char* argv[]);
int SNPs2SQLite(int argc, char* argv[]);
int DE2SQLite(int argc, char* argv[]);

// inplace text cleaning; any leading/trailing or internal quote characters are removed; excessive whitespace is reduced to single
char *
SQLRemoveQuotes(char *pszRawText)
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

// worker functions doing parsing from CSV and populating SQLite database
int
ProcessCSV2SQLite(int PMode,					// currently just the one mode...default is to parse from CSV and create/populate SQLite database
				  bool bSafe,					// if true then use indexing on all tables whilst inserting... much slower but perhaps safer if multiple threads ...
			      int CSVtype,					// input CSV file has this format (0: markers, 1: SNPs)
				  char *pszExprName,			// name by which this experiment is identified
				  char *pszExprDescr,			// describes experiment
				  char *pTargAssemb,			// assembly against which aligments for SNP discovery
				  int NumSpecies,				// number of species used in alignments
				  char *pszSpeciesNames[],		// names of species
				  char *pszInFile,				// parse from this input CSV file
				  char *pszDatabase);			// SQLite database file

int
ProcessDECSV2SQLite(int PMode,	// currently just the one mode...default is to parse from CSV and create/populate SQLite database
				  bool bSafe,					// if true then use indexing on all tables whilst inserting... much slower but perhaps safer if multiple threads ...
			      char *pszExprName,			// name by which this experiment is identified
				  char *pszExprDescr,			// describes experiment
				  char *pszCtrlConditions,		// control conditions
				  char *pszExprConditions,		// experiment conditions
				  char *pszInFile,				// parse from this input CSV file
				  char *pszDatabase);			// SQLite database file

#ifdef _WIN32
int Markers2SQLite(int argc, char* argv[])
{
// determine my process name
_splitpath(argv[0],NULL,NULL,gszProcName,NULL);
#else
int
Markers2SQLite(int argc, char** argv)
{
// determine my process name
CUtility::splitpath((char *)argv[0],NULL,gszProcName);
#endif

int iScreenLogLevel;		// level of screen diagnostics
int iFileLogLevel;			// level of file diagnostics
char szLogFile[_MAX_PATH];	// write diagnostics to this file
int Rslt = 0;   			// function result code >= 0 represents success, < 0 on failure
int Idx;					// used as iterator
int PMode;					// processing mode

char szName[cMaxIdntNameLen+1];			// name to identify the experiment from which markers were identified
char szDescr[cMaxIdntDescrLen+1];		// describes experimental conditions

char szAssemb[cMaxIdntNameLen+1];		// targeted assembly against which reads were aligned

char szSpeciesList[(cMaxIdntNameLen + 10) * cMaxExprCultivars];		// individual species parsed from this comma/tab/space separated list
int NumSNPSpecies;						 // number of SNP species used when aligning
char *pszSNPSpecies[cMaxExprCultivars];  // SNP called species or cultivars

char szInFile[_MAX_PATH];	// parse markers from this CSV file
char szOutFile[_MAX_PATH];	// write markers to this created SQLite database

// command line args
struct arg_lit  *help    = arg_lit0("h","help",                 "Print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"Print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to screen and logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"Diagnostics log file");

struct arg_int *mode = arg_int0("m","mode","<int>",				"Processing mode: 0 parse 'biokanga snpmarkers' process generated SNPs into SQLite database  (default 0)");
struct arg_str *name = arg_str1("n","name","<str>",				"Name by which experiment is identified");
struct arg_str *descr = arg_str0("N","descr","<str>",			"Description of experimental conditions");
struct arg_str *assemb = arg_str1("a","assemb","<str>",			"Cultivar/species used as target assembly when originally aligning reads");
struct arg_str *snpspecies = arg_str1("s","snpspecies","<str>", "Cultivar/species names, comma/tab/space separated, order specific, SNP called against targeted assembly");
struct arg_file *infile = arg_file1("i","in","<file>",		    "Input CSV file containing SNPs");
struct arg_file *outfile = arg_file1("o","out","<file>",		"Output SNPs to this SQLite database");

struct arg_end *end = arg_end(200);

void *argtable[] = {help,version,FileLogLevel,LogFile,
	                mode,name,descr,assemb,snpspecies,infile,outfile,
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
		return(1);
		}

	iScreenLogLevel = iFileLogLevel = FileLogLevel->count ? FileLogLevel->ival[0] : eDLInfo;
	if(iFileLogLevel < eDLNone || iFileLogLevel > eDLDebug)
		{
		printf("\nError: FileLogLevel '-l%d' specified outside of range %d..%d\n",iFileLogLevel,eDLNone,eDLDebug);
		return(1);
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
		return(1);
		}
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Subprocess %s Version %s starting",gpszSubProcess->pszName,cpszProgVer);
	PMode = mode->count ? mode->ival[0] : 0;
	if(PMode < 0 || PMode > 1)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Processing mode '-m%d' must be in range 0..1",PMode);
		return(1);
		}

	strncpy(szName,name->sval[0],cMaxIdntNameLen);
	szName[cMaxIdntNameLen]= '\0';
	SQLRemoveQuotes(szName);
	if(strlen(szName) < 1)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Expected experiment name '-n<name>' is empty");
		return(1);
		}
	if(descr->count > 0)
		{
		strncpy(szDescr,descr->sval[0],cMaxIdntDescrLen);
		szDescr[cMaxIdntDescrLen]='\0';
		SQLRemoveQuotes(szDescr);
		}
	if(descr->count == 0 || strlen(szDescr) < 1)
		sprintf(szDescr,"No description of experiment %s supplied",szName);

	strncpy(szAssemb,assemb->sval[0],cMaxIdntNameLen);
	szAssemb[cMaxIdntNameLen]= '\0';
	SQLRemoveQuotes(szAssemb);
	if(strlen(szAssemb) < 1)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Expected assembly name '-a<name>' is empty");
		return(1);
		}

	strncpy(szSpeciesList,snpspecies->sval[0],sizeof(szSpeciesList));
	szSpeciesList[sizeof(szSpeciesList)-1] = '\0';
	SQLRemoveQuotes(szSpeciesList);
	char *pChr = szSpeciesList;
	char *pStartChr;
	char Chr;
	int CurSpeciesLen;
	NumSNPSpecies=0;
	CurSpeciesLen = 0;
	pStartChr = pChr;
	while((Chr = *pChr++) != '\0')
		{
		if(Chr == ' ' || Chr == '\t' || Chr == ',')	// treat these as delimiters
			{
			pChr[-1] = '\0';
			if(CurSpeciesLen != 0)
				{
				pszSNPSpecies[NumSNPSpecies++] = pStartChr;
				CurSpeciesLen = 0;
				}
			pStartChr = pChr;
			continue;
			}
		CurSpeciesLen += 1;
		}
	if(CurSpeciesLen)
		pszSNPSpecies[NumSNPSpecies++] = pStartChr;

	if(!NumSNPSpecies)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Expected at least one SNP species");
		return(1);
		}

	strcpy(szInFile,infile->filename[0]);
	CUtility::TrimQuotedWhitespcExtd(szInFile);

	strcpy(szOutFile,outfile->filename[0]);
	CUtility::TrimQuotedWhitespcExtd(szOutFile);

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing parameters:");
	const char *pszDescr;
	switch(PMode) {
		case 0:
			pszDescr = "Load and parse marker SNPs CSV file and populate SQLite database";
			break;
		}

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Processing mode is : '%s'",pszDescr);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Experiment name : '%s'",szName);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Experiment description : '%s'",szDescr);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Assembly targeted : '%s'",szAssemb);
	for(Idx = 0; Idx < NumSNPSpecies; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"SNP Species (%d) : '%s'",Idx+1,pszSNPSpecies[Idx]);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Marker CSV file: '%s'",szInFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"SQLite database file: '%s'",szOutFile);

#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	gStopWatch.Start();
	Rslt = ProcessCSV2SQLite(PMode,false,0,szName,szDescr,szAssemb,NumSNPSpecies,pszSNPSpecies,szInFile,szOutFile);
	gStopWatch.Stop();
	Rslt = Rslt >=0 ? 0 : 1;
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Exit code: %d Total processing time: %s",Rslt,gStopWatch.Read());
	return(Rslt);
	}
else
	{
    printf("\n%s %s %s, Version %s\n", gszProcName,gpszSubProcess->pszName,gpszSubProcess->pszFullDescr,cpszProgVer);
	arg_print_errors(stdout,end,gszProcName);
	arg_print_syntax(stdout,argtable,"\nUse '-h' to view option and parameter usage\n");
	return(1);
	}
return 0;
}

#ifdef _WIN32
int SNPs2SQLite(int argc, char* argv[])
{
// determine my process name
_splitpath(argv[0],NULL,NULL,gszProcName,NULL);
#else
int
SNPs2SQLite(int argc, char** argv)
{
// determine my process name
CUtility::splitpath((char *)argv[0],NULL,gszProcName);
#endif
int iScreenLogLevel;		// level of screen diagnostics
int iFileLogLevel;			// level of file diagnostics
char szLogFile[_MAX_PATH];	// write diagnostics to this file
int Rslt = 0;   			// function result code >= 0 represents success, < 0 on failure
int Idx;					// used as iterator
int PMode;					// processing mode

char szName[cMaxIdntNameLen+1];			// name to identify the experiment from which markers were identified
char szDescr[cMaxIdntDescrLen+1];		// describes experimental conditions

char szAssemb[cMaxIdntNameLen+1];		// targeted assembly against which reads were aligned

char szSpeciesList[(cMaxIdntNameLen + 10) * cMaxExprCultivars];		// individual species parsed from this comma/tab/space separated list
int NumSNPSpecies;						 // number of SNP species used when aligning
char *pszSNPSpecies[cMaxExprCultivars];  // SNP called species

char szInFile[_MAX_PATH];	// parse markers from this CSV file
char szOutFile[_MAX_PATH];	// write markers to this created SQLite database

// command line args
struct arg_lit  *help    = arg_lit0("h","help",                 "Print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"Print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to screen and logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"Diagnostics log file");

struct arg_int *mode = arg_int0("m","mode","<int>",				"Processing mode: 0 parse 'biokanga align' process generated SNPs into SQLite database  (default 0)");
struct arg_str *name = arg_str1("n","name","<str>",				"Name by which experiment is identified");
struct arg_str *descr = arg_str0("N","descr","<str>",			"Description of experimental conditions");
struct arg_str *assemb = arg_str1("a","assemb","<str>",			"Cultivar/species used as target assembly when aligning reads");
struct arg_str *snpspecies = arg_str1("s","snpspecies","<str>", "Cultivar/species name SNPs called against targeted assembly");
struct arg_file *infile = arg_file1("i","in","<file>",		    "Input CSV file containing SNPs");
struct arg_file *outfile = arg_file1("o","out","<file>",		"Output SNPs to this SQLite database");

struct arg_end *end = arg_end(200);

void *argtable[] = {help,version,FileLogLevel,LogFile,
	                mode,name,descr,assemb,snpspecies,infile,outfile,
					end};

char **pAllArgs;
int argerrors;
argerrors = CUtility::arg_parsefromfile(argc,(char **)argv,&pAllArgs);
if(argerrors >= 0)
	argerrors = arg_parse(argerrors,pAllArgs,argtable);

/* special case: '--help' takes precedence over error reporting */
if (help->count > 0)
        {
		printf("\n%s SNPs to SQLite database, Version %s\nOptions ---\n", gszProcName,cpszProgVer);
        arg_print_syntax(stdout,argtable,"\n");
        arg_print_glossary(stdout,argtable,"  %-25s %s\n");
		printf("\nNote: Parameters can be entered into a parameter file, one parameter per line.");
		printf("\n      To invoke this parameter file then precede it's name with '@'");
		printf("\n      e.g. %s @myparams.txt\n",gszProcName);
		printf("\nPlease report any issues regarding usage of %s to stuart.stephen@csiro.au\n\n",gszProcName);
		return(1);
        }

    /* special case: '--version' takes precedence error reporting */
if (version->count > 0)
        {
		printf("\n%s Version %s\n",gszProcName,cpszProgVer);
		return(1);
        }


if (!argerrors)
	{
	if(FileLogLevel->count && !LogFile->count)
		{
		printf("\nError: FileLogLevel '-f%d' specified but no logfile '-F<logfile>\n'",FileLogLevel->ival[0]);
		return(1);
		}

	iScreenLogLevel = iFileLogLevel = FileLogLevel->count ? FileLogLevel->ival[0] : eDLInfo;
	if(iFileLogLevel < eDLNone || iFileLogLevel > eDLDebug)
		{
		printf("\nError: FileLogLevel '-l%d' specified outside of range %d..%d\n",iFileLogLevel,eDLNone,eDLDebug);
		return(1);
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
		return(1);
		}
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Version: %s",cpszProgVer);
	PMode = mode->count ? mode->ival[0] : 0;
	if(PMode < 0 || PMode > 1)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Processing mode '-m%d' must be in range 0..1",PMode);
		return(1);
		}

	strncpy(szName,name->sval[0],cMaxIdntNameLen);
	szName[cMaxIdntNameLen]= '\0';
	SQLRemoveQuotes(szName);
	if(strlen(szName) < 1)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Expected experiment name '-n<name>' is empty");
		return(1);
		}
	if(descr->count > 0)
		{
		strncpy(szDescr,descr->sval[0],cMaxIdntDescrLen);
		szDescr[cMaxIdntDescrLen]='\0';
		SQLRemoveQuotes(szDescr);
		}
	if(descr->count == 0 || strlen(szDescr) < 1)
		sprintf(szDescr,"No description of experiment %s supplied",szName);

	strncpy(szAssemb,assemb->sval[0],cMaxIdntNameLen);
	szAssemb[cMaxIdntNameLen]= '\0';
	SQLRemoveQuotes(szAssemb);
	if(strlen(szAssemb) < 1)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Expected assembly name '-a<name>' is empty");
		return(1);
		}

	strncpy(szSpeciesList,snpspecies->sval[0],sizeof(szSpeciesList));
	szSpeciesList[sizeof(szSpeciesList)-1] = '\0';
	SQLRemoveQuotes(szAssemb);
	char *pChr = szSpeciesList;
	char *pStartChr;
	char Chr;
	int CurSpeciesLen;
	NumSNPSpecies=0;
	CurSpeciesLen = 0;
	pStartChr = pChr;
	while((Chr = *pChr++) != '\0')
		{
		if(Chr == ' ' || Chr == '\t' || Chr == ',')	// treat these as delimiters
			{
			pChr[-1] = '\0';
			if(CurSpeciesLen != 0)
				{
				pszSNPSpecies[NumSNPSpecies++] = pStartChr;
				CurSpeciesLen = 0;
				}
			pStartChr = pChr;
			continue;
			}
		CurSpeciesLen += 1;
		}
	if(CurSpeciesLen)
		pszSNPSpecies[NumSNPSpecies++] = pStartChr;

	if(NumSNPSpecies != 1)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Expected one SNP species");
		return(1);
		}

	strcpy(szInFile,infile->filename[0]);
	CUtility::TrimQuotedWhitespcExtd(szInFile);

	strcpy(szOutFile,outfile->filename[0]);
	CUtility::TrimQuotedWhitespcExtd(szOutFile);

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing parameters:");
	const char *pszDescr;
	switch(PMode) {
		case 0:
			pszDescr = "Load and parse SNP CSV file and populate SQLite database";
			break;
		}

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Processing mode is : '%s'",pszDescr);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Experiment name : '%s'",szName);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Experiment description : '%s'",szDescr);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Assembly targeted : '%s'",szAssemb);
	for(Idx = 0; Idx < NumSNPSpecies; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"SNP Species (%d) : '%s'",Idx+1,pszSNPSpecies[Idx]);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"SNP CSV file: '%s'",szInFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"SQLite database file: '%s'",szOutFile);

#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	gStopWatch.Start();
	Rslt = ProcessCSV2SQLite(PMode,false,1,szName,szDescr,szAssemb,NumSNPSpecies,pszSNPSpecies,szInFile,szOutFile);
	gStopWatch.Stop();
	Rslt = Rslt >=0 ? 0 : 1;
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Exit code: %d Total processing time: %s",Rslt,gStopWatch.Read());
	return(Rslt);
	}
else
	{
	printf("\n%s SNPs %s, Version %s\n",gszProcName,"Generate SQLite SNP Database",cpszProgVer);
	arg_print_errors(stdout,end,gszProcName);
	arg_print_syntax(stdout,argtable,"\nUse '-h' to view option and parameter usage\n");
	return(1);
	}
return 0;
}


#ifdef _WIN32
int DE2SQLite(int argc, char* argv[])
{
// determine my process name
_splitpath(argv[0],NULL,NULL,gszProcName,NULL);
#else
int
DE2SQLite(int argc, char** argv)
{
// determine my process name
CUtility::splitpath((char *)argv[0],NULL,gszProcName);
#endif

int iScreenLogLevel;		// level of screen diagnostics
int iFileLogLevel;			// level of file diagnostics
char szLogFile[_MAX_PATH];	// write diagnostics to this file
int Rslt = 0;   			// function result code >= 0 represents success, < 0 on failure
int PMode;					// processing mode

char szName[cMaxTransIdntNameLen+1];			// name to identify the experiment from which DE was identified
char szDescr[cMaxTransIdntDescrLen+1];		// describes experiment

char szControlDescr[cMaxTransIdntDescrLen+1];		// describes control conditions
char szExperimentDescr[cMaxTransIdntDescrLen+1];		// describes experiment conditions

char szInFile[_MAX_PATH];	// parse transcripts from this CSV file
char szOutFile[_MAX_PATH];	// write markers to this created SQLite database

// command line args
struct arg_lit  *help    = arg_lit0("h","help",                 "Print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"Print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to screen and logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"Diagnostics log file");

struct arg_int *mode = arg_int0("m","mode","<int>",				"Processing mode: 0 parse 'kangade' process generated DE transcripts into SQLite database  (default 0)");
struct arg_str *name = arg_str1("n","name","<str>",				"Name by which experiment is identified");
struct arg_str *descr = arg_str0("N","descr","<str>",			"Description of experiment");

struct arg_str *control = arg_str1("c","control","<str>",		"Control experimental conditions");
struct arg_str *exper = arg_str1("e","experiment","<str>",		"Experiment experimental conditions");
struct arg_file *infile = arg_file1("i","in","<file>",		    "Input CSV file containing transcripts");
struct arg_file *outfile = arg_file1("o","out","<file>",		"Output transcripts to this SQLite database");

struct arg_end *end = arg_end(200);

void *argtable[] = {help,version,FileLogLevel,LogFile,
	                mode,name,descr,control,exper,infile,outfile,
					end};

char **pAllArgs;
int argerrors;
argerrors = CUtility::arg_parsefromfile(argc,(char **)argv,&pAllArgs);
if(argerrors >= 0)
	argerrors = arg_parse(argerrors,pAllArgs,argtable);

/* special case: '--help' takes precedence over error reporting */
if (help->count > 0)
        {
		printf("\n%s DE to SQLite database, Version %s\nOptions ---\n", gszProcName,cpszProgVer);
        arg_print_syntax(stdout,argtable,"\n");
        arg_print_glossary(stdout,argtable,"  %-25s %s\n");
		printf("\nNote: Parameters can be entered into a parameter file, one parameter per line.");
		printf("\n      To invoke this parameter file then precede it's name with '@'");
		printf("\n      e.g. %s @myparams.txt\n",gszProcName);
		printf("\nPlease report any issues regarding usage of %s to stuart.stephen@csiro.au\n\n",gszProcName);
		return(1);
        }

    /* special case: '--version' takes precedence error reporting */
if (version->count > 0)
        {
		printf("\n%s Version %s\n",gszProcName,cpszProgVer);
		return(1);
        }


if (!argerrors)
	{
	if(FileLogLevel->count && !LogFile->count)
		{
		printf("\nError: FileLogLevel '-f%d' specified but no logfile '-F<logfile>\n'",FileLogLevel->ival[0]);
		return(1);
		}

	iScreenLogLevel = iFileLogLevel = FileLogLevel->count ? FileLogLevel->ival[0] : eDLInfo;
	if(iFileLogLevel < eDLNone || iFileLogLevel > eDLDebug)
		{
		printf("\nError: FileLogLevel '-l%d' specified outside of range %d..%d\n",iFileLogLevel,eDLNone,eDLDebug);
		return(1);
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
		return(1);
		}
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Version: %s",cpszProgVer);
	PMode = mode->count ? mode->ival[0] : 0;
	if(PMode < 0 || PMode > 1)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Processing mode '-m%d' must be in range 0..1",PMode);
		return(1);
		}

	strncpy(szName,name->sval[0],cMaxIdntNameLen);
	szName[cMaxIdntNameLen]= '\0';
	SQLRemoveQuotes(szName);
	if(strlen(szName) < 1)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Expected experiment name '-n<name>' is empty");
		return(1);
		}
	if(descr->count > 0)
		{
		strncpy(szDescr,descr->sval[0],cMaxTransIdntDescrLen);
		szDescr[cMaxTransIdntDescrLen]='\0';
		SQLRemoveQuotes(szDescr);
		}
	if(descr->count == 0 || strlen(szDescr) < 1)
		sprintf(szDescr,"No description of experiment %s supplied",szName);

	strncpy(szControlDescr,control->sval[0],cMaxTransIdntDescrLen);
	szControlDescr[cMaxTransIdntDescrLen]= '\0';
	SQLRemoveQuotes(szControlDescr);
	if(strlen(szControlDescr) < 1)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Expected control conditions '-c<descr>' is empty");
		return(1);
		}

	strncpy(szExperimentDescr,exper->sval[0],cMaxTransIdntDescrLen);
	szExperimentDescr[cMaxTransIdntDescrLen]= '\0';
	SQLRemoveQuotes(szExperimentDescr);
	if(strlen(szExperimentDescr) < 1)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Expected experiment conditions '-e<descr>' is empty");
		return(1);
		}

	strcpy(szInFile,infile->filename[0]);
	CUtility::TrimQuotedWhitespcExtd(szInFile);

	strcpy(szOutFile,outfile->filename[0]);
	CUtility::TrimQuotedWhitespcExtd(szOutFile);

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing parameters:");
	const char *pszDescr;
	switch(PMode) {
		case 0:
			pszDescr = "Load and parse transcript DE CSV file and populate SQLite database";
			break;
		}

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Processing mode is : '%s'",pszDescr);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Experiment name : '%s'",szName);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Experiment description : '%s'",szDescr);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Control conditions : '%s'",szControlDescr);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Experiment conditions : '%s'",szExperimentDescr);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Kangade transcript CSV file: '%s'",szInFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"SQLite database file: '%s'",szOutFile);

#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	gStopWatch.Start();
	Rslt = ProcessDECSV2SQLite(PMode,false,szName,szDescr,szControlDescr,szExperimentDescr,szInFile,szOutFile);
	gStopWatch.Stop();
	Rslt = Rslt >=0 ? 0 : 1;
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Exit code: %d Total processing time: %s",Rslt,gStopWatch.Read());
	return(Rslt);
	}
else
	{
	printf("\n%s SNPs %s, Version %s\n",gszProcName,"Generate SQLite DE Database",cpszProgVer);
	arg_print_errors(stdout,end,gszProcName);
	arg_print_syntax(stdout,argtable,"\nUse '-h' to view option and parameter usage\n");
	return(1);
	}
return 0;
}



// worker function doing parsing from CSV and populating SQLite database
int
ProcessCSV2SQLite(int PMode,					// currently just the one mode...default is to parse from CSV and create/populate SQLite database
				  bool bSafe,					// if true then use indexing on all tables whilst inserting... much slower but perhaps safer if multiple threads ...
			      int CSVtype,					// input CSV file has this format (0: markers, 1: SNPs)
				  char *pszExprName,			// name by which this experiment is identified
				  char *pszExprDescr,			// describes experiment
				  char *pTargAssemb,			// assembly against which aligments for SNP discovery
				  int NumSpecies,				// number of species used in alignments
				  char *pszSpeciesNames[],		// names of species
				  char *pszInFile,				// parse from this input CSV file
				  char *pszDatabase)			// SQLite database file
{
int Rslt;
CSQLiteMarkers *pSQL = NULL;

if((pSQL = new CSQLiteMarkers()) == NULL)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Unable to instantiate instance of CSQLiteMarkers");
	return(eBSFerrObj);
	}

Rslt = pSQL->ProcessCSV2SQLite(PMode,bSafe,CSVtype,pszExprName,pszExprDescr,pTargAssemb,NumSpecies,pszSpeciesNames,pszInFile,pszDatabase);

if(pSQL != NULL)
	delete pSQL;
return(Rslt);
}

// worker function doing DE transcript parsing from CSV and populating SQLite database
int
ProcessDECSV2SQLite(int PMode,	// currently just the one mode...default is to parse from CSV and create/populate SQLite database
				  bool bSafe,					// if true then use indexing on all tables whilst inserting... much slower but perhaps safer if multiple threads ...
			      char *pszExprName,			// name by which this experiment is identified
				  char *pszExprDescr,			// describes experiment
				  char *pszCtrlConditions,		// control conditions
				  char *pszExprConditions,		// experiment conditions
				  char *pszInFile,				// parse from this input CSV file
				  char *pszDatabase)			// SQLite database file
{
int Rslt;
CSQLiteDE *pSQL = NULL;

if((pSQL = new CSQLiteDE()) == NULL)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Unable to instantiate instance of CSQLiteDE");
	return(eBSFerrObj);
	}

Rslt = pSQL->ProcessCSV2SQLite(PMode,bSafe,pszExprName,pszExprDescr,pszCtrlConditions,pszExprConditions,pszInFile,pszDatabase);

if(pSQL != NULL)
	delete pSQL;
return(Rslt);
}

