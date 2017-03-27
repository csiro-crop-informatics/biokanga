/*
 * CSIRO Open Source Software License Agreement (variation of the BSD / MIT License)
 * Copyright (c) 2017, Commonwealth Scientific and Industrial Research Organisation (CSIRO) ABN 41 687 119 230.
 * See LICENCE for the complete licence information (https://github.com/csiro-crop-informatics/biokanga/LICENSE)
 * Contact: Alex Whan <alex.whan@csiro.au>
 */

// genpseudogenome.cpp : Defines the entry point for the console application.
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

const int cMinNSeps = 5;				// allow a minium of 5 N's to be specified as being the length of sequence separators
const int cDfltNSeps = 100;				// default to 100 N's to be specified as being the length of sequence separators
const int cMaxNSeps = 100;				// allow at most 100 N's to be specified as being the length of sequence separators

// processing modes
typedef enum TAG_ePMode {
	ePMdefault,					// default processing mode
	ePMplaceholder				// used to set the enumeration range
	} etPMode;

// output format modes
typedef enum TAG_eFMode {
	eFMdefault,					// default
	eFMplaceholder				// used to set the enumeration rangeP
	} etFMode;

const int cMaxInFileSpecs = 100; // can handle up to this many input files

int
GPGProcess(etPMode PMode,					// processing mode
		etFMode FMode,					// output format mode
		int LenNSeps,					// generate with this number of 'N' bases separating concatenated sequences 
		char *pszTrackTitle,			// track title for output UCSC BED
		int NumInputFileSpecs,		  	// number of input file specs
		char *pszInfileSpecs[],		  	// names of inputs files containing multifasta
		char *pszOutGenomeFile,			// output pseudo genome file
		char *pszOutGeneFile);			// output pseudo gene file

#ifdef _WIN32
int genpseudogenome(int argc, char* argv[])
{
// determine my process name
_splitpath(argv[0],NULL,NULL,gszProcName,NULL);
#else
int
genpseudogenome(int argc, char** argv)
{
// determine my process name
CUtility::splitpath((char *)argv[0],NULL,gszProcName);
#endif

int iFileLogLevel;			// level of file diagnostics
int iScreenLogLevel;		// level of file diagnostics
char szLogFile[_MAX_PATH];	// write diagnostics to this file
int Rslt = 0;   			// function result code >= 0 represents success, < 0 on failure
int Idx;

int LenNSeps;				// generate with this number of 'N' bases separating concatenated sequences 
int PMode;				// processing mode
int FMode;				// format output mode

char szTrackTitle[cMaxDatasetSpeciesChrom];		// track title if output format is UCSC BED

char szOutGenomeFile[_MAX_PATH];			// pseudo genome to this file
char szOutGeneFile[_MAX_PATH];			// pseudo genes to this file
int NumInputFiles;						// number of input files
char *pszInfileSpecs[cMaxInFileSpecs];  // input (wildcards allowed) multifasta files

char szSQLiteDatabase[_MAX_PATH];	// results summaries to this SQLite file
char szExperimentName[cMaxDatasetSpeciesChrom+1];			// experiment name
char szExperimentDescr[1000];		// describes experiment

// command line args
struct arg_lit  *help    = arg_lit0("h","help",                 "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to screen and logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"diagnostics log file");

struct arg_int *pmode = arg_int0("m","mode","<int>",		    "generate processing mode: 0 - standard");
struct arg_int *format = arg_int0("M","format","<int>",		    "output format: 0 - output genome as single concatenated fasta)");
struct arg_int *lennseps = arg_int0("n","lennseps","<int>",		    "generate with this number of 'N' bases separating concatenated sequences (default 100, range 5..100)");
struct arg_file *infiles = arg_filen("i","in","<file>",0,cMaxInFileSpecs,"input from these multifasta files, wildcards allowed");
struct arg_file *outgenomefile = arg_file1("o","out","<file>",	"output pseudo genome to this file");
struct arg_file *outgenefile = arg_file1("O","outbed","<file>",	"Output pseudo gene (BED) file");

struct arg_str  *title = arg_str0("t","title","<string>",       "track title");
struct arg_file *summrslts = arg_file0("q","sumrslts","<file>",		"Output results summary to this SQLite3 database file");
struct arg_str *experimentname = arg_str0("w","experimentname","<str>",		"experiment name SQLite3 database file");
struct arg_str *experimentdescr = arg_str0("W","experimentdescr","<str>",	"experiment description SQLite3 database file");

struct arg_end *end = arg_end(200);

void *argtable[] = {help,version,FileLogLevel,LogFile,
					summrslts,experimentname,experimentdescr,
					pmode,lennseps,format,title,infiles,	outgenomefile,outgenefile,
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
		printf("\nPlease report any issues regarding usage of %s at https://github.com/csiro-crop-informatics/biokanga/issues\n\n",gszProcName);
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

	PMode = (etPMode)(pmode->count ? pmode->ival[0] : ePMdefault);
	if(PMode < ePMdefault || PMode >= ePMplaceholder)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Processing mode '-m%d' specified outside of range %d..%d\n",PMode,ePMdefault,(int)ePMplaceholder-1);
		exit(1);
		}

	FMode = (etFMode)(format->count ? format->ival[0] : eFMdefault);

	LenNSeps = lennseps->count ? lennseps->ival[0] : cDfltNSeps;
	if(LenNSeps < cMinNSeps || LenNSeps > cMaxNSeps)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Length of 'N' sequence separators '-n%d' specified outside of range %d..%d\n",LenNSeps,cMinNSeps,cMaxNSeps);
		exit(1);
		}


	FMode = (etFMode)(format->count ? format->ival[0] : eFMdefault);
	if(FMode < eFMdefault || FMode >= eFMplaceholder)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Output format mode '-m%d' specified outside of range %d..%d\n",FMode,eFMdefault,(int)eFMplaceholder-1);
		exit(1);
		}


	szTrackTitle[0] = '\0';
	if(title->count)
		{
		strncpy(szTrackTitle,title->sval[0],sizeof(szTrackTitle));
		szTrackTitle[sizeof(szTrackTitle)-1] = '\0';
		CUtility::TrimQuotedWhitespcExtd(szTrackTitle);
		CUtility::ReduceWhitespace(szTrackTitle);
		}
	if(szTrackTitle[0] == '\0')
		{
		gDiagnostics.DiagOut(eDLWarn,gszProcName,"Warning: output format requested to be UCSC BED but no track title with '-t<title' specified, defaulting to 'pseudo'\n");
		strcpy(szTrackTitle,"pseudo");
		}


	strcpy(szOutGeneFile,outgenefile->filename[0]);
	strcpy(szOutGenomeFile,outgenomefile->filename[0]);

	for(NumInputFiles=Idx=0;NumInputFiles < cMaxInFileSpecs && Idx < infiles->count; Idx++)
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
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: After removal of whitespace, no input file(s) specified with '-i<filespec>' option)\n");
		exit(1);
		}

// show user current resource limits
#ifndef _WIN32
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Resources: %s",CUtility::ReportResourceLimits());
#endif

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Version: %s Processing parameters:",cpszProgVer);

	const char *pszDescr;
	switch(PMode) {
		case ePMdefault:
			pszDescr = "Standard processing";
			break;
		default:
			pszDescr = "Default processing";
			break;
		}
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Processing mode is : '%s'",pszDescr);

	switch(FMode) {
		case eFMdefault:
			pszDescr = "output as single pseudo chrom genome fasta and associated gene BED file";
			break;
		}

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Length of inter-sequence separator N's: %d",LenNSeps);
	for(Idx = 0; Idx < NumInputFiles; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Process from input file (%d): '%s'",Idx+1,pszInfileSpecs[Idx]);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Output pseudo genome file : '%s'",szOutGenomeFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Output pseudo gene (BED) file : '%s'",szOutGeneFile);

	if(szExperimentName[0] != '\0')
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"This processing reference: %s",szExperimentName);

	if(gExperimentID > 0)
		{
		int ParamID;
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szLogFile),"log",szLogFile);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,sizeof(INT32),"mode",&PMode);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,sizeof(INT32),"format",&FMode);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,sizeof(INT32),"lennseps",&LenNSeps);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,sizeof(INT32),"NumInputFiles",&NumInputFiles);
		for(Idx = 0; Idx < NumInputFiles; Idx++)
			ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(pszInfileSpecs[Idx]),"in",pszInfileSpecs[Idx]);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szOutGenomeFile),"out",szOutGenomeFile);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szOutGeneFile),"outbed",szOutGeneFile);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szTrackTitle),"title",szTrackTitle);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szSQLiteDatabase),"sumrslts",szSQLiteDatabase);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szExperimentName),"experimentname",szExperimentName);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szExperimentDescr),"experimentdescr",szExperimentDescr);
		}


#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	gStopWatch.Start();
	Rslt = GPGProcess((etPMode)PMode,(etFMode)FMode,LenNSeps,szTrackTitle,NumInputFiles,pszInfileSpecs,szOutGenomeFile,szOutGeneFile);
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


int m_hInSrcFile;		// currently opened input file
int m_hOutGenomeFile; // currently opened output pseudo genome file
int m_hOutGeneFile;		// currently opened output bed file


const int cAllocInFasta = 300000000;			// buffer fasta input as being the maximal sized expected single contig or scaffold
const int cAllocOutFasta = 1000000;				// buffer output gene (BED) file
const int cAllocInBED = 1000000;				// buffer output gene (BED) file

UINT8 *m_pInFastaBuff;							// allocd to hold input fasta
char *m_pOutFastaBuff;							// allocd to hold buffered output fasta
char *m_pOutBEDBuff;							// allocd to hold buffred output BED 
int m_OutFastaOfs;								// bytes used m_pOutFastaBuff
int m_OutBEDOfs;								// bytes used m_pOutFastaBuff

UINT32 m_GenomeLen;								// total pseudo genome length
int m_TotEntries;								// total number of contigs or scaffolds in pseudo genome
etSeqBase m_100Ns[cMaxNSeps];					// to hold upto 100 eBaseNs used to separate contigs/scaffolds in output pseudo genome
int m_CurFastaCol;								// next fasta col to write into


void
GPGReset(void)
{
if(m_hInSrcFile != -1)
	{
	close(m_hInSrcFile);
	m_hInSrcFile = -1;
	}
if(m_hOutGenomeFile != -1)
	{
	close(m_hOutGenomeFile);
	m_hOutGenomeFile = -1;
	}
if(m_hOutGeneFile != -1)
	{
	close(m_hOutGeneFile);
	m_hOutGeneFile = -1;
	}
if(m_pInFastaBuff != NULL)
	{
	delete m_pInFastaBuff;
	m_pInFastaBuff = NULL;
	}
if(m_pOutFastaBuff != NULL)
	{
	delete m_pOutFastaBuff;
	m_pOutFastaBuff = NULL;
	}
if(m_pOutBEDBuff != NULL)
	{
	delete m_pOutBEDBuff;
	m_pOutBEDBuff = NULL;
	}

m_OutFastaOfs = 0;
m_OutBEDOfs = 0;
m_TotEntries = 0;
m_CurFastaCol = 0;
m_GenomeLen = 0;
memset(m_100Ns,eBaseN,sizeof(m_100Ns));
}

void
GPGInit(void)
{
m_hInSrcFile = -1;
m_hOutGenomeFile = -1;
m_hOutGeneFile = -1;
m_pInFastaBuff = NULL;
m_pOutFastaBuff = NULL;
m_pOutBEDBuff = NULL;
GPGReset();
}

int
OutputFasta(int SeqLen,			// sequence length
			etSeqBase *pSeq)	// sequence
{
int NumCols;
int ReadOfs = 0;
int	ReadLenRem = SeqLen;
while(SeqLen)
	{
	if(m_CurFastaCol == 70)
		{
		m_OutFastaOfs += sprintf(&m_pOutFastaBuff[m_OutFastaOfs],"\n");
		m_CurFastaCol = 0;
		}
	NumCols = min(SeqLen, 70 - m_CurFastaCol);
	CSeqTrans::MapSeq2Ascii(&pSeq[ReadOfs],NumCols,&m_pOutFastaBuff[m_OutFastaOfs],'N','U','I',true);
	m_OutFastaOfs += NumCols;
	m_CurFastaCol += NumCols;
	ReadOfs += NumCols;
	SeqLen -= NumCols;
	if((m_OutFastaOfs + 1000) > cAllocOutFasta)
		{
		CUtility::SafeWrite(m_hOutGenomeFile,m_pOutFastaBuff,m_OutFastaOfs);
		m_OutFastaOfs = 0;
		}
	}
return(eBSFSuccess);
}

int
OutputGene(char *pszGenome,
			char *pszName,
           char Strand,
		   UINT32 GeneStart,
		   UINT32 GeneLen)
{

if(m_OutBEDOfs + 1000  > cAllocInBED)
	{
	CUtility::SafeWrite(m_hOutGeneFile,m_pOutBEDBuff,m_OutBEDOfs);
	m_OutBEDOfs = 0;
	}

m_OutBEDOfs += sprintf(&m_pOutBEDBuff[m_OutBEDOfs],"%s\t%d\t%d\t%s\t%d\t%c",pszGenome,GeneStart,GeneStart+GeneLen,pszName,0,Strand);
m_OutBEDOfs += sprintf(&m_pOutBEDBuff[m_OutBEDOfs],"\t%d\t%d\t%d\t%d\t%d,\t%d,\n",GeneStart,GeneStart+GeneLen,0,1,GeneLen,0);
return(eBSFSuccess);
}



int 
LoadFasta(int LenNSeps,					// generate with this number of 'N' bases separating concatenated sequences
			char *pszGenome,char *pszFastaFile)
{
CFasta Fasta;
char szName[cBSFSourceSize];
char szDescription[cBSFDescriptionSize];
UINT32 GeneStart;

int SeqLen;
int Descrlen;
bool bFirstEntry;
bool bEntryCreated;
int Rslt;
int SeqID;

int ChromID;
size_t TotLen;

TotLen = 0;
ChromID = 0;

if((Rslt=Fasta.Open(pszFastaFile,true))!=eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open '%s' [%s] %s",pszFastaFile,Fasta.ErrText((teBSFrsltCodes)Rslt),Fasta.GetErrMsg());
	GPGReset();
	return(Rslt);
	}
bFirstEntry = true;
bEntryCreated = false;
SeqID = 0;
ChromID = 0;
TotLen = 0;
while((Rslt = SeqLen = Fasta.ReadSequence(m_pInFastaBuff,cAllocInFasta,true,false)) > eBSFSuccess)
	{
	if(!m_TotEntries || SeqLen == eBSFFastaDescr)		// just read a descriptor line or else the first entry without a descriptor
		{
		m_TotEntries += 1;
		if(SeqLen == eBSFFastaDescr)
			Descrlen = Fasta.ReadDescriptor(szDescription,cBSFDescriptionSize);
		if(SeqLen != eBSFFastaDescr || sscanf(szDescription," %s[ ,]",szName)!=1)
			sprintf(szName,"%s.%d",pszFastaFile,m_TotEntries);
		if(SeqLen == eBSFFastaDescr)
			continue;
		}

	// if not the first entry then separate with SepLen 'N's
	if(m_TotEntries > 1)
		{
		if((Rslt=OutputFasta(LenNSeps,m_100Ns)) < eBSFSuccess)
			break;
		m_GenomeLen += LenNSeps;
		}

	GeneStart = m_GenomeLen;
	// output sequence to fasta here
	if((Rslt=OutputFasta(SeqLen,m_pInFastaBuff)) < eBSFSuccess)
		break;
	m_GenomeLen += SeqLen;
	OutputGene(pszGenome,szName,'+',GeneStart,m_GenomeLen - GeneStart);
	}

if(Rslt < eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadFasta [%s] %s",Fasta.ErrText((teBSFrsltCodes)Rslt),Fasta.GetErrMsg());
	GPGReset();
	return(Rslt);
	}

return(eBSFSuccess);
}

int
GPGProcess(etPMode PMode,					// processing mode
		etFMode FMode,					// output format mode
		int LenNSeps,					// generate with this number of 'N' bases separating concatenated sequences 
		char *pszGenomeName,			// track title for output UCSC BED and psudeo genome fasta descriptor
		int NumInputFileSpecs,		  	// number of input file specs
		char *pszInfileSpecs[],		  	// names of inputs files containing multifasta
		char *pszOutGenomeFile,			// output pseudo genome file
		char *pszOutGeneFile)			// output pseudo gene file
{
int Rslt;
int Idx;
char *pszInfile;

GPGInit();

#ifdef _WIN32
if((m_hOutGenomeFile = open(pszOutGenomeFile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
if((m_hOutGenomeFile = open(pszOutGenomeFile,O_RDWR | O_CREAT |O_TRUNC, S_IREAD | S_IWRITE ))==-1)
#endif
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"OpenDumpFile:Unable to create or truncate  output file %s error: %s",pszOutGenomeFile,strerror(errno));
	GPGReset();
	return(false);
	}

#ifdef _WIN32
if((m_hOutGeneFile = open(pszOutGeneFile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
if((m_hOutGeneFile = open(pszOutGeneFile,O_RDWR | O_CREAT |O_TRUNC, S_IREAD | S_IWRITE ))==-1)
#endif
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: Unable to create or truncate  output file %s error: %s",pszOutGeneFile,strerror(errno));
	GPGReset();
	return(false);
	}

if((m_pInFastaBuff = new unsigned char [cAllocInFasta]) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessFastaFile:- Unable to allocate memory (%d bytes) for sequence buffer",cAllocInFasta);
	GPGReset();
	return(eBSFerrMem);
	}

if((m_pOutFastaBuff = new char [cAllocOutFasta]) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessFastaFile:- Unable to allocate memory (%d bytes) for sequence buffer",cAllocOutFasta);
	GPGReset();
	return(eBSFerrMem);
	}

if((m_pOutBEDBuff = new char [cAllocInBED]) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessFastaFile:- Unable to allocate memory (%d bytes) for sequence buffer",cAllocInBED);
	GPGReset();
	return(eBSFerrMem);
	}

m_GenomeLen = 0;
m_TotEntries = 0;
m_OutFastaOfs = sprintf(m_pOutFastaBuff,">%s\n",pszGenomeName);
m_OutBEDOfs = sprintf(m_pOutBEDBuff,"track type=bed name=\"%s\" description=\"%s\"\n",pszGenomeName,pszGenomeName);


// iterate the input files
CSimpleGlob glob(SG_GLOB_FULLSORT);
for(Idx = 0; Idx < NumInputFileSpecs; Idx++)
	{
	glob.Init();
	if(glob.Add(pszInfileSpecs[Idx]) < SG_SUCCESS)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to glob '%s",pszInfileSpecs[Idx]);
		GPGReset();
		return(eBSFerrOpnFile);
		}

	if(glob.FileCount() <= 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to locate any input source contig or scaffold fasta files matching '%s",pszInfileSpecs[Idx]);
		continue;
		}

	Rslt = eBSFSuccess;
	for (int FileID = 0; Rslt >= eBSFSuccess &&  FileID < glob.FileCount(); ++FileID)
		{
		pszInfile = glob.File(FileID);
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Loading and parsing contig or scaffolds from file '%s'\n",pszInfile);
		Rslt = LoadFasta(LenNSeps,pszGenomeName,pszInfile);
		if(Rslt != eBSFSuccess)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: Load failed for input contig or scaffolds file '%s'\n",pszInfile);
			GPGReset();
			return(Rslt);
			}
		}
	}
if(m_CurFastaCol != 0)
	m_OutFastaOfs += sprintf(&m_pOutFastaBuff[m_OutFastaOfs],"\n");
if(m_hOutGenomeFile != -1)
	{
	if(m_OutFastaOfs > 0)
		{
		CUtility::SafeWrite(m_hOutGenomeFile,m_pOutFastaBuff,m_OutFastaOfs);
		m_OutFastaOfs = 0;
		}
	close(m_hOutGenomeFile);
	m_hOutGenomeFile = -1;
	}

if(m_hOutGeneFile != -1)
	{
	if(m_OutBEDOfs)
		{
		CUtility::SafeWrite(m_hOutGeneFile,m_pOutBEDBuff,m_OutBEDOfs);
		m_OutBEDOfs = 0;
		}
	close(m_hOutGeneFile);
	m_hOutGeneFile = -1;
	}

GPGReset();
return(eBSFSuccess);
}
