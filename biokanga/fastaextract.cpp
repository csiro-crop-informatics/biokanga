/*
 * CSIRO Open Source Software License Agreement (variation of the BSD / MIT License)
 * Copyright (c) 2017, Commonwealth Scientific and Industrial Research Organisation (CSIRO) ABN 41 687 119 230.
 * See LICENCE for the complete licence information (https://github.com/csiro-crop-informatics/biokanga/LICENSE)
 * Contact: Alex Whan <alex.whan@csiro.au>
 */
// Extracts fasta sequences from a multifasta file
// Sequences to be extracted are identified by their descriptors

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


const int cMaxExtractDescrs = 20;					// at most this number of extraction descriptors can be processed

int
fastaextractproc(int PMode,				// processing mode - 0 by matching descriptor, 1 subsample by selecting every Nth sequence
		int NthSel,						// if > 0 then select every Nth descriptor.sequence to be extracted
		int XSense,						// extraction sense - 0 - original sense, 1 - reverse only, 2 - complement only 3 - reverse complement
		int NumDescrs,					// number of descriptors
		char *pszDescrs[],				// extract sequences with these descriptors
		char *pszInFasta,				// extract from this multifasta file
        char *pszOutFasta);				// output extracted sequences to this multifasta file



#ifdef _WIN32
int fastaextract(int argc, char* argv[])
{
// determine my process name
_splitpath(argv[0],NULL,NULL,gszProcName,NULL);
#else
int 
fastaextract(int argc, char** argv)
{
// determine my process name
CUtility::splitpath((char *)argv[0],NULL,gszProcName);
#endif

int iScreenLogLevel;		// level of screen diagnostics
int iFileLogLevel;			// level of file diagnostics
char szLogFile[_MAX_PATH];	// write diagnostics to this file
int Rslt = 0;   			// function result code >= 0 represents success, < 0 on failure

int PMode;					// processing mode
int NthSel;					// select every Nth sequence
int XSense;					// extraction sense - 0 - original sense, 1 - reverse only, 2 - complement only 3 - reverse complement

char szInFile[_MAX_PATH];		// read from this file
char szOutFile[_MAX_PATH];		// write extracted sequences to this file

int NumExtractDescrs;
char *pszExtractDescrs[cMaxExtractDescrs+1];

char szSQLiteDatabase[_MAX_PATH];	// results summaries to this SQLite file
char szExperimentName[cMaxDatasetSpeciesChrom+1];			// experiment name
char szExperimentDescr[1000];		// describes experiment

// command line args
struct arg_lit  *help    = arg_lit0("h","help",                 "Print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"Print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to screen and logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"Diagnostics log file");

struct arg_int *mode = arg_int0("m","mode","<int>",				"Processing mode: 0 extract by descriptor, 1 extract every Nth sequence");
struct arg_int *nthsel = arg_int0("n","nth","<int>",			"Extract every Nth descriptor 1 .. Nth");
struct arg_str *extractdescrs = arg_strn("e","extract","<str>",0,cMaxExtractDescrs,		"Extract sequences with these descriptors (wildcards allowed)");
struct arg_file *inputfile = arg_file1("i","in","<file>",		"Input file containing sequences to extract");
struct arg_file *outfile = arg_file1("o","out","<file>",		"Output extracted sequences to this file");

struct arg_int *xsense = arg_int0("x","xsense","<int>",           "Extracted sequences as: 0 - original sense, 1 - reverse only, 2 - complement only 3 - reverse complement");

struct arg_file *summrslts = arg_file0("q","sumrslts","<file>",					"Output results summary to this SQLite3 database file");
struct arg_str *experimentname = arg_str0("w","experimentname","<str>",			"experiment name SQLite3 database file");
struct arg_str *experimentdescr = arg_str0("W","experimentdescr","<str>",		"experiment description SQLite3 database file");

struct arg_end *end = arg_end(200);

void *argtable[] = {help,version,FileLogLevel,LogFile,
					summrslts,experimentname,experimentdescr,
	                mode,nthsel,extractdescrs, xsense, inputfile, outfile,
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
		printf("\n      To invoke this parameter file then precede its name with '@'");
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
		gProcessID = gSQLiteSummaries.AddProcess((char *)gszProcName,(char *)gszProcName,(char *)szExperimentDescr);
		if(gProcessID < 1)
			return(1);
		gProcessingID = gSQLiteSummaries.StartProcessing(gExperimentID,gProcessID,(char *)cpszProgVer);
		if(gProcessingID < 1)
			return(1);
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Initialised SQLite database '%s' for extraction results summary",szSQLiteDatabase);
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"SQLite database experiment identifier for '%s' is %d",szExperimentName,gExperimentID);
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"SQLite database process identifier for '%s' is %d",(char *)gszProcName,gProcessID);
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"SQLite database processing instance identifier is %d",gProcessingID);
		}
	else
		{
		szSQLiteDatabase[0] = '\0';
		szExperimentDescr[0] = '\0';
		}

	PMode = mode->count ? mode->ival[0] : 0;
	if(PMode < 0 || PMode > 2)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Processing mode '-m%d' must be in range 0..1",PMode);
		return(1);
		}

	
	if(PMode == 1)
		{
		NthSel = nthsel->count ? nthsel->ival[0] : 0;
		if(NthSel <= 0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Sample every Nth sequence requested but '-n<th>' not specified as %d",NthSel);
			return(1);
			}
		}
	else
		NthSel = 0;

	if(PMode == 0 && extractdescrs->count)
		{
		int Idx;
		int LenExtractDescr;
		NumExtractDescrs = 0;
		for(Idx=0;Idx < extractdescrs->count; Idx++)
			{
			pszExtractDescrs[Idx] = NULL;
			LenExtractDescr = (int)strlen(extractdescrs->sval[Idx]);
			if(pszExtractDescrs[NumExtractDescrs] == NULL)
				pszExtractDescrs[Idx] = new char [LenExtractDescr+1];
			strcpy(pszExtractDescrs[NumExtractDescrs],extractdescrs->sval[Idx]);
			CUtility::TrimQuotes(pszExtractDescrs[NumExtractDescrs]);
			if(pszExtractDescrs[NumExtractDescrs][0] != '\0')
				NumExtractDescrs++;
			}

		if(!NumExtractDescrs)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: After removal of whitespace, descriptors specified with '-e<descriptor>' option)\n");
			exit(1);
			}
		}
	else
		{
		NumExtractDescrs = 0;
		pszExtractDescrs[0] = NULL;
		}

	XSense = xsense->count ? xsense->ival[0] : 0;
	if(XSense < 0 || XSense > 3)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Extraction sense '-x%d' must be in range 0..3",XSense);
		return(1);
		}
	
	strcpy(szInFile,inputfile->filename[0]);
	CUtility::TrimQuotedWhitespcExtd(szInFile);

	strcpy(szOutFile,outfile->filename[0]);
	CUtility::TrimQuotedWhitespcExtd(szOutFile);

// show user current resource limits
#ifndef _WIN32
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Resources: %s",CUtility::ReportResourceLimits());
#endif

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing parameters:");
	const char *pszDescr;
	switch(PMode) {
		case 0:
			pszDescr = "Extract fasta sequence(s)";
			break;
		case 1:
			pszDescr = "Sample every Nth sequence";
			break;
		}

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Processing mode is : '%s'",pszDescr);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Experiment name : '%s'",szExperimentName);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Experiment description : '%s'",szExperimentDescr);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"SQLite database file: '%s'",szSQLiteDatabase);

	if(PMode == 0)
		{
		for(int Idx = 0; Idx < NumExtractDescrs; Idx++)
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Descriptor (%d) : '%s'",Idx+1,pszExtractDescrs[Idx]);
		}
	else
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Sample every Nth sequence: %d",NthSel);

	switch(XSense) {
			case 0:
				pszDescr = "original";
				break; 
			case 1:
				pszDescr = "reverse";
				break; 
			case 2:
				pszDescr = "complement";
				break; 
			case 3:
				pszDescr = "reverse complement";
				break;
			}
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Sequences extracted will be in : %s sense",pszDescr);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Extract sequences from file : '%s'",szInFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Write extracted features to file: '%s'",szOutFile);

	if(gExperimentID > 0)
		{
		int ParamID;
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szLogFile),"log",szLogFile);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,sizeof(PMode),"mode",&PMode);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,sizeof(XSense),"xsense",&XSense);

	    ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,sizeof(NumExtractDescrs),"NumExtractDescrs",&NumExtractDescrs);
		if(PMode == 0)
			{
			for(int Idx=0; Idx < NumExtractDescrs; Idx++)
				ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(pszExtractDescrs[Idx]),"extractdescrs",pszExtractDescrs[Idx]);
			}
		else
			ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,sizeof(NthSel),"nthsel",&NthSel);

		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szOutFile),"out",szOutFile);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szInFile),"in",szInFile);

		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szSQLiteDatabase),"sumrslts",szSQLiteDatabase);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szExperimentName),"experimentname",szExperimentName);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szExperimentDescr),"experimentdescr",szExperimentDescr);
		}


#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	gStopWatch.Start();

	Rslt = fastaextractproc(PMode,NthSel,XSense,NumExtractDescrs,pszExtractDescrs,szInFile,szOutFile);
	Rslt = Rslt >=0 ? 0 : 1;
	if(gExperimentID > 0)
		{
		if(gProcessingID)
			gSQLiteSummaries.EndProcessing(gProcessingID,Rslt);
		gSQLiteSummaries.EndExperiment(gExperimentID);
		}
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


return 0;
}

const size_t cMaxAllocBuffChunk = 0x0ffffff;		// allocate for input sequences buffering in these sized chunks

int m_hFEInFile;	// file handle for opened multifasta input file
int m_hFEOutFile;	// file handle for opened multifasta output file

size_t m_SeqBuffLen;		// number of bases currently buffered in m_pSeqBuff
size_t m_AllocdSeqBuffMem;  // size of memory currently allocated to m_pSeqBuff
UINT8 *m_pSeqBuff;			// buffers sequences as read from file

int m_NumExtractDescrs;		// number of extract descriptors
#ifdef _WIN32
Regexp *m_ExtractDescrsRE[cMaxExtractDescrs];	// compiled regular expressions
#else
regex_t m_ExtractDescrsRE[cMaxExtractDescrs];	// compiled regular expressions
#endif

void
FEReset(void)
{
if(m_hFEInFile != -1)
	{
	close(m_hFEInFile);
	m_hFEInFile = -1;
	}
if(m_hFEOutFile != -1)
	{
#ifdef _WIN32
	_commit(m_hFEOutFile);
#else
	fsync(m_hFEOutFile);
#endif
	close(m_hFEOutFile);
	m_hFEOutFile = -1;
	}

if(m_pSeqBuff != NULL)
	{
#ifdef _WIN32
	free(m_pSeqBuff);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pSeqBuff != MAP_FAILED)
		munmap(m_pSeqBuff,m_AllocdSeqBuffMem);
#endif
	m_pSeqBuff = NULL;
	}

m_AllocdSeqBuffMem = 0;
m_SeqBuffLen = 0;

m_NumExtractDescrs = 0;
}

void
FEInit(void)
{
m_hFEInFile = -1;
m_hFEOutFile = -1;
m_pSeqBuff = NULL;
FEReset();
}

etSeqBase *
AllocSeqBuff(size_t SeqLen)				// allocate for at least this sequence length
{
size_t memreq;
etSeqBase *pTmp;

if(m_pSeqBuff != NULL && m_AllocdSeqBuffMem >= SeqLen)
	return(m_pSeqBuff);

if(m_pSeqBuff == NULL)
	{
	memreq = max(SeqLen,(size_t)cMaxAllocBuffChunk);
#ifdef _WIN32
	m_pSeqBuff = (etSeqBase *) malloc(SeqLen);	// initial and perhaps the only allocation
	if(m_pSeqBuff == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AllocSeqBuff: Memory allocation of %lld bytes - %s",(INT64)SeqLen,strerror(errno));
		return(NULL);
		}
#else
	// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
	m_pSeqBuff = (etSeqBase *)mmap(NULL,memreq, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
	if(m_pSeqBuff == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AllocSeqBuff: Memory allocation of %lld bytes through mmap()  failed - %s",(INT64)memreq,strerror(errno));
		m_pSeqBuff = NULL;
		return(NULL);
		}
#endif
	}
else
	{
	memreq = SeqLen + cMaxAllocBuffChunk;
#ifdef _WIN32
	pTmp = (etSeqBase *) realloc(m_pSeqBuff,memreq);
#else
	pTmp = (etSeqBase *)mremap(m_pSeqBuff,m_AllocdSeqBuffMem,memreq,MREMAP_MAYMOVE);
	if(pTmp == MAP_FAILED)
		pTmp = NULL;
#endif
	if(pTmp == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AllocSeqBuff: Memory re-allocation to %lld bytes - %s",(INT64)memreq,strerror(errno));
		return(NULL);
		}
	m_pSeqBuff = pTmp;
	}
m_AllocdSeqBuffMem = memreq;
return(m_pSeqBuff);
}


int
CompileExtractRegExprs(int	NumExtractDescrs,	// number of regular expressions
		char **ppszExtractDescrs)		// array of regular expressions
{
int Idx;
int Len;
char szExtract[500];

#ifndef _WIN32
int RegErr;				// regular expression parsing error
char szRegErr[128];			// to hold RegErr as textual representation ==> regerror();
#endif

#ifdef _WIN32
try {
	for(Idx=0;Idx < NumExtractDescrs;Idx++)
		{
		m_ExtractDescrsRE[Idx] = new Regexp();
		strcpy(szExtract,ppszExtractDescrs[Idx]);
		Len = (int)strlen(szExtract);
		if(szExtract[Len-1] != '$')
			{
			szExtract[Len] = '$';
			szExtract[Len+1] = '\0';
			}
		m_ExtractDescrsRE[Idx]->Parse(szExtract,false);	// note case insensitive
		}
	}
catch(...)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to process include regexpr chrom '%s'",szExtract);
	return(eBSFerrMem);
	}
#else
for(Idx=0;Idx < NumExtractDescrs;Idx++)
	{
	strcpy(szExtract,ppszExtractDescrs[Idx]);
	Len = strlen(szExtract);
	if(szExtract[Len-1] != '$')
		{
		szExtract[Len] = '$';
		szExtract[Len+1] = '\0';
		}
	RegErr=regcomp(&m_ExtractDescrsRE[Idx],szExtract,REG_EXTENDED | REG_ICASE);	// note case insensitive
	if(RegErr)
		{
		regerror(RegErr,&m_ExtractDescrsRE[Idx],szRegErr,sizeof(szRegErr));
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to process include regexpr '%s' error: %s",szExtract,szRegErr);
		return(eBSFerrMem);
		}
	}
#endif
m_NumExtractDescrs = NumExtractDescrs;
return(eBSFSuccess);
}

bool									// returns true if descriptor accepted for extraction
ExtractDescrMatches(char *pszSeqDescr)	// descriptor as parsed from input multifasta
{
char szDescrIdent[200];
char Chr;
bool bAcceptDescr = false;
int Idx;

#ifdef _WIN32
RegexpMatch mc;
#else
regmatch_t mc;
int RegErr;					// regular expression parsing error
char szRegErr[128];			// to hold RegErr as textual representation ==> regerror();
#endif

for(Idx = 0; Idx < (sizeof(szDescrIdent) - 1); Idx++)
	{
	switch(Chr = *pszSeqDescr++) {
		case '\0': case ' ': case '\t':
			break;
		default:
			szDescrIdent[Idx] = Chr;
			continue;
		}
	break;
	}
szDescrIdent[Idx] = '\0';

for(Idx = 0; Idx < m_NumExtractDescrs; Idx++)
		{
#ifdef _WIN32
		if(m_ExtractDescrsRE[Idx]->Match(szDescrIdent,&mc))
#else
		if(!regexec(&m_ExtractDescrsRE[Idx],szDescrIdent,1,&mc,0))
#endif
			{
			bAcceptDescr = true;
			break;
			}
		}

return(bAcceptDescr);
}

int
WriteSeqFile(int XSense,						// extraction sense - 0 - original sense, 1 - reverse only, 2 - complement only 3 - reverse complement
			char *pszDescr,int SeqLen,etSeqBase *pSeq)
{
char Chr;
int LineLen;
int BuffIdx;
int SeqIdx;
UINT8 szOutBuff[0x07fff];

switch(XSense) {
	case 0:		// original sense
		break;
	case 1:     // reverse only
		CSeqTrans::ReverseSeq(SeqLen,pSeq);
		break;
	case 2:		// complement only
		CSeqTrans::ComplementStrand(SeqLen,pSeq);
		break;
	case 3:		// reverse and complement
		CSeqTrans::ReverseComplement(SeqLen,pSeq);
		break;
	}

LineLen = 0;
		// write out descriptor to file
BuffIdx = sprintf((char *)szOutBuff,">%s\n",pszDescr);
for(SeqIdx = 0; SeqIdx < SeqLen; SeqIdx++,pSeq++,LineLen++)
	{
	if(LineLen > 75)
		{
		szOutBuff[BuffIdx++] = '\n';
		LineLen = 0;
		}
	switch(*pSeq & 0x07) {
		case eBaseA:
			Chr = 'A';
			break;
		case eBaseC:
			Chr = 'C';
			break;
		case eBaseG:
			Chr = 'G';
			break;
		case eBaseT:
			Chr = 'T';
			break;
		default:
			Chr = 'N';
			break;
		}
	szOutBuff[BuffIdx++] = Chr;
	if(BuffIdx > (sizeof(szOutBuff) - 100))
		{
		CUtility::SafeWrite(m_hFEOutFile,szOutBuff,BuffIdx);
		BuffIdx = 0;
		}
	}
szOutBuff[BuffIdx++] = '\n';
CUtility::SafeWrite(m_hFEOutFile,szOutBuff,BuffIdx);
return(SeqLen);
}

int
fastaextractproc(int PMode,				// processing mode - 0 by matching descriptor, 1 subsample by selecting every Nth sequence
		int NthSel,						// if > 0 then select every Nth descriptor.sequence to be extracted
		int XSense,						// extraction sense - 0 - original sense, 1 - reverse only, 2 - complement only 3 - reverse complement
		int NumDescrs,					// number of descriptor reqular expressions (in pszDescrs[]) to match on
		char *pszDescrs[],				// extract sequences with these descriptors
		char *pszInFasta,				// extract from this multifasta file
        char *pszOutFasta)				// output extracted sequences to this multifasta file
{
int Rslt;
int NumMatches;
int SeqID;
int Descrlen;
int SeqLen;
size_t AvailBuffSize;
bool bExtract;
char szInDescription[cBSFDescriptionSize];

CFasta Fasta;

FEInit();

if(PMode == 0 && NumDescrs >= 1)
	{
	if((Rslt = CompileExtractRegExprs(NumDescrs,pszDescrs))!=eBSFSuccess)
		return(Rslt);
	}

#ifdef _WIN32
if((m_hFEOutFile = open(pszOutFasta, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
if((m_hFEOutFile = open(pszOutFasta,O_RDWR | O_CREAT |O_TRUNC, S_IREAD | S_IWRITE))==-1)
#endif
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create or truncate %s - %s",pszOutFasta,strerror(errno));
	FEReset();
	return(eBSFerrCreateFile);
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Output results file created/truncated: '%s'",pszOutFasta);

if((Rslt=Fasta.Open(pszInFasta,true))!=eBSFSuccess)
	{
	if(Rslt != eBSFerrNotFasta)
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessFastaFile: Unable to open '%s' [%s] %s",pszInFasta,Fasta.ErrText((teBSFrsltCodes)Rslt),Fasta.GetErrMsg());
	return(Rslt);
	}

if(m_pSeqBuff == NULL)				// if not already allocated then allocate to hold cMaxAllocBuffChunk bases 
	{
	SeqLen = cMaxAllocBuffChunk;
	if(AllocSeqBuff(SeqLen) == NULL)
		{
		Rslt = eBSFerrMem;
		Fasta.Close();
		return(Rslt);
		}
	}

NumMatches = 0;
SeqID = 0;
m_SeqBuffLen = 0;
AvailBuffSize = m_AllocdSeqBuffMem;
bExtract = false;
while((Rslt = SeqLen = Fasta.ReadSequence(&m_pSeqBuff[m_SeqBuffLen],(int)min(AvailBuffSize,(size_t)cMaxAllocBuffChunk),true,false)) > eBSFSuccess)
	{
	if(SeqLen == eBSFFastaDescr)		// just read a descriptor line
		{
		if(bExtract && m_SeqBuffLen)
			WriteSeqFile(XSense,szInDescription,(int)m_SeqBuffLen,m_pSeqBuff);
		else
			if(bExtract)
				gDiagnostics.DiagOut(eDLWarn,gszProcName,"Matching descriptor '%s' but no sequence in file - '%s'",szInDescription,pszOutFasta);
		bExtract = false;
		m_SeqBuffLen = 0;
		AvailBuffSize = m_AllocdSeqBuffMem;
		bExtract = false;
		SeqID++;

		Descrlen = Fasta.ReadDescriptor(szInDescription,cBSFDescriptionSize-10);

		// simply subsampling by selecting every Nth sequence?
		if(NumDescrs == 0 || PMode == 1)
			{
			if(PMode != 0 && NthSel > 1 && (SeqID % NthSel))
				continue;
			bExtract = true;
			}
		else	// is this a descriptor of interest?
			if(!ExtractDescrMatches(szInDescription))
				continue;
		NumMatches += 1;
		bExtract = true;				// flag that the sequence is to be extracted
		m_SeqBuffLen = 0;
		AvailBuffSize = m_AllocdSeqBuffMem;
		continue;
		}

	if(!bExtract)						// slough this sequence?
		continue;

	// write to file
	m_SeqBuffLen += SeqLen;
	AvailBuffSize -= SeqLen;
	if(AvailBuffSize < (size_t)(cMaxAllocBuffChunk / 8))
		{
		if(AllocSeqBuff(m_AllocdSeqBuffMem + SeqLen) == NULL)
			{
			Rslt = eBSFerrMem;
			Fasta.Close();
			return(Rslt);
			}
		AvailBuffSize = m_AllocdSeqBuffMem - m_SeqBuffLen;
		}
	}
if(Rslt < eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Encountered parse errors");
	while(Fasta.NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,Fasta.GetErrMsg());
	}
if(bExtract && m_SeqBuffLen)					// was previous sequence to be extracted?
	WriteSeqFile(XSense,szInDescription,(int)m_SeqBuffLen,m_pSeqBuff);

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Extracted %d sequences",NumMatches);
FEReset();
return(Rslt);
}