// PEscaffold.cpp : Defines the entry point for the console application.
// This process will attempt to scaffold contigs using paired end reads
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

#include "SQLiteSummaries.h"

const char *cpszProgVer = "0.5.0";		// increment with each release

const int cMaxChromNameLen = 50;		// handle chroms with names up to this length
const int cMaxPENameLen = 50;			// handle paired end names up to this length

const int cMinSeqLength = 10;			// minimal accepted sequence length
const int cMaxSeqLength = 2000;			// maximal accepted sequence length

const int cAllocChromNames	= 250000;	// alloc/realloc chromosome names in increments of this many
const int cAllocPENames		= 25000000;	// alloc/realloc paired end identifier names in increments of this many
const int cAllocScafolds	= 25000000;	// alloc/realloc scaffolds in increments of this many

const int cHashSize = 0x0ffffff;			// use this sized hash (24bits) 

#pragma pack(1)
typedef struct TAG_sPEScaffoldChrom {
	INT32 ChromID;						// uniquely identifies this chrom/contig
	char szChrom[cMaxChromNameLen+1];	// which has this name
	INT32 HashNext;						// if non-zero then identifier of next chrom with same hash
	} tsPEScaffoldChrom;

typedef struct TAG_sPEIdent {
	INT32 IdentID;						// uniquely identifies this PE identifier
	char szIdent[cMaxPENameLen+1];		// which has this name
	INT32 PEScafoldID;					// associated scaffold (invalid after scaffolds sorted)
	INT32 HashNext;						// if non-zero then identifier of next PEIdent with same hash
	} tsPEIdent;

typedef struct TAG_sPEScaffold {
	INT32 PEScafoldID;					// uniquely identifies this scaffold
	INT32 PE12SeqID;					// PE1/PE2 sequence identifier
	INT32 PE1ChromID;					// PE1 is aligned onto this chromosome/contig (0 if alignment unknown)
	INT32 PE2ChromID;					// PE2 is aligned onto this chromosome/contig (0 if alignment unknown)
	UINT8 PE1Sense:1;					// 1 if PE1 aligned sense onto PE1ChromID; 0 if aligned antisense
	UINT8 PE2Sense:1;					// 1 if PE2 aligned sense onto PE2ChromID; 0 if aligned antisense
} tsPEScaffold; 
#pragma pack()

int
Process(int PMode,					// processing mode
		char *pszInPE1File,			// input PE1 file
		char *pszInPE2File,			// input PE2 file
		char *pszOutFile);			// output corelations file


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
	return _T("CSVFilter");
}
// end of str library required code
#endif

CStopWatch gStopWatch;
CDiagnostics gDiagnostics;				// for writing diagnostics messages to log file
char gszProcName[_MAX_FNAME];			// process name

CSQLiteSummaries gSQLiteSummaries;		// for writing processing result summaries to SQLite database
int	gExperimentID = 0;					// SQLite experiment identifier
int gProcessID = 0;						// SQLite process identifier
int	gProcessingID = 0;					// SQLite processing identifier

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

int PMode;					// processing mode

char szInPE1File[_MAX_PATH];	// parse PE1 alignments from this file
char szInPE2File[_MAX_PATH];	// parse PE1 alignments from this file

char szOutFile[_MAX_PATH];		// write corelations to this file

char szSQLiteDatabase[_MAX_PATH];	// results summaries to this SQLite file
char szExperimentName[cMaxDatasetSpeciesChrom+1];			// experiment name
char szExperimentDescr[1000];		// describes experiment

// command line args
struct arg_lit  *help    = arg_lit0("h","help",                 "Print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"Print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to screen and logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"Diagnostics log file");

struct arg_int *mode = arg_int0("m","mode","<int>",				"Processing mode: 0 default");

struct arg_file *inpe1file = arg_file1("i","in","<file>",		"Input SAM file containing PE1 alignments");
struct arg_file *inpe2file = arg_file1("I","in","<file>",		"Input SAM file containing PE2 alignments");

struct arg_file *outfile = arg_file1("o","out","<file>",		"Output corelations to this file");

struct arg_file *summrslts = arg_file0("q","sumrslts","<file>",		"Output results summary to this SQLite3 database file");
struct arg_str *experimentname = arg_str0("w","experimentname","<str>",		"experiment name SQLite3 database file");
struct arg_str *experimentdescr = arg_str0("W","experimentdescr","<str>",	"experiment description SQLite3 database file");

struct arg_end *end = arg_end(20);

void *argtable[] = {help,version,FileLogLevel,LogFile,
					summrslts,experimentname,experimentdescr,
	                mode,inpe1file,inpe2file, outfile,
					end};

char **pAllArgs;
int argerrors;
argerrors = CUtility::arg_parsefromfile(argc,(char **)argv,&pAllArgs);
if(argerrors >= 0)
	argerrors = arg_parse(argerrors,pAllArgs,argtable);

/* special case: '--help' takes precedence over error reporting */
if (help->count > 0)
        {
		printf("\n%s PE1/PE2 Corelation, Version %s\nOptions ---\n", gszProcName,cpszProgVer);
        arg_print_syntax(stdout,argtable,"\n");
        arg_print_glossary(stdout,argtable,"  %-25s %s\n");
		printf("\nNote: Parameters can be entered into a parameter file, one parameter per line.");
		printf("\n      To invoke this parameter file then precede its name with '@'");
		printf("\n      e.g. %s @myparams.txt\n",gszProcName);
		printf("\nPlease report any issues regarding usage of %s at https://github.com/csiro-crop-informatics/biokanga/issues\n\n",gszProcName);
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
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Process %s Version %s starting",gszProcName,cpszProgVer);

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
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Initialised SQLite database '%s' for results summary collection",szSQLiteDatabase);
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
	if(PMode < 0 || PMode > 1)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Processing mode '-m%d' must be in range 0..0",PMode);
		return(1);
		}

	strcpy(szInPE1File,inpe1file->filename[0]);
	CUtility::TrimQuotedWhitespcExtd(szInPE1File);

	strcpy(szInPE2File,inpe2file->filename[0]);
	CUtility::TrimQuotedWhitespcExtd(szInPE2File);

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
			pszDescr = "Corelate PE1 and PE2 alignments";
			break;
		}

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Processing mode is : '%s'",pszDescr);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Experiment name : '%s'",szExperimentName);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Experiment description : '%s'",szExperimentDescr);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"SQLite database file: '%s'",szOutFile);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"PE1 SAM file: '%s'",szInPE1File);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"PE2 SAM file: '%s'",szInPE2File);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Corelations to file: '%s'",szOutFile);

#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	gStopWatch.Start();

	Rslt = Process(PMode,szInPE1File,szInPE2File,szOutFile);
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
    printf("\n%s, Version %s\n", gszProcName,cpszProgVer);
	arg_print_errors(stdout,end,gszProcName);
	arg_print_syntax(stdout,argtable,"\nUse '-h' to view option and parameter usage\n");
	exit(1);
	}
return 0;
}

// Following should one day, real soon, be made into a class

CMTqsort mtqsort;						// muti-threaded qsort
static int SortScaffolds(const void *arg1, const void *arg2);

int m_AllocdNumScaffoldChroms;			// current allocation can hold at most this many chroms
int m_NumScaffoldChroms;				// this many scaffold chroms are currently in m_pScaffoldChroms
size_t m_AllocdScaffoldChromsMem;		// m_pScaffoldChroms current memory allocation size
tsPEScaffoldChrom *m_pScaffoldChroms;	// allocated to hold scaffold chromosome names
int *m_pHashChroms;						// holds hashes mapping to chrom identifers 

int m_AllocdNumPEIdents;				// current allocation can hold at most this many PE identifiers
int m_NumPEIdents;						// this many PE identifiers are currently in m_pNumPEIdents
size_t m_AllocdPEIdentsMem;				// m_pNumPEIdents current memory allocation size
tsPEIdent *m_pPEIdents;					// allocated to hold PE identifier names
int *m_pHashPEIdents;					// holds hashes mapping to PEIdent identifers 

int m_AllocdNumScaffolds;				// current allocation can hold at most this many scaffolds
int m_NumScaffolds;						// this many scaffolds are currently in m_pScaffoldChroms
size_t m_AllocdScaffoldsMem;			// m_pScaffolds current memory allocation size
tsPEScaffold *m_pScaffolds;				// allocated to hold scaffolds

int m_hOutFile;							// corelations to this file

void Reset(void)
{
if(m_hOutFile != -1)
	{
	close(m_hOutFile);
	m_hOutFile = -1;
	}
if(m_pPEIdents != NULL)
	{
#ifdef _WIN32
	free(m_pPEIdents);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pPEIdents != MAP_FAILED)
		munmap(m_pPEIdents,m_AllocdPEIdentsMem);
#endif
	m_pPEIdents = NULL;
	}

if(m_pHashPEIdents)
	{
	delete m_pHashPEIdents;
	m_pHashPEIdents = NULL;
	}

if(m_pScaffoldChroms != NULL)
	{
#ifdef _WIN32
	free(m_pScaffoldChroms);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pScaffoldChroms != MAP_FAILED)
		munmap(m_pScaffoldChroms,m_AllocdScaffoldChromsMem);
#endif
	m_pScaffoldChroms = NULL;
	}

if(m_pHashChroms)
	{
	delete m_pHashChroms;
	m_pHashChroms = NULL;
	}

if(m_pScaffolds != NULL)
	{
#ifdef _WIN32
	free(m_pScaffolds);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pScaffolds != MAP_FAILED)
		munmap(m_pScaffolds,m_AllocdScaffoldsMem);
#endif
	m_pScaffolds = NULL;
	}

m_AllocdNumPEIdents = 0;
m_NumPEIdents = 0;
m_AllocdPEIdentsMem = 0;

m_AllocdNumScaffoldChroms = 0;
m_AllocdScaffoldChromsMem = 0;
m_NumScaffoldChroms = 0;

m_AllocdScaffoldsMem = 0;
m_AllocdNumScaffolds = 0;
m_NumScaffolds = 0;
}

void Init(void)
{
m_pPEIdents = NULL;
m_pHashPEIdents = NULL;
m_pScaffoldChroms = NULL;
m_pHashChroms = NULL;
m_pScaffolds = NULL;
m_hOutFile = -1;
Reset();
}

char *
TrimWhitespace(char *pTxt)
{
char *pStart;
char Chr;
	// strip leading whitespace
while(Chr = *pTxt++)
	if(!isspace(Chr))
			break;
if(Chr == '\0')					// empty line?
	return(pTxt-1);
pStart = pTxt-1;
while(Chr = *pTxt)			// fast forward to line terminator
	pTxt++;
pTxt-=1;
while(Chr = *pTxt--)
	if(!isspace(Chr))
		break;
pTxt[2] = '\0';
return(pStart);
}

// generate a cHashSize hash over pszName
int						// returned hash - hash will be in the range 0..cHashSize
GenNameHash(char *pszName)
{
int Hash;
char Chr;
Hash = 37199;			// circular prime as hash seed

while(Chr = *pszName++)
	{
	Hash = (Hash ^ (int)Chr) * 19937;	// a circular prime
	Hash ^= (Hash >> 21);
	Hash &= cHashSize;
	}
return(Hash & cHashSize);
}

// AddChromName
// If chrom already known then return existing chrom identifier otherwise add to m_pScaffoldChroms
INT32
AddChromName(char *pszChromName)
{
static INT32 PrevChromID = 0;
int Hash;
int HashIdx;
tsPEScaffoldChrom *pScaffoldChrom;

// check to see if chromosome name already known
if(PrevChromID != 0)
	{
	pScaffoldChrom = &m_pScaffoldChroms[PrevChromID-1];
	if(!stricmp(pszChromName,pScaffoldChrom->szChrom))
		return(pScaffoldChrom->ChromID);
	}

Hash = GenNameHash(pszChromName);
if((HashIdx = m_pHashChroms[Hash]) != 0)
	{
	do {
		pScaffoldChrom = &m_pScaffoldChroms[HashIdx-1];
		if(!stricmp(pszChromName,pScaffoldChrom->szChrom))
			{
			PrevChromID = pScaffoldChrom->ChromID;
			return(pScaffoldChrom->ChromID);
			}
		HashIdx = pScaffoldChrom->HashNext;
		}
	while(HashIdx > 0);
	}

// its a new chrom not previously seen
// realloc as may be required to hold this new chrom
if(m_NumScaffoldChroms == m_AllocdNumScaffoldChroms)
	{
	size_t memreq = m_AllocdScaffoldChromsMem + (cAllocChromNames * sizeof(tsPEScaffoldChrom));
#ifdef _WIN32
	pScaffoldChrom = (tsPEScaffoldChrom *) realloc(m_pScaffoldChroms,memreq);
#else
	pScaffoldChrom = (tsPEScaffoldChrom *)mremap(m_pScaffoldChroms,m_AllocdScaffoldChromsMem,memreq,MREMAP_MAYMOVE);
	if(pScaffoldChrom == MAP_FAILED)
		pScaffoldChrom = NULL;
#endif
	if(pScaffoldChrom == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddChromName: Memory re-allocation to %d bytes - %s",memreq,strerror(errno));
		return(eBSFerrMem);
		}
	m_pScaffoldChroms = pScaffoldChrom;
	m_AllocdScaffoldChromsMem = memreq;
	m_AllocdNumScaffoldChroms += cAllocChromNames;
	}

pScaffoldChrom = &m_pScaffoldChroms[m_NumScaffoldChroms++];
pScaffoldChrom->ChromID = m_NumScaffoldChroms;
pScaffoldChrom->HashNext = m_pHashChroms[Hash];
m_pHashChroms[Hash] = m_NumScaffoldChroms;
strncpy(pScaffoldChrom->szChrom,pszChromName,cMaxChromNameLen);
pScaffoldChrom->szChrom[cMaxChromNameLen] = '\0';
PrevChromID = m_NumScaffoldChroms;
return(m_NumScaffoldChroms);
}

char *
GetChromName(int ChromID)
{
tsPEScaffoldChrom *pScaffoldChrom;
pScaffoldChrom = &m_pScaffoldChroms[ChromID-1];
return(pScaffoldChrom->szChrom);
}

// AddPEIdent
// If PE name already known then return existing identifier otherwise add to m_pScaffoldChroms
INT32
AddPEIdent(char *pszIdentName)
{
static INT32 PrevPEIdentID = 0;
tsPEIdent *pPEIdent;
int Hash;
int HashIdx;

// check to see if PEIdent name already known
if(PrevPEIdentID != 0)
	{
	pPEIdent = &m_pPEIdents[PrevPEIdentID-1];
	if(!stricmp(pszIdentName,pPEIdent->szIdent))
		return(pPEIdent->IdentID);
	}

Hash = GenNameHash(pszIdentName);
if((HashIdx = m_pHashPEIdents[Hash]) != 0)
	{
	do {
		pPEIdent = &m_pPEIdents[HashIdx-1];
		if(!stricmp(pszIdentName,pPEIdent->szIdent))
			{
			PrevPEIdentID = pPEIdent->IdentID;
			return( pPEIdent->IdentID);
			}
		HashIdx = pPEIdent->HashNext;
		}
	while(HashIdx > 0);
	}

// its a new PE identifier not previously seen
// realloc as may be required to hold this new chrom
if(m_NumPEIdents == m_AllocdNumPEIdents)
	{
	size_t memreq = m_AllocdPEIdentsMem + (cAllocPENames * sizeof(tsPEIdent));
#ifdef _WIN32
	pPEIdent = (tsPEIdent *) realloc(m_pPEIdents,memreq);
#else
	pPEIdent = (tsPEIdent *)mremap(m_pPEIdents,m_AllocdPEIdentsMem,memreq,MREMAP_MAYMOVE);
	if(pPEIdent == MAP_FAILED)
		pPEIdent = NULL;
#endif
	if(pPEIdent == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddPEIdent: Memory re-allocation to %d bytes - %s",memreq,strerror(errno));
		return(eBSFerrMem);
		}
	m_pPEIdents = pPEIdent;
	m_AllocdPEIdentsMem = memreq;
	m_AllocdNumPEIdents += cAllocPENames;
	}

pPEIdent = &m_pPEIdents[m_NumPEIdents++];
pPEIdent->IdentID = m_NumPEIdents;
pPEIdent->PEScafoldID = 0;
pPEIdent->HashNext = m_pHashPEIdents[Hash];
m_pHashPEIdents[Hash] = m_NumPEIdents;
strncpy(pPEIdent->szIdent,pszIdentName,cMaxChromNameLen);
pPEIdent->szIdent[cMaxChromNameLen] = '\0';
PrevPEIdentID = m_NumPEIdents;
return(m_NumPEIdents);
}

char *
GetSeqName(int SeqID)
{
tsPEIdent *pIdent;
pIdent = &m_pPEIdents[SeqID-1];
return(pIdent->szIdent);
}

int
AddScaffold(bool bPE2,						// if false then PE1, if true then PE2
			char *pszPEIdent,				// paired end indentifier used to corelate paired ends
			char *pszChrom,					// PE aligns onto this chromosome
			char Strand)					// '+' or '-'
{
static int PrevScafoldID = 0;
tsPEScaffold *pPEScaffold;
tsPEIdent *pPEIdent;
int ChromID;
int PEIdentID;

ChromID = AddChromName(pszChrom);
PEIdentID = AddPEIdent(pszPEIdent);

// if processing PE2 then check to see if already have scaffold with same PEIdent
if(bPE2)
	{
	if(PrevScafoldID != 0)
		{
		pPEScaffold = &m_pScaffolds[PrevScafoldID-1];
		if(pPEScaffold->PE12SeqID == PEIdentID)
			{
			pPEScaffold->PE2ChromID = ChromID;
			pPEScaffold->PE2Sense = Strand == '+' ? 1 : 0;
			return(PrevScafoldID);
			}
		}

	pPEIdent = &m_pPEIdents[PEIdentID-1];
	if((PrevScafoldID = pPEIdent->PEScafoldID) != 0)
		{
		pPEScaffold = &m_pScaffolds[PrevScafoldID-1];
		pPEScaffold->PE2ChromID = ChromID;
		pPEScaffold->PE2Sense = Strand == '+' ? 1 : 0;
		return(PrevScafoldID);
		}
	}
else	// PE1 processing, normally expect that the PE1 identifiers are unique, here we are just confirming that they are unique
	{   // when satisfied PE1 identifiers are always unique then could skip this check...
	pPEIdent = &m_pPEIdents[PEIdentID-1];
	if((PrevScafoldID = pPEIdent->PEScafoldID) != 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddScaffold: duplicate PE1 identifer - %s onto %s",pszPEIdent,pszChrom);
		return(0);
		}
	}

// new scaffold required
// realloc as may be required to hold this new scaffold
if(m_NumScaffolds == m_AllocdNumScaffolds)
	{
	size_t memreq = m_AllocdScaffoldsMem + (cAllocScafolds * sizeof(tsPEScaffold));
#ifdef _WIN32
	pPEScaffold = (tsPEScaffold *) realloc(m_pScaffolds,memreq);
#else
	pPEScaffold = (tsPEScaffold *)mremap(m_pScaffolds,m_AllocdScaffoldsMem,memreq,MREMAP_MAYMOVE);
	if(pPEScaffold == MAP_FAILED)
		pPEScaffold = NULL;
#endif
	if(pPEScaffold == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddScaffold: Memory re-allocation to %d bytes - %s",memreq,strerror(errno));
		return(eBSFerrMem);
		}
	m_pScaffolds = pPEScaffold;
	m_AllocdScaffoldsMem = memreq;
	m_AllocdNumScaffolds += cAllocScafolds;
	}
pPEScaffold = &m_pScaffolds[m_NumScaffolds++];
memset(pPEScaffold,0,sizeof(tsPEScaffold));
pPEScaffold->PEScafoldID = m_NumScaffolds;
pPEIdent->PEScafoldID = m_NumScaffolds;
pPEScaffold->PE12SeqID = PEIdentID;
if(!bPE2)
	{
	pPEScaffold->PE1ChromID = ChromID;
	pPEScaffold->PE1Sense = Strand == '+' ? 1 : 0;
	}
else
	{
	pPEScaffold->PE2ChromID = ChromID;
	pPEScaffold->PE2Sense = Strand == '+' ? 1 : 0;
	}
return(m_NumScaffolds);
}

int
LoadSAM(bool bPE2,			// false if loading PE1, true if loading PE2
		char *pszSAMFile)	// load alignments from this SAM file
{
etClassifyFileType FileType;
FILE *pSAMStream;
int NumParsedElLines;
int NumAcceptedEls;
char szLine[16000];				// buffer input lines
char *pTxt;
char szDescriptor[128];			// parsed out descriptor
int Flags;						// parsed out flags
char szChrom[128];				// parsed out chrom
int StartLoci;					// start loci
int NumUnmappedEls;
int ScaffoldID;

	// open SAM for reading
if(pszSAMFile == NULL || *pszSAMFile == '\0')
	return(eBSFerrParams);

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Loading alignments for %s from: '%s'",bPE2 ? "PE2" : "PE1", pszSAMFile);

FileType = CUtility::ClassifyFileType(pszSAMFile);
if(FileType != eCFTSAM)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to classify file as SAM formated: '%s'",pszSAMFile);
	return(eBSFerrOpnFile);
	}

if((pSAMStream = fopen(pszSAMFile,"r"))==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"ParseSAMFileElements: Unable to fopen SAM format file %s error: %s",pszSAMFile,strerror(errno));
	return(eBSFerrOpnFile);
	}

NumParsedElLines = 0;
NumAcceptedEls = 0;
NumUnmappedEls = 0;
while(fgets(szLine,sizeof(szLine)-1,pSAMStream)!= NULL)
	{
	NumParsedElLines += 1;
	if(!(NumParsedElLines % 1000000) || NumParsedElLines == 1)
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Parsed %d SAM lines",NumParsedElLines);

	szLine[sizeof(szLine)-1] = '\0';
	pTxt = TrimWhitespace(szLine);
	if(*pTxt=='\0' || *pTxt=='@')	// simply slough lines which were just whitespace or start with '@'
		continue;
	
	// expecting to parse as "%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t", szDescriptor, Flags, m_szSAMTargChromName, StartLoci+1,MAPQ,szCigar,pszRNext,PNext,TLen);
	// interest is in the descriptor,chromname,flags
	sscanf(szLine,"%s\t%d\t%s\t%d\t",szDescriptor, &Flags, szChrom, &StartLoci);
		// check if element has been mapped, if not then slough ...
	if(Flags & 0x04)	// will be set if unmapped
		{
		NumUnmappedEls += 1;
	    continue;
		}

	if((ScaffoldID = AddScaffold(bPE2,szDescriptor,szChrom,Flags & 0x10 ? '+' : '-')) < 1)
		return(ScaffoldID);
	NumAcceptedEls += 1;
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Loading alignments (%d) for %s from: '%s' completed",NumAcceptedEls,bPE2 ? "PE2" : "PE1", pszSAMFile);
return(NumAcceptedEls);
}

// locate scaffold having matching PE1ChromID and PE2ChromID
// assumes scaffolds have been sorted PE1ChromID.PE2ChromID ascending
tsPEScaffold *
LocateMateScaffold(int PE1ChromID,int PE2ChromID)
{
int Mid;
int Hi;
int Lo;
tsPEScaffold *pEl1;
tsPEScaffold *pEl2;
Lo = 0;
Hi = m_NumScaffolds-1;
do {
	Mid = (Hi + Lo) / 2;
	pEl1 = &m_pScaffolds[Mid];
	if(pEl1->PE1ChromID == PE1ChromID)
		{
		if(pEl1->PE2ChromID == PE2ChromID)
			{
			// drill down to locate 1st instance - linear drill down as not expecting too many instances
			pEl2 = pEl1--;
			while(Mid-- && (pEl1->PE1ChromID == PE1ChromID && pEl1->PE2ChromID == PE2ChromID))
					pEl2 = pEl1--;
			return(pEl2);
			}
		if(pEl1->PE2ChromID > PE2ChromID)
			Hi = Mid - 1;
		else
			Lo = Mid + 1;
		}
	else
		{
		if(pEl1->PE1ChromID > PE1ChromID)
			Hi = Mid - 1;
		else
			Lo = Mid + 1;
		}
	}
while(Hi >= Lo);
return(NULL);		// unable to match any
}

int
ReportCorelationships(char *pszOutFile)
{
char *pszPE1Chrom;
char *pszPE2Chrom;
char *pszSeqID;
int PrevPE1ChromID;
int PrevPE2ChromID;
int NumPEAligned;
int NumSEAligned;
int NumSenseSense;
int NumSenseAnti;
int	RevNumSenseSense;
int	RevNumSenseAnti;
int NumUnpaired;
int NumPaired;
tsPEScaffold *pPEScaffold;
bool bRevMate;
tsPEScaffold *pMateScaffold;
int Idx;
int BuffIdx;
char szBuff[16000];

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Identifying and writing corelations to file: '%s'",pszOutFile);

PrevPE1ChromID = -1;
PrevPE2ChromID = -1;
NumPEAligned = 0;
NumSEAligned = 0;
NumSenseSense = 0;
NumSenseAnti = 0;
RevNumSenseSense = 0;
RevNumSenseAnti = 0;
NumUnpaired = 0;
NumPaired = 0;
BuffIdx = 0;
BuffIdx = sprintf(&szBuff[0],"\"PE1\",\"PE2\",\"NumAligned\",\"NumSenseSense\",\"NumSenseAnti\",\"Paired\",\"Self\",\"RevMate\",\"RevNumSenseSense\",\"RevNumSenseAnti\"\n");
pPEScaffold = m_pScaffolds;
for(Idx = 0; Idx < m_NumScaffolds; Idx++,pPEScaffold++)
	{
	if(pPEScaffold->PE1ChromID != PrevPE1ChromID || pPEScaffold->PE2ChromID != PrevPE2ChromID) // starting a different scaffold?
		{
		if(PrevPE1ChromID > 0 || PrevPE2ChromID > 0)
			{
			// write out any prev scaffold info here
			if(BuffIdx > (sizeof(szBuff)-500))
				{
				CUtility::SafeWrite(m_hOutFile,szBuff,BuffIdx);
				BuffIdx = 0;
				}
			if(PrevPE1ChromID > 0 && PrevPE2ChromID > 0)
				BuffIdx += sprintf(&szBuff[BuffIdx],"\"%s\",\"%s\",%d,%d,%d,\"Y\",\"%s\",\"%s\",%d,%d\n",
							pszPE1Chrom,pszPE2Chrom,NumPEAligned,NumSenseSense,NumSenseAnti,PrevPE1ChromID==PrevPE2ChromID ? "Y" : "N",bRevMate ? "Y" : "N",RevNumSenseSense,RevNumSenseAnti);
			else
				{
				if(PrevPE1ChromID > 0)
					BuffIdx += sprintf(&szBuff[BuffIdx],"\"%s\",\"%s\",%d,%d,%d,\"N\",\"N\",\"N\",0,0\n",pszPE1Chrom,"N/A",NumSEAligned,0,0);
				else
					BuffIdx += sprintf(&szBuff[BuffIdx],"\"%s\",\"%s\",%d,%d,%d,\"N\",\"N\",\"N\",0,0\n","N/A",pszPE2Chrom,NumSEAligned,0,0);
				}
			}

		NumPEAligned = 0;					// number aligning as paired ends
		NumSEAligned = 0;					// number aligning as single ended only
		NumSenseSense = 0;					// number aligning sense to sense
		NumSenseAnti = 0;					// number aligning sense to antisense (or antisense to sense)
		if(pPEScaffold->PE1ChromID > 0 && pPEScaffold->PE2ChromID > 0)
			{
			NumPEAligned = 1;
			if(pPEScaffold->PE1Sense == pPEScaffold->PE2Sense)
				NumSenseSense = 1;
			else
				NumSenseAnti = 1;
			NumPaired += 1;
			}
		else
			{
			NumSEAligned = 1;
			NumUnpaired += 1;
			}
		pszSeqID = GetSeqName(pPEScaffold->PE12SeqID);
		if(pPEScaffold->PE1ChromID)
			pszPE1Chrom = GetChromName(pPEScaffold->PE1ChromID);
		else
			pszPE1Chrom = (char *)"N/A";

		if(pPEScaffold->PE2ChromID)
			pszPE2Chrom = GetChromName(pPEScaffold->PE2ChromID);
		else
			pszPE2Chrom = (char *)"N/A";

		PrevPE1ChromID = pPEScaffold->PE1ChromID;
		PrevPE2ChromID = pPEScaffold->PE2ChromID;
	
		// if not to self then check if the PE2 chrom has any PE2 linking back to this PE1 and count linking sense/antisense
		bRevMate = false;
		RevNumSenseSense = 0;
		RevNumSenseAnti = 0;
		if(PrevPE1ChromID != PrevPE2ChromID && (PrevPE1ChromID > 0 && PrevPE2ChromID > 0))
			{
			pMateScaffold = LocateMateScaffold(PrevPE2ChromID,PrevPE1ChromID);
			if(pMateScaffold != NULL)
				{
				do {
					if(pMateScaffold->PE1Sense == pMateScaffold->PE2Sense)
						RevNumSenseSense += 1;
					else
						RevNumSenseAnti += 1;
					if(pMateScaffold->PEScafoldID == m_NumScaffolds)
						break;
					pMateScaffold += 1;
					}
				while(pMateScaffold->PE1ChromID == PrevPE2ChromID && pMateScaffold->PE2ChromID == PrevPE1ChromID);

				bRevMate = true;
				}
			}
		continue;
		}

	// same pair of chromosomes so accumulate counts
	if(pPEScaffold->PE1ChromID > 0 && pPEScaffold->PE2ChromID > 0)
		{
		NumPEAligned += 1;
		if(pPEScaffold->PE1Sense == pPEScaffold->PE2Sense)
			NumSenseSense += 1;
		else
			NumSenseAnti += 1;
		NumPaired += 1;
		}
	else
		{
		NumSEAligned += 1;
		NumUnpaired += 1;
		}
	}

if(BuffIdx > 0 || NumPEAligned > 0 || NumSEAligned > 0)
	{
	if(PrevPE1ChromID > 0 && PrevPE2ChromID > 0)
		BuffIdx += sprintf(&szBuff[BuffIdx],"\"%s\",\"%s\",%d,%d,%d,\"Y\",\"%s\",\"%s\",%d,%d\n",
							pszPE1Chrom,pszPE2Chrom,NumPEAligned,NumSenseSense,NumSenseAnti,PrevPE1ChromID==PrevPE2ChromID ? "Y" : "N",bRevMate ? "Y" : "N",RevNumSenseSense,RevNumSenseAnti);
	else
		{
		if(PrevPE1ChromID > 0)
			BuffIdx += sprintf(&szBuff[BuffIdx],"\"%s\",\"%s\",%d,%d,%d,\"N\",\"N\",\"N\",0,0\n",pszPE1Chrom,"N/A",NumSEAligned,0,0);
		else
			BuffIdx += sprintf(&szBuff[BuffIdx],"\"%s\",\"%s\",%d,%d,%d,\"N\",\"N\",\"N\",0,0\n","N/A",pszPE2Chrom,NumSEAligned,0,0);
		}

	CUtility::SafeWrite(m_hOutFile,szBuff,BuffIdx);
	}

close(m_hOutFile);
m_hOutFile = -1;

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Completed writing corelations (%d paired, %d orphaned) to file",NumPaired,NumUnpaired);
return(0);
}



int
Process(int PMode,					// processing mode
		char *pszInPE1File,			// input PE1 file
		char *pszInPE2File,			// input PE2 file
		char *pszOutFile)			// output corelations file
{
int Rslt;
size_t memreq;
Init();


#ifdef _WIN32
if((m_hOutFile = open(pszOutFile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
if((m_hOutFile = open(pszOutFile,O_RDWR | O_CREAT |O_TRUNC, S_IREAD | S_IWRITE))==-1)
#endif
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create or truncate %s - %s",pszOutFile,strerror(errno));
	Reset();
	return(eBSFerrCreateFile);
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Output results file created/truncated: '%s'",pszOutFile);

memreq = (size_t)(sizeof(tsPEScaffoldChrom) * cAllocChromNames);	
#ifdef _WIN32
m_pScaffoldChroms = (tsPEScaffoldChrom *) malloc(memreq);	// initial and perhaps the only allocation

if(m_pScaffoldChroms == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadReads: Memory allocation of %lld bytes for m_pScaffoldChroms - %s",(INT64)memreq,strerror(errno));
	Reset();
	return(eBSFerrMem);
	}
#else
	// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
m_pScaffoldChroms = (tsPEScaffoldChrom *)mmap(NULL,memreq, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
if(m_pScaffoldChroms == MAP_FAILED)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadReads: Memory allocation of %lld bytes for m_pScaffoldChroms through mmap()  failed - %s",(INT64)memreq,strerror(errno));
	m_pScaffoldChroms = NULL;
	Reset();
	return(eBSFerrMem);
	}
#endif
m_AllocdScaffoldChromsMem = memreq;
m_AllocdNumScaffoldChroms = cAllocChromNames;
m_NumScaffoldChroms = 0;

m_pHashChroms = new int [cHashSize];
if(m_pHashChroms == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadReads: Memory allocation of %d bytes for m_pHashChroms - %s",cHashSize * sizeof(int),strerror(errno));
	Reset();
	return(eBSFerrMem);
	}
memset(m_pHashChroms,0,sizeof(int) * cHashSize);

memreq = (size_t)(sizeof(tsPEScaffold) * cAllocScafolds);	
#ifdef _WIN32
m_pScaffolds = (tsPEScaffold *) malloc(memreq);	// initial and perhaps the only allocation

if(m_pScaffolds == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadReads: Memory allocation of %lld bytes for m_pScaffolds - %s",(INT64)memreq,strerror(errno));
	Reset();
	return(eBSFerrMem);
	}
#else
	// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
m_pScaffolds = (tsPEScaffold *)mmap(NULL,memreq, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
if(m_pScaffolds == MAP_FAILED)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadReads: Memory allocation of %lld bytes for m_pScaffolds through mmap()  failed - %s",(INT64)memreq,strerror(errno));
	m_pScaffolds = NULL;
	Reset();
	return(eBSFerrMem);
	}
#endif
m_AllocdScaffoldsMem = memreq;
m_AllocdNumScaffolds = cAllocScafolds;
m_NumScaffolds = 0;

memreq = (size_t)(sizeof(tsPEIdent) * cAllocPENames);	
#ifdef _WIN32
m_pPEIdents = (tsPEIdent *) malloc(memreq);	// initial and perhaps the only allocation

if(m_pPEIdents == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadReads: Memory allocation of %lld bytes for m_pPEIdents - %s",(INT64)memreq,strerror(errno));
	Reset();
	return(eBSFerrMem);
	}
#else
	// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
m_pPEIdents = (tsPEIdent *)mmap(NULL,memreq, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
if(m_pPEIdents == MAP_FAILED)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadReads: Memory allocation of %lld bytes for m_pPEIdents through mmap()  failed - %s",(INT64)memreq,strerror(errno));
	m_pPEIdents = NULL;
	Reset();
	return(eBSFerrMem);
	}
#endif
m_AllocdPEIdentsMem = memreq;
m_AllocdNumPEIdents = cAllocPENames;
m_NumPEIdents = 0;

m_pHashPEIdents = new int [cHashSize];
if(m_pHashPEIdents == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadReads: Memory allocation of %d bytes for m_pHashPEIdents - %s",cHashSize * sizeof(int),strerror(errno));
	Reset();
	return(eBSFerrMem);
	}
memset(m_pHashPEIdents,0,sizeof(int) * cHashSize);

Rslt = LoadSAM(false,pszInPE1File);
if(Rslt < 1)
	return(Rslt);

Rslt = LoadSAM(true,pszInPE2File);
if(Rslt < 1)
	return(Rslt);

if(m_NumScaffolds > 1)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Sorting %d scaffolds by PE1ChromID.PE2ChromID ascending",m_NumScaffolds);
	mtqsort.qsort(m_pScaffolds,m_NumScaffolds,sizeof(tsPEScaffold),SortScaffolds);
	tsPEScaffold *pScaffold;
	int Idx;
	pScaffold = m_pScaffolds;
	for(Idx=1;Idx <= m_NumScaffolds; Idx++,pScaffold++)
		pScaffold->PEScafoldID = Idx;
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Sorting scaffolds completed");
	}

ReportCorelationships(pszOutFile);

Reset();
return(0);
}


// Sort scaffolds by PE1ChromID.PE2ChromID ascending
static int
SortScaffolds(const void *arg1, const void *arg2)
{
tsPEScaffold *pEl1 = (tsPEScaffold *)arg1;
tsPEScaffold *pEl2 = (tsPEScaffold *)arg2;

if(pEl1->PE1ChromID > pEl2->PE1ChromID)
	return(1);
if(pEl1->PE1ChromID < pEl2->PE1ChromID)
	return(-1);
if(pEl1->PE2ChromID > pEl2->PE2ChromID)
	return(1);
if(pEl1->PE2ChromID < pEl2->PE2ChromID)
	return(-1);
return(0);
}