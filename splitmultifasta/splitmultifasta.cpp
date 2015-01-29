// splitmultifasta.cpp : Defines the entry point for the console application.
// Reads in specified fasta file containing multiple entries and writes out one fasta file per
// entry
// Each generated multifasta file will be number sequentially according to a user specified prefix
// The generated multifasta files can be limited by size or by number of entries


#include "stdafx.h"
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#if _WIN32

#include "../libbiokanga/commhdrs.h"
#else
#include "../libbiokanga/commhdrs.h"
#endif

const unsigned int cProgVer = 202;		// increment with each release


const int cMaxInBuffSize  = 10000000;	// read in chunks of this size from source fasta file
const int cMaxOutBuffSize = 10000000;	// and write in chunks of this size to output fasta file

const int cDfltMaxTotSeqLen = 100000000; // default max total sequence length
const int cDfltMaxFastaEntries = 100;	 // default max fasta entries in generated output fasta file



typedef struct TAG_sProcParams 
	{
	int MaxTotSeqLen;	// don't start new entries in current output file if total sequence length >= MaxTotSeqLen	
	int MaxFastaEntries; // limit number of entries in current output file to no more than MaxFastaEntries
	int hInFile;		// opened file handle for source multifasta file
	int hOutFile;		// opened file handle for current multifasta file being written
	int OutFileIdx;	// current output file sequential prefix
	int PushedBack;	// last pushed back char (only 1 level of pushback supported!)
	unsigned char *pInBuffer; // mem allocd to buffer chars being read from fasta
	int NumInBuffer;		   // num of chars currently in pInBuffer
	int InBuffIdx;			   // index of next char to read from pInBuffer[]
	unsigned char *pOutBuffer; // mem allocd to buffer chars being written to fasta
	int OutBuffIdx;			// number of chars in pOutBuffer[] ready to be written
	char szOutPath[_MAX_PATH]; // path on which to create output files (contains format char '%d' ready for seq number)
	char szOutFasta[_MAX_PATH];	// path + output file
	char szInFasta[_MAX_PATH];   // path + input file

	int TotSeqLen;			// current total seq length in output file
	int NumFastaEntries;	// current number of entries in output file
	char szGenome[128];		// descriptor line genome prefix
	} tsProcParams;


int ProcessFastaEntry(tsProcParams *pParams);
int WriteFasta(int Chr,tsProcParams *pParams);
int PushBac(int Chr,tsProcParams *pParams);
int GetNext(tsProcParams *pParams);
int Process(char *pszInFasta,char *pszOutPath,char *pszOutPrefix,char *pszGenome,int MaxTotSeqLen,int MaxFastaEntries);

CStopWatch gStopWatch;
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
int iScreenLogLevel;		// level of screen diagnostics
int iFileLogLevel;			// level of file diagnostics
char szLogFile[_MAX_PATH];	// write diagnostics to this file
int Rslt;

int iMaxTotSeqLen;
int iMaxFastaEntries;
char szInFile[_MAX_PATH];	// process multifasta from this file
char szOutPath[_MAX_PATH];	// output path and suffix of files containing entries
char szPrefix[128];			// output files generated will have this prefix
char szGenome[128];			// descriptor line prefix

// command line args
struct arg_lit  *help    = arg_lit0("hH","help",                "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_int *ScreenLogLevel=arg_int0("S", "ScreenLogLevel",	"<int>","Level of diagnostics written to logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"diagnostics log file");

struct arg_file *InFile = arg_file1("i","infile","<file>",	"process from these multifasta files");
struct arg_file *OutPath = arg_file1("o","outfile","<file>","write entries to this multifasta file (name prefixed by -p<prefix>id)");
struct arg_int  *MaxTotSeqLen = arg_int0("l","maxtotlen","<int>", "write max total seq length (default 100000000)");
struct arg_int  *MaxFastaEntries = arg_int0("e","maxentries","<int>","write max entries (default 100)");
struct arg_str  *Prefix = arg_str0("p","prefix","<string>","generated output file name prefix (default 'fsplit')");
struct arg_str  *Genome = arg_str0("g","genome","<string>","genome prefix to insert into descriptor lines - e.g. hg18");
struct arg_end *end = arg_end(20);

void *argtable[] = {help,version,FileLogLevel,ScreenLogLevel,LogFile,
					InFile,OutPath,MaxTotSeqLen,MaxFastaEntries,
					Prefix,Genome,
					end};

char **pAllArgs;
int argerrors;
argerrors = CUtility::arg_parsefromfile(argc,(char **)argv,&pAllArgs);
if(argerrors >= 0)
	argerrors = arg_parse(argerrors,pAllArgs,argtable);

    /* special case: '--help' takes precedence over error reporting */
if (help->count > 0)
        {
		printf("Usage: %s ", gszProcName);
        arg_print_syntax(stdout,argtable,"\n");
        arg_print_glossary(stdout,argtable,"  %-25s %s\n");
		exit(1);
        }

    /* special case: '--version' takes precedence error reporting */
if (version->count > 0)
        {
			printf("\n%s Version: %d.%2.2d\n",gszProcName,cProgVer/100,cProgVer%100);
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

	iMaxTotSeqLen = MaxTotSeqLen->count ? MaxTotSeqLen->ival[0] : cDfltMaxTotSeqLen;
	if(iMaxTotSeqLen < 1)
		{
		printf("\nError: Requested MaxTotSeqLen '-l%d' not supported",iMaxTotSeqLen);
		exit(1);
		}
	iMaxFastaEntries = MaxFastaEntries->count ? MaxFastaEntries->ival[0] : cDfltMaxTotSeqLen;
	if(iMaxTotSeqLen < 1)
		{
		printf("\nError: Requested MaxFastaEntries '-e%d' not supported",iMaxFastaEntries);
		exit(1);
		}

	if(!Prefix->count)
		strcpy(szPrefix,"fsplit");
	else
		strcpy(szPrefix,Prefix->sval[0]);

	if(!Genome->count)
		szGenome[0] = '\0';
	else
		strcpy(szGenome,Genome->sval[0]);

	strcpy(szInFile,InFile->filename[0]);
	strcpy(szOutPath,OutPath->filename[0]);

		// now that command parameters have been parsed then initialise diagnostics log system
	if(!gDiagnostics.Open(szLogFile,(etDiagLevel)iScreenLogLevel,(etDiagLevel)iFileLogLevel,true))
		{
		printf("\nError: Unable to start diagnostics subsystem.");
		if(szLogFile[0] != '\0')
			printf(" Most likely cause is that logfile '%s' can't be opened/created",szLogFile);
		exit(1);
		}

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Version: %d.%2.2d Processing parameters:",cProgVer/100,cProgVer%100);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"input multifasta files: '%s'",szInFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"output multifasta files:   '%s'",szOutPath);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"output file prefix: '%s'",szPrefix);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"descriptor line prefix: '%s'",szGenome[0] == '\0' ? "NONE SUPPLIED" : szGenome);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"max number of entries in output files: %d",iMaxFastaEntries);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"max total sequence lengths in output files: %d",iMaxTotSeqLen);

	gStopWatch.Start();
#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	Rslt = Process(szInFile,szOutPath,szPrefix,szGenome,iMaxTotSeqLen,iMaxFastaEntries);
	gStopWatch.Stop();
	Rslt = Rslt < 0 ? 1 : 0;
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Exit Code: %d Total processing time: %s",Rslt,gStopWatch.Read());
	exit(Rslt);
	}
else
	{
	arg_print_errors(stdout,end,gszProcName);
	arg_print_syntax(stdout,argtable,"\nend of help\n");
	exit(1);
	}
}

int
ProcessThisFile(char *pszInFasta,void *pHandlerParams)
{
int Rslt;
tsProcParams *pParams = (tsProcParams *)pHandlerParams;
strcpy(pParams->szInFasta,pszInFasta);
#ifdef _WIN32
if((pParams->hInFile = _open(pszInFasta,_O_RDWR | _O_BINARY | _O_SEQUENTIAL))==-1)
#else
if((pParams->hInFile = open(pszInFasta,O_RDWR))==-1)
#endif
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: Unable to open input file for processing - '%s' - %s", 
				pszInFasta,strerror(errno));
	return(eBSFerrOpnFile);
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"ProcessThisFile: Processing input file '%s'",pszInFasta); 

pParams->NumInBuffer = 0;
pParams->InBuffIdx = 0;

Rslt = 1;
while(Rslt == 1) {
	if(pParams->hOutFile == -1)
		{
		sprintf(pParams->szOutFasta,pParams->szOutPath,pParams->OutFileIdx++);
#ifdef _WIN32
		if((pParams->hOutFile = open(pParams->szOutFasta,( _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC),(_S_IREAD | _S_IWRITE)))==-1)
#else
		if((pParams->hOutFile = open(pParams->szOutFasta,O_RDWR | O_CREAT |O_TRUNC,S_IREAD | S_IWRITE))==-1)

#endif
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: Unable to open output file for processing - '%s' - %s", 
					pParams->szOutFasta,strerror(errno));
			close(pParams->hInFile);
			pParams->hInFile = -1;
			return(eBSFerrOpnFile);
			}
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"ProcessThisFile: Generating output file '%s'",pParams->szOutFasta); 
		pParams->OutBuffIdx = 0;
		pParams->TotSeqLen = 0;
		pParams->NumFastaEntries = 0;
		}
	Rslt = ProcessFastaEntry(pParams);
	switch(Rslt) {
		case 0:			// input file fully processed
			close(pParams->hInFile);
			pParams->hInFile = -1;
			break;

		case 1:			// new output file required
			if(pParams->OutBuffIdx > 0)
				CUtility::SafeWrite(pParams->hOutFile,pParams->pOutBuffer,pParams->OutBuffIdx);
			close(pParams->hOutFile);
			pParams->hOutFile = -1;
			break;

		default:		// error
			if(pParams->hInFile != -1)
				{
				close(pParams->hInFile);
				pParams->hInFile = -1;
				}
			if(pParams->hOutFile != -1)
				{
				close(pParams->hOutFile);
				pParams->hOutFile = -1;
				}
			break;
		}
	}
return(Rslt);
}


int
Process(char *pszInFasta,char *pszOutPath,char *pszOutPrefix,char *pszGenome,int MaxTotSeqLen,int MaxFastaEntries)
{
int Rslt;
char szDirPath[_MAX_PATH];
char szFileSpec[_MAX_PATH];
tsProcParams ProcParams;		// initialised to hold processing parameters
memset(&ProcParams,0,sizeof(ProcParams));

#ifdef _WIN32
char szDrive[_MAX_DRIVE];
char szDir[_MAX_DIR];
char szFname[_MAX_FNAME];
char szExt[_MAX_EXT];
_splitpath(pszOutPath,szDrive,szDir,szFname,szExt);
sprintf(szFileSpec,"%s%%d%s",pszOutPrefix,szFname);
_makepath(ProcParams.szOutPath,szDrive,szDir,szFileSpec,szExt);

_splitpath(pszInFasta,szDrive,szDir,szFname,szExt);
_makepath(szDirPath,szDrive,szDir,"","");
if(!szFname[0])
	strcpy_s(szFname,sizeof(szFname),".*");
if(!szExt[0])
	strcpy_s(szExt,sizeof(szExt),"*");
sprintf(szFileSpec,"%s%s",szFname,szExt);
#else

CUtility::splitpath(pszOutPath,ProcParams.szOutPath,szFileSpec);
sprintf(szDirPath,"%s%%d%s",pszOutPrefix,szFileSpec);
strcat(ProcParams.szOutPath,szDirPath);
CUtility::splitpath(pszInFasta,szDirPath,szFileSpec);
if(!szFileSpec[0])
	strcpy(szFileSpec,"*");
#endif

ProcParams.OutFileIdx = 1;				// generated files start from entry 1
ProcParams.hInFile = -1;
ProcParams.hOutFile = -1;
ProcParams.PushedBack = 0;
ProcParams.InBuffIdx = 0;
ProcParams.OutBuffIdx = 0;
ProcParams.NumInBuffer = 0;
ProcParams.MaxFastaEntries = MaxFastaEntries;
ProcParams.MaxTotSeqLen = MaxTotSeqLen;
if(pszGenome != NULL && pszGenome[0] != '\0')
	strcpy(ProcParams.szGenome,pszGenome);
else
	ProcParams.szGenome[0] = '\0';

if((ProcParams.pInBuffer = new unsigned char [cMaxInBuffSize])==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: Unable to allocate memory (%d bytes) for input buffering", 
				cMaxInBuffSize);
	return(eBSFerrMem);
	}

if((ProcParams.pOutBuffer = new unsigned char [cMaxOutBuffSize])==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: Unable to allocate memory (%d bytes) for output buffering", 
				cMaxOutBuffSize);
	delete ProcParams.pInBuffer;
	return(eBSFerrMem);
	}

Rslt = eBSFSuccess;
CSimpleGlob glob(SG_GLOB_FULLSORT);
if (glob.Add(szDirPath) >= SG_SUCCESS)
	{
	for (int n = 0; Rslt >= eBSFSuccess &&  n < glob.FileCount(); ++n)
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Will process in this order: %d '%s", n+1,glob.File(n));

    for (int n = 0; Rslt >= eBSFSuccess &&  n < glob.FileCount(); ++n)
		Rslt = ProcessThisFile(glob.File(n),&ProcParams);
	}
else
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to glob '%s",szDirPath);
	Rslt = eBSFerrOpnFile;	// treat as though unable to open file
    }

if(ProcParams.hOutFile != -1)			// close output file if still opened
	{
	if(ProcParams.OutBuffIdx && Rslt >= 0)
		CUtility::SafeWrite(ProcParams.hOutFile,ProcParams.pOutBuffer,ProcParams.OutBuffIdx);
	close(ProcParams.hOutFile);
	}
if(ProcParams.hInFile != -1)
	close(ProcParams.hInFile);
if(ProcParams.pInBuffer != NULL)
	delete ProcParams.pInBuffer;
if(ProcParams.pOutBuffer != NULL)
	delete ProcParams.pOutBuffer;
return(Rslt);
}


int		// 0: EOF -1: error >0 chr
GetNext(tsProcParams *pParams)
{
int Chr;
if((Chr = pParams->PushedBack) > 0)
	{
	pParams->PushedBack = 0;
	return(Chr);
	}
if(pParams->InBuffIdx == -1 || pParams->InBuffIdx >= pParams->NumInBuffer)
	{
	pParams->NumInBuffer = read(pParams->hInFile,pParams->pInBuffer,cMaxInBuffSize);
	if(pParams->NumInBuffer <= 0)
		{
		pParams->InBuffIdx = -1;
		return(pParams->NumInBuffer);
		}
	pParams->InBuffIdx = 0;
	}
return(pParams->pInBuffer[pParams->InBuffIdx++]);
}

int
PushBac(int Chr,tsProcParams *pParams)
{
pParams->PushedBack = Chr;
return(0);
}

int
WriteFasta(int Chr,tsProcParams *pParams)
{
if(pParams->OutBuffIdx < cMaxOutBuffSize)
	{
	pParams->pOutBuffer[pParams->OutBuffIdx++] = (char)Chr;
	return(0);
	}
if(pParams->hOutFile != -1)			
	{
	CUtility::SafeWrite(pParams->hOutFile,pParams->pOutBuffer,pParams->OutBuffIdx);
	pParams->OutBuffIdx = 0;
	pParams->pOutBuffer[pParams->OutBuffIdx++] = (char)Chr;
	return(0);
	}
return(-1);
}

// processing loop
// copy from current file input to output file until:
// a) number entries written == MaxFastaEntries
// OR
// b) new entry starting and total sequence written to file is >= MaxTotSeqLen
// NOTE: if MaxTotSeqLen <= 0 then MaxTotSeqLen is set to cDfltMaxTotSeqLen
// NOTE: if MaxFastaEntries <= 0 then MaxFastaEntries is set to cDfltMaxFastaEntries
//
int						// <0 error, 0 = end of input file, 1 new output file required
ProcessFastaEntry(tsProcParams *pParams)
{
int Rslt;
bool bInDescr = false;
char *pGenome;

while((Rslt = GetNext(pParams)) > 0) {
	switch((char)Rslt) {
		case '\r':			// silently slough CR - must have been generated on windows/msdos machine
			continue;

		case '\n':			// accept linefeeds - both Linux and windows are happy with fasta lines terminated by NL
			if(bInDescr)
				bInDescr = false;
			break;
		
		case '>':
			if(!bInDescr)
				{
				bInDescr = true;
				if(pParams->NumFastaEntries == pParams->MaxFastaEntries || pParams->TotSeqLen >= pParams->MaxTotSeqLen)
					{
					PushBac(Rslt,pParams);
					return(1);
					}
				pParams->NumFastaEntries++;
				// following processing is to handle inserting szGenome prefix between
				// the '>' and existing descriptor
				// idea is that if szGenome was specified as say 'hg18' and descriptor is '>chr10' then
				// the descriptor line is updated to be '>hg18.chr10'
				if((Rslt = WriteFasta((char)Rslt,pParams)) < 0) // error writing fasta?
					return(Rslt);
				if(pParams->szGenome[0] == '\0')
					continue;
				pGenome = pParams->szGenome;
				while((Rslt = *pGenome++)!='\0')
					if((Rslt = WriteFasta((char)Rslt,pParams)) < 0) // error writing fasta?
						return(Rslt);
				if((Rslt = WriteFasta('.',pParams)) < 0) // error writing fasta?
						return(Rslt);
				continue;
				}
			break;
			
		default:
			if(bInDescr)			// can accept whitespace within descriptors
				break;				// and descriptor text does not count against MaxTotSeqLen 
			if(isspace(Rslt))		// slough whitespace within sequence
				continue;
			pParams->TotSeqLen++;
			break;
		}
	if((Rslt = WriteFasta(Rslt,pParams)) < 0) // error writing fasta?
		return(Rslt);
	}
return(Rslt);
}



