// fastafilter.cpp : Defines the entry point for the console application.

/*
 * CSIRO Open Source Software License Agreement (GPLv3)
 * Copyright (c) 2017, Commonwealth Scientific and Industrial Research Organisation (CSIRO) ABN 41 687 119 230.
 * See LICENSE for the complete license information (https://github.com/csiro-crop-informatics/biokanga/LICENSE)
 * Contact: Alex Whan <alex.whan@csiro.au>
 */

// Filter fasta files
// Initially just has functionality for detecting duplicate sequence identifiers and forcing these to be unique by
// appending unique incrementing decimal values; additionally sequences containing more than specified numbers of consecutive indeterminates will have these runs truncated
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

const char *cpszProgVer = "0.0.2";		// increment with each release


const char *pszDfltIdent = "Seq";			// if no identifier then default to this sequence identifier
const size_t cAllocFastaIdents = 10000000;	// allocate for fasta identifiers in this sized increments
const int cMaxSepLen = 5;					// allow at most 5 chars to be used as unique identifier separators
const int cMaxIdentNameLen = 50;			// truncate identifiers to be at most this length excluding the initial '>'
const int cMaxDescrLineLen = 120;			// truncate fasta descriptor lines to be at most this length including the initial '>'


int Process(int PMode,
			int MaxNrun,				// truncate runs of 'N' to be no more than this length
			char *pszSepUnique,			// use these as the unique separator text
			char *pszInFasta,			// load from this fasta file 
			char *pszOutFile);			// output to this fasta file

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
int MaxNrun;				// truncate runs of 'N' to be no more than this length
char szSepUnique[cMaxSepLen+1];	// use these as the unique separator text
char szInFile[_MAX_PATH];	// load fasta from this file
char szOutFile[_MAX_PATH];	// write fasta to this file

// command line args
struct arg_lit  *help    = arg_lit0("h","help",                 "Print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"Print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to screen and logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"Diagnostics log file");

struct arg_int *mode = arg_int0("m","mode","<int>",				"Processing mode: 0 filter, 1 reverse complement all sequences (default 0)");
struct arg_int *maxnrun = arg_int0("n","maxnrun","<int>",		"Limit runs of indeterminates 'N's to be no more than this (defaults to 10)");
struct arg_str *sepunique = arg_str0("s","sepunique","<str>",   "Separator to use when suffixing duplicate identifiers (defaults to '.')");
struct arg_file *infile = arg_file1("i","in","<file>",			"Input fasta file");
struct arg_file *outfile = arg_file1("o","out","<file>",		"Output fasta file");

struct arg_end *end = arg_end(20);

void *argtable[] = {help,version,FileLogLevel,LogFile,
					mode,maxnrun,sepunique,infile,outfile,
					end};

char **pAllArgs;
int argerrors;
argerrors = CUtility::arg_parsefromfile(argc,(char **)argv,&pAllArgs);
if(argerrors >= 0)
	argerrors = arg_parse(argerrors,pAllArgs,argtable);

/* special case: '--help' takes precedence over error reporting */
if (help->count > 0)
        {
		printf("\n%s Duplicate fasta identifier rename, Version %s\nOptions ---\n", gszProcName,cpszProgVer);
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

	PMode = mode->count ? mode->ival[0] : 0;
	if(PMode < 0 || PMode > 1)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Processing mode '-m%d' must be in range 0..0",PMode);
		return(1);
		}

	MaxNrun = maxnrun->count ? maxnrun->ival[0] : 10;
	if(MaxNrun < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Max indeterminate run of Ns '-n%d' must be >= 0",MaxNrun);
		return(1);
		}

	szSepUnique[0] = '.';
	szSepUnique[1] = '\0';
	if(sepunique->count)
		{
		strncpy(szSepUnique,sepunique->sval[0],cMaxSepLen);
		szSepUnique[cMaxSepLen] = '\0';
		CUtility::TrimQuotedWhitespcExtd(szSepUnique);
		if(szSepUnique[0] == '\0')
			{
			szSepUnique[0] = '.';
			szSepUnique[1] = '\0';
			}
		}
	
	strcpy(szInFile,infile->filename[0]);
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
			pszDescr = "Rename duplicate sequence identifiers";
			break;
		case 1:
			pszDescr = "Reverse complement all sequences";
		}

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Processing mode is : '%s'",pszDescr);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Limit indeterminate run lengths to: %d",MaxNrun);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Unique separator : '%s'",szSepUnique);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Input file: '%s'",szInFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Output file: '%s'",szOutFile);

#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	gStopWatch.Start();

	Rslt = Process(PMode,MaxNrun,szSepUnique,szInFile,szOutFile);
	Rslt = Rslt >=0 ? 0 : 1;
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


typedef struct TAG_sFastaIdentifier {
	int IdentID;						// uniquely identifies this identifier
	int HashNext;					// offset in m_pFastaIdentifiers[] of next fasta identifier with same identifer hash
	int NumInstances;					// number of instances of this identifier 
	char szIdentifier[cMaxIdentNameLen+1];
	} tsFastaIdentifier;

int m_NumFastaIdents;					// currently this many fasta identifiers have been loaded into m_pFastaIdentifiers
int m_NumFastaIdentAlloc;				// allocation is for this many fasta identifiers
size_t m_FastaIdentAllocMem;			// how much memory has been allocated to hold m_pFastaIdentifiers
tsFastaIdentifier *m_pFastaIdentifiers;

int m_IdentHashes[0x0ffff];				// hashed indexes into m_pFastaIdentifiers

char m_szLineBuff[0x0ffffff];			// used to buffer fasta lines
char m_szSeqBuff[0x0fffffff];			// hold sequences with truncated runs of Ns

// AddFastaIdentifier
// If Identifier already known then returns number of instances this identifier has been returned
// otherwise add as being the first instance to m_pFastaIdentifiers
int
AddFastaIdentifier(char *pszFastaIdentifier)
{
static int PrevIdentID = 0;
UINT16 Hash;
int HashIdx;
tsFastaIdentifier *pFastaIdentifier;

// check to see if identifier already known
if(PrevIdentID != 0)
	{
	pFastaIdentifier = &m_pFastaIdentifiers[PrevIdentID-1];
	if(!stricmp(pszFastaIdentifier,pFastaIdentifier->szIdentifier))
		{
		pFastaIdentifier->NumInstances += 1;
		return(pFastaIdentifier->NumInstances);
		}
	}

Hash = CUtility::GenHash16(pszFastaIdentifier);
if((HashIdx = m_IdentHashes[Hash]) != 0)
	{
	do {
		pFastaIdentifier = &m_pFastaIdentifiers[HashIdx-1];
		if(!stricmp(pszFastaIdentifier,pFastaIdentifier->szIdentifier))
			{
			PrevIdentID = pFastaIdentifier->IdentID;
			pFastaIdentifier->NumInstances += 1;
			return(pFastaIdentifier->NumInstances);
			}
		HashIdx = pFastaIdentifier->HashNext;
		}
	while(HashIdx > 0);
	}

// its a new chrom not previously seen
// realloc as may be required to hold this new chrom
if(m_NumFastaIdents == m_NumFastaIdentAlloc)
	{
	size_t memreq = m_FastaIdentAllocMem + (cAllocFastaIdents * sizeof(tsFastaIdentifier));
#ifdef _WIN32
	pFastaIdentifier = (tsFastaIdentifier *) realloc(m_pFastaIdentifiers,memreq);
#else
	pFastaIdentifier = (tsFastaIdentifier *)mremap(m_pFastaIdentifiers,m_FastaIdentAllocMem,memreq,MREMAP_MAYMOVE);
	if(pFastaIdentifier == MAP_FAILED)
		pFastaIdentifier = NULL;
#endif
	if(pFastaIdentifier == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddFastaIdentifier: Memory re-allocation to %d bytes - %s",memreq,strerror(errno));
		return(eBSFerrMem);
		}
	m_pFastaIdentifiers = pFastaIdentifier;
	m_FastaIdentAllocMem = memreq;
	m_NumFastaIdentAlloc += cAllocFastaIdents;
	}

pFastaIdentifier = &m_pFastaIdentifiers[m_NumFastaIdents++];
pFastaIdentifier->IdentID = m_NumFastaIdents;
pFastaIdentifier->NumInstances = 1;
pFastaIdentifier->HashNext = m_IdentHashes[Hash];
m_IdentHashes[Hash] = m_NumFastaIdents;
strncpy(pFastaIdentifier->szIdentifier,pszFastaIdentifier,cMaxIdentNameLen);
pFastaIdentifier->szIdentifier[cMaxIdentNameLen] = '\0';
PrevIdentID = m_NumFastaIdents;
return(1);
}


void
Reset(void)
{
if(m_pFastaIdentifiers != NULL)
	{
#ifdef _WIN32
	free(m_pFastaIdentifiers);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pFastaIdentifiers != MAP_FAILED)
		munmap(m_pFastaIdentifiers,m_FastaIdentAllocMem);
#endif
	m_pFastaIdentifiers = NULL;
	}
m_NumFastaIdents = 0;
m_NumFastaIdentAlloc = 0;
m_FastaIdentAllocMem = 0;
memset(m_IdentHashes,0,sizeof(m_IdentHashes));
}

void Init(void)
{
m_pFastaIdentifiers = NULL;
Reset();
}


const UINT32 cMaxFastaLen = 500000000;  // allowing for fasta sequences of up to this max length
int
RevCplSeqs(char *pszInFasta,			// reverse complement fasta sequences in this file
		   char *pszOutFasta)			// and write out to this file
{
int Rslt;
int hOutFile;
bool bFirstEntry;
unsigned char *pSeqBuff;
unsigned char *pszSeqBuff;
char szDescription[cBSFDescriptionSize];
int SeqLen;
int Descrlen;
int ChrIdx;
CFasta InFasta;

hOutFile = -1;
if((pSeqBuff = new unsigned char [cMaxFastaLen]) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessFastaFile:- Unable to allocate memory (%d bytes) for raw sequence buffer",cMaxFastaLen);
	Reset();
	return(eBSFerrMem);
	}
if((pszSeqBuff = new unsigned char [cMaxFastaLen]) == NULL)
	{
	delete pSeqBuff;
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessFastaFile:- Unable to allocate memory (%d bytes) for char sequence buffer",cMaxFastaLen);
	Reset();
	return(eBSFerrMem);
	}



if((Rslt=InFasta.Open(pszInFasta,true))!=eBSFSuccess)
	{
	delete pSeqBuff;
	delete pszSeqBuff;
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open '%s' [%s] %s",pszInFasta,InFasta.ErrText((teBSFrsltCodes)Rslt),InFasta.GetErrMsg());
	Reset();
	return(Rslt);
	}

#ifdef _WIN32
hOutFile = open(pszOutFasta,O_CREATETRUNC );
#else
if((hOutFile = open(pszOutFasta,O_RDWR | O_CREAT,S_IREAD | S_IWRITE))!=-1)
    if(ftruncate(hOutFile,0)!=0)
			{
			delete pSeqBuff;
			delete pszSeqBuff;
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate %s - %s",pszOutFasta,strerror(errno));
			Reset();
			return(eBSFerrCreateFile);
			}
#endif

if(hOutFile < 0)
	{
	delete pSeqBuff;
	delete pszSeqBuff;
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate output file '%s'",pszOutFasta);
	Reset();
	return(eBSFerrCreateFile);
	}

szDescription[0] = '\0';
ChrIdx = 0;
bFirstEntry = true;
while((Rslt = SeqLen = InFasta.ReadSequence(pSeqBuff,cMaxFastaLen-1,true,false)) > eBSFSuccess)
	{
	if(SeqLen == eBSFFastaDescr)		// just read a descriptor line
		{
		ChrIdx = 0;
		szDescription[0] = '>';
		Descrlen = 1 + InFasta.ReadDescriptor(&szDescription[1],cBSFDescriptionSize-2);
		szDescription[Descrlen++] = '\n';
		CUtility::SafeWrite(hOutFile,szDescription,Descrlen);
		ChrIdx = 0;
		bFirstEntry = false;
		continue;
		}

	CSeqTrans::ReverseComplement(SeqLen,pSeqBuff);

	int BaseIdx = 0;
	int BasesLine;
	ChrIdx = 0;
	do {
		BasesLine = min(80,SeqLen - BaseIdx);
		CSeqTrans::MapSeq2Ascii(&pSeqBuff[BaseIdx],BasesLine,(char *)&pszSeqBuff[ChrIdx]);
		BaseIdx += BasesLine;
		ChrIdx += BasesLine;
		pszSeqBuff[ChrIdx] = '\n';
		ChrIdx += 1;
		}
	while(BaseIdx < SeqLen);
	if(pszSeqBuff[ChrIdx-1]!='\n')
		pszSeqBuff[ChrIdx++] ='\n';	
	CUtility::SafeWrite(hOutFile,pszSeqBuff,ChrIdx);
	ChrIdx = 0;
	}

if(hOutFile != -1)
	{
#ifdef _WIN32
	_commit(hOutFile);
#else
	fsync(hOutFile);
#endif
	close(hOutFile);
	hOutFile = -1;
	}

return(eBSFSuccess);
}

int Process(int PMode,
			int MaxNrun,				// truncate runs of 'N' to be no more than this length
			char *pszSepUnique,			// use these as the unique separator text
			char *pszInFasta,			// load from this fasta file 
			char *pszOutFasta)			// output to this fasta file
{
UINT32 LineNum;
int CurNrun;
int CurSeqLen;
int IdentLen;
char Chr;
char *pSrc;
char *pDst;
size_t memreq;
CFasta Fasta;
int NumInstances;
char szFastaDescr[cMaxIdentNameLen+10+cMaxDescrLineLen+1];
FILE *pInStream;
FILE *pOutStream;

if(PMode == 1)               // reverse complement sequences; no changes to sequence descriptors
	return(RevCplSeqs(pszInFasta,pszOutFasta));

memreq = (size_t)(sizeof(tsFastaIdentifier) * cAllocFastaIdents);	
#ifdef _WIN32
m_pFastaIdentifiers = (tsFastaIdentifier *) malloc(memreq);	// initial and perhaps the only allocation

if(m_pFastaIdentifiers == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: Memory allocation of %lld bytes for m_pFastaIdentifiers - %s",(INT64)memreq,strerror(errno));
	Reset();
	return(eBSFerrMem);
	}
#else
	// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
m_pFastaIdentifiers = (tsFastaIdentifier *)mmap(NULL,memreq, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
if(m_pFastaIdentifiers == MAP_FAILED)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: Memory allocation of %lld bytes for m_pFastaIdentifiers through mmap()  failed - %s",(INT64)memreq,strerror(errno));
	m_pFastaIdentifiers = NULL;
	Reset();
	return(eBSFerrMem);
	}
#endif
m_FastaIdentAllocMem = memreq;
m_NumFastaIdentAlloc = cAllocFastaIdents;
m_NumFastaIdents = 0;

if((pInStream = fopen(pszInFasta,"r"))==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to file %s for reading error: %s",pszInFasta,strerror(errno));
	Reset();
	return(eBSFerrOpnFile);
	}

if((pOutStream = fopen(pszOutFasta,"w"))==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create/truncate file %s for writing error: %s",pszInFasta,strerror(errno));
	Reset();
	return(eBSFerrOpnFile);
	}

int NumDups;
int CurLineLen;
int DescrLineLen;
int NumEmptyDescrs;
int NumOverLenIdents;
int NumOverLenDescr;

NumOverLenIdents = 0;
NumEmptyDescrs = 0;
NumOverLenDescr = 0;
NumDups = 0;
CurNrun = 0;
CurSeqLen = 0;
LineNum = 0;
while(fgets(m_szLineBuff,sizeof(m_szLineBuff)-1,pInStream)!= NULL)
	{
	LineNum += 1;
	CurLineLen = (int)strlen(m_szLineBuff);
	if(CurLineLen > 0)
		{
		if(m_szLineBuff[CurLineLen-1] == '\n')
			{
			m_szLineBuff[CurLineLen-1] = '\0';
			CurLineLen -= 1;
			}
		}

	if(m_szLineBuff[0] != '>')			// if not descriptor then must be parsing sequence
		{
		pSrc = m_szLineBuff;
		pDst = &m_szSeqBuff[CurSeqLen];
		while(Chr = *pSrc++)
			{
			if(isspace(Chr))			// slough any whitespace within sequences
				continue;
			if(Chr == 'a' || Chr =='A' ||		// any nucleotide is accepted 
				Chr == 'c' || Chr == 'C' || 
				Chr == 'g' || Chr == 'G' || 
				Chr == 't' || Chr == 'T' || 
				Chr == 'u' || Chr == 'U')
				{
				*pDst++ = Chr;					
				CurNrun = 0;			// terminates any run of Ns
				CurSeqLen += 1;
				continue;
				}
			if(Chr == 'n' || Chr == 'N')	// check if run of Ns
				{
				CurNrun+=1;
				if(CurNrun <= MaxNrun)		// accept upto MaxNrun 
					{
					*pDst++ = Chr;
					CurSeqLen += 1;
					}
				continue;
				}
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unrecognised base %c at line number %d in file '%s'",Chr,LineNum,pszInFasta);
			Reset();
			return(-1);
			}
		continue;
		}

	if(CurSeqLen > 0)						// descriptor line, write any previously accepted sequence out
		{
		pDst = m_szSeqBuff;					// write out no more than 75 bases per line
		Chr = *pDst;
		do	{
			*pDst = Chr;
			pSrc = pDst;
			pDst += min(75,CurSeqLen);
			Chr = *pDst;
			*pDst = '\0';
			fputs(pSrc,pOutStream);
			fputs("\n",pOutStream);
			CurSeqLen -= min(75,CurSeqLen);
			}
		while(CurSeqLen > 0);
		}

		DescrLineLen = (int)strlen(m_szLineBuff);
		if(DescrLineLen > cMaxDescrLineLen)
			if(NumOverLenDescr++ < 10)
				gDiagnostics.DiagOut(eDLWarn,gszProcName,"Truncating over length descriptor at line %d in file '%s'",LineNum,pszInFasta);
		m_szLineBuff[cMaxDescrLineLen] = '\0';
		pSrc = m_szLineBuff;
		IdentLen = 0;
		while(IdentLen < cMaxDescrLineLen && (Chr = *pSrc++) && !isspace(Chr))
			szFastaDescr[IdentLen++]=Chr;
		szFastaDescr[IdentLen]='\0';
		if(IdentLen == 1)
			{
			if(NumEmptyDescrs++ < 10)
				gDiagnostics.DiagOut(eDLWarn,gszProcName,"Empty descriptor at line %d in file '%s'",LineNum,pszInFasta);
			strcpy(&szFastaDescr[1],pszDfltIdent);
			IdentLen = (int)strlen(szFastaDescr);
			}
		else
			{
			if(IdentLen > cMaxIdentNameLen)
				{
				if(NumOverLenIdents++ < 10)
					gDiagnostics.DiagOut(eDLWarn,gszProcName,"Truncating overlength identifer at line %d in file '%s'",LineNum,pszInFasta);
				szFastaDescr[cMaxIdentNameLen] = '\0';
				IdentLen = cMaxIdentNameLen;
				}
			}
		NumInstances = AddFastaIdentifier(&szFastaDescr[1]);
		if(NumInstances > 1)
			{
			if(NumDups < 10)
				gDiagnostics.DiagOut(eDLWarn,gszProcName,"Nonunique identifier (%s) at line %d in file '%s'",&szFastaDescr[1],LineNum,pszInFasta);
			NumDups += 1;
			sprintf(&szFastaDescr[IdentLen],"%s%d",pszSepUnique,NumInstances-1);
			IdentLen = (int)strlen(szFastaDescr);
			char *pDst;
			char *pSrc;
			pDst = &szFastaDescr[IdentLen];
			pSrc = &m_szLineBuff[1];
			while((Chr = *pSrc) && !isspace(Chr))
				pSrc++;
			if(Chr != '\0')
				while(*pDst++ = *pSrc++);
			m_szLineBuff[cMaxDescrLineLen] = '\0';
			fputs(szFastaDescr,pOutStream);
			fputs("\n",pOutStream);
			continue;
			}
		else
			{
			fputs(m_szLineBuff,pOutStream);				// write out original descriptor, no changes except perhaps excessive length truncation
			fputs("\n",pOutStream);
			}
	}

if(CurSeqLen > 0)						// write any previously accepted sequence out
	{
	pDst = m_szSeqBuff;					// write out no more than 75 bases per line
	Chr = *pDst;
	do	{
		*pDst = Chr;
		pSrc = pDst;
		pDst += min(75,CurSeqLen);
		Chr = *pDst;
		*pDst = '\0';
		fputs(pSrc,pOutStream);
		fputs("\n",pOutStream);
		CurSeqLen -= min(75,CurSeqLen);
		}
	while(CurSeqLen > 0);
	}

fclose(pInStream);
fclose(pOutStream);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing completed, replaced %d identifiers in file '%s'",NumDups,pszInFasta);
Reset();
return(0);
}

