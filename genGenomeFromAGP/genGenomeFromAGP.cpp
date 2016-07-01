// genGenomeFromAGP.cpp : Defines the entry point for the console application.
// This program uses a AGP file and source contig files to build a genome fasta file containing entries for each chromosome/region defined in the AGP file
// Initially to be used with contstruction of the assemblies for fusarium oxy genome ready for import into the UCSC browser

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

#include "Contigs.h"
#include "AGPs.h"

const char *cpszProgVer = "1.1.0";		// increment with each release

// processing modes
typedef enum TAG_ePMode {
	ePMdefault,					// default processing mode
	ePMSupercontig,				// generate supercontig to chrom mappings
	ePMplaceholder				// used to set the enumeration range
	} etPMode;

int Process(etPMode PMode,		// processing mode
			char *pszContigFile,	// input from these multifasta (wildcards allowed) contig files
			char *pszAGPFile,		// input from this AGP file
			char *pszAssembFile);	// output to this genome assembly file

int 
MapSupercontigs(etPMode PMode,		// processing mode
			char *pszSuperContigFile,	// input from this AGP file containing supercontig contig mappings
			char *pszContigFile,		// input from this AGP file containing chrom contig mappings
			char *pszOutFile);			// output to mapping file

CStopWatch gStopWatch;
CDiagnostics gDiagnostics;				// for writing diagnostics messages to log file
char gszProcName[_MAX_FNAME];			// process name

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
	return _T("genGenomeFromAGP");
}
// end of str library required code
#endif


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

int iFileLogLevel;			// level of file diagnostics
int iScreenLogLevel;			// level of file diagnostics
char szLogFile[_MAX_PATH];	// write diagnostics to this file
int Rslt = 0;   			// function result code >= 0 represents success, < 0 on failure

etPMode PMode;				// processing mode

char szContigFile[_MAX_PATH];		// input from these multifasta (wildcards allowed) contig files
char szAGPFile[_MAX_PATH];			// input from this AGP file
char szAssembFile[_MAX_PATH];		// output to this genome assembly file


// command line args
struct arg_lit  *help    = arg_lit0("h","help",                 "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to screen and logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"diagnostics log file");

struct arg_int *pmode = arg_int0("m","mode","<int>",		    "processing mode: 0 - standard, 1 - generate supercontig to chrom mappings");
struct arg_file *contigfile = arg_file1("i","in","<file>",		"input from these multifasta (wildcards allowed) contig files, or supercontig AGP file in mode 1");
struct arg_file *agpfile = arg_file1("I","in","<file>",			"input from this AGP file");
struct arg_file *assembfile = arg_file1("o","out","<file>",		"output to this multifasta assembly file");

struct arg_end *end = arg_end(20);

void *argtable[] = {help,version,FileLogLevel,LogFile,
					pmode,contigfile,agpfile,assembfile,
					end};

char **pAllArgs;
int argerrors;
argerrors = CUtility::arg_parsefromfile(argc,(char **)argv,&pAllArgs);
if(argerrors >= 0)
	argerrors = arg_parse(argerrors,pAllArgs,argtable);

/* special case: '--help' takes precedence over error reporting */
if (help->count > 0)
        {
		printf("\n%s AGP file processing, Version %s\nOptions ---\n", gszProcName,cpszProgVer);
        arg_print_syntax(stdout,argtable,"\n");
        arg_print_glossary(stdout,argtable,"  %-25s %s\n");
		printf("\nNote: Parameters can be entered into a parameter file, one parameter per line.");
		printf("\n      To invoke this parameter file then precede it's name with '@'");
		printf("\n      e.g. %s @myparams.txt\n",gszProcName);
		printf("\nPlease report any issues regarding usage of %s to stuart.stephen@csiro.au\n\n",gszProcName);
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
	if(FileLogLevel->count && !LogFile->count)
		{
		printf("\nError: FileLogLevel '-f%d' specified but no logfile '-F<logfile>'",FileLogLevel->ival[0]);
		exit(1);
		}

	iScreenLogLevel = iFileLogLevel = FileLogLevel->count ? FileLogLevel->ival[0] : eDLInfo;
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
			// now that command parameters have been parsed then initialise diagnostics log system
	if(!gDiagnostics.Open(szLogFile,(etDiagLevel)iScreenLogLevel,(etDiagLevel)iFileLogLevel,true))
		{
		printf("\nError: Unable to start diagnostics subsystem.");
		if(szLogFile[0] != '\0')
			printf(" Most likely cause is that logfile '%s' can't be opened/created",szLogFile);
		exit(1);
		}

	PMode = (etPMode)(pmode->count ? pmode->ival[0] : ePMdefault);
	if(PMode < ePMdefault || PMode >= ePMplaceholder)
		{
		printf("\nError: Processing mode '-m%d' specified outside of range %d..%d",PMode,ePMdefault,(int)ePMplaceholder-1);
		exit(1);
		}

	strcpy(szContigFile,contigfile->filename[0]);
	strcpy(szAGPFile,agpfile->filename[0]);
	strcpy(szAssembFile,assembfile->filename[0]);

// show user current resource limits
#ifndef _WIN32
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Resources: %s",CUtility::ReportResourceLimits());
#endif

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Version: %s Processing parameters:",cpszProgVer);

	const char *pszDescr;
	switch(PMode) {
		case ePMdefault:
			pszDescr = "standard genome AGP generation";
			break;
		case ePMSupercontig:
			pszDescr = "generate supercontig to chrom mappings";
			break;
		}
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Processing mode is : '%s'",pszDescr);

#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	gStopWatch.Start();
	Rslt = Process(PMode,szContigFile,szAGPFile,szAssembFile);
	gStopWatch.Stop();
	Rslt = Rslt >=0 ? 0 : 1;
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Exit code: %d Total processing time: %s",Rslt,gStopWatch.Read());
	exit(Rslt);
	}
else
	{
	printf("\n%s  AGP file processing, Version %s\n",gszProcName,cpszProgVer);
	arg_print_errors(stdout,end,gszProcName);
	arg_print_syntax(stdout,argtable,"\nUse '-h' to view option and parameter usage\n");
	exit(1);
	}
return 0;
}



typedef struct TAG_sChromSupercontig {
	int ChromSuperID;				// used to unquely identify this instance (1..m_NumChromSupers)
	char szChrom[cMaxOjCompLen];
	char szSupercontig[cMaxOjCompLen];
	UINT32 Start;
	UINT32 End;
	} tsChromSupercontig;

const int cAllocChromSupers = 10000;
int m_AllocdChromSupers;	// currently allocated number of chrom supercontig mappings
int m_NumChromSupers;		// number of chrom supercontigs used
tsChromSupercontig *m_pChromSupers;

const int cOutBuffSize = 0x0fffff;	// buffer output to this size before write to file


CContigs *m_pContigs;
CAGPs *m_pAGPs;
CAGPs *m_pSuperAGPs;
int m_hAssembFile;

int m_ColPsn;			// next output column position
int m_BuffOfs;			// offset in m_OutputBuff[] to next write
char m_OutputBuff[cOutBuffSize]; // used to buffer output



void
Reset(void)
{
if(m_hAssembFile != -1)
	{
	if(m_BuffOfs > 0)
		CUtility::SafeWrite(m_hAssembFile,m_OutputBuff,m_BuffOfs);
	close(m_hAssembFile);
	m_hAssembFile = -1;
	}
if(m_pContigs != NULL)
	{
	delete m_pContigs;
	m_pContigs = NULL;
	}
if(m_pAGPs != NULL)
	{
	delete m_pAGPs;
	m_pAGPs = NULL;
	}
if(m_pSuperAGPs != NULL)
	{
	delete m_pSuperAGPs;
	m_pSuperAGPs = NULL;
	}
if(m_pChromSupers != NULL)
	{
	delete m_pChromSupers;
	m_pChromSupers = NULL;
	}
m_ColPsn = 0;
m_BuffOfs = 0;
m_AllocdChromSupers = 0;
m_NumChromSupers = 0;
}

void
Init(void)
{
m_pContigs = NULL;
m_pAGPs = NULL;
m_pSuperAGPs = NULL;
m_pChromSupers = NULL;
m_hAssembFile = -1;
Reset();
}

bool
OutputNs(int Ncnt)
{
if(Ncnt <= 0 || m_hAssembFile == -1)
	return(false);
while(Ncnt--)
	{
	if((m_BuffOfs + 100) > cOutBuffSize)
		{
		CUtility::SafeWrite(m_hAssembFile,m_OutputBuff,m_BuffOfs);
		m_BuffOfs = 0;
		}
	if(m_ColPsn >= 80)
		{
		m_BuffOfs += sprintf(&m_OutputBuff[m_BuffOfs],"\n");
		m_ColPsn = 0;
		}
	
	m_OutputBuff[m_BuffOfs++] = 'N';
	m_ColPsn += 1;
	}
return(true);
}

bool
OutputChrom(char *pszChrom)
{
if(pszChrom == NULL || pszChrom[0] == '\0' || m_hAssembFile == -1)
	return(false);
if((m_BuffOfs + 1000) > cOutBuffSize)
	{
	CUtility::SafeWrite(m_hAssembFile,m_OutputBuff,m_BuffOfs);
	m_BuffOfs = 0;
	}
if(m_ColPsn != 0)
	{
	m_BuffOfs += sprintf(&m_OutputBuff[m_BuffOfs],"\n");
	m_ColPsn = 0;
	}
m_BuffOfs += sprintf(&m_OutputBuff[m_BuffOfs],">%s\n",pszChrom);
return(true);
}

bool
OutputSeq(int SeqLen,UINT8 *pSeq)
{
char Base;

if(SeqLen < 1 || pSeq == NULL || m_hAssembFile == -1)
	return(false);
if((m_BuffOfs + 1000) > cOutBuffSize)
	{
	CUtility::SafeWrite(m_hAssembFile,m_OutputBuff,m_BuffOfs);
	m_BuffOfs = 0;
	}

while(SeqLen--)
	{
	if((m_BuffOfs + 100) > cOutBuffSize)
		{
		CUtility::SafeWrite(m_hAssembFile,m_OutputBuff,m_BuffOfs);
		m_BuffOfs = 0;
		}
	if(m_ColPsn >= 80)
		{
		m_BuffOfs += sprintf(&m_OutputBuff[m_BuffOfs],"\n");
		m_ColPsn = 0;
		}
	switch(*pSeq++) {
		case eBaseA:
			Base = 'a';
			break;
		case eBaseC:
			Base = 'c';
			break;
		case eBaseG:
			Base = 'g';
			break;
		case eBaseT:
			Base = 't';
			break;
		default:			// treat any other base as being indeterminate
			Base = 'n';
			break;
		}
	
	m_OutputBuff[m_BuffOfs++] = Base;
	m_ColPsn += 1;
	}
return(true);
}


int 
Process(etPMode PMode,		// processing mode
			char *pszContigFile,	// input from these multifasta (wildcards allowed) contig files
			char *pszAGPFile,		// input from this AGP file
			char *pszAssembFile)	// output to this genome assembly file
{
int Rslt;
tsAGPentry *pAGPentry;
Init();
if(PMode == ePMSupercontig)
	return(MapSupercontigs(PMode,pszContigFile,pszAGPFile,pszAssembFile));

m_pAGPs = new CAGPs;

// load and parse the AGP file
if((Rslt = m_pAGPs->LoadAGPs(pszAGPFile)) != eBSFSuccess)
	{
	Reset();
	return(Rslt);
	}

// load and parse the contig containing fasta file(s)
m_pContigs = new CContigs;
if((Rslt = m_pContigs->LoadContigs(pszContigFile)) != eBSFSuccess)
	{
	Reset();
	return(Rslt);
	}

// create file to write assembly sequences out to
#ifdef _WIN32
if((m_hAssembFile = open(pszAssembFile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
if((m_hAssembFile = open(pszAssembFile, O_RDWR | O_CREAT,S_IREAD | S_IWRITE))==-1)
#endif
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create %s - %s",pszAssembFile,strerror(errno));
	Reset();
	return(eBSFerrCreateFile);
	}

// now check that all contigs referenced in the AGP file have been loaded from the multifasta contig files
pAGPentry = NULL;
int NbaseCnt;
UINT8 *pSeq;
UINT32 Start = 1;
char szCurObject[cMaxOjCompLen];
szCurObject[0] = '\0';
while((pAGPentry = m_pAGPs->Next(pAGPentry))!=NULL)
	{
	if(szCurObject[0] == '\0' || stricmp(pAGPentry->szObjIdent,szCurObject))
		{
		strcpy(szCurObject,pAGPentry->szObjIdent);				// new chrom starting
		OutputChrom(szCurObject);
		Start = 1;
		if(pAGPentry->Start > Start)
			{
			NbaseCnt = pAGPentry->Start - Start;
			OutputNs(NbaseCnt);			
			Start = pAGPentry->Start;
			}
		}
	if(pAGPentry->Type == eAGPSSN)
		{
		NbaseCnt = 1 + pAGPentry->End - pAGPentry->Start;
		OutputNs(NbaseCnt);			
		Start = pAGPentry->End + 1;
		}
	else
		{
		if((pSeq = m_pContigs->LocateSeq(pAGPentry->Comp.szCompID,pAGPentry->Comp.Start))==NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: Unable to locate sequence for contig named '%s'",pAGPentry->Comp.szCompID);
			Reset();
			return(eBSFerrFastaDescr);
			}
		OutputSeq(1 + pAGPentry->Comp.End - pAGPentry->Comp.Start,pSeq);
		Start = pAGPentry->End + 1;
		}
	}
if(m_hAssembFile != -1 && m_BuffOfs)
	{
	CUtility::SafeWrite(m_hAssembFile,m_OutputBuff,m_BuffOfs);
	m_BuffOfs = 0;
	}
close(m_hAssembFile);
m_hAssembFile = -1;
Reset();
return(0);
}



int 
MapSupercontigs(etPMode PMode,		// processing mode
			char *pszSuperContigFile,	// input from this AGP file containing supercontig contig mappings
			char *pszContigFile,		// input from this AGP file containing chrom contig mappings
			char *pszOutFile)			// output to mapping file
{
int Rslt;
tsAGPentry *pAGPentry;
tsAGPentry *pSupercontig;
tsChromSupercontig *pChromSuper;
tsChromSupercontig *pCurChromSuper;
int Idx;

Init();

// load and parse the AGP files
m_pAGPs = new CAGPs;
if((Rslt = m_pAGPs->LoadAGPs(pszContigFile)) != eBSFSuccess)
	{
	Reset();
	return(Rslt);
	}
m_pSuperAGPs = new CAGPs;
if((Rslt =  m_pSuperAGPs->LoadAGPs(pszSuperContigFile)) != eBSFSuccess)
	{
	Reset();
	return(Rslt);
	}
if((m_pChromSupers = new tsChromSupercontig[cAllocChromSupers])==NULL)
	{
	Reset();
	return(eBSFerrMem);
	}
m_AllocdChromSupers = cAllocChromSupers;
m_NumChromSupers = 0;

// create file to write chrom supercontig mappings out to
#ifdef _WIN32
if((m_hAssembFile = open(pszOutFile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
if((m_hAssembFile = open(pszOutFile, O_RDWR | O_CREAT,S_IREAD | S_IWRITE))==-1)
#endif
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create %s - %s",pszOutFile,strerror(errno));
	Reset();
	return(eBSFerrCreateFile);
	}

// iterate over each contig, and get it's chrom and start/end
// now check that all contigs referenced in the AGP file have been loaded from the multifasta contig files
pCurChromSuper = NULL;
pAGPentry = NULL;
pSupercontig = NULL;
char szCurObject[cMaxOjCompLen];
szCurObject[0] = '\0';
while((pAGPentry = m_pAGPs->Next(pAGPentry))!=NULL)
	{
	if(szCurObject[0] == '\0' || stricmp(pAGPentry->szObjIdent,szCurObject))
		strcpy(szCurObject,pAGPentry->szObjIdent);				// new chrom starting
	if(pAGPentry->Type == eAGPSSN)
		continue;
	else
		{
		// locate contig in supercontigs
		if((pSupercontig = m_pSuperAGPs->LocateCompID(pAGPentry->Comp.szCompID))==NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: Unable to locate supercontig for contig named '%s'",pAGPentry->Comp.szCompID);
			Reset();
			return(eBSFerrFastaDescr);
			}

		if(pCurChromSuper == NULL || stricmp(pSupercontig->szObjIdent,pCurChromSuper->szSupercontig))
			{
			Idx = 0;

			if(m_NumChromSupers > 0)
				{
				pChromSuper = m_pChromSupers;
				for(Idx = 0; Idx < m_NumChromSupers; Idx++,pChromSuper++)
					{
					if(!stricmp(pSupercontig->szObjIdent,pChromSuper->szSupercontig))
						break;
					}
				}
			if(Idx == m_NumChromSupers)					// new supercontig starting
				{
				pChromSuper = &m_pChromSupers[Idx];
				memset(pChromSuper,0,sizeof(tsChromSupercontig));
				strcpy(pChromSuper->szChrom,pAGPentry->szObjIdent);
				strcpy(pChromSuper->szSupercontig,pSupercontig->szObjIdent);
				pChromSuper->Start = pAGPentry->Start;
				pChromSuper->End = pAGPentry->End;
				pChromSuper->ChromSuperID = ++m_NumChromSupers;
				pCurChromSuper = pChromSuper;
				continue;
				}
			}
		else
			pChromSuper = pCurChromSuper;
		if(pAGPentry->Start < pChromSuper->Start)
			pChromSuper->Start = pAGPentry->Start;
		if(pAGPentry->End > pChromSuper->End)
			pChromSuper->End = pAGPentry->End;
		}
	}
pChromSuper = m_pChromSupers;
m_BuffOfs = 0;
for(Idx = 0; Idx < m_NumChromSupers; Idx++,pChromSuper++)
	{
	m_BuffOfs += sprintf(&m_OutputBuff[m_BuffOfs],"%s,%s,%d,%d\n",pChromSuper->szSupercontig,pChromSuper->szChrom,pChromSuper->Start,pChromSuper->End);
	if((m_BuffOfs + 500) > cOutBuffSize)
		{
		CUtility::SafeWrite(m_hAssembFile,m_OutputBuff,m_BuffOfs);
		m_BuffOfs = 0;
		}
	}
if(m_hAssembFile != -1 && m_BuffOfs)
	{
	CUtility::SafeWrite(m_hAssembFile,m_OutputBuff,m_BuffOfs);
	m_BuffOfs = 0;
	}
close(m_hAssembFile);
m_hAssembFile = -1;
Reset();
return(0);
}


