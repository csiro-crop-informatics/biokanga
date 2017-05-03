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

const char *cpszProgVer = "0.9.0";		// increment with each release

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

const int cMaxInFileSpecs = 20; // can handle up to this many input files

int
Process(etPMode PMode,					// processing mode
		etFMode FMode,					// output format mode
		int LenNSeps,					// generate with this number of 'N' bases separating concatenated sequences 
		char *pszTrackTitle,			// track title for output UCSC BED
		int NumInputFileSpecs,		  	// number of input file specs
		char *pszInfileSpecs[],		  	// names of inputs files containing multifasta
		char *pszOutGenomeFile,			// output pseudo genome file
		char *pszOutGeneFile);			// output pseudo gene file

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
	return _T("kanga");
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
int iScreenLogLevel;		// level of file diagnostics
char szLogFile[_MAX_PATH];	// write diagnostics to this file
int Rslt = 0;   			// function result code >= 0 represents success, < 0 on failure
int Idx;

int LenNSeps;				// generate with this number of 'N' bases separating concatenated sequences 
etPMode PMode;				// processing mode
etFMode FMode;				// format output mode

char szTrackTitle[cMaxDatasetSpeciesChrom];		// track title if output format is UCSC BED

char szOutGenomeFile[_MAX_PATH];			// pseudo genome to this file
char szOutGeneFile[_MAX_PATH];			// pseudo genes to this file
int NumInputFiles;						// number of input files
char *pszInfileSpecs[cMaxInFileSpecs];  // input (wildcards allowed) multifasta files

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
struct arg_file *outgenefile = arg_file1("O","out","<file>",	"Output pseudo gene (BED) file");

struct arg_str  *title = arg_str0("t","title","<string>",       "track title");

struct arg_end *end = arg_end(20);

void *argtable[] = {help,version,FileLogLevel,LogFile,
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
		printf("\n%s the pseudo genome generator, Version %s\nOptions ---\n", gszProcName,cpszProgVer);
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
		printf("\n%s Version %s\n",gszProcName,cpszProgVer);
		exit(1);
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


#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	gStopWatch.Start();
	Rslt = Process(PMode,FMode,LenNSeps,szTrackTitle,NumInputFiles,pszInfileSpecs,szOutGenomeFile,szOutGeneFile);
	gStopWatch.Stop();
	Rslt = Rslt >=0 ? 0 : 1;
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Exit code: %d Total processing time: %s",Rslt,gStopWatch.Read());
	exit(Rslt);
	}
else
	{
	printf("\n%s the pseudo genome generator, Version %s\n",gszProcName,cpszProgVer);
	arg_print_errors(stdout,end,gszProcName);
	arg_print_syntax(stdout,argtable,"\nUse '-h' to view option and parameter usage\n");
	exit(1);
	}
return 0;
}


int m_hInFile;		// currently opened input file
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
Reset(void)
{
if(m_hInFile != -1)
	{
	close(m_hInFile);
	m_hInFile = -1;
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
Init(void)
{
m_hInFile = -1;
m_hOutGenomeFile = -1;
m_hOutGeneFile = -1;
m_pInFastaBuff = NULL;
m_pOutFastaBuff = NULL;
m_pOutBEDBuff = NULL;
Reset();
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
	Reset();
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
	Reset();
	return(Rslt);
	}

return(eBSFSuccess);
}

int
Process(etPMode PMode,					// processing mode
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

Init();

#ifdef _WIN32
if((m_hOutGenomeFile = open(pszOutGenomeFile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
if((m_hOutGenomeFile = open(pszOutGenomeFile,O_RDWR | O_CREAT |O_TRUNC, S_IREAD | S_IWRITE ))==-1)
#endif
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"OpenDumpFile:Unable to create or truncate  output file %s error: %s",pszOutGenomeFile,strerror(errno));
	Reset();
	return(false);
	}

#ifdef _WIN32
if((m_hOutGeneFile = open(pszOutGeneFile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
if((m_hOutGeneFile = open(pszOutGeneFile,O_RDWR | O_CREAT |O_TRUNC, S_IREAD | S_IWRITE ))==-1)
#endif
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: Unable to create or truncate  output file %s error: %s",pszOutGeneFile,strerror(errno));
	Reset();
	return(false);
	}

if((m_pInFastaBuff = new unsigned char [cAllocInFasta]) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessFastaFile:- Unable to allocate memory (%d bytes) for sequence buffer",cAllocInFasta);
	Reset();
	return(eBSFerrMem);
	}

if((m_pOutFastaBuff = new char [cAllocOutFasta]) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessFastaFile:- Unable to allocate memory (%d bytes) for sequence buffer",cAllocOutFasta);
	Reset();
	return(eBSFerrMem);
	}

if((m_pOutBEDBuff = new char [cAllocInBED]) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessFastaFile:- Unable to allocate memory (%d bytes) for sequence buffer",cAllocInBED);
	Reset();
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
		Reset();
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
			Reset();
			return(Rslt);
			}
		}
	}

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

Reset();
return(eBSFSuccess);
}
