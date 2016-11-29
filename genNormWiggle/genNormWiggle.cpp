// genNormWiggle.cpp : Defines the entry point for the console application.
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

const char *cpszProgVer = "0.0.3";		// increment with each release

const int cDfltKMerLen = 5;			// default K-mer length
const int cMaxKMerLen = 13;			// max K-mer length (restricted due to memory requirements)

const int cMaxDistLog2 = 20;		// can bin cnt distributions with max of 2^cMaxDistLog2

const int cChromSeqReAlloc = 5000000;	// realloc chrom sequence size

// processing modes
typedef enum TAG_ePMode {		
	ePMdefault,					// default is to process for read starts only
	ePMcoverage,				// or processing can be for read coverage
	ePMplaceholder				// used to set the enumeration range
	} etPMode;

// input format mode
typedef enum TAG_eIFormat {		
	eIFdefault,					// UCSC BED format
	eIFCSV,						//  CSV loci 
	eIFplaceholder				// used to set the enumeration range
	} etIFormat;

int
Process(etPMode PMode,					// processing mode
		etIFormat IFormat,				// input aligned reads file format 
		char Strand,					// process for this strand only
		int KMerLen,					// process for this length K-Mers
		int KMerInstsThres,				// report KMers with instances >= this threshold
		char *pszInFile,				// input BED file or csv file containing aligned read loci
		char *pszProcElsFile,			// optional input BED file containing elements over which K-mer deviations to be reported
		char *pszInBioseqFile,			// bioseq genome file
		char *pszRsltsFile,				// optional output UCSC wiggle file
	    char *pszKMerRsltsFile);		// output K-mer distributions file

char *Conf2Txt(teOctStructStats ConfParam); // returns short description of specified conformational characteristic

CStopWatch gStopWatch;
CDiagnostics gDiagnostics;				// for writing diagnostics messages to log file
char gszProcName[_MAX_FNAME];			// process name

#pragma pack(1)
typedef struct TAG_sKMerDist {
	UINT32 Insts;		// number of instances of this K-mer with more than 1 count
	UINT32 Cnts;		// total number of read counts observed for this K-mer
	UINT32 Dist[cMaxDistLog2];	// distribution of instance counts
} tsKMerDist;

// element distributions 
const int cMaxKMerInstCnts = 25;	   // allow upto this many K-mer instances in any element

typedef struct TAG_sKMerElDist {
	int FeatID;			// element identifier from the biobed file
	UINT64 KMerID;		// identifies this specific K-mer (UINT64 allows for K <= 31)
	UINT32 Insts;		// number of instances of this K-mer with at least 1 count within element
	UINT32 Cnts;		// total number of read counts observed for this K-mer within element
	UINT32 Min;			// minimum count observed for this K-mer within element
	UINT32 Max;			// maximum count observed for this K-mer within element
	double Mean;		// mean number of counts
	double Variance;	// count variance
	double StdDev;		// stddev of this K-mer within element
	struct {
		int ChromOfs;	// chrom offset at which K-mer starts
		int Cnt;		// count for this K-mer instance
		} Dist[cMaxKMerInstCnts];	// to hold K-mer cnts for each K-mer instance within this element when read start processing
} tsKMerElDist;

#pragma pack()

tsKMerDist *m_pKMerDist;
size_t m_KMerDists;

int m_hRsltsFile;					// handle for opened results file

CBEDfile *m_pBEDElFile;				// to hold instantiated BED element file
CBioSeqFile *m_pBioSeqFile;			// holds instantiated genome assembly sequences
etSeqBase *m_pChromSeq;				// holds current assembly chromosome sequence
int m_AllocdChromSeq;				// allocated assembly chromosome sequence length
int m_ChromSeqLen;					// current assembly chromosome sequence length
char m_szCurChrom[cMaxDatasetSpeciesChrom]; // for this chromosome
int m_CurChromID;					// for this chromosome

const int cAllocKMerElDistInsts = 10000;	// realloc in this many instance increments

int m_NumKMerElDistInsts;
int m_AllocdKMerElDistInsts;				// currently allocated for this many instances
size_t m_AllocdKMerElDistMem; 
tsKMerElDist *m_pKMerElDist;

CCSVFile *m_pCSVFile;			// used if processing input CSV file for coverage
CBEDfile *m_pBEDFile;			// used if processing input BED files for coverage

const int cMaxChromCov = 200;			// can handle at most this many chromosomes
const int cAllocCovCnts = 0x07fffff;	// allocate for chrom coverage cnts in this sized increments

typedef struct TAG_sChromCnts {
	char szChrom[cMaxDatasetSpeciesChrom+1];	// coverage is on this chromosome
	int AllocCovCnts;							// allocated to hold cnts for this sized chromosome
	int StartOfs;								// pCovCnts[offset] of first coverage cnt
	int EndOfs;									// pCovCnts[offset] of last coverage cnt
	UINT16 *pCovCnts;							// coverage counts
} tsChromCnts;

tsChromCnts m_ChromCnts[cMaxChromCov];
int m_NumChromsCov = 0;		

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
	return _T("genNormWiggle");
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
int iScreenLogLevel;		// level of screen diagnostics
int iFileLogLevel;			// level of file diagnostics
char szLogFile[_MAX_PATH];	// write diagnostics to this file

etPMode PMode;				// processing mode
etIFormat IFormat;			// format of input reads alignment file
int Rslt = 0;   			// function result code >= 0 represents success, < 0 on failure
int Strand;					// filter for this strand

char szRsltsFile[_MAX_PATH];	// output wiggle file
char szProcElsFile[_MAX_PATH];		// retain only loci covered by elements in this BED file
char szKMerRsltsFile[_MAX_PATH]; // output K-mer distribution file
char szInFile[_MAX_PATH];		// input file 
char szGenomeFile[_MAX_PATH];	// bioseq genome file

int KMerLen;				// process for this length K-mer
int KMerInstThres;			// report K-mers with instance counts >= this threshold

// command line args
struct arg_lit  *help    = arg_lit0("hH","help",                "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"diagnostics log file");
struct arg_int  *strand = arg_int0("s","strand","<int>",		"filter for this strand: 0 - any, 1 - Watson '+', 2 - Crick '-' (default is any)");
struct arg_int *pmode = arg_int0("m","mode","<int>",		    "processing mode: 0 - process read starts, 1 - process coverage (default = 0)");
struct arg_int *iformat = arg_int0("M","mode","<int>",		    "input loci file format: 0 - input BED, 1 - input CSV (default = 0)");
struct arg_file *infile = arg_file1("i","in","<file>",			"input from this CSV loci or BED alignment file");
struct arg_file *gfile = arg_file1("I","seq","<file>",			"input bioseq genome file");
struct arg_file *procelsfile = arg_file0("p","procels","<file>","optional input BED file containing elements over which K-mer stddev are to be reported");
struct arg_file *outfile = arg_file0("O","out","<file>",		"output optional UCSC wiggle format to this file");
struct arg_file *outkmerfile = arg_file1("o","out","<file>",	"output K-mer distributions to this CSV file");
struct arg_int *kmerlen = arg_int0("l","kmerlen","<int>",		"process for this length K-mer (default is 7)");
struct arg_int *kmerinstthres = arg_int0("L","kmerinstthres","<int>",	"report K-mers having instances >= this threshold (default is 50)");

struct arg_end *end = arg_end(20);

void *argtable[] = {help,version,FileLogLevel,LogFile,
					pmode,iformat,strand,kmerlen,kmerinstthres,infile,procelsfile,gfile,outfile,outkmerfile,
					end};

char **pAllArgs;
int argerrors;
argerrors = CUtility::arg_parsefromfile(argc,(char **)argv,&pAllArgs);
if(argerrors >= 0)
	argerrors = arg_parse(argerrors,pAllArgs,argtable);

    /* special case: '--help' takes precedence over error reporting */
if (help->count > 0)
        {
		printf("\n%s a component of the K-mer Adaptive Next Generation Aligner toolset, Version %s\nOptions ---\n", gszProcName,cpszProgVer);
        arg_print_syntax(stdout,argtable,"\n");
        arg_print_glossary(stdout,argtable,"  %-25s %s\n");
		printf("\nNote: Parameters can be entered into a parameter file, one parameter per line.");
		printf("\n      To invoke this parameter file then precede it's name with '@'");
		printf("\n      e.g. %s @myparams.txt\n",gszProcName);
		printf("\nPlease report any issues regarding usage of %s to stuart.stephen@csiro.au\n\n",gszProcName);
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

	PMode = (etPMode)(pmode->count ? pmode->ival[0] : 0);
	if(PMode < 0 || PMode >= ePMplaceholder)
		{
		printf("\nError: Processing mode '-m%d' specified outside of range %d..%d",PMode,0,(int)ePMplaceholder-1);
		exit(1);
		}

	IFormat = (etIFormat)(iformat->count ? iformat->ival[0] : 0);
	if(IFormat < 0 || IFormat >= eIFplaceholder)
		{
		printf("\nError: Input format '-M%d' specified outside of range %d..%d",IFormat,0,(int)eIFplaceholder-1);
		exit(1);
		}

	if(PMode == ePMdefault)
		{
		Strand = (int)'+';
		printf("\nWarning: forcing strand to be '+' only because read starts are being processed");
		}
	else
		{
		Strand = strand->count ? strand->ival[0] : 0;
		if(Strand < 0 || Strand > 2)
			{
			printf("\nError: Strand '-s%d' specified outside of range 0..2",Strand);
			exit(1);
			}

		switch(Strand) {
			case 1: Strand = (int)'+'; break;
			case 2: Strand = (int)'-'; break;
			case 0: Strand = (int)'*'; break;
			}
		}


	KMerLen = kmerlen->count ? kmerlen->ival[0] : cDfltKMerLen;
	if(KMerLen < 1 || KMerLen > cMaxKMerLen)
		{
		printf("\nError: KMerLen '-l%d' specified outside of range 1..%d",KMerLen,cMaxKMerLen);
		exit(1);
		}

	KMerInstThres = kmerinstthres->count ? kmerinstthres->ival[0] : 50;
	if(KMerInstThres < 0 || KMerInstThres > 100000)
		{
		printf("\nError: KMerInstThres '-L%d' specified outside of range 0..%d",KMerInstThres,100000);
		exit(1);
		}

	strcpy(szInFile,infile->filename[0]);
	if(procelsfile->count)
		strcpy(szProcElsFile,procelsfile->filename[0]);
	else
		szProcElsFile[0] = '\0';
	if(outfile->count)
		strcpy(szRsltsFile,outfile->filename[0]);
	else
		szRsltsFile[0] = '\0';
	strcpy(szKMerRsltsFile,outkmerfile->filename[0]);
	strcpy(szGenomeFile,gfile->filename[0]);

		// now that command parameters have been parsed then initialise diagnostics log system
	if(!gDiagnostics.Open(szLogFile,(etDiagLevel)iScreenLogLevel,(etDiagLevel)iFileLogLevel,true))
		{
		printf("\nError: Unable to start diagnostics subsystem.");
		if(szLogFile[0] != '\0')
			printf(" Most likely cause is that logfile '%s' can't be opened/created",szLogFile);
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
			pszDescr = "read starts";
			break;
		case ePMcoverage:			
			pszDescr = "read coverage";
			break;
		}
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"processing mode is : '%s'",pszDescr);

	switch(IFormat) {
		case eIFdefault:		// UCSC BED format
			pszDescr =  "UCSC BED format";
			break;
		case eIFCSV:			// 
			pszDescr = "CSV loci";
			break;
		}
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"input alignment file format : '%s'",pszDescr);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"process for this strand only: '%c'",(char)Strand);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"input %s file: '%s'",pszDescr, szInFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"input bioseq genome file: '%s'",szGenomeFile);
	if(szProcElsFile[0] != '\0')
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"BED file containing elements over which K-mer deviations to be reported : '%s'",szProcElsFile);
	if(szRsltsFile[0] != '\0')
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"output UCSC Wiggle format file: '%s'",szRsltsFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"output K-mer distribution file: '%s'",szKMerRsltsFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"process for these length K-mers: %d",KMerLen);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"report K-mers with at least this number of instances: %d",KMerInstThres);
#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	gStopWatch.Start();
	// processing here...
	Rslt = Process(PMode,IFormat,Strand,KMerLen,KMerInstThres,szInFile,szProcElsFile,szGenomeFile,szRsltsFile,szKMerRsltsFile);


	gStopWatch.Stop();
	Rslt = Rslt >=0 ? 0 : 1;
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Exit code: %d Total processing time: %s",Rslt,gStopWatch.Read());
	exit(Rslt);
	}
else
	{
	printf("\n%s component of the the K-mer Adaptive Next Generation Aligner toolkit, Version %s\n",gszProcName,cpszProgVer);
	arg_print_errors(stdout,end,gszProcName);
	arg_print_syntax(stdout,argtable,"\nUse '-h' to view option and parameter usage\n");
	exit(1);
	}
}


void
Init(void)
{
m_pBioSeqFile = NULL;		// holds instantiated genome assembly sequences
m_pChromSeq = NULL;			// holds current assembly chromosome sequence
m_pKMerDist = NULL;
m_KMerDists = 0;
m_pCSVFile = NULL;
m_pBEDFile = NULL;
m_pBEDElFile = NULL;
m_pKMerElDist = NULL;

m_hRsltsFile = -1;			// handle for opened results file
m_AllocdChromSeq = 0;		// allocated assembly chromosome sequence length
m_ChromSeqLen = 0;			// current assembly chromosome sequence length
m_CurChromID = 0;			// for this chromosome
m_szCurChrom[0] = '\0';

for(int ChromIdx = 0; ChromIdx < cMaxChromCov; ChromIdx++)
	{
	m_ChromCnts[ChromIdx].AllocCovCnts = 0;
	m_ChromCnts[ChromIdx].StartOfs = 0;
	m_ChromCnts[ChromIdx].EndOfs = 0;
	m_ChromCnts[ChromIdx].pCovCnts = NULL;
	m_ChromCnts[ChromIdx].szChrom[0] = '\0';
	}
m_NumChromsCov = 0;	

}

void
Reset(void)
{
if(m_hRsltsFile != -1)
	{
	close(m_hRsltsFile);
	m_hRsltsFile = -1;
	}
if(m_pBioSeqFile != NULL)
	{
	delete m_pBioSeqFile;
	m_pBioSeqFile = NULL;
	}

if(m_pChromSeq != NULL)
	{
	delete m_pChromSeq;
	m_pChromSeq = NULL;
	}

if(m_pKMerDist != NULL)
	{
    delete m_pKMerDist;
	m_pKMerDist = NULL;
	}
m_KMerDists = 0;
if(m_pBEDFile != NULL)
	{
	delete m_pBEDFile;
	m_pBEDFile = NULL;
	}
if(m_pCSVFile != NULL)
	{
	delete m_pCSVFile;
	m_pCSVFile = NULL;
	}
if(m_pBEDElFile != NULL)
	{
	delete m_pBEDElFile;
	m_pBEDElFile = NULL;
	}

if(m_pKMerElDist != NULL)
	{
#ifdef _WIN32
	free(m_pKMerElDist);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pKMerElDist != MAP_FAILED)
		munmap(m_pKMerElDist,m_AllocdKMerElDistMem);
#endif
	m_pKMerElDist = NULL;
	}
m_NumKMerElDistInsts = 0;
m_AllocdKMerElDistInsts = 0;
m_AllocdKMerElDistMem = 0;
m_AllocdChromSeq = 0;				// allocated assembly chromosome sequence length
m_ChromSeqLen = 0;					// current assembly chromosome sequence length
m_CurChromID = 0;						// for this chromosome
m_szCurChrom[0] = '\0';
if(m_NumChromsCov > 0)
	{
	for(int ChromIdx = 0; ChromIdx < m_NumChromsCov; ChromIdx++)
		{
		m_ChromCnts[ChromIdx].StartOfs = 0;
		m_ChromCnts[ChromIdx].EndOfs = 0;
		if(m_ChromCnts[ChromIdx].pCovCnts != NULL)
			{
#ifdef _WIN32
			free(m_ChromCnts[ChromIdx].pCovCnts);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
			if(m_ChromCnts[ChromIdx].pCovCnts != MAP_FAILED)
				munmap(m_ChromCnts[ChromIdx].pCovCnts,m_ChromCnts[ChromIdx].AllocCovCnts * sizeof(UINT16));
#endif
			m_ChromCnts[ChromIdx].pCovCnts = NULL;
			}
		m_ChromCnts[ChromIdx].AllocCovCnts = 0;
		m_ChromCnts[ChromIdx].szChrom[0] = '\0';
		}
	m_NumChromsCov = 0;	
	}
}

// GenSeqIdx
// Generates a sequence index for specified sequence of length SeqLen
// Max length sequence is 31
// Returns -1 if sequence contains any indeterminate or non a,c,g,t bases
INT64
GenSeqIdx(int SeqLen,etSeqBase *pSeq)
{
int Idx;
UINT64 SeqIdx;
int Base;
if(SeqLen > 31)
	return(-1);
SeqIdx=0;
for(Idx=0; Idx < SeqLen; Idx++,pSeq++)
	{
	Base = *pSeq & ~cRptMskFlg;
	if(Base > eBaseT)
		return(-1);
	SeqIdx <<= 2;
	SeqIdx |= Base;
	}
return(SeqIdx);
}

char *
StepIdx2Seq(int SeqLen,INT64 SeqIdx)
{
static char szSeqBuff[256];
char *pChr;
int Base;
int Idx;

szSeqBuff[SeqLen] = '\0';
pChr = &szSeqBuff[SeqLen-1];
for(Idx=0; Idx < SeqLen; Idx++,pChr--)
	{
	Base = SeqIdx & 0x03;
	switch(Base) {
		case 0: *pChr = 'a'; break;
		case 1: *pChr = 'c'; break;
		case 2: *pChr = 'g'; break;
		case 3: *pChr = 't'; break;
		}
	SeqIdx >>= 2;
	}
return(szSeqBuff);
}


int
LoadNxtChrom(void)
{
int Rslt;

if(m_pBioSeqFile == NULL)
	return(0);
if((m_CurChromID = m_pBioSeqFile->Next(m_CurChromID)) < 1)
	{
	m_CurChromID = 0;						
	m_szCurChrom[0] = '\0';
	return(eBSFerrEntry);
	}

m_pBioSeqFile->GetName(m_CurChromID,sizeof(m_szCurChrom),m_szCurChrom);
m_ChromSeqLen = m_pBioSeqFile->GetDataLen(m_CurChromID);
if(m_pChromSeq == NULL || m_ChromSeqLen >= m_AllocdChromSeq)
	{
	if(m_pChromSeq != NULL)
		{
		delete m_pChromSeq;
		m_pChromSeq = NULL;
		m_AllocdChromSeq = 0;
		}
	if((m_pChromSeq = new etSeqBase [m_ChromSeqLen + cChromSeqReAlloc]) == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate memory for chromosome sequence of %d bytes",m_ChromSeqLen + cChromSeqReAlloc);
		return(eBSFerrMem);
		}
	m_AllocdChromSeq = m_ChromSeqLen + cChromSeqReAlloc;
	}
if((Rslt=m_pBioSeqFile->GetData(m_CurChromID,eSeqBaseType,0,m_pChromSeq,m_ChromSeqLen))!=m_ChromSeqLen)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to load chromosome sequence for '%s' of expected length",m_szCurChrom,m_ChromSeqLen);
	return(Rslt);
	}
return(m_ChromSeqLen);
}


int
BuildReadCoverage(bool bStartOnly,		// if true then only process for read starts, not complete read overlaps
				char *pszChrom,			// coverage is onto this chrom
			  int StartOfs,				// coverage start at this offset 
			  int EndOfs,				// and ends at this offset inclusive
			  int Cnt)					// increment coverage by this
{
tsChromCnts *pChrom;
int ChromIdx;
int AllocCovCnts;
UINT16 *pCovCnts;
size_t ReallocTo;

if(pszChrom == NULL || pszChrom[0] == '\0')
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"BuildReadCoverage: No chromosome specified");
	Reset();
	return(eBSFerrChrom);
	}

// arbitary count clamping to be in range 1..20000 inclusive
if(Cnt < 1)
	Cnt = 1;
else
	if(Cnt > 20000)
		Cnt = 20000;

// ensure StartOfs and EndOfs are both >= 0
if(StartOfs < 0)
	StartOfs = 0;
if(EndOfs < 0)
	EndOfs = 0;

// ensure StartOfs <= EndOfs
if(StartOfs > EndOfs)
	{
	int TmpOfs = EndOfs;
	EndOfs = StartOfs;
	StartOfs = TmpOfs;
	}

// check if this is a new chrom or if coverage is onto an existing chrom
pChrom = &m_ChromCnts[0];
ChromIdx = 0;
if(m_NumChromsCov > 0)
	{
	for(ChromIdx = 0; ChromIdx < m_NumChromsCov; ChromIdx++,pChrom++)
		{
		if(!strnicmp(pszChrom,pChrom->szChrom,cMaxDatasetSpeciesChrom))
			break;
		}
	}
if(ChromIdx == m_NumChromsCov)	// if a new or first chrom
	{
	strncpy(pChrom->szChrom,pszChrom,cMaxDatasetSpeciesChrom);
	pChrom->szChrom[cMaxDatasetSpeciesChrom-1] = 0;
	pChrom->StartOfs = StartOfs;
	pChrom->EndOfs = EndOfs;
	AllocCovCnts = EndOfs + cAllocCovCnts;
	ReallocTo =  AllocCovCnts * sizeof(UINT16);
#ifdef _WIN32
	pChrom->pCovCnts = (UINT16 *) malloc(ReallocTo);	// initial and perhaps the only allocation
	if(pChrom->pCovCnts == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"BuildReadCoverage: Memory allocation of %lld bytes - %s",(INT64)ReallocTo,strerror(errno));
		Reset();
		return(eBSFerrMem);
		}
#else
	// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
	pChrom->pCovCnts = (UINT16 *)mmap(NULL,ReallocTo, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
	if(pChrom->pCovCnts == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"BuildReadCoverage: Memory allocation of %lld bytes through mmap()  failed - %s",(INT64)ReallocTo,strerror(errno));
		pChrom->pCovCnts = NULL;
		Reset();
		return(eBSFerrMem);
		}
#endif
	pChrom->AllocCovCnts = AllocCovCnts;
	memset(pChrom->pCovCnts,0,ReallocTo);
	m_NumChromsCov += 1;
	}

// check if chrom coverage cnts needs to be extended
if(EndOfs >= pChrom->AllocCovCnts)
	{
	AllocCovCnts = EndOfs + cAllocCovCnts;
	ReallocTo = AllocCovCnts * sizeof(UINT16);
#ifdef _WIN32
	pCovCnts = (UINT16 *) realloc(pChrom->pCovCnts,ReallocTo);
#else
	pCovCnts = (UINT16 *)mremap(pChrom->pCovCnts,pChrom->AllocCovCnts * sizeof(UINT16),ReallocTo,MREMAP_MAYMOVE);
	if(pCovCnts == MAP_FAILED)
		pCovCnts = NULL;
#endif
	if(pCovCnts == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"BuildReadCoverage: Memory re-allocation to %d bytes - %s",ReallocTo,strerror(errno));
		Reset();
		return(eBSFerrMem);
		}
	pChrom->pCovCnts = pCovCnts;
	memset(&pChrom->pCovCnts[pChrom->AllocCovCnts],0,(AllocCovCnts - pChrom->AllocCovCnts) * sizeof(UINT16));
	pChrom->AllocCovCnts = AllocCovCnts;
	}

if(EndOfs > pChrom->EndOfs)
	pChrom->EndOfs = EndOfs;
if(StartOfs < pChrom->StartOfs)
	pChrom->StartOfs = StartOfs;

pCovCnts = &pChrom->pCovCnts[StartOfs];

if(bStartOnly)
	{
	if((0x0fffe - *pCovCnts) > Cnt)
		*pCovCnts += Cnt;
	else
		*pCovCnts = 0x0fffe;
	}
else
	{
	while(StartOfs++ <= EndOfs)
		{
		// clamp accumulated cnts to be no greater than 0x0fffe
		if((0x0fffe - *pCovCnts) > Cnt)
			*pCovCnts += Cnt;
		else
			*pCovCnts = 0x0fffe;
		pCovCnts += 1;
		}
	}
return(eBSFSuccess);
}

int
WriteReadsWig(char *pszSrcFile,char *pszRsltsFile)
{
int BuffIdx;
char szLineBuff[8096];
tsChromCnts *pChrom;
UINT16 *pCnts;
int ChromIdx;
int SeqIdx;
bool bStartRegion;
int EmptyRegionLen;

#ifdef _WIN32
if((m_hRsltsFile = open(pszRsltsFile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
if((m_hRsltsFile = open(pszRsltsFile, O_RDWR | O_CREAT | O_TRUNC, S_IREAD | S_IWRITE ))==-1)
#endif
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create UCSC Wiggle file: %s - %s",pszRsltsFile,strerror(errno));
	Reset();
	return(eBSFerrCreateFile);
	}

BuffIdx = sprintf(szLineBuff,"track type=wiggle_0 color=50,150,255 autoScale=off maxHeightPixels=128:32:8 name=\"Reads - %s\" description=\"Reads distribution for %s\"\n",pszSrcFile,pszSrcFile);
CUtility::SafeWrite(m_hRsltsFile,szLineBuff,BuffIdx);
BuffIdx = 0;
pChrom = &m_ChromCnts[0];
for(ChromIdx = 0 ; ChromIdx < m_NumChromsCov; ChromIdx++,pChrom++)
	{
	EmptyRegionLen = 0;
	bStartRegion = true;
	pCnts = &pChrom->pCovCnts[pChrom->StartOfs];
	for(SeqIdx = pChrom->StartOfs; SeqIdx < pChrom->EndOfs; SeqIdx++,pCnts++)
		{
		if(*pCnts == 0)
			{
			EmptyRegionLen += 1;
			if(EmptyRegionLen > 10)
				bStartRegion = true;
			continue;
			}
		else
			{
			if(!bStartRegion && EmptyRegionLen > 0 && EmptyRegionLen <= 10)
				{
				while(EmptyRegionLen--)
					BuffIdx += sprintf(&szLineBuff[BuffIdx],"0\n");
				}
			EmptyRegionLen = 0;
			}
		if(bStartRegion)
			{
			BuffIdx += sprintf(&szLineBuff[BuffIdx],"fixedStep chrom=%s start=%d step=1\n",pChrom->szChrom,SeqIdx+1);
			CUtility::SafeWrite(m_hRsltsFile,szLineBuff,BuffIdx);
			BuffIdx = 0;
			bStartRegion = false;
			}
		BuffIdx += sprintf(&szLineBuff[BuffIdx],"%d\n",*pCnts);
		if(BuffIdx + 100 > sizeof(szLineBuff))
			{
			CUtility::SafeWrite(m_hRsltsFile,szLineBuff,BuffIdx);
			BuffIdx = 0;
			}

		}
	}
if(BuffIdx)
	CUtility::SafeWrite(m_hRsltsFile,szLineBuff,BuffIdx);
close(m_hRsltsFile);
m_hRsltsFile = -1;
return(eBSFSuccess);
}




int
AddFeatCounts(int FeatID,		// cnts are for this feature
			  int ChromOfs,		// cnts are associated with this chrom offset
			  int Cnts,
			  int KMerLen,etSeqBase *pSeq)
{
int Idx;
INT64 KMerID;
tsKMerElDist *pElDist; 
size_t ReallocTo;

// determine KMer identifier
if((KMerID = GenSeqIdx(KMerLen,pSeq)) < 0)
	return((int)KMerID);

if(m_pKMerElDist == NULL)
	{
	ReallocTo =  cAllocKMerElDistInsts * sizeof(tsKMerElDist);
#ifdef _WIN32
	pElDist = (tsKMerElDist *) malloc(ReallocTo);	// initial and perhaps the only allocation
	if(pElDist == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"BuildReadCoverage: Memory allocation of %lld bytes - %s",(INT64)ReallocTo,strerror(errno));
		Reset();
		return(eBSFerrMem);
		}
#else
	// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
	pElDist = (tsKMerElDist *)mmap(NULL,ReallocTo, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
	if(pElDist == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"BuildReadCoverage: Memory allocation of %lld bytes through mmap()  failed - %s",(INT64)ReallocTo,strerror(errno));
		pElDist = NULL;
		Reset();
		return(eBSFerrMem);
		}
#endif
	m_pKMerElDist = pElDist;
	m_AllocdKMerElDistInsts = cAllocKMerElDistInsts;
	m_AllocdKMerElDistMem = ReallocTo;
	m_NumKMerElDistInsts = 0;
	Idx = 0;
	}
else
	{
	// locate instance, currently a linear scan but this can be optimised if needed
	pElDist = m_pKMerElDist;
	for(Idx = 0; Idx < m_NumKMerElDistInsts; Idx++,pElDist++)
		if(pElDist->KMerID == KMerID && pElDist->FeatID == FeatID)
			break;
	}

if(!m_NumKMerElDistInsts || Idx == m_NumKMerElDistInsts)			// need to extend and start new instance?
	{
	if(m_NumKMerElDistInsts == m_AllocdKMerElDistInsts)
		{
		ReallocTo = (m_NumKMerElDistInsts + cAllocKMerElDistInsts) * sizeof(tsKMerElDist);
#ifdef _WIN32
		pElDist = (tsKMerElDist *) realloc(m_pKMerElDist,ReallocTo);
#else
		pElDist = (tsKMerElDist *)mremap(m_pKMerElDist,m_AllocdKMerElDistMem,ReallocTo,MREMAP_MAYMOVE);
		if(pElDist == MAP_FAILED)
			pElDist = NULL;
#endif
		if(pElDist == NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddFeatCounts: Memory re-allocation to %d bytes - %s",ReallocTo,strerror(errno));
			Reset();
			return(eBSFerrMem);
			}
		m_pKMerElDist = pElDist;
		}

	pElDist = &m_pKMerElDist[m_NumKMerElDistInsts++];
	memset(pElDist,0,sizeof(tsKMerElDist));
	pElDist->KMerID = KMerID;
	pElDist->FeatID = FeatID;
	pElDist->Min = Cnts;
	pElDist->Max = Cnts;
	}

pElDist->Insts += 1;
pElDist->Cnts += Cnts;
if(pElDist->Min > (UINT32)Cnts)
	pElDist->Min = Cnts;
if(pElDist->Max < (UINT32)Cnts)
	pElDist->Max = Cnts;

if(pElDist->Insts <= cMaxKMerInstCnts)
	{
	pElDist->Dist[pElDist->Insts-1].Cnt = Cnts;
	pElDist->Dist[pElDist->Insts-1].ChromOfs = ChromOfs;
	}
return(pElDist->Insts);
}


int
GenFeatCountVariance(int FeatID,		// cnts are for this feature
					int Cnts,int KMerLen,etSeqBase *pSeq)
{
int Idx;
INT64 KMerID;
tsKMerElDist *pElDist; 

// determine KMer identifier
if((KMerID = GenSeqIdx(KMerLen,pSeq)) < 0)
	return((int)KMerID);

	// locate instance, currently a linear scan but this can be optimised if needed
pElDist = m_pKMerElDist;
for(Idx = 0; Idx < m_NumKMerElDistInsts; Idx++,pElDist++)
	if(pElDist->KMerID == KMerID && pElDist->FeatID == FeatID)
		break;
if(Idx == m_NumKMerElDistInsts)		// should never occur but bugs are bugs...
	return(-1);

if(pElDist->Mean == 0.0)
	pElDist->Mean = (double)pElDist->Cnts / pElDist->Insts;
pElDist->Variance += (pElDist->Mean - Cnts) * (pElDist->Mean - Cnts);
return(pElDist->Insts);
}

// When feature has been processed then generate all population standard deviations for feature K-Mers
int
GenFeatCountStdDevs(int FeatID)		// stddev for this feature
{
int Idx;
int NumFeatKMers;
tsKMerElDist *pElDist; 

if(!m_NumKMerElDistInsts || m_pKMerElDist == NULL)
	return(0);

pElDist = m_pKMerElDist;
NumFeatKMers = 0;
for(Idx = 0; Idx < m_NumKMerElDistInsts; Idx++,pElDist++)
	if(pElDist->FeatID == FeatID)
		{
		pElDist->Variance /= pElDist->Insts;			// using population variance not sample
		pElDist->StdDev = sqrt(pElDist->Variance);
		NumFeatKMers += 1;
		}
return(NumFeatKMers);
}

int
SetDistCounts(int Cnts,int KMerLen,etSeqBase *pSeq)
{
int Idx;
int KMerIdx;
tsKMerDist *pKMer;

// determine index
if((KMerIdx = (int)GenSeqIdx(KMerLen,pSeq)) < 0)
	return(KMerIdx);

pKMer = &m_pKMerDist[KMerIdx];
pKMer->Insts += 1;
pKMer->Cnts += Cnts;

Idx = 0;
while(Idx < cMaxDistLog2 && Cnts > 0) {
	Idx += 1;
	Cnts >>= 1;
	}
pKMer->Dist[Idx-1] += 1;
return(pKMer->Insts);
}


int
GenKMerElDistsZ(int KMerLen,					// process for KMers of this length
			int KMerInstsThres,				// must be at least this number of instances to be reported
			char cStrand,					// only process elements on this strand
			char *pszProcElsFile,			// BED file defining elements of interest start/end loci
			char *pszKMerElRsltsFile)		// results into this file
{
int Rslt;
int CurEntryID;

int ChromSeqLen;
tsChromCnts *pChrom;
UINT16 *pCnts;
int ChromIdx;
int SeqIdx;
etSeqBase *pSeq;
int BuffIdx;
char szLineBuff[0x03fff];


int NumFeats2Proc;
int CurFeatID;
char szFeatName[128];
char szFeatChrom[128];
char szPrevFeatChrom[128];
int FeatStart;
int FeatEnd;
char cFeatStrand;

int NumKMerInsts;
int NumKMerInsts0;
int NumKMerInstsN;
int MinKMerCnts;
int MaxKMerCnts;
double MeanKMerCnts;
double SumDevSqrd;
double Variance;
double PopStdDev;
int KMerIdx;
char *pszSeq;
int KMerInstIdx;
int CntsDist100[cMaxKMerInstCnts];					// to hold K-mer cnts for each feature if read start processing

size_t TotKMerInstCnts;

if((m_pBEDElFile = new CBEDfile)==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CBEDfile");
	Reset();
	return(eBSFerrObj);
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Starting to load BED elements %s",pszProcElsFile);
if((Rslt=m_pBEDElFile->Open(pszProcElsFile,eBTAnyBed)) !=eBSFSuccess)
	{
	while(m_pBEDElFile->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,m_pBEDElFile->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open '%s' for processing",pszProcElsFile);
	Reset();
	return(eBSFerrOpnFile);
	}

NumFeats2Proc = m_pBEDElFile->GetNumFeatures();
if(NumFeats2Proc < 1)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"%s is empty, contains no elements...");
	return(0);
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"%d elements loaded, now processing...",NumFeats2Proc);

#ifdef _WIN32
if((m_hRsltsFile = open(pszKMerElRsltsFile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
if((m_hRsltsFile = open(pszKMerElRsltsFile, O_RDWR | O_CREAT | O_TRUNC, S_IREAD | S_IWRITE ))==-1)
#endif
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create K-mer distribution file: %s - %s",pszKMerElRsltsFile,strerror(errno));
	Reset();
	return(eBSFerrCreateFile);
	}

BuffIdx = 0;
ChromSeqLen = 0;
for(KMerIdx = 0; KMerIdx <  (int)m_KMerDists; KMerIdx++)
	{
	pszSeq = StepIdx2Seq(KMerLen,KMerIdx);
	if(!(KMerIdx % 500))
		printf("Processing K-mer '%s'\r",pszSeq);

	CurFeatID = 0;
	szPrevFeatChrom[0] = '\0';
	while((CurFeatID = m_pBEDElFile->GetNextFeatureID(CurFeatID)) > 0)
		{
		m_pBEDElFile->GetFeature(CurFeatID,szFeatName,szFeatChrom,&FeatStart,&FeatEnd,NULL,&cFeatStrand);
		if(cStrand != '*' && cFeatStrand != cStrand)
			continue;

		// if a feature is on a new chrom then need to load sequence for that chrom
		if(szPrevFeatChrom[0] == '\0' || stricmp(szPrevFeatChrom,szFeatChrom))
			{
			printf("\nProcessing features on chrom %s...",szFeatChrom);
			// get sequence for this feature chrom
			if((CurEntryID = m_pBioSeqFile->LocateEntryIDbyName(szFeatChrom)) < 1)
				{
				// unable to locate any sequences for this 
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to locate sequence for feature chrom %s",szFeatChrom);
				Reset();
				return(CurEntryID);
				}

			// get genome sequence for chrom
			// ensure allocation can hold sequence 
			ChromSeqLen =m_pBioSeqFile->GetDataLen(CurEntryID);
			if(ChromSeqLen <= FeatEnd)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Located sequence for %s but it is too short",szFeatChrom);
				Reset();
				return(eBSFerrFastqSeq);
				}
			
			if(m_pChromSeq == NULL || m_AllocdChromSeq < ChromSeqLen)
				{
				if(m_pChromSeq != NULL)
					{
					delete m_pChromSeq;
					m_pChromSeq = NULL;
					}
				m_AllocdChromSeq = 0;
				if((m_pChromSeq = new UINT8 [ChromSeqLen])==NULL)
					{
					gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate memory '%d' for chrom sequence",ChromSeqLen);
					Reset();
					return(eBSFerrMem);
					}
				m_AllocdChromSeq = ChromSeqLen;
				}

			if((Rslt = m_pBioSeqFile->GetData(CurEntryID,eSeqBaseType,0,m_pChromSeq,ChromSeqLen)) < ChromSeqLen)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"GetData failed for chrom %s",pChrom->szChrom);
				Reset();
				return(Rslt);
				}
	

			pChrom = &m_ChromCnts[0];
			for(ChromIdx = 0 ; ChromIdx < m_NumChromsCov; ChromIdx++,pChrom++)
				if(!stricmp(pChrom->szChrom,szFeatChrom))
					break;
			if(ChromIdx == m_NumChromsCov)
				pChrom = NULL;
			strcpy(szPrevFeatChrom,szFeatChrom);
			}

		// pChrom will be NULL if no reads aligned which are within this current chrom

		// iterate over element sequence, and process counts at KMer instance of interest
		NumKMerInsts = 0;
		NumKMerInsts0 = 0;
		NumKMerInstsN = 0;
		MinKMerCnts = -1;
		MaxKMerCnts = -1;
		TotKMerInstCnts = 0;
		pSeq = &m_pChromSeq[FeatStart];
		for(SeqIdx = FeatStart; SeqIdx < (FeatEnd - KMerLen); SeqIdx++,pSeq++)
			{
			if(KMerIdx == (int)GenSeqIdx(KMerLen,pSeq))
				{
				NumKMerInsts += 1;
				if(pChrom != NULL && (pChrom->StartOfs <= SeqIdx && pChrom->EndOfs >= SeqIdx))
					pCnts = &pChrom->pCovCnts[SeqIdx - pChrom->StartOfs];
				else
					pCnts = NULL;
				if(pCnts == NULL || *pCnts == 0)
					{
					NumKMerInsts0 += 1;
					continue;
					}
				TotKMerInstCnts += *pCnts;
				if(MinKMerCnts == -1 || MinKMerCnts >  *pCnts)
					MinKMerCnts = *pCnts;
				if(MaxKMerCnts == -1 || MaxKMerCnts <  *pCnts)
					MaxKMerCnts = *pCnts;
				NumKMerInstsN += 1;
				}
			}

		
		if(NumKMerInstsN >= KMerInstsThres)
			{
			memset(CntsDist100,0,sizeof(CntsDist100));
			KMerInstIdx = 0;
			MeanKMerCnts = (double)TotKMerInstCnts / NumKMerInstsN;
			// total and mean known, can now calc Variance and population standard deviations
			SumDevSqrd = 0.0;
			pSeq = &m_pChromSeq[FeatStart];
			for(SeqIdx = FeatStart; SeqIdx < (FeatEnd - KMerLen); SeqIdx++,pSeq++)
				{
				if(KMerIdx == (int)GenSeqIdx(KMerLen,pSeq))
					{
					if(pChrom != NULL && (pChrom->StartOfs <= SeqIdx && pChrom->EndOfs >= SeqIdx))
						pCnts = &pChrom->pCovCnts[SeqIdx - pChrom->StartOfs];
					else
						pCnts = NULL;
					if(pCnts == NULL || *pCnts == 0)
						continue;
					if(KMerInstIdx < cMaxKMerInstCnts)
						CntsDist100[KMerInstIdx++] = *pCnts;
					SumDevSqrd += (MeanKMerCnts - (double)*pCnts) * (MeanKMerCnts - (double)*pCnts);
					}
				}

			Variance = SumDevSqrd / NumKMerInstsN;
			PopStdDev = sqrt(Variance);
			}
		else
			continue;
	
		// can now report the stats for this element with the requested KMer
		BuffIdx += sprintf(&szLineBuff[BuffIdx],"%d,\"%s\",\"%s\",%d,%d,%d,%d,%d,%f,%f,%f",
			KMerIdx + 1,pszSeq,szFeatName,NumKMerInsts,NumKMerInstsN,(int)TotKMerInstCnts,MinKMerCnts,MaxKMerCnts,MeanKMerCnts,Variance,PopStdDev);
		for(KMerInstIdx = 0; KMerInstIdx < cMaxKMerInstCnts; KMerInstIdx += 1)
			BuffIdx += sprintf(&szLineBuff[BuffIdx],",%d",CntsDist100[KMerInstIdx]);
		BuffIdx += sprintf(&szLineBuff[BuffIdx],"\n");

		if((BuffIdx + 1000) > sizeof(szLineBuff))
			{
			CUtility::SafeWrite(m_hRsltsFile,szLineBuff,BuffIdx);
			BuffIdx = 0;
			}
		}
	}

if(BuffIdx > 0)
	CUtility::SafeWrite(m_hRsltsFile,szLineBuff,BuffIdx);

close(m_hRsltsFile);
m_hRsltsFile = -1;
return(eBSFSuccess);
}


int
GenKMerElDists(int KMerLen,					// process for KMers of this length
			int KMerInstsThres,				// must be at least this number of instances to be reported
			char cStrand,					// only process elements on this strand
			char *pszProcElsFile,			// BED file defining elements of interest start/end loci
			char *pszKMerElRsltsFile)		// results into this file
{
int Rslt;
int CurEntryID;

int ChromSeqLen;
tsChromCnts *pChrom;
UINT16 *pCnts;
int ChromIdx;
int SeqIdx;
etSeqBase *pSeq;
int BuffIdx;
char szLineBuff[0x03fff];


int NumFeats2Proc;
int CurFeatID;
char szFeatName[128];
char szFeatChrom[128];
char szPrevFeatChrom[128];
int FeatStart;
int FeatEnd;
char cFeatStrand;

int NumKMerInsts;
int NumKMerInsts0;
int NumKMerInstsN;
int MinKMerCnts;
int MaxKMerCnts;

int KMerIdx;
char *pszSeq;
int KMerInstIdx;

size_t TotKMerInstCnts;

if((m_pBEDElFile = new CBEDfile)==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CBEDfile");
	Reset();
	return(eBSFerrObj);
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Starting to load BED elements %s",pszProcElsFile);
if((Rslt=m_pBEDElFile->Open(pszProcElsFile,eBTAnyBed)) !=eBSFSuccess)
	{
	while(m_pBEDElFile->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,m_pBEDElFile->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open '%s' for processing",pszProcElsFile);
	Reset();
	return(eBSFerrOpnFile);
	}

NumFeats2Proc = m_pBEDElFile->GetNumFeatures();
if(NumFeats2Proc < 1)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"%s is empty, contains no elements...");
	return(0);
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"%d elements loaded, now processing...",NumFeats2Proc);

#ifdef _WIN32
if((m_hRsltsFile = open(pszKMerElRsltsFile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
if((m_hRsltsFile = open(pszKMerElRsltsFile, O_RDWR | O_CREAT | O_TRUNC, S_IREAD | S_IWRITE ))==-1)
#endif
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create K-mer distribution file: %s - %s",pszKMerElRsltsFile,strerror(errno));
	Reset();
	return(eBSFerrCreateFile);
	}

BuffIdx = 0;
ChromSeqLen = 0;

CurFeatID = 0;
szPrevFeatChrom[0] = '\0';
printf("\nfeat %d\r",CurFeatID);
while((CurFeatID = m_pBEDElFile->GetNextFeatureID(CurFeatID)) > 0)
	{
	if(!(CurFeatID % 1000))
		printf("feat %d\r",CurFeatID);

	m_pBEDElFile->GetFeature(CurFeatID,szFeatName,szFeatChrom,&FeatStart,&FeatEnd,NULL,&cFeatStrand);
	if(cStrand != '*' && cFeatStrand != cStrand)
		continue;

	// if a feature is on a new chrom then need to load sequence for that chrom
	if(szPrevFeatChrom[0] == '\0' || stricmp(szPrevFeatChrom,szFeatChrom))
		{
		printf("\nProcessing features on chrom %s...",szFeatChrom);
		// get sequence for this feature chrom
		if((CurEntryID = m_pBioSeqFile->LocateEntryIDbyName(szFeatChrom)) < 1)
			{
			// unable to locate any sequences for this 
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to locate sequence for feature chrom %s",szFeatChrom);
			Reset();
			return(CurEntryID);
			}

		// get genome sequence for chrom
		// ensure allocation can hold sequence 
		ChromSeqLen =m_pBioSeqFile->GetDataLen(CurEntryID);
		if(ChromSeqLen <= FeatEnd)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Located sequence for %s but it is too short",szFeatChrom);
			Reset();
			return(eBSFerrFastqSeq);
			}
			
		if(m_pChromSeq == NULL || m_AllocdChromSeq < ChromSeqLen)
			{
			if(m_pChromSeq != NULL)
				{
				delete m_pChromSeq;
				m_pChromSeq = NULL;
				}
			m_AllocdChromSeq = 0;
			if((m_pChromSeq = new UINT8 [ChromSeqLen])==NULL)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate memory '%d' for chrom sequence",ChromSeqLen);
				Reset();
				return(eBSFerrMem);
				}
			m_AllocdChromSeq = ChromSeqLen;
			}

		if((Rslt = m_pBioSeqFile->GetData(CurEntryID,eSeqBaseType,0,m_pChromSeq,ChromSeqLen)) < ChromSeqLen)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"GetData failed for chrom %s",pChrom->szChrom);
			Reset();
			return(Rslt);
			}
	

		pChrom = &m_ChromCnts[0];
		for(ChromIdx = 0 ; ChromIdx < m_NumChromsCov; ChromIdx++,pChrom++)
			if(!stricmp(pChrom->szChrom,szFeatChrom))
				break;
		if(ChromIdx == m_NumChromsCov)
			pChrom = NULL;
		strcpy(szPrevFeatChrom,szFeatChrom);
		}

	// pChrom will be NULL if no reads aligned which are within this current chrom
	// iterate over element sequence, and process counts at KMer instance of interest
	NumKMerInsts = 0;
	NumKMerInsts0 = 0;
	NumKMerInstsN = 0;
	MinKMerCnts = -1;
	MaxKMerCnts = -1;
	TotKMerInstCnts = 0;

	m_NumKMerElDistInsts = 0;
	pSeq = &m_pChromSeq[FeatStart];
	for(SeqIdx = FeatStart; SeqIdx < (FeatEnd - KMerLen); SeqIdx++,pSeq++)
		{
		if(pChrom != NULL && (pChrom->StartOfs <= SeqIdx && pChrom->EndOfs >= SeqIdx))
			pCnts = &pChrom->pCovCnts[SeqIdx - pChrom->StartOfs];
		else
			pCnts = NULL;
		if(pCnts == NULL || *pCnts == 0)
			{
			NumKMerInsts0 += 1;
			continue;
			}

		AddFeatCounts(CurFeatID,SeqIdx,*pCnts,KMerLen,pSeq);
		}

	if(m_NumKMerElDistInsts)
		{
		pSeq = &m_pChromSeq[FeatStart];
		for(SeqIdx = FeatStart; SeqIdx < (FeatEnd - KMerLen); SeqIdx++,pSeq++)
			{
			if(pChrom != NULL && (pChrom->StartOfs <= SeqIdx && pChrom->EndOfs >= SeqIdx))
				pCnts = &pChrom->pCovCnts[SeqIdx - pChrom->StartOfs];
			else
				pCnts = NULL;
			if(pCnts == NULL || *pCnts == 0)
				continue;

			GenFeatCountVariance(CurFeatID,*pCnts,KMerLen,pSeq);
			}
		GenFeatCountStdDevs(CurFeatID);

		tsKMerElDist *pElDist = m_pKMerElDist; 
		for(KMerIdx = 0; KMerIdx < m_NumKMerElDistInsts; KMerIdx++,pElDist++)
			{
			if(pElDist->FeatID != CurFeatID || pElDist->Insts < (UINT32)KMerInstsThres)
				continue;

			// can now report the stats for this element with the requested KMer
			pszSeq = StepIdx2Seq(KMerLen,pElDist->KMerID);
			BuffIdx += sprintf(&szLineBuff[BuffIdx],"%d,\"%s\",\"%s\",\"%s\",%d,%d,%d,%d,%d,%f,%f,%f",
				KMerIdx + 1,pszSeq,szFeatName,szFeatChrom,FeatStart,FeatEnd,pElDist->Insts,pElDist->Min,pElDist->Max,pElDist->Mean,pElDist->Variance,pElDist->StdDev);
				
			for(KMerInstIdx = 0; KMerInstIdx < cMaxKMerInstCnts; KMerInstIdx += 1)
				{
				BuffIdx += sprintf(&szLineBuff[BuffIdx],",%d,%d",pElDist->Dist[KMerInstIdx].Cnt,pElDist->Dist[KMerInstIdx].ChromOfs);
				if(pElDist->Dist[KMerInstIdx].Cnt)
					{
					if(pElDist->Dist[KMerInstIdx].ChromOfs < 10)
						pSeq = &m_pChromSeq[0];
					else
						pSeq = &m_pChromSeq[pElDist->Dist[KMerInstIdx].ChromOfs-10];
					pszSeq = CSeqTrans::MapSeq2Ascii(pSeq,10);
					BuffIdx += sprintf(&szLineBuff[BuffIdx],",\"%s\"",pszSeq);
					if((pElDist->Dist[KMerInstIdx].ChromOfs + KMerLen) > ChromSeqLen)
						pSeq = &m_pChromSeq[ChromSeqLen - 10];
					else
						pSeq = &m_pChromSeq[pElDist->Dist[KMerInstIdx].ChromOfs+KMerLen];
					pszSeq = CSeqTrans::MapSeq2Ascii(pSeq,10);
					BuffIdx += sprintf(&szLineBuff[BuffIdx],",\"%s\"",pszSeq);
					}
				else
					BuffIdx += sprintf(&szLineBuff[BuffIdx],",\"\",\"\"");
				}
			BuffIdx += sprintf(&szLineBuff[BuffIdx],"\n");

			if((BuffIdx + 1000) > sizeof(szLineBuff))
				{
				CUtility::SafeWrite(m_hRsltsFile,szLineBuff,BuffIdx);
				BuffIdx = 0;
				}
			}
		}
	}


if(BuffIdx > 0)
	CUtility::SafeWrite(m_hRsltsFile,szLineBuff,BuffIdx);

close(m_hRsltsFile);
m_hRsltsFile = -1;
return(eBSFSuccess);
}



int
GenKMerDists(int KMerLen,
			int KMerInstsThres,				// only report k-mers with cnts of at least this threshold
			char *pszKMerRsltsFile)			// output K-mer distributions file
{
int Rslt;
int Idx;
int CurEntryID;
int ChromSeqLen;
tsChromCnts *pChrom;
UINT16 *pCnts;
int ChromIdx;
int SeqIdx;
etSeqBase *pSeq;
int BuffIdx;
char szLineBuff[8096];

#ifdef _WIN32
if((m_hRsltsFile = open(pszKMerRsltsFile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
if((m_hRsltsFile = open(pszKMerRsltsFile, O_RDWR | O_CREAT | O_TRUNC, S_IREAD | S_IWRITE ))==-1)
#endif
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create K-mer distribution file: %s - %s",pszKMerRsltsFile,strerror(errno));
	Reset();
	return(eBSFerrCreateFile);
	}

pChrom = &m_ChromCnts[0];
for(ChromIdx = 0 ; ChromIdx < m_NumChromsCov; ChromIdx++,pChrom++)
	{
	if((CurEntryID = m_pBioSeqFile->LocateEntryIDbyName(pChrom->szChrom)) < 1)
		{
		// unable to locate any sequences for this 
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to locate sequence for %s",pChrom->szChrom);
		Reset();
		return(CurEntryID);
		}
	

	// ensure that chrom sequence is at least that length as covered by reads
	ChromSeqLen =m_pBioSeqFile->GetDataLen(CurEntryID);
	if(ChromSeqLen <= pChrom->EndOfs)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Located sequence for %s but it is too short",pChrom->szChrom);
		Reset();
		return(eBSFerrFastqSeq);
		}

	// get genome sequence for current chrom
	// ensure allocation can hold sequence
	if(m_pChromSeq == NULL || m_AllocdChromSeq < ChromSeqLen)
		{
		if(m_pChromSeq != NULL)
			{
			delete m_pChromSeq;
			m_pChromSeq = NULL;
			}
		m_AllocdChromSeq = 0;
		if((m_pChromSeq = new UINT8 [ChromSeqLen])==NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate memory '%d' for chrom sequence",ChromSeqLen);
			Reset();
			return(eBSFerrMem);
			}
		m_AllocdChromSeq = ChromSeqLen;
		}

	if((Rslt = m_pBioSeqFile->GetData(CurEntryID,eSeqBaseType,0,m_pChromSeq,ChromSeqLen)) < ChromSeqLen)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"GetData failed for chrom %s",pChrom->szChrom);
		Reset();
		return(Rslt);
		}
	
	pCnts = &pChrom->pCovCnts[pChrom->StartOfs];
	pSeq = &m_pChromSeq[pChrom->StartOfs];
	for(SeqIdx = pChrom->StartOfs; SeqIdx < (pChrom->EndOfs - KMerLen); SeqIdx++,pCnts++,pSeq++)
		{
		if(*pCnts == 0)
			continue;

		SetDistCounts(*pCnts,KMerLen,pSeq);
		}
	}

// now write out the distributions
tsKMerDist *pKMer;
int KMerIdx;
char *pszSeq;

pKMer = m_pKMerDist;
BuffIdx = 0;
for(KMerIdx = 0; KMerIdx < (int)m_KMerDists; KMerIdx++,pKMer++)
	{
	if(pKMer->Insts < (UINT32)KMerInstsThres)
		continue;
	pszSeq = StepIdx2Seq(KMerLen,KMerIdx);
	BuffIdx += sprintf(&szLineBuff[BuffIdx],"%d,\"%s\",%d,%d,%f",KMerIdx + 1,pszSeq,pKMer->Insts, pKMer->Cnts,
									pKMer->Insts > 0 ? (double)pKMer->Cnts / pKMer->Insts : 0.0f );
	for(Idx = 0; Idx < cMaxDistLog2; Idx++)
		BuffIdx += sprintf(&szLineBuff[BuffIdx],",%d",pKMer->Dist[Idx]);
	BuffIdx += sprintf(&szLineBuff[BuffIdx],"\n");
	if((BuffIdx + 1000) > sizeof(szLineBuff))
		{
		CUtility::SafeWrite(m_hRsltsFile,szLineBuff,BuffIdx);
		BuffIdx = 0;
		}
	}
if(BuffIdx > 0)
	CUtility::SafeWrite(m_hRsltsFile,szLineBuff,BuffIdx);
close(m_hRsltsFile);
m_hRsltsFile = -1;
return(eBSFSuccess);
}

int GenCSVWiggle(bool bStartOnly,				// if true then only process for read starts, not complete read overlaps
				 char FiltStrand,				// process for alignments on this strand only
				 char *pszInFile)				// CSV loci input file
{
int Rslt;
int NumEls;
int NumFields;
char *pszStrand;
char *pszChrom;
int StartLoci;
int EndLoci;

if((m_pCSVFile = new CCSVFile)==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CCSVFile");
	Reset();
	return(eBSFerrObj);
	}

if((Rslt=m_pCSVFile->Open(pszInFile)) !=eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open '%s' for processing",pszInFile);
	Reset();
	return(eBSFerrOpnFile);
	}
NumEls = 0;
while((Rslt=m_pCSVFile->NextLine()) > 0)	// onto next line containing fields
	{
	NumFields = m_pCSVFile->GetCurFields();
	if(NumFields < 7)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expected at least 7 fields in '%s', GetCurFields() returned '%d'",pszInFile,NumFields);
		return(eBSFerrFieldCnt);
		}

	if(FiltStrand != '*')
		{
		if(NumFields >= 8)					// check if strand has been specified
			{
			m_pCSVFile->GetText(8,&pszStrand);
			if(pszStrand[0] != '-')			// assume anything other than '-' is on the plus strand
				pszStrand[0] = '+';
			if(FiltStrand != pszStrand[0])
				continue;
			}
		else
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Strand filtering requested, expected at least 8 fields in '%s', GetCurFields() returned '%d'",pszInFile,NumFields);
			return(eBSFerrFieldCnt);
			}
		}
	NumEls += 1;
	m_pCSVFile->GetText(4,&pszChrom);
	m_pCSVFile->GetInt(5,&StartLoci);
	m_pCSVFile->GetInt(6,&EndLoci);

	if((Rslt=BuildReadCoverage(bStartOnly,pszChrom,StartLoci,EndLoci,1))!=eBSFSuccess)
		break;
	}

return(Rslt);
}

int GenBEDWiggle(bool bStartOnly,				// if true then only process for read starts, not complete read overlaps
				 char FiltStrand,				// process for alignments on this strand only
				 char *pszInFile)				// UCSC BED input file
{
int Rslt;
int CurFeatureID;
int StartLoci;
int EndLoci;
int Score;
char szChrom[128];
char szPrevChrom[128];
char szFeatName[128];
char Strand;
int RefID;
int IntergenicStart;
int NumEls;

if((m_pBEDFile = new CBEDfile)==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CBEDfile");
	Reset();
	return(eBSFerrObj);
	}

if((Rslt=m_pBEDFile->Open(pszInFile,eBTAnyBed)) !=eBSFSuccess)
	{
	while(m_pBEDFile->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,m_pBEDFile->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open '%s' for processing",pszInFile);
	Reset();
	return(eBSFerrOpnFile);
	}

// now iterate over the features, filtering as may be appropriate
szPrevChrom[0] = '\0';
CurFeatureID = 0;
IntergenicStart = 0;
RefID = 0;
Rslt = eBSFSuccess;
NumEls = 0;
while(Rslt == eBSFSuccess && (CurFeatureID = m_pBEDFile->GetNextFeatureID(CurFeatureID)) > 0)
	{
	m_pBEDFile->GetFeature(CurFeatureID,// feature instance identifier
				szFeatName,				// where to return feature name
				szChrom,				// where to return chromosome name
				&StartLoci,				// where to return feature start on chromosome (0..n) 
				&EndLoci,				// where to return feature end on chromosome
 				&Score,					// where to return score
 				&Strand);				// where to return strand

	if(CurFeatureID == 1 || stricmp(szChrom,szPrevChrom))	// if new chromosome then reset IntergenicStart
		{
		strcpy(szPrevChrom,szChrom);
		IntergenicStart = 0;
		}

	if(FiltStrand != '*' && FiltStrand != Strand)
		continue;

	NumEls += 1;
	Rslt=BuildReadCoverage(bStartOnly,szChrom,StartLoci,EndLoci,1);
	}

return(Rslt);
}


int
Process(etPMode PMode,					// processing mode
		etIFormat IMode,					// input aligned reads file format 
		char Strand,					// process for alignments on this strand only
		int KMerLen,					// process for this length K-Mers
		int KMerInstsThres,				// report KMers with instances >= this threshold
		char *pszInFile,				// input BED file or CSV file containing aligned read loci
		char *pszProcElsFile,			// optional input BED file containing elements over which K-mer deviations to be reported
		char *pszInBioseqFile,			// bioseq genome file
		char *pszRsltsFile,				// optional output UCSC wiggle file
	    char *pszKMerRsltsFile)			// output K-mer distributions file
{
int Rslt;
int Pwr;
bool bStartOnly;
char szKMerElRsltsFile[_MAX_PATH];

Init();

bStartOnly = PMode == ePMdefault ? true : false;

if((m_pBioSeqFile = new CBioSeqFile) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CBioSeqFile object");
	Reset();
	return(eBSFerrObj);
	}

if((Rslt = m_pBioSeqFile->Open(pszInBioseqFile))!=eBSFSuccess)
	{
	while(m_pBioSeqFile->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,m_pBioSeqFile->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open genome assembly sequence file '%s'",pszInBioseqFile);
	Reset();
	return(Rslt);
	}

m_KMerDists = 4;
for(Pwr = 1; Pwr < KMerLen; Pwr++)
	m_KMerDists *= 4;
if((m_pKMerDist = new tsKMerDist [m_KMerDists])==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate memory '%d' for KMerDists",m_KMerDists * sizeof(tsKMerDist));
	Reset();
	return(eBSFerrMem);
	}
memset(m_pKMerDist,0,m_KMerDists * sizeof(tsKMerDist));

switch(IMode) {
	case eIFCSV:		// default is for CSV loci processing
		Rslt = GenCSVWiggle(bStartOnly,Strand,pszInFile);
		break;

	case eIFdefault:			// UCSC BED processing
		Rslt = GenBEDWiggle(bStartOnly,Strand,pszInFile);
		break;

	default:
		break;
	}
if(Rslt == eBSFSuccess)
	{
	if((Rslt = GenKMerDists(KMerLen,KMerInstsThres,pszKMerRsltsFile)) >= 0)
		{
		// genome wide alignment KMer distributions known, now if element file specified then look at the variation intra-element
		if(pszProcElsFile != NULL && pszProcElsFile[0] != '\0')
			{
			sprintf(szKMerElRsltsFile,"%s.e",pszKMerRsltsFile);
			Rslt = GenKMerElDists(KMerLen,KMerInstsThres,Strand,pszProcElsFile,szKMerElRsltsFile);
			}
		if(Rslt >= 0 && pszRsltsFile != NULL && pszRsltsFile[0] != '\0')
			Rslt = WriteReadsWig(pszInFile,pszRsltsFile);
		}
	}
Reset();
return(Rslt);
}



