// fasta2dist.cpp : Defines the entry point for the console application.
// Processes probes against targets
// Probes and targets can be either either multifasta *.fa) or bioseq (.seq)
// 

#include "stdafx.h"
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#if _WIN32

#include "../libbiokanga/commhdrs.h"

#else
#include "../libbiokanga/commhdrs.h"
#endif

const unsigned int cProgVer = 101;			// increment with each release

CStopWatch gStopWatch;
CDiagnostics gDiagnostics;					// for writing diagnostics messages to log file
char gszProcName[_MAX_FNAME];				// process name

const int cAllocMatchLoci = 100000;			// allocate match loci blocks holding this many match loci

const int cMinProbeLen = 2;				    // minimum length probes
const int cDfltMinProbeLen = 4;				// default minimum length probe
const int cMaxProbeLen = 10000;				// maximum length probes

const int cMinBinWidth  = 10000;		    // minimum width bins
const int cDfltBinWidth = 1000000;			// default width bins
const int cMaxBinWidth  = 1000000000;		// maximum width bins

const int cDfltMismatches = 0;				// default number of allowed mismatches
const int cMaxMismatches = 100;				// max number of allowed mismatches

typedef enum TAG_eProcMode {				// processing modes
	ePMHChromStats = 0,						// default is to report COV, stdev and mean for each chromosome
	ePMHChromBins							// otherwise report for bins along length of each chromosome	
	} teProcMode;

typedef struct TAG_sFSAMProbe {
	struct TAG_sFSAMProbe *pNext;			// probes are linked
	int ProbID;								// uniquely identifies this probe
	char szDescr[80];						// as parsed from fasta descriptor
	int TotHits;							// total hits over all chroms
	int CurChromHits;						// total hits on current chromosome

	double CurChromMean;					// mean number of hits over bins in current chromosome
	double CurChromStdDev;					// stddev of hits over bins in current chromosome
	double CurChromCOV;						// COV of hits over bins in current chromosome
	int ProbeLen;							// number of bases in Bases[]
	etSeqBase Bases[1];						// to hold probe sequence
} sFSAMProbe;

typedef struct TAG_sProcCtx {
	teProcMode ProcMode;	// processing mode
	char *pszProbeFile;		// probes from this multfasta file
	char *pszTargFile;		// targets from this bioseq file
	char *pszRsltsFile;		// output to this results file
	int MinLength;			// process probes which are of at least this length
	int MaxLength;			// process probes which are no longer than this length
	char *pszChroms;		// only process these target chroms
	char Strand;			// target strand to be processed - '+', '-' or '*' for both
	int hRsltsFile;			// handle for opened results file
	int BinWidth;			// chromosome histogram bin width
	int MaxMismatches;		// maximum number of mismatches accepted

	int CurTargID;			// target identifier
	char *pszCurChrom;		// current chrom being processed
	int NumBins;			// number of bins for current chromosome
	int AllocBins;			// number of bins currently allocated (each bin contains slots for all probes)
	unsigned int *pBins;	// ptr to allocated bins

	int NumProbes;			// number of probes loaded
	int MaxProbeLen;		// max length of any probe
	int MinProbeLen;		// min length of any probe
	sFSAMProbe *pProbes;	// pts to linked list of sAMProbes
} tsProcCtx;




int Process(teProcMode Mode,
			char *pszProbeFile,			// probes from this multfasta file
			char *pszTargFile,			// targets from this bioseq file
			char *pszRsltsFile,			// hits into this results file
			int MinLength,				// process probes which are of at least this length
			int MaxLength,				// process probes which are no longer than this length
			char *pszChroms,			// only process these target chroms
			char Strand,				// process against this target strand - '+','-' or '*' for both
			int MaxMismatches,			// accept up to at most this many mismatches in any hit
			int BinWidth);				// use this sized histogram bins when accumulating chromosome probe hits

int ReportMatches(CBioSeqFile *pBioSeq,int hRsltsFile);
int ReportCounts(int hRsltsFile,char *pszChrom);
int ReportMisses(CBioSeqFile *pBioSeq,int hRsltsFile,char *pszChrom);
int AddMatch(sFSAMProbe *pProbe,char Strand,int Loci,tsProcCtx *pProcCtx);
int LoadProbes(tsProcCtx *pProcCtx);
void DeleteProbes(tsProcCtx *pProcCtx);
int										// returns the number of hits by probes against the chromosome sequence
FindMatches(tsProcCtx *pProcCtx,			// processing parameters
		   tBSFEntryID TargEntryID,		// current target entry identifier
		   unsigned char *pTargSeq,		// pts to target sequence
		   int TargSeqLen,				// target sequence length
		   char *pszChrom);				// target is on this chromosome 
int	CleanNameList(char *pszTxt);
bool NameInList(char *pszNameList, char *pszName);

int CalcChromDensity(tsProcCtx *pProcCtx);
int ReportChromDensities(tsProcCtx *pProcCtx);
int ReportAllChromDensities(tsProcCtx *pProcCtx);

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
int iMode;					// Processing mode: 0 Target hit loci, 1: Probe hit counts
int iMaxMismatches;			//accept up to at most this many mismatches in any hit (0..100), default is 0");
char cStrand;				//Strand '+' or '-', default is '*' for both");
int iMinLength;				//Filter out probes of less than this length (default: 4, min: 2)");
int iMaxLength;				//Filter out probes of longer than this length (default: 10000, max: 10000)");
int iBinWidth;				// chromosome histogram bin width
char szChroms[2048];		//Comma/space separated list of chromosomes to process (default: all)");

char szRsltsFile[_MAX_PATH];
char szTargFile[_MAX_PATH];
char szProbeFile[_MAX_PATH];

// command line args
struct arg_lit  *help    = arg_lit0("hH","help",                "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_int *ScreenLogLevel=arg_int0("S", "ScreenLogLevel",	"<int>","Level of diagnostics written to screen 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"diagnostics log file");

struct arg_int *Mode=arg_int0("m", "mode","<int>",				"Processing mode: 0 report COV, stdev, mean for each chromosome, 1: hit density distribution along each chromosome");
struct arg_file *ProbeFile = arg_file1("i","probes","<file>",	"input from multifasta (.fa) probe file");
struct arg_file *TargFile = arg_file1("I","targets","<file>",	"input from bioseq (.seq) target file");
struct arg_file *RsltsFile = arg_file1("o","allresults","<file>",	"output results to this file as CSV");
struct arg_int *BinWidth=arg_int0("b", "binwidth","<int>",		"chromosome histogram bin width (default: 1000000 nt)");
struct arg_int *MaxMismatches=arg_int0("t", "maxmissmatches","<int>",	"accept up to at most this many mismatches in any hit (0..100), default is 0");
struct arg_str *Strand=arg_str0("s", "strand",	"<string>",		"Strand '+' or '-', default is '*' for both");
struct arg_int *MinLength=arg_int0("l", "Minlen","<int>",		"Filter out probes of less than this length (default: 4, min: 2)");
struct arg_int *MaxLength=arg_int0("L", "Maxlen","<int>",		"Filter out probes of longer than this length (default: 10000, max: 10000)");
struct arg_str *Chroms=arg_strn("c", "chroms",	"<string>",	0,50,	"Comma/space separated list of chromosomes to process (default: all)");


struct arg_end *end = arg_end(20);

void *argtable[] = {help,version,FileLogLevel,ScreenLogLevel,LogFile,
					Mode,ProbeFile,TargFile,RsltsFile,MaxMismatches,BinWidth,
					Strand,MinLength,MaxLength,Chroms,
					end};
char **pAllArgs;
int argerrors;
argerrors = CUtility::arg_parsefromfile(argc,(char **)argv,&pAllArgs);
if(argerrors >= 0)
	argerrors = arg_parse(argerrors,pAllArgs,argtable);

/* special case: '--help' takes precedence over error reporting */
if (help->count > 0)
        {
		printf("\n%s ", gszProcName);
        arg_print_syntax(stdout,argtable,"\n");
        arg_print_glossary(stdout,argtable,"  %-25s %s\n");
		printf("\nNote: Parameters can be entered into a parameter file, one parameter per line.");
		printf("\n      To invoke this parameter file then precede its name with '@'");
		printf("\n      e.g. %s @myparams.txt\n\n",gszProcName);
		exit(1);
        }

    /* special case: '--version' takes precedence error reporting */
if (version->count > 0)
        {
		printf("\n%s Version %d.%2.2d\n",gszProcName,cProgVer/100,cProgVer%100);
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

	iMode = Mode->count ? Mode->ival[0] : ePMHChromStats;
	if(iMode < ePMHChromStats || iMode > ePMHChromBins)
		{
		printf("\nError: Processing mode '-m%d' not in range %d..%d\n",iMode,ePMHChromStats,ePMHChromBins);
		exit(1);
		}

	iMinLength = MinLength->count ? MinLength->ival[0] : cDfltMinProbeLen;
	if(iMinLength < cMinProbeLen || iMinLength > cMaxProbeLen)
		{
		printf("\nError: Minimum length probe '-l%d' not in range %d..%d\n",iMinLength,cMinProbeLen,cMaxProbeLen);
		exit(1);
		}
	iMaxLength = MaxLength->count ? MaxLength->ival[0] : cMaxProbeLen;
	if(iMaxLength < iMinLength || iMaxLength > cMaxProbeLen)
		{
		printf("\nError: Maximum length probe '-L%d' not in range %d..%d\n",iMaxLength,iMinLength,cMaxProbeLen);
		exit(1);
		}

	iMaxMismatches = MaxMismatches->count ? MaxMismatches->ival[0] : cDfltMismatches;
	if(iMaxMismatches < 0 || iMaxMismatches > cMaxMismatches)
		{
		printf("\nError: Maximum number of mismatches '-t%d' not in range 0..%d\n",iMaxMismatches,cDfltMismatches);
		exit(1);
		}

	iBinWidth = BinWidth->count ? BinWidth->ival[0] : cDfltBinWidth;
	if(iBinWidth < cMinBinWidth || iBinWidth > cMaxBinWidth)
		{
		printf("\nError: bin width '-b%d' not in range %d..%d\n",iBinWidth,cMinBinWidth,cMaxBinWidth);
		exit(1);
		}

	szChroms[0] = '\0';
	if(Chroms->count)
		{
		for(int Idx = 0; Idx < Chroms->count; Idx++)
			{
			if(Idx)
				strcat(szChroms,",");
			strcat(szChroms,Chroms->sval[Idx]);
			}
		CleanNameList(szChroms);
		}
	
	if(Strand->count)
		{
		cStrand = Strand->sval[0][0];
		if(cStrand != '*' && cStrand != '+' && cStrand != '-' && cStrand != '1' && cStrand != '0')
			{
			printf("\nError: Requested strand '-s%c' is invalid, must be '-s+' or '-s-' or '-s*'",cStrand);
			exit(1);
			}
		if(cStrand == '0')
			cStrand = '+';
		else
			if(cStrand == '0')
				cStrand = '-';
		}
	else
		cStrand = '*';

	strncpy(szProbeFile,ProbeFile->filename[0],sizeof(szProbeFile));
	szProbeFile[sizeof(szProbeFile)-1] = '\0';
	
	strncpy(szTargFile,TargFile->filename[0],sizeof(szTargFile));
	szTargFile[sizeof(szTargFile)-1] = '\0';

	strncpy(szRsltsFile,RsltsFile->filename[0],sizeof(szRsltsFile));
	szRsltsFile[sizeof(szRsltsFile)-1] = '\0';

		// now that command parameters have been parsed then initialise diagnostics log system
	if(!gDiagnostics.Open(szLogFile,(etDiagLevel)iScreenLogLevel,(etDiagLevel)iFileLogLevel,true))
		{
		printf("\nError: Unable to start diagnostics subsystem.");
		if(szLogFile[0] != '\0')
			printf(" Most likely cause is that logfile '%s' can't be opened/created",szLogFile);
		exit(1);
		}

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Version: %d.%2.2d Processing parameters:",cProgVer/100,cProgVer%100);
	switch(iMode) {
		case ePMHChromStats:					// default is to report densities over all chromosomes
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Processing mode: %s","Report COV,mean and stdev for each chromosome");
			break;
		case ePMHChromBins:					// otherwise report for each chromosome
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Processing mode: %s","Report hit density along each chromosome");
			break;
		}
	
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Probes from multifasta file: %s",szProbeFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Targets from bioseq sequence file: %s",szTargFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Results to file: '%s'", szRsltsFile);

	if(szChroms[0] != '\0')
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Only process these target chromosomes: %s",szChroms);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Filter out probes with length less than: %d",iMinLength);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Filter out probes with length longer than: %d",iMaxLength);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Target strand: '%c'",cStrand);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Max mismatches: %d",iMaxMismatches);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Use chromosome histogram bin widths of: %d",iBinWidth);

#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif

	gStopWatch.Start();
	Rslt = Process((teProcMode)iMode,
				szProbeFile,			// probes from this multfasta file
				szTargFile,				// targets from this bioseq file
				szRsltsFile,			// output to this results file
				iMinLength,				// process probes which are of at least this length
				iMaxLength,				// process probes which are no longer than this length
				szChroms,				// only process these target chroms
				cStrand,				// process against this target strand - '+','-' or '*' for both
				iMaxMismatches,			// accept up to at most this many mismatches in any hit
				iBinWidth);				// use this sized histogram bins when accumulating chromosome probe hits


	gStopWatch.Stop();
	Rslt = Rslt >=0 ? 0 : 1;
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Exit code: %d Total processing time: %s",Rslt,gStopWatch.Read());
	exit(Rslt);
	}
else
	{
	arg_print_errors(stdout,end,gszProcName);
	arg_print_syntax(stdout,argtable,"\n [@<parameterfile>]\nend of help\n");
	exit(1);
	}
}

// CleanNameList
// a) Changes all occurances of one or more whitespace to single comma's
// b) Removes all whitespace, characters < 0x020, characters > 0x07f, single and double quotes
int				// returns number of names
CleanNameList(char *pszTxt)
{
char *pDst;
char Chr;
bool bInName;
bool bComma;
int NumNames;

if(pszTxt == NULL || pszTxt[0] == '\0')
	return(0);
pDst = pszTxt;
bInName = false;
bComma = false;
NumNames = 0;
while((Chr = *pszTxt++))	
	{
	switch(Chr) {
		case ',':
			if(!bComma && NumNames)
				{
				*pDst++ = ',';
				bComma = true;
				bInName = false;
				}
			continue;

		case '"':				// simply slough all double quotes
		case '\'':				// simply slough all single quotes
			continue;

		default:
			if(isspace(Chr))
				{
				if(bComma || !NumNames)
				   continue;
				*pDst++ = ',';
				bComma = true;
				bInName = false;
				continue;
				}
			if(Chr < 0x021 || Chr > 0x07e)
				continue;
			*pDst++ = Chr;
			if(!bInName)
				{
				bInName = true;
				bComma = false;
				NumNames += 1;
				}
			continue;
		}
	}
*pDst = '\0';
if(NumNames && pDst[-1] == ',') // ensure list not terminated by comma
	pDst[-1] = '\0';
return(NumNames);
}

bool			// true if name in CSV name list, otherwise returns false
NameInList(char *pszNameList, char *pszName)
{
int NameLen;

if(pszNameList == NULL || pszNameList[0] == '\0' ||
   pszName == NULL || pszName[0] == '\0')
   return(false);

NameLen = (int)strlen(pszName);
while(*pszNameList != '\0')
	{
	if(!strnicmp(pszNameList,pszName,NameLen))
		return(true);
	while(*pszNameList != '\0' && *pszNameList++ != ',');
	}
return(false);
}


int Process(teProcMode ProcMode,
			char *pszProbeFile,			// probes from this multfasta file
			char *pszTargFile,			// targets from this bioseq file
			char *pszRsltsFile,			// hits into this results file
			int MinLength,				// process probes which are of at least this length
			int MaxLength,				// process probes which are no longer than this length
			char *pszChroms,			// only process these target chroms
			char Strand,				// process against this target strand - '+','-' or '*' for both
			int MaxMismatches,			// accept up to at most this many mismatches in any hit
			int BinWidth)				// use this sized histogram bins when accumulating chromosome probe hits
{
CBioSeqFile *pBioSeq;
unsigned char *pTargSeq;			// to hold sequence for each chromosome
int AllocdTargSeqLen;				// length of currently allocated pTargSeq
int hRsltsFile;
char szChromName[cBSFSourceSize+1];	// to hold current chromosome name
tBSFEntryID CurEntry;				// current entry being processed
int TargSeqLen;						// length of chromosome sequence being processed
int Rslt;							// used to hold processing result

tsProcCtx ProcCtx;
memset(&ProcCtx,0,sizeof(ProcCtx));
ProcCtx.ProcMode = ProcMode;
ProcCtx.pszProbeFile = pszProbeFile;
ProcCtx.pszTargFile = pszTargFile;
ProcCtx.pszRsltsFile = pszRsltsFile;
ProcCtx.MinLength = MinLength;
ProcCtx.MaxLength = MaxLength;
ProcCtx.pszChroms = pszChroms;
ProcCtx.hRsltsFile = -1;
ProcCtx.Strand = Strand;
ProcCtx.BinWidth = BinWidth;
ProcCtx.MaxMismatches = MaxMismatches;

if((Rslt=LoadProbes(&ProcCtx))<=eBSFSuccess)
	return(Rslt);

// open bioseq file containing chromosome sequences
if((pBioSeq = new CBioSeqFile()) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CBioSeqFile");
	DeleteProbes(&ProcCtx);
	return(eBSFerrObj);
	}
if((Rslt=pBioSeq->Open(pszTargFile,cBSFTypeSeq))!=eBSFSuccess)
	{
	while(pBioSeq->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,pBioSeq->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open input bioseq sequence file '%s'",pszTargFile);
	DeleteProbes(&ProcCtx);
	delete pBioSeq;
	return(Rslt);
	}

// open and truncate CSV output file
#ifdef _WIN32
hRsltsFile = open(pszRsltsFile, _O_WRONLY | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE );
#else
hRsltsFile = open(pszRsltsFile, O_WRONLY | O_CREAT | O_TRUNC, S_IREAD | S_IWRITE );
#endif
if(hRsltsFile == -1)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open output results file '%s', error: %s",pszRsltsFile,strerror(errno));
	DeleteProbes(&ProcCtx);
	delete pBioSeq;
	return(eBSFerrOpnFile);
	}
ProcCtx.hRsltsFile = hRsltsFile;

	// iterate over each entry (chromosome) in target sequence file
pTargSeq = NULL;
TargSeqLen = 0;
AllocdTargSeqLen = 0;
CurEntry = 0;
Rslt = eBSFSuccess;
while((CurEntry = pBioSeq->Next(CurEntry)) > 0)
	{
		// get entry name - assume it will be a chromosome/contig name
	pBioSeq->GetName(CurEntry,sizeof(szChromName),szChromName);
		// check only processing specific chromosomes
	if(pszChroms != NULL && pszChroms[0] != '\0' && !NameInList(pszChroms,szChromName))
		continue;

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing: %s",szChromName);
	ProcCtx.pszCurChrom = szChromName;
	ProcCtx.CurTargID = CurEntry;

		// need to ensure buffer can hold the sequence
	TargSeqLen = pBioSeq->GetDataLen(CurEntry);
	if(pTargSeq == NULL || (TargSeqLen+1) > AllocdTargSeqLen)
		{
		if(pTargSeq != NULL)
			delete pTargSeq;
		AllocdTargSeqLen = TargSeqLen + TargSeqLen/10; // alloc more than actually required, next sequence may not then need to be realloc'd
		pTargSeq = (unsigned char *)new unsigned char [AllocdTargSeqLen];
		if(pTargSeq == NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate memory (%d bytes) for holding chromosome sequence",AllocdTargSeqLen);
			Rslt = eBSFerrMem;
			break;
			}
		}

	// determine and ensure memory allocated for histogram bins to cover chromosome
	// note that each bin contains slots for all probes to be processed
	ProcCtx.NumBins = ((TargSeqLen + BinWidth - 1) / BinWidth);
	if(ProcCtx.pBins == NULL || ProcCtx.AllocBins < ProcCtx.NumBins)
		{
		if(ProcCtx.pBins != NULL)
			delete ProcCtx.pBins;
		ProcCtx.AllocBins = ProcCtx.NumBins + ProcCtx.NumBins/20; // alloc more than actually required, next may not need to be realloc'd
		ProcCtx.pBins = (unsigned int *)new unsigned int [ProcCtx.AllocBins * ProcCtx.NumProbes];
		if(ProcCtx.pBins == NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate memory (%d bytes) for holding histogram bins",ProcCtx.AllocBins * ProcCtx.NumProbes * sizeof(unsigned int));
			Rslt = eBSFerrMem;
			break;
			}
		}
	memset(ProcCtx.pBins,0,ProcCtx.AllocBins * ProcCtx.NumProbes * sizeof(unsigned int));

	pBioSeq->GetData(CurEntry,cBSFTypeSeq,0,pTargSeq,TargSeqLen);
	pTargSeq[TargSeqLen] = eBaseEOS;
	FindMatches(&ProcCtx,		// processing context
		   CurEntry,		// current target entry identifier
		   pTargSeq,		// pts to target sequence
		   TargSeqLen,		// target sequence length
		   szChromName);	// target is on this chromosome 

	CalcChromDensity(&ProcCtx);

	switch(ProcMode) {
		case ePMHChromStats:
			ReportAllChromDensities(&ProcCtx);
			break;
		case ePMHChromBins:
			ReportChromDensities(&ProcCtx);
			break;
		}
	}

DeleteProbes(&ProcCtx);
if(hRsltsFile!=-1)
	close(hRsltsFile);

if(ProcCtx.pBins!=NULL)
	delete ProcCtx.pBins;
if(pTargSeq!=NULL)
	delete pTargSeq;
if(pBioSeq != NULL)
	delete pBioSeq;
return(Rslt);
}


unsigned int *
GetpSlot(int StartLoci,sFSAMProbe *pProbe,tsProcCtx *pProcCtx)
{
return(&pProcCtx->pBins[((StartLoci / pProcCtx->BinWidth) * pProcCtx->NumProbes) + pProbe->ProbID - 1]);
}

// CalcChromDensity
// Calculates densitities for all probes on current chromosome
int
CalcChromDensity(tsProcCtx *pProcCtx)
{
unsigned int *pSlot;
double SumCnts;
double SumCntsSqrd;
double Variance;

sFSAMProbe *pProbe;

SumCnts = 0.0;
SumCntsSqrd = 0.0;
pProbe = pProcCtx->pProbes;
do {
	pSlot = &pProcCtx->pBins[pProbe->ProbID - 1];
	if(pProcCtx->NumBins == 1)
		{
		pProbe->CurChromMean = (double)*pSlot;
		pProbe->CurChromStdDev = 0.0;
		if(pProbe->CurChromMean <= 0.0)
			pProbe->CurChromCOV = 100.0;
		else
			pProbe->CurChromCOV = 0.0;
		}
	else
		{
		for(int BinIdx = 0; BinIdx < pProcCtx->NumBins; BinIdx++, pSlot += pProcCtx->NumProbes)
			{
			SumCnts += *pSlot;
			SumCntsSqrd += *pSlot * *pSlot;
			}
		Variance = (SumCntsSqrd - (((SumCnts * SumCnts)) / pProcCtx->NumBins))/(pProcCtx->NumBins-1);
		pProbe->CurChromMean = SumCnts/pProcCtx->NumBins;
		pProbe->CurChromStdDev = sqrt(Variance);
		if(pProbe->CurChromMean <= 0.0)
			pProbe->CurChromCOV = 100.0;
		else
			pProbe->CurChromCOV = (100.0 * pProbe->CurChromStdDev) / pProbe->CurChromMean;
		}
	}
while((pProbe=pProbe->pNext)!=NULL);
return(eBSFSuccess);
}


int
AddMatch(int StartLoci,sFSAMProbe *pProbe,tsProcCtx *pProcCtx)
{
unsigned int *pSlot = GetpSlot(StartLoci,pProbe,pProcCtx);
*pSlot += 1;
pProbe->TotHits += 1;
pProbe->CurChromHits += 1;
return(eBSFSuccess);
}

void
DeleteProbes(tsProcCtx *pProcCtx)
{
sFSAMProbe *pCurProbe;
sFSAMProbe *pNxtProbe;
if((pNxtProbe = pProcCtx->pProbes) != NULL)
	while((pCurProbe = pNxtProbe) != NULL) 
		{
		pNxtProbe = pCurProbe->pNext;
		delete pCurProbe;
		}
pProcCtx->pProbes = NULL;
pProcCtx->NumProbes = 0;
}


int
LoadProbes(tsProcCtx *pProcCtx)
{
int Rslt;
int UnderLength;
int OverLength;
sFSAMProbe *pCurProbe;
sFSAMProbe *pNxtProbe;
CFasta *pFasta;

char szDescription[cBSFDescriptionSize];
etSeqBase SeqBuff[cMaxProbeLen+1];			// sequences will be terminated by eBaseEOS
int SeqLen;
bool bDescriptor;
bool bFirst;
int TruncCnt;


pProcCtx->NumProbes = 0;
pProcCtx->pProbes = NULL;

// open fasta file containing probes
if((pFasta = new CFasta()) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CFasta");
	return(eBSFerrObj);
	}
if((Rslt=pFasta->Open(pProcCtx->pszProbeFile))!=eBSFSuccess)
	{
	while(pFasta->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,pFasta->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open input fasta file containing probes '%s'",pProcCtx->pszProbeFile);
	delete pFasta;
	return(Rslt);
	}

// read each probe sequence into memory array
// note that only the first cMaxProbeLen of each probe is processed
bDescriptor = false;
bFirst = true;
TruncCnt = 0;
UnderLength = 0;
OverLength = 0;
while((Rslt = SeqLen = pFasta->ReadSequence(SeqBuff,cMaxProbeLen)) > eBSFSuccess)
	{
	if(SeqLen == eBSFFastaDescr)		// just read a descriptor line
		{
		pFasta->ReadDescriptor(szDescription,cBSFDescriptionSize);
		TruncCnt = 0;
		bDescriptor = true;
		bFirst = true;					
		continue;
		}
	else
		if(!bDescriptor)					// if there was no descriptor then dummy up one...
			{
			if(!bFirst)						// only process 1st cMaxProbeLen
				{
				if(!TruncCnt++)
					gDiagnostics.DiagOut(eDLInfo,gszProcName,"LoadProbes [%s] truncated probe %s at %d length",pProcCtx->pszProbeFile,szDescription, cMaxProbeLen);
				continue;
				}
			sprintf(szDescription,"Probe%d",pProcCtx->NumProbes+1);
			TruncCnt = 0;
			bDescriptor = true;
			}
	if(SeqLen < pProcCtx->MinLength)
		{
		UnderLength += 1;
		continue;
		}
	if(SeqLen > pProcCtx->MaxLength)
		{
		OverLength += 1;
		continue;
		}
	bFirst = false;

	SeqBuff[SeqLen] = eBaseEOS;

	if(pProcCtx->MaxProbeLen == 0 || SeqLen > pProcCtx->MaxProbeLen)
		pProcCtx->MaxProbeLen = SeqLen;
	if(pProcCtx->MinProbeLen == 0 || SeqLen < pProcCtx->MinProbeLen)
			pProcCtx->MinProbeLen = SeqLen;
	pNxtProbe = (sFSAMProbe *)new unsigned char [ sizeof(sFSAMProbe) + SeqLen + 1];
	memset(pNxtProbe,0,sizeof(sFSAMProbe));
	pNxtProbe->ProbID = ++pProcCtx->NumProbes;
	pNxtProbe->ProbeLen = SeqLen;
	strncpy(pNxtProbe->szDescr,szDescription,sizeof(pNxtProbe->szDescr)-1);
	pNxtProbe->szDescr[sizeof(pNxtProbe->szDescr)-1] = '\0';
	CSeqTrans::RemoveMasking(SeqBuff,SeqLen);
	memcpy(pNxtProbe->Bases,SeqBuff,SeqLen+1);
	if(pProcCtx->pProbes == NULL)
		pProcCtx->pProbes = pCurProbe = pNxtProbe;
	else
		{
		pCurProbe->pNext = pNxtProbe;
		pCurProbe = pNxtProbe;
		}
	bDescriptor = false;
	}
if(Rslt < eBSFSuccess)
	{
	if(pProcCtx->pProbes != NULL)
		DeleteProbes(pProcCtx);
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessFastaFile [%s] %s",pFasta->ErrText((teBSFrsltCodes)Rslt),pFasta->GetErrMsg());
	return(Rslt);
	}
delete pFasta;

if(UnderLength)
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Filtered out %d underlength probes",UnderLength);
if(OverLength)
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Filtered out %d overlength probes",OverLength);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Number of probes to process %d",pProcCtx->NumProbes);
pProcCtx->NumProbes = pProcCtx->NumProbes;
return(pProcCtx->NumProbes);
}




int										// returns the number of hits by probes against the chromosome sequence
FindMatches(tsProcCtx *pProcCtx,			// processing parameters
		   tBSFEntryID TargEntryID,		// current target entry identifier
		   unsigned char *pTargSeq,		// pts to target sequence
		   int TargSeqLen,				// target sequence length
		   char *pszChrom)				// target is on this chromosome 
{
int Rslt;
sFSAMProbe *pPutProbe;					// putative matching probe
int CurTargIdx;
int TargEndIdx;
int TargStartIdx;
etSeqBase *pProbe;
etSeqBase *pTarg;
etSeqBase ProbeBase;
etSeqBase TargBase;
int	TotMatches;
int	TotMismatches;

int PlusHits;
int MinHits;
int NumCompared;

int StartLoci;

PlusHits = 0;
MinHits = 0;

// reset probe chromosome stats
if((pPutProbe = pProcCtx->pProbes)==NULL)
	return(0);
do {
	pPutProbe->CurChromCOV = 0.0;
	pPutProbe->CurChromHits = 0;
	pPutProbe->CurChromMean = 0.0;
	pPutProbe->CurChromStdDev = 0.0;
	}
while((pPutProbe = pPutProbe->pNext)!=NULL);

// first process the '+' strand
if((pProcCtx->Strand == '*' || pProcCtx->Strand == '+'))
	{
	for(CurTargIdx = 0; CurTargIdx < TargSeqLen; CurTargIdx++)
		{
		if(!(CurTargIdx & 0x0fffff))
			printf("\rTarget Loci '+' strand %s: %9.9d",pszChrom,CurTargIdx);		// let user know we haven't crashed...
		pPutProbe = pProcCtx->pProbes;
		do
			{
			TargStartIdx = CurTargIdx;
			TargEndIdx = CurTargIdx;
			TotMatches = 0;
			TotMismatches = 0;


			pProbe = pPutProbe->Bases;
			pTarg = (etSeqBase *)&pTargSeq[CurTargIdx];
			NumCompared = 0;
			while((ProbeBase = (*pProbe++ & ~cRptMskFlg)) != eBaseEOS && (TargBase = (*pTarg++ & ~cRptMskFlg)) != eBaseEOS)
				{
				NumCompared += 1;
				if(ProbeBase == TargBase)
					TotMatches++;
				else			
					if(++TotMismatches > pProcCtx->MaxMismatches)
						break;
				TargEndIdx += 1;
				}

			if(NumCompared == pPutProbe->ProbeLen && 
					TotMismatches <= pProcCtx->MaxMismatches)
				{
				if((Rslt=AddMatch(TargStartIdx,pPutProbe,pProcCtx))>=eBSFSuccess)
					PlusHits += 1;
				}
			}
		while((pPutProbe = pPutProbe->pNext)!=NULL);
		}
	printf("\rTarget Loci '+' strand %s: %9.9d",pszChrom,TargSeqLen);
	}


// and now the '-' strand
if((pProcCtx->Strand == '*' || pProcCtx->Strand == '-'))
	{
	printf("\n");
	CSeqTrans::ReverseComplement(TargSeqLen,pTargSeq);
	for(CurTargIdx = 0; CurTargIdx < TargSeqLen; CurTargIdx++)
		{
		if(!(CurTargIdx & 0x0fffff))
			printf("\rTarget Loci '-' strand %s: %9.9d",pszChrom,CurTargIdx);		// let user know we haven't crashed...
		pPutProbe = pProcCtx->pProbes;
		do
			{
			TargStartIdx = CurTargIdx;
			TargEndIdx = CurTargIdx;
			TotMatches = 0;
			TotMismatches = 0;

			pProbe = pPutProbe->Bases;
			pTarg = (etSeqBase *)&pTargSeq[CurTargIdx];
			NumCompared = 0;
			while((ProbeBase = (*pProbe++ & ~cRptMskFlg)) != eBaseEOS && (TargBase = (*pTarg++ & ~cRptMskFlg)) != eBaseEOS)
				{
				NumCompared += 1;
				if(ProbeBase == TargBase)
					TotMatches++;
				else			
					if(++TotMismatches > pProcCtx->MaxMismatches)
						break;
				TargEndIdx += 1;
				}

			if(NumCompared == pPutProbe->ProbeLen && 
				TotMismatches <= pProcCtx->MaxMismatches)
				{
				StartLoci = TargSeqLen - TargEndIdx;
				if((Rslt=AddMatch(StartLoci,pPutProbe,pProcCtx))>=eBSFSuccess)
					MinHits += 1;
				}
	
			}
		while((pPutProbe = pPutProbe->pNext)!=NULL);
		}
	printf("\rTarget Loci '-' strand %s: %9.9d",pszChrom,TargSeqLen);
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Hits on plus strand: %d, on minus strand: %d",PlusHits,MinHits);
return(PlusHits + MinHits);
}

int 
ReportChromDensities(tsProcCtx *pProcCtx)
{
unsigned int *pSlot;
sFSAMProbe *pProbe;
int Len;
char szLineBuff[16000];
static bool bFirst = true;

Len = 0;
if(bFirst)
	{
	Len = sprintf(szLineBuff,"\"Bin\",\"Chrom\",\"TargID\",\"Probe\",\"ProbeID\",\"Hits\"\n");
	bFirst = false;
	}

pProbe = pProcCtx->pProbes;
do {
	pSlot = &pProcCtx->pBins[pProbe->ProbID - 1];
	for(int BinIdx = 0; BinIdx < pProcCtx->NumBins; BinIdx++, pSlot += pProcCtx->NumProbes)
		{
		Len += sprintf(&szLineBuff[Len],"%d,\"%s\",%d,\"%s\",%d,%d\n",
			BinIdx+1,pProcCtx->pszCurChrom,pProcCtx->CurTargID,pProbe->szDescr,pProbe->ProbID,*pSlot);
		if(Len > (sizeof(szLineBuff) * 9) / 10)
			{
			CUtility::SafeWrite(pProcCtx->hRsltsFile,szLineBuff,Len);
			Len = 0;
			}		
		}
	}
while((pProbe=pProbe->pNext)!=NULL);
if(Len)
	CUtility::SafeWrite(pProcCtx->hRsltsFile,szLineBuff,Len);
#ifdef _WIN32
	_commit(pProcCtx->hRsltsFile);
#else
	fsync(pProcCtx->hRsltsFile);
#endif
return(eBSFSuccess);
}

int 
ReportAllChromDensities(tsProcCtx *pProcCtx)
{
sFSAMProbe *pProbe;
int Len;
char szLineBuff[16000];
static bool bFirst = true;

Len = 0;
if(bFirst)
	{
	Len = sprintf(szLineBuff,"\"Chrom\",\"TargID\",\"Probe\",\"ProbeID\",\"COV\",\"StdDev\",\"Mean\",\"Hits\"\n");
	bFirst = false;
	}

pProbe = pProcCtx->pProbes;
do {
	
   Len += sprintf(&szLineBuff[Len],"\"%s\",%d,\"%s\",%d,%1.8f,%1.8f,%1.8f,%d\n",
			pProcCtx->pszCurChrom,pProcCtx->CurTargID,pProbe->szDescr,pProbe->ProbID,pProbe->CurChromCOV,pProbe->CurChromStdDev,pProbe->CurChromMean,pProbe->CurChromHits);
	if(Len > (sizeof(szLineBuff) * 9) / 10)
		{
		CUtility::SafeWrite(pProcCtx->hRsltsFile,szLineBuff,Len);
		Len = 0;
		}		
	}
while((pProbe=pProbe->pNext)!=NULL);
if(Len)
	CUtility::SafeWrite(pProcCtx->hRsltsFile,szLineBuff,Len);
#ifdef _WIN32
	_commit(pProcCtx->hRsltsFile);
#else
	fsync(pProcCtx->hRsltsFile);
#endif
return(eBSFSuccess);
}

