// csv2fasta.cpp : Defines the entry point for the console application.
// Generates either a multifasta or concatenated fasta file containing all sequences specified in a CSV loci file


#include "stdafx.h"
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#if _WIN32

#include "../libbiokanga/commhdrs.h"
#else
#include "../libbiokanga/commhdrs.h"
#endif
#include "IncExclChroms.h"

const unsigned int cProgVer = 315;		// increment with each release

const int cMaxTruncLength = 50000000;	// maximal element sequence truncation length
const int cDfltTruncLength= 1000000;	// by default truncate sequences to this length
const int cDfltMinLengthRange = 10;		// default minimum element length
const int cMaxLenCSVSeq = 32000;		// maximal element length when sequence is from adhoc CSV file 

const int cMaxAllocBuffChunk = 0x007ffffff;	// buffer for fasta sequences is realloc'd in this sized chunks

typedef struct TAG_sProcParams 
	{
	bool bCapsSoftMask;		// treat uppercase as softmasked bases
	CBioSeqFile *pSeqFile;	// preopened bioseq file
	} tsProcParams; 

typedef enum TAG_eRsltsFormat {
	eRsltsMultifasta = 0,				// output sequences as multifasta
	eRsltsFasta,						// output sequences concatenated as single fasta record with no sequence separators
	eRsltsRMESFasta						// output sequences concatenated in single fasta record but with 'Z' as sequence separator char
	} teRsltsFormat;

int 
Process(teCSVFormat CSVFormat,		// expected input CSV format
		teRsltsFormat RsltsFormat,	// results output mode (default 0) 0: multifasta, 1: concatenated fasta 
		char Strand,				// specifies '+'/'-'/'*' strand to process
		int OfsLoci,				// offset start loci by this many bases
		int DeltaLen,				// delta length by this many bases
		bool bSkipFirst,			// true if first line contains header and should be skipped
		int MinLength,				// core elements must be of at least this length
		int TruncLength,			// truncate sequences to be no longer than this length
		int	NumIncludeChroms,		// number of chromosome regular expressions to include
		char **ppszIncludeChroms,	// array of include chromosome regular expressions
		int	NumExcludeChroms,		// number of chromosome expressions to exclude
		char **ppszExcludeChroms,	// array of exclude chromosome regular expressions
		int MinIdentity,			// minimum identity (0..100)
		int MaxIdentity,			// maximum identity (MinIdentity..100)
		int	DescrIDX,				// which field to use as the descriptor
		int SeqIDX,					// which field contains the sequence
		char *pszInLociFile,		// CSV file containing elements
		char *pszInSeqFile,			// file containing genome assembly sequences
		char *pszRsltsFile);		// file to write results into

// ProcessFastaFile
// Parse input fasta format file into a biosequence file
bool ProcessFastaFile(char *pszFile,
				 void *pParams);

int
CreateBioseqFastaFile(bool bTargDeps,			// true if process only if any independent src files newer than target
					  char *pszSrcDirPath,
					  char *pszDestBioseqFile,
					  char *pszRefSpecies,
					  char *pszDescr,char *pszTitle,
					  bool bCapsSoftMask);

char *RsltsFormat2Text(teRsltsFormat Format);
char *CSVFormat2Text(teCSVFormat Format);


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
	return _T("csv2fasta");
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

int Rslt;
int Idx;
int LenChromList;
int CSVFormat;
int RsltsFormat;
char Strand;				// if '-' or '+' then only process that specific strand
bool bSkipFirst;			// true if first line contains header and should be skipped
int MinLength;				// core elements must be of at least this length
int TruncLength;			// and truncate to be no longer than this length
int NumIncludeChroms;
char *pszIncludeChroms[cMaxIncludeChroms];
int NumExcludeChroms;
char *pszExcludeChroms[cMaxExcludeChroms];

int MinIdentity;			// out species must align with at least this identity
int MaxIdentity;			// and no more than this identity
int	DescrIDX;				// which field to use as the descriptor
int SeqIDX;				// which field contains the sequence
int OfsLoci;				// offset loci by this many bases
int DeltaLen;				// change element lengths by this many bases
char szInLociFile[_MAX_PATH];	// input element loci from this file
char szInSeqFile[_MAX_PATH];	// input bioseq file containing assembly
char szRsltsFile[_MAX_PATH];	// output loci + sequences or gene identifiers to this file


// command line args
struct arg_lit  *help    = arg_lit0("hH","help",                "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_int *ScreenLogLevel=arg_int0("S", "ScreenLogLevel",	"<int>","Level of diagnostics written to logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"diagnostics log file");

struct arg_int  *csvformat = arg_int0("m","inpformat","<int>",	"input CSV file type 0:Loci 1:locsfx probe 2:locsfx target 3:quickcount 4: Genhyperconserved x2 5: GenMAlignScore m3 6:arbitary seq");
struct arg_file *inlocifile = arg_file1("i","inloci","<file>",	"sequence loci CSV file");
struct arg_file *inseqfile = arg_file0("I","assembly","<file>",	"genome assembly bioseq or fasta file");
struct arg_file *rsltsfile = arg_file1("o","output","<file>",	"output file");
struct arg_int  *rsltsformat = arg_int0("p","outformat","<int>","results output as 0: multifasta,  1: concatenated fasta, 2: RMES format fasta");
struct arg_lit  *skipfirst = arg_lit0("x","skipfirst",          "skip first line of CSV - header line");
struct arg_str *strand=arg_str0("s", "strand","<str>",          "process for this strand '+' or '-' only (default is to process both)");

struct arg_int  *minlength = arg_int0("l","minlength","<int>",	"minimum element length (default 10)");
struct arg_int  *trunclength = arg_int0("T","truncatelength","<int>","truncate elements to be no  longer than this length (default 1000000)");
struct arg_int  *ofsloci   = arg_int0("u","offset","<int>",	    "offset element loci by this many bases, -2048 to +2048 (default 0)");
struct arg_int  *deltalen  = arg_int0("U","deltalen","<int>",	"delta element lengths by this many bases, -2048 to +2048 (default 0)");
struct arg_str  *ExcludeChroms = arg_strn("Z","chromexclude","<string>",0,cMaxExcludeChroms,"high priority - regular expressions defining species.chromosomes to exclude");
struct arg_str  *IncludeChroms = arg_strn("z","chromeinclude","<string>",0,cMaxIncludeChroms,"low priority - regular expressions defining species.chromosomes to include");

struct arg_int *minidentity = arg_int0("k", "minident","<int>",	"input type 3 minimum out species alignment identity (default 0)");
struct arg_int *maxidentity = arg_int0("K", "maxident","<int>",	"input type 3 maximum out species alignment identity (default 100)");
struct arg_int  *descridx = arg_int0("d","descridx","<int>",	"input type 4 index (1..100) of field in CSV containing value to use in fasta descriptor line");
struct arg_int  *seqidx = arg_int0("q","seqidx","<int>",		"input type 4 index (1..100) of field in CSV containing sequence");

struct arg_end *end = arg_end(20);

void *argtable[] = {help,version,FileLogLevel,ScreenLogLevel,LogFile,
					csvformat,inlocifile,inseqfile,rsltsfile,rsltsformat,skipfirst,strand,ofsloci,deltalen,minlength,trunclength,ExcludeChroms,IncludeChroms,minidentity,maxidentity,descridx,seqidx,
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

	CSVFormat = csvformat->count ? csvformat->ival[0] : eCSVFdefault;
	if(CSVFormat < eCSVFdefault || CSVFormat > eCSVFSeq)
		{
		printf("\nError: expected input CSV format specified '-m%d' must be in range 0..6",CSVFormat);
		exit(1);
		}

	RsltsFormat = rsltsformat->count ? rsltsformat->ival[0] : eRsltsMultifasta;
	if(RsltsFormat < eRsltsMultifasta || RsltsFormat > eRsltsRMESFasta)
		{
		printf("\nError: RsltsFormat specified '-p%d' must be in range 0..1",RsltsFormat);
		exit(1);
		}

	Strand = '*';
	if(CSVFormat == eCSVFdefault)
		{
		Strand = strand->count ? *(char *)strand->sval[0] : '*';
		if(!(Strand == '+' || Strand == '-' || Strand == '*'))
			{
			printf("\nError: Strand specified '-s%c' must be one of '+', '-' or '*'",Strand);
			exit(1);
			}
		}

	bSkipFirst = skipfirst->count ? true : false;

	OfsLoci = ofsloci->count ? ofsloci->ival[0] : 0;
	if(abs(OfsLoci) > 1024)
		{
		printf("Error: loci offset '-u%d' must be in the range -2048 to +2048",OfsLoci);
		exit(1);
		}

	DeltaLen = deltalen->count ? deltalen->ival[0] : 0;
	if(abs(DeltaLen) > 1024)
		{
		printf("Error: delta length '-U%d' must be in the range -2048 to +2048",DeltaLen);
		exit(1);
		}

	MinLength = minlength->count ? minlength->ival[0] : cDfltMinLengthRange;
	if(MinLength < 1 || MinLength > cMaxTruncLength)
		{
		printf("Error: Mininum element length '-l%d' is not in range 1..%d",MinLength,cMaxTruncLength);
		exit(1);
		}

	TruncLength = trunclength->count ? trunclength->ival[0] : cDfltTruncLength;
	if(TruncLength < MinLength || TruncLength > cMaxTruncLength)
		{
		printf("Error: Maximum element length '-T%d' is not in range %d..%d",TruncLength,MinLength,cMaxTruncLength);
		exit(1);
		}

	NumIncludeChroms = IncludeChroms->count;
	for(Idx=0;Idx < IncludeChroms->count; Idx++)
		{
		LenChromList = (int)strlen(IncludeChroms->sval[Idx]);
		pszIncludeChroms[Idx] = new char [LenChromList+1];
		strcpy(pszIncludeChroms[Idx],IncludeChroms->sval[Idx]);
		CUtility::TrimQuotes(pszIncludeChroms[Idx]);
		}

	NumExcludeChroms = ExcludeChroms->count;
	for(Idx=0;Idx < ExcludeChroms->count; Idx++)
		{
		LenChromList = (int)strlen(ExcludeChroms->sval[Idx]);
		pszExcludeChroms[Idx] = new char [LenChromList+1];
		strcpy(pszExcludeChroms[Idx],ExcludeChroms->sval[Idx]);
		CUtility::TrimQuotes(pszExcludeChroms[Idx]);
		}

	if(CSVFormat == eCSVFhyperX2)
		{
		MinIdentity = minidentity->count ? minidentity->ival[0] : 0;
		if(MinIdentity < 0 || MinIdentity > 100)
			{
			printf("\nError: Minimum identity specified as '-d%d' must be in range 0..100\n",MinIdentity);
			exit(1);
			}

		MaxIdentity = maxidentity->count ? maxidentity->ival[0] : 100;
		if(MaxIdentity < MinIdentity || MaxIdentity > 100)
			{
			printf("\nError: Maximum identity specified as '-D%d' must be in range %d..100\n",MaxIdentity,MinIdentity);
			exit(1);
			}
		}
	else
		{
		MinIdentity = 0;
		MaxIdentity = 100;
		}

	if(CSVFormat == eCSVFSeq)
		{
		DescrIDX = descridx->ival[0];
		if(DescrIDX < 1 || DescrIDX > 100)
			{
			printf("\nError: Descriptor field index '-d%d' specified outside of range 1..100",DescrIDX);
			exit(1);
			}
		SeqIDX = seqidx->ival[0];
		if(SeqIDX < 1 || SeqIDX > 100)
			{
			printf("\nError: Sequence field index '-d%d' specified outside of range 1..100",SeqIDX);
			exit(1);
			}
		}
	else
		{
		DescrIDX = 0;
		SeqIDX = 0;
		}

	strncpy(szInLociFile,inlocifile->filename[0],_MAX_PATH);
	szInLociFile[_MAX_PATH-1] = '\0';

	if(CSVFormat == eCSVFSeq || CSVFormat == eCSVFQuickCount)
		szInSeqFile[0] = '\0';
	else
		{
		if(inseqfile->count == 0)
			{
			printf("\nError: No input bioseq assembly file specified with '-I<file'");
			exit(1);
			}
		strncpy(szInSeqFile,inseqfile->filename[0],_MAX_PATH);
		szInSeqFile[_MAX_PATH-1] = '\0';
		}
	strncpy(szRsltsFile,rsltsfile->filename[0],_MAX_PATH);
	szRsltsFile[_MAX_PATH-1] = '\0';

		// now that command parameters have been parsed then initialise diagnostics log system
	if(!gDiagnostics.Open(szLogFile,(etDiagLevel)iScreenLogLevel,(etDiagLevel)iFileLogLevel,true))
		{
		printf("\nError: Unable to start diagnostics subsystem.");
		if(szLogFile[0] != '\0')
			printf(" Most likely cause is that logfile '%s' can't be opened/created",szLogFile);
		exit(1);
		}
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Version: %d.%2.2d Processing parameters:",cProgVer/100,cProgVer%100);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Input CSV element loci file: '%s'",szInLociFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Expecting CSV format as: %s",CSVFormat2Text((teCSVFormat)CSVFormat));
	if(CSVFormat == eCSVFdefault)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Process for this strand : '%c'",Strand);
	if(!(CSVFormat == eCSVFSeq || CSVFormat == eCSVFQuickCount))
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Input bioseq genome assembly file: '%s'",szInSeqFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Output to file: '%s'",szRsltsFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Generating Results as: %s",RsltsFormat2Text((teRsltsFormat)RsltsFormat));
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"First line contains header: %s",bSkipFirst ? "yes" : "no");
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Offset element start loci by : %d",OfsLoci);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Delta element length by : %d",DeltaLen);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Minimum element length: %d",MinLength);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Truncate element length: %d",TruncLength);

	for(Idx = 0; Idx < NumIncludeChroms; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"reg expressions defining chroms to include: '%s'",pszIncludeChroms[Idx]);
	for(Idx = 0; Idx < NumExcludeChroms; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"reg expressions defining chroms to exclude: '%s'",pszExcludeChroms[Idx]); 

	if(CSVFormat == eCSVFhyperX2)
		{
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Minimum out species identity: %d",MinIdentity);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Maximum out species identity: %d",MaxIdentity);
		}
	if(CSVFormat == eCSVFSeq)
		{
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Descriptor field index: %d",DescrIDX);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Sequence field index: %d",SeqIDX);
		}

#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	gStopWatch.Start();
	// processing here...
	Rslt = Process((teCSVFormat)CSVFormat,(teRsltsFormat)RsltsFormat,Strand,OfsLoci,DeltaLen,bSkipFirst,MinLength,TruncLength,NumIncludeChroms,pszIncludeChroms,NumExcludeChroms,pszExcludeChroms,MinIdentity,MaxIdentity,DescrIDX,SeqIDX,szInLociFile,szInSeqFile,szRsltsFile);

	gStopWatch.Stop();
	Rslt = Rslt >=0 ? 0 : 1;
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Exit code: %d Total processing time: %s",Rslt,gStopWatch.Read());
	exit(Rslt);
	}
else
	{
	arg_print_errors(stdout,end,gszProcName);
	arg_print_syntax(stdout,argtable,"\nend of help\n");
	exit(1);
	}
}

char *
RsltsFormat2Text(teRsltsFormat Format)
{
switch(Format) {
	case eRsltsMultifasta:	
		return((char *)"Sequences as multiple fasta records");
	case eRsltsFasta:	
		return((char *)"Sequences concatenated into single Fasta record with no sequence separators");
	case eRsltsRMESFasta:
		return((char *)"Sequences concatenated into single Fasta record with 'Z' as a sequence separator (for R'MES package input)");
	default:
		break;
	}
return((char *)"Unsupported");
}

char *
CSVFormat2Text(teCSVFormat Format)
{
switch(Format) {
	case eCSVFdefault:
		return((char *)"Default 9+ field loci");
	case eCSVFprobe:
		return((char *)"locsfx probe loci");
	case eCSVFtarget:
		return((char *)"locsfx target loci");
	case eCSVFQuickCount:
		return((char *)"quickcount COV density");
	case eCSVFhyperX2:
		return((char *)"Genhyperconserved with -x2 option");
	case eCSVFAlignM3:
		return((char *)"GenMAlignScore with -m3 option");
	case eCSVFSeq:
		return((char *)"CSV adhoc format in which one of the fields contains a sequence and another contains a fasta descriptor");
	default:
		break;
	}
return((char *)"Unsupported");
}

int			// actual length, or if under MinLen then -1, or if start would be < 0 or end < start then returns -2
AdjLoci(bool bOnMinStrand,		// true if element is on '-' strand
		int *pStartLoci,		// starting loci on '+' strand
		int *pEndLoci,			// ending loci on '+' strand
		int OfsLoci,			// offset loci by
		int DeltaLen,			// adjust length by
		int MinLen,				// must be at least this length
		int TruncLen,			// truncate to this length
		int *pTruncCnts)		// will be incremented if truncation was needed
{
int Ofs;
int DLen;

if(*pStartLoci < 0 || *pStartLoci > *pEndLoci)
	return(-2);

if(OfsLoci != 0 || DeltaLen != 0)
	{
	if(bOnMinStrand)
		{
		Ofs = OfsLoci * -1;
		DLen = DeltaLen * -1;
		}
	else
		{
		Ofs = OfsLoci;
		DLen = DeltaLen;
		}

	*pStartLoci += Ofs;
	*pEndLoci += Ofs;
	if(bOnMinStrand)
		*pStartLoci += DLen;
	else
		*pEndLoci += DLen;
	}

if(*pStartLoci > *pEndLoci || (!bOnMinStrand && *pStartLoci < 0))
	return(-2);
DLen = 1 + *pEndLoci - *pStartLoci;
if(DLen < MinLen)
	return(-1);
if(DLen > TruncLen)
	{
	if(bOnMinStrand)
		{
		*pStartLoci = 1 + *pEndLoci - TruncLen;
		if(*pStartLoci < 0)
			return(-2);
		}
	else
		*pEndLoci = *pStartLoci + TruncLen - 1;
	DLen= TruncLen;
	*pTruncCnts += 1;
	}
return(DLen);
}

int 
Process(teCSVFormat CSVFormat,		// expected input CSV format
		teRsltsFormat RsltsFormat,	// output format (default 0) 0: multifasta, 1: concatenated fasta
		char Strand,				// specifies '+'/'-'/'*' strand to process
		int OfsLoci,				// offset start loci by this many bases
		int DeltaLen,				// delta length by this many bases
		bool bSkipFirst,			// true if first line contains header and should be skipped
		int MinLength,				// core elements must be of at least this length
		int TruncLength,			// truncate elements to at most this length
		int	NumIncludeChroms,		// number of chromosome regular expressions to include
		char **ppszIncludeChroms,	// array of include chromosome regular expressions
		int	NumExcludeChroms,		// number of chromosome expressions to exclude
		char **ppszExcludeChroms,	// array of exclude chromosome regular expressions
		int MinIdentity,			// minimum identity (0..100)
		int MaxIdentity,			// maximum identity (MinIdentity..100)
		int	DescrIDX,				// which field to use as the descriptor
		int SeqIDX,					// which field contains the sequence
		char *pszInLociFile,		// CSV file containing elements
		char *pszInSeqFile,			// file containing genome assembly sequences
		char *pszRsltsFile)			// file to write fasta into
{
int NumFields;

int Rslt;
int SrcID;
char *pszChrom;
int ChromID;
char *pszElType;
char *pszRefSpecies;
char *pszRelSpecies;
char *pszStrand;
int StartLoci;
int EndLoci;
int Len;
bool bOnMinStrand;
int NumReads;

int OGUnalignedBases;
int OGMatchCnt;
int OGMismatchCnt;
int OGInDelCnt;

int Identity;
int NumCols;
int SeqOfs;
int NxtFastaCol;
bool bFastaStart;

char szDescr[cMaxGenFastaLineLen+1];		// to hold fasta descriptors
char *pszDescr;
char *pszSeq;
int ReqAllocLen;
int AllocSeqBuffLen = 0;
int AllocLineBuffLen = 0;

etSeqBase *pElSeq;
char Chr;
int NumElsProcessed;
int NumElsAccepted;
int NumUnderLen;
int NumTruncated;
int NumUnderIdentity;
int NumOverIdentity;
int NumLoadSeqFails;
int NumStrandSkipped;
int NumExclChroms;

CIncExclChroms IncExclChroms;

if((Rslt=IncExclChroms.InitChromExclusion(NumIncludeChroms,ppszIncludeChroms,NumExcludeChroms,ppszExcludeChroms))!=eBSFSuccess)
	return(Rslt);

char *pszLineBuff = NULL;
int BuffLen;
etSeqBase *pElSeqBuff = NULL;
int hRsltFile = -1;
CBioSeqFile *pBioseq = NULL;

CCSVFile *pCSV = new CCSVFile;
if(pCSV == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CCSVfile");
	return(eBSFerrObj);
	}

if(CSVFormat == eCSVFSeq || CSVFormat == eCSVFQuickCount)
	{
	pCSV->SetMaxFieldLen(cMaxLenCSVSeq);
	pCSV->SetMaxFields(max(DescrIDX,SeqIDX));
	}

if((Rslt=pCSV->Open(pszInLociFile))!=eBSFSuccess)
	{
	while(pCSV->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,pCSV->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open file: %s",pszInLociFile);
	delete pCSV;
	return(Rslt);
	}

if(CSVFormat != eCSVFSeq && CSVFormat != eCSVFQuickCount)
	{
	// assume existing bioseq file, if not then try to create from what must have been a fasta file
	if((pBioseq = new CBioSeqFile) == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CBioSeqFile object");
		delete pCSV;
		return(eBSFerrObj);
		}

	Rslt = pBioseq->Open(pszInSeqFile);
	if(Rslt != eBSFSuccess)
		{
		char szBioseqFile[_MAX_PATH];
		strcpy(szBioseqFile,pszInSeqFile);
		strcat(szBioseqFile,".bsb");
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open BSB file, creating %s from input assumed fasta file %s",szBioseqFile,pszInSeqFile);
		if((Rslt = CreateBioseqFastaFile(false,pszInSeqFile,szBioseqFile,(char *)"test",(char *)"test",(char *)"test",(char *)false))>0)
			Rslt = pBioseq->Open(szBioseqFile);
		}

	if(Rslt != eBSFSuccess)
		{
		while(pBioseq->NumErrMsgs())
			gDiagnostics.DiagOut(eDLFatal,gszProcName,pBioseq->GetErrMsg());
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open assembly sequence file '%s'",pszInSeqFile);
		delete pCSV;
		delete pBioseq;
		return(Rslt);
		}
	}
else
	pBioseq = NULL;

#ifdef _WIN32
if((hRsltFile = open(pszRsltsFile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
if((hRsltFile = open(pszRsltsFile, O_RDWR | O_CREAT |O_TRUNC,S_IREAD | S_IWRITE))==-1)
#endif
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create %s - %s",pszRsltsFile,strerror(errno));
	delete pCSV;
	if(pBioseq != NULL)
		delete pBioseq;
	return(eBSFerrCreateFile);
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Output file created/truncated: '%s'",pszRsltsFile);
bFastaStart = true;
NumExclChroms = 0;
NumUnderLen = 0;
NumTruncated = 0;
NumUnderIdentity = 0;
NumOverIdentity = 0;
NumLoadSeqFails = 0;
NumElsProcessed = 0;
NumElsAccepted = 0;
NumStrandSkipped = 0;
bOnMinStrand = false;				// assume always on plus strand unless strand is present in CSV
while((Rslt=pCSV->NextLine()) > 0)	// onto next line containing fields
	{
	if(bSkipFirst || (!NumElsProcessed && pCSV->IsLikelyHeaderLine()))
		{
		bSkipFirst = false;
		continue;
		}

	NumElsProcessed++;
	NumFields = pCSV->GetCurFields();
	switch(CSVFormat) {
		case eCSVFdefault:
		case eCSVFprobe:
			if(NumFields < 7)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expected 7+ fields in '%s', GetCurFields() returned '%d'",pszInLociFile,NumFields);
				Rslt = eBSFerrFieldCnt;
				break;
				}

			pCSV->GetText(4,&pszChrom);
			if((Rslt = IncExclChroms.IncludeThisChrom(pszChrom)) < 1)
				{
				if(Rslt < 0)
					break;
				NumExclChroms += 1;
				continue;
				}


			pCSV->GetInt(1,&SrcID);
			pCSV->GetText(2,&pszElType);
			pCSV->GetText(3,&pszRefSpecies);

			pCSV->GetInt(5,&StartLoci);
			pCSV->GetInt(6,&EndLoci);
			pCSV->GetInt(7,&Len);
			if(StartLoci > EndLoci)
				{
				int tmp = StartLoci;
				StartLoci = EndLoci;
				EndLoci = tmp;
				}
			bOnMinStrand = false;
			if(NumFields >= 8)					// check if strand has been specified
				{
				pCSV->GetText(8,&pszStrand);
				if(pszStrand[0] == '-')		// assume anything other than '-' is on the plus strand
					bOnMinStrand = true;
				}

			if((Strand == '+' && bOnMinStrand) || (Strand == '-' && !bOnMinStrand))
				{
				NumStrandSkipped += 1;
				continue;
				}

			if((Len = AdjLoci(bOnMinStrand,&StartLoci,&EndLoci,OfsLoci,DeltaLen,MinLength,TruncLength,&NumTruncated)) < 0)
				{
				NumUnderLen += 1;
				continue;
				}

			if(NumFields >= 11)					// check if number of reads has been specified
				pCSV->GetInt(11,&NumReads);
			else
				NumReads = 1;

			sprintf(szDescr,"%d %s|%s|%d|%d|%d|%c|%d",
				SrcID,pszRefSpecies,pszChrom,StartLoci,EndLoci,Len,bOnMinStrand ? '-' : '+',NumReads);
			break;


		case eCSVFtarget:				// use target loci
			if(NumFields < 12)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expected 12+ fields in '%s', GetCurFields() returned '%d'",pszInLociFile,NumFields);
				Rslt = eBSFerrFieldCnt;
				break;
				}

			pCSV->GetText(9,&pszChrom);
			if((Rslt = IncExclChroms.IncludeThisChrom(pszChrom))<1)
				{
				if(Rslt < 0)
					break;
				NumExclChroms += 1;
				continue;
				}


			pCSV->GetInt(7,&Len);
			pCSV->GetInt(1,&SrcID);
			pCSV->GetText(2,&pszElType);
			pCSV->GetText(8,&pszRefSpecies);

			pCSV->GetInt(10,&StartLoci);
			pCSV->GetInt(11,&EndLoci);
			pCSV->GetText(3,&pszRelSpecies);
			if((Len = AdjLoci(bOnMinStrand,&StartLoci,&EndLoci,OfsLoci,DeltaLen,MinLength,TruncLength,&NumTruncated)) < 0)
				{
				NumUnderLen += 1;
				continue;
				}

			sprintf(szDescr,"%d %s|%s|%d|%d|%d",
					SrcID,pszRefSpecies,pszChrom,StartLoci,EndLoci,Len);
			break;
	
		case eCSVFAlignM3:			// GenMAlignScore with -m3 option
			if(NumFields < 7)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expected 7+ fields in '%s', GetCurFields() returned '%d'",pszInLociFile,NumFields);
				Rslt = eBSFerrFieldCnt;
				break;
         		}

			pCSV->GetText(4,&pszChrom);
			if((Rslt = IncExclChroms.IncludeThisChrom(pszChrom))<1)
				{
				if(Rslt < 0)
					break;
				NumExclChroms += 1;
				continue;
				}


			pCSV->GetInt(2,&Len);
			pCSV->GetInt(1,&SrcID);
			pszElType = (char *)"BlockID";
			pCSV->GetText(3,&pszRefSpecies);

			pCSV->GetInt(6,&StartLoci);
			pCSV->GetInt(7,&EndLoci);
			pszRelSpecies = pszRefSpecies;
			if((Len = AdjLoci(bOnMinStrand,&StartLoci,&EndLoci,OfsLoci,DeltaLen,MinLength,TruncLength,&NumTruncated)) < 0)
				{
				NumUnderLen += 1;
				continue;
				}
			sprintf(szDescr,"%d %s|%s|%d|%d|%d",
					SrcID,pszRefSpecies,pszChrom,StartLoci,EndLoci,Len);
			break;

		case eCSVFhyperX2:
			if(NumFields < 13)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expected at least 13 fields in '%s', GetCurFields() returned '%d'",pszInLociFile,NumFields);
				Rslt = eBSFerrFieldCnt;
				break;
				}

			pCSV->GetText(4,&pszChrom);
			if((Rslt = IncExclChroms.IncludeThisChrom(pszChrom))<1)
				{
				if(Rslt < 0)
					break;
				NumExclChroms += 1;
				continue;
				}

			pCSV->GetInt(7,&Len);

			pCSV->GetInt(11,&OGMatchCnt);
			Identity = (int)(((100.0 * OGMatchCnt) / (double)Len) + 0.001);

			if(Identity < MinIdentity)
				{
				NumUnderIdentity += 1;
				continue;
				}
		
			if(Identity > MaxIdentity)
				{
				NumOverIdentity += 1;
				continue;
				}

			pCSV->GetInt(1,&SrcID);
			pCSV->GetText(2,&pszElType);
			pCSV->GetText(3,&pszRefSpecies);

			pCSV->GetInt(5,&StartLoci);
			pCSV->GetInt(6,&EndLoci);
			if((Len = AdjLoci(bOnMinStrand,&StartLoci,&EndLoci,OfsLoci,DeltaLen,MinLength,TruncLength,&NumTruncated)) < 0)
				{
				NumUnderLen += 1;
				continue;
				}
			pCSV->GetText(8,&pszRelSpecies);

			pCSV->GetInt(10,&OGUnalignedBases);
			pCSV->GetInt(12,&OGMismatchCnt);
			pCSV->GetInt(13,&OGInDelCnt);
			sprintf(szDescr,"%d %s|%s|%d|%d|%d",
					SrcID,pszRefSpecies,pszChrom,StartLoci,EndLoci,Len);
			break;

		case eCSVFQuickCount:
			if(NumFields < 6)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expected 6+ fields in '%s', GetCurFields() returned '%d'",pszInLociFile,NumFields);
				Rslt = eBSFerrFieldCnt;
				break;
				}

			Rslt = pCSV->GetInt(1,&SrcID);
			if(Rslt < eBSFSuccess)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process","Unable to read sequence identifier from field 1 at line %d from %s",pCSV->GetLineNumber(),pszInLociFile);
				break;
				}
			sprintf(szDescr,"Seq%d",SrcID);
			Rslt = pCSV->GetText(6,&pszSeq);
			if(Rslt < eBSFSuccess)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process","Unable to read sequence field 6 at line %d from %s",pCSV->GetLineNumber(),pszInLociFile);
				break;
				}
			Len = (int)strlen(pszSeq);
			if(Len < MinLength)
				{
				NumUnderLen += 1;
				continue;
				}
	
			if(Len >TruncLength)
				{
				Len = TruncLength;
				NumTruncated += 1;
				}
			break;

		case eCSVFSeq:
			Rslt = pCSV->GetText(DescrIDX,&pszDescr);
			if(Rslt < eBSFSuccess)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process","Unable to read descriptor field %d at line %d from %s",DescrIDX,pCSV->GetLineNumber(),pszInLociFile);
				break;
				}
			strncpy(szDescr,pszDescr,sizeof(szDescr));
			szDescr[sizeof(szDescr)-1] = '\0';
			Rslt = pCSV->GetText(SeqIDX,&pszSeq);
			if(Rslt < eBSFSuccess)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process","Unable to read sequence field %d at line %d from %s",SeqIDX,pCSV->GetLineNumber(),pszInLociFile);
				break;
				}
			Len = (int)strlen(pszSeq);
			if(Len < MinLength)
				{
				NumUnderLen += 1;
				continue;
				}
	
			if(Len > TruncLength)
				{
				NumTruncated += 1;
				TruncLength = Len;
				}
			break;
		}

	if(Rslt < eBSFSuccess)
		break;

	ReqAllocLen = Len;
	if(pElSeqBuff == NULL || ReqAllocLen > AllocSeqBuffLen)
		{
		if(pElSeqBuff != NULL)
			delete pElSeqBuff;
		ReqAllocLen += 10000;		  // a little extra to reduce number of potential subsequent reallocations
		if((pElSeqBuff = (etSeqBase *)new unsigned char[ReqAllocLen])==NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate %s bytes as a sequence buffer",ReqAllocLen);
			Rslt = eBSFerrMem;
			break;
			}
		AllocSeqBuffLen = ReqAllocLen;
		}

	ReqAllocLen = (Len * 3)/2;			// need to account for EOLs plus fasta descriptors etc
	if(pszLineBuff == NULL || ReqAllocLen > AllocLineBuffLen)	
		{
		if(pszLineBuff != NULL)
			delete pszLineBuff;
		ReqAllocLen += 10000;			// a little extra to reduce number of potential subsequent reallocations
		if((pszLineBuff = (char *)new char[ReqAllocLen])==NULL) 
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate %s bytes as a line buffer",ReqAllocLen);
			Rslt = eBSFerrMem;
			break;
			}
		AllocLineBuffLen = ReqAllocLen;
		}

	// now get the sequence
	if(CSVFormat != eCSVFSeq && CSVFormat != eCSVFQuickCount)
		{
		if((Rslt= ChromID = pBioseq->LocateEntryIDbyName(pszChrom))<=0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to locate chrom '%s' in assembly file '%s'",pszChrom,pszInSeqFile);
			break;
			}

		if((Rslt=pBioseq->GetData(ChromID,eSeqBaseType,StartLoci,pElSeqBuff,Len)) != Len)
			{
			if(NumLoadSeqFails++ < 25)
				gDiagnostics.DiagOut(eDLInfo,gszProcName,"Loading sequence failed from chrom: %s loci %d len: %d file: '%s'",pszChrom,StartLoci,Len,pszInSeqFile);
			continue;
			}
		}
	else
		{
		// copy in sequence removing any embedded whitespace, quotes and commas
		Len = 0;
		pElSeq = pElSeqBuff;
		while(Chr = *pszSeq++) {
			if(isspace(Chr))
				continue;
			switch(Chr) {
				case '"':
				case '\'':
				case ',':
					continue;
				case 'a':
					*pElSeq++ = eBaseA | cRptMskFlg;
					break;
				case 'A':
					*pElSeq++ = eBaseA;
					break;
				case 'c':
					*pElSeq++ = eBaseC  | cRptMskFlg;
					break;
				case 'C':
					*pElSeq++ =  eBaseC;
					break;
				case 'g':
					*pElSeq++ =  eBaseG  | cRptMskFlg;
					break;
				case 'G':
					*pElSeq++ =  eBaseG;
					break;
				case 't': case 'u': 
					*pElSeq++ =  eBaseT  | cRptMskFlg;
					break;
				case 'T': case 'U':
					*pElSeq++ =  eBaseT;
					break;
				case '-':
					*pElSeq++ =  eBaseInDel;
					break;
				default: // 'N'
					*pElSeq++ = eBaseN;
					break;
				}
			Len++;
			}
		// check on min/max length
		if(Len < MinLength)
			{
			NumUnderLen += 1;
			continue;
			}
	
		if(Len > TruncLength)
			{
			NumTruncated += 1;
			Len = TruncLength;
			}
		}

	if(bOnMinStrand) // if on minus strand then reverse complement..
		CSeqTrans::ReverseComplement((unsigned int)Len,pElSeqBuff);

	// processing mode specific
	switch(RsltsFormat) {
		case eRsltsFasta:
		case eRsltsRMESFasta:
			if(bFastaStart)
				{
				BuffLen = sprintf(pszLineBuff,">lcl|%s\n",RsltsFormat== eRsltsFasta ? "concatenated" : "RMESconcatenated" );
				bFastaStart = false;
				NxtFastaCol = 0;
				}
			else
				if(RsltsFormat == eRsltsRMESFasta)			// R'MES format is for 'z' to be used as a separator between concatenated sequences 
					BuffLen += sprintf(&pszLineBuff[BuffLen],"Z");
			SeqOfs = 0;
			while(Len)
				{
				NumCols = Len > 70 ? 70 : Len;
				if((NumCols + NxtFastaCol) > 70)
					NumCols = 70 - NxtFastaCol;
				CSeqTrans::MapSeq2Ascii(&pElSeqBuff[SeqOfs],NumCols,&pszLineBuff[BuffLen]);
				BuffLen += NumCols;
				NxtFastaCol += NumCols;
				if(NxtFastaCol >= 70)
					{
					BuffLen += sprintf(&pszLineBuff[BuffLen],"\n");
					NxtFastaCol = 0;
					}
				Len -= NumCols;
				SeqOfs += NumCols;
				}
			break;

		case eRsltsMultifasta:
			BuffLen = sprintf(pszLineBuff,">lcl|%s\n",szDescr);
			SeqOfs = 0;
			while(Len)
				{
				NumCols = Len > 70 ? 70 : Len;
				CSeqTrans::MapSeq2Ascii(&pElSeqBuff[SeqOfs],NumCols,&pszLineBuff[BuffLen]);
				BuffLen += NumCols;
				BuffLen += sprintf(&pszLineBuff[BuffLen],"\n");
				Len -= NumCols;
				SeqOfs += NumCols;
				}
			break;
		}
	if(BuffLen)
		{
		CUtility::SafeWrite(hRsltFile,pszLineBuff,BuffLen);
		BuffLen = 0;
		}
	NumElsAccepted += 1;
	}
close(hRsltFile);
delete pCSV;
if(pBioseq != NULL)
	delete pBioseq;
if(pElSeqBuff != NULL)
	delete pElSeqBuff;
if(pszLineBuff != NULL)
	delete pszLineBuff;
if(Rslt >= eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Elements processed: %d, Accepted: %d",NumElsProcessed,NumElsAccepted);
	gDiagnostics.DiagOut(eDLInfo,gszProcName," Truncated %d and filtered out %d because of Strand: %d, Chrom: %d, UnderLen: %d, UnderIdentity: %d, OverIdentity: %d, LoadSeqFails: %d",
					NumTruncated,NumElsProcessed - NumElsAccepted,NumStrandSkipped,NumExclChroms,NumUnderLen,NumUnderIdentity,NumOverIdentity,NumLoadSeqFails);
	}

return(Rslt > eBSFSuccess ? NumElsAccepted : Rslt);
}


// ProcessFastaFile
// Parse input fasta format file into a biosequence file
bool ProcessFastaFile(char *pszFile,
				 void *pParams)			
{
CFasta Fasta;
unsigned char *pSeqBuff;
char szName[cBSFSourceSize];
char szDescription[cBSFDescriptionSize];
int SeqLen;
int Descrlen;
bool bFirstEntry;
bool bEntryCreated;
tBSFEntryID EntryID;
int Rslt;
int SeqID;
bool bCapsSoftMask;

if((pSeqBuff = new unsigned char [cMaxAllocBuffChunk]) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessFastaFile:- Unable to allocate memory (%d bytes) for sequence buffer",cMaxAllocBuffChunk);
	return(false);
	}

CBioSeqFile *pSeqFile = ((tsProcParams *)pParams)->pSeqFile;
bCapsSoftMask = ((tsProcParams *)pParams)->bCapsSoftMask;

gDiagnostics.DiagOut(eDLInfo,gszProcName,"ProcessFastaFile:- Adding %s..",pszFile);
if((Rslt=Fasta.Open(pszFile,true))!=eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open '%s' [%s] %s",pszFile,Fasta.ErrText((teBSFrsltCodes)Rslt),Fasta.GetErrMsg());
	delete pSeqBuff;
	return(false);
	}
bFirstEntry = true;
bEntryCreated = false;
SeqID = 0;
while((Rslt = SeqLen = Fasta.ReadSequence(pSeqBuff,cMaxAllocBuffChunk,true,bCapsSoftMask)) > eBSFSuccess)
	{
	if(SeqLen == eBSFFastaDescr)		// just read a descriptor line
		{
		SeqID++;
		if(bEntryCreated)	// close any previous entry
			pSeqFile->SealEntry();
		Descrlen = Fasta.ReadDescriptor(szDescription,cBSFDescriptionSize);
		// An assumption - will one day bite real hard - is that the
		// fasta descriptor line starts with some form of unique identifier.
		// Use this identifier as the entry name.
		if(sscanf(szDescription," %s[ ,]",szName)!=1)
			sprintf(szName,"%s.%d",pszFile,++SeqID);
		EntryID = pSeqFile->CreateEntry(szName, szDescription,eSeqBaseType,bCapsSoftMask);
		if(EntryID < eBSFSuccess)
			{
			delete pSeqBuff;
			return(false);
			}
		bFirstEntry = false;
		bEntryCreated = true;
		continue;
		}
	else
		if(bFirstEntry)	// if there was no descriptor then dummy up one...
			{
			SeqID++;
			sprintf(szName,"%s.%d",pszFile,SeqID);
			EntryID =  pSeqFile->CreateEntry(szName, (char *)"No Description provided",eSeqBaseType,bCapsSoftMask);
			if(EntryID < eBSFSuccess)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unble to create entry '%s' [%s] %s",szName,Fasta.ErrText((teBSFrsltCodes)EntryID),Fasta.GetErrMsg());
				delete pSeqBuff;
				return(false);
				}
			bFirstEntry = false;
			bEntryCreated = true;
			}
	if((Rslt=pSeqFile->AddData(SeqLen,pSeqBuff)) < eBSFSuccess) // save sequence
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to add data [%s] %s",Fasta.ErrText((teBSFrsltCodes)Rslt),Fasta.GetErrMsg());
		delete pSeqBuff;
		return(false);
		}
	}
if(Rslt < eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessFastaFile [%s] %s",Fasta.ErrText((teBSFrsltCodes)Rslt),Fasta.GetErrMsg());
	delete pSeqBuff;
	return(false);
	}
if(bEntryCreated)
	pSeqFile->SealEntry();
delete pSeqBuff;
return(true);
}



int
CreateBioseqFastaFile(bool bTargDeps,			// true if process only if any independent src files newer than target
					  char *pszSrcDirPath,char *pszDestBioseqFile,char *pszRefSpecies,
						char *pszDescr,char *pszTitle,
						bool bCapsSoftMask)
{
int Rslt;
tsProcParams ProcParams;
bool bCreate;

CSimpleGlob glob(SG_GLOB_FULLSORT);
if (glob.Add(pszSrcDirPath) >= SG_SUCCESS)
	{
	bCreate = bTargDeps ? false : true;
	for (int n = 0; n < glob.FileCount(); ++n)
		{
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Will process in this order: %d '%s", n+1,glob.File(n));
		if(!bCreate && CUtility::ChkTargDepend(NULL,0,pszDestBioseqFile,glob.File(n),NULL)!= 0)
			bCreate = true;
		}
	}
else
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to glob '%s",pszSrcDirPath);
	return(eBSFerrOpnFile);	
    }
if(!bCreate)
	{
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Target '%s' is already updated",pszDestBioseqFile);
	return(eBSFSuccess);
	}

CBioSeqFile *pSeqFile = new CBioSeqFile;
if(pSeqFile == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"CreateBioseqFastaFile: Unable to instantiate CBioSeqFile");
	return(eBSFerrObj);
	}
if((Rslt=pSeqFile->Open(pszDestBioseqFile,cBSFTypeSeq,true)) != eBSFSuccess)
	{
	while(pSeqFile->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,pSeqFile->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"CreateBioseqFastaFile: Unable to create/truncate %s",pszDestBioseqFile);
	delete pSeqFile;
	return(Rslt);
	}
if((Rslt = pSeqFile->SetDescription(pszDescr)) != eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"CreateBioseqFastaFile: Unable to set description '%s' into %s",pszDescr,pszDestBioseqFile);
	delete pSeqFile;
	return(Rslt);
	}
if((Rslt=pSeqFile->SetTitle(pszTitle)) != eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"CreateBioseqFastaFile: Unable to set title '%s' into %s",pszTitle,pszDestBioseqFile);
	delete pSeqFile;
	return(Rslt);
	}

if((Rslt=pSeqFile->SetDatasetName(pszRefSpecies)) != eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"CreateBioseqFastaFile: Unable to set dataset %s",pszRefSpecies);
	delete pSeqFile;
	return(Rslt);
	}
ProcParams.bCapsSoftMask = bCapsSoftMask;
ProcParams.pSeqFile = pSeqFile;

srand(0);	// rand() used when processing for FMIndex generation to randomise long sequences of eBaseN otherwise ds_sort takes too long

for (int n = 0; Rslt >= eBSFSuccess &&  n < glob.FileCount(); ++n)
	Rslt=ProcessFastaFile(glob.File(n),&ProcParams);

pSeqFile->Close();
delete pSeqFile;
return(Rslt);
}



