// locseqs.cpp : Defines the entry point for the console application.
// Dynamically creates a suffix array from concatenated target sequences (with markers separating each target sequence)
// Then processes probes against this suffix array, locating unanchored exact matches of at least a specified minimum length and then maximally
// extends these matches left and right until a specified number of mismatches.
// Probes and targets can be from either multifasta (.fa) or bioseq (.seq) files.
// This processing is designed for situations wherein there are large numbers of small probes and targets (all vs. all for example), not
// for processing probes against chromosomes or a complete assembly - use locsfx in that instance!



#include "stdafx.h"
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#if _WIN32

#include "../libbiokanga/commhdrs.h"

#else
#include "../libbiokanga/commhdrs.h"
#endif

#include "SeqSfx.h"


CStopWatch gStopWatch;
CDiagnostics gDiagnostics;				// for writing diagnostics messages to log file
char gszProcName[_MAX_FNAME];			// process name

const unsigned int cProgVer = 301;		// increment with each release

// processing modes
typedef enum eProcMode {
	ePMHitLoci = 0,						// default processing for hit loci
	ePMHitFreq							// generate hit frequency summary
} etProcMode;

typedef enum eFileType {
	eFTFasta = 0,				// fasta
	eFTBioSeq					// bioseq
} etFileType;

typedef struct TAG_sSeqDescrs {
	int MaxSeqLen;				// max length of any sequence
	int MinSeqLen;				// minimum length of any sequence
	unsigned int TotSeqLen;		// total length of all sequences
	int NumDescrs;				// number of descriptors actually loaded
	int AllocDescrs;			// descriptors allocated
	sAMProbe *pDescrs;			// pts to array of allocated sequence descriptors

	int CurCacheLen;			// current total cached seq length
	int AllocCacheLen;			// allocated cache sequence length
	etSeqBase *pCacheSeq;		// to hold cached sequences

	int CurSeqLen;				// currently loaded sequence length
	int AllocSeqBases;			// currently allocd to hold bases ptd to by gpProbeBases;
	etSeqBase *pSeq;			// to hold current sequence(s)
} tsSeqDescrs;


const int cDfltProbeHits = 1000;		// report upto this many hits for any probe
const int cDfltRptMsk = eRPTHignor;		// default is to still process matches containing repeats
const int cDfltRptThres = 80;			// if repeat masking then default is to slough matches which contain more than 80% repeats
const int cAllocAMDescrs = 10000;		// allocate for sequence descriptors in this many instance increments
const int cMaxAllocCacheSeq = 200000000; // allocate at most this much for holding cached sequences

const int cMinCoreLen = 50;				// minimum match length as proportion of probe length if proportional core length specified
const int cMaxMatchProp = 100;			// maximum match length as proportion of probe length if proportional core length specified
const int cDfltMatchProp = cMaxMatchProp;	// default

const int cMinAbsLen  = 8;				// min core length if absolute core length specified
const int cDfltAbsLen = 25;				// default min core length if absolute core length specified

const int cMaxMismatches = 100;			// shouldn't really be a need for an upper limit on max mismatches but... 

etProcMode gProcMode = ePMHitLoci;	// processing mode
etFileType gProbeFileType = eFTFasta; // probe file type
void *gpProbeFile = NULL;			// will be either CFasta or CBioSeqFile when probe file opened
etFileType gTargFileType = eFTFasta; // target file type
void *gpTargFile = NULL;			// will be either CFasta or CBioSeqFile when target file opened


int gMaxProbeHits = cDfltProbeHits;	// maximum number of hits to for any probe
int gMatchCnt = 0;					// total number of matches

tsSeqDescrs gProbeDescrs;			// probe descriptors
tsSeqDescrs gTargDescrs;			// target descriptors

char *gpszProbeSpecies;				// probe species
char *gpszTargSpecies;				// target species
int ghRsltsFile;					// where to write results into
char *gpGenomeChrom;				// current genome chromosome being processed

char *gpszRsltsFile;				// results file name

int TrimQuotes(char *pszTxt);
char *RptMask2Txt(etRPTMasking RptMask,int RptThres);
int 
Process(etProcMode ProcMode,	// processing mode
		int MaxProbeHits,		// Maximum number of hits per chromosome per probe
		bool bAbsHitLen,		// if true then MinCoreLen and MinMatchLen are absolute lengths
		int MinCoreLen,			// otherwise min exact core length  a percentage of actual probe length
		int MinMatchLen,		// otherwise min match length is a percentage of actual probe length
	    int MaxMismatches,		// hits to have at most this many mismatches in the left or right flank around exactly matching core
		etFileType ProbeInType,	//input probe file type
		char *pszProbesFile,	// path+name of probes file
		etFileType TargInType,	//input target file type
		char *pszTargFile,		// path+name of target sequence
		char *pszRsltsFile,		// results into this file
		char *pszProbeSpecies,  // probe species
		char *pszTargSpecies,	// target species
		char Strand,			// which probe strand '*', '+' or '-'
	   etRPTMasking RptMask,	// repeat mask processing mode
	   int RptThres);			// filter out matches with repeats above this percentage threshold

int
ReportHitFreqs(bool bAbsHitLen,		// if true then pctMinCoreLen and pctMinMatchLen are absolute length
		   int pctMinCoreLen,		// otherwise core length is a percentage of actual probe length
		   int pctMinMatchLen,		// otherwise minimum match length is a percentage of actual probe length
	       int hRsltsFile);			// write hits to this file)
		

etSeqBase *LoadSeq(tsSeqDescrs *pDescrs,int DescrIdx);	// load probe sequence
		
int LoadFastaDescriptors(tsSeqDescrs *pSeqDescrs,   // where to load descriptors into
						 CFasta *pFasta,		// preopened fasta file containing probe sequences
			   char Strand);		// treat sequences as being on this strand

int LoadBioSeqDescriptors(tsSeqDescrs *pSeqDescrs,   // where to load descriptors into
						  CBioSeqFile *pBioSeqFile, // preopened bioseq file containing probe sequences
			   char Strand);			 // treat sequences as being on this strand	

void DeleteDescrs(tsSeqDescrs *pDescrs);

void DeleteProbeFile(void);

int											// < 0 to terminate, 0 to stop processing current probe, > 0 to continue
HandleApproxMatches(sAMProbe *pProbe,		// probe
					unsigned int ProbeOfs,		// probe offset at which probe hit starts 
					unsigned int MatchLen,		// number of bases hit
					unsigned int TargetEntryID,	// target suffix array entry 
					unsigned int TargOfs,		// target loci hit	
					etSeqBase *pSeq);			// target sequence hit

int											// < 0 to terminate, 0 to stop processing current probe, > 0 to continue
HandleExactMatches(sAMProbe *pProbe,		// probe
					unsigned int ProbeOfs,		// probe offset at which probe hit starts 
					unsigned int MatchLen,		// number of bases hit
					unsigned int TargetEntryID,	// target suffix array entry 
					unsigned int TargOfs,		// target loci hit	
					etSeqBase *pSeq);			// target sequence hit

int 
ExactMatch(bool bAbsHitLen,			// if true then MinCoreLen is an absolute length
		   int pctMinCoreLen,		// otherwise it's a percentage of actual probe length
		   etRPTMasking RptMask,	// repeat mask processing mode
		   int RptThres,			// filter out matches with repeats above this percentage threshold
		   int hRsltsFile,			// write hits to this file
		   CSeqSfx *pSfxArray,	// hits against sequences in this multi-suffix array
		   char Strand);				// '*':hits against both strands, '+' strand only, '-' strand only

int 
ApproxMatch(bool bAbsHitLen,		// if true then pctMinCoreLen and pctMinMatchLen are absolute length
		   int pctMinCoreLen,		// otherwise core length is a percentage of actual probe length
		   int pctMinMatchLen,		// otherwise minimum match length is a percentage of actual probe length
	       int MaxMismatches,		// hits to have at most this many mismatches in the left or right flank around exactly matching core
		   etRPTMasking RptMask,	// repeat mask processing mode
		   int RptThres,			// filter out matches with repeats above this percentage threshold
			int hRsltsFile,			// write hits to this file
			CSeqSfx *pSfxArray,	// hits against sequences in this multi-suffix array
			char Strand);			// '*':hits against both strands, '+' strand only, '-' strand only

CSeqSfx *LoadTargSfx(void);

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

char szProbeSpecies[80];		// probe species
char szTargSpecies[80];		// target species

int Rslt;

etProcMode eMode;			// processing mode
char cStrand;				// which strand '+' or '-' or '*'
int iProbeInType;			// input probe file type: 0=multifasta, 1=bioseq
int iTargInType;			// input target file type: 0=multifasta, 1=bioseq
int iMaxProbeHits;			// Maximum number of hits per chromosome per probe
bool bAbsHitLen;			// if true then MinCoreLen is an absolute length
int iMinCoreLen;			// otherwise it's a percentage of actual probe length
int iMinMatchLen;			// Minimum required match length as either absolute or proportion of probe length (10..100)
int iMaxMismatches;			// hits to have at most this many mismatches in the left or right flank around exactly matching core
int iRptMask;				// target repeat mask interpretation		
int iRptThres;				// filter out matches with more than this percentage (0-100) repeats

char szRsltsFile[_MAX_PATH];
char szTargFile[_MAX_PATH];
char szProbesFile[_MAX_PATH];

// command line args
struct arg_lit  *help    = arg_lit0("hH","help",                "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_int *ScreenLogLevel=arg_int0("S", "ScreenLogLevel",	"<int>","Level of diagnostics written to screen 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"diagnostics log file");

struct arg_int *ProcMode = arg_int0("m","mode","<int>",			"processing mode 0: output hit loci, 1: output hit frequency");

struct arg_file *ProbesFile = arg_file1("i","inprobes","<file>","input from probes multifasta or bioseq  file");
struct arg_int *ProbeInType=arg_int0("t", "probeintype","<int>","input probe file type: 0=multifasta, 1=bioseq");
struct arg_file *TargFile = arg_file1("I","inindex","<file>",	"input from target multifasta or bioseq file");
struct arg_int *TargInType=arg_int0("T", "targintype","<int>",	"input target file type: 0=multifasta, 1=bioseq");

struct arg_file *RsltsFile = arg_file1("o","outmatches","<file>","output matches to this file as CSV");

struct arg_int *MaxProbeHits=arg_int0("p", "MaxProbeHits","<int>","Maximum number of probe loci hits to report in mode 0 (1..1000000, default is 1000)");

struct arg_lit  *AbsHitLen = arg_lit0("a","abscorelen",				"if specified then MinCoreLen is an absolute length (10..N) and not a proportion");
struct arg_int *MinCoreLen=arg_int0("l", "MinCoreLen","<int>",		"Minimum required exact core length as either absolute or proportion of probe length (10..100)");
struct arg_int *MinMatchLen=arg_int0("L", "MinMatchLen","<int>",	"Minimum required match length as either absolute or proportion of probe length (10..100)");
struct arg_int *MaxMismatches=arg_int0("y", "MaxMismatches","<int>","Maximum number of mismatches allowed in exact core flanking regions");
struct arg_str *Strand=arg_str0("s", "strand",	"<string>",			"Strand '+' or '-', default is '*' for both");
struct arg_str *ProbeSpecies = arg_str1("r","ref","<string>",		"probe species");
struct arg_str *TargSpecies = arg_str1("R","rel","<string>",		"target species");
struct arg_int *RptMask=arg_int0("x", "rptmask","<int>",			"probe sequence repeat masking 0:ignor, 1: repeats, 2: invert repeat sense");
struct arg_int *RptThres=arg_int0("X", "rptthres","<int>",			"filter out exact matches with repeats which are greater than this percentage (0..100)");

struct arg_end *end = arg_end(20);

void *argtable[] = {help,version,FileLogLevel,ScreenLogLevel,LogFile,
					ProcMode,ProbeInType,ProbesFile,TargInType,TargFile,RsltsFile,MaxProbeHits,AbsHitLen,MinCoreLen,MinMatchLen,MaxMismatches,RptMask,RptThres,Strand,
					ProbeSpecies,TargSpecies,
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
		printf("\n      To invoke this parameter file then precede it's name with '@'");
		printf("\n      e.g. %s @myparams.txt\n\n",gszProcName);
		exit(1);
        }

    /* special case: '--version' takes precedence error reporting */
if (version->count > 0)
        {
		printf("\n%s Version %d.%2.2d",gszProcName,cProgVer/100,cProgVer%100);
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

	eMode = (etProcMode)(ProcMode->count ? ProcMode->ival[0] : ePMHitLoci);
	if(eMode < ePMHitLoci || eMode > ePMHitFreq)
		{
		printf("\nError: Requested processing mode '-m%d' not supported",eMode);
		exit(1);
		}

	if(eMode == ePMHitLoci)
		{
		iMaxProbeHits = MaxProbeHits->count ? MaxProbeHits->ival[0] : cDfltProbeHits;			
		if(iMaxProbeHits < 1 || iMaxProbeHits > 1000000)
			{
			printf("\nError: Requested max hits per probe '-p%d' is invalid, must be 1..1000000",iMaxProbeHits);
			exit(1);
			}
		}
	else
		iMaxProbeHits = 0;

	iMaxMismatches = MaxMismatches->count ? MaxMismatches->ival[0] : 0;			// no mismatches allowed is the default
	if(iMaxMismatches < 0 || iMaxMismatches > cMaxMismatches)
		{
		printf("\nError: Requested max mismatches '-y%d' is invalid, must be %d..%d",iMaxMismatches,cMaxMismatches);
		exit(1);
		}

	bAbsHitLen = AbsHitLen->count ? true : false;
	if(!bAbsHitLen)
		{
		iMinCoreLen = MinCoreLen->count ? MinCoreLen->ival[0] : cDfltMatchProp;
		if(iMinCoreLen < cMinCoreLen || iMinCoreLen > cMaxMatchProp)
			{
			printf("\nError: Requested hit length proportion '-l%d' is invalid, must be %d..%d",iMinCoreLen,cMinCoreLen,cMaxMatchProp);
			exit(1);
			}

		if(!iMaxMismatches)
			iMinMatchLen = iMinCoreLen;
		else
			{
			iMinMatchLen = MinMatchLen->count ? MinMatchLen->ival[0] : iMinCoreLen;
			if(iMinMatchLen < iMinCoreLen || iMinMatchLen > cMaxMatchProp)
				{
				printf("\nError: Requested min match length proportion '-L%d' is invalid, must be %d..%d",iMinMatchLen,iMinCoreLen,cMaxMatchProp);
				exit(1);
				}
			}
		}
	else
		{
		iMinCoreLen = MinCoreLen->count ? MinCoreLen->ival[0] : cDfltAbsLen;
		if(iMinCoreLen < cMinAbsLen)
			{
			printf("\nError: Requested absolute exact core length '-n%d' is invalid, must be at least %d",iMinCoreLen,cMinAbsLen);
			exit(1);
			}
		if(!iMaxMismatches)
			iMinMatchLen = iMinCoreLen;
		else
			{
			iMinMatchLen = MinMatchLen->count ? MinMatchLen->ival[0] : iMinCoreLen;
			if(iMinMatchLen < iMinCoreLen)
				{
				printf("\nError: Requested absolute min match length '-L%d' is invalid, must be at least %d",iMinMatchLen,iMinCoreLen);
				exit(1);
				}
			}
		}

	if(RptMask->count)
		{
		iRptMask = RptMask->ival[0];
		if(iRptMask < eRPTHignor || iRptMask > eRPTHunmasked)
			{
			printf("\nError: Requested repeat mask processing '-X%d' is invalid, must be 0..2",iRptMask);
			exit(1);
			}
		}
	else
		iRptMask = cDfltRptMsk;



	if(iRptMask != cDfltRptMsk)
		{
		if(RptThres->count)
			{
			iRptThres = RptThres->ival[0];
			if(iRptThres < 0 || iRptThres > 100)
				{
				printf("\nError: Requested repeat threshold '-x%d' is invalid, must be 0..100",iRptThres);
				exit(1);
				}
			}
		else
			iRptThres = cDfltRptThres;
		}
	else
		iRptThres = 0;


	iProbeInType = ProbeInType->count ? ProbeInType->ival[0] : 0;
	if(iProbeInType < 0 || iProbeInType > 1)
		{
		printf("\nError: Requested probe file format '-t%d' is invalid, must be 0 or 1",iProbeInType);
		exit(1);
		}

	iTargInType = TargInType->count ? TargInType->ival[0] : 0;
	if(iTargInType < 0 || iTargInType > 1)
		{
		printf("\nError: Requested target file format '-T%d' is invalid, must be 0 or 1",iTargInType);
		exit(1);
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

	if(TargFile->count)
		strcpy(szTargFile,TargFile->filename[0]);
	else
		strcpy(szTargFile,"in.fmi");


	if(RsltsFile->count)
		strcpy(szRsltsFile,RsltsFile->filename[0]);
	else
		strcpy(szRsltsFile,"out.csv");

	if(ProbesFile->count)
		strcpy(szProbesFile,ProbesFile->filename[0]);
	else
		strcpy(szProbesFile,"probes.fa");

	if(ProbeSpecies->count)
		{
		strcpy(szProbeSpecies,ProbeSpecies->sval[0]);
		TrimQuotes(szProbeSpecies);
		}
	else
		strcpy(szProbeSpecies,"Probe");

	if(TargSpecies->count)
		{
		strcpy(szTargSpecies,TargSpecies->sval[0]);
		TrimQuotes(szTargSpecies);
		}
	else
		strcpy(szTargSpecies,"Target");

		// now that command parameters have been parsed then initialise diagnostics log system
	if(!gDiagnostics.Open(szLogFile,(etDiagLevel)iScreenLogLevel,(etDiagLevel)iFileLogLevel,true))
		{
		printf("\nError: Unable to start diagnostics subsystem.");
		if(szLogFile[0] != '\0')
			printf(" Most likely cause is that logfile '%s' can't be opened/created",szLogFile);
		exit(1);
		}

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Version: %d.%2.2d Processing parameters:",cProgVer/100,cProgVer%100);

	switch(eMode) {
		case ePMHitLoci:
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Processing mode: report hit loci for each probe");
			break;

		case ePMHitFreq:
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Processing mode: report hit frequency for each probe");
			break;
		}

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Probe sequences from %s file: %s",iProbeInType > 0 ? "bioseq" : "multifasta",szProbesFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Probe species: '%s'",szProbeSpecies);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Target sequences from %s file: %s",iTargInType > 0 ? "bioseq" : "multifasta",szTargFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Target species: '%s'",szTargSpecies);
	
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Output as CVS file: '%s'",szRsltsFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"All hit loci will be generated");
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Target strand: '%c'",cStrand);
	if(eMode == ePMHitLoci)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Max hits per probe to report: %d",iMaxProbeHits);
	if(!bAbsHitLen)
		{
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Minimum exact match core length as proportion of probe length: %d%%",iMinCoreLen);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Minimum matching length as proportion of probe length: %d%%",iMinMatchLen);
		}
	else
		{
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Minimum exact match core length as absolute length: %d",iMinCoreLen);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Minimum matching length as absolute length: %d",iMinMatchLen);
		}

	gDiagnostics.DiagOutMsgOnly(eDLInfo,RptMask2Txt((etRPTMasking)iRptMask,iRptThres));
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Maximum number of mismatches allowed in exact core flanks: %d",iMaxMismatches);
#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	gStopWatch.Start();
	Rslt = Process(eMode,iMaxProbeHits,bAbsHitLen,iMinCoreLen,iMinMatchLen,iMaxMismatches,(etFileType)iProbeInType,szProbesFile,(etFileType)iTargInType,szTargFile,szRsltsFile,szProbeSpecies,szTargSpecies,cStrand,(etRPTMasking)iRptMask,iRptThres);
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

int
TrimWhitespace(char *pTxt)
{
int Len = 0;
char *pStart = pTxt;
char Chr;
if(pTxt == NULL || *pTxt == '\0')
	return(0);

	// strip any leading whitespace
while((Chr = *pTxt) && isspace(Chr))
	pTxt++;
if(Chr == '\0')					// empty line?
	{
	*pStart = '\0';
	return(0);
	}

while(*pTxt)			// fast forward to string terminator '\0' 
	{
	*pStart++ = *pTxt++;
	Len++;
	}
pStart--;				// backup to last chr 
while(isspace(*pStart))
	{
	pStart--;
	Len--;
	}
pStart[1] = '\0';
return(Len);
}

// TrimQuotes
// Removes any whitespace, leading and trailing quotes from pszTxt
int
TrimQuotes(char *pszTxt)
{
char *pDst;
char *pChr;
char Chr;
int Len = 0;
TrimWhitespace(pszTxt);
pChr = pszTxt;
pDst = pszTxt;
while((Chr = *pChr++))	
	{
	if((!Len || *pChr == '\0') && (Chr == '"' || Chr == '\''))
		continue;
	*pDst++ = Chr;
	Len++;
	}
*pDst = Chr;
return(TrimWhitespace(pszTxt));
}


char *
RptMask2Txt(etRPTMasking RptMask,int RptThres)
{
static char szTxt[128];
switch(RptMask) {
	case eRPTHignor:	// treat all bases as not being a repeat (ignore any cRptMskFlg)
		return((char *)"No repeats processing");
	case eRPTHmasked:	// treat any base cRptMskFlg as a repeat masked base
		sprintf(szTxt,"Probe matches with repeats above %d%%%% threshold will be sloughed",RptThres);
		return(szTxt);
	case eRPTHunmasked:	// treat any base without cRptMskFlg as a repeat masked base
		sprintf(szTxt,"Probe matches with non-repeats above %d%%%% threshold will be sloughed",RptThres);
		return(szTxt);
	default:
		break;
	}
return((char *)"Unsupported masking mode");
}

int
Process(etProcMode ProcMode,	// processing mode
		int MaxProbeHits,		// Maximum number of hits per chromosome per probe
		bool bAbsHitLen,		// if true then MinCoreLen and MinMatchLen are absolute lengths
		int MinCoreLen,			// otherwise min exact core length  a percentage of actual probe length
		int MinMatchLen,		// otherwise min match length is a percentage of actual probe length
	    int MaxMismatches,		// hits to have at most this many mismatches in the left or right flank around exactly matching core
		etFileType ProbeFileType,	// input probe file type
		char *pszProbesFile,	// path+name of probes file
		etFileType TargFileType,	// input target file type
		char *pszTargFile,		// path+name of target sequence
		char *pszRsltsFile,		// results into this file
		char *pszProbeSpecies,  // probe species
		char *pszTargSpecies,	// target species
		char Strand,			// which probe strand '*', '+' or '-'
	   etRPTMasking RptMask,	// repeat mask processing mode
	   int RptThres)
{
CSeqSfx *pTargSfx = NULL;
CFasta *pFasta = NULL;
CBioSeqFile *pBioSeqFile = NULL;

char szHeader[200];					// to hold output CSV header row text
int LineLen;
int Rslt=0;

memset(&gProbeDescrs,0,sizeof(gProbeDescrs));
memset(&gTargDescrs,0,sizeof(gTargDescrs));

gProcMode = ProcMode;
gProbeFileType = ProbeFileType;
gpProbeFile = NULL;
gTargFileType = TargFileType;
gpTargFile = NULL;

switch(ProbeFileType) {
	case eFTFasta:						// input probe(s) from multifasta file
		// open fasta file containing probes
		if((pFasta = new CFasta()) == NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CFasta");
			return(eBSFerrObj);
			}

		if((Rslt=pFasta->Open(pszProbesFile))!=eBSFSuccess)
			{
			while(pFasta->NumErrMsgs())
				gDiagnostics.DiagOut(eDLFatal,gszProcName,pFasta->GetErrMsg());
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open input fasta file containing probes '%s'",pszProbesFile);
			delete pFasta;
			return(Rslt);
			}
		
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Loading probes from fasta file '%s'",pszProbesFile);
		if((Rslt = LoadFastaDescriptors(&gProbeDescrs,pFasta,'+')) <= eBSFSuccess)	// expect at least one probe to be loaded!
			{
			delete pFasta;
			DeleteDescrs(&gProbeDescrs);
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"No probes loaded");
			return(Rslt);
			}
		gpProbeFile = pFasta;
		break;

	case eFTBioSeq:							    // input probes from bioseq file
		if((pBioSeqFile = new CBioSeqFile()) == NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CBioSeqFile");
			return(eBSFerrObj);
			}

		if((Rslt=pBioSeqFile->Open(pszProbesFile))!=eBSFSuccess)
			{
			while(pBioSeqFile->NumErrMsgs())
					gDiagnostics.DiagOut(eDLFatal,gszProcName,pBioSeqFile->GetErrMsg());
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open input probe bioseq file '%s'",pszProbesFile);
			delete pBioSeqFile;
			return(Rslt);
			}

		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Loading probes from bioseq file '%s'",pszProbesFile);
		if((Rslt = LoadBioSeqDescriptors(&gProbeDescrs,pBioSeqFile,'+')) <= eBSFSuccess)	// expect at least one probe to be loaded!
			{
			delete pBioSeqFile;
			DeleteDescrs(&gProbeDescrs);
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"No probes loaded");
			return(Rslt);
			}
		gpProbeFile = pBioSeqFile;
		break;
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Total Probes: %d",gProbeDescrs.NumDescrs);

switch(TargFileType) {
	case eFTFasta:						// input target(s) from multifasta file
		// open fasta file containing probes
		if((pFasta = new CFasta()) == NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CFasta");
			DeleteDescrs(&gProbeDescrs);
			return(eBSFerrObj);
			}

		if((Rslt=pFasta->Open(pszTargFile))!=eBSFSuccess)
			{
			while(pFasta->NumErrMsgs())
				gDiagnostics.DiagOut(eDLFatal,gszProcName,pFasta->GetErrMsg());
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open fasta file containing targets '%s'",pszTargFile);
			delete pFasta;
			DeleteDescrs(&gProbeDescrs);
			return(Rslt);
			}
		
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Loading targets from fasta file '%s'",pszTargFile);
		if((Rslt = LoadFastaDescriptors(&gTargDescrs,pFasta,'+')) <= eBSFSuccess)	// expect at least one sequence descr to be loaded!
			{
			delete pFasta;
			DeleteDescrs(&gProbeDescrs);
			DeleteDescrs(&gTargDescrs);
			return(Rslt);
			}
		gpTargFile = pFasta;
		break;

	case eFTBioSeq:							    // input targs from bioseq file
		if((pBioSeqFile = new CBioSeqFile()) == NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CBioSeqFile");
			DeleteDescrs(&gProbeDescrs);
			DeleteDescrs(&gTargDescrs);
			return(eBSFerrObj);
			}

		if((Rslt=pBioSeqFile->Open(pszTargFile))!=eBSFSuccess)
			{
			while(pBioSeqFile->NumErrMsgs())
					gDiagnostics.DiagOut(eDLFatal,gszProcName,pBioSeqFile->GetErrMsg());
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open targets bioseq file '%s'",pszTargFile);
			delete pBioSeqFile;
			DeleteDescrs(&gProbeDescrs);
			DeleteDescrs(&gTargDescrs);
			return(Rslt);
			}

		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Loading targets from bioseq file '%s'",pszTargFile);
		if((Rslt = LoadBioSeqDescriptors(&gTargDescrs,pBioSeqFile,'+')) <= eBSFSuccess)	// expect at least one seq descr to be loaded!
			{
			delete pBioSeqFile;
			DeleteDescrs(&gProbeDescrs);
			DeleteDescrs(&gTargDescrs);
			return(Rslt);
			}
		gpTargFile = pBioSeqFile;
		break;
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Total targets: %d",gTargDescrs.NumDescrs);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Generating multiple target suffix array on concatenated targets");
if((pTargSfx = LoadTargSfx())==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create sequence suffix over target sequences in '%s'",pszTargFile);
	DeleteDescrs(&gProbeDescrs);
	DeleteDescrs(&gTargDescrs);
	return(-1);
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Generated multiple target suffix array");

// open and truncate CSV output file
#ifdef _WIN32
ghRsltsFile = open(pszRsltsFile, _O_WRONLY | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE );
#else
ghRsltsFile = open(pszRsltsFile, O_WRONLY | O_CREAT | O_TRUNC, S_IREAD | S_IWRITE );
#endif
if(ghRsltsFile == -1)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open output CSV file '%s', error: %s",pszRsltsFile,strerror(errno));
	DeleteDescrs(&gProbeDescrs);
	DeleteDescrs(&gTargDescrs);
	delete pTargSfx;
	return(eBSFerrOpnFile);
	}

// output header
if(ProcMode == ePMHitFreq)				// if frequencies to be reported
	LineLen = sprintf(szHeader,"\"RefID\",\"SrcRefID\",\"ProbeRef\",\"MinHits\",\"PlusHits\",\"TotHits\"\n");
else
	LineLen = sprintf(szHeader,"\"RefID\",\"Type\",\"TargSpecies\",\"Chrom\",\"StartLoci\",\"EndLoci\",\"Len\",\"ProbeSpecies\",\"Region\",\"ProbeStrand\",\"ProbeID\",\"ProbeDescr\",\"ProbeOfs\",\"ProbeLen\"\n");

if(write(ghRsltsFile,szHeader,LineLen)!=LineLen)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Fatal: write of header len %d to '%s' failed, error was - ",LineLen,pszRsltsFile,strerror(errno));
	DeleteDescrs(&gProbeDescrs);
	DeleteDescrs(&gTargDescrs);
	delete pTargSfx;
	close(ghRsltsFile);
	return(eBSFerrFileAccess);
	}
_commit(ghRsltsFile);

gpszRsltsFile = pszRsltsFile;
gMaxProbeHits = MaxProbeHits;

gMatchCnt = 0;
Rslt = eBSFSuccess;
gpszProbeSpecies = pszProbeSpecies;
gpszTargSpecies = pszTargSpecies;

if((gProbeDescrs.pSeq = new etSeqBase[gProbeDescrs.MaxSeqLen + 1])==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate memory for maximal sized probe of %d",gProbeDescrs.MaxSeqLen);
	DeleteDescrs(&gProbeDescrs);
	DeleteDescrs(&gTargDescrs);
	DeleteProbeFile();
	return(eBSFerrMem);
	}
gProbeDescrs.AllocSeqBases = gProbeDescrs.MaxSeqLen+1;	
gProbeDescrs.CurSeqLen = 0;

if(gProbeDescrs.MinSeqLen < cMaxAllocCacheSeq)
	{
	gProbeDescrs.AllocCacheLen = min(gProbeDescrs.TotSeqLen,cMaxAllocCacheSeq) + 1;
	if((gProbeDescrs.pCacheSeq = new etSeqBase[gProbeDescrs.AllocCacheLen])==NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate %d bytes memory for caching probe sequnences ",gProbeDescrs.AllocCacheLen);
		DeleteDescrs(&gProbeDescrs);
		DeleteDescrs(&gTargDescrs);
		DeleteProbeFile();
		return(eBSFerrMem);
		}
	gProbeDescrs.CurCacheLen = 0;
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Now processing probes against targets");

if(MaxMismatches == 0)			// if hits must be exact with no mismatches allowed
		Rslt = ExactMatch(bAbsHitLen,MinCoreLen,RptMask,RptThres,ghRsltsFile,pTargSfx,Strand);
	else								 // mismatches are allowed
		Rslt = ApproxMatch(bAbsHitLen,
			MinCoreLen, 
			MinMatchLen,
			MaxMismatches,
			RptMask,
			RptThres,
			ghRsltsFile,				// write hits to this file
			pTargSfx,					// hits against sequences in this multi-suffix array
			Strand);					// '*':hits against both strands, '+' strand only, '-' strand only
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing probes against targets completed");

if(gProcMode == ePMHitFreq)				// if frequencies to be reported
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing hit frequencies");
	Rslt = ReportHitFreqs(bAbsHitLen,		// if true then pctMinCoreLen and pctMinMatchLen are absolute length
						MinCoreLen,		// otherwise core length is a percentage of actual probe length
						MinMatchLen,		// otherwise minimum match length is a percentage of actual probe length
						ghRsltsFile);		// write hits to this file)
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing hit frequencies completed");
	}

DeleteDescrs(&gProbeDescrs);
DeleteDescrs(&gTargDescrs);

if(ghRsltsFile!=-1)
	{
	_commit(ghRsltsFile);
	close(ghRsltsFile);
	ghRsltsFile = -1;
	}
return(Rslt);
}

int
ReportHitFreqs(bool bAbsHitLen,		// if true then pctMinCoreLen and pctMinMatchLen are absolute length
		   int pctMinCoreLen,		// otherwise core length is a percentage of actual probe length
		   int pctMinMatchLen,		// otherwise minimum match length is a percentage of actual probe length
	       int hRsltsFile)			// write hits to this file)
{
char szBuff[4096];
int BuffLen;
int MinCoreLen;
int CoreID;
int ProbeIdx;
int WrtLen;
sAMProbe *pCurDescr;
pCurDescr = gProbeDescrs.pDescrs;
BuffLen = 0;
CoreID = 0;
MinCoreLen = pctMinCoreLen;
for(ProbeIdx = 0; ProbeIdx < gProbeDescrs.NumDescrs; ProbeIdx++,pCurDescr++)
	{
	if(bAbsHitLen) 
		{
		if(pCurDescr->ProbeLen < pctMinCoreLen)
			continue;
		if(pCurDescr->ProbeLen < pctMinMatchLen)
			continue;
		}
	else
		{
		if(pctMinCoreLen == 100)
			MinCoreLen = pCurDescr->ProbeLen;
		else
			MinCoreLen = (pCurDescr->ProbeLen * pctMinCoreLen) / 100;
		}
	if(MinCoreLen < cMinAbsLen)		// ensure min core length is always >= cMinAbsLen
		MinCoreLen = cMinAbsLen;
	if(MinCoreLen > pCurDescr->ProbeLen)
		continue;

	BuffLen += sprintf(&szBuff[BuffLen],"%d,%d,\"%s\",%d,%d,%d\n",
			++CoreID,pCurDescr->ProbID,pCurDescr->szDescr,pCurDescr->TotMinHits,pCurDescr->TotPlusHits,pCurDescr->TotMinHits + pCurDescr->TotPlusHits);
	if(BuffLen > (sizeof(szBuff) - 500))
		{
		if((WrtLen = write(hRsltsFile,szBuff,BuffLen)) != BuffLen)
			{
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"Fatal: write of len %d to '%s' failed, error was - ",BuffLen,gpszRsltsFile,strerror(errno));
			return(eBSFerrFileAccess);
			}
		_commit(hRsltsFile);
		BuffLen = 0;
		}
	}
if(BuffLen)
	{
	if((WrtLen = write(hRsltsFile,szBuff,BuffLen)) != BuffLen)
		{
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Fatal: write of len %d to '%s' failed, error was - ",BuffLen,gpszRsltsFile,strerror(errno));
		return(eBSFerrFileAccess);
		}
	_commit(hRsltsFile);
	}
return(CoreID);
}


void 
DeleteProbeFile(void)
{
if(gpProbeFile != NULL)
	{
	switch(gProbeFileType) {
		case eFTFasta:
			delete (CFasta *)gpProbeFile;
			break;
		case eFTBioSeq:
			delete (CBioSeqFile *)gpProbeFile;
			break;
		}
	gpProbeFile = NULL;
	}
}

void
DeleteDescrs(tsSeqDescrs *pDescrs)
{
if(pDescrs->pDescrs != NULL)
	delete pDescrs->pDescrs;

if(pDescrs->pSeq != NULL)
	delete pDescrs->pSeq;
memset(pDescrs,0,sizeof(*pDescrs));
}


int								// < 0 to terminate, 0 to stop processing current probe, > 0 to continue
FilterApproxMatches(sAMProbe *pProbe,	// probe
			unsigned int ProbeOfs,		// probe offset at which probe hit starts 
			unsigned int HitLen,		// number of bases hit
			unsigned int ProbeCoreOfs,	// probe offset at which exact core starts
			unsigned int ProbeCoreLen,	// length of exact core
			unsigned int NumLeftMismatches, // number of mismatches in left flank
			unsigned int NumRightMismatches, // number of mismatches in right flank
			unsigned int TargetRef,			// user assigned target which was to be associated with any hit to this sequence
			unsigned int TargOfs,		// target loci hit	
			etSeqBase *pSeq)			// target sequence, hit starts at &pSeq[TargOfs]
{
static int ProbeID = -1;
static int TargID = -1;
static int NumHits = -1;
if(gProcMode == ePMHitLoci)
	{
	if(NumHits < 0 || ProbeID != pProbe->ProbID || TargID != TargetRef) 
		{
		ProbeID = pProbe->ProbID;
		TargID = TargetRef;
		NumHits = 1;
		return(1);
		}
	if(NumHits++ > (gMaxProbeHits * 1000))	// this is simply to ensure that if a low limit was set then we don't bother continueing with probe when it will be filtered out later
		return(0);
	}
return(1);			// currently no real filtering is implemented
}

int								// < 0 to terminate, 0 to stop processing current probe, > 0 to continue 
FilterExactMatches(sAMProbe *pProbe,	// probe
			unsigned int ProbeOfs,		// probe offset at which exact core starts 
			unsigned int ProbeCoreLen,	// length of exact core
			unsigned int TargetRef,		// user assigned target which was to be associated with any hit to this sequence
			unsigned int TargOfs,		// target loci hit	
			etSeqBase *pSeq)			// target sequence, hit starts at &pSeq[TargOfs]
{
static int ProbeID = -1;
static int TargID = -1;
static int NumHits = -1;
if(gProcMode == ePMHitLoci)
	{
	if(NumHits < 0 || ProbeID != pProbe->ProbID || TargID != TargetRef) 
		{
		ProbeID = pProbe->ProbID;
		TargID = TargetRef;
		NumHits = 1;
		return(1);
		}
	if(NumHits++ > (gMaxProbeHits * 1000))	// this is simply to ensure that if a low limit was set then we don't bother continueing with probe when it will be filtered out later
		return(0);
	}
return(1);			// currently no real filtering is implemented
}


int											// < 0 to terminate, 0 to stop processing current probe, > 0 to continue
HandleApproxMatches(sAMProbe *pProbe,		// probe
					unsigned int ProbeOfs,		// probe offset at which probe hit starts 
					unsigned int MatchLen,		// number of bases hit
					unsigned int TargetRef,		// user assigned target which was to be associated with any hit to this sequence
					unsigned int TargOfs,		// target loci hit	
					etSeqBase *pSeq)			// target sequence hit
{
static int CoreID = 0;
char szLineBuff[1024];
int Len;
int WrtLen;
sAMProbe *pTargDescr;
pTargDescr = &gTargDescrs.pDescrs[TargetRef-1];

if(pProbe->CurStrand == '-')
	{
	ProbeOfs = pProbe->ProbeLen - (ProbeOfs + MatchLen); 
	pProbe->TotMinHits += 1;
	}
else
	pProbe->TotPlusHits += 1;

if(gProcMode == ePMHitLoci)
	{
	if((pProbe->TotMinHits + pProbe->TotPlusHits) > gMaxProbeHits)
		return(0); 

	Len = sprintf(szLineBuff,"%d,\"%s\",\"%s\",\"%s\",%d,%d,%d,\"%s\",0,\"%c\",%d,\"%s\",%d,%d\n",
			++CoreID,"exactcore",gpszTargSpecies,pTargDescr->szDescr,
			TargOfs,TargOfs+MatchLen-1,MatchLen,gpszProbeSpecies,pProbe->CurStrand,
			pProbe->ProbID,pProbe->szDescr,ProbeOfs,pProbe->ProbeLen);

	if((WrtLen = write(ghRsltsFile,szLineBuff,Len)) != Len)
		{
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Fatal: write of len %d to '%s' failed, error was - ",Len,gpszRsltsFile,strerror(errno));
		return(eBSFerrFileAccess);
		}
	_commit(ghRsltsFile);
	}
return(1);
}


int											// < 0 to terminate, 0 to stop processing current probe, > 0 to continue
HandleExactMatches(sAMProbe *pProbe,		// probe
					unsigned int ProbeOfs,		// probe offset at which probe hit starts 
					unsigned int MatchLen,		// number of bases hit
					unsigned int TargetRef,			// user assigned target which was to be associated with any hit to this sequence
					unsigned int TargOfs,		// target loci hit	
					etSeqBase *pSeq)			// target sequence hit
{
static int CoreID = 0;
char szLineBuff[1024];
int Len;
int WrtLen;
sAMProbe *pTargDescr;

if(pProbe->CurStrand == '-')
	{
	ProbeOfs = pProbe->ProbeLen - (ProbeOfs + MatchLen); 
	pProbe->TotMinHits += 1;
	}
else
	pProbe->TotPlusHits += 1;

pTargDescr = &gTargDescrs.pDescrs[TargetRef-1];
if(gProcMode == ePMHitLoci)
	{
	if((pProbe->TotMinHits + pProbe->TotPlusHits) > gMaxProbeHits)
		return(0);

	Len = sprintf(szLineBuff,"%d,\"%s\",\"%s\",\"%s\",%d,%d,%d,\"%s\",0,\"%c\",%d,\"%s\",%d,%d\n",
			++CoreID,"exactcore",gpszTargSpecies,pTargDescr->szDescr,
			TargOfs,TargOfs+MatchLen-1,MatchLen,gpszProbeSpecies,pProbe->CurStrand,
			pProbe->ProbID,pProbe->szDescr,ProbeOfs,pProbe->ProbeLen);

	if((WrtLen = write(ghRsltsFile,szLineBuff,Len)) != Len)
		{
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Fatal: write of len %d to '%s' failed, error was - ",Len,gpszRsltsFile,strerror(errno));
		return(eBSFerrFileAccess);
		}
	}
return(1);
}

// ExactMatch
// Locate all exact matches against probes in gpProbes
int 
ExactMatch(bool bAbsHitLen,			// if true then MinCoreLen is an absolute length
		   int pctMinCoreLen,		// otherwise it's a percentage of actual probe length
		   etRPTMasking RptMask,	// repeat mask processing mode
		   int RptThres,			// filter out matches with repeats above this percentage threshold
		   int hRsltsFile,			// write hits to this file
		   CSeqSfx *pSfxArray,		// hits against sequences in this multi-suffix array
		   char Strand)				// '*':hits against both strands, '+' strand only, '-' strand only
{
sAMProbe *pCurDescr;
etSeqBase *pProbeBases;
CStopWatch CurTime;
unsigned long Elapsed;
unsigned long PrvTime;

int ProbeIdx;
int MinHitLen;
int Rslt=0;

pCurDescr = gProbeDescrs.pDescrs;
CurTime.Start();
PrvTime = 0;
for(ProbeIdx = 0; ProbeIdx < gProbeDescrs.NumDescrs; ProbeIdx++, pCurDescr++)
	{
	Elapsed = CurTime.ReadUSecs();
	if(!(pCurDescr->ProbID % 1000) || pCurDescr->ProbID == 1 || Elapsed > PrvTime)
		{
		PrvTime = Elapsed;
		// NOTE: tabs in following printf are simply a poor mans clear to end of line hack!
		printf("\r\t\t\t\t\t\t\t\rProbe %d '%s'",pCurDescr->ProbID,pCurDescr->szDescr);
		}

	// no point in continuing with this probe if already had more hits than max requested
	if(gProcMode == ePMHitLoci && ((pCurDescr->TotMinHits + pCurDescr->TotPlusHits) > gMaxProbeHits))
		continue;

	if(bAbsHitLen) 
		{
		if(pCurDescr->ProbeLen >= pctMinCoreLen)
			MinHitLen = pctMinCoreLen;
		else
			continue;
		}
	else
		if(pctMinCoreLen == 100)
			MinHitLen = pCurDescr->ProbeLen;
		else
			MinHitLen = (pCurDescr->ProbeLen * pctMinCoreLen) / 100;

	if(MinHitLen < cMinAbsLen)		// ensure min core length is always >= cMinAbsLen
		MinHitLen = cMinAbsLen;
	if(MinHitLen > pCurDescr->ProbeLen)
		continue;
		
	if((pProbeBases = LoadSeq(&gProbeDescrs,ProbeIdx))==NULL)
		return(0);

	if(Strand == '*' || Strand == '+')
		{
		if(pCurDescr->CurStrand != '+')
			{
			pCurDescr->CurStrand = '+';
			CSeqTrans::ReverseComplement(pCurDescr->ProbeLen,pProbeBases);
			}

		if((Rslt = pSfxArray->LocateExacts(pCurDescr,pCurDescr->ProbeLen,pProbeBases,MinHitLen,FilterExactMatches,HandleExactMatches,RptMask,RptThres))>=0)
			gMatchCnt += Rslt;
		else
			return(Rslt);
	}


	if(Strand == '*' || Strand == '-')
		{
		if(pCurDescr->CurStrand != '-')
			{
			pCurDescr->CurStrand = '-';
			CSeqTrans::ReverseComplement(pCurDescr->ProbeLen,pProbeBases);
			}
		// no point in continuing with this probe if already had more hits than max requested
		if(gProcMode == ePMHitLoci && ((pCurDescr->TotMinHits + pCurDescr->TotPlusHits) > gMaxProbeHits))
			continue;
		Rslt = pSfxArray->LocateExacts(pCurDescr,pCurDescr->ProbeLen,pProbeBases,MinHitLen,FilterExactMatches,HandleExactMatches,RptMask,RptThres);
		if(Rslt >= 0)
			gMatchCnt += Rslt;
		else
			return(Rslt);
		}
	}

return(gMatchCnt);
}

// ApproxMatch
// Locate all matches against probes in gpProbes with at most this percentage mismatches
int 
ApproxMatch(bool bAbsHitLen,		// if true then pctMinCoreLen and pctMinMatchLen are absolute length
		   int pctMinCoreLen,		// otherwise core length is a percentage of actual probe length
		   int pctMinMatchLen,			// otherwise minimum match length is a percentage of actual probe length
	       int MaxMismatches,		// hits to have at most this many mismatches in the left or right flank around exactly matching core
		   etRPTMasking RptMask,	// repeat mask processing mode
		   int RptThres,			// filter out matches with repeats above this percentage threshold
			int hRsltsFile,			// write hits to this file
			CSeqSfx *pSfxArray,	// hits against sequences in this multi-suffix array
			char Strand)			// '*':hits against both strands, '+' strand only, '-' strand only
{
unsigned int LeftMaxExtend;
unsigned int RightMaxExtend;
int MinCoreLen;
int MinMatchLen;
CStopWatch CurTime;
unsigned long Elapsed;
unsigned long PrvTime;

sAMProbe *pCurDescr;
etSeqBase *pProbeBases;
int ProbeIdx;
int Rslt=0;
int NumMatched;
NumMatched = 0;

if(pctMinCoreLen > pctMinMatchLen)
	return(0);

pCurDescr = gProbeDescrs.pDescrs;
for(ProbeIdx = 0; ProbeIdx < gProbeDescrs.NumDescrs; ProbeIdx++, pCurDescr++)
	{
	Elapsed = CurTime.ReadUSecs();
	if(!(pCurDescr->ProbID % 1000) || pCurDescr->ProbID == 1 || Elapsed > PrvTime)
		{
		PrvTime = Elapsed;
		// NOTE: tabs in following printf are simply a poor mans clear to end of line hack!
		printf("\r\t\t\t\t\t\t\t\rProbe %d '%s'",pCurDescr->ProbID,pCurDescr->szDescr);
		}

	if(gProcMode == ePMHitLoci && ((pCurDescr->TotMinHits + pCurDescr->TotPlusHits) > gMaxProbeHits))
		continue;

	if(bAbsHitLen) 
		{
		if(pCurDescr->ProbeLen < pctMinCoreLen)
			continue;
		if(pCurDescr->ProbeLen < pctMinMatchLen)
			continue;
		MinCoreLen = pctMinCoreLen;
		MinMatchLen = pctMinMatchLen;
		}
	else
		{
		if(pctMinCoreLen == 100)
			MinCoreLen = pCurDescr->ProbeLen;
		else
			MinCoreLen = (pCurDescr->ProbeLen * pctMinCoreLen) / 100;

		if(pctMinMatchLen == 100)
			MinMatchLen = pCurDescr->ProbeLen;
		else
			MinMatchLen = (pCurDescr->ProbeLen * pctMinMatchLen) / 100;
		}
	if(MinCoreLen < cMinAbsLen)		// ensure min core length is always >= cMinAbsLen
		MinCoreLen = cMinAbsLen;
	if(MinMatchLen < MinCoreLen)		// ensure min match length is always >= MinCoreLen
		MinMatchLen = MinCoreLen;

	if(MinCoreLen > pCurDescr->ProbeLen)
		continue;

	LeftMaxExtend = pCurDescr->ProbeLen - MinCoreLen;
	RightMaxExtend = LeftMaxExtend;

	if((pProbeBases = LoadSeq(&gProbeDescrs,ProbeIdx))==NULL)
		return(0);

	if(Strand == '*' || Strand == '+')
		{
		if(pCurDescr->CurStrand != '+')
			{
			pCurDescr->CurStrand = '+';
			CSeqTrans::ReverseComplement(pCurDescr->ProbeLen,pProbeBases);
			}

		if((Rslt = pSfxArray->LocateNearExacts(pCurDescr,
							pCurDescr->ProbeLen, // probe length
							pProbeBases,	//probe sequence
							MinCoreLen,		// exactly matching core must be of at least this length
							LeftMaxExtend,	// try extending left for at most this many bases
							RightMaxExtend,	// try extending right for at most this many bases
							MaxMismatches,	// allow this many mismatches in the left flank
							MaxMismatches,	// allow this many mismatches in the right flank
							MinMatchLen,		// matches must of at least this length
							FilterApproxMatches,	// filter handler function for any matches
							HandleApproxMatches,	// handler for accepted matches
							RptMask,
							RptThres))>=0) // suffix array entry
			NumMatched += Rslt;
		else
			return(Rslt);

		}
	if(Strand == '*' || Strand == '-')
		{
		if(pCurDescr->CurStrand != '-')
			{
			pCurDescr->CurStrand = '-';
			CSeqTrans::ReverseComplement(pCurDescr->ProbeLen,pProbeBases);
			}
		// no point in continuing with this probe if already had more hits than max requested
		if(gProcMode == ePMHitLoci && ((pCurDescr->TotMinHits + pCurDescr->TotPlusHits) > gMaxProbeHits))
			continue;

		Rslt = pSfxArray->LocateNearExacts(pCurDescr,
							pCurDescr->ProbeLen, // probe length
							pProbeBases,	//probe sequence
							MinCoreLen,		// exactly matching core must be of at least this length
							LeftMaxExtend,	// try extending left for at most this many bases
							RightMaxExtend,	// try extending right for at most this many bases
							MaxMismatches,	// allow this many mismatches in the left flank
							MaxMismatches,	// allow this many mismatches in the right flank
							MinMatchLen,		// matches must of at least this length
							FilterApproxMatches,	// filter handler function for any matches
							HandleApproxMatches,	// handler for accepted matches
							RptMask,
							RptThres);
		if(Rslt >= 0)		
			NumMatched += Rslt;
		else
			return(Rslt);

		}
	}
return(NumMatched);
}

int
LoadFastaDescriptors(tsSeqDescrs *pSeqDescrs,   // where to load descriptors into
						 CFasta *pFasta,		// preopened fasta file containing probe sequences
			   char Strand)		// treat sequences as being on this strand
{
int Rslt;
sAMProbe *pCurDescr;

char szDescription[cBSFDescriptionSize];
INT64 DescrFileOfs;
int SeqLen;
bool bDescriptor;

memset(pSeqDescrs,0,sizeof(tsSeqDescrs));
if((pSeqDescrs->pDescrs = new sAMProbe [cAllocAMDescrs])==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate memory for %d probe descriptors",cAllocAMDescrs);
	return(eBSFerrMem);
	}
pSeqDescrs->AllocDescrs = cAllocAMDescrs;
pSeqDescrs->MinSeqLen = -1;

// iterate over each probe descriptor
bDescriptor = false;
while((Rslt = SeqLen = pFasta->ReadSequence()) > eBSFSuccess)
	{
	if(SeqLen == eBSFFastaDescr)		// just read a descriptor line
		{
		pFasta->ReadDescriptor(szDescription,cBSFDescriptionSize);
		DescrFileOfs = pFasta->GetDescrFileOfs();
		bDescriptor = true;
		continue;
		}
	else									// just read sequence
		if(!bDescriptor)					// if there was no descriptor then dummy up one...
			{
			sprintf(szDescription,"Probe%d",gProbeDescrs.NumDescrs+1);
			DescrFileOfs = 0;
			bDescriptor = true;
			}
	
	if(pSeqDescrs->MaxSeqLen == 0 || SeqLen > pSeqDescrs->MaxSeqLen)
		pSeqDescrs->MaxSeqLen = SeqLen;
	if(pSeqDescrs->MinSeqLen == -1 || SeqLen < pSeqDescrs->MinSeqLen)
		pSeqDescrs->MinSeqLen = SeqLen;

	if(pSeqDescrs->NumDescrs == pSeqDescrs->AllocDescrs)		// need to allocate for more?
		{
		if((pCurDescr = (sAMProbe *)new sAMProbe [pSeqDescrs->NumDescrs + cAllocAMDescrs])==NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate memory for %d probe descriptors",pSeqDescrs->NumDescrs + cAllocAMDescrs);
			return(eBSFerrMem);
			}
		memcpy(pCurDescr,pSeqDescrs->pDescrs,sizeof(sAMProbe) * pSeqDescrs->NumDescrs);
		delete pSeqDescrs->pDescrs;
		pSeqDescrs->pDescrs = pCurDescr;
		pSeqDescrs->AllocDescrs += cAllocAMDescrs;
		}

	pCurDescr = &pSeqDescrs->pDescrs[pSeqDescrs->NumDescrs];
	pCurDescr->ProbID = ++pSeqDescrs->NumDescrs;
	pCurDescr->ProbeLen = SeqLen;
	pSeqDescrs->TotSeqLen += SeqLen;				
	pCurDescr->TotMinHits = 0;
	pCurDescr->TotPlusHits = 0;
	pCurDescr->Strand = Strand;
	pCurDescr->CurStrand = Strand;
	TrimWhitespace(szDescription);
	strncpy(pCurDescr->szDescr,szDescription,sizeof(pCurDescr->szDescr)-1);
	pCurDescr->szDescr[sizeof(pCurDescr->szDescr)-1] = '\0';
	pCurDescr->SeqRefLoc = DescrFileOfs;
	pCurDescr->CacheSeqOfs = -1;
	bDescriptor = false;
	}
if(Rslt < eBSFSuccess)
	{
	if(pSeqDescrs->pDescrs != NULL)
		{
		delete pSeqDescrs->pDescrs;
		memset(pSeqDescrs,0,sizeof(tsSeqDescrs));
		return(Rslt);
		}
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessFastaFile [%s] %s",pFasta->ErrText((teBSFrsltCodes)Rslt),pFasta->GetErrMsg());
	return(Rslt);
	}
if(Rslt < eBSFSuccess)
	return(Rslt);
return(pSeqDescrs->NumDescrs);
}



// LoadBioSeqDescriptors
// Loads bioseq file probe descriptors, sequences not loaded
int 
LoadBioSeqDescriptors(tsSeqDescrs *pSeqDescrs,   // where to load descriptors into
		    CBioSeqFile *pBioSeqFile, // preopened bioseq file containing probes
		   char Strand)				// assume sequences are on this strand
{
sAMProbe *pCurDescr;
int AllocProbeLen;
int CurProbeEntry;
int Rslt;
int SeqLen;
int NumEntries;

char szProbeName[cMaxDatasetSpeciesChrom+1];
char szProbeDescr[cMaxDatasetSpeciesChrom+1];

memset(pSeqDescrs,0,sizeof(*pSeqDescrs));

if((Rslt = NumEntries = pBioSeqFile->NumEntries()) < 1)
	return(Rslt);
if((pSeqDescrs->pDescrs = new sAMProbe [NumEntries])==NULL)
	return(eBSFerrMem);
pSeqDescrs->AllocDescrs = NumEntries;
pSeqDescrs->MinSeqLen = -1;
pCurDescr = pSeqDescrs->pDescrs;

AllocProbeLen = 0;
Rslt = 0;
CurProbeEntry = 0;
while((Rslt = CurProbeEntry = pBioSeqFile->Next(CurProbeEntry)) > 0)
	{
	pBioSeqFile->GetNameDescription(CurProbeEntry,sizeof(szProbeName),szProbeName,sizeof(szProbeDescr),szProbeDescr);
	SeqLen = pBioSeqFile->GetDataLen(CurProbeEntry);
	if(pSeqDescrs->MaxSeqLen == 0 || SeqLen > pSeqDescrs->MaxSeqLen)
		pSeqDescrs->MaxSeqLen = SeqLen;
	if(pSeqDescrs->MinSeqLen == -1 || SeqLen < pSeqDescrs->MinSeqLen)
		pSeqDescrs->MinSeqLen = SeqLen;
	pCurDescr->ProbID = ++pSeqDescrs->NumDescrs;
	pCurDescr->SeqRefLoc = (INT64)CurProbeEntry;
	pCurDescr->CacheSeqOfs = -1;
	pCurDescr->ProbeLen = SeqLen;
	pSeqDescrs->TotSeqLen += (unsigned int)SeqLen;				
	pCurDescr->TotMinHits = 0;
	pCurDescr->TotPlusHits = 0;
	pCurDescr->Strand = Strand;
	pCurDescr->CurStrand = Strand;
	TrimWhitespace(szProbeName);
	strncpy(pCurDescr->szDescr,szProbeName,sizeof(pCurDescr->szDescr)-1);
	pCurDescr->szDescr[sizeof(pCurDescr->szDescr)-1] = '\0';
	pCurDescr += 1;
	}

if(Rslt < eBSFSuccess && Rslt != eBSFerrEntry)
	{
	if(pSeqDescrs->pDescrs != NULL)
		{
		delete pSeqDescrs->pDescrs;
		memset(pSeqDescrs,0,sizeof(*pSeqDescrs));
		}
	return(Rslt);
	}
return(pSeqDescrs->NumDescrs);
}

// load sequence
etSeqBase *
LoadSeq(tsSeqDescrs *pDescrs,int DescrIdx)
{
int Rslt;
sAMProbe *pCurDescr;
CFasta *pFasta = NULL;
CBioSeqFile *pBioSeqFile = NULL;

if(pDescrs == NULL || gpProbeFile == NULL || pDescrs->pSeq == NULL || DescrIdx >= pDescrs->NumDescrs)
	return(NULL);

pCurDescr = &pDescrs->pDescrs[DescrIdx];
if(pDescrs->pCacheSeq != NULL && pCurDescr->CacheSeqOfs >= 0)		// with any luck sequence may already have been cached
	{
	memcpy(pDescrs->pSeq,&pDescrs->pCacheSeq[pCurDescr->CacheSeqOfs],pCurDescr->ProbeLen);
	pDescrs->CurSeqLen = pCurDescr->ProbeLen;
	pCurDescr->CurStrand = pCurDescr->Strand;
	return(pDescrs->pSeq);
	}

switch(gProbeFileType) {
	case eFTFasta:		// Fasta
		pFasta = (CFasta *)gpProbeFile;
		if((Rslt=pFasta->Reset(pCurDescr->SeqRefLoc))<eBSFSuccess)
			return(NULL);
		if((Rslt = pDescrs->CurSeqLen = pFasta->ReadSequence(pDescrs->pSeq,pCurDescr->ProbeLen)) <= eBSFSuccess)
			return(NULL);
		if(Rslt == eBSFFastaDescr)		// just read a descriptor line, slough and read sequence
			if((Rslt = pDescrs->CurSeqLen = pFasta->ReadSequence(pDescrs->pSeq,pCurDescr->ProbeLen)) <= eBSFSuccess || Rslt == eBSFFastaDescr)
				return(NULL);
		break;

	case eFTBioSeq:
		pBioSeqFile = (CBioSeqFile *)gpProbeFile;
		pDescrs->CurSeqLen = pBioSeqFile->GetDataLen((tBSFEntryID)pCurDescr->SeqRefLoc);
		pBioSeqFile->GetData((tBSFEntryID)pCurDescr->SeqRefLoc,eSeqBaseType,0,pDescrs->pSeq,pCurDescr->ProbeLen);
		break;
	}
CSeqTrans::RemoveMasking(pDescrs->pSeq,pDescrs->CurSeqLen);
pCurDescr->CurStrand = pCurDescr->Strand;

if(pDescrs->pCacheSeq != NULL)
	{
	if((pDescrs->CurCacheLen + pDescrs->CurSeqLen) <= pDescrs->AllocCacheLen)
		{
		memcpy(&pDescrs->pCacheSeq[pDescrs->CurCacheLen],pDescrs->pSeq,pDescrs->CurSeqLen);
		pCurDescr->CacheSeqOfs = pDescrs->CurCacheLen;
		pDescrs->CurCacheLen += pDescrs->CurSeqLen;
		}
	}
return(pDescrs->pSeq);
}


CSeqSfx *
LoadTargSfx(void)
{
int Rslt;
int TargIdx;
sAMProbe *pTargDescr;
CFasta *pFasta = NULL;
CBioSeqFile *pBioSeqFile = NULL;

CSeqSfx * pTargSfx = new CSeqSfx();
if(pTargSfx == NULL)
	return(NULL);

	// now to generate memory suffix array containing all targets
if((gTargDescrs.pSeq = new etSeqBase[gTargDescrs.MaxSeqLen + 1])==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate memory for maximal sized target of %d",gTargDescrs.MaxSeqLen);
	return(NULL);
	}
gTargDescrs.AllocSeqBases = gTargDescrs.MaxSeqLen+1;	
gTargDescrs.CurSeqLen = 0;
pTargDescr = gTargDescrs.pDescrs; 
for(TargIdx = 0; TargIdx < gTargDescrs.NumDescrs; TargIdx++,pTargDescr++)
	{
	switch(gTargFileType) {
		case eFTFasta:		// Fasta
			pFasta = (CFasta *)gpTargFile;
			if((Rslt=pFasta->Reset(pTargDescr->SeqRefLoc))<eBSFSuccess)
				return(NULL);
			if((Rslt = gTargDescrs.CurSeqLen = pFasta->ReadSequence(gTargDescrs.pSeq,gTargDescrs.AllocSeqBases)) <= eBSFSuccess)
				return(NULL);
			if(Rslt == eBSFFastaDescr)		// just read a descriptor line, slough and read sequence
				if((Rslt = gTargDescrs.CurSeqLen = pFasta->ReadSequence(gTargDescrs.pSeq,gTargDescrs.AllocSeqBases)) <= eBSFSuccess || Rslt == eBSFFastaDescr)
					return(NULL);
			break;

		case eFTBioSeq:
			pBioSeqFile = (CBioSeqFile *)gpTargFile;
			gTargDescrs.CurSeqLen = pBioSeqFile->GetDataLen((tBSFEntryID)pTargDescr->SeqRefLoc);
			pBioSeqFile->GetData((tBSFEntryID)pTargDescr->SeqRefLoc,eSeqBaseType,0,gTargDescrs.pSeq,gTargDescrs.AllocSeqBases);
			break;
		}

	if((Rslt=pTargSfx->Add(pTargDescr->ProbID,pTargDescr->ProbeLen,gTargDescrs.pSeq)) < eBSFSuccess)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Failed added suffix entry at target index %d",TargIdx);
		return(NULL);
		}
	}
delete gTargDescrs.pSeq;
gTargDescrs.pSeq = NULL;
gTargDescrs.AllocSeqBases = 0;
gTargDescrs.CurSeqLen = 0;
pTargSfx->GenSfx();
return(pTargSfx);
}
