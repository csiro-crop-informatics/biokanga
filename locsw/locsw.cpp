// locsw.cpp : Defines the entry point for the console application.
// Processes probes against targets using using optimised smith/waterman
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

const unsigned int cProgVer = 100;		// increment with each release

CStopWatch gStopWatch;
CDiagnostics gDiagnostics;				// for writing diagnostics messages to log file
char gszProcName[_MAX_FNAME];			// process name

const int cAllocMatchLoci = 100000;			// allocate match loci blocks holding this many match loci

const int cMinProbeLen = 4;				    // minimum length probes
const int cDfltMinProbeLen = 20;			// default minimum length probe
const int cMaxProbeLen = 10000;				// maximum length probes

const int cDfltMismatches = 5;				// default number of mismatches
const int cMaxMismatches = 100;				// max number of mismatches
const int cMaxInDels = 10;					// max number of InDels
const int cMaxInDelExtns = 100;				// max total length of InDels

const int cMaxEvtScore = 20;				// max score for any match,mismatch,InDel and InDel extension
const double cDfltMinAccumScore = 1.0;		// default minimim acculated normalised score required before reporting hit
const int cDfltMtchScore = 3;				// default score for a match (added to accumulated score)
const int cDfltMismtchScore = 1;			// default score for a mismatch (subtracted from accumulated score)
const int cDfltInDelScore = 5;				// default score for InDel opening (subtracted from accumulated score)
const int cDfltInDelExtnScore = 1;			// default score for InDel extension (subtracted from accumulated score)


typedef enum TAG_eProcMode {				// processing modes
	ePMHitLoci = 0,
	ePMHitCnts
	} teProcMode;

typedef struct TAG_sScores {
	char Strand;			// target strand to be processed - '+', '-' or '*' for both
	int hRsltsFile;			// handle for opened results file

	double ThresAccumScore;	// report all normalised accumulative scores which are >= this threshold and if
	int ThresAccumMismatches; // no more than this number of mismatches and
	int ThresAccumInDels;	// no more than this number of InDels and
	int ThresAccumInDelLen; // also total length of all InDels is less than this value

	int Match;				// (added to AccumScore) score for exact matches
	int Mismatch;			// (subtracted from AccumScore) score for mismatches
	int InDelOpen;			// (subtracted from AccumScore) score for InDel open
	int InDelExtn;			// (subtracted from AccumScore) score for InDel extensions

	int MinScore;			// terminate current search branch if random walk accumulative score drops below this score
	int MaxInDelLen;		// terminate search branch if any single InDel length is longer than this value

	int MinHits;			// report if number of hits to probe is >= MinHits
	int MaxHits;				// report if number of hits to probe is <= MaxHits
} tsScores;



typedef struct TAG_sFSAMProbe {
	struct TAG_sFSAMProbe *pNext;			// probes are linked
	int ProbID;								// uniquely identifies this probe
	char szDescr[80];						// as parsed from fasta descriptor
	bool bAtMaxHits;						// true if this probe has already at max allowed hits
	int TotPlusHits;						// total hits on '+' strand
	int TotMinHits;							// total hits on '-' strand
	int ProbeLen;							// number of bases in Bases[]
	etSeqBase Bases[1];						// to hold probe sequence
} sFSAMProbe;

typedef struct TAG_sAccumWalkScore {
	char Strand;			// target strand being processed
	int PeakAccumScore;		// peak accumulative score in random walk
	int ProbeIdx;			// inital index into probe sequence
	int TargStartIdx;		// initial index into targ sequence at which matching started
	int TargEndIdx;			// index into targ sequence at which match ended
	
	sFSAMProbe *pCurProbe;	// current probe
	int ProbeLen;			// probe length
	etSeqBase *pProbe;		// probe sequence (must be eBaseEOS terminated)
	int TargEntryID;		// target entry identifier
	int TargLen;			// target length
	etSeqBase *pTarg;		// targ sequence (must be eBaseEOS terminated)
	int CurInDelLen;		// current length of any InDel
	int AccumScore;			// accumulative walk score
	int TotMatches;			// total number of matches
	int TotMismatches;		// total number of mismatches
	int TotInDels;			// total number of InDels
	int TotInDelLen;		// total length of all InDels
} tsAccumWalkScore;

// holds match loci 
typedef struct TAG_sMatchLoci {
	sFSAMProbe *pProbe;					// probe which matched
	char Strand;						// was on to this strand
	tBSFEntryID EntryID;				// matches on to this chromosome/sequence
	int StartLoci;						// and at this loci
	int MatchLen;						// and is this length
	double NormScore;					// with this normalised score
} sMatchLoci;

typedef struct TAG_sMatchLociBlock {
	struct TAG_sMatchLociBlock *pNext;	// pts to next allocated block of match loci
	int NumLoci;						// number of loci currently in this block
	sMatchLoci Loci[cAllocMatchLoci];
	} sMatchLociBlock;

int gNumProbes = 0;						// number of probes loaded
int gNumProbesAtMaxHits = 0;			// number of probes which have MaxHits to target
int gMaxProbeLen = 0;					// max length of any probe
int gMinProbeLen = 0;					// min length of any probe
sFSAMProbe *gpProbes = NULL;			// pts to linked list of sAMProbes

int gMatchCnt;							// global match count
sMatchLociBlock *gpMatchBlocks = NULL;	// pts to linked list of match blocks

int Process(teProcMode Mode,
			char *pszProbeFile,			// probes from this multfasta file
			char *pszTargFile,			// targets from this bioseq file
			char *pszRsltsFile,			// hits into this results file
			int MinLength,				// process probes which are of at least this length
			int MaxLength,				// process probes which are no longer than this length
			int MinHits,				// only report probes which have at least this many hits
			int MaxHits,				// only report probes which have at most this many hits
			char *pszChroms,			// only process these target chroms
			char Strand,				// process against this target strand - '+','-' or '*' for both
			int MaxMismatches,			// accept up to at most this many mismatches in any hit
			int MaxInDels,				// accept up to at most this many InDels
			int MaxInDelLen,			// accept total InDel lengths of upto this many bases
			double NormAccumScore,		// report hits with normalised scores which are at least this value 
			int MatchScore,				// score each match this value
			int MismatchScore,			// score each mismatch this value
			int InDelOpen,				// score each InDel opening this value
			int InDelExtn);				// score each InDel extension this value

int ReportMatches(CBioSeqFile *pBioSeq,int hRsltsFile,int MinHits,int MaxHits);
int ReportCounts(int hRsltsFile,int MinHits,int MaxHits,char *pszChrom);
int ReportMisses(CBioSeqFile *pBioSeq,int hRsltsFile,char *pszChrom);

int AddMatch(sFSAMProbe *pProbe,char Strand,tBSFEntryID EntryID,int Loci,double NormScore);
void DeleteMatches(void);
int LoadProbes(char *pszProbeFile,int MaxMismatches,int MinLenght, int MaxLength);
void DeleteProbes(void);
int										// returns the number of hits by probes against the chromosome sequence
FindMatches(teProcMode Mode,
		   tsScores *pSWScores,			// scoring to use
		   tBSFEntryID CurEntry,		// current target entry identifier
		   unsigned char *pTargSeq,		// pts to target sequence
		   int TargSeqLen,				// target sequence length
		   char *pszChrom);				// target is on this chromosome 

int	CleanNameList(char *pszTxt);
bool NameInList(char *pszNameList, char *pszName);

int ProcessAlignment(double NormScore,tsAccumWalkScore *pAccumScore,	// current accmulative walk score
			  tsScores *pScores);
int FindSWAlignment(tsAccumWalkScore *pAccumScore,	// current accmulative walk score
			  tsScores *pScores);


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
int iMaxMismatches;			//accept up to at most this many mismatches in any hit (0..100), default is 5");
int iMaxInDels;				//accept up to at most this many InDels (0..10), default is 1");
int iMaxInDelLen;			//accept total InDel lengths of upto this many bases (1..100), default is 5");
double dNormAccumScore;		//report hits with normalised scores which are at least this value ( >= 0.0), default is 1.0");
int iMatchScore;			//score add each match this value (0..20), default is 3");
int iMismatchScore;			//score subtract each mismatch this value (0..20), default is 1");
int iInDelOpenScore;		//score subtract InDel opening this value (0..20), default is 5");
int iInDelExtnScore;		//score subtract InDel extension this value (0..20), default is 1");
char cStrand;				//Strand '+' or '-', default is '*' for both");
int iMinLength;				//Filter out probes of less than this length (default: 20, min: 4)");
int iMaxLength;				//Filter out probes of longer than this length (default: 10000, max: 10000)");
int iMinHits;				//Only report probes with hit counts of at least this many (default: 1)");
int iMaxHits;				//Only report probes with hit counts of no more than this many (default: 1000)");
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

struct arg_int *Mode=arg_int0("m", "mode","<int>","Processing mode: 0 Target hit loci, 1: Probe hit counts");
struct arg_file *ProbeFile = arg_file1("i","probes","<file>",	"input from multifasta (.fa) probe file");
struct arg_file *TargFile = arg_file1("I","targets","<file>",	"input from bioseq (.seq) target file");
struct arg_file *RsltsFile = arg_file1("o","results","<file>",	"output match results to this file as CSV");
struct arg_int *MaxMismatches=arg_int0("t", "maxmissmatches","<int>",	"accept up to at most this many mismatches in any hit (0..100), default is 5");
struct arg_int *MaxInDels=arg_int0("j", "maxindels","<int>",			"accept up to at most this many InDels (0..10), default is 0");
struct arg_int *MaxInDelLen=arg_int0("J", "maxtotindellen","<int>",		"accept total InDel lengths of upto this many bases (1..100), default is 5");
struct arg_dbl *NormAccumScore=arg_dbl0("k", "mintotscore","<int>",		"only report hits with normalised scores which are at least this value ( >= 0.0), default is 1.0");
struct arg_int *MatchScore=arg_int0("K", "matchscore","<int>",			"score add each match this value (0..20), default is 3");
struct arg_int *MismatchScore=arg_int0("T", "mismatchscore","<int>",	"score subtract each mismatch this value (0..20), default is 1");
struct arg_int *InDelOpenScore=arg_int0("q", "indelscore","<int>",		"score subtract InDel opening this value (0..20), default is 5");
struct arg_int *InDelExtnScore=arg_int0("Q", "indelextnscore","<int>",	"score subtract InDel extension this value (0..20), default is 1");
struct arg_str *Strand=arg_str0("s", "strand",	"<string>",		"Strand '+' or '-', default is '*' for both");
struct arg_int *MinLength=arg_int0("l", "Minlen","<int>",		"Filter out probes of less than this length (default: 20, min: 4)");
struct arg_int *MaxLength=arg_int0("L", "Maxlen","<int>",		"Filter out probes of longer than this length (default: 10000, max: 10000)");
struct arg_int *MinHits=arg_int0("r", "minhits","<int>",		"Only report probes with hit counts of at least this many (default: 1)");
struct arg_int *MaxHits=arg_int0("R", "maxhits","<int>",		"Only report probes with hit counts of no more than this many (default: 1000)");
struct arg_str *Chroms=arg_strn("c", "chroms",	"<string>",	0,50,	"Comma/space separated list of chromosomes to process (default: all)");


struct arg_end *end = arg_end(20);

void *argtable[] = {help,version,FileLogLevel,ScreenLogLevel,LogFile,
					Mode,ProbeFile,TargFile,RsltsFile,MaxMismatches,MaxInDels,MaxInDelLen,NormAccumScore,MatchScore,MismatchScore,InDelOpenScore,
					InDelExtnScore,Strand,MinLength,MaxLength,MinHits,MaxHits,Chroms,
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

	iMode = Mode->count ? Mode->ival[0] : ePMHitLoci;
	if(iMode < ePMHitLoci || iMode > ePMHitCnts)
		{
		printf("\nError: Processing mode '-m%d' not in range %d..%d\n",iMode,ePMHitLoci,ePMHitCnts);
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

	iMaxInDels = MaxInDels->count ? MaxInDels->ival[0] : 0;
	if(iMaxInDels < 0 || iMaxInDels > cMaxInDels)
		{
		printf("\nError: Maximum number of InDels '-j%d' not in range 0..%d\n",iMaxInDels,cMaxInDels);
		exit(1);
		}

	iMaxInDelLen = MaxInDelLen->count ? MaxInDelLen->ival[0] : 5;
	if(iMaxInDelLen < 1 || iMaxInDelLen > cMaxInDelExtns)
		{
		printf("\nError: Maximum total length of InDels '-J%d' not in range 0..%d\n",iMaxInDelLen,cMaxInDelExtns);
		exit(1);
		}


	dNormAccumScore = NormAccumScore->count ? NormAccumScore->dval[0] : cDfltMinAccumScore;
	if(dNormAccumScore < 0.0)
		{
		printf("\nError: Minimum normalised accumulative score '-k%f' must be >= 0.0\n",dNormAccumScore);
		exit(1);
		}

	iMatchScore = MatchScore->count ? MatchScore->ival[0] : cDfltMtchScore;
	if(iMatchScore < 0 || iMatchScore > cMaxEvtScore)
		{
		printf("\nError: Match score '-K%d' not in range 0..%d\n",iMatchScore,cMaxEvtScore);
		exit(1);
		}

	iMismatchScore = MismatchScore->count ? MismatchScore->ival[0] : cDfltMismtchScore;
	if(iMismatchScore < 0 || iMismatchScore > cMaxEvtScore)
		{
		printf("\nError: Mismatch score '-T%d' not in range 0..%d\n",iMismatchScore,cMaxEvtScore);
		exit(1);
		}

	iInDelOpenScore = InDelOpenScore->count ? InDelOpenScore->ival[0] : cDfltInDelScore;
	if(iInDelOpenScore < 0 || iInDelOpenScore > cMaxEvtScore)
		{
		printf("\nError: InDel open score '-q%d' not in range 0..%d\n",iInDelOpenScore,cMaxEvtScore);
		exit(1);
		}

	iInDelExtnScore = InDelExtnScore->count ? InDelExtnScore->ival[0] : cDfltInDelExtnScore;
	if(iInDelExtnScore < 0 || iInDelExtnScore > cMaxEvtScore)
		{
		printf("\nError: InDel extension '-Q%d' not in range 0..%d\n",iInDelExtnScore,cMaxEvtScore);
		exit(1);
		}

	iMinHits = MinHits->count ? MinHits->ival[0] : 1;
	if(iMinHits < 1)
		{
		printf("\nError: Minimum hits '-r%d' must be at least 1\n",iMinHits);
		exit(1);
		}
	iMaxHits = MaxHits->count ? MaxHits->ival[0] : 1000;
	if(iMaxHits < iMinHits)
		{
		printf("\nError: Maximum hits '-R%d' must be at least minimum hits %d\n",iMaxHits,iMinHits);
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
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Probes from multifasta file: %s",szProbeFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Targets from bioseq sequence file: %s",szTargFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Results to file: '%s'", szRsltsFile);
	if(szChroms[0] != '\0')
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Only process these target chromosomes: %s",szChroms);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"accept up to at most this many InDels: %d",iMaxInDels);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"accept total InDel lengths of upto: %d",iMaxInDelLen);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"report hits with normalised scores which are at least: %f",dNormAccumScore);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"score each base match: %d",iMatchScore);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"score subtract each base mismatch: %d",iMismatchScore);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"core subtract InDel opening: %d",iInDelOpenScore);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"score subtract InDel extension: %d",iInDelExtnScore);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Filter out probes with length less than: %d",iMinLength);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Filter out probes with length longer than: %d",iMaxLength);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Filter out probes with less than %d hits",iMinHits);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Filter out probes with more than %d hits",iMaxHits);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Target strand: '%c'",cStrand);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Max mismatches: %d",iMaxMismatches);

#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif

	gStopWatch.Start();
	Rslt = Process((teProcMode)iMode,
				szProbeFile,			// probes from this multfasta file
				szTargFile,				// targets from this bioseq file
				szRsltsFile,			// hits into this results file
				iMinLength,				// process probes which are of at least this length
				iMaxLength,				// process probes which are no longer than this length
				iMinHits,				// only report probes which have at least this many hits
				iMaxHits,				// only report probes which have at most this many hits
				szChroms,				// only process these target chroms
				cStrand,				// process against this target strand - '+','-' or '*' for both
				iMaxMismatches,			// accept up to at most this many mismatches in any hit
				iMaxInDels,				// accept up to at most this many InDels
				iMaxInDelLen,			// accept total InDel lengths of upto this many bases
				dNormAccumScore,		// report hits with normalised scores which are at least this value 
				iMatchScore,			// score each match this value
				iMismatchScore,			// score each mismatch this value
				iInDelOpenScore,		// score each InDel opening this value
				iInDelExtnScore);		// score each InDel extension this value

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


int Process(teProcMode Mode,
			char *pszProbeFile,			// probes from this multfasta file
			char *pszTargFile,			// targets from this bioseq file
			char *pszRsltsFile,			// hits into this results file
			int MinLength,				// process probes which are of at least this length
			int MaxLength,				// process probes which are no longer than this length
			int MinHits,				// only report probes which have at least this many hits
			int MaxHits,				// only report probes which have at most this many hits
			char *pszChroms,			// only process these target chroms
			char Strand,				// process against this target strand - '+','-' or '*' for both
			int MaxMismatches,			// accept up to at most this many mismatches in any hit
			int MaxInDels,				// accept up to at most this many InDels
			int MaxInDelLen,			// accept total InDel lengths of upto this many bases
			double NormAccumScore,		// report hits with normalised scores which are at least this value 
			int MatchScore,				// score each match this value
			int MismatchScore,			// score each mismatch this value
			int InDelOpen,				// score each InDel opening this value
			int InDelExtn)				// score each InDel extension this value
{
CBioSeqFile *pBioSeq;
unsigned char *pTargSeq;			// to hold sequence for each chromsome
int AllocdTargSeqLen;				// length of currently allocated pTargSeq
int hRsltsFile;
char szChromName[cBSFSourceSize+1];	// to hold current chromosome name
tBSFEntryID CurEntry;				// current entry being processed
int TargSeqLen;						// length of chromosome sequence being processed
int Rslt;							// used to hold processing result

tsScores SWScores;
memset(&SWScores,0,sizeof(SWScores));

SWScores.hRsltsFile = -1;
SWScores.Strand = Strand;

SWScores.ThresAccumScore = NormAccumScore,		// report hits with normalised scores which are at least this value;				
SWScores.ThresAccumMismatches = MaxMismatches;	// no more than this number of mismatches and
SWScores.ThresAccumInDels = MaxInDels;		// no more than this number of InDels and
SWScores.ThresAccumInDelLen = MaxInDelLen;	// also total length of all InDels is less than this value

SWScores.Match = MatchScore;			// (added to AccumScore) score for exact matches
SWScores.Mismatch = MismatchScore;		// (subtracted from AccumScore) score for mismatches
SWScores.InDelOpen = InDelOpen;			// (subtracted from AccumScore) score for InDel open
SWScores.InDelExtn = InDelExtn;			// (subtracted from AccumScore) score for InDel extensions

SWScores.MinScore = -100;					// terminate search branch if random walk accumulative score drops below this score
SWScores.MaxInDelLen = MaxInDelLen;			// terminate search branch if any single InDel length is longer than this value

SWScores.MinHits = MinHits;			// report if number of hits to probe is >= MinHits
SWScores.MaxHits = MaxHits;			// report if number of hits to probe is <= MaxHits


if((Rslt=LoadProbes(pszProbeFile,MaxMismatches,MinLength,MaxLength))<=eBSFSuccess)
	return(Rslt);

// open bioseq file containing chromosome sequences
if((pBioSeq = new CBioSeqFile()) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CBioSeqFile");
	DeleteProbes();
	return(eBSFerrObj);
	}
if((Rslt=pBioSeq->Open(pszTargFile,cBSFTypeSeq))!=eBSFSuccess)
	{
	while(pBioSeq->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,pBioSeq->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open input bioseq sequence file '%s'",pszTargFile);
	DeleteProbes();
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
	DeleteProbes();
	delete pBioSeq;
	return(eBSFerrOpnFile);
	}
SWScores.hRsltsFile = hRsltsFile;

	// iterate over each entry (chromosome) in target sequence file
pTargSeq = NULL;
TargSeqLen = 0;
AllocdTargSeqLen = 0;
CurEntry = 0;
gMatchCnt = 0;
Rslt = eBSFSuccess;
while((CurEntry = pBioSeq->Next(CurEntry)) > 0)
	{
		// get entry name - assume it will be a chromosome/contig name
	pBioSeq->GetName(CurEntry,sizeof(szChromName),szChromName);
		// check only processing specific chromosomes
	if(pszChroms != NULL && pszChroms[0] != '\0' && !NameInList(pszChroms,szChromName))
		continue;

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing: %s",szChromName);

		// need to ensure buffer can hold the sequence
	TargSeqLen = pBioSeq->GetDataLen(CurEntry);
	if(pTargSeq == NULL || (TargSeqLen+1) > AllocdTargSeqLen)
		{
		if(pTargSeq != NULL)
			delete pTargSeq;
		AllocdTargSeqLen = TargSeqLen + TargSeqLen/10; // alloc more than actually required, next sequence may not need to be realloc'd
		pTargSeq = (unsigned char *)new unsigned char [AllocdTargSeqLen];
		if(pTargSeq == NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate memory (%d bytes) for holding chromosome sequence",AllocdTargSeqLen);
			Rslt = eBSFerrMem;
			break;
			}
		}

	pBioSeq->GetData(CurEntry,cBSFTypeSeq,0,pTargSeq,TargSeqLen);
	pTargSeq[TargSeqLen] = eBaseEOS;
	FindMatches(Mode,
		   &SWScores,		// scoring to use
		   CurEntry,		// current target entry identifier
		   pTargSeq,		// pts to target sequence
		   TargSeqLen,		// target sequence length
		   szChromName);	// target is on this chromosome 

	if(gNumProbesAtMaxHits >=  gNumProbes)
		{
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"All probes have had at least %d hits, search terminating",MaxHits);
		break;
		}
	}

if(Mode == ePMHitLoci)
	ReportMatches(pBioSeq,hRsltsFile,MinHits,MaxHits);
else
	ReportCounts(hRsltsFile,MinHits,MaxHits,(char *)"ALLCHROMES");
		
DeleteMatches();
DeleteProbes();
if(hRsltsFile!=-1)
	close(hRsltsFile);

if(pTargSeq!=NULL)
	delete pTargSeq;
if(pBioSeq != NULL)
	delete pBioSeq;
return(Rslt);
}

void
DeleteMatches(void)
{
sMatchLociBlock *pNxtBlock;
sMatchLociBlock *pCurBlock;
if((pNxtBlock = gpMatchBlocks) != NULL)
	while((pCurBlock = pNxtBlock) != NULL) 
		{
		pNxtBlock = pCurBlock->pNext;
		delete pCurBlock;
		}
gpMatchBlocks = NULL;
gMatchCnt = 0;
}


int
AddMatch(sFSAMProbe *pProbe,char Strand,tBSFEntryID EntryID,int StartLoci,int MatchLen,double NormScore)
{
sMatchLoci *pMatch;
sMatchLociBlock *pNewBlock;

if(gpMatchBlocks == NULL ||
   gpMatchBlocks->NumLoci == cAllocMatchLoci)
	{
	if((pNewBlock = new sMatchLociBlock)==NULL)
		return(eBSFerrMem);
	pNewBlock->NumLoci = 0;
	pNewBlock->pNext = gpMatchBlocks;
	gpMatchBlocks = pNewBlock;
	}
pMatch = &gpMatchBlocks->Loci[gpMatchBlocks->NumLoci++];
pMatch->Strand = Strand;
pMatch->EntryID = EntryID;
pMatch->StartLoci = StartLoci;
pMatch->MatchLen = MatchLen;
pMatch->NormScore = NormScore;
pMatch->pProbe = pProbe;
gMatchCnt+=1;
return(eBSFSuccess);
}

void
DeleteProbes(void)
{
sFSAMProbe *pCurProbe;
sFSAMProbe *pNxtProbe;
if((pNxtProbe = gpProbes) != NULL)
	while((pCurProbe = pNxtProbe) != NULL) 
		{
		pNxtProbe = pCurProbe->pNext;
		delete pCurProbe;
		}
gpProbes = NULL;
}


int
LoadProbes(char *pszProbeFile,int MaxMismatches,int MinLength,int MaxLength)
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


gNumProbes = 0;
gNumProbesAtMaxHits = 0;
gpProbes = NULL;

// open fasta file containing probes
if((pFasta = new CFasta()) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CFasta");
	return(eBSFerrObj);
	}
if((Rslt=pFasta->Open(pszProbeFile))!=eBSFSuccess)
	{
	while(pFasta->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,pFasta->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open input fasta file containing probes '%s'",pszProbeFile);
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
					gDiagnostics.DiagOut(eDLInfo,gszProcName,"LoadProbes [%s] truncated probe %s at %d length",pszProbeFile,szDescription, cMaxProbeLen);
				continue;
				}
			sprintf(szDescription,"Probe%d",gNumProbes+1);
			TruncCnt = 0;
			bDescriptor = true;
			}
	if(SeqLen < MinLength)
		{
		UnderLength += 1;
		continue;
		}
	if(SeqLen > MaxLength)
		{
		OverLength += 1;
		continue;
		}
	bFirst = false;

	SeqBuff[SeqLen] = eBaseEOS;

	if(gMaxProbeLen == 0 || SeqLen > gMaxProbeLen)
		gMaxProbeLen = SeqLen;
	if(gMinProbeLen == 0 || SeqLen < gMinProbeLen)
			gMinProbeLen = SeqLen;
	pNxtProbe = (sFSAMProbe *)new unsigned char [ sizeof(sFSAMProbe) + SeqLen + 1];
	pNxtProbe->pNext = NULL;
	pNxtProbe->ProbID = ++gNumProbes;
	pNxtProbe->ProbeLen = SeqLen;
	pNxtProbe->TotMinHits = 0;
	pNxtProbe->TotPlusHits = 0;
	pNxtProbe->bAtMaxHits = false;
	strncpy(pNxtProbe->szDescr,szDescription,sizeof(pNxtProbe->szDescr)-1);
	pNxtProbe->szDescr[sizeof(pNxtProbe->szDescr)-1] = '\0';
	CSeqTrans::RemoveMasking(SeqBuff,SeqLen);
	memcpy(pNxtProbe->Bases,SeqBuff,SeqLen+1);
	if(gpProbes == NULL)
		gpProbes = pCurProbe = pNxtProbe;
	else
		{
		pCurProbe->pNext = pNxtProbe;
		pCurProbe = pNxtProbe;
		}
	bDescriptor = false;
	}
if(Rslt < eBSFSuccess)
	{
	if(gpProbes != NULL)
		DeleteProbes();
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessFastaFile [%s] %s",pFasta->ErrText((teBSFrsltCodes)Rslt),pFasta->GetErrMsg());
	return(Rslt);
	}
delete pFasta;

if(UnderLength)
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Filtered out %d underlength probes",UnderLength);
if(OverLength)
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Filtered out %d overlength probes",OverLength);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Number of probes to process %d",gNumProbes);
return(gNumProbes);
}


int
ReportCounts(int hRsltsFile,int MinHits,int MaxHits,char *pszChrom)
{
int TotHits;
int RsltsLen;
char szOutBuff[0x3fff];
sFSAMProbe *pCurProbe = gpProbes;
while(pCurProbe != NULL)
	{
	TotHits = pCurProbe->TotPlusHits + pCurProbe->TotMinHits;
	if(TotHits >= MinHits && TotHits <= MaxHits)
		{
		RsltsLen = sprintf((char *)szOutBuff,"%d,%d,\"%s\",%d,%d,\"%s\"\n",++gMatchCnt,pCurProbe->ProbID,pCurProbe->szDescr,pCurProbe->TotPlusHits,pCurProbe->TotMinHits,pszChrom);
		CUtility::SafeWrite(hRsltsFile,szOutBuff,RsltsLen);
		}
	pCurProbe->TotPlusHits = 0;
	pCurProbe->TotMinHits = 0;
	pCurProbe = pCurProbe->pNext;
	}
return(gMatchCnt);
}


int ReportMatches(CBioSeqFile *pBioSeq,int hRsltsFile,int MinHits,int MaxHits)
{
sMatchLociBlock *pCurBlock;
sMatchLoci *pCurMatch;
sFSAMProbe *pProbe;
char szChromName[cBSFSourceSize+1];		// to hold cached chromosome name
char szTitle[cBSFSourceSize+1];			// to hold bioseq title used for target assembly name
tBSFEntryID CurEntryID;					// entry cached

char szRsltsBuff[0x7fff];		// used to hold results before writing to disk
int	RsltsLen;							// current length of results in buffer
int TargHits;
int ProbeHits;
int LociIdx;

if((pCurBlock = gpMatchBlocks) == NULL)
	return(0);

pBioSeq->GetTitle(sizeof(szTitle),szTitle);

CurEntryID = 0;
TargHits = 0;
RsltsLen = 0;
do {
	pCurMatch = &pCurBlock->Loci[0];
	for(LociIdx = 0; LociIdx < pCurBlock->NumLoci; LociIdx++,pCurMatch++)
		{
		pProbe = pCurMatch->pProbe;
		ProbeHits = pProbe->TotPlusHits + pProbe->TotMinHits;
		if(ProbeHits < MinHits ||  ProbeHits > MaxHits)
			continue;
		
		if(CurEntryID == 0 || pCurMatch->EntryID != CurEntryID)
			{
			pBioSeq->GetName(pCurMatch->EntryID,sizeof(szChromName),szChromName);
			CurEntryID = pCurMatch->EntryID;
			}

		RsltsLen += sprintf(&szRsltsBuff[RsltsLen],"%d,\"SWaligned\",\"%s\",\"%s\",%u,%u,%d,\"%s\",%d,\"%c\",%.3f\n",
								++TargHits,szTitle,szChromName,pCurMatch->StartLoci,pCurMatch->StartLoci + pCurMatch->MatchLen - 1,
								pCurMatch->MatchLen,pProbe->szDescr,0,pCurMatch->Strand,pCurMatch->NormScore);
		if(RsltsLen > (sizeof(szRsltsBuff)/10))
			{
			CUtility::SafeWrite(hRsltsFile,szRsltsBuff,RsltsLen);
			RsltsLen = 0;
			}

		}
	}
while(pCurBlock = pCurBlock->pNext);
if(RsltsLen)
	CUtility::SafeWrite(hRsltsFile,szRsltsBuff,RsltsLen);
return(TargHits);
}

int										// returns the number of hits by probes against the chromosome sequence
FindMatches(teProcMode Mode,
		   tsScores *pSWScores,			// scoring and processing parameters
		   tBSFEntryID CurEntry,		// current target entry identifier
		   unsigned char *pTargSeq,		// pts to target sequence
		   int TargSeqLen,				// target sequence length
		   char *pszChrom)				// target is on this chromosome 
{
int Rslt;
tsAccumWalkScore CurAccumScore;
sFSAMProbe *pPutProbe;					// putative matching probe
int CurTargIdx;
etSeqBase *pProbe;
etSeqBase *pTarg;
etSeqBase ProbeBase;
etSeqBase TargBase;
int PlusHits;
int MinHits;
int NumCompared;
bool bInDelsAllowed;					// will be true if matches can contain InDels
double NormScore;

PlusHits = 0;
MinHits = 0;
bInDelsAllowed = pSWScores->ThresAccumInDels > 0 ? true : false;
	
// first process the '+' strand
if((pSWScores->Strand == '*' || pSWScores->Strand == '+'))
	{
	for(CurTargIdx = 0; (CurTargIdx < TargSeqLen) && (gNumProbesAtMaxHits <  gNumProbes); CurTargIdx++)
		{
		if(!(CurTargIdx & 0x0fffff))
			printf("\rTarget Loci '+' strand %s: %9.9d",pszChrom,CurTargIdx);		// let user know we haven't crashed...
		pPutProbe = gpProbes;
		while(pPutProbe != NULL)
			{
			if((pPutProbe->TotMinHits + pPutProbe->TotPlusHits) <= pSWScores->MaxHits)
				{
				memset(&CurAccumScore,0,sizeof(CurAccumScore));
				CurAccumScore.Strand = '+';
				CurAccumScore.pCurProbe = pPutProbe;
				CurAccumScore.ProbeLen = pPutProbe->ProbeLen;
				CurAccumScore.pProbe = pPutProbe->Bases;
				CurAccumScore.TargEntryID = CurEntry;
				CurAccumScore.TargLen = TargSeqLen;
				CurAccumScore.pTarg = pTargSeq;
				CurAccumScore.TargStartIdx = CurTargIdx;
				CurAccumScore.TargEndIdx = CurTargIdx;

				if(!bInDelsAllowed)				// if matches not to contain InDels 
					{
					pProbe = pPutProbe->Bases;
					pTarg = (etSeqBase *)&pTargSeq[CurTargIdx];
					NumCompared = 0;
					while((ProbeBase = (*pProbe++ & ~cRptMskFlg)) != eBaseEOS && (TargBase = (*pTarg++ & ~cRptMskFlg)) != eBaseEOS)
						{
						NumCompared += 1;
						if(ProbeBase == TargBase)
							{
							CurAccumScore.TotMatches++;
							CurAccumScore.AccumScore += pSWScores->Match;
							if(CurAccumScore.PeakAccumScore < CurAccumScore.AccumScore)
								CurAccumScore.PeakAccumScore = CurAccumScore.AccumScore;
							}
						else			
							{
							if(++CurAccumScore.TotMismatches > pSWScores->ThresAccumMismatches)
								break;
							CurAccumScore.AccumScore -= pSWScores->Mismatch;
							if(CurAccumScore.AccumScore < pSWScores->MinScore)
								break;
							}
						CurAccumScore.TargEndIdx += 1;
						}

					NormScore = (double)CurAccumScore.PeakAccumScore / pPutProbe->ProbeLen;
					if(NumCompared == CurAccumScore.ProbeLen && 
						CurAccumScore.TotMismatches <= pSWScores->ThresAccumMismatches &&
						NormScore >= pSWScores->ThresAccumScore)
						if((Rslt=ProcessAlignment(NormScore,&CurAccumScore,pSWScores))>=eBSFSuccess)
							{
							pPutProbe->TotPlusHits += 1;
							PlusHits += 1;
							}
					}
				else							// else matches may contain InDels 
					{
					if((Rslt=FindSWAlignment(&CurAccumScore,pSWScores))>=eBSFSuccess)
						{
						pPutProbe->TotPlusHits += Rslt;
						PlusHits += Rslt;
						}
					}
				}
			else
				{
				if(!pPutProbe->bAtMaxHits)
					{
					pPutProbe->bAtMaxHits = true;
					gNumProbesAtMaxHits+=1;
					}
				}
			pPutProbe = pPutProbe->pNext;
			}
		}
	printf("\rTarget Loci '+' strand %s: %9.9d",pszChrom,TargSeqLen);
	}


// and now the '-' strand
if((pSWScores->Strand == '*' || pSWScores->Strand == '-'))
	{
	printf("\n");
	CSeqTrans::ReverseComplement(TargSeqLen,pTargSeq);
	for(CurTargIdx = 0; (CurTargIdx < TargSeqLen) && (gNumProbesAtMaxHits <  gNumProbes); CurTargIdx++)
		{
		if(!(CurTargIdx & 0x0fffff))
			printf("\rTarget Loci '-' strand %s: %9.9d",pszChrom,CurTargIdx);		// let user know we haven't crashed...
		pPutProbe = gpProbes;
		while(pPutProbe != NULL)
			{
			if((pPutProbe->TotMinHits + pPutProbe->TotPlusHits) <= pSWScores->MaxHits)
				{
				memset(&CurAccumScore,0,sizeof(CurAccumScore));
				CurAccumScore.Strand = '-';
				CurAccumScore.pCurProbe = pPutProbe;
				CurAccumScore.ProbeLen = pPutProbe->ProbeLen;
				CurAccumScore.pProbe = pPutProbe->Bases;
				CurAccumScore.ProbeIdx = 0;
				CurAccumScore.TargEntryID = CurEntry;
				CurAccumScore.TargLen = TargSeqLen;
				CurAccumScore.pTarg = pTargSeq;
				CurAccumScore.TargStartIdx = CurTargIdx;
				CurAccumScore.TargEndIdx = CurTargIdx;

				if(!bInDelsAllowed)				// if matches not to contain InDels 
					{
					pProbe = pPutProbe->Bases;
					pTarg = (etSeqBase *)&pTargSeq[CurTargIdx];
					NumCompared = 0;
					while((ProbeBase = (*pProbe++ & ~cRptMskFlg)) != eBaseEOS && (TargBase = (*pTarg++ & ~cRptMskFlg)) != eBaseEOS)
						{
						NumCompared += 1;
						if(ProbeBase == TargBase)
							{
							CurAccumScore.TotMatches++;
							CurAccumScore.AccumScore += pSWScores->Match;
							if(CurAccumScore.PeakAccumScore < CurAccumScore.AccumScore)
								CurAccumScore.PeakAccumScore = CurAccumScore.AccumScore;
							}
						else			
							{
							if(++CurAccumScore.TotMismatches > pSWScores->ThresAccumMismatches)
								break;
							CurAccumScore.AccumScore -= pSWScores->Mismatch;
							if(CurAccumScore.AccumScore < pSWScores->MinScore)
								break;
							}
						CurAccumScore.TargEndIdx += 1;
						}

					NormScore = (double)CurAccumScore.PeakAccumScore / pPutProbe->ProbeLen;
					if(NumCompared == CurAccumScore.ProbeLen && 
						CurAccumScore.TotMismatches <= pSWScores->ThresAccumMismatches &&
						NormScore >= (double)pSWScores->ThresAccumScore)
							if((Rslt=ProcessAlignment(NormScore,&CurAccumScore,pSWScores))>=eBSFSuccess)
								{
								pPutProbe->TotMinHits += 1;
								MinHits += 1;
								}
					}
				else							// else matches allowed to contain InDels 
					{
					if((Rslt=FindSWAlignment(&CurAccumScore,pSWScores))>=eBSFSuccess)
						{
						pPutProbe->TotMinHits += Rslt;
						MinHits += Rslt;
						}
					}
				}
			else
				{
				if(!pPutProbe->bAtMaxHits)
					{
					pPutProbe->bAtMaxHits = true;
					gNumProbesAtMaxHits+=1;
					}
				}
			pPutProbe = pPutProbe->pNext;
			}
		}
	printf("\rTarget Loci '-' strand %s: %9.9d",pszChrom,TargSeqLen);
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Hits on plus strand: %d, on minus strand: %d",PlusHits,MinHits);
return(PlusHits + MinHits);
}


int
ProcessAlignment(double NormScore,tsAccumWalkScore *pAccumScore,	// current accmulative walk score
			  tsScores *pScores)
{
int StartLoci;
int TargMatchLen = pAccumScore->TargEndIdx - pAccumScore->TargStartIdx;
if(pAccumScore->Strand == '+')
	StartLoci = pAccumScore->TargStartIdx;
else
	StartLoci = pAccumScore->TargLen - pAccumScore->TargEndIdx;
return(AddMatch(pAccumScore->pCurProbe,pAccumScore->Strand,pAccumScore->TargEntryID,StartLoci,TargMatchLen,NormScore));
}

// FindSWAlignment
// Find alignment using optimised smith-waterman
int												// < 0 if errors, otherwise total number of hits
FindSWAlignment(tsAccumWalkScore *pAccumScore,	// current accmulative walk score
			  tsScores *pScores)
{
int Rslt;
int TotHits;
etSeqBase ProbeBase;
etSeqBase TargBase;
etSeqBase *pProbe;
etSeqBase *pTarg;
tsAccumWalkScore CurAccumWalkScore;
tsAccumWalkScore TmpAccumWalkScore;
double NormScore;

memcpy(&CurAccumWalkScore,pAccumScore,sizeof(tsAccumWalkScore));
pProbe = &CurAccumWalkScore.pProbe[CurAccumWalkScore.ProbeIdx];
pTarg = &CurAccumWalkScore.pTarg[CurAccumWalkScore.TargEndIdx];

Rslt = eBSFSuccess;
TotHits = 0;
while((ProbeBase = (*pProbe++ & ~cRptMskFlg))!=eBaseEOS && (TargBase = (*pTarg++ & ~cRptMskFlg))!=eBaseEOS)
	{
	if(ProbeBase > eBaseT || TargBase > eBaseT ||	// treat indeterminate bases as mismatches
		TargBase != ProbeBase)	// if mismatch then could be because either a true mismatch, or InDel
		{
		if((!CurAccumWalkScore.CurInDelLen && CurAccumWalkScore.TotInDels < pScores->ThresAccumInDels) &&
			CurAccumWalkScore.CurInDelLen < pScores->MaxInDelLen && CurAccumWalkScore.TotInDelLen < pScores->ThresAccumInDelLen)
			{
			memcpy(&TmpAccumWalkScore,&CurAccumWalkScore,sizeof(tsAccumWalkScore));
			if((TmpAccumWalkScore.CurInDelLen += 1) == 1)	// if this is starting an InDel
				{
				TmpAccumWalkScore.TotInDels += 1;
				TmpAccumWalkScore.AccumScore -= pScores->InDelOpen;
				}
			else											// must be an InDel extension
				TmpAccumWalkScore.AccumScore -= pScores->InDelExtn;

			TmpAccumWalkScore.TotInDelLen += 1;
			TmpAccumWalkScore.ProbeIdx += 1;				// explore as InDel in probe
			if((Rslt=FindSWAlignment(&TmpAccumWalkScore,pScores)) <  eBSFSuccess)
				return(Rslt);
			TotHits += Rslt;

			TmpAccumWalkScore.ProbeIdx -= 1;				// explore as InDel in target
			TmpAccumWalkScore.TargEndIdx += 1;
			if((Rslt=FindSWAlignment(&TmpAccumWalkScore,pScores)) <  eBSFSuccess)
				return(Rslt);
			TotHits += Rslt;
			}

		// finished exploring as InDel, treat as mismatch
		if(CurAccumWalkScore.TotMismatches > pScores->ThresAccumMismatches)
			break;
		
		CurAccumWalkScore.CurInDelLen = 0;				
		CurAccumWalkScore.AccumScore -= pScores->Mismatch; 
		CurAccumWalkScore.TotMismatches += 1;		
		if(CurAccumWalkScore.AccumScore < pScores->MinScore) 
			break;
		}
	else				// else an exact match
		{
		CurAccumWalkScore.AccumScore += pScores->Match;
		CurAccumWalkScore.TotMatches += 1;
		if(CurAccumWalkScore.AccumScore > CurAccumWalkScore.PeakAccumScore)		// note peak score
			CurAccumWalkScore.PeakAccumScore = CurAccumWalkScore.AccumScore;
		}

	CurAccumWalkScore.CurInDelLen = 0;
	CurAccumWalkScore.ProbeIdx += 1;
	CurAccumWalkScore.TargEndIdx += 1;
	}

NormScore = (double)CurAccumWalkScore.PeakAccumScore / CurAccumWalkScore.ProbeLen;
if(NormScore >= pScores->ThresAccumScore &&
	CurAccumWalkScore.TotMismatches <= pScores->ThresAccumMismatches &&
	CurAccumWalkScore.TotInDels <= pScores->ThresAccumInDels &&
	CurAccumWalkScore.TotInDelLen <= pScores->ThresAccumInDelLen)
	if((Rslt=ProcessAlignment(NormScore,&CurAccumWalkScore,pScores))>=eBSFSuccess)
		TotHits += 1;
return(Rslt >= eBSFSuccess ? TotHits : Rslt);
}


