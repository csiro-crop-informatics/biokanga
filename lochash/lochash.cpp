// lochash.cpp : Defines the entry point for the console application
// Processes probes against targets using hashes
// Locates unanchored matches in each probe with at most the specified number of mismatches at a specified minimum length, and then extends
// left and right until there are at most the specified number of mismatches.
// Probes and targets can be either either multifasta *.fa) or bioseq (.seq)
// 
// 
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "stdafx.h"

#if _WIN32

#include "../conservlib/commhdrs.h"

#else
#include "../libbiokanga/commhdrs.h"
#endif

const unsigned int cProgVer = 303;		// increment with each release

CStopWatch gStopWatch;
CDiagnostics gDiagnostics;				// for writing diagnostics messages to log file
char gszProcName[_MAX_FNAME];			// process name

const int cMaxFSAMProbeLen = 10000;		   // maximum allowed probe length
const int cMaxHashLen = 12;				   // max number of bases over which to calculate hashes
const int cAllocMatchLoci = 100000;		   // allocate match loci blocks holding this many match loci

const int cMinProbeLen = 8;				    // minimum length probes
const int cDfltMinProbeLen = 16;			// default minimum length probe
const int cMaxProbeLen = 10000;				// maximum length probes

typedef enum TAG_eProcMode {				// processing modes
	ePMGlobalCnts = 0,
	ePMChromCnts,
	ePMHitLoci
	} teProcMode;

typedef struct TAG_sFSAMProbe {
	struct TAG_sFSAMProbe *pNext;			// probes are linked by probe identifier
	struct TAG_sFSAMProbe *pNextPlusHashed;	// and by HashPlusIdx for exacts and 1 mismatches allowed
	struct TAG_sFSAMProbe *pNextMinHashed;	// and by HashMinIdx

	struct TAG_sFSAMProbe *pNextPlus2Hashed;// and by HashPlusIdx for 1 mismatches allowed
	struct TAG_sFSAMProbe *pNextMin2Hashed;	// and by HashMinIdx

	unsigned int HashPlusIdx;				// hash over plus strand initial bases of sequence
	unsigned int HashMinIdx;				// hash over minus strand initial bases of sequence

	unsigned int HashPlus2Idx;				// hash over plus strand initial bases of sequence
	unsigned int HashMin2Idx;				// hash over minus strand initial bases of sequence


	int ProbID;								// uniquely identifies this probe
	char szDescr[80];						// as parsed from fasta descriptor

	tBSFEntryID PrvEntryID;					// used to check if previously hit detected to same chromosome
	int PrvLociHit;							// used to check if previously hit detected to same loci
	int TotPlusHits;						// total hits on '+' strand
	int TotMinHits;							// total hits on '-' strand
	int ProbeLen;							// number of bases in Bases[]
	etSeqBase Bases[1];						// to hold probe sequence
} sFSAMProbe;


// holds match loci 
typedef struct TAG_sMatchLoci {
	sFSAMProbe *pProbe;					// probe which matched
	char Strand;						// was on to this strand
	tBSFEntryID EntryID;				// matches on to this chromosome/sequence
	int MatchPsn;						// and at this loci
} sMatchLoci;

typedef struct TAG_sMatchLociBlock {
	struct TAG_sMatchLociBlock *pNext;	// pts to next allocated block of match loci
	int NumLoci;						// number of loci currently in this block
	sMatchLoci Loci[cAllocMatchLoci];
	} sMatchLociBlock;

// holds approximate matches as linked list prior to sorting when deduping
typedef struct TAG_sFSAMMatch {
	struct TAG_sFSAMMatch *pNext;		// matches are linked
	int MatchID;						// uniquely identifies this match
	int MatchLen;						// number of bases in Bases[]
	etSeqBase Bases[1];					// to hold match sequence
} sFSAMMatch;

int gNumProbes = 0;					// number of probes loaded
int gMaxProbeLen = 0;				// max length of any probe
int gMinProbeLen = 0;				// min length of any probe
sFSAMProbe *gpProbes = NULL;		// pts to linked list of sAMProbes

int	gCurHashLen;					// hash length in use
int gCurHashMsk;					// mask to keep hash in range
int gCurHashEls;					// number of hash elements 

sFSAMProbe **gppPlusHashes = NULL;	// pts to allocated array of ptrs to sAMProbes on plus strand indexed by HashIdx

sFSAMProbe **gppMinHashes = NULL;	// pts to allocated array of ptrs to sAMProbes on minus strand indexed by HashIdx

sFSAMProbe **gppPlus2Hashes = NULL;	// pts to allocated array of ptrs to sAMProbes on plus strand indexed by HashIdx - used if 1 mismatch allowed
sFSAMProbe **gppMin2Hashes = NULL;	// pts to allocated array of ptrs to sAMProbes on minus strand indexed by HashIdx - used if 1 mismatch allowed

int gMatchCnt;							// global match count
sMatchLociBlock *gpMatchBlocks = NULL;	// pts to linked list of match blocks

int Process(teProcMode Mode,char *pszInputFile,char *pszProbeFile,char *pszOutputFile,char *pszChroms,bool bMisses,bool bPerChrom,char Strand,int iMaxMismatches,bool bBED,bool bBedCompat,int MinRange,int MaxRange,int MinLength, int MaxLength);
int ReportMatches(CBioSeqFile *pBioSeq,int hRsltsFile,int MinCounts,int MaxCounts,bool bBED,bool bBedCompat);
int ReportCounts(int hOutFile,int MinRange,int MaxRange,char *pszChrom);
int ReportMisses(CBioSeqFile *pBioSeq,int hRsltsFile,char *pszChrom);

sMatchLoci *AddMatch(sFSAMProbe *pProbe,char Strand,tBSFEntryID EntryID,int Loci);
void DeleteMatches(void);
int LoadProbes(char *pszProbesFile,int MaxMismatches,int MinLenght, int MaxLength);
void DeleteProbes(void);
int ExactMatch(teProcMode Mode,int hRsltsFile,tBSFEntryID CurEntry,unsigned char *pChromSeqBuff,int ChromSeqLen,char *pszChrom, char Strand,bool bBED,int MinCounts, int MaxCounts,bool bMisses); 
int ApproxMatch(teProcMode Mode,int hRsltsFile,tBSFEntryID CurEntry,unsigned char *pChromSeqBuff,int ChromSeqLen,char *pszChrom, char Strand,int MaxMismatches,bool bBED,bool bMisses);
int	CleanNameList(char *pszTxt);
bool NameInList(char *pszNameList, char *pszName);


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
int iMaxMismatches;			// maximum mismatches allowed
char cStrand;				// which strand '+' or '-' or '*'
int iMode;					// processing mode: 0 global counts, 1 chrom counts, 2 hit loci
bool bBED;					// generate loci (mode=2) as BED instead of default CSV
int iMinRange;				// filter out scores of less than this value
int iMaxRange;				// filter out scores of more than this value
bool bBedCompat;			// limit scores to 0...999 to maintain UCSC BED score compatiblity
char szChromsList[2048];	// CSV list of chromosomes to process or '\0' for all
bool bPerChrom;				// Range filtering is per chromosome (default: over all chromosomes)
bool bMisses;				// output probes which are not matched

int iMinLength;				// Filter out probes of less than this length (default: 16, min: 8)");
int iMaxLength;				// Filter out probes of longer than this length (default: 10000, max: 10000)");

char szOutputFile[_MAX_PATH];
char szInputFile[_MAX_PATH];
char szProbesFile[_MAX_PATH];

// command line args
struct arg_lit  *help    = arg_lit0("hH","help",                "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_int *ScreenLogLevel=arg_int0("S", "ScreenLogLevel",	"<int>","Level of diagnostics written to screen 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"diagnostics log file");

struct arg_int *Mode=arg_int0("M", "mode","<int>","Results mode: 0=global counts, 1=chrom counts, 2=probe hit loci");
struct arg_file *InFile = arg_file1("i","inindex","<file>",		"input from bioseq (.seq) target file");
struct arg_file *ProbesFile = arg_file1("I","inprobes","<file>","input from multifasta (.fa) probe file");
struct arg_file *OutFile = arg_file1("o","outmatches","<file>",	"output (misses if '-x')matches to this file as CSV");
struct arg_int *MaxMismatches=arg_int0("m", "maxmissmatches","<int>","Maximum number of mismatches (0..1)");
struct arg_str *Strand=arg_str0("s", "strand",	"<string>",		"Strand '+' or '-', default is '*' for both");
struct arg_lit  *Bed    = arg_lit0("b","bed",					"generate loci (mode=2) as BED instead of default CSV");
struct arg_lit  *BedCompat    = arg_lit0("B","bedcompat",		"limit scores to 0...999 to maintain UCSC BED score compatiblity");

struct arg_int *MinLength=arg_int0("l", "minrange","<int>","Filter out probes of less than this length (default: 16, min: 8)");
struct arg_int *MaxLength=arg_int0("L", "maxrange","<int>","Filter out probes of longer than this length (default: 10000, max: 10000)");


struct arg_int *MinRange=arg_int0("r", "minrange","<int>","Filter out hit counts of less than this value (default: 1)");
struct arg_int *MaxRange=arg_int0("R", "maxrange","<int>","Filter out hit counts of more than this value (default: 999)");

struct arg_str *Chroms=arg_strn("c", "chroms",	"<string>",	0,50,	"Comma/space separated list of chromosomes to process (default: all)");
struct arg_lit *PerChrom=arg_lit0("C", "chromorange",			"Range filtering is per chromosome (default: over all chromosomes)");

struct arg_lit  *Misses  = arg_lit0("x","misses",				"output probes which are not matched");

struct arg_end *end = arg_end(20);

void *argtable[] = {help,version,FileLogLevel,ScreenLogLevel,LogFile,
					Mode,InFile,OutFile,ProbesFile,Misses,
					MaxMismatches,
					Strand,Bed,BedCompat,
					MinRange,MaxRange,Chroms,PerChrom,
					MinLength,MaxLength,
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
	if(MaxMismatches->count)
		{
		iMaxMismatches = MaxMismatches->ival[0];
		if(iMaxMismatches < 0)
			iMaxMismatches = 0;
		else
			if(iMaxMismatches > 1)
				iMaxMismatches = 1;
		}
	else
		iMaxMismatches = 0;

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

	bMisses = Misses->count ? true : false;
	if(Mode->count)
		{
		iMode = Mode->ival[0];
		if(iMode < ePMGlobalCnts)
			iMode = ePMGlobalCnts;
		else
			if(iMode > ePMHitLoci)
				iMode = ePMHitLoci;
		}
	else
		iMode = ePMGlobalCnts;
	if(bMisses && iMode == ePMHitLoci)
		{
		printf("\nError: Can't process hit loci ('-m2') at same time as reporting missing probes ('-x')\n");
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
		printf("\nError: Maximum length probe '-l%d' not in range %d..%d\n",iMaxLength,iMinLength,cMaxProbeLen);
		exit(1);
		}


	if(!bMisses)
		{
		if(MinRange->count)
			{
			iMinRange = MinRange->ival[0];
			if(iMinRange < 0)
				iMinRange = 0;
			else
				if(iMinRange > 10000000)
					iMinRange = 10000000;
			}
		else
			iMinRange = 1;

		if(MaxRange->count)
			{
			iMaxRange = MaxRange->ival[0];
			if(iMaxRange < iMinRange)
				iMaxRange = iMinRange;
			else
				if(iMaxRange > 10000000)
					iMaxRange = 10000000;
			}
		else
			iMaxRange = iMinRange < 999 ? 999 : iMinRange;
		}
	else
		{
		iMinRange = 1;
		iMaxRange = 1;
		}

	if(!bMisses)
		bBedCompat = BedCompat->count ? true : false;
	else
		bBedCompat = false;

	bPerChrom = PerChrom->count ? true : false;

	szChromsList[0] = '\0';
	if(Chroms->count)
		{
		for(int Idx = 0; Idx < Chroms->count; Idx++)
			{
			if(Idx)
				strcat(szChromsList,",");
			strcat(szChromsList,Chroms->sval[Idx]);
			}
		CleanNameList(szChromsList);
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

	if(iMode == ePMHitLoci)
		bBED = Bed->count ? true : false;
	else
		bBED = false;

	if(InFile->count)
		strcpy(szInputFile,InFile->filename[0]);
	else
		strcpy(szInputFile,"in.fmi");


	if(OutFile->count)
		strcpy(szOutputFile,OutFile->filename[0]);
	else
		strcpy(szOutputFile,bBED ? "out.bed" : "out.csv");

	if(ProbesFile->count)
		strcpy(szProbesFile,ProbesFile->filename[0]);
	else
		strcpy(szProbesFile,"probes.fa");

		// now that command parameters have been parsed then initialise diagnostics log system
	if(!gDiagnostics.Open(szLogFile,(etDiagLevel)iScreenLogLevel,(etDiagLevel)iFileLogLevel,true))
		{
		printf("\nError: Unable to start diagnostics subsystem.");
		if(szLogFile[0] != '\0')
			printf(" Most likely cause is that logfile '%s' can't be opened/created",szLogFile);
		exit(1);
		}

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Version: %d.%2.2d Processing parameters:",cProgVer/100,cProgVer%100);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Input from bioseq sequence file: %s",szInputFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Probes from multifasta file: %s",szProbesFile);
	if(!bMisses)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Output as %s file: '%s'",bBED ? "BED" : "CSV", szOutputFile);
	else
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Output as non matched probes file: '%s'",szOutputFile);

	if(szChromsList[0] != '\0')
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Only process these chromosomes: %s",szChromsList);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Probes from multifasta file: %s",szProbesFile);
	switch(iMode) {
		case ePMGlobalCnts:
			gDiagnostics.DiagOutMsgOnly(eDLInfo,(char *)(bMisses ? "non matched probes per genome will be generated" : "Hit counts per genome will be generated"));
			break;
		case ePMChromCnts:
			gDiagnostics.DiagOutMsgOnly(eDLInfo,(char *)(bMisses ? "non matched probes per chromosome will be generated" : "Hit counts per chromosome will be generated"));
			break;
		case ePMHitLoci:
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"All filtered hit loci will be generated");
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Scores %s",bBedCompat ? "are capped at 999" : "No cap on scores");
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Range filtering is %s",bPerChrom ?  "per chrom" : "over whole genome");
			break;
		}

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Filter out probes with length less than: %d",iMinLength);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Filter out probes with length longer than: %d",iMaxLength);

	if(!bMisses)
		{
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Filter out probes with less than %d hits",iMinRange);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Filter out probes with more than %d hits",iMaxRange);
		}

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Target strand: '%c'",cStrand);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Max mismatches: %d",iMaxMismatches);
#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	gStopWatch.Start();
	Rslt = Process((teProcMode)iMode,szInputFile,szProbesFile,szOutputFile,szChromsList,bMisses,bPerChrom,cStrand,iMaxMismatches,bBED,bBedCompat,iMinRange,iMaxRange,iMinLength,iMaxLength);
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


int 
Process(teProcMode Mode,char *pszInputFile,char *pszProbesFile,char *pszOutputFile,char *pszChroms,bool bMisses, bool bPerChrom,char Strand,int MaxMismatches,bool bBED,bool bBedCompat,int MinRange,int MaxRange,int MinLength,int MaxLength)
{
CBioSeqFile *pBioSeq;
unsigned char *pChromSeqBuff;			// to hold sequence for each chromsome
int AllocdChromSeqLen;					// length of currently allocated pChromSeqBuff
int hOutFile;
char szChromName[cBSFSourceSize+1];		// to hold current chromosome name
tBSFEntryID CurEntry;					// current entry being processed
int ChromSeqLen;						// length of chromosome sequence being processed
int Rslt;								// used to hold processing result

if((Rslt=LoadProbes(pszProbesFile,MaxMismatches,MinLength,MaxLength))<=eBSFSuccess)
	return(Rslt);

// open bioseq file containing chromosome sequences
if((pBioSeq = new CBioSeqFile()) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CBioSeqFile");
	DeleteProbes();
	return(eBSFerrObj);
	}
if((Rslt=pBioSeq->Open(pszInputFile,cBSFTypeSeq))!=eBSFSuccess)
	{
	while(pBioSeq->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,pBioSeq->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open input bioseq sequence file '%s'",pszInputFile);
	DeleteProbes();
	delete pBioSeq;
	return(Rslt);
	}

// open and truncate CSV output file
#ifdef _WIN32
hOutFile = open(pszOutputFile, _O_WRONLY | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE );
#else
hOutFile = open(pszOutputFile, O_WRONLY | O_CREAT | O_TRUNC, S_IREAD | S_IWRITE );
#endif
if(hOutFile == -1)
	{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open output %s file '%s', error: %s",bBED ? "BED" : "CSV", pszOutputFile,strerror(errno));
	DeleteProbes();
	delete pBioSeq;
	return(eBSFerrOpnFile);
	}

	// iterate over each entry (chromosome) in sequence file
pChromSeqBuff = NULL;
ChromSeqLen = 0;
AllocdChromSeqLen = 0;
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
	ChromSeqLen = pBioSeq->GetDataLen(CurEntry);
	if(pChromSeqBuff == NULL || ChromSeqLen > AllocdChromSeqLen)
		{
		if(pChromSeqBuff != NULL)
			delete pChromSeqBuff;
		AllocdChromSeqLen = ChromSeqLen + ChromSeqLen/10; // alloc more than actually required, next sequence may not need to be realloc'd
		pChromSeqBuff = (unsigned char *)new unsigned char [AllocdChromSeqLen];
		if(pChromSeqBuff == NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate memory (%d bytes) for holding chromosome sequence",AllocdChromSeqLen);
			Rslt = eBSFerrMem;
			break;
			}
		}

	pBioSeq->GetData(CurEntry,cBSFTypeSeq,0,pChromSeqBuff,ChromSeqLen);
	if(MaxMismatches == 0)
		ExactMatch(Mode,hOutFile,CurEntry,pChromSeqBuff,ChromSeqLen,szChromName,Strand,bBED,MinRange,MaxRange,bMisses); 
	else
		ApproxMatch(Mode,hOutFile,CurEntry,pChromSeqBuff,ChromSeqLen,szChromName,Strand,MaxMismatches,bBED,bMisses);
	
	if(Mode == ePMChromCnts)
		{
		if(bMisses)
			ReportMisses(pBioSeq,hOutFile,szChromName);
		else
			ReportCounts(hOutFile,MinRange,MaxRange,szChromName);
		}

	if(Mode == ePMHitLoci && bPerChrom)
		{
		ReportMatches(pBioSeq,hOutFile,MinRange,MaxRange,bBED,bBedCompat);
		DeleteMatches();
		}
	}

if(Mode == ePMGlobalCnts)
	{
	gMatchCnt = 0;
	if(bMisses)
		ReportMisses(pBioSeq,hOutFile,"ALLCHROMES");
	else
		ReportCounts(hOutFile,MinRange,MaxRange,"ALLCHROMES");
	}
else
	if(Mode == ePMHitLoci &&  !bPerChrom)
		ReportMatches(pBioSeq,hOutFile,MinRange,MaxRange,bBED,bBedCompat);
DeleteMatches();
DeleteProbes();
if(hOutFile!=-1)
	close(hOutFile);

if(pChromSeqBuff!=NULL)
	delete pChromSeqBuff;
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


sMatchLoci *
AddMatch(sFSAMProbe *pProbe,char Strand,tBSFEntryID EntryID,int Loci)
{
sMatchLoci *pMatch;
sMatchLociBlock *pNewBlock;

if(gpMatchBlocks == NULL ||
   gpMatchBlocks->NumLoci == cAllocMatchLoci)
	{
	if((pNewBlock = new sMatchLociBlock)==NULL)
		return(NULL);
	pNewBlock->NumLoci = 0;
	pNewBlock->pNext = gpMatchBlocks;
	gpMatchBlocks = pNewBlock;
	}
pMatch = &gpMatchBlocks->Loci[gpMatchBlocks->NumLoci++];
pMatch->Strand = Strand;
pMatch->EntryID = EntryID;
pMatch->MatchPsn = Loci;
pMatch->pProbe = pProbe;
gMatchCnt+=1;
return(pMatch);
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
if(gppPlusHashes != NULL)
	delete gppPlusHashes;
gppPlusHashes = NULL;
if(gppMinHashes != NULL)
	delete gppMinHashes;
gppMinHashes = NULL;
if(gppPlus2Hashes != NULL)
	delete gppPlus2Hashes;
if(gppMin2Hashes != NULL)
	delete gppMinHashes;
gppMin2Hashes = NULL;
}



// Generate hash over Len (max gCurHashLen) pBases
int 
GenPlusHash(etSeqBase *pBases)		// assumed to pt to first base
{
int HashIdx = 0;
int Len = gCurHashLen;
while(Len--)
	{
	HashIdx <<= 2;
	HashIdx |= (*pBases++ & 0x03);
	}
return(HashIdx & gCurHashMsk);
}

// Generate 1 mismatch hash over Len (max gCurHashLen) pBases
// Every 2nd base starting from 1st (0,2,4..) is used in hash
int 
GenPlus1Hash(etSeqBase *pBases)		// assumed to pt to first base
{
int HashIdx = 0;
int Len = (gCurHashLen+1)/2;
while(Len--)
	{
	HashIdx <<= 2;
	HashIdx |= (*pBases & 0x03);
	pBases += 2;
	}
return(HashIdx & gCurHashMsk);
}

// Generate 1 mismatch hash over Len (max gCurHashLen) pBases
// Every 2nd base starting from 2nd (1,3,5..) is used in hash
int 
GenPlus2Hash(etSeqBase *pBases)		// assumed to pt to first base
{
int HashIdx = 0;
int Len = gCurHashLen/2;
pBases++;
while(Len--)
	{
	HashIdx <<= 2;
	HashIdx |= (*pBases & 0x03);
	pBases += 2;
	}
return(HashIdx & gCurHashMsk);
}


// Generate hash over Len (max gCurHashLen) pBases reverse complemented
int 
GenMinHash(etSeqBase *pBases)		// assumed to pt to last base
{
int HashIdx = 0;
etSeqBase CplBase;
int	Len = gCurHashLen;
while(Len--)
	{
	HashIdx <<= 2;
	CplBase = *pBases-- & 0x03;
	switch(CplBase) {
		case eBaseA:
			CplBase = eBaseT;
			break;
		case eBaseC:
			CplBase = eBaseG;
			break;
		case eBaseG:
			CplBase = eBaseC;
			break;
		case eBaseT:
			CplBase = eBaseA;
			break;
		}
	HashIdx |= CplBase;
	}
return(HashIdx & gCurHashMsk);
}


// Generate hash over Len (max gCurHashLen) pBases reverse complemented
// Every 2nd base starting from last is used in hash
int 
GenMin1Hash(etSeqBase *pBases)		// assumed to pt to last base
{
int HashIdx = 0;
etSeqBase CplBase;
int	Len = (gCurHashLen+1)/2;
while(Len--)
	{
	HashIdx <<= 2;
	CplBase = *pBases & 0x03;
	pBases-=2;
	switch(CplBase) {
		case eBaseA:
			CplBase = eBaseT;
			break;
		case eBaseC:
			CplBase = eBaseG;
			break;
		case eBaseG:
			CplBase = eBaseC;
			break;
		case eBaseT:
			CplBase = eBaseA;
			break;
		}
	HashIdx |= CplBase;
	}
return(HashIdx & gCurHashMsk);
}


// Generate hash over Len (max gCurHashLen) pBases reverse complemented
// Every 2nd base starting from last-1 is used in hash
int 
GenMin2Hash(etSeqBase *pBases)		// assumed to pt to last base
{
int HashIdx = 0;
etSeqBase CplBase;
int	Len = gCurHashLen/2;
pBases-=1;
while(Len--)
	{
	HashIdx <<= 2;
	CplBase = *pBases & 0x03;
	pBases-=2;
	switch(CplBase) {
		case eBaseA:
			CplBase = eBaseT;
			break;
		case eBaseC:
			CplBase = eBaseG;
			break;
		case eBaseG:
			CplBase = eBaseC;
			break;
		case eBaseT:
			CplBase = eBaseA;
			break;
		}
	HashIdx |= CplBase;
	}
return(HashIdx & gCurHashMsk);
}


int
LoadProbes(char *pszProbesFile,int MaxMismatches,int MinLength,int MaxLength)
{
int Rslt;
int UnderLength;
int OverLength;
sFSAMProbe *pCurProbe;
sFSAMProbe *pNxtProbe;
CFasta *pFasta;

char szDescription[cBSFDescriptionSize];
etSeqBase SeqBuff[cMaxFSAMProbeLen];
int SeqLen;
bool bDescriptor;
bool bFirst;
int TruncCnt;


gNumProbes = 0;
gpProbes = NULL;
if(gppPlusHashes)
	delete gppPlusHashes;
gppPlusHashes = NULL;
if(gppMinHashes)
	delete gppMinHashes;
gppMinHashes = NULL;

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

// read each probe sequence into memory array
// note that only the first cMaxFSAMProbeLen of each probe is processed
bDescriptor = false;
bFirst = true;
TruncCnt = 0;
UnderLength = 0;
OverLength = 0;
while((Rslt = SeqLen = pFasta->ReadSequence(SeqBuff,cMaxFSAMProbeLen)) > eBSFSuccess)
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
			if(!bFirst)						// only process 1st cMaxFSAMProbeLen
				{
				if(!TruncCnt++)
					gDiagnostics.DiagOut(eDLInfo,gszProcName,"LoadProbes [%s] truncated probe %s at %d length",pszProbesFile,szDescription, cMaxFSAMProbeLen);
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

	if(gMaxProbeLen == 0 || SeqLen > gMaxProbeLen)
		gMaxProbeLen = SeqLen;
	if(gMinProbeLen == 0 || SeqLen < gMinProbeLen)
			gMinProbeLen = SeqLen;
	pNxtProbe = (sFSAMProbe *)new unsigned char [ sizeof(sFSAMProbe) + SeqLen];
	pNxtProbe->pNext = NULL;
	pNxtProbe->pNextPlusHashed = NULL;
	pNxtProbe->pNextMinHashed = NULL;
	pNxtProbe->pNextPlus2Hashed = NULL;
	pNxtProbe->pNextMin2Hashed = NULL;
	pNxtProbe->ProbID = ++gNumProbes;
	pNxtProbe->ProbeLen = SeqLen;
	pNxtProbe->TotMinHits = 0;
	pNxtProbe->TotPlusHits = 0;
	strncpy(pNxtProbe->szDescr,szDescription,sizeof(pNxtProbe->szDescr)-1);
	pNxtProbe->szDescr[sizeof(pNxtProbe->szDescr)-1] = '\0';
	CSeqTrans::RemoveMasking(SeqBuff,SeqLen);
	memcpy(pNxtProbe->Bases,SeqBuff,SeqLen);
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

if(MaxMismatches == 0)
	{
	gCurHashLen = min(gMinProbeLen,cMaxHashLen);
	gCurHashEls = 0x00001 << (gCurHashLen * 2);
	gCurHashMsk = ~(-1 << (gCurHashLen * 2));
	gppPlusHashes = new sFSAMProbe * [gCurHashEls];	// allows for hashes over gCurHashLen bases
	memset(gppPlusHashes,0,sizeof(sFSAMProbe *) * gCurHashEls);

	gppMinHashes = new sFSAMProbe * [gCurHashEls];	// allows for hashes over gCurHashLen bases
	memset(gppMinHashes,0,sizeof(sFSAMProbe *) * gCurHashEls);

	pCurProbe = gpProbes;
	// construct hashes over probes
	while(pCurProbe != NULL)
		{
		pCurProbe->HashPlusIdx = GenPlusHash(pCurProbe->Bases);

		if(gppPlusHashes[pCurProbe->HashPlusIdx] != NULL)
			pCurProbe->pNextPlusHashed = gppPlusHashes[pCurProbe->HashPlusIdx];
		else
			pCurProbe->pNextPlusHashed = NULL;
		gppPlusHashes[pCurProbe->HashPlusIdx] = pCurProbe;

		pCurProbe->HashMinIdx  = GenMinHash(&pCurProbe->Bases[pCurProbe->ProbeLen-1]);
		if(gppMinHashes[pCurProbe->HashMinIdx] != NULL)
			pCurProbe->pNextMinHashed = gppMinHashes[pCurProbe->HashMinIdx];
		else
			pCurProbe->pNextMinHashed = NULL;
		gppMinHashes[pCurProbe->HashMinIdx] = pCurProbe;

		pCurProbe = pCurProbe->pNext;
		}
	}
else
	{
	gCurHashLen = min(gMinProbeLen,cMaxHashLen * 2);
	gCurHashEls = 0x00001 << (gCurHashLen);
	gCurHashMsk = ~(-1 << (gCurHashLen));
	gppPlusHashes = new sFSAMProbe * [gCurHashEls];	// allows for hashes over gCurHashLen bases
	memset(gppPlusHashes,0,sizeof(sFSAMProbe *) * gCurHashEls);

	gppPlus2Hashes = new sFSAMProbe * [gCurHashEls];	// allows for hashes over gCurHashLen bases
	memset(gppPlus2Hashes,0,sizeof(sFSAMProbe *) * gCurHashEls);

	gppMinHashes = new sFSAMProbe * [gCurHashEls];	// allows for hashes over gCurHashLen bases
	memset(gppMinHashes,0,sizeof(sFSAMProbe *) * gCurHashEls);

	gppMin2Hashes = new sFSAMProbe * [gCurHashEls];	// allows for hashes over gCurHashLen bases
	memset(gppMin2Hashes,0,sizeof(sFSAMProbe *) * gCurHashEls);

	pCurProbe = gpProbes;

	// construct hashes over probes
	while(pCurProbe != NULL)
		{
		pCurProbe->HashPlusIdx = GenPlus1Hash(pCurProbe->Bases);
		if(gppPlusHashes[pCurProbe->HashPlusIdx] != NULL)
			pCurProbe->pNextPlusHashed = gppPlusHashes[pCurProbe->HashPlusIdx];
		else
			pCurProbe->pNextPlusHashed = NULL;
		gppPlusHashes[pCurProbe->HashPlusIdx] = pCurProbe;
		pCurProbe->HashPlus2Idx = GenPlus2Hash(pCurProbe->Bases);
		if(gppPlus2Hashes[pCurProbe->HashPlus2Idx] != NULL)
			pCurProbe->pNextPlus2Hashed = gppPlus2Hashes[pCurProbe->HashPlus2Idx];
		else
			pCurProbe->pNextPlus2Hashed = NULL;
		gppPlus2Hashes[pCurProbe->HashPlus2Idx] = pCurProbe;


		pCurProbe->HashMinIdx  = GenMin1Hash(&pCurProbe->Bases[pCurProbe->ProbeLen-1]);
		if(gppMinHashes[pCurProbe->HashMinIdx] != NULL)
			pCurProbe->pNextMinHashed = gppMinHashes[pCurProbe->HashMinIdx];
		else
			pCurProbe->pNextMinHashed = NULL;
		gppMinHashes[pCurProbe->HashMinIdx] = pCurProbe;

		pCurProbe->HashMin2Idx  = GenMin2Hash(&pCurProbe->Bases[pCurProbe->ProbeLen-1]);
		if(gppMin2Hashes[pCurProbe->HashMin2Idx] != NULL)
			pCurProbe->pNextMin2Hashed = gppMin2Hashes[pCurProbe->HashMin2Idx];
		else
			pCurProbe->pNextMin2Hashed = NULL;
		gppMin2Hashes[pCurProbe->HashMin2Idx] = pCurProbe;

		pCurProbe = pCurProbe->pNext;
		}
	}

if(UnderLength)
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Filtered out %d underlength probes",UnderLength);
if(OverLength)
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Filtered out %d overlength probes",OverLength);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Number of probes to process %d",gNumProbes);
return(gNumProbes);
}


int
ReportCounts(int hOutFile,int MinRange,int MaxRange,char *pszChrom)
{
int TotHits;
int RsltsLen;
char szOutBuff[0x3fff];
sFSAMProbe *pCurProbe = gpProbes;
while(pCurProbe != NULL)
	{
	TotHits = pCurProbe->TotPlusHits + pCurProbe->TotMinHits;
	if(TotHits >= MinRange && TotHits <= MaxRange)
		{
		RsltsLen = sprintf((char *)szOutBuff,"%d,%d,\"%s\",%d,%d,\"%s\"\n",++gMatchCnt,pCurProbe->ProbID,pCurProbe->szDescr,pCurProbe->TotPlusHits,pCurProbe->TotMinHits,pszChrom);
		write(hOutFile,szOutBuff,RsltsLen);
		}
	pCurProbe->TotPlusHits = 0;
	pCurProbe->TotMinHits = 0;
	pCurProbe = pCurProbe->pNext;
	}
return(gMatchCnt);
}

// ReportMisses
// Reports all probes which had 0 matches
int 
ReportMisses(CBioSeqFile *pBioSeq,int hRsltsFile,char *pszChrom)
{
unsigned char szRsltsBuff[0x7fff];		// used to hold results before writing to disk
unsigned char szSeqsBuff[0x7fff];		// used to hold results before writing to disk

int	RsltsLen;							// current length of results in buffer
char *pszAsciiBases;
sFSAMProbe *pCurProbe = gpProbes;
while(pCurProbe != NULL)
	{
	if(!(pCurProbe->TotPlusHits + pCurProbe->TotMinHits))
		{
		pszAsciiBases = CSeqTrans::MapSeq2Ascii(pCurProbe->Bases,pCurProbe->ProbeLen,(char *)szSeqsBuff);
		RsltsLen = sprintf((char *)szRsltsBuff,"%d,\"%s\",%d,\"%s\",\"%s\"\n",++gMatchCnt,pszChrom,pCurProbe->ProbID,pCurProbe->szDescr,pszAsciiBases);
		write(hRsltsFile,szRsltsBuff,RsltsLen);
		}
	else
		{
		pCurProbe->TotPlusHits = 0;
		pCurProbe->TotMinHits = 0;
		}
	pCurProbe = pCurProbe->pNext;
	}
return(gMatchCnt);
}

int ReportMatches(CBioSeqFile *pBioSeq,int hRsltsFile,int MinCounts,int MaxCounts,bool bBED,bool bBedCompat)
{
sMatchLociBlock *pCurBlock;
sMatchLoci *pCurMatch;
sFSAMProbe *pProbe;
char szChromName[cBSFSourceSize+1];		// to hold cached chromosome name
tBSFEntryID CurEntryID;					// entry cached

unsigned char szRsltsBuff[0x7fff];		// used to hold results before writing to disk
int	RsltsLen;							// current length of results in buffer
char *pszAsciiBases;
int MatchCnt;
int TotHits;
int LociIdx;

if((pCurBlock = gpMatchBlocks) == NULL)
	return(0);
CurEntryID = 0;
MatchCnt = 0;
do {
	pCurMatch = &pCurBlock->Loci[0];
	for(LociIdx = 0; LociIdx < pCurBlock->NumLoci; LociIdx++,pCurMatch++)
		{
		pProbe = pCurMatch->pProbe;
		TotHits = pProbe->TotPlusHits + pProbe->TotMinHits;
		if(TotHits < MinCounts ||  TotHits > MaxCounts)
			continue;
		
		if(CurEntryID == 0 || pCurMatch->EntryID != CurEntryID)
			{
			pBioSeq->GetName(pCurMatch->EntryID,sizeof(szChromName),szChromName);
			CurEntryID = pCurMatch->EntryID;
			}

		if(bBED && TotHits > 0)
			RsltsLen = sprintf((char *)szRsltsBuff,"%s\t%d\t%d\t%s\t%d\t%c\n",szChromName,pCurMatch->MatchPsn,pCurMatch->MatchPsn+pProbe->ProbeLen,pProbe->szDescr,(bBedCompat && TotHits > 999) ? 999 : TotHits, pCurMatch->Strand);
		else
			{
			pszAsciiBases = CSeqTrans::MapSeq2Ascii(pProbe->Bases,pProbe->ProbeLen);
			RsltsLen = sprintf((char *)szRsltsBuff,"%d,\"%s\",\"%s\",\"%c\",%lu,%d,\"%s\",%d\n",++MatchCnt,pProbe->szDescr,szChromName,pCurMatch->Strand,pCurMatch->MatchPsn,pProbe->ProbeLen,pszAsciiBases,TotHits);
			}
		write(hRsltsFile,szRsltsBuff,RsltsLen);
		MatchCnt+=1;
		}
	}
while(pCurBlock = pCurBlock->pNext);
return(MatchCnt);
}

int										// returns the number of hits by probes against the chromosome sequence
ExactMatch(teProcMode Mode,int hRsltsFile,tBSFEntryID CurEntry,unsigned char *pChromSeqBuff,int ChromSeqLen,char *pszChrom, char Strand,bool bBED,int MinCounts, int MaxCounts,bool bMisses) 
{
sFSAMProbe *pPutProbe;					// putative matching probe
unsigned int HashIdx;
int NumWindows;
unsigned char *pLeft;
unsigned char *pRight;
unsigned char *pPutProbeSeq;
unsigned char *pPutSeq;
int PutSeqIdx;
int Loci;
int PlusHits;
int MinHits;
int TotHits;

etSeqBase ProbeBase;
etSeqBase SeqBase;

if(ChromSeqLen < gCurHashLen)					// don't bother if chrom seq is smaller than hash seq size
	return(0);

Loci = 0;
NumWindows = ChromSeqLen - gCurHashLen - 1;
pLeft = pChromSeqBuff;
pRight = &pChromSeqBuff[gCurHashLen];

PlusHits = 0;
MinHits = 0;
printf("\n   000000000");

HashIdx = GenPlusHash(pLeft);				// initial hash
while(NumWindows--)
	{
	if(!(Loci & 0x0fffff))
		printf("\b\b\b\b\b\b\b\b\b%9.9d",Loci);		// let user know we haven't crashed...
	if((*pLeft & ~ cRptMskFlg) <= eBaseT)
		{
		// any probes on plus strand with same hash?
		if((Strand == '*' || Strand == '+') && (pPutProbe = gppPlusHashes[HashIdx])!=NULL)
			{
			while(pPutProbe != NULL)
				{
				pPutProbe->PrvEntryID = -1;		// no hit at current loci for this probe
				pPutProbe->PrvLociHit = -1;
				if(pPutProbe->ProbeLen <= (ChromSeqLen - Loci) && (!bMisses || !pPutProbe->TotPlusHits))
					{
					pPutProbeSeq = pPutProbe->Bases;
					pPutSeq = pLeft;
					for(PutSeqIdx = 0; PutSeqIdx < pPutProbe->ProbeLen; PutSeqIdx++,pPutSeq++,pPutProbeSeq++)
						{
						ProbeBase = *pPutProbeSeq & ~ cRptMskFlg;
						SeqBase = *pPutSeq & ~ cRptMskFlg;
						if(SeqBase > eBaseT || ProbeBase > eBaseT)
							break;
						if(SeqBase != ProbeBase)
							break;
						}
					if(PutSeqIdx == pPutProbe->ProbeLen)
						{
						PlusHits += 1;
						pPutProbe->TotPlusHits++;
						TotHits = pPutProbe->TotPlusHits + pPutProbe->TotMinHits;
						pPutProbe->PrvEntryID = CurEntry;
						pPutProbe->PrvLociHit = Loci;
						if(Mode==2 && TotHits <= MaxCounts && TotHits >= MinCounts)		
							AddMatch(pPutProbe,'+',CurEntry,Loci);
						}
					}
				pPutProbe = pPutProbe->pNextPlusHashed;
				}
			}

		// any probes on minus strand with same hash?
		if((Strand == '*' || Strand == '-') &&(pPutProbe = gppMinHashes[HashIdx])!=NULL)
			{
			while(pPutProbe != NULL)
				{
				if(!(pPutProbe->PrvEntryID == CurEntry && pPutProbe->PrvLociHit == Loci) &&	// only one hit per probe per loci allowed 
						pPutProbe->ProbeLen <= (ChromSeqLen - Loci)  && (!bMisses || !pPutProbe->TotMinHits))
					{
					pPutProbeSeq = &pPutProbe->Bases[pPutProbe->ProbeLen-1];
					pPutSeq = pLeft;
					for(PutSeqIdx = 0; PutSeqIdx < pPutProbe->ProbeLen; PutSeqIdx++,pPutSeq++,pPutProbeSeq--)
						{
						ProbeBase = *pPutProbeSeq & ~ cRptMskFlg;
						SeqBase = *pPutSeq & ~ cRptMskFlg;
						if(SeqBase > eBaseT || ProbeBase > eBaseT)
							break;
						switch(ProbeBase) {
							case eBaseA:
								if(SeqBase == eBaseT)
									continue;
								break;
							case eBaseC:
								if(SeqBase == eBaseG)
									continue;
								break;
							case eBaseG:
								if(SeqBase == eBaseC)
									continue;
								break;
							case eBaseT:
								if(SeqBase == eBaseA)
									continue;
								break;
							}					
						break;
						}
					if(PutSeqIdx == pPutProbe->ProbeLen)
						{
						MinHits += 1;
						pPutProbe->TotMinHits++;
						TotHits = pPutProbe->TotPlusHits + pPutProbe->TotMinHits;
						pPutProbe->PrvEntryID = CurEntry;
						pPutProbe->PrvLociHit = Loci;
						
						if(Mode==ePMHitLoci && TotHits <= MaxCounts && TotHits >= MinCounts)
							AddMatch(pPutProbe,'-',CurEntry,Loci);
						}
					}
				pPutProbe = pPutProbe->pNextMinHashed;
				}
			}
		}
		// gen new hash
	HashIdx <<= 2;
	HashIdx &= gCurHashMsk;
	HashIdx |= (*pRight++ & 0x03); 
	pLeft++;
	Loci++;
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Hits on plus strand: %d on minus strand: %d",PlusHits,MinHits);
return(PlusHits + MinHits);
}


// ApproxMatch
// Currently only allows for 1 mismatches
int										// returns the number of hits by probes against the chromosome sequence
ApproxMatch(teProcMode Mode,
			int hRsltsFile,
			tBSFEntryID CurEntry,
			unsigned char *pChromSeqBuff,
			int ChromSeqLen,
			char *pszChrom, 
			char Strand,				// which strand '*', '+' or '-'	
			int MaxMismatches,			// allowed mismatches, currently limited to a single mismatch
			bool bBED,
			bool bMisses)				// true if interested in probes which do not have at least one hit 
{
sFSAMProbe *pPutProbe;					// putative matching probe
unsigned char szRsltsBuff[0x7fff];		// used to hold results before writing to disk
int	RsltsLen;							// current length of results in buffer
char *pszAsciiBases;
unsigned int Hash1Idx;
unsigned int Hash2Idx;
int NumWindows;
unsigned char *pLeft;
unsigned char *pPutProbeSeq;
unsigned char *pPutSeq;
int PutSeqIdx;
int Loci;
int CurMismatches;
int PlusHits;
int MinHits;
etSeqBase ProbeBase;
etSeqBase SeqBase;

if(ChromSeqLen < gCurHashLen)					// don't bother if chrom seq is smaller than hash seq size
	return(0);

Loci = 0;
PlusHits = 0;
MinHits = 0;
NumWindows = ChromSeqLen - gCurHashLen - 1;
pLeft = pChromSeqBuff;
printf("\n   000000000");

while(NumWindows--)
	{
	if(!(Loci & 0x0fffff))
		printf("\b\b\b\b\b\b\b\b\b%9.9d",Loci);		// let user know we haven't crashed...

	Hash1Idx = GenPlus1Hash(pLeft);			// hash window on even bases 0,2,4...
	Hash2Idx = GenPlus2Hash(pLeft);			// hash window on odd bases 1,3,5
	
	// any probes on plus strand with same hash over even bases?
	if((Strand == '*' || Strand == '+') && (pPutProbe = gppPlusHashes[Hash1Idx])!=NULL)
		{
		while(pPutProbe != NULL)	// iterate over all probes with same hash
			{
			pPutProbe->PrvEntryID = -1;		// no hit at current loci for this probe
			pPutProbe->PrvLociHit = -1;
			if(pPutProbe->ProbeLen <= (ChromSeqLen - Loci)  && (!bMisses || !pPutProbe->TotPlusHits))
				{
				pPutProbeSeq = pPutProbe->Bases;
				pPutSeq = pLeft;
				CurMismatches = 0;
				for(PutSeqIdx = 0; PutSeqIdx < pPutProbe->ProbeLen; PutSeqIdx++,pPutSeq++,pPutProbeSeq++)
					{
					ProbeBase = *pPutProbeSeq & ~ cRptMskFlg;
					SeqBase = *pPutSeq & ~ cRptMskFlg;
					if(SeqBase > eBaseT || ProbeBase > eBaseT)	// only interested in a,c,g,t
						break;
					if(SeqBase == ProbeBase)					
						continue;
					if(++CurMismatches > MaxMismatches)			// too many mismatches?
						break;
					}

				if(PutSeqIdx == pPutProbe->ProbeLen)
					{
					pPutProbe->TotPlusHits++;
					pPutProbe->PrvEntryID = CurEntry;
					pPutProbe->PrvLociHit = Loci;
					PlusHits++;
					if(Mode==ePMHitLoci)
						{
						pszAsciiBases = CSeqTrans::MapSeq2Ascii(pPutProbe->Bases,pPutProbe->ProbeLen);
						RsltsLen = sprintf((char *)szRsltsBuff,"%d,\"%s\",\"%s\",\"%c\",%lu,%d,\"%s\"\n",++gMatchCnt,pPutProbe->szDescr,pszChrom,'+',Loci,pPutProbe->ProbeLen,pszAsciiBases);
						write(hRsltsFile,szRsltsBuff,RsltsLen);
						}
					}
				}
			pPutProbe = pPutProbe->pNextPlusHashed;
			}
		}

	// any probes on plus strand with same hash over odd bases?
	if((Strand == '*' || Strand == '+') && (pPutProbe = gppPlus2Hashes[Hash2Idx])!=NULL)
		{
		while(pPutProbe != NULL)
			{
			if(!(pPutProbe->PrvEntryID == CurEntry && pPutProbe->PrvLociHit == Loci) &&	// only one hit per probe per loci allowed 
				pPutProbe->ProbeLen <= (ChromSeqLen - Loci)  && (!bMisses || !pPutProbe->TotPlusHits))
				{
				pPutProbeSeq = pPutProbe->Bases;
				pPutSeq = pLeft;
				CurMismatches = 0;
				for(PutSeqIdx = 0; PutSeqIdx < pPutProbe->ProbeLen; PutSeqIdx++,pPutSeq++,pPutProbeSeq++)
					{
					ProbeBase = *pPutProbeSeq & ~ cRptMskFlg;
					SeqBase = *pPutSeq & ~ cRptMskFlg;
					if(SeqBase > eBaseT || ProbeBase > eBaseT)
						break;
					if(SeqBase == ProbeBase)					
						continue;
					if(++CurMismatches > MaxMismatches)			// too many mismatches?
						break;
					}
				if(PutSeqIdx == pPutProbe->ProbeLen)
					{
					pPutProbe->TotPlusHits++;
					pPutProbe->PrvEntryID = CurEntry;
					pPutProbe->PrvLociHit = Loci;
					PlusHits++;
					if(Mode==ePMHitLoci)
						{
						pszAsciiBases = CSeqTrans::MapSeq2Ascii(pPutProbe->Bases,pPutProbe->ProbeLen);
						RsltsLen = sprintf((char *)szRsltsBuff,"%d,\"%s\",\"%s\",\"%c\",%lu,%d,\"%s\"\n",++gMatchCnt,pPutProbe->szDescr,pszChrom,'+',Loci,pPutProbe->ProbeLen,pszAsciiBases);
						write(hRsltsFile,szRsltsBuff,RsltsLen);
						}
					}
				}
			pPutProbe = pPutProbe->pNextPlus2Hashed;
			}
		}


	// any probes on minus strand with same hash over even bases?
	if((Strand == '*' || Strand == '-') &&(pPutProbe = gppMinHashes[Hash1Idx])!=NULL)
		{
		while(pPutProbe != NULL)
			{
			if(!(pPutProbe->PrvEntryID == CurEntry && pPutProbe->PrvLociHit == Loci) &&	// only one hit per probe per loci allowed 
				pPutProbe->ProbeLen <= (ChromSeqLen - Loci)  && (!bMisses || !pPutProbe->TotMinHits))
				{
				pPutProbeSeq = &pPutProbe->Bases[pPutProbe->ProbeLen-1];
				pPutSeq = pLeft;
				CurMismatches = 0;
				for(PutSeqIdx = 0; PutSeqIdx < pPutProbe->ProbeLen; PutSeqIdx++,pPutSeq++,pPutProbeSeq--)
					{
					ProbeBase = *pPutProbeSeq & ~ cRptMskFlg;
					SeqBase = *pPutSeq & ~ cRptMskFlg;
					if(SeqBase > eBaseT || ProbeBase > eBaseT)
						break;
					switch(ProbeBase) {
						case eBaseA:
							if(SeqBase == eBaseT)
								continue;
							break;
						case eBaseC:
							if(SeqBase == eBaseG)
								continue;
							break;
						case eBaseG:
							if(SeqBase == eBaseC)
								continue;
							break;
						case eBaseT:
							if(SeqBase == eBaseA)
								continue;
							break;
						}					
					if(++CurMismatches > MaxMismatches)			// too many mismatches?
						break;
					}
				if(PutSeqIdx == pPutProbe->ProbeLen)
					{
					pPutProbe->TotMinHits++;
					pPutProbe->PrvEntryID = CurEntry;
					pPutProbe->PrvLociHit = Loci;
					MinHits++;
					if(Mode==ePMHitLoci)
						{
						pszAsciiBases = CSeqTrans::MapSeq2Ascii(pPutProbe->Bases,pPutProbe->ProbeLen);
						RsltsLen = sprintf((char *)szRsltsBuff,"%d,\"%s\",\"%s\",\"%c\",%lu,%d,\"%s\"\n",++gMatchCnt,pPutProbe->szDescr,pszChrom,'-',Loci,pPutProbe->ProbeLen,pszAsciiBases);
						write(hRsltsFile,szRsltsBuff,RsltsLen);
						}
					}
				}
			pPutProbe = pPutProbe->pNextMinHashed;
			}
		}

	// any probes on minus strand with same hash over odd bases?
	if((Strand == '*' || Strand == '-') &&(pPutProbe = gppMin2Hashes[Hash2Idx])!=NULL)
		{
		while(pPutProbe != NULL)
			{
			if(!(pPutProbe->PrvEntryID == CurEntry && pPutProbe->PrvLociHit == Loci) &&	// only one hit per probe per loci allowed 
				pPutProbe->ProbeLen <= (ChromSeqLen - Loci)  && (!bMisses || !pPutProbe->TotMinHits))
				{
				pPutProbeSeq = &pPutProbe->Bases[pPutProbe->ProbeLen-1];
				pPutSeq = pLeft;
				CurMismatches = 0;
				for(PutSeqIdx = 0; PutSeqIdx < pPutProbe->ProbeLen; PutSeqIdx++,pPutSeq++,pPutProbeSeq--)
					{
					ProbeBase = *pPutProbeSeq & ~ cRptMskFlg;
					SeqBase = *pPutSeq & ~ cRptMskFlg;
					if(SeqBase > eBaseT || ProbeBase > eBaseT)
						break;
					switch(ProbeBase) {
						case eBaseA:
							if(SeqBase == eBaseT)
								continue;
							break;
						case eBaseC:
							if(SeqBase == eBaseG)
								continue;
							break;
						case eBaseG:
							if(SeqBase == eBaseC)
								continue;
							break;
						case eBaseT:
							if(SeqBase == eBaseA)
								continue;
							break;
						}					
					if(++CurMismatches > MaxMismatches)			// too many mismatches?
						break;
					}
				if(PutSeqIdx == pPutProbe->ProbeLen)
					{
					pPutProbe->TotMinHits++;
					pPutProbe->PrvEntryID = CurEntry;
					pPutProbe->PrvLociHit = Loci;
					MinHits++;
					if(Mode==ePMHitLoci)
						{
						pszAsciiBases = CSeqTrans::MapSeq2Ascii(pPutProbe->Bases,pPutProbe->ProbeLen);
						RsltsLen = sprintf((char *)szRsltsBuff,"%d,\"%s\",\"%s\",\"%c\",%lu,%d,\"%s\"\n",++gMatchCnt,pPutProbe->szDescr,pszChrom,'-',Loci,pPutProbe->ProbeLen,pszAsciiBases);
						write(hRsltsFile,szRsltsBuff,RsltsLen);
						}
					}
				}
			pPutProbe = pPutProbe->pNextMin2Hashed;
			}
		}
	pLeft++;
	Loci++;
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Hits on plus strand: %d on minus strand: %d",PlusHits,MinHits);
return(PlusHits + MinHits);
}
