// genhyperdropouts.cpp : Defines the entry point for the console application.
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

const unsigned int cProgVer = 219;		// increment with each release

const int cDfltJoinOverlap = 4;			// join cores which overlap and have starts which at most differ by cMaxJoinOverlap
const int cMaxJoinOverlap  = 1000;		// max allowed join overlap

// processing mode
typedef enum eProcMode {
	eProcModeDropouts,				// output elements in Ref which are not in Rel
	eProcModeIntersect,				// output elements in Ref which are in Rel
	eProcModeRefUnique,				// output elements including partial elements in Ref which are not in Rel
	eProcModeCombine				// output the elements in either Ref or Rel plus union of elements in both Ref and Rel
} etProcMode;

const int cMaxLengthRange = 10000000;	// maximal element length

typedef struct TAG_sLenRangeClass {
	int ID;					// uniquely identifies this range
	int Min;				// minimum length in this range
	int Max;				// maximum length in this range
	const char *pszDescr;	// descriptive text
	}tsLenRangeClass;

// length range classes
tsLenRangeClass LenRangeClasses[] = {
	{1,0,4,"0-4"},
	{2,5,9,"5-9"},
	{3,10,14,"10-14"},
	{4,15,19,"15-19"},
	{5,20,29,"20-29"},
    {6,30,49,"30-49"},
	{7,50,74,"50-74"},
	{8,75,99,"75-99"},
	{9,100,124,"100-124"},
	{10,125,149,"125-149"},
	{11,150,174,"150-174"},
	{12,175,199,"175-199"},
	{13,200,249,"200-249"},
	{14,250,299,"250-299"},
	{15,300,349,"300-349"},
	{16,350,399,"350-399"},
	{17,400,449,"400-449"},
	{18,450,499,"450-499"},
	{19,500,599,"500-599"},
	{20,600,699,"600-699"},
	{21,700,799,"700-799"},
	{22,800,899,"800-899"},
	{23,900,999,"900-999"},
	{24,1000,1249,"1000-1249"},
	{25,1250,1499,"1250-1499"},
	{26,1500,1749,"1500-1749"},
	{27,1750,1999,"1750-1999"},
	{28,2000,INT_MAX,"2000+"}
};
const int cLenRanges = sizeof(LenRangeClasses)/sizeof(tsLenRangeClass);		  // number of length range classes


typedef struct TAG_sProcParams 
	{
	int hLociFile;
	int hStatsFile;
	etProcMode ProcMode;	// processing mode
	int OverlapBases;		// minimum required absolute overlap between ref and rel element
	int OverlapPercentage;  // percentage overlap of ref element (1..100) 0 to ignore this param
	int MinLength;		// only process ref elements of at least this length
	int MaxLength;		// only process ref elements which are no longer than this length
	int JoinOverlap;	// join elements which overlap and have starts which only differ by this many bases
	char szRefFile[_MAX_PATH];	// CSV file containing ref elements
	char szRelFile[_MAX_PATH];	// CSV file containing rel elements
	char szOutStatsFile[_MAX_PATH]; // write counts/stats into this CSV file 
	char szOutLociFile[_MAX_PATH];  // write element loci into this CSV file

	int NumRefEls;
	CHyperEls *pRefHypers;
	int NumRelEls;
	CHyperEls *pRelHypers;

	char szRefSpecies[cMaxDatasetSpeciesChrom];	// replace element's ref species with this in OutLociFile
	char szRelSpecies[cMaxDatasetSpeciesChrom];	// replace element's rel species with this in OutLociFile
	char szElType[cMaxDatasetSpeciesChrom];		// replace element's type with this in OutLociFile
} tsProcParams; 



bool CleanupResources(tsProcParams *pProcParams);

int 
Process(etProcMode ProcMode,	// processing mode
		int OverlapBases,		// minimum required overlap between ref and rel element
		int OverlapPercentage, // minimum required percentage overlap of ref element (1..100) 0 to ignore this param
		int MinLength,		// only process ref elements of at least this length
		int MaxLength,		// only process ref elements which are no longer than this length
		int JoinOverlap,	// join elements which overlap and have starts which only differ by this many bases
		char *pszRefFile,	// CSV file containing ref elements
		char *pszRelFile,	// CSV file containing rel elements
		char *pszOutStatsFile, // write counts/stats into this CSV file 
		char *pszOutLociFile,  // write element loci into this CSV file
		char *pszRefSpecies,	// replacement ref species in output loci file
		char *pszRelSpecies,	// replacement rel species in output loci file
		char *pszElType);		// replacement element type in output loci file

int GenCombined(tsProcParams *pProcParams);
int FindDropouts(tsProcParams *pProcParams);
int FindIntersects(tsProcParams *pProcParams);
bool OutputDropoutResults(int *pCntStepCnts,tsProcParams *pProcParams);
bool OutputRefUniques(int CoreID,char *pszElType,char *pszSpecies,char *pszSpeciesList,char *pszChrom,int StartLoci, int EndLoci,tsProcParams *pProcParams);
int GenRefUniques(tsProcParams *pProcParams);
char *ProcMode2Txt(etProcMode ProcMode);

int				// returns number of elements parsed
ParseFileElements(int SrcID,		// source identifier to associate with elements parsed from this file
				  char *pszFile);	// file (in CSV format) containing element loci


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
	return _T("genhyperdropouts");
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
int iProcMode;
int iOverlapBases;
int iOverlapPercentage;
int iMinLength;
int iMaxLength;
int iJoinOverlap;
char szRefFile[_MAX_PATH];	// process ref hypers from this file
char szRelFile[_MAX_PATH];	// process rel hypers from this file
char szOutStatsFile[_MAX_PATH];	// write stats to this file
char szOutLociFile[_MAX_PATH];	// write loci to this file
char szRefSpecies[cMaxDatasetSpeciesChrom];	// use this species as the ref species in generated szOutLociFile
char szRelSpecies[cMaxDatasetSpeciesChrom];	// use this species/list as the rel species in generated szOutLociFile
char szElType[cMaxDatasetSpeciesChrom];		// use this as the element type in generated szOutLociFile

// command line args
struct arg_lit  *help    = arg_lit0("hH","help",                "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_int *ScreenLogLevel=arg_int0("S", "ScreenLogLevel",	"<int>","Level of diagnostics written to logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"diagnostics log file");

struct arg_file *RefFile = arg_file1("i","reffile","<file>",	"reference hyper element CSV file");
struct arg_file *RelFile = arg_file1("I","relfile","<file>",	"relative hyper element CSV file");
struct arg_file *OutStatsFile = arg_file0("o",NULL,"<file>",	"output to counts statistics file as CSV");
struct arg_file *OutLociFile = arg_file0("O",NULL,"<file>",		"output loci to file as CSV");

struct arg_str  *RefSpecies = arg_str0("r","refspecies","<string>","output loci file ref species");
struct arg_str  *RelSpecies = arg_str0("R","relspecies","<string>","output loci file rel species");
struct arg_str  *ElType = arg_str0("t","eltype","<string>","output loci file element type");

struct arg_int  *ProcMode = arg_int0("p","mode","<int>",		 "processing mode: 0:Dropouts, 1:Intersect, 2:RefUnique, 3:Combined (default 0)");
struct arg_int  *OverlapBases = arg_int0("l","OverlapBases","<int>", "minimum number of overlap bases required (default 10)");
struct arg_int  *OverlapPercentage = arg_int0("L","minpercent","<int>","minimum overlap as a percentage 0..100 (default 50) of ref element length");
struct arg_int  *MinLength = arg_int0("m","minlength","<int>",   "minimum ref element length (default 0)");
struct arg_int  *MaxLength = arg_int0("M","maxlength","<int>",   "maximum ref element length (default 1000000)");
struct arg_int  *JoinOverlap = arg_int0("j","joinoverlap","<int>","joins cores which only differ by this many bases in start loci (default 4)");



struct arg_end *end = arg_end(20);

void *argtable[] = {help,version,FileLogLevel,ScreenLogLevel,LogFile,
					ProcMode,
					RefFile,RelFile,OutStatsFile,OutLociFile,
					OverlapBases,OverlapPercentage,MinLength,MaxLength,JoinOverlap,
					RefSpecies,RelSpecies,ElType,
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


	iProcMode = ProcMode->count ? ProcMode->ival[0] : eProcModeDropouts;
	if(iProcMode < eProcModeDropouts || iProcMode > eProcModeCombine)
		{
		printf("Error: Processing mode '-p%d' is not in range 0..3",iProcMode);
		exit(1);
		}
	switch(iProcMode) {
		case eProcModeDropouts:				// output elements in Ref which are not in Rel
		case eProcModeIntersect:			// output elements in Ref which are in Rel
			if(!OutStatsFile->count && !OutLociFile->count)
				{
				printf("\nError: No output file(s) specified with '-o<statsfile>; or '-O<locifile>'.");
				exit(1);
				}
			if(OutStatsFile->count)
				{
				strncpy(szOutStatsFile,OutStatsFile->filename[0],_MAX_PATH);
				szOutStatsFile[_MAX_PATH-1] = '\0';
				}
			else
				szOutStatsFile[0] = '\0';
			break;

		case eProcModeCombine:				// output the elements in either Ref or Rel plus union of elements in both Ref and Rel
		case eProcModeRefUnique:			// output elements including partial elements in Ref which are not in Rel
			if(!OutLociFile->count)
				{
				printf("\nError: No output file specified with '-O<locifile>'.");
				exit(1);
				}
			szOutStatsFile[0] = '\0';
			break;
		
		}
	if(OutLociFile->count)
		{
		strncpy(szOutLociFile,OutLociFile->filename[0],_MAX_PATH);
		szOutLociFile[_MAX_PATH-1] = '\0';

		if(RefSpecies->count)
			{
			strncpy(szRefSpecies,RefSpecies->sval[0],sizeof(szRefSpecies));
			szRefSpecies[sizeof(szRefSpecies)-1] = '\0';
			}
		else
			szRefSpecies[0] = '\0';

		if(RelSpecies->count)
			{
			strncpy(szRelSpecies,RelSpecies->sval[0],sizeof(szRelSpecies));
			szRelSpecies[sizeof(szRelSpecies)-1] = '\0';
			}
		else
			szRelSpecies[0] = '\0';

		if(ElType->count)
			{
			strncpy(szElType,ElType->sval[0],sizeof(szElType));
			szElType[sizeof(szElType)-1] = '\0';
			}
		else
			szElType[0] = '\0';
		}
	else
		{
		szOutLociFile[0] = '\0';
		szRefSpecies[0] = '\0';
		szRelSpecies[0] = '\0';
		szElType[0] = '\0';
		}

	iMinLength = MinLength->count ? MinLength->ival[0] : 0;
	if(iMinLength < 0 || iMinLength > cMaxLengthRange)
		{
		printf("Error: Minimum element length '-m%d' is not in range 0..%d",iMinLength,cMaxLengthRange);
		exit(1);
		}

	iMaxLength = MaxLength->count ? MaxLength->ival[0] : cMaxLengthRange;
	if(iMaxLength < iMinLength || iMaxLength > cMaxLengthRange)
		{
		printf("Error: Maximum element length '-M%d' is not in range %d..%d",iMaxLength,iMinLength,cMaxLengthRange);
		exit(1);
		}

	if(iProcMode != eProcModeCombine)
		{
		iOverlapBases = OverlapBases->count ? OverlapBases->ival[0] : 5;
		if(iOverlapBases < 0 || iOverlapBases > iMaxLength)
			{
			printf("Error: Minimum overlap length '-l%d' is not in range 1..%d",iOverlapBases,iMaxLength);
			exit(1);
			}
		iOverlapPercentage = OverlapPercentage->count ? OverlapPercentage->ival[0] : 0;
		if(iOverlapPercentage < 0 || iOverlapPercentage > 100)
			{
			printf("Error: Minimum percentage overlap '-L%d' is not in range 0..100",iOverlapPercentage);
			exit(1);
			}

		iJoinOverlap = JoinOverlap->count ? JoinOverlap->ival[0] : cDfltJoinOverlap;
		if(iJoinOverlap < 0 || iJoinOverlap > cMaxJoinOverlap)
			{
			printf("Error: join overlap length '-j%d' is not in range 0..%d",iJoinOverlap,cMaxJoinOverlap);
			exit(1);
			}
		}
	else
		{
		iOverlapBases = 0;
		iOverlapPercentage = 0;
		iJoinOverlap = 0;
		}

	strncpy(szRefFile,RefFile->filename[0],_MAX_PATH);
	szRefFile[_MAX_PATH-1] = '\0';
	strncpy(szRelFile,RelFile->filename[0],_MAX_PATH);
	szRelFile[_MAX_PATH-1] = '\0';



		// now that command parameters have been parsed then initialise diagnostics log system
	if(!gDiagnostics.Open(szLogFile,(etDiagLevel)iScreenLogLevel,(etDiagLevel)iFileLogLevel,true))
		{
		printf("\nError: Unable to start diagnostics subsystem.");
		if(szLogFile[0] != '\0')
			printf(" Most likely cause is that logfile '%s' can't be opened/created",szLogFile);
		exit(1);
		}
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Version: %d.%2.2d Processing parameters:",cProgVer/100,cProgVer%100);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Processing Mode: %d (%s)",iProcMode,ProcMode2Txt((etProcMode)iProcMode));
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Reference CSV file: '%s'",szRefFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Relative CSV file: '%s'",szRelFile);
	if(szOutStatsFile[0] != '\0')
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Write stats into CSV file: '%s'",szOutStatsFile);
	if(szOutLociFile[0] != '\0')
		{
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Output processed loci into CSV file: '%s'",szOutLociFile);
		if(szRefSpecies[0] != '\0')
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Replacement output loci file ref species: '%s'",szRefSpecies);
		if(szRelSpecies[0] != '\0')
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Replacement output loci file rel species: '%s'",szRelSpecies);
		if(szElType[0] != '\0')
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Replacement output loci file element type: '%s'",szElType);
		}

	if(iProcMode != eProcModeCombine)
		{
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Minimum bases overlap: %d",iOverlapBases);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Minimum percentage overlap: %d",iOverlapPercentage);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Join overlap: %d",iJoinOverlap);
		}

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Minimum element length: %d",iMinLength);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Maximum element length: %d",iMaxLength);


	// processing here...
	gStopWatch.Start();
#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif	
	Rslt = Process((etProcMode)iProcMode,iOverlapBases,iOverlapPercentage,iMinLength,iMaxLength,iJoinOverlap,szRefFile,szRelFile,szOutStatsFile,szOutLociFile,
		szRefSpecies,szRelSpecies,szElType);
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
ProcMode2Txt(etProcMode ProcMode)
{
switch(ProcMode) {
	case eProcModeDropouts:			// output elements in Ref which are not in Rel
		return((char *)"Dropouts");
	case eProcModeIntersect:		// output elements in Ref which are in Rel
		return((char *)"Intersects");
	case eProcModeRefUnique:		// output elements including partial elements in Ref which are not in Rel
		return((char *)"Uniques");
	case eProcModeCombine:			// output the elements in either Ref or Rel plus union of elements in both Ref and Rel
		return((char *)"Combined");

	default:
		break;
	}
return((char *)"Unrecognised");
}

// returns highest priority functional region from feature bits
teFuncRegion 
MapFeatureBits2Idx(int FeatureBits)
{
if(FeatureBits & cFeatBitCDS)	// CDS has the highest priority
	return(eFRCDS);
if(FeatureBits & cFeatBit5UTR)
	return(eFR5UTR);
if(FeatureBits & cFeatBit3UTR)
	return(eFR3UTR);
if(FeatureBits & cFeatBitIntrons)
	return(eFRIntronic);
if(FeatureBits & cFeatBitUpstream)
	return(eFRUpstream);
if(FeatureBits & cFeatBitDnstream)
	return(eFRDnstream);
return(eFRIntergenic);			// and intergenic the lowest
}

// GetLengthRangeClass
// Returns ptr to length range class for specified Length, or NULL if can't classify into a range
tsLenRangeClass *GetLengthRangeClass(int Length)
{
int Idx;
tsLenRangeClass *pRange = LenRangeClasses;
for(Idx = 0; Idx < cLenRanges; Idx++,pRange++)
	if(Length >= pRange->Min && Length <= pRange->Max)
		return(pRange);
return(NULL);
}

// GetRangeClass
// Returns ptr to length range class for specified range identifier, or NULL if can't classify into a range
tsLenRangeClass *GetRangeClass(int RangeID)
{
if(RangeID < 1 || RangeID > cLenRanges)
	return(NULL);
return(&LenRangeClasses[RangeID-1]);
}

int 
Process(etProcMode ProcMode,		// processing mode
		int OverlapBases,		// minimum required absolute overlap between ref and rel element
		int OverlapPercentage, // minimum required percentage overlap of ref element (1..100) 0 to ignore this param
		int MinLength,			  // only process ref elements of at least this length
		int MaxLength,			  // only process ref elements which are no longer than this length
		int JoinOverlap,		  // join elements which overlap and have starts which only differ by this many bases
		char *pszRefFile,		  // CSV file containing ref elements
		char *pszRelFile,		  // CSV file containing rel elements
		char *pszOutStatsFile,	  // write counts/stats into this CSV file 
		char *pszOutLociFile,	  // write element loci into this CSV file
		char *pszRefSpecies,		// replacement ref species in output loci file
		char *pszRelSpecies,	// replacement rel species in output loci file
		char *pszElType)		// replacement element type in output loci file
{
int Rslt;
int TmpNumRefEls;
int TmpNumRelEls;

tsProcParams ProcParams;

memset(&ProcParams,0,sizeof(ProcParams));

ProcParams.MaxLength = MaxLength;
ProcParams.MinLength = MinLength;
ProcParams.OverlapBases = OverlapBases;
ProcParams.JoinOverlap = JoinOverlap;
ProcParams.OverlapPercentage = OverlapPercentage;
ProcParams.ProcMode = ProcMode;
if(pszOutStatsFile[0] != '\0')
	strcpy(ProcParams.szOutStatsFile,pszOutStatsFile);
if(pszOutLociFile[0] != '\0')
	strcpy(ProcParams.szOutLociFile,pszOutLociFile);
if(pszRefFile[0] != '\0')
	strcpy(ProcParams.szRefFile,pszRefFile);
if(pszRelFile[0] != '\0')
	strcpy(ProcParams.szRelFile,pszRelFile);
if(pszRefSpecies[0] != '\0')
	strcpy(ProcParams.szRefSpecies,pszRefSpecies);
if(pszRelSpecies[0] != '\0')
	strcpy(ProcParams.szRelSpecies,pszRelSpecies);
if(pszElType[0] != '\0')
	strcpy(ProcParams.szElType,pszElType);

if(pszOutStatsFile[0] != '\0')
	{

#ifdef _WIN32
	if((ProcParams.hStatsFile = open(pszOutStatsFile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
	if((ProcParams.hStatsFile = open(pszOutStatsFile,O_RDWR | O_CREAT |O_TRUNC, S_IREAD | S_IWRITE ))==-1)
#endif	
	{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create stats file: %s - %s",pszOutStatsFile,strerror(errno));
		CleanupResources(&ProcParams);
		return(eBSFerrCreateFile);
		}
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Stats CSV file created/truncated: '%s'",pszOutStatsFile);
	}


if(pszOutLociFile[0] != '\0')
	{
#ifdef _WIN32
	if((ProcParams.hLociFile = open(pszOutLociFile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
	if((ProcParams.hLociFile = open(pszOutLociFile, O_RDWR | O_CREAT | O_TRUNC, S_IREAD | S_IWRITE ))==-1)
#endif
	{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create loci file: %s - %s",pszOutLociFile,strerror(errno));
		CleanupResources(&ProcParams);
		return(eBSFerrCreateFile);
		}
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Loci CSV file created/truncated: '%s'",pszOutLociFile);
	}

// parse reference element file
ProcParams.pRefHypers = new CHyperEls;

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Parsing ref CSV file: '%s'",pszRefFile);
if((Rslt = ProcParams.pRefHypers->ParseCSVFileElements(pszRefFile,MinLength,MaxLength)) < 0)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Parse errors in CSV file: '%s'",pszRefFile);
	CleanupResources(&ProcParams);
	return(Rslt);
	}

TmpNumRefEls = ProcParams.pRefHypers->NumEls();
if(TmpNumRefEls == 0)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"No reference elements with length range %d..%d in CSV file: '%s'",MinLength,MaxLength,pszRefFile);
	CleanupResources(&ProcParams);
	return(Rslt);
	}

ProcParams.NumRefEls = ProcParams.pRefHypers->DedupeSort(JoinOverlap);				// dedupe and sort elements
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Ref Elements %d - Parsed: %d FiltLen: %d FiltDup: %d",
											ProcParams.NumRefEls,
											ProcParams.pRefHypers->NumElsParsed(),
											ProcParams.pRefHypers->NumElsFiltLen(),
											ProcParams.pRefHypers->NumElsDeduped());

if(ProcMode == eProcModeCombine)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Parsing rel CSV file: '%s'",pszRelFile);
	if((Rslt=ProcParams.pRefHypers->ParseCSVFileElements(pszRelFile,MinLength,MaxLength)) < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Parse errors in CSV file: '%s'",pszRelFile);
		CleanupResources(&ProcParams);	
		return(Rslt);
		}
	TmpNumRelEls = ProcParams.pRefHypers->NumEls() - TmpNumRefEls;
	if((TmpNumRelEls +  TmpNumRefEls) == 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"No elements in length range %d..%d in CSV files: '%s' or '%s'",
				MinLength,MaxLength,pszRefFile,pszRelFile);
		CleanupResources(&ProcParams);
		return(Rslt);
		}
	if(TmpNumRelEls == 0)
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"No relative elements in length range %d..%d in CSV file: '%s'",
							MinLength,MaxLength,pszRelFile);

	ProcParams.NumRefEls = ProcParams.pRefHypers->DedupeSort(JoinOverlap);		// dedupe and sort elements
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Ref Elements %d  - Rel Parsed: %d FiltLen: %d FiltDup: %d",
						    ProcParams.NumRefEls,
							ProcParams.pRefHypers->NumElsParsed(),
							ProcParams.pRefHypers->NumElsFiltLen(),
							ProcParams.pRefHypers->NumElsDeduped());
	}
else
	{
	// parse relative element file
	// if a minimum overlay was specified then relative elements must
	// be of at least that length
	ProcParams.pRelHypers = new CHyperEls;
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Parsing rel CSV file (%d..%d): '%s'",MinLength,MaxLength,pszRelFile);
	if((Rslt=ProcParams.pRelHypers->ParseCSVFileElements(pszRelFile,MinLength,MaxLength)) < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Parse errors in CSV file: '%s'",pszRelFile);
		CleanupResources(&ProcParams);	
		return(Rslt);
		}
	TmpNumRelEls = ProcParams.pRelHypers->NumEls();
	if(TmpNumRelEls == 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"No relative elements with length range %d..%d in CSV file: '%s'",MinLength,cMaxLengthRange,pszRelFile);
		CleanupResources(&ProcParams);
		return(Rslt);
		}

	ProcParams.NumRelEls = ProcParams.pRelHypers->DedupeSort(JoinOverlap);		// dedupe and sort elements
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Rel Elements %d  - Parsed: %d FiltLen: %d FiltDup: %d",
						    ProcParams.NumRelEls,
							ProcParams.pRelHypers->NumElsParsed(),
							ProcParams.pRelHypers->NumElsFiltLen(),
							ProcParams.pRelHypers->NumElsDeduped());
	}

switch(ProcMode) {
	case eProcModeDropouts:				// output elements in Ref which are not at least partially overlapped by an element in Rel
		Rslt = FindDropouts(&ProcParams);
		break;
	case eProcModeIntersect:			// output elements in Ref which are also at least partially overlapped by an element in Rel
		Rslt = FindIntersects(&ProcParams);
		break;

	case eProcModeRefUnique:				// output elements including partial elements in Ref which are not in Rel
		Rslt = GenRefUniques(&ProcParams);
		break;

	case eProcModeCombine:				// output the elements in either Ref or Rel plus union of elements in both Ref and Rel
		Rslt = GenCombined(&ProcParams);
		break;
	}


CleanupResources(&ProcParams);
return(Rslt > 0 ? 0 : Rslt);
}


bool 
CleanupResources(tsProcParams *pProcParams)
{
if(pProcParams->hLociFile > 0)
	{
	close(pProcParams->hLociFile);
	pProcParams->hLociFile = 0;
	}
if(pProcParams->hStatsFile > 0)
	{
	close(pProcParams->hStatsFile);
	pProcParams->hStatsFile = 0;
	}
if(pProcParams->pRefHypers != NULL)
	{
	delete pProcParams->pRefHypers;
	pProcParams->pRefHypers = NULL;
	}
pProcParams->NumRefEls = 0;
if(pProcParams->pRelHypers != NULL)
	{
	delete pProcParams->pRelHypers;
	pProcParams->pRelHypers = NULL;
	}
pProcParams->NumRelEls = 0;

return(true);
}


bool
UpdateSpeciesElType(tsHyperElement *pCurEl,
					int *pCurRefSpeciesID,
					int *pCurRelSpeciesID,
					int *pCurElTypeID,
					char **ppszRefSpecies,
					char **ppszRelSpecies,
					char **ppszElType,
					tsProcParams *pProcParams)
{
bool bUpdated = false;
if(*ppszRefSpecies == NULL || *pCurRefSpeciesID != pCurEl->RefSpeciesID)
	{
	if(pProcParams->szRefSpecies[0] != '\0')
		*ppszRefSpecies = pProcParams->szRefSpecies;
	else
		*ppszRefSpecies = pProcParams->pRefHypers->GetRefSpecies(pCurEl->RefSpeciesID);
	*pCurRefSpeciesID = pCurEl->RefSpeciesID;
	bUpdated = true;
	}

if(*ppszRelSpecies == NULL || *pCurRelSpeciesID != pCurEl->RelSpeciesID)
	{
	if(pProcParams->szRelSpecies[0] != '\0')
		*ppszRelSpecies = pProcParams->szRelSpecies;
	else
		*ppszRelSpecies = pProcParams->pRefHypers->GetRelSpecies(pCurEl->RelSpeciesID);
	*pCurRelSpeciesID = pCurEl->RelSpeciesID;
	bUpdated = true;
	}

if(*ppszElType == NULL || *pCurElTypeID != pCurEl->ElTypeID)
	{
	if(pProcParams->szElType[0] != '\0')
		*ppszElType = pProcParams->szElType;
	else
		*ppszElType = pProcParams->pRefHypers->GetType(pCurEl->ElTypeID);
	*pCurElTypeID = pCurEl->ElTypeID;
	bUpdated = true;
	}
return(bUpdated);
}

int
GenCombined(tsProcParams *pProcParams)
{
int CurElID;
int CurChromID;
int CurStartLoci;
int CurEndLoci;
int CurLength;
int NxtEndLoci;
int NumCombinedEls;
int RefChromID;
char *pszRefSpecies;
char *pszRelSpecies;
char *pszElType;
char szLineBuff[cMaxReadLen+1];
int LineLen;
tsHyperElement *pCurEl;
int CurRefSpeciesID;
int CurRelSpeciesID;
int CurElTypeID;
char *pszRefChrom = NULL;

// now iterate over all elements combining those that overlap
NumCombinedEls = 0;
RefChromID = -1;
CurRefSpeciesID = -1;
CurRelSpeciesID = -1;
CurElTypeID = -1;
CurLength = 0;
for(CurElID = 1; CurElID <= pProcParams->NumRefEls; CurElID++)
	{
	pCurEl=pProcParams->pRefHypers->GetElement(CurElID);
	NxtEndLoci = pCurEl->StartLoci + pCurEl->Len - 1;
	if(CurElID == 1)
		{
		CurChromID = pCurEl->ChromID;
		CurStartLoci = pCurEl->StartLoci;
		CurEndLoci = NxtEndLoci;

		UpdateSpeciesElType(pCurEl,
					&CurRefSpeciesID,
					&CurRelSpeciesID,
					&CurElTypeID,
					&pszRefSpecies,
					&pszRelSpecies,
					&pszElType,
					pProcParams);
		}
	else
		{
		if(pCurEl->ChromID != CurChromID ||
			pCurEl->StartLoci > (CurEndLoci + 1))
			{
			if(pszRefChrom == NULL || RefChromID != CurChromID)
				{
				pszRefChrom = pProcParams->pRefHypers->GetChrom(CurChromID);
				RefChromID = CurChromID;
				}
			
			CurLength = 1 + CurEndLoci - CurStartLoci;

			if(CurLength >= pProcParams->MinLength && 
				CurLength <= pProcParams->MaxLength && 
				pProcParams->hLociFile > 0)
				{
				LineLen=sprintf(szLineBuff,"%d,\"%s\",\"%s\",\"%s\",%d,%d,%d,\"%s\",0\n",
								++NumCombinedEls,
								pszElType,
								pszRefSpecies,
								pszRefChrom,
								CurStartLoci,CurEndLoci,CurLength,
								pszRelSpecies);
				CUtility::SafeWrite(pProcParams->hLociFile,szLineBuff,LineLen);
				}
			CurChromID = pCurEl->ChromID;
			CurStartLoci = pCurEl->StartLoci;
			CurEndLoci = NxtEndLoci;

			UpdateSpeciesElType(pCurEl,
					&CurRefSpeciesID,
					&CurRelSpeciesID,
					&CurElTypeID,
					&pszRefSpecies,
					&pszRelSpecies,
					&pszElType,
					pProcParams);
			}
		else
			{
			if(NxtEndLoci > CurEndLoci)
				CurEndLoci = NxtEndLoci;
			}
		}
	}

CurLength = 1 + CurEndLoci - CurStartLoci;
if(CurLength >= pProcParams->MinLength && 
	CurLength <= pProcParams->MaxLength && 
	pProcParams->hLociFile > 0)
	{
	LineLen=sprintf(szLineBuff,"%d,\"%s\",\"%s\",\"%s\",%d,%d,%d,\"%s\",0\n",
								++NumCombinedEls,
								pszElType,
								pszRefSpecies,
								pszRefChrom,
								CurStartLoci,CurEndLoci,CurLength,
								pszRelSpecies);
	CUtility::SafeWrite(pProcParams->hLociFile,szLineBuff,LineLen);
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Total original elements: %d, combined elements: %d", pProcParams->NumRefEls,NumCombinedEls);
return(NumCombinedEls);
}

int
FindDropouts(tsProcParams *pProcParams)
{
int RefElID;
int RelElID;
int RefOverlapBases;
tsHyperElement *pRefEl;
tsHyperElement *pRelEl;
char *pszRefChrom;
char *pszRefSpecies;
char *pszRelSpecies;
char *pszElType;
int RelChromID;
int CurRefSpeciesID;
int CurRelSpeciesID;
int CurElTypeID;

int RefEndLoci;
int RefStartLociMin;
int RefStartLociMax;
int RefEndLociMin;
int RefEndLociMax;
int RelEndLoci;
int nThInstance;

int NumElsLocated;
char szLineBuff[cMaxReadLen+1];
int LineLen;
int StepCnts[7*cLenRanges];
// now iterate over all ref elements locating rel elements
NumElsLocated = 0;
if(pProcParams->OverlapPercentage < 1)
	RefOverlapBases = pProcParams->OverlapBases;
memset(StepCnts,0,sizeof(StepCnts));

CurRefSpeciesID = -1;
CurRelSpeciesID = -1;
CurElTypeID = -1;

for(RefElID = 1; RefElID <= pProcParams->NumRefEls; RefElID++)
	{
	pRefEl=pProcParams->pRefHypers->GetElement(RefElID);
	pszRefChrom = pProcParams->pRefHypers->GetChrom(pRefEl->ChromID);
	RelChromID = pProcParams->pRelHypers->GetChromID(pszRefChrom);
	RefEndLoci = pRefEl->StartLoci+pRefEl->Len -1;

	UpdateSpeciesElType(pRefEl,
					&CurRefSpeciesID,
					&CurRelSpeciesID,
					&CurElTypeID,
					&pszRefSpecies,
					&pszRelSpecies,
					&pszElType,
					pProcParams);

	pRelEl = NULL;
	if(RelChromID >= 1)	
		{
		if(pProcParams->OverlapPercentage >= 1)
			{
			RefOverlapBases = (pRefEl->Len * pProcParams->OverlapPercentage)/100;
			if(RefOverlapBases < pProcParams->OverlapBases)
				RefOverlapBases = pProcParams->OverlapBases;
			}

		if((RefStartLociMin = pRefEl->StartLoci - RefOverlapBases) < 0)
			RefStartLociMin = 0;
		RefStartLociMax = pRefEl->StartLoci + RefOverlapBases;

		if((RefEndLociMin = RefEndLoci - RefOverlapBases) < RefStartLociMin)
			RefEndLociMin = RefStartLociMin;
		RefEndLociMax = RefEndLoci + RefOverlapBases;
		nThInstance = 1;
		
			// see if any rel elements overlap ref element
		while((RelElID = pProcParams->pRelHypers->Locate(RelChromID,pRefEl->StartLoci,pRefEl->Len,RefOverlapBases > 0 ? RefOverlapBases : 1,nThInstance++))>0)
			{
			pRelEl = pProcParams->pRelHypers->GetElement(RelElID);
			RelEndLoci = pRelEl->StartLoci+pRelEl->Len - 1;
			if(pRelEl->StartLoci <= RefStartLociMax &&
     			pRelEl->StartLoci >= RefStartLociMin &&
				RelEndLoci <= RefEndLociMax &&
				RelEndLoci >= RefEndLociMin)
					break;
			}
		if(RelElID > 0)
			continue;
		}

	// else - no rel element overlaps ref element
	NumElsLocated++; // at least one element of interest

	// classify range
	tsLenRangeClass *pRange;
	if((pRange = GetLengthRangeClass(pRefEl->Len))==NULL)
		{
		gDiagnostics.DiagOut(eDLWarn,gszProcName,"Unexpected length range classification error - length = %d",pRefEl->Len);
		break;
		}

	int FeatIdx = MapFeatureBits2Idx(pRefEl->Features);
	StepCnts[((pRange->ID-1) * 7) + FeatIdx] += 1;


	if(pProcParams->hLociFile > 0)
		{
		LineLen=sprintf(szLineBuff,"%d,\"%s\",\"%s\",\"%s\",%d,%d,%d,\"%s\",%d,%d\n",
			pRefEl->SrcID,
			pszElType,
			pszRefSpecies,
			pszRefChrom,
			pRefEl->StartLoci,
			RefEndLoci,
			pRefEl->Len,
			pszRelSpecies,
			pRefEl->Features,
			pRelEl == NULL ? 0 : pRelEl->SrcID);
		CUtility::SafeWrite(pProcParams->hLociFile,szLineBuff,LineLen);
		}
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Total ref elements: %d, matching elements: %d", pProcParams->NumRefEls,NumElsLocated);
OutputDropoutResults(StepCnts,pProcParams);
return(NumElsLocated);
}

// FindIntersects
int
FindIntersects(tsProcParams *pProcParams)
{
int FeatIdx;
int RelIdx;
int RefElID;
int RelElID;
int RefOverlapBases;
int OverlapLoci;
int OverlapLen;
tsHyperElement *pRefEl;
char *pszRefChrom;
char *pszRefSpecies;
char *pszRelSpecies;
char *pszElType;
int RelChromID;
int CurRefSpeciesID;
int CurRelSpeciesID;
int CurElTypeID;
int NumElsLocated;
char szLineBuff[cMaxReadLen+1];
int LineLen;
int StepCnts[7*cLenRanges];
// now iterate over all ref elements locating rel elements
NumElsLocated = 0;
if(pProcParams->OverlapPercentage < 1)
	RefOverlapBases = pProcParams->OverlapBases;
memset(StepCnts,0,sizeof(StepCnts));
CurRefSpeciesID = -1;
CurRelSpeciesID = -1;
CurElTypeID = -1;
for(RefElID = 1; RefElID <= pProcParams->NumRefEls; RefElID++)
	{
	pRefEl=pProcParams->pRefHypers->GetElement(RefElID);
	pszRefChrom = pProcParams->pRefHypers->GetChrom(pRefEl->ChromID);
	RelChromID = pProcParams->pRelHypers->GetChromID(pszRefChrom);

	if(RelChromID < 1)
		continue;
	
	UpdateSpeciesElType(pRefEl,
					&CurRefSpeciesID,
					&CurRelSpeciesID,
					&CurElTypeID,
					&pszRefSpecies,
					&pszRelSpecies,
					&pszElType,
					pProcParams);


	if(pProcParams->OverlapPercentage >= 1)
		{
		RefOverlapBases = (pRefEl->Len * pProcParams->OverlapPercentage)/100;
		if(RefOverlapBases < pProcParams->OverlapBases)
			RefOverlapBases = pProcParams->OverlapBases;
		}

	// see if any rel elements overlap ref element by at least RefOverlapBases bases
	RelIdx=1;
	while((RelElID = pProcParams->pRelHypers->Locate(RelChromID,pRefEl->StartLoci,pRefEl->Len,RefOverlapBases,RelIdx++,&OverlapLoci,&OverlapLen)) > 0)
		{
		NumElsLocated++; // at least one element of interest

		// classify range
		tsLenRangeClass *pRange;
		if((pRange = GetLengthRangeClass(OverlapLen))==NULL)
			{
			gDiagnostics.DiagOut(eDLWarn,gszProcName,"Unexpected length range classification error - length = %d",OverlapLen);
			break;
			}

		FeatIdx = MapFeatureBits2Idx(pRefEl->Features);
		StepCnts[((pRange->ID-1) * 7) + FeatIdx] += 1;

		if(pProcParams->hLociFile > 0)
			{
			LineLen=sprintf(szLineBuff,"%d,\"%s\",\"%s\",\"%s\",%d,%d,%d,\"%s\",%d\n",
				pRefEl->SrcID,
				pszElType,
				pszRefSpecies,
				pszRefChrom,
				OverlapLoci,OverlapLoci + OverlapLen - 1,OverlapLen,
				pszRelSpecies,
				pRefEl->Features);
			CUtility::SafeWrite(pProcParams->hLociFile,szLineBuff,LineLen);
			}
		}
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Total ref elements: %d, matching elements: %d", pProcParams->NumRefEls,NumElsLocated);
OutputDropoutResults(StepCnts,pProcParams);
return(NumElsLocated);
}


// GenRefUniques
// Finds all reference subelements which are not overlapped by a relative element
int
GenRefUniques(tsProcParams *pProcParams)
{
int RefElID;
int RelElID;
int RefStartLoci;
int RefEndLoci;
tsHyperElement *pRefEl;
tsHyperElement *pRelEl;
char *pszRefChrom;
int RelChromID;
int NumElsOutput;
int nthRelEl;
int RefLen;

int RefSpeciesID;
int RelSpeciesID;
char *pszRefSpecies;
char *pszRelSpecies;
char *pszElType;
int CurRefSpeciesID;
int CurRelSpeciesID;
int CurElTypeID;

// iterate over all ref elements
RefSpeciesID = -1;
RelSpeciesID = -1;
CurRefSpeciesID = -1;
CurRelSpeciesID = -1;
CurElTypeID = -1;
NumElsOutput = 0;
for(RefElID = 1; RefElID <= pProcParams->NumRefEls; RefElID++)
	{
	pRefEl=pProcParams->pRefHypers->GetElement(RefElID);
	pszRefChrom = pProcParams->pRefHypers->GetChrom(pRefEl->ChromID);
	UpdateSpeciesElType(pRefEl,
					&CurRefSpeciesID,
					&CurRelSpeciesID,
					&CurElTypeID,
					&pszRefSpecies,
					&pszRelSpecies,
					&pszElType,
					pProcParams);

	RelChromID = pProcParams->pRelHypers->GetChromID(pszRefChrom);
	if(pRefEl->Len >= pProcParams->MinLength && RelChromID >= 1)	// at least one relative element on same chromosome
		{	
			// see if any rel elements partially overlap ref element
		nthRelEl = 1;
		RefStartLoci = pRefEl->StartLoci;
		RefEndLoci = pRefEl->Len + RefStartLoci - 1;
		RefLen = pRefEl->Len;
		while(RefLen > 0 && RefLen >= pProcParams->MinLength && (RelElID = pProcParams->pRelHypers->Locate(RelChromID,pRefEl->StartLoci,pRefEl->Len,1,nthRelEl++))>0)
			{
			pRelEl=pProcParams->pRelHypers->GetElement(RelElID);
			
			if(pRelEl->StartLoci > RefStartLoci &&  (pRelEl->StartLoci - RefStartLoci) >= pProcParams->MinLength)
				OutputRefUniques(++NumElsOutput,pszElType,pszRefSpecies,pszRelSpecies,pszRefChrom,RefStartLoci,pRelEl->StartLoci-1,pProcParams);
			
			RefStartLoci = pRelEl->StartLoci + pRelEl->Len;
			RefLen = RefEndLoci - RefStartLoci;
			if(RefLen < 0)
				RefLen = 0;
			else
				RefLen += 1;
			}

		if(RefLen > 0 && RefLen  >= pProcParams->MinLength)
			OutputRefUniques(++NumElsOutput,pszElType,pszRefSpecies,pszRelSpecies,pszRefChrom,RefStartLoci,RefEndLoci,pProcParams);
		}
	}
return(NumElsOutput);
}

bool
OutputRefUniques(int CoreID,char *pszElType,char *pszSpecies,char *pszSpeciesList,char *pszChrom,int StartLoci, int EndLoci,tsProcParams *pProcParams)
{
char szLineBuff[2058];
int Len;
if(pProcParams->hLociFile <= 0)
	return(false);

Len = sprintf(szLineBuff,"%d,\"%s\",\"%s\",\"%s\",%d,%d,%d,\"%s\",%d\n",
		CoreID,pszElType,
		pszSpecies,
		pszChrom,StartLoci,EndLoci,EndLoci-StartLoci+1,
		pszSpeciesList,0);

CUtility::SafeWrite(pProcParams->hLociFile,szLineBuff,Len);
return(true);
}


bool	
OutputDropoutResults(int *pCntStepCnts,tsProcParams *pProcParams)
{
tsLenRangeClass *pRange;
static bool bOutputHdrFirst = true;
char szLineBuff[2058];
int Idx;
int Steps;
int Instances;
int Len;
int *pStep;
pStep = pCntStepCnts;

if(pProcParams->hStatsFile <= 0)
	return(false);

if(bOutputHdrFirst)
	{
	bOutputHdrFirst = false;
	Len = sprintf(szLineBuff,"\"LenRange\",\"Mismatches\",\"TotInstances\",\"IG\",\"US\",\"5'UTR\",\"CDS\",\"INTRON\",\"3'UTR\",\"DS\",\"5'ExSplice\",\"3'ExSplice\"");
	CUtility::SafeWrite(pProcParams->hStatsFile,szLineBuff,Len);
	}
pStep = pCntStepCnts;
for(Idx = 0; Idx < cLenRanges; Idx++, pStep += 7)
	{
	for(Instances = Steps = 0; Steps < 7; Steps++)
		Instances += pStep[Steps];
	pRange = GetRangeClass(Idx+1);
	Len = sprintf(szLineBuff,"\n\"%s\",%d,%d",
						pRange->pszDescr,0,Instances);
	for(Steps = 0; Steps < 7; Steps++)
			Len += sprintf(&szLineBuff[Len],",%d",pStep[Steps]);
	Len += sprintf(&szLineBuff[Len],",0,0");
	CUtility::SafeWrite(pProcParams->hStatsFile,szLineBuff,Len);
	}
return(true);
}
