// genrollups.cpp : Defines the entry point for the console application.
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


const unsigned int cProgVer = 304;		// increment with each release

const int cRUmaxSpecies = 50;			// max number of species handled

// Caution: check genhyperconserved before changing the following....
const int cMaxNumRangeClasses = 28;		// max number of range classes
const int cMaxMatchDistSegments  = 100;	// max match distribution profile segments (applies if eProcModeOutspecies)
const int cMaxNumIdentityClasses = 101;	// 0-100% identity classes

CStopWatch gStopWatch;
CDiagnostics gDiagnostics;				// for writing diagnostics messages to log file
char gszProcName[_MAX_FNAME];			// process name

// user selects one of following bin class range sets with '-c<rangeclass>'
typedef enum eRCClass {
	eRCCFull = 0,			// 20-29,30-49,50-74,75-99,100-124,125-149,150-174,175-199,200-249,250-299,300-349,350-399,400-449,450-499,500-599,600-699,700-799,800-899,900-999,1000-1249,1250-1499,1500-1749,1750-1999,2000+
	eRCCReduced,			// 20-49,50-99,100-149,150-199,200-249,250-299,300+
	eRCCMinimal,			// 20-49,50-99,100-199,200-299,300+
	eRCCMinimalA,			// 20-99,100-199,200+
	eRCCMinimalB,			// 20-49,50-99,100+
	eRCCMinimalC			// 100+

} etRCClass;

// user selects one of following processing modes with '-m<mode>'
typedef enum eRPMode {
	eRPMTotals = 0,
	eRPMRegional,
	eRPMLociBases,
	eRPMLociRegional,
	eRPMOutspeciesSummary
} etRPMode;

// user selects one of following generated results layout with '-l<layout>'
// only applies if processing mode is eRPMOutspeciesSummary 
typedef enum eRPRsltLayout {
	eRPRsltStandard = 0,
	eRPRsltTable,
	eRPRsltSegDistRows,
	eRPRsltSegDistTables,
	eRPRsltSegDistCols,
	eRPRsltIdentDistCols,
	eRPRsltRegionTableA,
	eRPRsltRegionTable
} etRPRsltLayout;

// processing parameters passsed between functions
typedef struct TAG_sProcParams 
	{
	etRPMode Mode;			// processing mode
	teFuncRegion Region;	// genomic functional region to filter on
	bool bSegsExcludeMatch; // true if cores which 100% match outspecies are to be excluded from segment distributions
	etRCClass RangeClass;	// into which range class characterisation 
	bool bPercentages;		// output as percentages
	int hRsltsFile;			// results file handle
	char *pszRsltsFile;		// results file name
	bool bHeaderOut;		// true if results header has been written
	etRPRsltLayout RsltsLayout; // layout format to use if processing eRPMOutspeciesSummary
	int NumDistSegs;		// number of distribution segments 
	int NumSpecies;			// number of species in szSpecies following
	char szSpecies[cRUmaxSpecies][cMaxDatasetSpeciesChrom];	// species names of interest
	int Align2Core;			// must be at least this number of bases in alignment to count as aligned to a core
	double PCAlign2Core;	// minimum percentage (1..100) aligned to core to count as as aligned to core (default is 0%)");
	double IDAlign2Core;	// minimum identity (1..100) aligned to core to count as as aligned to core (default is 0%)");
	double IDOSIdentity;	// minimum outspecies (matches/matchesPlusmismatches) (1..100) aligned to core to count as as aligned to core (default is 1%)");
} tsProcParams; 

int	InitRangeClasses(etRCClass RangeClass);
char *RangeClass2Txt(etRCClass RangeClass);
int Process(etRPMode Mode,						// processing mode
				etRPRsltLayout RsltsLayout,		// table layout to generate results into
				bool bSegsExcludeMatch,			// true if cores which 100% match outspecies are to be excluded from segment distributions
				etRCClass RangeClass,			// length range class
				teFuncRegion Region,			// process this genomic functional region
				bool bPercentages,				// results as percentages instead of counts
				int Align2Core,					// must be at least this number of bases in alignment to count as aligned to a core
				double PCAlign2Core,			// minimum percentage (1..100) aligned to core to count as as aligned to core (default is 1%)");
				double IDAlign2Core,			// minimum identity (1..100) aligned to core to count as as aligned to core (default is 1%)");
				double dIDOSIdentity,			// minimum outspecies (matches/matchesPlusmismatches) (1..100) aligned to core to count as as aligned to core (default is 1%)");
				char *pszInputFiles,			// process from this/these inputfiles - can contain wildcards
				char *pszOutputFile);			// results into this file

bool ProcessThisFile(char *pszFile,void *pParams);	// will be ptr to tsProcParams
int WriteRsltsHeader(tsProcParams *pParams);
int ParseFileNameSpecies(char *pszFile,	// filename to parse for ref+rel species
				 char *pszClass,		// user supplied buffer to hold element class (NULL allowed)
				 char *pszRefSpecies,	// user supplied buffer to hold ref species name
				 char *pszRelSpecies,	// user supplied buffer to hold rel species names
 				 char *pszType);			// user supplied buffer to hold element type (NULL allowed)
char *Mode2Text(etRPMode Mode);
char *Region2Text(int Region);
int ProcessLengthSummary(char *pszFile,CCSVFile *pCSV,tsProcParams *pParams);
int ProcessLociLengthSummary(char *pszFile,CCSVFile *pCSV,tsProcParams *pParams);
int ProcessRegionalLengthSummary(char *pszFile,CCSVFile *pCSV,tsProcParams *pParams);
int ProcessLociRegionalLengthSummary(char *pszFile,CCSVFile *pCSV,tsProcParams *pParams);
int ProcessOutspeciesSummary(char *pszFile,CCSVFile *pCSV,tsProcParams *pParams);


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

char szOutputFile[_MAX_PATH];
char szInputFile[_MAX_PATH];

etRPMode iMode;
int iRegion;
bool bPercentages;
etRCClass RangeClass;
bool bSegsExcludeMatch;
etRPRsltLayout iRsltsLayout;		// out species processing file format to generate
int iAlign2Core;					// at least this many bases aligned in outspecies for it to count as aligned to core
double dPCAlign2Core;				// minimum percentage (1..100) aligned to core to count as as aligned to core (default is 1%)");
double dIDAlign2Core;				// minimum identity (1..100) aligned to core to count as as aligned to core (default is 1%)");
double dIDOSIdentity;				// minimum outspecies (matches/matchesPlusmismatches) (1..100) aligned to core to count as as aligned to core (default is 1%)");


// command line args
struct arg_lit  *help    = arg_lit0("hH","help",                "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_int *ScreenLogLevel=arg_int0("S", "ScreenLogLevel",	"<int>","Level of diagnostics written to logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"diagnostics log file");

struct arg_file *InFile =  arg_file1("i",NULL,"<file>",			"input from .csv files");
struct arg_file *OutFile = arg_file1("o",NULL,"<file>",			"output to rollup statistics file as CSV");
struct arg_int  *Mode = arg_int0("m","mode","<int>",			"processing mode - 0:Totals,1:Regional,2:Loci bases Totals, 3: Loci bases Regional, 4: Outspecies totals");
struct arg_int  *Region = arg_int0("r","region","<int>",		"which region, 0:IG,1:5'US,2:5'UTR,3:CDS,4:Intronic,5:3'UTR,6:3'DS,7:All (default is 7:All)");
struct arg_lit  *SegsExcludeMatch = arg_lit0("s","segsmatchprof","exclude 100% match cores in segment distribution profiles");
struct arg_lit  *Percentages = arg_lit0("p","percent",			"generate results as percentages");
struct arg_int  *Range = arg_int0("c","binclass","<int>",		"range bin class characterisation - 0: full, 1: reduced, 2: minimal, 3: UCSCUltra 200+");

struct arg_int  *RsltsLayout = arg_int0("l","layout","<int>",	"File layout for results if out species processing - 0: row values, 1: column values, 2: Segment Distribution rows, 3: Segment Distribution tables, 4: Segment Distribution columns, 5: core identity columns, 6: regions");
struct arg_int  *Align2Core = arg_int0("a","align2core","<int>","minimum bases in outspecies for alignment to count as aligned to core (default is 1)");
struct arg_dbl  *PCAlign2Core = arg_dbl0("P","pcalign2core","<dbl>","minimum percentage bases aligned to core to count as as aligned to core (default is 0%)");
struct arg_dbl  *IDAlign2Core = arg_dbl0("A","idalign2core","<dbl>","minimum identity to count as as aligned to core (default is 0%)");
struct arg_dbl  *IDOSIdentity = arg_dbl0("k","osidentity","<dbl>","minimum outspecies identity ,matches/(matches+mismatches), to count as as aligned to core (default is 0%)");

struct arg_end *end = arg_end(20);

void *argtable[] = {help,version,FileLogLevel,ScreenLogLevel,LogFile,
					InFile,OutFile,
					Mode,Region,SegsExcludeMatch,Percentages,Range,RsltsLayout,
					Align2Core,PCAlign2Core,IDAlign2Core,IDOSIdentity,
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

	iMode = (etRPMode)(Mode->count ? Mode->ival[0] : eRPMTotals);
	if(iMode < eRPMTotals || iMode > eRPMOutspeciesSummary)
		{
		printf("\nError: Processing mode '-m%d' is not supported\n",(int)iMode);
		exit(1);
		}
	if(iMode == eRPMRegional || iMode == eRPMLociRegional)
		{
		iRegion = Region->count ? Region->ival[0] : eFRAuto;
		if(iRegion < eFRIntergenic || iRegion > eFRAuto)
			{
			printf("\nError: Region '-r%d' is not supported\n",iRegion);
			exit(1);
			}
		}
	else
		iRegion = eFRAuto;

	if(iMode == eRPMOutspeciesSummary)
		{
		iRsltsLayout = RsltsLayout->count ? (etRPRsltLayout)RsltsLayout->ival[0] : eRPRsltStandard;
		if(iRsltsLayout < eRPRsltStandard || iRsltsLayout > eRPRsltRegionTable)
			{
			printf("\nError: Results layout format '-x%d' is not supported\n",iRsltsLayout);
			exit(1);
			}

		iAlign2Core = Align2Core->count ? Align2Core->ival[0] : 1;
		if(iAlign2Core < 1)
			{
			printf("\nError: Aligned bases to core '-a%d' must be >= 1",iAlign2Core);
			exit(1);
			}

		dPCAlign2Core = PCAlign2Core->count ? PCAlign2Core->dval[0] : 0.0;
		if(dPCAlign2Core < 0.0 || dPCAlign2Core > 100.0)
			{
			printf("\nError: Minimum percentage aligned bases to core '-P%f' must be in range 0.0 to 100.0",dPCAlign2Core);
			exit(1);
			}
		
		dIDAlign2Core = IDAlign2Core->count ? IDAlign2Core->dval[0] : 0.0;
		if(dIDAlign2Core < 0.0 || dIDAlign2Core > 100.0)
			{
			printf("\nError: Minimum identity to core '-A%f' must be in range 0.0 to 100.0",dIDAlign2Core);
			exit(1);
			}

		dIDOSIdentity = IDOSIdentity->count ? IDOSIdentity->dval[0] : 0.0;
		if(dIDOSIdentity < 0.0 || dIDOSIdentity > 100.0)
			{
			printf("\nError: Minimum outspecies identity '-k%f' must be in range 0.0 to 100.0",dIDOSIdentity);
			exit(1);
			}
		}
	else
		{
		iRsltsLayout = eRPRsltStandard;
		iAlign2Core = 1;
		dIDAlign2Core = 0.0;
		dPCAlign2Core = 0.0;
		dIDOSIdentity = 0.0;
		}
	 

	bSegsExcludeMatch = SegsExcludeMatch->count ? true : false;
	RangeClass = Range->count ? (etRCClass)Range->ival[0] : eRCCFull;
	if(RangeClass < eRCCFull || RangeClass > eRCCMinimalC)
		{
		printf("\nError: Range class characterisation '-c%d' is not supported\n",RangeClass);
		exit(1);
		}

	if(!InFile->count)
		{
		printf("\nError: Input files to process '-i<inputfilespec>' is not specified\n");
		exit(1);
		}
	strcpy(szInputFile,InFile->filename[0]);
	if(!OutFile->count)
		{
		printf("\nError: Output file to generate '-o<outputfilespec>' is not specified\n");
		exit(1);
		}
	strcpy(szOutputFile,OutFile->filename[0]);
	bPercentages = Percentages->count ? true : false;

	// now that command parameters have been parsed then initialise diagnostics log system
	if(!gDiagnostics.Open(szLogFile,(etDiagLevel)iScreenLogLevel,(etDiagLevel)iFileLogLevel,true))
		{
		printf("\nError: Unable to start diagnostics subsystem.");
		if(szLogFile[0] != '\0')
			printf(" Most likely cause is that logfile '%s' can't be opened/created",szLogFile);
		exit(1);
		}
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Version: %d.%2.2d Processing parameters:",cProgVer/100,cProgVer%100);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Processing Mode: '%s'",Mode2Text(iMode));
	if(iMode == eRPMOutspeciesSummary)
		{
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Segment distribution profiles to exclude 100%% identity cores: %s",
				bSegsExcludeMatch ? "yes" : "no");
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Generated processed results layout format: %d",iRsltsLayout);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Alignments must have at least : %d bases",iAlign2Core);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Minimum percentage bases must be least : %1.2f%% ",dPCAlign2Core);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Minimum identity at least : %1.2f%%",dIDAlign2Core);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Minimum outspecies identity at least : %1.2f%%",dIDOSIdentity);
		}

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Range Classes: %s",RangeClass2Txt(RangeClass));

	if(iMode == eRPMRegional || iMode == eRPMLociRegional)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Region: %d - '%s'",iMode,Region2Text(iRegion));

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Results as percentages: '%s'",bPercentages ? "yes" : "no");
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Input filespec: '%s'",szInputFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Output file: '%s'",szOutputFile);
	gStopWatch.Start();
#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	Rslt = Process(iMode,iRsltsLayout,bSegsExcludeMatch,RangeClass,(teFuncRegion)iRegion,bPercentages,iAlign2Core,dPCAlign2Core,dIDAlign2Core,dIDOSIdentity,szInputFile,szOutputFile);
	gStopWatch.Stop();
	Rslt = Rslt < 0 ? 1 : 0;
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Exit Code: %d Total processing time: %s",Rslt,gStopWatch.Read());
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
RangeClass2Txt(etRCClass RangeClass)
{
switch(RangeClass) {
	case eRCCFull:     return((char *)"20-29,30-49,50-74,75-99,100-124,125-149,150-174,175-199,200-249,250-299,300-349,350-399,400-449,450-499,500-599,600-699,700-799,800-899,900-999,1000-1249,1250-1499,1500-1749,1750-1999,2000+");
	case eRCCReduced:  return((char *)"20-49,50-99,100-149,150-199,200-249,250-299,300+");
	case eRCCMinimal:  return((char *)"20-49,50-99,100-199,200-299,300+");
	case eRCCMinimalA: return((char *)"20-99,100-199,200+");
	case eRCCMinimalB: return((char *)"20-49,50-99,100+");
	case eRCCMinimalC: return((char *)"100+");
	default: return((char *)"Undefined");
	}
}

int Process(etRPMode Mode,						// processing mode
				etRPRsltLayout RsltsLayout,		// table layout to generate results into
				bool bSegsExcludeMatch,			// true if cores which 100% match outspecies are to be excluded from segment distributions
				etRCClass RangeClass,			// length range class
				teFuncRegion Region,			// process this genomic functional region
				bool bPercentages,				// results as percentages instead of counts
				int Align2Core,					// must be at least this many bases in alignment to count as an alignment
				double PCAlign2Core,			// minimum percentage (0..100) aligned to core to count as as aligned to core (default is 0%)");
				double IDAlign2Core,			// minimum identity (0..100) aligned to core to count as as aligned to core (default is 0%)");
				double IDOSIdentity,			// minimum outspecies (matches/matchesPlusmismatches) (1..100) aligned to core to count as as aligned to core (default is 0%)");
				char *pszInputFiles,			// process from this/these inputfiles - can contain wildcards
				char *pszOutputFile)			// results into this file
{
int Rslt;
char szDirPath[_MAX_PATH];
char szFileSpec[_MAX_PATH];

#ifndef _WIN32
struct stat FileStat;
#endif

tsProcParams ProcParams;
memset(&ProcParams,0,sizeof(tsProcParams));
ProcParams.Mode = Mode;
ProcParams.RsltsLayout = RsltsLayout;
ProcParams.bSegsExcludeMatch = bSegsExcludeMatch;
ProcParams.RangeClass = RangeClass;
ProcParams.pszRsltsFile = pszOutputFile;
ProcParams.bPercentages = bPercentages;
ProcParams.Align2Core = Align2Core;
ProcParams.PCAlign2Core = PCAlign2Core;
ProcParams.IDAlign2Core= IDAlign2Core;
ProcParams.IDOSIdentity = IDOSIdentity;

#ifdef _WIN32
char szDrive[_MAX_DRIVE];
char szDir[_MAX_DIR];
char szFname[_MAX_FNAME];
char szExt[_MAX_EXT];
_splitpath(pszInputFiles,szDrive,szDir,szFname,szExt);
_makepath(szDirPath,szDrive,szDir,"","");
if(!szFname[0])
	strcpy(szFname,"*");
if(!szExt[0])
	strcpy(szExt,".*");
sprintf(szFileSpec,"%s%s",szFname,szExt);
#else
CUtility::splitpath(pszInputFiles,szDirPath,szFileSpec);
if(!szFileSpec[0])
	strcpy(szFileSpec,"*");
#endif


#ifdef _WIN32
if((ProcParams.hRsltsFile = open(pszOutputFile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
if((ProcParams.hRsltsFile = open(pszOutputFile,O_RDWR | O_CREAT |O_TRUNC, S_IREAD | S_IWRITE ))==-1)
#endif
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create %s - %s",pszOutputFile,strerror(errno));
	return(eBSFerrCreateFile);
	}

ProcParams.Region = Region;
if((Mode == eRPMRegional || Mode == eRPMLociRegional) && Region == eFRAuto)				// process all regions
	{
	Rslt = eBSFSuccess;
	for(int iRegion = eFRIntergenic; iRegion <= eFRDnstream; iRegion++)
		{
		ProcParams.Region = (teFuncRegion)iRegion;

		CSimpleGlob glob(SG_GLOB_FULLSORT);
		if (glob.Add(pszInputFiles) >= SG_SUCCESS)
			{
			for (int n = 0; Rslt >= eBSFSuccess &&  n < glob.FileCount(); ++n)
				gDiagnostics.DiagOut(eDLInfo,gszProcName,"Will process in this order: %d '%s", n+1,glob.File(n));

			for (int n = 0; Rslt >= eBSFSuccess &&  n < glob.FileCount(); ++n)
				Rslt = ProcessThisFile(glob.File(n),&ProcParams);
			}
		else
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to glob '%s",pszInputFiles);
			Rslt = eBSFerrOpnFile;	// treat as though unable to open file
			}

		if(Rslt < eBSFSuccess)
			break;
		}
	}
else
	{
	Rslt = eBSFSuccess;
	CSimpleGlob glob(SG_GLOB_FULLSORT);
	if (glob.Add(pszInputFiles) >= SG_SUCCESS)
		{
		for (int n = 0; Rslt >= eBSFSuccess &&  n < glob.FileCount(); ++n)
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"Will process in this order: %d '%s", n+1,glob.File(n));
		for (int n = 0; Rslt >= eBSFSuccess &&  n < glob.FileCount(); ++n)
			Rslt = ProcessThisFile(glob.File(n),&ProcParams);
		}
	else
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to glob '%s",pszInputFiles);
		Rslt = eBSFerrOpnFile;	// treat as though unable to open file
		}
	}
if(ProcParams.hRsltsFile != -1)
	close(ProcParams.hRsltsFile);
return(Rslt);
}

bool 
ProcessThisFile(char *pszFile,void *pParams)
{
tsProcParams *pProcParams = (tsProcParams *)pParams;
CCSVFile *pCSV = new CCSVFile;
int Rslt;

InitRangeClasses(pProcParams->RangeClass);

pCSV->SetMaxFields(15 + (cMaxMatchDistSegments*4));
if((Rslt=pCSV->Open(pszFile))!=eBSFSuccess)
	{
	while(pCSV->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,pCSV->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open file: %s",pszFile);
	return(false);
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing File: %s",pszFile);


switch(pProcParams->Mode) {
	case eRPMTotals:			// process length summary files
		Rslt = ProcessLengthSummary(pszFile,pCSV,pProcParams);
		break;
	case eRPMRegional:
		Rslt = ProcessRegionalLengthSummary(pszFile,pCSV,pProcParams);
		break;
	case eRPMLociBases:			// process loci length summary files
		Rslt = ProcessLociLengthSummary(pszFile,pCSV,pProcParams);
		break;
	case eRPMLociRegional:
		Rslt = ProcessLociRegionalLengthSummary(pszFile,pCSV,pProcParams);
		break;
	case eRPMOutspeciesSummary:
		Rslt = ProcessOutspeciesSummary(pszFile,pCSV,pProcParams);
		break;
	}
delete pCSV;
return(Rslt >= eBSFSuccess ? true : false);
}


// ParseFileNameSpecies
// parse filename to determine reference and relative species --- UGH!!!
// Expected filename format is -
// <Class><RefSpecies><RelSpecies>_<Type>*.csv
// <Class> 'L' or 'S' but could be some other single char
// <RefSpecies> species name with the last one or two chars digits e.g. hg17
// <RelSpecies> concatenated species names, last species name is terminated by an underscore '_'
// <Type>	currently 'u' for ultra or 'h' for hyper
int
ParseFileNameSpecies(char *pszFile,		// filename to parse for ref+rel species
				 char *pszClass,		// user supplied buffer to hold element class
				 char *pszRefSpecies,	// user supplied buffer to hold ref species name
				 char *pszRelSpecies,	// user supplied buffer to hold rel species names
 				 char *pszType)			// user supplied buffer to hold element type
{
char *pSrcChr;
char *pDstChr;
char Chr;
bool bParseRef;
bool bParsedRef;
bool bParseRel;
bool bParsedRel;
bool bDigit;

char szFname[256];
char szDir[256];

#ifdef _WIN32
char szDrive[_MAX_DRIVE];
char szExt[_MAX_EXT];
_splitpath(pszFile,szDrive,szDir,szFname,szExt);
#else
CUtility::splitpath(pszFile,szDir,szFname);
#endif

if(pszClass != NULL)
	switch(szFname[0]) {
		case 's': case 'S':
			strcpy(pszClass,"summary");
			break;
		case 'l': case 'L':
			strcpy(pszClass,"loci");
			break;
		default:
			*pszClass = szFname[0];
			pszClass[1] = '\0';
			break;
		}

pSrcChr = &szFname[1];
pDstChr = pszRefSpecies;
bParseRef = true;
bParsedRef = false;
bParseRel = false;
bParsedRel = false;
bDigit = false;
*pszRefSpecies = '\0';
*pszRelSpecies = '\0';
while(Chr = *pSrcChr++)
	{
	if(Chr == '_')
		{
		bParsedRel = true;
		*pDstChr = '\0';
		if(pszType != NULL)
			{
			if((Chr = *pSrcChr) == '\0')
				{
				*pszType = '?';
				pszType[1] = '\0';
				}
			else
				switch(Chr) {
					case 'u': case 'U':
						strcpy(pszType,"ultra");
						break;
					case 'h': case 'H':
						strcpy(pszType,"hyper");
						break;
					default:
						*pszType = Chr;
						pszType[1] = '\0';
						break;
					}
			}
		break;
		}
	if(bParseRef)
		{
		if(Chr >= '0' && Chr <= '9')
			bDigit = true;
		else
			{
			if(bDigit)		// was previously parsing digits but not now then that means rel species starting
				{
				*pDstChr = '\0';
				bParseRef = false;
				bParsedRef = true;
				bParseRel = true;
				pDstChr = pszRelSpecies;
				}
			}
		}
	*pDstChr++ = Chr;
	}
if(!(bParsedRef && bParsedRel))
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to parse filename '%s' into ref+rel species names",szFname);
	return(eBSFerrFileName);
	}
return(eBSFSuccess);
}


// ParseNumSpecies
// Initialises pProcParams with parsed species names in space or comma delimited list ptd at by pszSpeciesList
// Returns number of species
int
ParseNumSpecies(char *pszSpeciesList,tsProcParams *pProcParams)
{
// parse out species list
char Chr;
char *pSpecies;
int NumSpecies = 0;
bool InToken = false;
if(pszSpeciesList == NULL || *pszSpeciesList == '\0')
	return(0);

while(Chr = *pszSpeciesList++)
	{
	if(Chr == '"' || Chr == '\'') // change any single or double quotes into spaces
		Chr = ' ';
	if(isspace(Chr) || Chr==',')
		{
		if(!InToken)			// slough whitespace or ',' if not inside a token parse
			continue;
		InToken = false;
		pszSpeciesList[-1] = '\0';
		if(pProcParams != NULL)
			{
			strncpy(pProcParams->szSpecies[NumSpecies],pSpecies,cMaxDatasetSpeciesChrom);
			pProcParams->szSpecies[NumSpecies][cMaxDatasetSpeciesChrom-1] = '\0';
			}
		pszSpeciesList[-1] = Chr;
		NumSpecies++;
		if(NumSpecies >= cRUmaxSpecies)
			break;
		continue;
		}
	if(!InToken)			// if not already inside token then start token 
		{
		pSpecies = pszSpeciesList-1;
		InToken = true;
		}
	}
if(InToken)
	{
	pszSpeciesList[-1] = '\0';
	if(pProcParams != NULL)
		{
		strncpy(pProcParams->szSpecies[NumSpecies],pSpecies,cMaxDatasetSpeciesChrom);
		pProcParams->szSpecies[NumSpecies][cMaxDatasetSpeciesChrom-1] = '\0';
		}
	pszSpeciesList[-1] = Chr;
	NumSpecies++;
	}
return(NumSpecies);
}



char *
Mode2Text(etRPMode Mode)
{
switch(Mode) {
	case eRPMTotals:
		return((char *)"Totals summary");
	case eRPMRegional:
		return((char *)"Regional summary");
	case eRPMLociBases:
		return((char *)"Loci bases Totals summary");
	case eRPMLociRegional:
		return((char *)"Loci bases Regional summary");
	case eRPMOutspeciesSummary:
		return((char *)"Out species summary");
	default:
		break;
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Unsupported processing mode '%d' requested",(int)Mode);
return((char *)"Unsupported");
}


// Regions can be overlapping
// feature bits - note that features can concurrently overlap multiple gene components
// bit 0 (0x01) part of feature overlaps CDS
// bit 1 (0x02) part of feature overlaps 5'UTR
// bit 2 (0x04) part of feature overlaps 3'UTR
// bit 3 (0x08) part of feature overlaps Intron
// bit 4 (0x10) part of feature overlaps 5'upstream regulatory 	
// bit 5 (0x20) part of feature overlaps 3'downstream regulatory
// bit 6 (x040) part of feature overlaps intron5'/3'exon splice site
// bit 7 (0x80) part of feature overlaps exon5'/3'intron splice site
// Normalise regions maps regions into a logical index
// Intergenic	-> 0
// US			-> 1
// 5'UTR		-> 2
// CDS			-> 3
// INTRON		-> 4
// 3'UTR		-> 5
// DS			-> 6
int
NormaliseRegion(int Region)
{
if(Region & cFeatBitCDS)
	return(eFRCDS);			// CDS
if(Region & cFeatBit5UTR)
	return(eFR5UTR);		// 5'UTR
if(Region & cFeatBit3UTR)
	return(eFR3UTR);		// 3'UTR
if(Region & cFeatBitIntrons)
	return(eFRIntronic);	// Intron
if(Region & cFeatBitUpstream)
	return(eFRUpstream);	// 5'US
if(Region & cFeatBitDnstream)
	return(eFRDnstream);	// 3'DS
return(eFRIntergenic);		// IG
}

char *
Region2Text(int Region)
{
switch(Region) {
	case 0:
		return((char *)"IG");
	case 1:
		return((char *)"US");
	case 2:
		return((char *)"5'UTR");
	case 3:
		return((char *)"CDS");
	case 4:
		return((char *)"INTRON");
	case 5:
		return((char *)"3'UTR");
	case 6:
		return((char *)"DS");
	case 7:
		return((char *)"ANY");
	default:
		break;
	}
gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unsupported region '%d' requested",Region);
return((char *)"Unsupported");
}

int 
WriteRsltsHeader(tsProcParams *pParams)
{
char szLineBuff[2048];
int Len;
int SegIdx;

if(pParams->bHeaderOut)					// if header already written?
	return(eBSFSuccess);

if(pParams->hRsltsFile == -1)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Results file '%s' was closed when called to write headers",pParams->pszRsltsFile);
	return(eBSFerrFileClosed);
	}
switch(pParams->Mode) {
	case eRPMTotals: case eRPMLociBases:
		Len = sprintf(szLineBuff,"\"ElType\",\"Ref Species\",\"Rel Species\",\"Total\"");
		break;
	case eRPMRegional: case eRPMLociRegional:
		Len = sprintf(szLineBuff,"\"ElType\",\"Ref Species\",\"Rel Species\",\"Region\",\"Total\"");
		break;

	case eRPMOutspeciesSummary:
		switch(pParams->RsltsLayout) {
			case eRPRsltStandard:
				Len = sprintf(szLineBuff,"\"ElType\",\"RefSpecies\",\"CoreSpecies\",\"OutSpecies\",\"Region\",\"Qualifier\",\"Name\",\"Value\"\n");
				break;
			
			case eRPRsltTable:
				Len = sprintf(szLineBuff,"\"Core Set\",\"Out Species\",\"Ref Core Len\",\"Ref Core Cnt\",\"Aligned Cnt\",\"Aligned Core %%\",\"Aligned Bases\",\"Aligned Mismatches\",\"Aligned Identical\",\"Aligned Identity\"\n");
				break;

			case eRPRsltSegDistRows:
				Len = sprintf(szLineBuff,"\"ElType\",\"RefSpecies\",\"CoreSpecies\",\"OutSpecies\",\"LengthRange\"");
				for(SegIdx = 0; SegIdx < pParams->NumDistSegs; SegIdx++)
					Len += sprintf(&szLineBuff[Len],",\"Seg%dMatches\",\"Seg%dMismatches\",\"Seg%dInDels\",\"Seg%dUnaligns\"",
								SegIdx+1,SegIdx+1,SegIdx+1,SegIdx+1);
				Len += sprintf(&szLineBuff[Len],"\n");
				break;

			case eRPRsltSegDistTables:
				Len = sprintf(szLineBuff,"\"ElType\",\"RefSpecies\",\"CoreSpecies\",\"OutSpecies\",\"LengthRange\",\"Category\"");
				for(SegIdx = 1; SegIdx <= pParams->NumDistSegs; SegIdx++)
					Len += sprintf(&szLineBuff[Len],",\"Segment%d\"",SegIdx);
				Len += sprintf(&szLineBuff[Len],"\n");
				break;

			case eRPRsltSegDistCols:
				Len = sprintf(szLineBuff,"\"ElType\",\"RefSpecies\",\"CoreSpecies\",\"OutSpecies\",\"LengthRange\",\"Segment\",\"Matches\",\"Mismatches\",\"InDels\",\"Unaligned\"\n");
				break;

			case eRPRsltIdentDistCols:
				Len = sprintf(szLineBuff,"\"ElType\",\"RefSpecies\",\"CoreSpecies\",\"OutSpecies\",\"LengthRange\",\"Core Instances\",\"OG Instances\",\"Identity\",\"OG Identity Instances\"\n");
				break;

			case eRPRsltRegionTable:
				Len = sprintf(szLineBuff,"\"Core Set\",\"Out Species\",\"Ref Core Len\",\"Ref Core Cnt\",\"Aligned Cnt\",\"IG\",\"US\",\"5'UTR\",\"CDS\",\"INTRON\",\"3'UTR\",\"DS\",\"5'ExSplice\",\"3'ExSplice\"\n");
				break;

			case eRPRsltRegionTableA:
				Len = sprintf(szLineBuff,"\"Core Set\",\"Out Species\",\"Ref Core Len\",\"Ref Core Cnt\",\"Aligned Cnt\",\"IG\",\"US\",\"5'UTR\",\"CDS\",\"INTRON\",\"3'UTR\",\"DS\",\"5'ExSplice\",\"3'ExSplice\",\"Excl IG\",\"Excl US\",\"Excl 5'UTR\",\"Excl CDS\",\"Excl INTRON\",\"Excl 3'UTR\",\"Excl DS\",\"Splice Sites\"\n");
				break;

			default:
				Len = sprintf(szLineBuff,",\"Illegal Layout\"\n");
				break;
			}
		break;

	default:
		Len = sprintf(szLineBuff,",\"Illegal Mode\"");

	}

if(pParams->Mode != eRPMOutspeciesSummary)
	switch(pParams->RangeClass) {
		case eRCCFull:
			Len += sprintf(&szLineBuff[Len],",\"20-29\",\"30-49\",\"50-74\",\"75-99\",\"100-124\",\"125-149\",\"150-174\",\"175-199\"");
			Len += sprintf(&szLineBuff[Len],",\"200-249\",\"250-299\",\"300-349\",\"350-399\",\"400-449\",\"450-499\"");
			Len += sprintf(&szLineBuff[Len],",\"500-599\",\"600-699\",\"700-799\",\"800-899\",\"900-999\"");
			Len += sprintf(&szLineBuff[Len],",\"1000-1249\",\"1250-1499\",\"1500-1749\",\"1750-1999\",\"2000+\"");
			break;

		case eRCCReduced:
			Len += sprintf(&szLineBuff[Len],",\"20-49\",\"50-99\",\"100-149\",\"150-199\",\"200-249\",\"250-299\",\"300+\"");
			break;

		case eRCCMinimal:
			Len += sprintf(&szLineBuff[Len],",\"20-49\",\"50-99\",\"100-199\",\"200-299\",\"300+\"");
			break;

		case eRCCMinimalA:
			Len += sprintf(&szLineBuff[Len],",\"20-99\",\"100-199\",\"200+\"");
			break;

		case eRCCMinimalB:
			Len += sprintf(&szLineBuff[Len],",\"20-49\",\"59-99\",\"100+\"");
			break;

		case eRCCMinimalC:
			Len += sprintf(&szLineBuff[Len],",\"100+\"");
			break;

		default:
			Len += sprintf(&szLineBuff[Len],",\"Illegal Range Class\"");
			break;
		}	

if(Len && (write(pParams->hRsltsFile,szLineBuff,Len)!=Len))
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to write results file '%s' - '%s'",pParams->pszRsltsFile,strerror(errno));
	return(eBSFerrWrite);
	}
pParams->bHeaderOut = true;
return(eBSFSuccess);
}


int 
ProcessLengthSummary(char *pszFile,CCSVFile *pCSV,tsProcParams *pParams)
{
int Rslt;
int NumRows;
int NumFields;
int Len;
char *pszRange;
int RangeCnt;
int Total;
int Totals[cMaxNumRangeClasses];			// to hold all totals
int Idx;
char szLineBuff[2048];
char szRefSpecies[200];
char szRelSpecies[200];
char szElementType[20];


if((Rslt=ParseFileNameSpecies(pszFile,NULL,szRefSpecies,szRelSpecies,szElementType))!=eBSFSuccess)
	return(Rslt);

// Skip first row as that contains titles
if((Rslt=pCSV->NextLine())<12)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expected at least 12 header fields in 1st row of '%s', NextLine returned '%d'",pszFile,Rslt);
	return(eBSFerrFieldCnt);
	}

// subsequent rows contain fields of interest
memset(Totals,0,sizeof(Totals));
NumRows = 0;
Total = 0;
while((Rslt=pCSV->NextLine()) == 12)	// onto next line containing fields
	{
	if(NumRows == cMaxNumRangeClasses)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expected only 28 rows (plus header) in '%s', file has more rows",pszFile);
		return(eBSFerrRowCnt);
		}
	NumFields = pCSV->GetCurFields();
	if(NumFields != 12)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expected 12 fields in '%s', GetCurFields() returned '%d'",pszFile,NumFields);
		return(eBSFerrFieldCnt);
		}
	pCSV->GetText(1,&pszRange);
	pCSV->GetInt(3,&RangeCnt);

	switch(pParams->RangeClass) {
		case eRCCFull:
			Totals[NumRows] = RangeCnt;
			break;
		case eRCCReduced:
			switch(NumRows) {
				case 0:						// 0-4
					Totals[0] = RangeCnt;
					break;
				case 1:						// 5-9
					Totals[0] += RangeCnt;
					break;
				case 2:						// 10-14
					Totals[1] = RangeCnt;
					break;
				case 3:						// 15-19
					Totals[1] += RangeCnt;
					break;
				case 4:						// 20-29
					Totals[2] = RangeCnt;
					break;
				case 5:						// 30-49
					Totals[2] += RangeCnt;
					break;
				case 6:						// 50-74
					Totals[3] = RangeCnt;
					break;
				case 7:						// 75-99
					Totals[3] += RangeCnt;
					break;
				case 8:						// 100-124
					Totals[4] = RangeCnt;
					break;
				case 9:						// 125-149
					Totals[4] += RangeCnt;
					break;
				case 10:					// 150-174
					Totals[5] = RangeCnt;
					break;
				case 11:					// 175-199
					Totals[5] += RangeCnt;
					break;
				case 12:					// 200-249
					Totals[6] = RangeCnt;
					break;
				case 13:					// 250-300
					Totals[7] = RangeCnt;
					Totals[8] = 0;
					break;
				default:
					Totals[8] += RangeCnt;
					break;
				}
			break;
		case eRCCMinimal:
			switch(NumRows) {
				case 0:						// 0-4
					Totals[0] += RangeCnt;
					break;
				case 1: case 2: case 3:		// 5-9, 10-14,15-19
					Totals[0] += RangeCnt;
					break;
				case 4:						// 20-29
					Totals[1] += RangeCnt;
					break;
				case 5:						// 30-49
					Totals[1] += RangeCnt;
					break;
				case 6:						// 50-74
					Totals[2] += RangeCnt;
					break;
				case 7:						// 75-99
					Totals[2] += RangeCnt;
					break;
				case 8:						// 100-124
					Totals[3] += RangeCnt;
					break;
				case 9:	case 10: case 11:	// 125-149, 150-174,175-199
					Totals[3] += RangeCnt;
					break;
				case 12:					// 200-249
					Totals[4] = RangeCnt;
					break;
				case 13:					// 250-300
					Totals[4] += RangeCnt;
					break;
				default:
					Totals[5] += RangeCnt;
					break;
				}
			break;

		case eRCCMinimalA:
			switch(NumRows) {
				case 0:	case 1: case 2: case 3:		// 0-20
					Totals[0] += RangeCnt;
					break;
				case 4:	case 5: case 6: case 7:	// 20-99
					Totals[1] += RangeCnt;
					break;
				case 8:	case 9:	case 10: case 11: // 100-199
					Totals[2] += RangeCnt;
					break;
				default:
					Totals[3] += RangeCnt;
					break;
				}
			break;

		case eRCCMinimalB:
			switch(NumRows) {
				case 0:	case 1: case 2: case 3:		// 0-20
					Totals[0] += RangeCnt;
					break;
				case 4:	case 5: 					// 20-49
					Totals[1] += RangeCnt;
					break;
				case 6: case 7:						// 50-99
					Totals[2] += RangeCnt;
					break;
				default:							// 100+
					Totals[3] += RangeCnt;
					break;
				}
			break;

		case eRCCMinimalC:
			switch(NumRows) {
				case 0:	case 1: case 2: case 3:		// 0-20
				case 4:	case 5: 					// 20-49
				case 6: case 7:						// 50-99
					break;
				default:							// 100+
					Totals[0] += RangeCnt;
					break;
				}
			break;
		}
	NumRows++;
	}
if(NumRows != cMaxNumRangeClasses)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expected 28 rows (plus header) in '%s', only '%d' parsed",pszFile,NumRows);
	return(eBSFerrRowCnt);
	}

switch(pParams->RangeClass) {
	case eRCCFull:
		for(Idx = 4; Idx < cMaxNumRangeClasses; Idx++)
			Total += Totals[Idx];
		Len=sprintf(szLineBuff,"\n\"%s\",\"%s\",\"%s\",%d",szElementType,szRefSpecies,szRelSpecies,Total);
		for(Idx = 4; Idx < cMaxNumRangeClasses; Idx++)
			if(pParams->bPercentages)
				Len+=sprintf(&szLineBuff[Len],",%1.2f",Total > 0 ? ((double)Totals[Idx]*100.0)/(double)Total : 0.0);
			else
				Len+=sprintf(&szLineBuff[Len],",%d",Totals[Idx]);
		break;

	case eRCCReduced:
		for(Idx = 2; Idx < 9; Idx++)
			Total += Totals[Idx];
		Len=sprintf(szLineBuff,"\n\"%s\",\"%s\",\"%s\",%d",szElementType,szRefSpecies,szRelSpecies,Total);
		for(Idx = 2; Idx < 9; Idx++)
			if(pParams->bPercentages)
				Len+=sprintf(&szLineBuff[Len],",%1.2f",Total > 0 ? ((double)Totals[Idx]*100.0)/(double)Total : 0.0);
			else
				Len+=sprintf(&szLineBuff[Len],",%d",Totals[Idx]);
		break;

	case eRCCMinimal:
		for(Idx = 1; Idx < 6; Idx++)
			Total += Totals[Idx];
		Len=sprintf(szLineBuff,"\n\"%s\",\"%s\",\"%s\",%d",szElementType,szRefSpecies,szRelSpecies,Total);
		for(Idx = 1; Idx < 6; Idx++)
			if(pParams->bPercentages)
				Len+=sprintf(&szLineBuff[Len],",%1.2f",Total > 0 ? ((double)Totals[Idx]*100.0)/(double)Total : 0.0);
			else
				Len+=sprintf(&szLineBuff[Len],",%d",Totals[Idx]);
		break;

	case eRCCMinimalA:
		for(Idx = 1; Idx < 4; Idx++)
			Total += Totals[Idx];
		Len=sprintf(szLineBuff,"\n\"%s\",\"%s\",\"%s\",%d",szElementType,szRefSpecies,szRelSpecies,Total);
		for(Idx = 1; Idx < 4; Idx++)
			if(pParams->bPercentages)
				Len+=sprintf(&szLineBuff[Len],",%1.2f",Total > 0 ? ((double)Totals[Idx]*100.0)/(double)Total : 0.0);
			else
				Len+=sprintf(&szLineBuff[Len],",%d",Totals[Idx]);
		break;

	case eRCCMinimalB:
		for(Idx = 1; Idx < 4; Idx++)
			Total += Totals[Idx];
		Len=sprintf(szLineBuff,"\n\"%s\",\"%s\",\"%s\",%d",szElementType,szRefSpecies,szRelSpecies,Total);
		for(Idx = 1; Idx < 4; Idx++)
			if(pParams->bPercentages)
				Len+=sprintf(&szLineBuff[Len],",%1.2f",Total > 0 ? ((double)Totals[Idx]*100.0)/(double)Total : 0.0);
			else
				Len+=sprintf(&szLineBuff[Len],",%d",Totals[Idx]);
		break;

	case eRCCMinimalC:
		Total += Totals[0];
		Len=sprintf(szLineBuff,"\n\"%s\",\"%s\",\"%s\",%d",szElementType,szRefSpecies,szRelSpecies,Total);
		if(pParams->bPercentages)
			Len+=sprintf(&szLineBuff[Len],",%1.2f",Total > 0 ? ((double)Totals[0]*100.0)/(double)Total : 0.0);
		else
			Len+=sprintf(&szLineBuff[Len],",%d",Totals[0]);
		break;

	default:
		Len=sprintf(szLineBuff,",-1");
		break;
	}

WriteRsltsHeader(pParams);

if(write(pParams->hRsltsFile,szLineBuff,Len)!=Len)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to write results file '%s' - '%s'",pParams->pszRsltsFile,strerror(errno));
	return(eBSFerrWrite);
	}

return(eBSFSuccess);
}


typedef struct TAG_sLenRangeClass {
	int ID;					// uniquely identifies this range
	int Min;				// minimum length in this range
	int Max;				// maximum length in this range
	const char *pszDescr;			// descriptive text
	}tsLenRangeClass;

// length range classes
tsLenRangeClass LenRangeClassesFull[] = {
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
const int cLenRangesFull = sizeof(LenRangeClassesFull)/sizeof(tsLenRangeClass);		  // number of length range classes


tsLenRangeClass LenRangeClassesReduced[] = {
	{1,0,9,"0-9"},
	{2,10,19,"10-19"},
	{3,20,49,"20-49"},
    {4,50,99,"50-99"},
	{5,100,149,"100-149"},
	{6,150,199,"150-199"},
	{7,200,249,"200-249"},
	{8,250,299,"250-299"},
	{9,300,INT_MAX,"300+"}
};
const int cLenRangesReduced = sizeof(LenRangeClassesReduced)/sizeof(tsLenRangeClass);		  // number of length range classes


tsLenRangeClass LenRangeClassesMinimal[] = {
	{1,0,19,"0-19"},
	{2,20,49,"20-49"},
    {3,50,99,"50-99"},
	{4,100,199,"100-199"},
	{5,200,299,"200-299"},
	{6,300,INT_MAX,"300+"}
};
const int cLenRangesMinimal = sizeof(LenRangeClassesMinimal)/sizeof(tsLenRangeClass);	

tsLenRangeClass LenRangeClassesMinimalA[] = {
	{1,0,19,"0-19"},
	{2,20,99,"20-99"},
    {3,100,199,"100-199"},
	{4,200,INT_MAX,"200+"}
};
const int cLenRangesMinimalA = sizeof(LenRangeClassesMinimalA)/sizeof(tsLenRangeClass);	

tsLenRangeClass LenRangeClassesMinimalB[] = {
	{1,0,19,"0-19"},
	{2,20,49,"20-49"},
    {3,50,99,"50-99"},
	{4,100,INT_MAX,"100+"}
};
const int cLenRangesMinimalB = sizeof(LenRangeClassesMinimalB)/sizeof(tsLenRangeClass);	

tsLenRangeClass LenRangeClassesMinimalC[] = {
	{1,100,INT_MAX,"100+"}
};

const int cLenRangesMinimalC = sizeof(LenRangeClassesMinimalC)/sizeof(tsLenRangeClass);

tsLenRangeClass *pLenRangeClasses = LenRangeClassesFull;
int NumLenRanges = cLenRangesFull;

int					// returns number of ranges
InitRangeClasses(etRCClass RangeClass)
{
switch(RangeClass) {
	case eRCCFull:			// 20-29,30-49,50-74,75-99,100-124,125-149,150-174,175-199,200-249,250-299,300-349,350-399,400-449,450-499,500-599,600-699,700-799,800-899,900-999,1000-1249,1250-1499,1500-1749,1750-1999,2000+
		pLenRangeClasses = LenRangeClassesFull;
		NumLenRanges = cLenRangesFull;
		break;

	case eRCCReduced:		// 20-49,50-99,100-149,150-199,200-249,250-299,300+
		pLenRangeClasses = LenRangeClassesReduced;
		NumLenRanges = cLenRangesReduced;
		break;

	case eRCCMinimal:		// 20-49,50-99,100-199,200-299,300+
		pLenRangeClasses = LenRangeClassesMinimal;
		NumLenRanges = cLenRangesMinimal;
		break;

	case eRCCMinimalA:		// 20-99,100-199,200+
		pLenRangeClasses = LenRangeClassesMinimalA;
		NumLenRanges = cLenRangesMinimalA;
		break;

	case eRCCMinimalB:		// 20-49,50-99,100+
		pLenRangeClasses = LenRangeClassesMinimalB;
		NumLenRanges = cLenRangesMinimalB;
		break;

	case eRCCMinimalC:		// 100+
		pLenRangeClasses = LenRangeClassesMinimalC;
		NumLenRanges = cLenRangesMinimalC;
		break;

	default:
		return(-1);
	}
return(NumLenRanges);
}

// GetLengthRangeClass
// Returns ptr to length range class for specified Length, or NULL if can't classify into a range
tsLenRangeClass *GetLengthRangeClass(int Length)
{
int Idx;
tsLenRangeClass *pRange = pLenRangeClasses;
for(Idx = 0; Idx < NumLenRanges; Idx++,pRange++)
	if(Length >= pRange->Min && Length <= pRange->Max)
		return(pRange);
return(NULL);
}


// Expected CSV format
// ElementID	Uniquely identifies this element (1..n)
// CoreType		Currently either 'hypercore' or 'ultracore'
// Chromosome	Chromosome on reference species
// StartLoci	Starting offset (0..n) on Chromosome
// EndLoci		Ending offset (StartLoci + length -1) on Chromosome
// Length		Element length
// SpeciesList	List of all species sharing this element, 1st species is the reference
// Region		In which region core is located
//
// There are an optional 4 fields present if CSV file was generated with outspecies processing (-x2 option)
// Unaligned	Count of bases in outspecies not aligned
// Matches		Count of bases in outspecies which align and match
// Mismatches	Count of bases in outspecies which align but mismatch
// InDels		Count of bases in reference or outspecies which are InDels
int 
ProcessLociLengthSummary(char *pszFile,CCSVFile *pCSV,tsProcParams *pParams)
{
int Rslt;
int NumRows;
int NumFields;
int Len;
int Total;
int Totals[cMaxNumRangeClasses];			// to hold all totals
int Idx;
char szLineBuff[2048];
char szRefSpecies[200];
char szRelSpecies[200];
char szElementType[20];
tsLenRangeClass *pRangeClass;
int SeqLen;


if((Rslt=ParseFileNameSpecies(pszFile,NULL,szRefSpecies,szRelSpecies,szElementType))!=eBSFSuccess)
	return(Rslt);

// no headers, straight into the rows which are expected to contain 9 fields
NumRows = 0;
Total = 0;
memset(Totals,0,sizeof(Totals));

while((Rslt=pCSV->NextLine()) > 0)	// onto next line containing fields
	{
	NumFields = pCSV->GetCurFields();
	if(NumFields < 9)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expected at least 9 fields in '%s', GetCurFields() returned '%d'",pszFile,NumFields);
		return(eBSFerrFieldCnt);
		}
	pCSV->GetInt(7,&SeqLen);
	pRangeClass = GetLengthRangeClass(SeqLen);
	Totals[pRangeClass->ID-1] += SeqLen;
	NumRows++;
	}

switch(pParams->RangeClass) {
	case eRCCFull:
		for(Idx = 4; Idx < cMaxNumRangeClasses; Idx++)
			Total += Totals[Idx];
		Len=sprintf(szLineBuff,"\n\"%s\",\"%s\",\"%s\",%d",szElementType,szRefSpecies,szRelSpecies,Total);
		for(Idx = 4; Idx < cMaxNumRangeClasses; Idx++)
			if(pParams->bPercentages)
				{
					Len+=sprintf(&szLineBuff[Len],",%1.2f",Total > 0 ? ((double)Totals[Idx]*100.0)/(double)Total : 0.0);
				}
			else
				Len+=sprintf(&szLineBuff[Len],",%d",Totals[Idx]);
		break;

	case eRCCReduced:
		for(Idx = 2; Idx < 9; Idx++)
			Total += Totals[Idx];
		Len=sprintf(szLineBuff,"\n\"%s\",\"%s\",\"%s\",%d",szElementType,szRefSpecies,szRelSpecies,Total);
		for(Idx = 2; Idx < 9; Idx++)
		if(pParams->bPercentages)
			{
			Len+=sprintf(&szLineBuff[Len],",%1.2f",Total > 0 ? ((double)Totals[Idx]*100.0)/(double)Total : 0.0);
			}
		else
			Len+=sprintf(&szLineBuff[Len],",%d",Totals[Idx]);
		break;

	case eRCCMinimal:
		for(Idx = 1; Idx < 6; Idx++)
			Total += Totals[Idx];
		Len=sprintf(szLineBuff,"\n\"%s\",\"%s\",\"%s\",%d",szElementType,szRefSpecies,szRelSpecies,Total);
		for(Idx = 1; Idx < 6; Idx++)
			if(pParams->bPercentages)
				{
					Len+=sprintf(&szLineBuff[Len],",%1.2f",Total > 0 ? ((double)Totals[Idx]*100.0)/(double)Total : 0.0);
				}
			else
				Len+=sprintf(&szLineBuff[Len],",%d",Totals[Idx]);		
		break;

	case eRCCMinimalA:
		for(Idx = 1; Idx < 4; Idx++)
			Total += Totals[Idx];
		Len=sprintf(szLineBuff,"\n\"%s\",\"%s\",\"%s\",%d",szElementType,szRefSpecies,szRelSpecies,Total);
		for(Idx = 1; Idx < 4; Idx++)
			if(pParams->bPercentages)
				{
					Len+=sprintf(&szLineBuff[Len],",%1.2f",Total > 0 ? ((double)Totals[Idx]*100.0)/(double)Total : 0.0);
				}
			else
				Len+=sprintf(&szLineBuff[Len],",%d",Totals[Idx]);		
		break;

	case eRCCMinimalB:
		for(Idx = 1; Idx < 4; Idx++)
			Total += Totals[Idx];
		Len=sprintf(szLineBuff,"\n\"%s\",\"%s\",\"%s\",%d",szElementType,szRefSpecies,szRelSpecies,Total);
		for(Idx = 1; Idx < 4; Idx++)
			if(pParams->bPercentages)
				{
					Len+=sprintf(&szLineBuff[Len],",%1.2f",Total > 0 ? ((double)Totals[Idx]*100.0)/(double)Total : 0.0);
				}
			else
				Len+=sprintf(&szLineBuff[Len],",%d",Totals[Idx]);		
		break;

	case eRCCMinimalC:
		Total += Totals[0];
		Len=sprintf(szLineBuff,"\n\"%s\",\"%s\",\"%s\",%d",szElementType,szRefSpecies,szRelSpecies,Total);
		if(pParams->bPercentages)
			{
			Len+=sprintf(&szLineBuff[Len],",%1.2f",Total > 0 ? ((double)Totals[0]*100.0)/(double)Total : 0.0);
			}
		else
			Len+=sprintf(&szLineBuff[Len],",%d",Totals[0]);		
		break;

	default:
		Len=sprintf(szLineBuff,",-1");
		break;
	}

WriteRsltsHeader(pParams);

if(write(pParams->hRsltsFile,szLineBuff,Len)!=Len)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to write results file '%s' - '%s'",pParams->pszRsltsFile,strerror(errno));
	return(eBSFerrWrite);
	}

return(eBSFSuccess);
}



int 
ProcessRegionalLengthSummary(char *pszFile,CCSVFile *pCSV,tsProcParams *pParams)
{
int Rslt;
int NumRows;
int NumFields;
int Len;
char *pszRange;
int RangeCnt;
int Total;
int Totals[cMaxNumRangeClasses];			// to hold all totals
int Idx;
char szLineBuff[2048];
char szRefSpecies[200];
char szRelSpecies[200];
char szClass[20];
char szElementType[20];
char *pszRegion;

pszRegion = Region2Text(pParams->Region);

if((Rslt=ParseFileNameSpecies(pszFile,szClass,szRefSpecies,szRelSpecies,szElementType))!=eBSFSuccess)
	return(Rslt);

// Skip first row as that contains titles
if((Rslt=pCSV->NextLine())<12)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expected at least 12 header fields in 1st row of '%s', NextLine returned '%d'",pszFile,Rslt);
	return(eBSFerrFieldCnt);
	}
// subsequent rows contain fields of interest
memset(Totals,0,sizeof(Totals));
NumRows = 0;
Total = 0;
while((Rslt=pCSV->NextLine()) == 12)	// onto next line containing fields
	{
	if(NumRows == cMaxNumRangeClasses)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expected only 28 rows (plus header) in '%s', file has more rows",pszFile);
		return(eBSFerrRowCnt);
		}
	NumFields = pCSV->GetCurFields();
	if(NumFields < 12)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expected at least 12 fields in '%s', GetCurFields() returned '%d'",pszFile,NumFields);
		return(eBSFerrFieldCnt);
		}
	pCSV->GetText(1,&pszRange);
	pCSV->GetInt(pParams->Region+4,&RangeCnt);

	switch(pParams->RangeClass) {
		case eRCCFull:
			Totals[NumRows] = RangeCnt;
			break;
		case eRCCReduced:
			switch(NumRows) {
				case 0:						// 0-4
					Totals[0] = RangeCnt;
					break;
				case 1:						// 5-9
					Totals[0] += RangeCnt;
					break;
				case 2:						// 10-14
					Totals[1] = RangeCnt;
					break;
				case 3:						// 15-19
					Totals[1] += RangeCnt;
					break;
				case 4:						// 20-29
					Totals[2] = RangeCnt;
					break;
				case 5:						// 30-49
					Totals[2] += RangeCnt;
					break;
				case 6:						// 50-74
					Totals[3] = RangeCnt;
					break;
				case 7:						// 75-99
					Totals[3] += RangeCnt;
					break;
				case 8:						// 100-124
					Totals[4] = RangeCnt;
					break;
				case 9:						// 125-149
					Totals[4] += RangeCnt;
					break;
				case 10:					// 150-174
					Totals[5] = RangeCnt;
					break;
				case 11:					// 175-199
					Totals[5] += RangeCnt;
					break;
				case 12:					// 200-249
					Totals[6] = RangeCnt;
					break;
				case 13:					// 250-300
					Totals[7] = RangeCnt;
					Totals[8] = 0;
					break;
				default:
					Totals[8] += RangeCnt;
					break;
				}
			break;
		case eRCCMinimal:
			switch(NumRows) {
				case 0:						// 0-4
					Totals[0] = RangeCnt;
					break;
				case 1: case 2: case 3:		// 5-9, 10-14,15-19
					Totals[0] += RangeCnt;
					break;
				case 4:						// 20-29
					Totals[1] = RangeCnt;
					break;
				case 5:						// 30-49
					Totals[1] += RangeCnt;
					break;
				case 6:						// 50-74
					Totals[2] = RangeCnt;
					break;
				case 7:						// 75-99
					Totals[2] += RangeCnt;
					break;
				case 8:						// 100-124
					Totals[3] = RangeCnt;
					break;
				case 9:	case 10: case 11:	// 125-149, 150-174,175-199
					Totals[3] += RangeCnt;
					break;
				case 12:					// 200-249
					Totals[4] = RangeCnt;
					break;
				case 13:					// 250-300
					Totals[4] += RangeCnt;
					Totals[5] = 0;
					break;
				default:
					Totals[5] += RangeCnt;
					break;
			}
			break;
		case eRCCMinimalA:
			switch(NumRows) {
				case 0: case 1: case 2: case 3:		// < 20
					Totals[0] += RangeCnt;
					break;

				case 4:	case 5:	case 6:	case 7:		// 20-99
					Totals[1] += RangeCnt;
					break;
										
				case 8:	case 9:	case 10: case 11:	// 100-199
					Totals[2] += RangeCnt;
					break;

				default:
					Totals[3] += RangeCnt;
					break;
			}
			break;
		case eRCCMinimalB:
			switch(NumRows) {
				case 0: case 1: case 2: case 3:		// < 20
					Totals[0] += RangeCnt;
					break;

				case 4:	case 5:			// 20-49
					Totals[1] += RangeCnt;
					break;
										
				case 6:	case 7:			// 50-99
					Totals[2] += RangeCnt;
					break;

				default:
					Totals[3] += RangeCnt;
					break;
			}
			break;

		case eRCCMinimalC:
			switch(NumRows) {
				case 0: case 1: case 2: case 3:		// < 20
					break;

				case 4:	case 5:			// 20-49
					break;
										
				case 6:	case 7:			// 50-99
					break;

				default:
					Totals[0] += RangeCnt;
					break;
			}
			break;
		}
	NumRows++;
	}
if(NumRows != cMaxNumRangeClasses)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expected 28 rows (plus header) in '%s', only '%d' parsed",pszFile,NumRows);
	return(eBSFerrRowCnt);
	}

switch(pParams->RangeClass) {
	case eRCCFull:
		for(Idx = 4; Idx < cMaxNumRangeClasses; Idx++)
			Total += Totals[Idx];
		Len=sprintf(szLineBuff,"\n\"%s\",\"%s\",\"%s\",\"%s\",%d",szElementType,szRefSpecies,szRelSpecies,pszRegion,Total);
		for(Idx = 4; Idx < cMaxNumRangeClasses; Idx++)
			if(pParams->bPercentages)
				Len+=sprintf(&szLineBuff[Len],",%1.2f",Total > 0 ? ((double)Totals[Idx]*100.0)/(double)Total : 0.0);
			else
				Len+=sprintf(&szLineBuff[Len],",%d",Totals[Idx]);		
		break;

	case eRCCReduced:
		for(Idx = 2; Idx < 9; Idx++)
			Total += Totals[Idx];
		Len=sprintf(szLineBuff,"\n\"%s\",\"%s\",\"%s\",\"%s\",%d",szElementType,szRefSpecies,szRelSpecies,pszRegion,Total);
		for(Idx = 2; Idx < 9; Idx++)
			if(pParams->bPercentages)
				Len+=sprintf(&szLineBuff[Len],",%1.2f",Total > 0 ? ((double)Totals[Idx]*100.0)/(double)Total : 0.0);
			else
				Len+=sprintf(&szLineBuff[Len],",%d",Totals[Idx]);		
		break;

	case eRCCMinimal:
		for(Idx = 1; Idx < 6; Idx++)
			Total += Totals[Idx];
		Len=sprintf(szLineBuff,"\n\"%s\",\"%s\",\"%s\",\"%s\",%d",szElementType,szRefSpecies,szRelSpecies,pszRegion,Total);
		for(Idx = 1; Idx < 6; Idx++)
			if(pParams->bPercentages)
				Len+=sprintf(&szLineBuff[Len],",%1.2f",Total > 0 ? ((double)Totals[Idx]*100.0)/(double)Total : 0.0);
			else
				Len+=sprintf(&szLineBuff[Len],",%d",Totals[Idx]);		
		break;

	case eRCCMinimalA:
		for(Idx = 1; Idx < 4; Idx++)
			Total += Totals[Idx];
		Len=sprintf(szLineBuff,"\n\"%s\",\"%s\",\"%s\",\"%s\",%d",szElementType,szRefSpecies,szRelSpecies,pszRegion,Total);
		for(Idx = 1; Idx < 4; Idx++)
			if(pParams->bPercentages)
				Len+=sprintf(&szLineBuff[Len],",%1.2f",Total > 0 ? ((double)Totals[Idx]*100.0)/(double)Total : 0.0);
			else
				Len+=sprintf(&szLineBuff[Len],",%d",Totals[Idx]);		
		break;

	case eRCCMinimalB:
		for(Idx = 1; Idx < 4; Idx++)
			Total += Totals[Idx];
		Len=sprintf(szLineBuff,"\n\"%s\",\"%s\",\"%s\",\"%s\",%d",szElementType,szRefSpecies,szRelSpecies,pszRegion,Total);
		for(Idx = 1; Idx < 4; Idx++)
			if(pParams->bPercentages)
				Len+=sprintf(&szLineBuff[Len],",%1.2f",Total > 0 ? ((double)Totals[Idx]*100.0)/(double)Total : 0.0);
			else
				Len+=sprintf(&szLineBuff[Len],",%d",Totals[Idx]);		
		break;

	case eRCCMinimalC:
		Total += Totals[0];
		Len=sprintf(szLineBuff,"\n\"%s\",\"%s\",\"%s\",\"%s\",%d",szElementType,szRefSpecies,szRelSpecies,pszRegion,Total);
		if(pParams->bPercentages)
			Len+=sprintf(&szLineBuff[Len],",%1.2f",Total > 0 ? ((double)Totals[0]*100.0)/(double)Total : 0.0);
		else
			Len+=sprintf(&szLineBuff[Len],",%d",Totals[0]);		
		break;


	default:
		Len=sprintf(szLineBuff,",-1");
		break;
	}
WriteRsltsHeader(pParams);

if(write(pParams->hRsltsFile,szLineBuff,Len)!=Len)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to write results file '%s' - '%s'",pParams->pszRsltsFile,strerror(errno));
	return(eBSFerrWrite);
	}
return(eBSFSuccess);
}

// Expected CSV format
// ElementID	Uniquely identifies this element (1..n)
// CoreType		Currently either 'hypercore' or 'ultracore'
// Chromosome	Chromosome on reference species
// StartLoci	Starting offset (0..n) on Chromosome
// EndLoci		Ending offset (StartLoci + length -1) on Chromosome
// Length		Element length
// SpeciesList	List of all species sharing this element, 1st species is the reference
// Region		In which region core is located
//
// There are an optional 4 fields present if CSV file was generated with outspecies processing (-x2 option)
// Unaligned	Count of bases in outspecies not aligned
// Matches		Count of bases in outspecies which align and match
// Mismatches	Count of bases in outspecies which align but mismatch
// InDels		Count of bases in reference or outspecies which are InDels
int 
ProcessLociRegionalLengthSummary(char *pszFile,CCSVFile *pCSV,tsProcParams *pParams)
{
int Rslt;
int NumRows;
int NumFields;
int Len;
int Total;
int Totals[cMaxNumRangeClasses];			// to hold all totals
int Idx;
char szLineBuff[2048];
char szRefSpecies[200];
char szRelSpecies[200];
char szElementType[20];
tsLenRangeClass *pRangeClass;
int SeqLen;
int Region;
int NormRegion;
char *pszRegion;

pszRegion = Region2Text(pParams->Region);

if((Rslt=ParseFileNameSpecies(pszFile,NULL,szRefSpecies,szRelSpecies,szElementType))!=eBSFSuccess)
	return(Rslt);

// no headers, straight into the rows which are expected to contain at least 9 fields
NumRows = 0;
Total = 0;
memset(Totals,0,sizeof(Totals));
while((Rslt=pCSV->NextLine()) > 0)	// onto next line containing fields
	{
	NumFields = pCSV->GetCurFields();
	if(NumFields < 9)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expected at least 9 fields in '%s', GetCurFields() returned '%d'",pszFile,NumFields);
		return(eBSFerrFieldCnt);
		}
	pCSV->GetInt(9,&Region);
	NormRegion = NormaliseRegion(Region);
	if(pParams->Region != eFRAuto && pParams->Region != NormRegion)
		continue;

	pCSV->GetInt(7,&SeqLen);
	pRangeClass = GetLengthRangeClass(SeqLen);
	Totals[pRangeClass->ID-1] += SeqLen;
	NumRows++;
	}

switch(pParams->RangeClass) {
	case eRCCFull:
		for(Idx = 4; Idx < cMaxNumRangeClasses; Idx++)
			Total += Totals[Idx];
		Len=sprintf(szLineBuff,"\n\"%s\",\"%s\",\"%s\",\"%s\",%d",szElementType,szRefSpecies,szRelSpecies,pszRegion,Total);
		for(Idx = 4; Idx < cMaxNumRangeClasses; Idx++)
			if(pParams->bPercentages)
				Len+=sprintf(&szLineBuff[Len],",%1.2f",Total > 0 ? ((double)Totals[Idx]*100.0)/(double)Total : 0.0);
			else
				Len+=sprintf(&szLineBuff[Len],",%d",Totals[Idx]);		
		break;

	case eRCCReduced:
		for(Idx = 2; Idx < 9; Idx++)
			Total += Totals[Idx];
		Len=sprintf(szLineBuff,"\n\"%s\",\"%s\",\"%s\",\"%s\",%d",szElementType,szRefSpecies,szRelSpecies,pszRegion,Total);
		for(Idx = 2; Idx < 9; Idx++)
			if(pParams->bPercentages)
				Len+=sprintf(&szLineBuff[Len],",%1.2f",Total > 0 ? ((double)Totals[Idx]*100.0)/(double)Total : 0.0);
			else
				Len+=sprintf(&szLineBuff[Len],",%d",Totals[Idx]);		
		break;

	case eRCCMinimal:
		for(Idx = 1; Idx < 6; Idx++)
			Total += Totals[Idx];
		Len=sprintf(szLineBuff,"\n\"%s\",\"%s\",\"%s\",\"%s\",%d",szElementType,szRefSpecies,szRelSpecies,pszRegion,Total);
		for(Idx = 1; Idx < 6; Idx++)
		if(pParams->bPercentages)
			Len+=sprintf(&szLineBuff[Len],",%1.2f",Total > 0 ? ((double)Totals[Idx]*100.0)/(double)Total : 0.0);
		else
			Len+=sprintf(&szLineBuff[Len],",%d",Totals[Idx]);
		break;

	case eRCCMinimalA:
		for(Idx = 1; Idx < 4; Idx++)
			Total += Totals[Idx];
		Len=sprintf(szLineBuff,"\n\"%s\",\"%s\",\"%s\",\"%s\",%d",szElementType,szRefSpecies,szRelSpecies,pszRegion,Total);
		for(Idx = 1; Idx < 4; Idx++)
		if(pParams->bPercentages)
			Len+=sprintf(&szLineBuff[Len],",%1.2f",Total > 0 ? ((double)Totals[Idx]*100.0)/(double)Total : 0.0);
		else
			Len+=sprintf(&szLineBuff[Len],",%d",Totals[Idx]);
		break;

	case eRCCMinimalB:
		for(Idx = 1; Idx < 4; Idx++)
			Total += Totals[Idx];
		Len=sprintf(szLineBuff,"\n\"%s\",\"%s\",\"%s\",\"%s\",%d",szElementType,szRefSpecies,szRelSpecies,pszRegion,Total);
		for(Idx = 1; Idx < 4; Idx++)
		if(pParams->bPercentages)
			Len+=sprintf(&szLineBuff[Len],",%1.2f",Total > 0 ? ((double)Totals[Idx]*100.0)/(double)Total : 0.0);
		else
			Len+=sprintf(&szLineBuff[Len],",%d",Totals[Idx]);
		break;

	case eRCCMinimalC:
		Total += Totals[0];
		Len=sprintf(szLineBuff,"\n\"%s\",\"%s\",\"%s\",\"%s\",%d",szElementType,szRefSpecies,szRelSpecies,pszRegion,Total);
		if(pParams->bPercentages)
			Len+=sprintf(&szLineBuff[Len],",%1.2f",Total > 0 ? ((double)Totals[0]*100.0)/(double)Total : 0.0);
		else
			Len+=sprintf(&szLineBuff[Len],",%d",Totals[0]);
		break;

	default:
		Len=sprintf(szLineBuff,",-1");
		break;
	}
WriteRsltsHeader(pParams);

if(write(pParams->hRsltsFile,szLineBuff,Len)!=Len)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to write results file '%s' - '%s'",pParams->pszRsltsFile,strerror(errno));
	return(eBSFerrWrite);
	}

return(eBSFSuccess);
}


int
FormatSummaryLine(char *pszBuffer,const char *pszElementType,const char *pszRefSpecies,const char *pszCoreList,const char *pszOutSpecies,const char *pszRegion,const char *pszQual,const char *pszName,int Value)
{
int Len = sprintf(pszBuffer,"\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%d\n",
				  pszElementType,pszRefSpecies,pszCoreList,pszOutSpecies,pszRegion,pszQual,pszName,Value);
return(Len);
}


typedef struct TAG_sDistSegCnts {
	double Matches;			// number of exact matches in this segment
	double Mismatches;			// number of mismatches
	double InDels;				// total number of InDel bases
	double Unaligned;			// number of unaligned bases
	} tsDistSegCnts;

typedef struct TAG_sOutSpeciesSummary {
	int OGMatchBaseTotal;
	int OGMismatchBaseTotal;
	int OGInDelBaseTotal;
	int OGUnalignedBaseTotal;
	int OGUnalignedCoreBaseTotal;

	int InstTotal;								// total number of core or reference instances
	int AlignedInstTotal;						// total number of aligned cores
	int BaseTotal;								// total number of bases over all core or reference instances
	
	int RegionTotals[7];						// genomic region totals
	int IntronExonSpliceSiteTotals;				// 5' splice site totals
	int ExonIntronSpliceSiteTotals;				// 3' splice site totals


	int InstTotals[cMaxNumRangeClasses];		// to hold total core or reference instance counts
	int BaseTotals[cMaxNumRangeClasses];		// total number of bases over all core or reference instances in each length range
	int	OGUnalignedCores[cMaxNumRangeClasses];	// to hold outspecies number of instances with no alignment to core
	int	OGUnalignedCoreBases[cMaxNumRangeClasses]; // to hold total number of bases in instances of cores unaligned to out species

	int OGMatchBaseTotals[cMaxNumRangeClasses];	// total number of bases in out species matching core 
	int OGMismatchBaseTotals[cMaxNumRangeClasses]; // total number of bases in out species not matching core
	int OGInDelBaseTotals[cMaxNumRangeClasses];	 // total number of bases in out species which are InDels
	int OGUnalignedBaseTotals[cMaxNumRangeClasses]; // total number of bases in out species which are not aligned
	int OGAligned2Cores[cMaxNumRangeClasses];		// total number of alignments of at least MinAlignedBases in outspecies

	int OGRegions[cMaxNumRangeClasses][7];			// region counts
	int OGIntronExonSpliceSite[cMaxNumRangeClasses];	// 5'splice site
	int OGExonIntronSpliceSite[cMaxNumRangeClasses];	// 3'splice site

	int OGIntergenicRegion[cMaxNumRangeClasses];	// contained exclusively within intergenic region counts
	int OGIntergenicRegionTotals;
	int OGUSRegion[cMaxNumRangeClasses];			// contained exclusively within US counts
	int	OGUSRegionTotals;
	int OG5UTRRegion[cMaxNumRangeClasses];			// contained exclusively within 5'UTR counts
	int	OG5UTRRegionTotals;
	int OGCDSRegion[cMaxNumRangeClasses];			// contained exclusively within CDS counts
	int OGCDSRegionTotals;
	int OGIntronRegion[cMaxNumRangeClasses];		// contained exclusively within Intron counts
	int OGIntronRegionTotals;
	int OG3UTRRegion[cMaxNumRangeClasses];			// contained exclusively within 3'UTR counts
	int OG3UTRRegionTotals;
	int OGDSRegion[cMaxNumRangeClasses];			// contained exclusively within DS counts
	int	OGDSRegionTotals;
	int OGSpliceSites[cMaxNumRangeClasses];			// counts of either 5' or 3' splice sites
	int OGSpliceSitesTotals;

	int OGMatchIdentities[cMaxNumRangeClasses][cMaxNumIdentityClasses];	// number of ref core elements with identities 0..100 in each range class
	int OGTotMatchIdentities[cMaxNumIdentityClasses]; // total number of ref core elements with identities 0..100

	tsDistSegCnts DistSegCnts[cMaxNumRangeClasses][cMaxMatchDistSegments]; // to hold distribution segment counts for each length range
	tsDistSegCnts DistSegCntsTotals[cMaxNumRangeClasses];

	int OGUnalignedCoresTotal;

	char szRefSpecies[200];
	char szOutSpecies[200];
	char szCoreList[1024];
	char szElementType[120];

	tsLenRangeClass *pRangeClass;
	char *pszRegion;
	int NumDistSegs;
	int StartIdx;
	int EndIdx;
} tsOutSpeciesSummary;


int
RsltStandard(tsOutSpeciesSummary *pOSSummary,tsProcParams *pParams)
{
char szLineBuff[0x03fff];
char szSegBuff[128];
int SegIdx;
int Idx;

int Len;

WriteRsltsHeader(pParams);

Len = FormatSummaryLine(szLineBuff,pOSSummary->szElementType,pOSSummary->szRefSpecies,pOSSummary->szCoreList,pOSSummary->szOutSpecies,pOSSummary->pszRegion,"All Cores","Total Cores",pOSSummary->InstTotal);
Len += FormatSummaryLine(&szLineBuff[Len],pOSSummary->szElementType,pOSSummary->szRefSpecies,pOSSummary->szCoreList,pOSSummary->szOutSpecies,pOSSummary->pszRegion,"All Cores","Total Bases in Cores",pOSSummary->BaseTotal);
Len += FormatSummaryLine(&szLineBuff[Len],pOSSummary->szElementType,pOSSummary->szRefSpecies,pOSSummary->szCoreList,pOSSummary->szOutSpecies,pOSSummary->pszRegion,"All Cores","Total Cores to OG Unaligned",pOSSummary->OGUnalignedCoresTotal);
Len += FormatSummaryLine(&szLineBuff[Len],pOSSummary->szElementType,pOSSummary->szRefSpecies,pOSSummary->szCoreList,pOSSummary->szOutSpecies,pOSSummary->pszRegion,"All Cores","Total Bases in Unaligned Cores",pOSSummary->OGUnalignedCoreBaseTotal);

Len += FormatSummaryLine(&szLineBuff[Len],pOSSummary->szElementType,pOSSummary->szRefSpecies,pOSSummary->szCoreList,pOSSummary->szOutSpecies,pOSSummary->pszRegion,"All Cores","OG Matching Bases",pOSSummary->OGMatchBaseTotal);
Len += FormatSummaryLine(&szLineBuff[Len],pOSSummary->szElementType,pOSSummary->szRefSpecies,pOSSummary->szCoreList,pOSSummary->szOutSpecies,pOSSummary->pszRegion,"All Cores","OG Missmatch Bases",pOSSummary->OGMismatchBaseTotal);
Len += FormatSummaryLine(&szLineBuff[Len],pOSSummary->szElementType,pOSSummary->szRefSpecies,pOSSummary->szCoreList,pOSSummary->szOutSpecies,pOSSummary->pszRegion,"All Cores","OG InDel Bases",pOSSummary->OGInDelBaseTotal);
Len += FormatSummaryLine(&szLineBuff[Len],pOSSummary->szElementType,pOSSummary->szRefSpecies,pOSSummary->szCoreList,pOSSummary->szOutSpecies,pOSSummary->pszRegion,"All Cores","OG Unaligned Bases",pOSSummary->OGUnalignedBaseTotal);

for(SegIdx = 0; SegIdx < pOSSummary->NumDistSegs; SegIdx++)
	{
	sprintf(szSegBuff,"Segment %d",SegIdx+1);
	Len += FormatSummaryLine(&szLineBuff[Len],pOSSummary->szElementType,pOSSummary->szRefSpecies,
				pOSSummary->szCoreList,pOSSummary->szOutSpecies,
				pOSSummary->pszRegion,"All Cores (Matches)",szSegBuff,
				(int)(pOSSummary->DistSegCntsTotals[SegIdx].Matches + 0.5));
	Len += FormatSummaryLine(&szLineBuff[Len],pOSSummary->szElementType,pOSSummary->szRefSpecies,
				pOSSummary->szCoreList,pOSSummary->szOutSpecies,
				pOSSummary->pszRegion,"All Cores (Mismatches)",szSegBuff,
				(int)(pOSSummary->DistSegCntsTotals[SegIdx].Mismatches + 0.5));
	Len += FormatSummaryLine(&szLineBuff[Len],pOSSummary->szElementType,pOSSummary->szRefSpecies,
				pOSSummary->szCoreList,pOSSummary->szOutSpecies,
				pOSSummary->pszRegion,"All Cores (InDels)",szSegBuff,
				(int)(pOSSummary->DistSegCntsTotals[SegIdx].InDels + 0.5));
	Len += FormatSummaryLine(&szLineBuff[Len],pOSSummary->szElementType,pOSSummary->szRefSpecies,
				pOSSummary->szCoreList,pOSSummary->szOutSpecies,
				pOSSummary->pszRegion,"All Cores (Unaligned)",szSegBuff,
				(int)(pOSSummary->DistSegCntsTotals[SegIdx].Unaligned + 0.5));

	}

if(write(pParams->hRsltsFile,szLineBuff,Len)!=Len)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to write results file '%s' - '%s'",pParams->pszRsltsFile,strerror(errno));
	return(eBSFerrWrite);
	}
Len = 0;
for(Idx = pOSSummary->StartIdx; Idx < pOSSummary->EndIdx; Idx++)
	{
	char *pszRange = (char *)pOSSummary->pRangeClass[Idx].pszDescr;
	Len += FormatSummaryLine(&szLineBuff[Len],pOSSummary->szElementType,pOSSummary->szRefSpecies,pOSSummary->szCoreList,pOSSummary->szOutSpecies,pOSSummary->pszRegion,pszRange,"Cores",pOSSummary->InstTotals[Idx]);
	Len += FormatSummaryLine(&szLineBuff[Len],pOSSummary->szElementType,pOSSummary->szRefSpecies,pOSSummary->szCoreList,pOSSummary->szOutSpecies,pOSSummary->pszRegion,pszRange,"Bases in Cores",pOSSummary->BaseTotals[Idx]);
	Len += FormatSummaryLine(&szLineBuff[Len],pOSSummary->szElementType,pOSSummary->szRefSpecies,pOSSummary->szCoreList,pOSSummary->szOutSpecies,pOSSummary->pszRegion,pszRange,"Cores to OG Unaligned",pOSSummary->OGUnalignedCores[Idx]);
	Len += FormatSummaryLine(&szLineBuff[Len],pOSSummary->szElementType,pOSSummary->szRefSpecies,pOSSummary->szCoreList,pOSSummary->szOutSpecies,pOSSummary->pszRegion,pszRange,"Bases in Unaligned Cores",pOSSummary->OGUnalignedCoreBases[Idx]);

	Len += FormatSummaryLine(&szLineBuff[Len],pOSSummary->szElementType,pOSSummary->szRefSpecies,pOSSummary->szCoreList,pOSSummary->szOutSpecies,pOSSummary->pszRegion,pszRange,"OG Matching Bases",pOSSummary->OGMatchBaseTotals[Idx]);
	Len += FormatSummaryLine(&szLineBuff[Len],pOSSummary->szElementType,pOSSummary->szRefSpecies,pOSSummary->szCoreList,pOSSummary->szOutSpecies,pOSSummary->pszRegion,pszRange,"OG Missmatch Bases",pOSSummary->OGMismatchBaseTotals[Idx]);
	Len += FormatSummaryLine(&szLineBuff[Len],pOSSummary->szElementType,pOSSummary->szRefSpecies,pOSSummary->szCoreList,pOSSummary->szOutSpecies,pOSSummary->pszRegion,pszRange,"OG InDel Bases",pOSSummary->OGInDelBaseTotals[Idx]);
	Len += FormatSummaryLine(&szLineBuff[Len],pOSSummary->szElementType,pOSSummary->szRefSpecies,pOSSummary->szCoreList,pOSSummary->szOutSpecies,pOSSummary->pszRegion,pszRange,"OG Unaligned Bases",pOSSummary->OGUnalignedBaseTotals[Idx]);

	for(SegIdx = 0; SegIdx < pOSSummary->NumDistSegs; SegIdx++)
		{
		sprintf(szSegBuff,"Segment %d (Matches)",SegIdx+1);
		Len += FormatSummaryLine(&szLineBuff[Len],pOSSummary->szElementType,pOSSummary->szRefSpecies,
			pOSSummary->szCoreList,pOSSummary->szOutSpecies,pOSSummary->pszRegion,pszRange,szSegBuff,
			(int)(pOSSummary->DistSegCnts[Idx][SegIdx].Matches+0.5));

		sprintf(szSegBuff,"Segment %d (Mismatches)",SegIdx+1);
		Len += FormatSummaryLine(&szLineBuff[Len],pOSSummary->szElementType,pOSSummary->szRefSpecies,
			pOSSummary->szCoreList,pOSSummary->szOutSpecies,pOSSummary->pszRegion,pszRange,szSegBuff,
			(int)(pOSSummary->DistSegCnts[Idx][SegIdx].Mismatches+0.5));

		sprintf(szSegBuff,"Segment %d (InDels)",SegIdx+1);
		Len += FormatSummaryLine(&szLineBuff[Len],pOSSummary->szElementType,pOSSummary->szRefSpecies,
			pOSSummary->szCoreList,pOSSummary->szOutSpecies,pOSSummary->pszRegion,pszRange,szSegBuff,
			(int)(pOSSummary->DistSegCnts[Idx][SegIdx].InDels+0.5));

		sprintf(szSegBuff,"Segment %d (Unaligned)",SegIdx+1);
		Len += FormatSummaryLine(&szLineBuff[Len],pOSSummary->szElementType,pOSSummary->szRefSpecies,
			pOSSummary->szCoreList,pOSSummary->szOutSpecies,pOSSummary->pszRegion,pszRange,szSegBuff,
			(int)(pOSSummary->DistSegCnts[Idx][SegIdx].Unaligned+0.5));
		}

	if(write(pParams->hRsltsFile,szLineBuff,Len)!=Len)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to write results file '%s' - '%s'",pParams->pszRsltsFile,strerror(errno));
		return(eBSFerrWrite);
		}
	Len = 0;
	}

if(Len)
	if(write(pParams->hRsltsFile,szLineBuff,Len)!=Len)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to write results file '%s' - '%s'",pParams->pszRsltsFile,strerror(errno));
		return(eBSFerrWrite);
		}

return(eBSFSuccess);
}

// RsltTable
// Generates table showing matching base identities for each outspecies
// <refspecies>,<specieslist>,<outspecies>,<region>,<refcorelen>,<numrefcores>,<numoutcores>,<identity1>,<identity2>
// <identity1> Identity as NumMatchingBases/NumBasesInCoresWhichAlignToRefCores
// <identity2> Identity as NumMatchingBases/(NumMatchingBases + NumMissmatchingBases)

int
RsltTable(tsOutSpeciesSummary *pOSSummary,tsProcParams *pParams)
{
char szLineBuff[0x03fff];
int Len;

WriteRsltsHeader(pParams);

char szCoreSet[200];
strcpy(szCoreSet,pOSSummary->szRefSpecies);
strcat(szCoreSet,pOSSummary->szCoreList);
Len = sprintf(szLineBuff,"\"%s\",\"%s\",\"%s\",%d,%d,%2.2f,%d,%d,%d,%2.2f\n",
			szCoreSet,pOSSummary->szOutSpecies,"All",		// core length range
			pOSSummary->InstTotal,
			pOSSummary->AlignedInstTotal,
			pOSSummary->InstTotal == 0 ? 0.0 : (pOSSummary->AlignedInstTotal*100.0)/pOSSummary->InstTotal,
			pOSSummary->OGMatchBaseTotal + pOSSummary->OGMismatchBaseTotal,
			pOSSummary->OGMismatchBaseTotal,
			pOSSummary->OGMatchBaseTotal,
			(pOSSummary->OGMatchBaseTotal + pOSSummary->OGMismatchBaseTotal) == 0 ? 0.0 :
			(pOSSummary->OGMatchBaseTotal * 100.0) / (pOSSummary->OGMatchBaseTotal + pOSSummary->OGMismatchBaseTotal));


for(int Idx = pOSSummary->StartIdx; Idx < pOSSummary->EndIdx; Idx++)
	{
	const char *pszRange = pOSSummary->pRangeClass[Idx].pszDescr;
	Len += sprintf(&szLineBuff[Len],"\"%s\",\"%s\",\"%s\",%d,%d,%2.2f,%d,%d,%d,%2.2f\n",
			szCoreSet,pOSSummary->szOutSpecies,pszRange,								// core length range
			pOSSummary->InstTotals[Idx],
			pOSSummary->OGAligned2Cores[Idx],
			pOSSummary->InstTotals[Idx] == 0 ? 0.0 :
			((pOSSummary->InstTotals[Idx] - pOSSummary->OGUnalignedCores[Idx]) * 100.0) /pOSSummary->InstTotals[Idx],
			pOSSummary->OGMatchBaseTotals[Idx] + pOSSummary->OGMismatchBaseTotals[Idx],
			pOSSummary->OGMismatchBaseTotals[Idx],
			pOSSummary->OGMatchBaseTotals[Idx],
			(pOSSummary->OGMatchBaseTotals[Idx] + pOSSummary->OGMismatchBaseTotals[Idx]) == 0 ? 0.0 :
			(pOSSummary->OGMatchBaseTotals[Idx] * 100.0) / (pOSSummary->OGMatchBaseTotals[Idx] + pOSSummary->OGMismatchBaseTotals[Idx]));

	}

if(Len)
	if(write(pParams->hRsltsFile,szLineBuff,Len)!=Len)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to write results file '%s' - '%s'",pParams->pszRsltsFile,strerror(errno));
		return(eBSFerrWrite);
		}

return(eBSFSuccess);
}

// RsltSegDistRows
// Output the segment distributions as proportions in rows
int
RsltSegDistRows(tsOutSpeciesSummary *pOSSummary,tsProcParams *pParams)
{
char szLineBuff[0x03fff];
int SegIdx;
int Idx;
double SegTotal;

int Len = 0;
WriteRsltsHeader(pParams);
for(Idx = pOSSummary->StartIdx; Idx < pOSSummary->EndIdx; Idx++)
	{
	const char *pszRange = pOSSummary->pRangeClass[Idx].pszDescr;
	Len = sprintf(szLineBuff,"\"%s\",\"%s\",\"%s\",\"%s\",\"%s\"",
		pOSSummary->szElementType,pOSSummary->szRefSpecies,pOSSummary->szCoreList,pOSSummary->szOutSpecies,pszRange);
	for(SegIdx = 0; SegIdx < pOSSummary->NumDistSegs; SegIdx++)
		{
		if(pParams->bPercentages)
			{
			SegTotal = pOSSummary->DistSegCnts[Idx][SegIdx].Matches +
					   pOSSummary->DistSegCnts[Idx][SegIdx].Mismatches +
						pOSSummary->DistSegCnts[Idx][SegIdx].InDels +
						pOSSummary->DistSegCnts[Idx][SegIdx].Unaligned;
			if(SegTotal <= 0.0)
				Len += sprintf(&szLineBuff[Len],",0.00000,0.00000,0.00000,0.00000");
			else
				Len += sprintf(&szLineBuff[Len],",%1.8f,%1.8f,%1.8f,%1.8f",
							pOSSummary->DistSegCnts[Idx][SegIdx].Matches/SegTotal,
							pOSSummary->DistSegCnts[Idx][SegIdx].Mismatches/SegTotal,
							pOSSummary->DistSegCnts[Idx][SegIdx].InDels/SegTotal,
							pOSSummary->DistSegCnts[Idx][SegIdx].Unaligned/SegTotal);
			}
		else
			Len += sprintf(&szLineBuff[Len],",%d,%d,%d,%d",
				(int)(pOSSummary->DistSegCnts[Idx][SegIdx].Matches+0.5),
				(int)(pOSSummary->DistSegCnts[Idx][SegIdx].Mismatches+0.5),
				(int)(pOSSummary->DistSegCnts[Idx][SegIdx].InDels+0.5),
				(int)(pOSSummary->DistSegCnts[Idx][SegIdx].Unaligned+0.5));
		}
	Len += sprintf(&szLineBuff[Len],"\n");
	if(write(pParams->hRsltsFile,szLineBuff,Len)!=Len)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to write results file '%s' - '%s'",pParams->pszRsltsFile,strerror(errno));
		return(eBSFerrWrite);
		}
	Len = 0;
	}

if(Len)
	if(write(pParams->hRsltsFile,szLineBuff,Len)!=Len)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to write results file '%s' - '%s'",pParams->pszRsltsFile,strerror(errno));
		return(eBSFerrWrite);
		}

return(eBSFSuccess);
}

// RsltSegDistTables
// Output the segment distributions as a table with segments the columns and rows from matches/mismatches/indels/unaligned
int
RsltSegDistTables(tsOutSpeciesSummary *pOSSummary,tsProcParams *pParams)
{
char szLineBuff[0x03fff];
int SegIdx;
int Idx;
double SegTotals[500];
double CntCatVal;

int Len = 0;
WriteRsltsHeader(pParams);
for(Idx = pOSSummary->StartIdx; Idx < pOSSummary->EndIdx; Idx++)
	{
	const char *pszRange = pOSSummary->pRangeClass[Idx].pszDescr;

	if(pParams->bPercentages)
		{
		for(SegIdx = 0; SegIdx < pOSSummary->NumDistSegs; SegIdx++)
			{
			SegTotals[SegIdx] = pOSSummary->DistSegCnts[Idx][SegIdx].Matches +
						   pOSSummary->DistSegCnts[Idx][SegIdx].Mismatches +
							pOSSummary->DistSegCnts[Idx][SegIdx].InDels +
							pOSSummary->DistSegCnts[Idx][SegIdx].Unaligned;
			if(SegTotals[SegIdx] < 0.0)
				SegTotals[SegIdx] = 0.0;
			}
		}

	for(int CntCat=0;CntCat<4;CntCat++)
		{
		Len = sprintf(szLineBuff,"\"%s\",\"%s\",\"%s\",\"%s\",\"%s\"",
			pOSSummary->szElementType,pOSSummary->szRefSpecies,pOSSummary->szCoreList,pOSSummary->szOutSpecies,pszRange);

		switch(CntCat) {
			case 0:			// matches
				Len += sprintf(&szLineBuff[Len],",Matches");
				break;
			case 1:			// mismatches
				Len += sprintf(&szLineBuff[Len],",Mismatches");
				break;

			case 2:			// indels
				Len += sprintf(&szLineBuff[Len],",InDels");
				break;

			case 3:			// unaligned
				Len += sprintf(&szLineBuff[Len],",Unaligned");
				break;
			}

		for(SegIdx = 0; SegIdx < pOSSummary->NumDistSegs; SegIdx++)
			{
			switch(CntCat) {
				case 0:			// matches
					CntCatVal=pOSSummary->DistSegCnts[Idx][SegIdx].Matches;
					break;
				case 1:			// mismatches
					CntCatVal=pOSSummary->DistSegCnts[Idx][SegIdx].Mismatches;
					break;

				case 2:			// indels
					CntCatVal=pOSSummary->DistSegCnts[Idx][SegIdx].InDels;
					break;

				case 3:			// unaligned
					CntCatVal=pOSSummary->DistSegCnts[Idx][SegIdx].Unaligned;
					break;
				}

			if(pParams->bPercentages)
				{
				if(SegTotals[SegIdx] <= 0.0)
					Len += sprintf(&szLineBuff[Len],",0.00000");
				else
					Len += sprintf(&szLineBuff[Len],",%1.8f",CntCatVal/SegTotals[SegIdx]);
				}
			else
				Len += sprintf(&szLineBuff[Len],",%d",(int)(CntCatVal+0.5));
			}
		Len += sprintf(&szLineBuff[Len],"\n");
		if(write(pParams->hRsltsFile,szLineBuff,Len)!=Len)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to write results file '%s' - '%s'",pParams->pszRsltsFile,strerror(errno));
			return(eBSFerrWrite);
			}
		Len = 0;
		}
	}

if(Len)
	if(write(pParams->hRsltsFile,szLineBuff,Len)!=Len)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to write results file '%s' - '%s'",pParams->pszRsltsFile,strerror(errno));
		return(eBSFerrWrite);
		}

return(eBSFSuccess);
}

// RsltSegDistCols
// Output the segment distributions as proportions in columns
int
RsltSegDistCols(tsOutSpeciesSummary *pOSSummary,tsProcParams *pParams)
{
char szLineBuff[0x03fff];
int SegIdx;
int Idx;
double SegTotal;
char *pszRange;

int Len = 0;
WriteRsltsHeader(pParams);
for(Idx = pOSSummary->StartIdx; Idx < pOSSummary->EndIdx; Idx++)
	{
	pszRange = (char *)pOSSummary->pRangeClass[Idx].pszDescr;
	Len = 0;
	for(SegIdx = 0; SegIdx < pOSSummary->NumDistSegs; SegIdx++)
		{
		Len += sprintf(&szLineBuff[Len],"\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",\"Segment%d\"",
			pOSSummary->szElementType,pOSSummary->szRefSpecies,pOSSummary->szCoreList,pOSSummary->szOutSpecies,pszRange,SegIdx+1);

		if(pParams->bPercentages)
			{
			SegTotal = pOSSummary->DistSegCnts[Idx][SegIdx].Matches +
					   pOSSummary->DistSegCnts[Idx][SegIdx].Mismatches +
						pOSSummary->DistSegCnts[Idx][SegIdx].InDels +
						pOSSummary->DistSegCnts[Idx][SegIdx].Unaligned;
			if(SegTotal <= 0.0)
				Len += sprintf(&szLineBuff[Len],",0.00000,0.00000,0.00000,0.00000\n");
			else
				Len += sprintf(&szLineBuff[Len],",%1.8f,%1.8f,%1.8f,%1.8f\n",
							pOSSummary->DistSegCnts[Idx][SegIdx].Matches/SegTotal,
							pOSSummary->DistSegCnts[Idx][SegIdx].Mismatches/SegTotal,
							pOSSummary->DistSegCnts[Idx][SegIdx].InDels/SegTotal,
							pOSSummary->DistSegCnts[Idx][SegIdx].Unaligned/SegTotal);
			}
		else
			Len += sprintf(&szLineBuff[Len],",%d,%d,%d,%d\n",
				(int)(pOSSummary->DistSegCnts[Idx][SegIdx].Matches+0.5),
				(int)(pOSSummary->DistSegCnts[Idx][SegIdx].Mismatches+0.5),
				(int)(pOSSummary->DistSegCnts[Idx][SegIdx].InDels+0.5),
				(int)(pOSSummary->DistSegCnts[Idx][SegIdx].Unaligned+0.5));
		}
	if(write(pParams->hRsltsFile,szLineBuff,Len)!=Len)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to write results file '%s' - '%s'",pParams->pszRsltsFile,strerror(errno));
		return(eBSFerrWrite);
		}
	Len = 0;
	}

if(Len)
	if(write(pParams->hRsltsFile,szLineBuff,Len)!=Len)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to write results file '%s' - '%s'",pParams->pszRsltsFile,strerror(errno));
		return(eBSFerrWrite);
		}

return(eBSFSuccess);
}

// RsltIdentDistCols
// Output the identities and number of instances for each length range class in columns
int
RsltIdentDistCols(tsOutSpeciesSummary *pOSSummary,tsProcParams *pParams)
{
char szLineBuff[0x07fff];
int Idx;
char *pszRange;

int Len = 0;
WriteRsltsHeader(pParams);
for(Idx = pOSSummary->StartIdx; Idx < pOSSummary->EndIdx; Idx++)
	{
	pszRange = (char *)pOSSummary->pRangeClass[Idx].pszDescr;
	Len = 0;
	for(int IdentityIdx = 0; IdentityIdx < cMaxNumIdentityClasses; IdentityIdx++)
		{
		Len += sprintf(&szLineBuff[Len],"\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%d,%d,%2.1f,%d\n",
			pOSSummary->szElementType,pOSSummary->szRefSpecies,pOSSummary->szCoreList,pOSSummary->szOutSpecies,pszRange,
			pOSSummary->InstTotals[Idx],pOSSummary->InstTotals[Idx] - pOSSummary->OGUnalignedCores[Idx],
			(double)IdentityIdx, 
			pOSSummary->OGMatchIdentities[Idx][IdentityIdx]);
		}

	if(write(pParams->hRsltsFile,szLineBuff,Len)!=Len)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to write results file '%s' - '%s'",pParams->pszRsltsFile,strerror(errno));
		return(eBSFerrWrite);
		}
	Len = 0;
	}

if(Len)
	if(write(pParams->hRsltsFile,szLineBuff,Len)!=Len)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to write results file '%s' - '%s'",pParams->pszRsltsFile,strerror(errno));
		return(eBSFerrWrite);
		}

return(eBSFSuccess);
}

// RsltRegionTable
// Output the regional characterisation by length range
int
RsltRegionTable(tsOutSpeciesSummary *pOSSummary,tsProcParams *pParams)
{
char szLineBuff[0x03fff];
int Len;

WriteRsltsHeader(pParams);
char szCoreSet[200];
strcpy(szCoreSet,pOSSummary->szRefSpecies);
strcat(szCoreSet,pOSSummary->szCoreList);
Len = sprintf(szLineBuff,"\"%s\",\"%s\",\"%s\",%d,%d",
			szCoreSet,pOSSummary->szOutSpecies,"All",		// core length range
			pOSSummary->InstTotal,
			pOSSummary->AlignedInstTotal);

for(int RegionIdx = 0; RegionIdx < 7; RegionIdx++)
	Len += sprintf(&szLineBuff[Len],",%d",pOSSummary->RegionTotals[RegionIdx]);
Len += sprintf(&szLineBuff[Len],",%d,%d\n",pOSSummary->IntronExonSpliceSiteTotals,pOSSummary->ExonIntronSpliceSiteTotals);

for(int Idx = pOSSummary->StartIdx; Idx < pOSSummary->EndIdx; Idx++)
	{
	char *pszRange = (char *)pOSSummary->pRangeClass[Idx].pszDescr;
	Len += sprintf(&szLineBuff[Len],"\"%s\",\"%s\",\"%s\",%d,%d",
			szCoreSet,pOSSummary->szOutSpecies,pszRange,								// core length range
			pOSSummary->InstTotals[Idx],
			pOSSummary->OGAligned2Cores[Idx]);
	for(int RegionIdx = 0; RegionIdx < 7; RegionIdx++)
		Len += sprintf(&szLineBuff[Len],",%d",pOSSummary->OGRegions[Idx][RegionIdx]);
	Len += sprintf(&szLineBuff[Len],",%d,%d\n",pOSSummary->OGIntronExonSpliceSite[Idx],pOSSummary->OGExonIntronSpliceSite[Idx]);
	}

if(Len)
	if(write(pParams->hRsltsFile,szLineBuff,Len)!=Len)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to write results file '%s' - '%s'",pParams->pszRsltsFile,strerror(errno));
		return(eBSFerrWrite);
		}

return(eBSFSuccess);
}


// RsltRegionTableA
// Output the regional characterisation by length range
int
RsltRegionTableA(tsOutSpeciesSummary *pOSSummary,tsProcParams *pParams)
{
char szLineBuff[0x03fff];
int Len;

WriteRsltsHeader(pParams);
char szCoreSet[2048];
strcpy(szCoreSet,pOSSummary->szRefSpecies);
strcat(szCoreSet,pOSSummary->szCoreList);
Len = sprintf(szLineBuff,"\"%s\",\"%s\",\"%s\",%d,%d",
			szCoreSet,pOSSummary->szOutSpecies,"All",		// core length range
			pOSSummary->InstTotal,
			pOSSummary->AlignedInstTotal);

for(int RegionIdx = 0; RegionIdx < 7; RegionIdx++)
	Len += sprintf(&szLineBuff[Len],",%d",pOSSummary->RegionTotals[RegionIdx]);
Len += sprintf(&szLineBuff[Len],",%d,%d",
			   pOSSummary->IntronExonSpliceSiteTotals,pOSSummary->ExonIntronSpliceSiteTotals);
Len += sprintf(&szLineBuff[Len],",%d,%d,%d,%d,%d,%d,%d,%d\n",
			   pOSSummary->OGIntergenicRegionTotals,pOSSummary->OGUSRegionTotals,
			   pOSSummary->OG5UTRRegionTotals,pOSSummary->OGCDSRegionTotals,
			   pOSSummary->OGIntronRegionTotals,pOSSummary->OG3UTRRegionTotals,
			   pOSSummary->OGDSRegionTotals,pOSSummary->OGSpliceSitesTotals);

for(int Idx = pOSSummary->StartIdx; Idx < pOSSummary->EndIdx; Idx++)
	{
	char *pszRange = (char *)pOSSummary->pRangeClass[Idx].pszDescr;
	Len += sprintf(&szLineBuff[Len],"\"%s\",\"%s\",\"%s\",%d,%d",
			szCoreSet,pOSSummary->szOutSpecies,pszRange,								// core length range
			pOSSummary->InstTotals[Idx],
			pOSSummary->OGAligned2Cores[Idx]);
	for(int RegionIdx = 0; RegionIdx < 7; RegionIdx++)
		Len += sprintf(&szLineBuff[Len],",%d",pOSSummary->OGRegions[Idx][RegionIdx]);
	Len += sprintf(&szLineBuff[Len],",%d,%d",pOSSummary->OGIntronExonSpliceSite[Idx],pOSSummary->OGExonIntronSpliceSite[Idx]);
	Len += sprintf(&szLineBuff[Len],",%d,%d,%d,%d,%d,%d,%d,%d\n",
			   pOSSummary->OGIntergenicRegion[Idx],pOSSummary->OGUSRegion[Idx],
			   pOSSummary->OG5UTRRegion[Idx],pOSSummary->OGCDSRegion[Idx],
			   pOSSummary->OGIntronRegion[Idx],pOSSummary->OG3UTRRegion[Idx],
			   pOSSummary->OGDSRegion[Idx],pOSSummary->OGSpliceSites[Idx]);
	}

if(Len)
	if(write(pParams->hRsltsFile,szLineBuff,Len)!=Len)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to write results file '%s' - '%s'",pParams->pszRsltsFile,strerror(errno));
		return(eBSFerrWrite);
		}

return(eBSFSuccess);
}


// Expected CSV format
// ElementID	Uniquely identifies this element (1..n)
// CoreType		Currently either 'hypercore' or 'ultracore'
// Chromosome	Chromosome on reference species
// StartLoci	Starting offset (0..n) on Chromosome
// EndLoci		Ending offset (StartLoci + length -1) on Chromosome
// Length		Element length
// SpeciesList	List of all species sharing this element, 1st species is the reference
// Region		In which region core is located
//
// There are an optional 4 fields present if CSV file was generated with outspecies processing (-x2 option)
// Unaligned	Count of bases in outspecies not aligned
// Matches		Count of bases in outspecies which align and match
// Mismatches	Count of bases in outspecies which align but mismatch
// InDels		Count of bases in reference or outspecies which are InDels
int 
ProcessOutspeciesSummary(char *pszFile,CCSVFile *pCSV,tsProcParams *pParams)
{
int Rslt;
int NumRows;
int NumFields;

int AlignedBases;
int OGUnalignedBases;
int OGMatchCnt;
int OGMismatchCnt;
int OGInDelCnt;
int SegCnt;
int SegIdx;

char *pszElementType;
char *pszSpeciesList;
char *pszRefSpecies;
int SeqLen;
int Region;
int NormRegion;

tsLenRangeClass *pRangeClass;


tsOutSpeciesSummary OSSummary;


// no headers, straight into the rows which are expected to contain at least 13 fields
NumRows = 0;
memset(&OSSummary,0,sizeof(OSSummary));
OSSummary.pszRegion = Region2Text(pParams->Region);


switch(pParams->RangeClass) {
	case eRCCFull:
		OSSummary.StartIdx = 4;
		OSSummary.EndIdx = cMaxNumRangeClasses;
		break;

	case eRCCReduced:
		OSSummary.StartIdx = 2; 
		OSSummary.EndIdx = 9;
		break;

	case eRCCMinimal:
		OSSummary.StartIdx = 1; 
		OSSummary.EndIdx = 6;
		break;

	case eRCCMinimalA:
		OSSummary.StartIdx = 1; 
		OSSummary.EndIdx = 4;
		break;

	case eRCCMinimalB:
		OSSummary.StartIdx = 1; 
		OSSummary.EndIdx = 4;
		break;

	case eRCCMinimalC:
		OSSummary.StartIdx = 0; 
		OSSummary.EndIdx = 1;
		break;
	}


pParams->NumSpecies = 0;
while((Rslt=pCSV->NextLine()) > 0)	// onto next line containing fields
	{
	NumFields = pCSV->GetCurFields();

	if(NumFields < 9)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expected at least 9 fields in outspecies '%s', GetCurFields() returned '%d'",pszFile,NumFields);
		return(eBSFerrFieldCnt);
		}


	pCSV->GetText(2,(char **)&pszElementType);
	strncpy(OSSummary.szElementType,pszElementType,sizeof(OSSummary.szElementType));
	OSSummary.szElementType[sizeof(OSSummary.szElementType)-1] = '\0';
	pCSV->GetText(3,(char **)&pszRefSpecies);
	pCSV->GetText(8,(char **)&pszSpeciesList);

	if(!pParams->NumSpecies)
		{
		pParams->NumSpecies = ParseNumSpecies(pszSpeciesList,pParams);
		strcpy(OSSummary.szRefSpecies,pParams->szSpecies[0]);
		strcpy(OSSummary.szOutSpecies,pParams->szSpecies[pParams->NumSpecies-1]);
		OSSummary.szCoreList[0] = '\0';
		for(int Idx = 1; Idx < (pParams->NumSpecies-1); Idx++)
			{
			if(Idx > 1)
				strcat(OSSummary.szCoreList,",");
			strcat(OSSummary.szCoreList,pParams->szSpecies[Idx]);
			}
		}
	pCSV->GetInt(9,&Region);
	NormRegion = NormaliseRegion(Region);
	if(pParams->Region != eFRAuto && pParams->Region != NormRegion)
		continue;

	pCSV->GetInt(7,&SeqLen);
	pRangeClass = GetLengthRangeClass(SeqLen);

	if(NumFields > 9)		// if base alignment detail in file
		{
		pCSV->GetInt(10,&OGUnalignedBases);
		pCSV->GetInt(11,&OGMatchCnt);
		pCSV->GetInt(12,&OGMismatchCnt);
		pCSV->GetInt(13,&OGInDelCnt);
		AlignedBases = OGMatchCnt + OGMismatchCnt;
		if((AlignedBases >= pParams->Align2Core || AlignedBases == SeqLen) &&
			(((AlignedBases * 100.0) / SeqLen) >= pParams->PCAlign2Core) &&
			(((OGMatchCnt * 100.0) / SeqLen) >= pParams->IDAlign2Core) &&
			(((OGMatchCnt * 100.0) / (OGMatchCnt + OGMismatchCnt)) >= pParams->IDOSIdentity))
			OSSummary.OGAligned2Cores[pRangeClass->ID-1] += 1;
		else
			{
			OGMatchCnt = 0;
			OGMismatchCnt = 0;
			OGInDelCnt = 0;
			OGUnalignedBases = SeqLen;
			NormRegion = -1;
			Region = -1;
			}
		}
	else
		{
		OGUnalignedBases = 0;
		OGMatchCnt = SeqLen;
		OGMismatchCnt = 0;
		OGInDelCnt = 0;
		}
	
	OSSummary.BaseTotals[pRangeClass->ID-1] += SeqLen;
	OSSummary.InstTotals[pRangeClass->ID-1] += 1;
	OSSummary.OGMatchBaseTotals[pRangeClass->ID-1] += OGMatchCnt;
	OSSummary.OGMismatchBaseTotals[pRangeClass->ID-1] += OGMismatchCnt;
	OSSummary.OGInDelBaseTotals[pRangeClass->ID-1] += OGInDelCnt;
	if(NormRegion != -1)			// if regional mapping in source file
		{
		if(pParams->RsltsLayout == eRPRsltRegionTableA)
			{
			// determine exclusive regions (loci exclusive, totally contained, within the region)
			if(Region & cFeatBitCDS && !(Region & (cFeatBit5UTR | cFeatBit3UTR)))
				{
				OSSummary.OGCDSRegion[pRangeClass->ID-1] += 1;
				OSSummary.OGCDSRegionTotals += 1;
				}
			else
				if(Region & cFeatBit5UTR && !(Region & (cFeatBitCDS | cFeatBit3UTR)))
					{
					OSSummary.OG5UTRRegion[pRangeClass->ID-1] += 1;
					OSSummary.OG5UTRRegionTotals += 1;
					}
				else
					if(Region & cFeatBit3UTR && !(Region & (cFeatBitCDS | cFeatBit5UTR)))
						{
						OSSummary.OG3UTRRegion[pRangeClass->ID-1] += 1;
						OSSummary.OG3UTRRegionTotals += 1;
						}
					else
						if(Region & cFeatBitIntrons && !(Region & (cFeatBitCDS | cFeatBit5UTR | cFeatBit3UTR)))
							{
							OSSummary.OGIntronRegion[pRangeClass->ID-1] += 1;
							OSSummary.OGIntronRegionTotals += 1;
							}
						else
							if(Region & cFeatBitUpstream && !(Region & (cFeatBitCDS | cFeatBitIntrons | cFeatBit5UTR | cFeatBit3UTR)))
								{
								OSSummary.OGUSRegion[pRangeClass->ID-1] += 1;
								OSSummary.OGUSRegionTotals += 1;
								}
							else
								if(Region & cFeatBitDnstream && !(Region & (cFeatBitCDS | cFeatBitIntrons | cFeatBit5UTR | cFeatBit3UTR)))
									{
									OSSummary.OGDSRegion[pRangeClass->ID-1] += 1;
									OSSummary.OGDSRegionTotals += 1;
									}
								else
									if(!Region)
										{
										OSSummary.OGIntergenicRegion[pRangeClass->ID-1] += 1;
										OSSummary.OGIntergenicRegionTotals += 1;
										}

			if(Region & (cIntronExonSpliceSite | cExonIntronSpliceSite))
				{
				OSSummary.OGSpliceSites[pRangeClass->ID-1] += 1;
				OSSummary.OGSpliceSitesTotals += 1;
				}


			// now process region counts with region unnormalised
			if(Region & cFeatBitCDS)
				{
				OSSummary.OGRegions[pRangeClass->ID-1][eFRCDS] += 1;
				OSSummary.RegionTotals[eFRCDS] += 1;
				}
			if(Region & cFeatBit5UTR)
				{
				OSSummary.OGRegions[pRangeClass->ID-1][eFR5UTR] += 1;
				OSSummary.RegionTotals[eFR5UTR] += 1;
				}
			if(Region & cFeatBit3UTR)
				{
				OSSummary.OGRegions[pRangeClass->ID-1][eFR3UTR] += 1;
				OSSummary.RegionTotals[eFR3UTR] += 1;
				}
			if(Region & cFeatBitIntrons)
				{
				OSSummary.OGRegions[pRangeClass->ID-1][eFRIntronic] += 1;
				OSSummary.RegionTotals[eFRIntronic] += 1;
				}
			if(Region & cFeatBitUpstream)
				{
				OSSummary.OGRegions[pRangeClass->ID-1][eFRUpstream] += 1;
				OSSummary.RegionTotals[eFRUpstream] += 1;
				}
			if(Region & cFeatBitDnstream)
				{
				OSSummary.OGRegions[pRangeClass->ID-1][NormRegion] += 1;
				OSSummary.RegionTotals[eFRDnstream] += 1;
				}

			if(!Region)
				{
				OSSummary.OGRegions[pRangeClass->ID-1][eFRIntergenic] += 1;
				OSSummary.RegionTotals[eFRIntergenic] += 1;
				}

			if(Region & cIntronExonSpliceSite)
				{
				OSSummary.OGIntronExonSpliceSite[pRangeClass->ID-1] += 1;
				OSSummary.IntronExonSpliceSiteTotals += 1;
				}
			
			if(Region & cExonIntronSpliceSite)
				{		
				OSSummary.OGExonIntronSpliceSite[pRangeClass->ID-1] += 1;
				OSSummary.ExonIntronSpliceSiteTotals += 1;
				}
			}
		else
			{
			// now process region counts where region has been normalised
			OSSummary.OGRegions[pRangeClass->ID-1][NormRegion] += 1;
			OSSummary.RegionTotals[NormRegion] += 1;
			if(Region & cIntronExonSpliceSite)
				{
				OSSummary.OGIntronExonSpliceSite[pRangeClass->ID-1] += 1;
				OSSummary.IntronExonSpliceSiteTotals += 1;
				}
			if(Region & cExonIntronSpliceSite)
				{		
				OSSummary.OGExonIntronSpliceSite[pRangeClass->ID-1] += 1;
				OSSummary.ExonIntronSpliceSiteTotals += 1;
				}
			}
		}

	if(NumFields <= 15)
		OSSummary.NumDistSegs = 0;
	else
		OSSummary.NumDistSegs = (NumFields - 15)/4;		// 4 fields per segment (matches,mismatches,indels,unaligned)
	pParams->NumDistSegs = OSSummary.NumDistSegs;

		// need to normalise the counts because the orginal count binning
		// was on discrete base boundaries. For example -
		// Assume 10 bins, and sequence length was 20
		// All 10 bins contain counts for 2 bases per bin (BinWidth = 2)
		// Assume sequence length was 21
		// The first 9 bins contain counts for 2 bases per bin, but the last bin contains counts for 3 bases
		// Assume sequence length now 28
		// The first 2 bins contain counts for 2 bases per bin (BinWidth = 2), but the last 8 bins contains counts for 3 bases (BinWidth = 3)
		//
		// To normalise, scale the bases per bin to be BinCnt *= Sequencelength/(NumberBins*BinWidth)
		// Using same examples as above -
		// Assume 10 bins, and sequence length was 20
		// Scale all 10 bins by 20/(10*2) or 1.0
		// Assume sequence length was 21
		// Scale first 9 bins by 21/(10*2) or 1.05 , last bin by 21/(10*3) or 0.70
	if(OSSummary.NumDistSegs && (!pParams->bSegsExcludeMatch || SeqLen > OGMatchCnt))
		{
		int BinWidth = SeqLen / OSSummary.NumDistSegs;	// intial bin width or number of bases in that bin
		int BinWidthChg = OSSummary.NumDistSegs - (SeqLen % OSSummary.NumDistSegs); // bin index at which scaling changes
		double CntScale = (double)SeqLen / (OSSummary.NumDistSegs * BinWidth);
		for(SegIdx = 0; SegIdx < OSSummary.NumDistSegs; SegIdx++)
			{
			if(SegIdx == BinWidthChg)
				CntScale = (double)SeqLen / (OSSummary.NumDistSegs * (BinWidth+1));
			pCSV->GetInt((SegIdx * 4) + 16,&SegCnt);
			OSSummary.DistSegCnts[pRangeClass->ID-1][SegIdx].Matches += SegCnt * CntScale;
			pCSV->GetInt((SegIdx * 4) + 17,&SegCnt);
			OSSummary.DistSegCnts[pRangeClass->ID-1][SegIdx].Mismatches += SegCnt * CntScale;
			pCSV->GetInt((SegIdx * 4) + 18,&SegCnt);
			OSSummary.DistSegCnts[pRangeClass->ID-1][SegIdx].InDels += SegCnt * CntScale;
			pCSV->GetInt((SegIdx * 4) + 19,&SegCnt);
			OSSummary.DistSegCnts[pRangeClass->ID-1][SegIdx].Unaligned += SegCnt * CntScale;
			}
		}

	if(OGUnalignedBases == SeqLen)		// if no alignment at all...
		{
		OSSummary.OGUnalignedCores[pRangeClass->ID-1] += 1;
		OSSummary.OGUnalignedCoreBases[pRangeClass->ID-1] += OGUnalignedBases;
		}
	else			// else at least 1 base aligned!
		{
		OSSummary.OGUnalignedBaseTotals[pRangeClass->ID-1] += OGUnalignedBases;
		double MatchPercent = (OGMatchCnt * 100.0) / SeqLen;
		int MatchPercentIdx = (int)(MatchPercent + 0.001);	// 0.001 added in case of small round down errors in floating points
		OSSummary.OGMatchIdentities[pRangeClass->ID-1][MatchPercentIdx] += 1;
		}

	NumRows++;
	}

// all element rows have been processed

OSSummary.InstTotal = 0;
OSSummary.BaseTotal = 0;
for(int Idx = OSSummary.StartIdx; Idx < OSSummary.EndIdx; Idx++)
		{
		OSSummary.InstTotal += OSSummary.InstTotals[Idx];
		OSSummary.BaseTotal += OSSummary.BaseTotals[Idx];
		OSSummary.AlignedInstTotal += OSSummary.OGAligned2Cores[Idx];
		OSSummary.OGMatchBaseTotal += OSSummary.OGMatchBaseTotals[Idx];
		OSSummary.OGMismatchBaseTotal += OSSummary.OGMismatchBaseTotals[Idx];
		OSSummary.OGInDelBaseTotal += OSSummary.OGInDelBaseTotals[Idx];
		OSSummary.OGUnalignedBaseTotal += OSSummary.OGUnalignedBaseTotals[Idx];
		OSSummary.OGUnalignedCoreBaseTotal += OSSummary.OGUnalignedCoreBases[Idx];
		OSSummary.OGUnalignedCoresTotal += OSSummary.OGUnalignedCores[Idx];
		
		for(int IdentityIdx = 0; IdentityIdx < cMaxNumIdentityClasses; IdentityIdx++)
			OSSummary.OGTotMatchIdentities[IdentityIdx] += OSSummary.OGMatchIdentities[Idx][IdentityIdx];

		for(SegIdx = 0; SegIdx < OSSummary.NumDistSegs; SegIdx++)
			{
			OSSummary.DistSegCntsTotals[SegIdx].Matches += OSSummary.DistSegCnts[Idx][SegIdx].Matches;
			OSSummary.DistSegCntsTotals[SegIdx].Mismatches += OSSummary.DistSegCnts[Idx][SegIdx].Mismatches;
			OSSummary.DistSegCntsTotals[SegIdx].InDels += OSSummary.DistSegCnts[Idx][SegIdx].InDels;
			OSSummary.DistSegCntsTotals[SegIdx].Unaligned += OSSummary.DistSegCnts[Idx][SegIdx].Unaligned;
			}
		}
OSSummary.pRangeClass = pLenRangeClasses;

switch(pParams->RsltsLayout) {
	case eRPRsltStandard:
		Rslt = RsltStandard(&OSSummary,pParams);
		break;
	case eRPRsltTable:
		Rslt = RsltTable(&OSSummary,pParams);
		break;
	case eRPRsltSegDistRows:
		Rslt = RsltSegDistRows(&OSSummary,pParams);
		break;
	case eRPRsltSegDistTables:
		Rslt = RsltSegDistTables(&OSSummary,pParams);
		break;
	case eRPRsltSegDistCols:
		Rslt = RsltSegDistCols(&OSSummary,pParams);
		break;
	case eRPRsltIdentDistCols:
		Rslt = RsltIdentDistCols(&OSSummary,pParams);
		break;
	case eRPRsltRegionTable:
		Rslt = RsltRegionTable(&OSSummary,pParams);
		break;
	case eRPRsltRegionTableA:
		Rslt = RsltRegionTableA(&OSSummary,pParams);
		break;
	}
return(eBSFSuccess);
}



