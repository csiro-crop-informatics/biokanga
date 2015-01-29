// CSVFilter.cpp : Defines the entry point for the console application.


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

const char *cpszProgVer = "1.0.1";		// increment with each release

const int cMinElLen = 0;				// allow element lengths to be specified down to cMinElLen
const int cDfltElLen= 20;				// if element lengths not specified then default to cDfltMinElLen
const int cMaxElLen = 1000000;			// minimum element lengths can be specified upto this length

const int cMinDeltaLen = 0;				// allow 0% change in length between ref + outspecies
const int cDfltDeltaLen= 5;				// default is a 5% change in length between ref + outspecies
const int cMaxDeltaLen = 100;			// ignore any differences in length 

const int cMinIdentity = 100;			// require 100% identity between ref + outspecies
const int cDfltIdentity= 95;			// default is that there must be a 95% identity between ref + outspecies
const int cMaxIdentity = 0;				// ignore any differences in identity

const int cMaxIncludeLoci = 10;			// maximun number of include loci filter files
const int cMaxExcludeLoci = 10;			// maximun number of exclude loci filter files

const int cMaxIncludeChroms = 20;		// max number of include chromosomes regular expressions
const int cMaxExcludeChroms = 20;		// max number of exclude chromosomes regular expressions

const int cMaxOutputFields = 20;		// max number of field output specifiers

const int cMaxNumFields = 200;			// assume that there will never be more than 200 fields in CSV file

const int cElLociAllocChunk = 5000000;	// grow core loci array by this many elements
const int cElChrAllocChunk = 100000000;	// grow chr buffer by this size
const int cElChromAllocChunk = 10000;	// grow chrom array by this many elements
const int cElSpeciesAllocChunk = 100;	// grow species array by this many elements


// processing mode
typedef enum eProcMode {
	eProcStandard,						// standard processing
	eProcOutspecies						// process outspecies CSV files containing outspecies identities, indels etc
} etProcMode;


typedef struct TAG_sFiltState {
	bool bFiltOut;						// true if this element has been marked to be filtered out
	unsigned int fOverLen:1;			// element is overlength
	unsigned int fUnderLen:1;			// element is underlength
	unsigned int fExcludeLoci:1;		// element is on loci to be excluded
	unsigned int fIncludeLoci:1;		// element is not on loci to be included

	unsigned int fOverlap:1;			// element overlaps another element
	unsigned int fChrom:1;				// element is on filtered chrom
	unsigned int fOutRefID:1;			// element filtered RefID
	unsigned int fInRefID:1;			// element filtered RefID
	unsigned int fSpecies:1;			// element filtered species

	unsigned int fDeltaIdentity:1;		// element filtered identity
	unsigned int fDeltaLen:1;			// element filtered delta length
	unsigned int fExcludeRegion:1;		// element filtered because in regions to be filtered out
	unsigned int fIncludeRegion:1;		// element filtered because not in regions to be retained
	unsigned int fSelectN:1;			// element filtered to reduce number of rows output

	} tsFiltState;

typedef struct TAG_sCSVEl {
	tsFiltState Filter;
	int RefID;							// element reference identifier as parsed from source file
	int SpeciesID;						// species
	int ChromID;						// chromsome
	int StartLoci;						// loci start
	int EndLoci;						// loci end
	int LineLen;						// strlen() of line 
	size_t LineOfs;						// offset into pElChrBuff of line as it would be written, -1 if no line to write 
} tsCSVEl;

typedef struct TAG_sElChrom {
	int ChromID;							// unique identifier for this chromosome
	UINT16 Hash;							// hash over chromosome name
	char szChrom[cMaxDatasetSpeciesChrom];	// chromosome
} tsElChrom;

typedef struct TAG_sElSpecies {
	int SpeciesID;							// unique identifier for this species
	UINT16 Hash;							// hash over species name
	char szSpecies[cMaxDatasetSpeciesChrom];	// species
} tsElSpecies;

typedef struct TAG_sExcludeSpeciesChrom {
	char szSpecies[cMaxDatasetSpeciesChrom];			// which species
	char szChromID[cMaxDatasetSpeciesChrom];			// chromosome is not to be processed
	} tsExcludeSpeciesChrom;


const int cMaxExcludeHistory = 100;
typedef struct TAG_sExcludeEl {
	struct TAG_sExcludeEl *pNext;
	struct TAG_sExcludeEl *pPrev;
	char szSpecies[cMaxDatasetSpeciesChrom];			// which species
	char szChrom[cMaxDatasetSpeciesChrom];			// chromosome is not to be processed
	bool bExclude;		// true if to be excluded, false if not
	} tsExcludeEl;

typedef struct TAG_sProcParams 
	{
	etProcMode ProcMode;			// processing mode
	char *pszSpeciesList;			// comma separated species list

	int NumSpecies;					// number of species
	char szSpecies[cMaxAlignedSpecies][cMaxDatasetSpeciesChrom];	// species names of interest

	bool bOverlaps;				// filter out overlaps
	bool bNoOverlaps;			// filter out none overlaps

	int RegionsIn;				// exclusive regions to be retained
	int RegionsOut;				// regions to be filtered out

	char *pszInFile;				// input CSV loci file
	CCSVFile *pCSV;

	char *pszStatsFile;				// output stats file to create/write
	int hStatsFile;					// file handle for stats file

	char *pszOutFile;				// output filtered file to create/write
	int hOutFile;					// file handle for output filtered file

	int SelectN;					// max number of output rows into output file

	int NumIncludeChroms;			// number of chromosomes explicitly defined to be included
	char **ppszIncludeChroms;		// ptr to array of reg expressions defining chroms to include - overides exclude
	int NumExcludeChroms;			// number of chromosomes explicitly defined to be excluded
	char **ppszExcludeChroms;		// ptr to array of reg expressions defining chroms to include
#ifdef _WIN32
	Regexp *IncludeChromsRE[cMaxIncludeChroms];	// compiled regular expressions
	Regexp *ExcludeChromsRE[cMaxExcludeChroms];
#else
	regex_t IncludeChromsRE[cMaxIncludeChroms];	// compiled regular expressions
	regex_t ExcludeChromsRE[cMaxExcludeChroms];
#endif
	int	MinElLen;					// minimum element length required
	int	MaxElLen;					// maximum element length required

	int Align2Core;					// at least this many bases aligned in outspecies for it to count as aligned to core
	double PCAlign2Core;				// minimum percentage (0..100) aligned to core to count as as aligned to core (default is 0%)");
	double IDIdent2Core;				// minimum identity (0..100) aligned to core to count as as aligned to core (default is 0%)");
	double IDOSIdentity;				// minimum outspecies (matches/matchesPlusmismatches) (0..100) aligned to core to count as as aligned to core (default is 1%)");

	CFilterRefIDs *pFilterInRefIDs;	// used when filtering in by RefID
	CFilterRefIDs *pFilterOutRefIDs;// used when filtering out by RefID

	int gNumExcludeEls;		// current number of elements in gExcludeChroms
	tsExcludeEl *gpMRA;		// pts to most recently accessed or added
	tsExcludeEl *gpLRA;		// pts to least recently accessed
	tsExcludeEl gExcludeChroms[cMaxExcludeHistory];

	CFilterLoci *pExcludeFiltLoci;	// exclude these loci
	CFilterLoci *pIncludeFiltLoci;	// include these loci



	tsCSVEl *pElLocs;				// pts to mem allocd to hold array of loci (was allocated with malloc/realloc, or mmap/mremap, not c++'s new...)
	size_t MaxElLocs;					// max allocd locii
	size_t NumElLocs;					// number of loci ptd at by pElLocs

	char *pElChrBuff;				// pts to mem allocd (was allocated with malloc/realloc, or mmap/mremap, not c++'s new...) to hold buffered lines to output
	size_t MaxAllocdChrs;			// allocd number of chrs in pElChrBuff
	size_t NumUsedChrs;				// number of loci ptd at by pElLocs

	int CachElSpeciesID;			// identifier of last retrieved species from pElSpecies
	tsElSpecies *pElSpecies;			// pts to mem allocd to hold array of species names
	int MaxElSpecies;				// max allocd species
	int NumElSpecies;				// number of species ptd at by pElSpecies


	int CachElChromID;				// identifier of last retrieved chromosome from pElChroms
	tsElChrom *pElChroms;			// pts to mem allocd to hold array of chromosome names
	int MaxElChroms;				// max allocd chroms
	int NumElChroms;				// number of chroms ptd at by pElSpecies

	bool bFilterOutputFields;			// true if output fields to be filtered according to OutputFields[];
	bool OutputFields[cMaxNumFields];	// true if field is to be output otherwise false 
} tsProcParams; 




int Process(etProcMode ProcMode,		// processing mode
			char *pszInFile,			// file to filter
			char *pszOutFile,			// write filtered file
			char *pszStatsFile,			// write stats file
			int FilterRegionsIn,		// retain any of these (exclusive) regions
			int FilterRegionsOut,		// remove any of these regions
			int SelectN,				// output at most this many rows
			char *pszOutputFields,		// only output these fields
			char *pszSpeciesList,		// list of species to filter in ('\0' if no filtering)
			bool bOverlaps,				// filter out overlaps
			bool bNoOverlaps,			// filter out none overlaps
			char *pszFilterOutRefIDFile, // RefIDs to filter out
			char *pszFilterInRefIDFile,	 // RefIDs to filter in unless filtered out
			int MinElLen,				 // elements must be of at least this length
			int MaxElLen,				 // no longer than this length
			int Align2Core,			// at least this many bases aligned in outspecies for it to count as aligned to core
			double PCAlign2Core,		// minimum percentage (0..100) aligned to core to count as as aligned to core (default is 0%)");
			double IDIdent2Core,		// minimum identity (0..100) aligned to core to count as as aligned to core (default is 0%)");
			double IDOSIdentity,		// minimum outspecies (matches/matchesPlusmismatches) (0..100) aligned to core to count as as aligned to core (default is 1%)");
			int NumIncludeChroms,		 // number of chromosome regular expressions to include
			char **ppszIncludeChroms,	 // array of include chromosome regular expressions
			int NumExcludeChroms,		 // number of chromosome expressions to exclude
			char **ppszExcludeChroms,	 // array of exclude chromosome regular expressions
			int NumIncludeLoci,			 // number of include region files 
			char **ppszIncludeLoci,		 // array of include region files
			int NumExcludeLoci,			 // number of exclude region files
			char **ppszExcludeLoci);	 // array of exclude region files

char *ProcMode2Txt(etProcMode ProcMode);
int TrimQuotes(char *pszTxt);
int ParseNumFields(char *pszFieldList,tsProcParams *pProcParams);
int ParseNumSpecies(char *pszSpeciesList,tsProcParams *pProcParams);
int ParseRegions(char *pszRegionList);
char *Regions2Txt(int Regions);
bool ExcludeThisChrom(char *pszSpecies,char *pszChrom,tsProcParams *pProcParams);
tsExcludeEl *LocateExclude(char *pszSpecies,char *pszChrom,tsProcParams *pParams);
bool AddExcludeHistory(char *pszSpecies,char *pszChrom,bool bExclude,tsProcParams *pParams);
int FilterCSV(tsProcParams *pProcParams);

UINT16 GenNameHash(char *pszName);
tsElChrom *LocateElChrom(char *pszChrom,tsProcParams *pProcParams);
int AddChrom(char *pszChrom,tsProcParams *pParams);

tsElSpecies *LocateElSpecies(char *pszSpecies,tsProcParams *pProcParams);
int AddSpecies(char *pszSpecies,tsProcParams *pParams);

tsCSVEl *
AddCSVEl(tsFiltState FiltState,	// current filter state for this element
		 int RefID,				// reference identifier from input file
		 char *pszSpecies,		// species
		 char *pszChrom,		// chrom
		 int StartLoci,			// start loci
		 int EndLoci,			// end loci
		 int LineLen,			// num of chars (excl trailing '\0') in pszLineBuff
		 char *pszLine,			// line to output if not filtered
		 tsProcParams *pParams);

static int CompareCSVEls( const void *arg1, const void *arg2);


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
	return _T("CSVFilter");
}
// end of str library required code
#endif

CStopWatch gStopWatch;
CDiagnostics gDiagnostics;				// for writing diagnostics messages to log file
char gszProcName[_MAX_FNAME];			// process name

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
int LenFileList;
int iNumSpecies;

int iProcMode;

char szInFile[_MAX_PATH];
char szOutFile[_MAX_PATH];
char szStatsFile[_MAX_PATH];
char szFilterOutRefIDFile[_MAX_PATH];
char szFilterInRefIDFile[_MAX_PATH];
char szSpeciesList[cMaxAlignedSpecies*cMaxDatasetSpeciesChrom];
char szRegionsIn[128];
int iRegionsIn;
int iRegionsOut;
char szRegionsOut[128];
int iMinElLen;
int iMaxElLen;
int NumIncludeChroms;
char *pszIncludeChroms[cMaxIncludeChroms];
int NumExcludeChroms;
char *pszExcludeChroms[cMaxExcludeChroms];
int NumIncludeLoci;
char *pszIncludeLoci[cMaxIncludeLoci];
int NumExcludeLoci;
char *pszExcludeLoci[cMaxExcludeLoci];
bool bNoOverlaps;
bool bOverlaps;

int iAlign2Core;					// at least this many bases aligned in outspecies for it to count as aligned to core
double dPCAlign2Core;				// minimum percentage (0..100) aligned to core to count as as aligned to core (default is 0%)");
double dIDIdent2Core;				// minimum identity (0..100) aligned to core to count as as aligned to core (default is 0%)");
double dIDOSIdentity;				// minimum outspecies (matches/matchesPlusmismatches) (0..100) aligned to core to count as as aligned to core (default is 1%)");

int iSelectN;						// randomly select at most this many rows to return
int NumOutputFields;
char szOutputFields[cMaxOutputFields * 10];

// command line args
struct arg_lit  *help    = arg_lit0("h","help",                 "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to screen and logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"diagnostics log file");

struct arg_int  *ProcMode = arg_int0("m","procmode","<int>",		"mode 0: filter on loci only, mode 1: filter on identity and delta length");
struct arg_file *InFile = arg_file1("i",NULL,"<file>",				"input csv file to be filtered");
struct arg_file *OutFile = arg_file0("o",NULL,"<file>",				"optional output filtered csv file");
struct arg_file *StatsFile = arg_file0("O",NULL,"<file>",			"optional output stats file");
struct arg_str  *RegionsOut = arg_str0("R","regionsout","<string>","Remove these regions (space or comma delimit), 1: Intergenic, 2: US, 3: 5'UTR, 4: CDS, 5: Intron, 6: 3'UTR, 7: DS, 8: 5'Splice, 9: 3'Splice");
struct arg_str  *RegionsIn = arg_str0("r","regionsin","<string>","Retain these exclusive regions (space or comma delimit), 1: Intergenic, 2: US, 3: 5'UTR, 4: CDS, 5: Intron, 6: 3'UTR, 7: DS");
struct arg_str  *SpeciesList = arg_str0("s","includespecies","<string>","Filter in only these species");
struct arg_lit  *Overlaps = arg_lit0("j","nooverlaps",				"filter out all loci which is overlapped by another loci");
struct arg_lit  *NoOverlaps = arg_lit0("J","overlaps",				"filter out all loci which is not overlapped by another loci");
struct arg_file *FilterOutRefIDFile = arg_file0("X",NULL,"<file>",	"First: filter out with RefIDs from this filter file");
struct arg_file *FilterInRefIDFile = arg_file0("x",NULL,"<file>",	"Second: filter in with RefIDs from this filter file");
struct arg_int *MinElLen=arg_int0("l", "MinElLen",	"<int>",		"filter on minimum absolute length");
struct arg_int *MaxElLen=arg_int0("L", "MaxElLen",	"<int>",		"filter on maximum absolute length");

struct arg_int  *Align2Core = arg_int0("a","align2core","<int>","minimum bases in outspecies for alignment to count as aligned to core (default is 0)");
struct arg_dbl  *IDAlign2Core = arg_dbl0("P","pcalign2core","<dbl>","minimum percentage bases aligned to core to count as as aligned to core (default is 0%)");
struct arg_dbl  *IDIdent2Core = arg_dbl0("A","IDIdent2Core","<dbl>","minimum identity core to count as as aligned to core (default is 0%)");
struct arg_dbl  *IDOSIdentity = arg_dbl0("k","osidentity","<dbl>","minimum outspecies identity ,matches/(matches+mismatches), to count as as aligned to core (default is 0%)");

struct arg_file *ExcludeLoci = arg_filen("E","exclude","<file>",0,cMaxExcludeLoci,	"First: filter out these Loci");
struct arg_file *IncludeLoci = arg_filen("I","include","<file>",0,cMaxIncludeLoci,	"Second: filter in these Loci");
struct arg_str  *ExcludeChroms = arg_strn("Z","chromexclude","<string>",0,cMaxExcludeChroms,"priority 1: regular expressions defining species.chromosomes to explicitly exclude");
struct arg_str  *IncludeChroms = arg_strn("z","chromeinclude","<string>",0,cMaxIncludeChroms,"priority 2: regular expressions defining species.chromosomes to explicitly include");

struct arg_int *SelectN=arg_int0("N", "SelectN","<int>",			"Randomly select at most N (without relacement) rows to return - default is to return all");
struct arg_str  *OutputFields = arg_str0("n","fieldsoutput","<string>","Output these fields only (default is to output all fields)");

struct arg_end *end = arg_end(20);

void *argtable[] = {help,version,FileLogLevel,LogFile,
					ProcMode,
					InFile,OutFile,StatsFile,RegionsOut,RegionsIn,SpeciesList,NoOverlaps,Overlaps,FilterOutRefIDFile,FilterInRefIDFile,MinElLen,MaxElLen,
					Align2Core,IDAlign2Core,IDIdent2Core,IDOSIdentity,
					ExcludeLoci,IncludeLoci,ExcludeChroms,IncludeChroms,SelectN,OutputFields,
					end};

char **pAllArgs;
int argerrors;
argerrors = CUtility::arg_parsefromfile(argc,(char **)argv,&pAllArgs);
if(argerrors >= 0)
	argerrors = arg_parse(argerrors,pAllArgs,argtable);

    /* special case: '--help' takes precedence over error reporting */
if (help->count > 0)
        {
		printf("\n%s CSV File Filter, Version %s\nOptions ---\n", gszProcName,cpszProgVer);
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

	iProcMode = ProcMode->count ? ProcMode->ival[0] : eProcStandard;
	if(iProcMode < eProcStandard || iProcMode > eProcOutspecies)
		{
		printf("\nError: Requested processing mode '-x%d' not supported",iProcMode);
		exit(1);
		}

	if(RegionsIn->count)
		{
		strcpy(szRegionsIn,RegionsIn->sval[0]);
		TrimQuotes(szRegionsIn);
		if((iRegionsIn = ParseRegions(szRegionsIn)) < 0)
			{
			printf("Error: unable to parse '-r%s' into regions to retain",szRegionsIn);
			exit(1);
			}
		if(iRegionsIn & (cIntronExonSpliceSite | cExonIntronSpliceSite))
			{
			printf("Error: Can't specify 5'Splice and 3'Splice to be exclusively retained with '-r%s'",szRegionsIn);
			exit(1);
			}
		}
	else
		{
		szRegionsIn[0] = '\0';
		iRegionsIn = 0;
		}

	if(RegionsOut->count)
		{
		strcpy(szRegionsOut,RegionsOut->sval[0]);
		TrimQuotes(szRegionsOut);
		if((iRegionsOut = ParseRegions(szRegionsOut)) < 0)
			{
			printf("Error: unable to parse '-R%s' into regions to remove",szRegionsOut);
			exit(1);
			}
		}
	else
		{
		szRegionsOut[0] = '\0';
		iRegionsOut = 0;
		}

	if(SpeciesList->count)
		{
		strcpy(szSpeciesList,SpeciesList->sval[0]);
		TrimQuotes(szSpeciesList);
		iNumSpecies = ParseNumSpecies(szSpeciesList,NULL);
		if(iNumSpecies < 1)
			{
			printf("Error: Species list %s",iNumSpecies < 1 ? "is unable to be parsed" : "must contain at least 1 species");
			exit(1);
			}
		else
			if(iNumSpecies > cMaxAlignedSpecies)
				{
				printf("Error: Species list '-s%s' specifies (%d) more than %d species",szSpeciesList,iNumSpecies,cMaxAlignedSpecies);
				exit(1);
				}
		}
	else
		szSpeciesList[0] = '\0';

	bOverlaps = Overlaps->count ? true : false;
	bNoOverlaps = NoOverlaps->count ? true : false;
	if(bOverlaps && bNoOverlaps)
		{
		printf("\nError: Overlaps '-j' and NoOverlaps '-J' are exclusive and can't both be specified");
		exit(1);
		}

	iSelectN = SelectN->count ? SelectN->ival[0] : INT_MAX;
	if(iSelectN < 1)
		{
		printf("\n '-N%d' < 1, assuming you meant '-N1",iSelectN);
		iSelectN = 1;
		}

	if(OutputFields->count)
		{
		strncpy(szOutputFields,OutputFields->sval[0],sizeof(szOutputFields));
		TrimQuotes(szOutputFields);
		NumOutputFields = ParseNumFields(szOutputFields,NULL);
		if(NumOutputFields < 1)
			{
			printf("\nError: No output fields specified with '-nFields'");
			exit(1);
			}
		}
	else
		{
		szOutputFields[0] = '\0';
		NumOutputFields = 0;
		}

	iMinElLen = MinElLen->count ? MinElLen->ival[0] : cDfltElLen;
	if(iMinElLen < cMinElLen)
		{
		printf("\nSpecified minimum element length '-l%d' < %d, assuming you meant '-l%d'",iMinElLen,cMinElLen,cMinElLen);
		iMinElLen = cMinElLen;
		}
	else
		{
		if(iMinElLen > cMaxElLen)
			{
			printf("\nSpecified minimum element length '-l%d' > %d, assuming you meant '-l%d'",iMinElLen,cMaxElLen,cMaxElLen);
			iMinElLen = cMaxElLen;
			}
		}


	iMaxElLen = MaxElLen->count ? MaxElLen->ival[0] : cMaxElLen;
	if(iMaxElLen < iMinElLen)
		{
		printf("\nSpecified maximum hyper length '-L%d' < %d, assuming you meant '-L%d'",iMaxElLen,iMinElLen,iMinElLen);
		iMaxElLen = iMinElLen;
		}
	else
		{
		if(iMaxElLen > cMaxElLen)
			{
			printf("\nSpecified maximum element length was '-L%d' > %d, assuming you meant '-L%d'",iMaxElLen,cMaxElLen,cMaxElLen);
			iMaxElLen = cMaxElLen;
			}
		}

	if(iProcMode == eProcOutspecies)
		{
		iAlign2Core = Align2Core->count ? Align2Core->ival[0] : 0;
		if(iAlign2Core < 0)
			{
			printf("\nError: Aligned bases to core '-a%d' must be >= 0",iAlign2Core);
			exit(1);
			}

		dPCAlign2Core = IDAlign2Core->count ? IDAlign2Core->dval[0] : 0.0;
		if(dPCAlign2Core < 0.0 || dPCAlign2Core > 100.0)
			{
			printf("\nError: Minimum percentage aligned bases to core '-P%f' must be in range 0.0 to 100.0",dPCAlign2Core);
			exit(1);
			}
			
		dIDIdent2Core = IDIdent2Core->count ? IDIdent2Core->dval[0] : 0.0;
		if(dIDIdent2Core < 0.0 || dIDIdent2Core > 100.0)
			{
			printf("\nError: Minimum identity to core '-A%f' must be in range 0.0 to 100.0",dIDIdent2Core);
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
		iAlign2Core = 0;
		dIDIdent2Core = 0.0;
		dPCAlign2Core = 0.0;
		dIDOSIdentity = 0.0;
		}

	if(FilterOutRefIDFile->count)
		{
		strncpy(szFilterOutRefIDFile,FilterOutRefIDFile->filename[0],sizeof(szFilterOutRefIDFile));
		szFilterOutRefIDFile[sizeof(szFilterOutRefIDFile)-1] = '\0';
		}
	else
		szFilterOutRefIDFile[0] = '\0';


	if(FilterInRefIDFile->count)
		{
		strncpy(szFilterInRefIDFile,FilterInRefIDFile->filename[0],sizeof(szFilterInRefIDFile));
		szFilterInRefIDFile[sizeof(szFilterInRefIDFile)-1] = '\0';
		}
	else
		szFilterInRefIDFile[0] = '\0';


	NumIncludeLoci = IncludeLoci->count;
	for(Idx=0;Idx < IncludeLoci->count; Idx++)
		{
		LenFileList = (int)strlen(IncludeLoci->filename[Idx]);
		pszIncludeLoci[Idx] = new char [LenFileList+1];
		strcpy(pszIncludeLoci[Idx],IncludeLoci->filename[Idx]);
		TrimQuotes(pszIncludeLoci[Idx]);
		}
	NumExcludeLoci = ExcludeLoci->count;
	for(Idx=0;Idx < ExcludeLoci->count; Idx++)
		{
		LenFileList = (int)strlen(ExcludeLoci->filename[Idx]);
		pszExcludeLoci[Idx] = new char [LenFileList+1];
		strcpy(pszExcludeLoci[Idx],ExcludeLoci->filename[Idx]);
		TrimQuotes(pszExcludeLoci[Idx]);
		}

	NumIncludeChroms = IncludeChroms->count;
	for(Idx=0;Idx < IncludeChroms->count; Idx++)
		{
		LenFileList = (int)strlen(IncludeChroms->sval[Idx]);
		pszIncludeChroms[Idx] = new char [LenFileList+1];
		strcpy(pszIncludeChroms[Idx],IncludeChroms->sval[Idx]);
		TrimQuotes(pszIncludeChroms[Idx]);
		}

	NumExcludeChroms = ExcludeChroms->count;
	for(Idx=0;Idx < ExcludeChroms->count; Idx++)
		{
		LenFileList = (int)strlen(ExcludeChroms->sval[Idx]);
		pszExcludeChroms[Idx] = new char [LenFileList+1];
		strcpy(pszExcludeChroms[Idx],ExcludeChroms->sval[Idx]);
		TrimQuotes(pszExcludeChroms[Idx]);
		}

	strncpy(szInFile,InFile->filename[0],sizeof(szInFile));
	if(OutFile->count)
		strncpy(szOutFile,OutFile->filename[0],sizeof(szOutFile));
	else
		szOutFile[0] = '\0';

	if(StatsFile->count)
		strncpy(szStatsFile,StatsFile->filename[0],sizeof(szStatsFile));
	else
		szStatsFile[0] = '\0';

			// now that command parameters have been parsed then initialise diagnostics log system
	if(!gDiagnostics.Open(szLogFile,(etDiagLevel)iScreenLogLevel,(etDiagLevel)iFileLogLevel,true))
		{
		printf("\nError: Unable to start diagnostics subsystem.");
		if(szLogFile[0] != '\0')
			printf(" Most likely cause is that logfile '%s' can't be opened/created",szLogFile);
		exit(1);
		}

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Version: %s Processing parameters:",cpszProgVer);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Processing mode: %s",ProcMode2Txt((etProcMode)iProcMode));

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"input loci csv file to be filtered: %s",szInFile);
	if(szOutFile[0] != '\0')
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"output filtered loci csv file: %s",szOutFile);
	if(szStatsFile[0] != '\0')
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"output stats csv file: %s",szStatsFile);

	if(iRegionsIn)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Filter in these regions: '%s'",Regions2Txt(iRegionsIn));
	if(iRegionsOut)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Filter out these regions: '%s'",Regions2Txt(iRegionsOut));
	
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"At most output this many rows: %d",iSelectN);

	if(NumOutputFields)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Output only these fields: %s",szOutputFields);
	else
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Output all fields");

	if(szSpeciesList[0] != '\0')
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Filter out any species not in: %s",szSpeciesList);
	else
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Accept all species");

	if(bOverlaps)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Filter out any overlapping loci");

	if(bNoOverlaps)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Filter out any none overlapping loci");

	if(szFilterOutRefIDFile[0])
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"filter out with RefIDs from this filter file: %s",szFilterOutRefIDFile);
	if(szFilterOutRefIDFile[0])
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"filter in with RefIDs from this filter file: %s",szFilterInRefIDFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"elements must be at least this length: %d",iMinElLen);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"elements must be no longer than this length: %d",iMaxElLen);

	if(iProcMode == eProcOutspecies)
		{
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Alignments must have at least : %d bases",iAlign2Core);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Minimum percentage bases must be least : %1.2f%% ",dPCAlign2Core);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Minimum identity to core at least : %1.2f%%",dIDIdent2Core);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Minimum outspecies identity at least : %1.2f%%",dIDOSIdentity);
		}

	for(Idx = 0; Idx < NumIncludeChroms; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"reg expressions defining chroms to include: '%s'",pszIncludeChroms[Idx]);
	for(Idx = 0; Idx < NumExcludeChroms; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"reg expressions defining chroms to exclude: '%s'",pszExcludeChroms[Idx]); 

	for(Idx = 0; Idx < NumIncludeLoci; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"files with loci to include: '%s'",pszIncludeLoci[Idx]);
	for(Idx = 0; Idx < NumExcludeLoci; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"files with loci to include: '%s'",pszExcludeLoci[Idx]); 


		gStopWatch.Start();
#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	Rslt = Process((etProcMode)iProcMode,				// processing mode
			szInFile,			// file to filter
			szOutFile,			// write filtered file
			szStatsFile,		// write stats file
			iRegionsIn,			// regions to retain
			iRegionsOut,		// regions to remove
			iSelectN,			// output at most this many rows
			szOutputFields,		// only output these fields
			szSpeciesList,		// list of species to accept
			bOverlaps,			// filter out overlaps
			bNoOverlaps,		// filter out none overlaps
			szFilterOutRefIDFile, // RefIDs to filter out
			szFilterInRefIDFile, // RefIDs to filter in unless filtered out
			iMinElLen,			 // elements must be of at least this length
			iMaxElLen,			// no longer than this length
			iAlign2Core,		// at least this many bases aligned in outspecies for it to count as aligned to core
			dPCAlign2Core,		// minimum percentage (0..100) aligned to core to count as as aligned to core (default is 0%)");
			dIDIdent2Core,		// minimum identity (0..100) aligned to core to count as as aligned to core (default is 0%)");	
			dIDOSIdentity,		// minimum outspecies (matches/matchesPlusmismatches) (0..100) aligned to core to count as as aligned to core (default is 1%)");
			NumIncludeChroms,	// number of chromosome regular expressions to include
			pszIncludeChroms,	// array of include chromosome regular expressions
			NumExcludeChroms,	// number of chromosome expressions to exclude
			pszExcludeChroms,	// array of exclude chromosome regular expressions
			NumIncludeLoci,		// number of include region files 
			pszIncludeLoci,	// array of include region files
			NumExcludeLoci,		// number of exclude region files
			pszExcludeLoci);	// array of exclude region files
		gStopWatch.Stop();
	Rslt = Rslt < 0 ? 1 : 0;
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Exit Code: %d Total processing time: %s",Rslt,gStopWatch.Read());
	exit(Rslt);
	}
else
	{
	printf("\n%s CSV File Filter, Version %s\n",gszProcName,cpszProgVer);
	arg_print_errors(stdout,end,gszProcName);
	arg_print_syntax(stdout,argtable,"\nUse '-h' to view option and parameter usage\n");
	exit(1);
	}
}

char *
ProcMode2Txt(etProcMode ProcMode)
{
switch(ProcMode) {
	case eProcStandard:
		return((char *)"Processing loci files");
	case eProcOutspecies:
		return((char *)"Processing outspecies files with identity counts");
	default:
		break;
	}
return((char *)"Unsupported");
}

// TrimQuotes
// Removes any leading and trailing quotes from pszTxt
int
TrimQuotes(char *pszTxt)
{
char *pDst = pszTxt;
char Chr;
int Len = 0;
while((Chr = *pszTxt++))	
	{
	if((!Len || *pszTxt == '\0') && (Chr == '"' || Chr == '\''))
		continue;
	*pDst++ = Chr;
	Len++;
	}
*pDst = Chr;
return(Len);
}


// ParseNumFields
// Initialises pProcParams with parsed fields which are to be output
// Expects field identifiers as space or comma delimited list ptd at by pszFieldList
// Returns number of fields
int
ParseNumFields(char *pszFieldList,tsProcParams *pProcParams)
{
// parse out field list
char Chr;
int FieldID = 0;
int NumFields = 0;
bool InToken = false;
if(pszFieldList == NULL || *pszFieldList == '\0')
	return(0);

if(pProcParams != NULL)
	memset(pProcParams->OutputFields,0,sizeof(bool) * cMaxNumFields);	// default is that no fields to be output

while(Chr = *pszFieldList++)
	{
	if(Chr == '"' || Chr == '\'') // change any single or double quotes into spaces
		Chr = ' ';
	if(isspace(Chr) || Chr==',')
		{
		if(!InToken)			// slough whitespace or ',' if not inside a token parse
			continue;
		InToken = false;		// was in token but no longer
		if(FieldID < 1 || FieldID > cMaxNumFields)
			{
			printf("\nError: Parsed field %d outside expected range 1..%d\n",FieldID,cMaxNumFields);
			return(-1);
			}

		if(pProcParams != NULL)
			pProcParams->OutputFields[FieldID-1] = true;
		FieldID = 0;
		NumFields++;
		if(NumFields > cMaxNumFields)
			break;
		continue;
		}
	
	// char not a space or ','
	if(!InToken)			// if not already inside token then start token 
		InToken = true;
	if(Chr >= '0' && Chr <= '9')
		FieldID = (FieldID * 10) + Chr - '0';
	else
		{
		printf("\nError: Illegal char '%c' in field list",Chr);
		return(-1);
		}
	}
if(InToken)
	{
	if(FieldID < 1 || FieldID > cMaxNumFields)
		{
		printf("\nError: Parsed field %d outside expected range 1..%d\n",FieldID,cMaxNumFields);
		return(-1);
		}
	if(pProcParams != NULL)
		pProcParams->OutputFields[FieldID-1] = true;
	NumFields++;
	}
if(pProcParams != NULL)
	pProcParams->bFilterOutputFields = NumFields ? true : false;
return(NumFields);
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
		if(NumSpecies >= cMaxAlignedSpecies)
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



// ParseRegions
// Parses space or comma delimited list of regions in which
// 1 == Intergenic, 2 == US, 3 == 5'UTR, 4 == CDS, 5 == Intron, 6 == 3'UTR, 7 == DS, 8 == 5'Splice, 9 == 3'Splice
//
// Returns bitmap of regions or -1 if parse errors
int
ParseRegions(char *pszRegionList)
{
// parse out region list
char Chr;
int Region = 0;
if(pszRegionList == NULL || *pszRegionList == '\0')
	return(-1);

while(Chr = *pszRegionList++) {
	if(isspace(Chr) || Chr == ',')		// accept spaces and commas as separators
		continue;
	switch(Chr) {
		case '1':						// intergenic to be filtered
			Region |= cFeatBitIG;
			break;
		case '2':						// 5'US to be filtered
			Region |= cFeatBitUpstream;
			break;
		case '3':						// 5'UTR to be filtered
			Region |= cFeatBit5UTR;
			break;
		case '4':
			Region |= cFeatBitCDS;		// CDS to be filtered
			break;
		case '5':
			Region |=  cFeatBitIntrons;	// any intronic to be filtered
			break;
		case '6':
			Region |=  cFeatBit3UTR;	// any 3'UTR to be filtered
			break;
		case '7':
			Region |=  cFeatBitDnstream;	// any 3'DS to be filtered 	
			break;
		case '8':
			Region |=  cIntronExonSpliceSite;	// any 5' splice site to be filtered
			break;
		case '9':
			Region |=  cExonIntronSpliceSite; // any 3' splice site to be filtered
			break;
		default:
			return(-1);
		}
	}
return(Region);
}


// Regions2Txt
// Returns textual representation of regions
char *
Regions2Txt(int Regions)
{
static char szRegions[200];
if(!Regions)
	return((char *)"None specified");
if(Regions & cFeatBitIG || Regions == 0)
	strcpy(szRegions,"Intergenic");
else
	szRegions[0] = '\0';
if(Regions & cFeatBitUpstream)
	{
	if(szRegions[0] != '\0')
		strcat(szRegions,", ");
	strcat(szRegions,"5'US");
	}
if(Regions & cFeatBit5UTR)
	{
	if(szRegions[0] != '\0')
		strcat(szRegions,", ");
	strcat(szRegions,"5'UTR");
	}
if(Regions & cFeatBitCDS)
	{
	if(szRegions[0] != '\0')
		strcat(szRegions,", ");
	strcat(szRegions,"CDS");
	}
if(Regions & cFeatBitIntrons)
	{
	if(szRegions[0] != '\0')
		strcat(szRegions,", ");
	strcat(szRegions,"Introns");
	}
if(Regions & cFeatBit3UTR)
	{
	if(szRegions[0] != '\0')
		strcat(szRegions,", ");
	strcat(szRegions,"3'UTR");
	}
if(Regions & cFeatBitDnstream)
	{
	if(szRegions[0] != '\0')
		strcat(szRegions,", ");
	strcat(szRegions,"3'DS");
	}
if(Regions & cIntronExonSpliceSite)
	{
	if(szRegions[0] != '\0')
		strcat(szRegions,", ");
	strcat(szRegions,"3'Splice");
	}
return(szRegions);
}


// AddExcludeHistory
// Adds a tsExcludeEl to the cached history
// The newly added element will be the MRA
bool
AddExcludeHistory(char *pszSpecies,char *pszChrom,bool bExclude,tsProcParams *pParams)
{
tsExcludeEl *pEl;
if(pParams->gNumExcludeEls < cMaxExcludeHistory)
	pEl = &pParams->gExcludeChroms[pParams->gNumExcludeEls++];
else
	{
	pEl = pParams->gpLRA;		// reuse the least recently accessed element
	pParams->gpLRA = pEl->pPrev;	// 2nd LRA now becomes the LRA
	pParams->gpLRA->pNext = NULL;
	}
if(pParams->gpMRA != NULL)
	pParams->gpMRA->pPrev = pEl;
pEl->pNext = pParams->gpMRA;
pEl->pPrev = NULL;
pParams->gpMRA = pEl;
if(pParams->gpLRA == NULL)
	pParams->gpLRA = pEl;
pEl->bExclude = bExclude;
strcpy(pEl->szSpecies,pszSpecies);
strcpy(pEl->szChrom,pszChrom);
return(bExclude);
}

// LocateExclude
// Locates - starting from the MRA - a tsExcludeEl which matches on szSpecies and szChrom
// If matches then this tsExcludeEl is made the MRA
// Returns ptr to matching tsExcludeEl if located or NULL
tsExcludeEl *
LocateExclude(char *pszSpecies,char *pszChrom,tsProcParams *pParams)
{
tsExcludeEl *pEl;
pEl = pParams->gpMRA;
while(pEl != NULL)
	{
	if(!stricmp(pEl->szSpecies,pszSpecies) && !stricmp(pEl->szChrom,pszChrom))
		{
		if(pParams->gNumExcludeEls ==1 || pParams->gpMRA == pEl)	// if only, or already the MRA then no need for any relinking 
			return(pEl);

		if(pParams->gpLRA == pEl)						// if was the LRA then the 2nd LRA becomes the LRA
			{
			pParams->gpLRA = pEl->pPrev;
			pParams->gpLRA->pNext = NULL;
			}
		else									// not the LRA, and not the MRA
			{
			pEl->pPrev->pNext = pEl->pNext;
			pEl->pNext->pPrev = pEl->pPrev;
			}
		pParams->gpMRA->pPrev = pEl;
		pEl->pNext = pParams->gpMRA;
		pEl->pPrev = NULL;
		pParams->gpMRA = pEl;
		return(pEl);
		}
	pEl = pEl->pNext;
	}
return(NULL);
}


// ExcludeThisChrom
// Returns true if szSpecies.szChromosome is to be excluded from processing
// ExcludeThisChrom
bool
ExcludeThisChrom(char *pszSpecies,char *pszChrom,tsProcParams *pProcParams)
{
#ifdef _WIN32
RegexpMatch mc;
#else
regmatch_t mc;
#endif

char szSpeciesChrom[200];
int Idx;
if(!pProcParams->NumExcludeChroms && !pProcParams->NumIncludeChroms)
	return(false);

// check if this species and chromosome are already known to be included/excluded
tsExcludeEl *pEl;
if((pEl = LocateExclude(pszSpecies,pszChrom,pProcParams))!=NULL)
	return(pEl->bExclude);
// haven't seen this species or chromosome before - or else they have been discarded from history...
sprintf(szSpeciesChrom,"%s.%s",pszSpecies,pszChrom);

// to be explicitly excluded?
for(Idx = 0; Idx < pProcParams->NumExcludeChroms; Idx++)
#ifdef _WIN32	
	if(pProcParams->ExcludeChromsRE[Idx]->Match(szSpeciesChrom,&mc))
#else
	if(!regexec(&pProcParams->ExcludeChromsRE[Idx],szSpeciesChrom,1,&mc,0))
#endif
		return(AddExcludeHistory(pszSpecies,pszChrom,true,pProcParams));

// to be included?
for(Idx = 0; Idx < pProcParams->NumIncludeChroms; Idx++)
	{
#ifdef _WIN32
	if(pProcParams->IncludeChromsRE[Idx]->Match(szSpeciesChrom,&mc))
#else
	if(!regexec(&pProcParams->IncludeChromsRE[Idx],szSpeciesChrom,1,&mc,0))
#endif
		return(AddExcludeHistory(pszSpecies,pszChrom,false,pProcParams));
	}

// if chroms to explicitly include specified, and this chrom not matched then this chrom is to be excluded
return(AddExcludeHistory(pszSpecies,pszChrom,pProcParams->NumIncludeChroms ? true : false,pProcParams));
}


void
Cleanup(tsProcParams *pProcParams)
{
if(pProcParams->pFilterInRefIDs != NULL)
	{
	delete pProcParams->pFilterInRefIDs;
	pProcParams->pFilterInRefIDs = NULL;
	}

if(pProcParams->pFilterOutRefIDs != NULL)
	{
	delete pProcParams->pFilterOutRefIDs;
	pProcParams->pFilterOutRefIDs = NULL;
	}

if(pProcParams->pIncludeFiltLoci != NULL)
	{
	delete pProcParams->pIncludeFiltLoci;
	pProcParams->pIncludeFiltLoci = NULL;
	}

if(pProcParams->pExcludeFiltLoci != NULL)
	{
	delete pProcParams->pExcludeFiltLoci;
	pProcParams->pExcludeFiltLoci = NULL;
	}

if(pProcParams->pCSV != NULL)
	{
	delete pProcParams->pCSV;
	pProcParams->pCSV = NULL;
	}

if(pProcParams->hOutFile != -1)
	{
	close(pProcParams->hOutFile);
	pProcParams->hOutFile = -1;
	}

if(pProcParams->hStatsFile != -1)
	{
	close(pProcParams->hStatsFile);
	pProcParams->hStatsFile = -1;
	}


if(pProcParams->pElLocs != NULL)
	{
#ifdef _WIN32
	free(pProcParams->pElLocs);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(pProcParams->pElLocs != MAP_FAILED)
		munmap(pProcParams->pElLocs,pProcParams->MaxElLocs * sizeof(tsCSVEl));
#endif
	pProcParams->pElChrBuff = NULL;
	}

if(pProcParams->pElChrBuff != NULL)
	{
#ifdef _WIN32
	free(pProcParams->pElChrBuff);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(pProcParams->pElChrBuff != MAP_FAILED)
		munmap(pProcParams->pElChrBuff,pProcParams->MaxAllocdChrs);
#endif
	pProcParams->pElChrBuff = NULL;
	}

if(pProcParams->pElSpecies != NULL)
	{
	delete pProcParams->pElSpecies;
	pProcParams->pElSpecies = NULL;
	}

if(pProcParams->pElChroms != NULL)
	{
	delete pProcParams->pElChroms;
	pProcParams->pElChroms = NULL;
	}

pProcParams->MaxElLocs = 0;
pProcParams->NumElLocs = 0;
pProcParams->MaxAllocdChrs = 0;
pProcParams->NumUsedChrs = 0;
}

int Process(etProcMode ProcMode,		// processing mode
			char *pszInFile,			// file to filter
			char *pszOutFile,			// write filtered file
			char *pszStatsFile,			// write stats file
			int RegionsIn,					// regions to retain
			int RegionsOut,				// regions to remove
			int SelectN,				// output at most this many rows
			char *pszOutputFields,		// only output these fields
			char *pszSpeciesList,		// list of species to filter in ('\0' if no filtering)
			bool bOverlaps,				// filter out overlaps
			bool bNoOverlaps,			// filter out none overlaps
			char *pszFilterOutRefIDFile, // RefIDs to filter out
			char *pszFilterInRefIDFile,	 // RefIDs to filter in unless filtered out
			int MinElLen,				 // elements must be of at least this length
			int MaxElLen,				// no longer than this length
			int Align2Core,				// at least this many bases aligned in outspecies for it to count as aligned to core
			double PCAlign2Core,		// minimum percentage (0..100) aligned to core to count as as aligned to core (default is 0%)");
			double IDIdent2Core,		// minimum identity (0..100) aligned to core to count as as aligned to core (default is 0%)");
			double IDOSIdentity,		// minimum outspecies (matches/matchesPlusmismatches) (0..100) aligned to core to count as as aligned to core (default is 1%)");
			int NumIncludeChroms,		// number of chromosome regular expressions to include
			char **ppszIncludeChroms,	// array of include chromosome regular expressions
			int NumExcludeChroms,		// number of chromosome expressions to exclude
			char **ppszExcludeChroms,	// array of exclude chromosome regular expressions
			int NumIncludeLoci,			// number of include region files 
			char **ppszIncludeLoci,		// array of include region files
			int NumExcludeLoci,			// number of exclude region files
			char **ppszExcludeLoci)		// array of exclude region files
{
int Idx;
int Rslt;
tsProcParams ProcParams;

#ifndef _WIN32
int RegErr;				// regular expression parsing error
char szRegErr[128];			// to hold RegErr as textual representation ==> regerror();
#endif

memset(&ProcParams,0,sizeof(ProcParams));
ProcParams.ProcMode = ProcMode;
ProcParams.hOutFile = -1;
ProcParams.hStatsFile = -1;
ProcParams.RegionsIn = RegionsIn;
ProcParams.RegionsOut= RegionsOut;
ProcParams.SelectN = SelectN;
ProcParams.MaxElLen = MaxElLen;
ProcParams.MinElLen = MinElLen;
ProcParams.pszOutFile = pszOutFile;
ProcParams.pszStatsFile = pszStatsFile;
ProcParams.pszInFile = pszInFile;
ProcParams.pszSpeciesList = pszSpeciesList;
ProcParams.bNoOverlaps = bNoOverlaps;
ProcParams.bOverlaps = bOverlaps;
ProcParams.bFilterOutputFields = false;
ProcParams.Align2Core = Align2Core;
ProcParams.PCAlign2Core = PCAlign2Core;
ProcParams.IDIdent2Core = IDIdent2Core;
ProcParams.IDOSIdentity = IDOSIdentity;

if(pszOutputFields != NULL && pszOutputFields[0] != '\0')
	if(ParseNumFields(pszOutputFields,&ProcParams) < 1)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to parse output field identifiers");
		return(-1);
		}

if(NumIncludeLoci)
	{
	if((ProcParams.pIncludeFiltLoci = new CFilterLoci) == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instance CFilterLoci");
		return(-1);
		}
	for(Idx=0;Idx<NumIncludeLoci; Idx++)
		{
		if((Rslt = ProcParams.pIncludeFiltLoci->Load(ppszIncludeLoci[Idx])) < 0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open/load %s",ppszIncludeLoci[Idx]);
			Cleanup(&ProcParams);
			return(Rslt);
			}
		}
	}
if(NumExcludeLoci)
	{
	if((ProcParams.pExcludeFiltLoci = new CFilterLoci) == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instance CFilterLoci");
		Cleanup(&ProcParams);
		return(-1);
		}
	for(Idx=0;Idx<NumExcludeLoci; Idx++)
		{
		if((Rslt = ProcParams.pExcludeFiltLoci->Load(ppszExcludeLoci[Idx]))<0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open/load %s",ppszExcludeLoci[Idx]);
			Cleanup(&ProcParams);
			return(Rslt);
			}
		}
	}

if(pszFilterInRefIDFile != NULL && pszFilterInRefIDFile[0])
	{
	ProcParams.pFilterInRefIDs = new CFilterRefIDs;
	if((Rslt=ProcParams.pFilterInRefIDs->Open(pszFilterInRefIDFile)) < 0)
		{
		while(ProcParams.pFilterInRefIDs->NumErrMsgs())
			gDiagnostics.DiagOut(eDLFatal,gszProcName,ProcParams.pFilterInRefIDs->GetErrMsg());
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open/parse input filter RefID file %s",pszFilterInRefIDFile);
		Cleanup(&ProcParams);
		return(Rslt);
		}
	}

if(pszFilterOutRefIDFile != NULL && pszFilterOutRefIDFile[0])
	{
	ProcParams.pFilterOutRefIDs = new CFilterRefIDs;
	if((Rslt=ProcParams.pFilterOutRefIDs->Open(pszFilterOutRefIDFile)) < 0)
		{
		while(ProcParams.pFilterOutRefIDs->NumErrMsgs())
			gDiagnostics.DiagOut(eDLFatal,gszProcName,ProcParams.pFilterOutRefIDs->GetErrMsg());
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open/parse input filter RefID file %s",pszFilterOutRefIDFile);
		Cleanup(&ProcParams);
		return(Rslt);
		}
	}

ProcParams.NumIncludeChroms = NumIncludeChroms;		// number of chromosomes explicitly defined to be included
ProcParams.ppszIncludeChroms = ppszIncludeChroms;	// ptr to array of reg expressions defining chroms to include - overides exclude
ProcParams.NumExcludeChroms = NumExcludeChroms;		// number of chromosomes explicitly defined to be excluded
ProcParams.ppszExcludeChroms = ppszExcludeChroms;	// ptr to array of reg expressions defining chroms to include

// parse out species list
if(pszSpeciesList == NULL || pszSpeciesList[0] == '\0')
	ProcParams.NumSpecies = 0;
else
	ProcParams.NumSpecies = ParseNumSpecies(pszSpeciesList,&ProcParams); 

#ifdef _WIN32
try {
	for(Idx=0;Idx < NumIncludeChroms;Idx++)
		{
		ProcParams.IncludeChromsRE[Idx] = new Regexp();
		ProcParams.IncludeChromsRE[Idx]->Parse(ppszIncludeChroms[Idx],false);	// note case insensitive
		}
	}
catch(...)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to process include regexpr chrom '%s'",ppszIncludeChroms[Idx]);
	Cleanup(&ProcParams);
	return(eBSFerrMem);
	}
try {
	for(Idx=0;Idx < NumExcludeChroms;Idx++)
		{
		ProcParams.ExcludeChromsRE[Idx] = new Regexp();
		ProcParams.ExcludeChromsRE[Idx]->Parse(ppszExcludeChroms[Idx],false);	// note case insensitive
		}
	}
catch(...)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to process exclude regexpr chrom '%s'",ppszExcludeChroms[Idx]);
	Cleanup(&ProcParams);
	return(eBSFerrMem);
	}

#else
for(Idx=0;Idx < NumIncludeChroms;Idx++)
	{

	RegErr=regcomp(&ProcParams.IncludeChromsRE[Idx],ppszIncludeChroms[Idx],REG_EXTENDED | REG_ICASE);	// note case insensitive
	if(RegErr)
		{
		regerror(RegErr,&ProcParams.IncludeChromsRE[Idx],szRegErr,sizeof(szRegErr));
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to process include chrom '%s' error: %s",ppszIncludeChroms[Idx],szRegErr);
		Cleanup(&ProcParams);
		return(eBSFerrMem);
		}
	}
for(Idx=0;Idx < NumExcludeChroms;Idx++)
	{
	RegErr = regcomp(&ProcParams.ExcludeChromsRE[Idx],ppszExcludeChroms[Idx],REG_EXTENDED | REG_ICASE);	// note case insensitive
	if(RegErr)
		{
		regerror(RegErr,&ProcParams.ExcludeChromsRE[Idx],szRegErr,sizeof(szRegErr));
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to process exclude chrom '%s' error: %s",ppszExcludeChroms[Idx],szRegErr);
		Cleanup(&ProcParams);
		return(eBSFerrMem);
		}
	}
#endif

ProcParams.pCSV = new CCSVFile;
if(ProcParams.pCSV == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CCSVfile");
	Cleanup(&ProcParams);
	return(eBSFerrObj);
	}
ProcParams.pCSV->SetMaxFields(cMaxNumFields);	
if((Rslt=ProcParams.pCSV->Open(pszInFile))!=eBSFSuccess)
	{
	while(ProcParams.pCSV->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,ProcParams.pCSV->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open file: %s",pszInFile);
	Cleanup(&ProcParams);
	return(Rslt);
	}

if(pszOutFile != NULL && pszOutFile[0] != '\0')
	{
#ifdef _WIN32
if((ProcParams.hOutFile = open(pszOutFile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
if((ProcParams.hOutFile = open(pszOutFile, O_RDWR | O_CREAT | O_TRUNC, S_IREAD | S_IWRITE ))==-1)
#endif
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create output filtered file: %s - %s",pszOutFile,strerror(errno));
		Cleanup(&ProcParams);
		return(eBSFerrCreateFile);
		}
	}

if(pszStatsFile != NULL && pszStatsFile[0] != '\0')
	{
	char szLineBuff[2048];
	int LineLen;


#ifdef _WIN32
	if((ProcParams.hStatsFile = open(pszStatsFile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
	if((ProcParams.hStatsFile = open(pszStatsFile, O_RDWR | O_CREAT | O_TRUNC, S_IREAD | S_IWRITE ))==-1)
#endif
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create output stats file: %s - %s",pszStatsFile,strerror(errno));
		Cleanup(&ProcParams);
		return(eBSFerrCreateFile);
		}

	LineLen = sprintf(szLineBuff,"\"RefID\",\"Filtered\",\"OverLen\",\"UnderLen\",\"ExcludeLoci\",\"IncludeLoci\",\"Overlaps\",\"Chrom\",\"OutRefID\",\"InRefID\",\"Species\"\n");
	CUtility::SafeWrite(ProcParams.hStatsFile,szLineBuff,LineLen);	
	}

Rslt = FilterCSV(&ProcParams);

Cleanup(&ProcParams);
return(Rslt);
}

// GenNameHash
// Generates a 16bit hash on specified lowercased name
UINT16 
GenNameHash(char *pszName)
{
unsigned long hash = 5381;
char Chr;
while (Chr = *pszName++)
	hash = ((hash << 5) + hash) + tolower(Chr);
return ((UINT16)hash);
}


tsElChrom *
LocateElChrom(char *pszChrom,
				tsProcParams *pProcParams) // global processing parameters
{
tsElChrom *pChrom;
UINT16 Hash;
int Idx;
if(pProcParams->pElChroms == NULL || pProcParams->NumElChroms < 1)
	return(NULL);
Hash = GenNameHash(pszChrom);

// quick check against cached chrom name
if(pProcParams->CachElChromID > 0)
	{
	pChrom = &pProcParams->pElChroms[pProcParams->CachElChromID-1];
	if(Hash == pChrom->Hash && !stricmp(pChrom->szChrom,pszChrom))
		return(pChrom);
	}
pChrom = pProcParams->pElChroms;
for(Idx = 0; Idx < pProcParams->NumElChroms; Idx++,pChrom++)
	{
	if(Hash == pChrom->Hash && !stricmp(pChrom->szChrom,pszChrom))
		{
		pProcParams->CachElChromID = pChrom->ChromID;
		return(pChrom);
		}
	}
return(NULL);
}



int
AddChrom(char *pszChrom,tsProcParams *pParams)
{
tsElChrom *pChrom;
int Idx;
UINT16 Hash;
Hash = GenNameHash(pszChrom);

if(pParams->pElChroms != NULL && pParams->NumElChroms)
	{
	if(pParams->CachElChromID > 0)
		{
		pChrom = &pParams->pElChroms[pParams->CachElChromID-1];
		if((pChrom->Hash == Hash) && !stricmp(pszChrom,pChrom->szChrom))
			return(pChrom->ChromID);
		}

	pChrom = pParams->pElChroms;
	for(Idx = 0; Idx < pParams->NumElChroms; Idx++,pChrom++)
		{
		if((pChrom->Hash == Hash) && !stricmp(pszChrom,pChrom->szChrom))
			return(pParams->CachElChromID = pChrom->ChromID);
		}
	}

if(pParams->pElChroms == NULL || pParams->NumElChroms == pParams->MaxElChroms)
	{
	if((pChrom = new tsElChrom[pParams->MaxElChroms + cElChromAllocChunk])==NULL)
		{
		return(eBSFerrMem);
		}
	if(pParams->pElChroms != NULL)
		{
		memcpy(pChrom,pParams->pElChroms,sizeof(tsElChrom) * pParams->NumElChroms);
		pParams->pElChroms = pChrom;
		pParams->MaxElChroms += cElChromAllocChunk;
		}
	else
		{
		pParams->MaxElChroms = cElChromAllocChunk;
		pParams->NumElLocs = 0;
		}
	pParams->pElChroms = pChrom;
	}
pChrom = &pParams->pElChroms[pParams->NumElChroms++];
pChrom->ChromID = pParams->NumElChroms;
pChrom->Hash = Hash;
strcpy(pChrom->szChrom,pszChrom);
return(pParams->CachElChromID = pChrom->ChromID);
}

int
AddSpecies(char *pszSpecies,tsProcParams *pParams)
{
tsElSpecies *pSpecies;
int Idx;
UINT16 Hash;
Hash = GenNameHash(pszSpecies);

if(pParams->pElSpecies != NULL && pParams->NumElSpecies)
	{
	if(pParams->CachElSpeciesID > 0)
		{
		pSpecies = &pParams->pElSpecies[pParams->CachElSpeciesID-1];
		if((pSpecies->Hash == Hash) && !stricmp(pszSpecies,pSpecies->szSpecies))
			return(pSpecies->SpeciesID);
		}

	pSpecies = pParams->pElSpecies;
	for(Idx = 0; Idx < pParams->NumElSpecies; Idx++,pSpecies++)
		{
		if((pSpecies->Hash == Hash) && !stricmp(pszSpecies,pSpecies->szSpecies))
			return(pParams->CachElSpeciesID = pSpecies->SpeciesID);
		}
	}

if(pParams->pElSpecies == NULL || pParams->NumElSpecies == pParams->MaxElSpecies)
	{
	if((pSpecies = new tsElSpecies[pParams->MaxElSpecies + cElSpeciesAllocChunk])==NULL)
		{
		return(eBSFerrMem);
		}
	if(pParams->pElSpecies != NULL)
		{
		memcpy(pSpecies,pParams->pElSpecies,sizeof(tsElSpecies) * pParams->NumElSpecies);
		pParams->pElSpecies = pSpecies;
		pParams->MaxElSpecies += cElSpeciesAllocChunk;
		}
	else
		{
		pParams->MaxElChroms = cElSpeciesAllocChunk;
		pParams->NumElLocs = 0;
		}
	pParams->pElSpecies = pSpecies;
	}
pSpecies = &pParams->pElSpecies[pParams->NumElSpecies++];
pSpecies->SpeciesID = pParams->NumElSpecies;
pSpecies->Hash = Hash;
strcpy(pSpecies->szSpecies,pszSpecies);
return(pParams->CachElSpeciesID = pSpecies->SpeciesID);
}



tsCSVEl *
AddCSVEl(tsFiltState FiltState,	// current filter state for this element
		 int RefID,				// reference identifier from input file
		 char *pszSpecies,		// species
		 char *pszChrom,		// chrom
		 int StartLoci,			// start loci
		 int EndLoci,			// end loci
		 int LineLen,			// num of chars (excl trailing '\0') in pszLineBuff
		 char *pszLine,			// line to output if not filtered
		 tsProcParams *pParams)
{
char *pTmp;
size_t memreq;
tsCSVEl *pEl;

int ChromID = AddChrom(pszChrom,pParams);
int SpeciesID = AddSpecies(pszSpecies,pParams);

if(LineLen && ((LineLen + pParams->NumUsedChrs + 1) >= pParams->MaxAllocdChrs))
	{
	memreq = pParams->MaxAllocdChrs + cElChrAllocChunk;

#ifdef _WIN32
	pTmp = (char *) realloc(pParams->pElChrBuff,memreq);
#else
	pTmp = (char *)mremap(pParams->pElChrBuff,pParams->MaxAllocdChrs,memreq,MREMAP_MAYMOVE);
	if(pTmp == MAP_FAILED)
		pTmp = NULL;
#endif
	if(pTmp == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddCSVEl: Memory re-allocation to %d bytes - %s",memreq,strerror(errno));
		return(NULL);
		}
	pParams->pElChrBuff = pTmp;
	pParams->MaxAllocdChrs = memreq;
	}

if(pParams->pElLocs == NULL || pParams->NumElLocs >= pParams->MaxElLocs)
	{
	memreq = sizeof(tsCSVEl) * (pParams->MaxElLocs + cElLociAllocChunk);

#ifdef _WIN32
	pEl = (tsCSVEl *) realloc(pParams->pElLocs,memreq);
#else
	pEl = (tsCSVEl *)mremap(pParams->pElLocs,sizeof(tsCSVEl) * pParams->MaxElLocs,memreq,MREMAP_MAYMOVE);
	if(pEl == MAP_FAILED)
		pEl = NULL;
#endif
	if(pEl == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddCSVEl: Memory re-allocation to %d bytes - %s",memreq,strerror(errno));
		return(NULL);
		}
	pParams->pElLocs = pEl;
	pParams->MaxElLocs += cElLociAllocChunk;
	}

pEl = &pParams->pElLocs[pParams->NumElLocs++];
pEl->RefID = RefID;
pEl->Filter = FiltState;
pEl->ChromID = ChromID;
pEl->SpeciesID = SpeciesID;
pEl->StartLoci = StartLoci;
pEl->EndLoci = EndLoci;
pEl->LineLen = LineLen;
if(LineLen)
	{
	strcpy(&pParams->pElChrBuff[pParams->NumUsedChrs],pszLine);
	pEl->LineOfs = pParams->NumUsedChrs;
	pParams->NumUsedChrs += LineLen + 1;
	}
else
	pEl->LineOfs = -1;
	
return(pEl);
}

int 
FilterCSV(tsProcParams *pParams)
{
int Rslt = 0;
int NumFields;
int RefID;
char *pszType;
char *pszSpecies;
char *pszChrom;
int StartLoci;
int EndLoci;
int Region;
char *pszFieldVal;
int Len;
int FieldIdx;
int BuffLen;
char szLineBuff[0x07fff];

int NumElsRead;		// number of elements before filtering
int NumElsAccepted;	// number of elements accepted after filtering
tsFiltState FiltState; // filter state for current element

tsCSVEl *pEl;
tsCSVEl *pEl1;
size_t Idx;
size_t Idy;
char *pszLine;
int ChkOverlaps;

// prealloc memory 
#ifdef _WIN32
	pParams->pElChrBuff = (char *) malloc(cElChrAllocChunk);	// initial and perhaps the only allocation

	if(pParams->pElChrBuff == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"FilterCSV: Memory allocation of %lld bytes - %s",(INT64)cElChrAllocChunk,strerror(errno));
		Cleanup(pParams);
		return(eBSFerrMem);
		}
#else
// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
	pParams->pElChrBuff = (char *)mmap(NULL,cElChrAllocChunk, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
	if(pParams->pElChrBuff == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"FilterCSV: Memory allocation of %lld bytes through mmap()  failed - %s",(INT64)cElChrAllocChunk,strerror(errno));
		pParams->pElChrBuff = NULL;
		Cleanup(pParams);
		return(eBSFerrMem);
		}
#endif
	pParams->MaxAllocdChrs = cElChrAllocChunk;
	pParams->NumUsedChrs = 0;

#ifdef _WIN32
	pParams->pElLocs = (tsCSVEl *) malloc(cElLociAllocChunk * sizeof(tsCSVEl));	// initial and perhaps the only allocation

	if(pParams->pElLocs == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"FilterCSV: Memory allocation of %lld bytes - %s",(INT64)cElLociAllocChunk * sizeof(tsCSVEl),strerror(errno));
		Cleanup(pParams);
		return(eBSFerrMem);
		}
#else
// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
	pParams->pElLocs = (tsCSVEl *)mmap(NULL,cElLociAllocChunk * sizeof(tsCSVEl), PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
	if(pParams->pElLocs == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"FilterCSV: Memory allocation of %lld bytes through mmap()  failed - %s",(INT64)cElLociAllocChunk * sizeof(tsCSVEl),strerror(errno));
		pParams->pElLocs = NULL;
		Cleanup(pParams);
		return(eBSFerrMem);
		}
	
#endif
	pParams->MaxElLocs = cElLociAllocChunk;
	pParams->NumElLocs = 0;

NumElsRead =0;		// number of elements before filtering
NumElsAccepted =0;	// number of elements accepted after filtering
printf("\n           Parsing CSV line: 0\r");
while((Rslt=pParams->pCSV->NextLine()) > 0)	// onto next line containing fields
	{
	if(!(NumElsRead % 100000))
		printf("           Parsing CSV line: %d\r",NumElsRead);

	NumFields = pParams->pCSV->GetCurFields();

	if(pParams->ProcMode == eProcStandard)
		{
		if(NumFields < 7)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expected at least 7 fields in '%s', GetCurFields() returned '%d'",pParams->pszInFile,NumFields);
			return(eBSFerrFieldCnt);
			}
		}
	else
		{
		if(NumFields < 13)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expected at least 13 fields in '%s', GetCurFields() returned '%d'",pParams->pszInFile,NumFields);
			return(eBSFerrFieldCnt);
			}
		}


	NumElsRead += 1;
	memset(&FiltState,0,sizeof(FiltState));
	pParams->pCSV->GetInt(1,&RefID);
	if(pParams->pFilterOutRefIDs != NULL && pParams->pFilterOutRefIDs->Locate(RefID))
		{
		FiltState.fOutRefID = 1;
		FiltState.bFiltOut = true;
		}

	if(pParams->pFilterInRefIDs != NULL && !pParams->pFilterInRefIDs->Locate(RefID))
		{	
		FiltState.fInRefID = 1;
		FiltState.bFiltOut = true;
		}

	pParams->pCSV->GetInt(7,&Len);
	if(Len < pParams->MinElLen)
		{
		FiltState.fUnderLen = 1;
		FiltState.bFiltOut = true;
		}
	else
		if(Len > pParams->MaxElLen)
			{
			FiltState.fOverLen = 1;
			FiltState.bFiltOut = true;
			}

	pParams->pCSV->GetText(2,&pszType);
	pParams->pCSV->GetText(3,&pszSpecies);
	if(pParams->NumSpecies)
		{
		for(Idx=0; Idx < (size_t)pParams->NumSpecies; Idx++)
			{
			if(!stricmp(pszSpecies,pParams->szSpecies[Idx]))
				break;
			}
		if(Idx == pParams->NumSpecies)
			{
			FiltState.fSpecies = 1;
			FiltState.bFiltOut = true;
			}
		}

	pParams->pCSV->GetText(4,&pszChrom);
	if(ExcludeThisChrom(pszSpecies,pszChrom,pParams))
		{
		FiltState.fChrom = 1;
		FiltState.bFiltOut = true;
		}

	pParams->pCSV->GetInt(5,&StartLoci);
	pParams->pCSV->GetInt(6,&EndLoci);

	if(pParams->pExcludeFiltLoci != NULL && pParams->pExcludeFiltLoci->Locate(pszChrom,StartLoci,EndLoci))
		{	
		FiltState.fExcludeLoci = 1;
		FiltState.bFiltOut = true;
		}

	if(pParams->pIncludeFiltLoci != NULL && !pParams->pIncludeFiltLoci->Locate(pszChrom,StartLoci,EndLoci))
		{	
		FiltState.fIncludeLoci = 1;
		FiltState.bFiltOut = true;
		}

	if(pParams->RegionsIn || pParams->RegionsOut)
		{
		if(NumFields < 9)				// need at least 9 fields to filter on regions
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expected at least 9 fields for region filtering in '%s', GetCurFields() returned '%d'",pParams->pszInFile,NumFields);
			return(eBSFerrFieldCnt);
			}
		pParams->pCSV->GetInt(9,&Region);
		if(pParams->RegionsOut)
			{
			if((Region == 0 && (pParams->RegionsOut & cFeatBitIG)) || Region & pParams->RegionsOut)	
				{
				FiltState.fExcludeRegion = 1;
				FiltState.bFiltOut = true;
				}
			}
		if(pParams->RegionsIn)
			{
			if(Region == 0 && !(pParams->RegionsIn & cFeatBitIG))	
				{
				FiltState.fIncludeRegion = 1;
				FiltState.bFiltOut = true;
				}
			else
				if(Region != 0)
					{
					int Msk = 0x01;
					for(Idx = 0; Idx < 6; Idx++,Msk <<= 1)
						{
						if((pParams->RegionsIn & Msk) && (Region & 0x03f) == Msk)
							break;
						}
					if(Idx == 6)
						{
						FiltState.fIncludeRegion = 1;
						FiltState.bFiltOut = true;
						}
					}
			}
		}

	if(pParams->ProcMode == eProcOutspecies)
		{
		int Matching;
		int Mismatching;
		pParams->pCSV->GetInt(11,&Matching);
		pParams->pCSV->GetInt(12,&Mismatching);
		
		if(pParams->Align2Core)				// filter on number of bases aligned to core?
			{
			if((Matching + Mismatching) < pParams->Align2Core)
				{
				FiltState.fDeltaIdentity = 1;
				FiltState.bFiltOut = true;
				}
			}

		if(pParams->PCAlign2Core > 0.0)			// filter on percentage aligned
			{
			if(((100.0 * (Matching + Mismatching))/Len) < pParams->PCAlign2Core)
				{
				FiltState.fDeltaIdentity = 1;
				FiltState.bFiltOut = true;
				}
			}
		if(pParams->IDIdent2Core > 0.0)			// filter on matches/corelen
			{
			if(((100.0 * Matching)/Len) < pParams->IDIdent2Core)
				{
				FiltState.fDeltaIdentity = 1;
				FiltState.bFiltOut = true;
				}
			}
		if(pParams->IDOSIdentity > 0.0)			// filter on outspecies identity
			{
			if(!(Matching + Mismatching) || ((100.0 * Matching)/(Matching + Mismatching)) < pParams->IDOSIdentity)
				{
				FiltState.fDeltaIdentity = 1;
				FiltState.bFiltOut = true;
				}
			}
		}

	BuffLen = 0;
	szLineBuff[0] = '\0';
	if(!FiltState.bFiltOut)	// if this element is thus far unfiltered then build output text
		{
		if(pParams->hOutFile != -1)
			{
			BuffLen = 0;
			for(FieldIdx = 1; FieldIdx <= NumFields; FieldIdx++)
				{
				if(pParams->bFilterOutputFields && !pParams->OutputFields[FieldIdx-1])
					continue;
				pParams->pCSV->GetText(FieldIdx,&pszFieldVal);
				if(BuffLen)
					BuffLen += sprintf(&szLineBuff[BuffLen],",");
				if(pParams->pCSV->GetQuoted(FieldIdx))
					BuffLen += sprintf(&szLineBuff[BuffLen],"\"%s\"",pszFieldVal);
				else
					BuffLen += sprintf(&szLineBuff[BuffLen],"%s",pszFieldVal);
				}
			BuffLen += sprintf(&szLineBuff[BuffLen],"\n");
			}
		NumElsAccepted += 1;
		}

	if(AddCSVEl(FiltState,RefID,pszSpecies,pszChrom,StartLoci,EndLoci,BuffLen,szLineBuff,pParams)==NULL)
		return(eBSFerrMem);
	}
printf("           Parsing CSV line: %d\n           Parsing CSV Completed\n",NumElsRead);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Parsing %d lines from CSV Completed of which %d elements were accepted for further processing",NumElsRead,NumElsAccepted);

int NumAccepted=0;
int NumRejected = 0;
int OverLens = 0;
int UnderLens = 0;
int Identities = 0;
int ExcludeLoci = 0;
int IncludeLoci = 0;
int IncludeRegions = 0;
int ExcludeRegions = 0;
int Overlaps= 0;
int Chroms = 0;
int Loci= 0;
int OutRefIDs = 0;
int InRefIDs = 0;
int Species = 0; 
int SelectNs = 0;

// now to check on overlaps
ChkOverlaps = 0;
if(pParams->NumElLocs > 0)
	{
	if(pParams->bOverlaps || pParams->bNoOverlaps)
		{
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing for overlaps started");
		qsort(pParams->pElLocs,pParams->NumElLocs,sizeof(tsCSVEl),CompareCSVEls);

		pEl = pParams->pElLocs;
		for(Idx = 0; Idx < (pParams->NumElLocs - 1); Idx++, pEl++)
			{
			pEl1 = pEl+1;
			for(Idy = Idx + 1; Idy < pParams->NumElLocs; Idy++, pEl1++)
				{
				if(pEl->SpeciesID != pEl1->SpeciesID ||
					pEl->ChromID != pEl1->ChromID ||
					pEl->EndLoci < pEl1->StartLoci)
					break;

				pEl->Filter.fOverlap = 1;
				pEl1->Filter.fOverlap = 1;
				}
			if(pEl->Filter.fOverlap)
				ChkOverlaps += 1;
			}
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing for overlaps completed");
		}
	}

pEl = pParams->pElLocs;
NumAccepted=0;
for(Idx = 0; Idx < pParams->NumElLocs; Idx++,pEl++)
	{
	if(pEl->Filter.fOverlap && pParams->bOverlaps ||
		!pEl->Filter.fOverlap && pParams->bNoOverlaps)
		pEl->Filter.bFiltOut = true;
	if(!pEl->Filter.bFiltOut)
		NumAccepted++;
	}

#ifdef USEOLDNONRANDOM
int NumSelected;
if(NumAccepted > pParams->SelectN)		// need to reduce number of rows output?
	{
	NumSelected = 0;
	pEl = pParams->pElLocs;
	int Num2Select = pParams->SelectN;
	SelectNs = 0;
	for(Idx = 0; Idx < pParams->NumElLocs; Idx++,pEl++)
		{
		if(pEl->Filter.bFiltOut)
			continue;
		if(!SelectNs)
			{
			SelectNs = NumAccepted  / Num2Select;
			NumAccepted -= 1;
			NumSelected += 1;
			Num2Select -= 1;
			continue;
			}
		pEl->Filter.fSelectN = true;
		pEl->Filter.bFiltOut = true;
		SelectNs -= 1;
		NumAccepted -= 1;
		}
	}
else
	NumSelected = NumAccepted;
#endif

if(NumAccepted > pParams->SelectN)		// need to reduce number of rows output?
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Limiting number of accepted elements by randomly selecting elements to report");
	TRandomCombined<TRanrotWGenerator,TRandomMersenne> RG((int)time(0));
	int ToRemove;
	while(NumAccepted > pParams->SelectN)
		{
		ToRemove = RG.IRandom(0,(long)pParams->NumElLocs);
		pEl = &pParams->pElLocs[ToRemove];
		for(Idx = ToRemove; Idx < pParams->NumElLocs; Idx++,pEl++)
			{
			if(!pEl->Filter.bFiltOut)
				{
				pEl->Filter.fSelectN = true;
				pEl->Filter.bFiltOut = true;
				NumAccepted -= 1;
				break;
				}
			}

		if(Idx != pParams->NumElLocs)
			continue;

		pEl = pParams->pElLocs;
		for(Idx = 0; Idx < (size_t)ToRemove; Idx++,pEl++)
			{
			if(!pEl->Filter.bFiltOut)
				{
				pEl->Filter.fSelectN = true;
				pEl->Filter.bFiltOut = true;
				NumAccepted -= 1;
				break;
				}
			}
		}
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Limiting number of accepted elements processing completed");
	}

if(pParams->hOutFile != -1)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Write accepted elements to output file '%s'",pParams->pszOutFile);
	pEl = pParams->pElLocs;
	for(Idx = 0; Idx < pParams->NumElLocs; Idx++,pEl++)
		{
		if(pEl->Filter.bFiltOut)
			continue;
		if(pParams->hOutFile != -1)
			{
			pszLine = &pParams->pElChrBuff[pEl->LineOfs];
			CUtility::SafeWrite(pParams->hOutFile,pszLine,pEl->LineLen);
			}
		}
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Write accepted elements to output file completed");
	}


// now for the stats part!
if(pParams->hStatsFile != -1)
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Write distribution stats to output file '%s'",pParams->pszStatsFile);
NumAccepted=0;
NumRejected = 0;
OverLens = 0;
UnderLens = 0;
Identities = 0;
ExcludeLoci = 0;
IncludeLoci = 0;
IncludeRegions = 0;
ExcludeRegions = 0;
Overlaps= 0;
Chroms = 0;
Loci= 0;
OutRefIDs = 0;
InRefIDs = 0;
Species = 0; 
SelectNs = 0;
pEl = pParams->pElLocs;
for(Idx = 0; Idx < pParams->NumElLocs; Idx++,pEl++)
	{
	if(pParams->hStatsFile != -1)
		{
		BuffLen = sprintf(szLineBuff,"%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d\n",
			pEl->RefID,pEl->Filter.bFiltOut ? 1 : 0,pEl->Filter.fOverLen,pEl->Filter.fUnderLen,
			pEl->Filter.fDeltaIdentity,
			pEl->Filter.fExcludeLoci,pEl->Filter.fIncludeLoci,
			(pEl->Filter.fOverlap && pParams->bOverlaps || !pEl->Filter.fOverlap && pParams->bNoOverlaps) ? 1 : 0,
				pEl->Filter.fChrom,pEl->Filter.fOutRefID,pEl->Filter.fInRefID,pEl->Filter.fSpecies,pEl->Filter.fSelectN,pEl->Filter.fIncludeRegion,pEl->Filter.fExcludeRegion);
		CUtility::SafeWrite(pParams->hStatsFile,szLineBuff,BuffLen);
		}

	if(!pEl->Filter.bFiltOut)
		{
		NumAccepted += 1;
		continue;
		}

	if(pEl->Filter.fOverLen)
		OverLens+=1;
	if(pEl->Filter.fUnderLen)
		UnderLens+=1;
	if(pEl->Filter.fExcludeLoci)
		ExcludeLoci+=1;
	if(pEl->Filter.fIncludeLoci)
		IncludeLoci+=1;
	if(pEl->Filter.fOverlap && pParams->bOverlaps ||
		!pEl->Filter.fOverlap && pParams->bNoOverlaps)
		Overlaps+=1;
	if(pEl->Filter.fChrom)
		Chroms+=1;
	if(pEl->Filter.fOutRefID)
		OutRefIDs+=1;
	if(pEl->Filter.fInRefID)
		InRefIDs+=1;
	if(pEl->Filter.fSpecies)
		Species+=1;
	if(pEl->Filter.fDeltaIdentity)
		Identities += 1;
	if(pEl->Filter.fSelectN)
		SelectNs += 1;

	if(pEl->Filter.fIncludeRegion)
		IncludeRegions += 1;
	if(pEl->Filter.fExcludeRegion)
		ExcludeRegions += 1;

	NumRejected += 1;
	}

if(pParams->hStatsFile != -1)
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Write distribution stats to output file completed");

BuffLen = sprintf(szLineBuff,"Processed: %d, Accepted: %d, Rejected: %d, Overlen: %d, UnderLen: %d, Identity: %d, ExcludeLoci: %d, Include Loci: %d, Overlaps: %d, Chroms: %d, OutRefIDs: %d, InRefIDs: %d, Species: %d, SelectNs: %d, IncludeRegions: %d, ExcludeRegions: %d",
				  (int)pParams->NumElLocs,NumAccepted,NumRejected,OverLens,UnderLens,Identities,ExcludeLoci,IncludeLoci,Overlaps,Chroms,OutRefIDs,InRefIDs,Species,SelectNs,IncludeRegions,ExcludeRegions); 

gDiagnostics.DiagOut(eDLInfo,gszProcName,szLineBuff);
return(Rslt);
}


// CompareCSVEls
// Used to sort CSV input loci
int 
CompareCSVEls( const void *arg1, const void *arg2)
{
tsCSVEl *pEl1 = (tsCSVEl *)arg1;
tsCSVEl *pEl2 = (tsCSVEl *)arg2;

if(pEl1->SpeciesID < pEl2->SpeciesID)
	return(-1);
if(pEl1->SpeciesID > pEl2->SpeciesID)
	return(1);
if(pEl1->ChromID == pEl2->ChromID)
	{
	if(pEl1->StartLoci < pEl2->StartLoci)
		return(-1);
	if(pEl1->StartLoci > pEl2->StartLoci)
		return(1);
	if(pEl1->EndLoci < pEl2->EndLoci) 
		return(-1);									
	if(pEl1->EndLoci > pEl2->EndLoci)
		return(1);
	return(0);
	}
return(pEl1->ChromID < pEl2->ChromID ? -1 : 1);
}

