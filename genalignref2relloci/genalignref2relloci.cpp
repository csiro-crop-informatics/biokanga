// GenAlignRef2RelLoci.cpp : Defines the entry point for the console application.
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

const unsigned int cProgVer = 301;		// increment with each release

const int cCoreLociAllocChunk = 10000;	// grow core loci array by this many elements
const int cCoreChromAllocChunk = 1000;	// grow core chrom array by this many elements
const int cMaxOverlapIDs = 2000;		// can handle this many overlaps before truncating list

typedef enum eProcMode {
	eProcModeStandard = 0,		// default processing
	eProcModeFilterRef			// filter loci by reference species
} etProcMode;

typedef struct TAG_sCoreLoci {
	int RefID;
	int ChromID;
	int StartLoci;
	int EndLoci;
} tsCoreLoci;

typedef struct TAG_sCoreChrom {
	int ChromID;							// unique identifier for this chromosome
	tsCoreLoci *pCoreLoci;					// pts to first core on this chromosome
	UINT16 Hash;								// hash over chromosome name
	char szChrom[cMaxDatasetSpeciesChrom];			// core chromosome
} tsCoreChrom;

typedef struct TAG_sSpeciesAlign {
	int SpeciesID;				// species identifier
	char szSpeciesName[cMaxDatasetSpeciesChrom]; // species name
	bool bIsAligned;			// true if species in current alignment block and ref loci was mapped
	int AlignLen;				// current alignment length including InDels
	int BlockStartLoci;			// loci at start of StartLoci alignment block
	int BlockEndLoci;			// loci at start of EndLoci block
	int NxtStartLoci;			// expected next block start loci
	int StartLoci;				// start loci 
	int EndLoci;				// end loci
	int ChromID;				// chromosome
	char szChromName[cMaxDatasetSpeciesChrom]; // chrom name
	char Strand;				// strand
	etSeqBase Seq[cMaxAlignSeqLen]; // sequence for this species - may contain InDels
} tsSpeciesAlign;

typedef struct TAG_sProcParams {
	etProcMode iProcMode;		// processing mode
	char *pszAlignFile;			// alignment file
	CMAlignFile *pAlignments;	// opened alignment file
	char *pszRefLociFile;		// file containing reference loci to map on to relative species
	CCSVFile *pCSVLoci;			// opened reference loci file
	char *pszMappedLociSpec;	// output file spec (can contain place holder '#' for species name replacement)
	int NumMappedLociFiles;		// number of mapped loci files
	char szMappedLociFiles[cMaxAlignedSpecies][_MAX_PATH];	// mapped file names after any species name replacement
	int hMappedLociFiles[cMaxAlignedSpecies];			// handles for opened map loci files
	char *pszSpeciesList;		// comma separated species list

	int MinLen;		// elements must be of at least this length
	int MaxLen;		// elements must be no longer than this length

	bool bFilterRefSpecies;		// true if reference species must match that parsed from ref loci
	bool bOverlapDetail;		// true if overlap detail to be generated

	char *pszRefSpecies;		// current ref species as parsed from ref loci
	char *pszElType;			// element type as parsed from ref loci
	int RefSrcID;				// source identifier as parsed from ref loci
	int RefFeatures;			// features as parsed from ref loci

	int RefSpeciesID;			// reference species identifier
	int RefSpeciesIdx;			// reference species index
	int NumSpecies;				// number of species being processed
	tsSpeciesAlign *pSpeciesAlign;	// pts to array of species alignments

	tsCoreLoci *pCoreLocs;			// pts to mem allocd to hold array of core loci
	int MaxCoreLocs;				// max allocd core locii
	int NumCoreLocs;				// number of core loci ptd at by pCoreLocs

	int CachCoreChromID;			// identifier of last retrieved core chromosome from pCoreChroms
	tsCoreChrom *pCoreChroms;		// pts to mem allocd to hold array of core chromosome names
	int MaxCoreChroms;				// max allocd core chroms
	int NumCoreChroms;				// number of core chroms ptd at by pCoreLocs

	CFilterRefIDs *pFilterRefIDs;	// used when filtering by RefID
} tsProcParams;

int TrimQuotes(char *pszTxt);
int ParseNumSpecies(char *pszSpeciesList,tsProcParams *pProcParams);

int 
Process(etProcMode ProcMode,		// processing mode 0: default, 1:filter by ref species in loci file
		char *pszRefLociFile,		// CSV file containing reference loci to be mapped
		char *pszFilterRefIDFile, // exclude any RefIDs in this filter file
		char *pszAlignFile,			// bio multialignment (.algn) file to process
			char *pszMappedLociSpec,// where to write out mapped loci (if contains place holder # then species name inserted)
			char *pszOverlapFile,	// where to write overlap detail
			char *pszSpeciesList,	// list of species to be processed
			int	MinLen,		// elements must be of at least this length
			int	MaxLen);		// elements must be no longer than this length

int MapLoci(tsProcParams *pProcParams);

int			// returned number of relative species loci mapped (ref species not counted!)					
ProcessStartBlock(int CurBlockID,int StartRefLoci2Map,tsProcParams *pProcParams);

int						// returned number of relative species loci mapped (ref species not counted!)				
ProcessNextBlock(int CurBlockID,tsProcParams *pProcParams);
int					// returned number of relative species loci mapped (ref species not counted!)
ProcessEndBlock(int CurBlockID,int StartBlockID,int RefLoci2Map,tsProcParams *pProcParams);

int										// number of replacements made or -1 if errors
ReplacePlaceHolders(int MaxLen,			// size of buffer ptd to by pszDest
				   char *pszDest,		// source is copied into this buffer with any replacements
				   char *pszSource,		// source with place holder chars
				   char *pszReplace,	// replacement text for any place holders in pszSource
				   char PlaceHolder);	// replacement place holder char - normally '#'

int
OutputLoci(tsProcParams *pProcParams);
int
OutputOverlaps(char *pszOverlapFile,tsProcParams *pProcParams);

int AddLoci(int RefID,char *pszChrom,int StartLoci,int Len,tsProcParams *pParams);
int AddChrom(char *pszChrom,tsProcParams *pParams);
tsCoreChrom *LocateCoreChrom(char *pszChrom,tsProcParams *pProcParams);
UINT16 GenNameHash(char *pszName);
tsCoreLoci *GetFirstLociOverlaps(char *pszChrom,	// chromosome
				int ChromOfs,		// start loci
				int ChromEndOfs,		// end loci
				tsProcParams *pProcParams); // global processing parameters

tsCoreLoci *GetLociOverlaps(char *pszChrom,	// chromosome
				int ChromOfs,		// start loci
				int ChromEndOfs,		// end loci
				tsCoreLoci *pPrevCoreLoci, // previous loci or NULL if none
				tsProcParams *pProcParams); // global processing parameters
static int CompareCoreEls( const void *arg1, const void *arg2);


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
int iProcMode;				// processing mode
char szRefLociFile[_MAX_PATH]; // CSV file containing reference loci to map
char szAlignFile[_MAX_PATH]; // .algn file containing alignment blocks
char szMappedFile[_MAX_PATH]; // mapped loci written to this file
char szOverlapFile[_MAX_PATH];		// overlap detail written to this file
char szFilterRefIDFile[_MAX_PATH];  // exclude any RefIDs in this filter file

int iMinLen;				// minimum length
int iMaxLen;				// maximum length

char szSpeciesList[2048];
int iNumSpecies;

// command line args
struct arg_lit  *help    = arg_lit0("hH","help",                "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_int *ScreenLogLevel=arg_int0("S", "ScreenLogLevel",	"<int>","Level of diagnostics written to logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",	"diagnostics log file");

struct arg_int  *ProcMode = arg_int0("m","procmode","<int>",	"processing mode 0: default, 1:filter by ref species in loci file");
struct arg_file *RefLociFile = arg_file1("i",NULL,"<file>",		"input from CSV reference species loci file");
struct arg_file *FilterRefIDFile = arg_file0("X",NULL,"<file>",	"exclude any RefIDs in this filter file");
struct arg_file *AlignFile = arg_file1("I",NULL,"<file>",		"input from .algn file");
struct arg_file *MappedFile = arg_file1("o",NULL,"<file>",	    "write mapped reference loci to this file (can contain place holder '#' for species name replacement)");
struct arg_file *OverlapFile = arg_file0("O",NULL,"<file>",	    "write overlap detail to this file");

struct arg_str  *SpeciesList = arg_str0("s","species","<string>","species list, ordered by processing priority (ref must be first)");
struct arg_int *MinLen = arg_int0("l", "minlen","<int>",		"minimum length element (default 20)");
struct arg_int *MaxLen = arg_int0("L", "maxlen","<int>",		"maximum length element (default 100000000)");

struct arg_end *end = arg_end(20);

void *argtable[] = {help,version,FileLogLevel,ScreenLogLevel,LogFile,
					ProcMode,RefLociFile,FilterRefIDFile,AlignFile,MappedFile,OverlapFile,SpeciesList,
					MinLen,MaxLen,
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

	iProcMode = ProcMode->count ? ProcMode->ival[0] : eProcModeStandard;
	if(iProcMode < eProcModeStandard || iProcMode > eProcModeFilterRef)
		{
		printf("\nError: Requested processing mode '-m%d' not supported",iProcMode);
		exit(1);
		}


	iMinLen = MinLen->count ? MinLen->ival[0] : 20;
	if(iMinLen < 0 || iMinLen > 100000000)
		{
		printf("\nError: Minimum length '-l%d' not supported",iMinLen);
		exit(1);
		}

	iMaxLen = MaxLen->count ? MaxLen->ival[0] : 100000000;
	if(iMaxLen < iMinLen || iMaxLen > 100000000)
		{
		if(iMaxLen < iMinLen)
			printf("\nError: Maximum length '-L%d' less than minimum '-l%d'",iMaxLen,iMinLen);
		else
			printf("\nError: Maximum length '-L%d' not supported",iMaxLen);
		exit(1);
		}

	if(OverlapFile->count)
		{
		strncpy(szOverlapFile,OverlapFile->filename[0],sizeof(szRefLociFile));
		szOverlapFile[sizeof(szOverlapFile)-1] = '\0';
		}
	else
		szOverlapFile[0] = '\0';

	if(FilterRefIDFile->count)
		{
		strncpy(szFilterRefIDFile,FilterRefIDFile->filename[0],sizeof(szFilterRefIDFile));
		szFilterRefIDFile[sizeof(szFilterRefIDFile)-1] = '\0';
		}
	else
		szFilterRefIDFile[0] = '\0';

	strncpy(szRefLociFile,RefLociFile->filename[0],sizeof(szRefLociFile));
	szRefLociFile[sizeof(szRefLociFile)-1] = '\0';

	strncpy(szAlignFile,AlignFile->filename[0],sizeof(szAlignFile));
	szAlignFile[sizeof(szAlignFile)-1] = '\0';

	strncpy(szMappedFile,MappedFile->filename[0],sizeof(szMappedFile));
	szMappedFile[sizeof(szMappedFile)-1] = '\0';

	if(!SpeciesList->count)
		{
		szSpeciesList[0] = '\0';
		iNumSpecies = 0;
		}
	else
		{
		strcpy(szSpeciesList,SpeciesList->sval[0]);
		TrimQuotes(szSpeciesList);
		iNumSpecies = ParseNumSpecies(szSpeciesList,NULL);
		if(iNumSpecies < 2)
			{
			if(iNumSpecies <= 0)
				printf("Error: Unable to parse species list '-s\"%s\"'",szSpeciesList);
			else
				printf("Error: At least 2 species (including reference) must be specified '-s\"%s\"'",szSpeciesList);
			exit(1);
			}
		}

	if(szOverlapFile[0] && iNumSpecies > 2)
		{
		printf("Error: Overlap detail requested, can only process reference plus single relative speciesn not '-s\"%s\"'",szSpeciesList);
		exit(1);
		}

			// now that command parameters have been parsed then initialise diagnostics log system
	if(!gDiagnostics.Open(szLogFile,(etDiagLevel)iScreenLogLevel,(etDiagLevel)iFileLogLevel,true))
		{
		printf("\nError: Unable to start diagnostics subsystem.");
		if(szLogFile[0] != '\0')
			printf(" Most likely cause is that logfile '%s' can't be opened/created",szLogFile);
		exit(1);
		}

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Version: %d.%2.2d Processing parameters:",cProgVer/100,cProgVer%100);
	switch(iProcMode) {
		case eProcModeStandard:		
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Processing mode: mapping loci for all species from alignments");
			break;

		case eProcModeFilterRef:		
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Processing mode: mapping loci for all species from alignments with ref species loci filtering");
			break;

		default:
			break;
	}

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"CSV file containing reference loci to map: '%s'",szRefLociFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"bio multialignment (.algn) file to process: '%s'",szAlignFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Write mapped loci to: '%s'",szMappedFile);
	if(iNumSpecies)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"species list: '%s'",	szSpeciesList);
	else
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"species list: * (all)");

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Minimum length elements: %d",iMinLen);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Maximum length elements: %d",iMaxLen);
	if(szOverlapFile[0])
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Generate overlap detail into: '%s'",szOverlapFile);

	if(szFilterRefIDFile[0])
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Exclude any RefIDs in this filter file: '%s'",szFilterRefIDFile);

	gStopWatch.Start();
#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	Rslt = Process((etProcMode)iProcMode,		// processing mode
					szRefLociFile,	// reference loci file
					szFilterRefIDFile, // exclude any RefIDs in this filter file
					szAlignFile,	// bio multialignment (.algn) file to process
					szMappedFile,	// where to write out mapped loci
					szOverlapFile,	// where to write out overlap detail
					szSpeciesList,	// species to process
					iMinLen,		// elements must be of at least this length
					iMaxLen);		// elements must be no longer than this length

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
		if(pProcParams != NULL && pProcParams->pSpeciesAlign != NULL)
			{
			strncpy(pProcParams->pSpeciesAlign[NumSpecies].szSpeciesName,pSpecies,cMaxDatasetSpeciesChrom);
			pProcParams->pSpeciesAlign[NumSpecies].szSpeciesName[cMaxDatasetSpeciesChrom-1] = '\0';
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
		strncpy(pProcParams->pSpeciesAlign[NumSpecies].szSpeciesName,pSpecies,cMaxDatasetSpeciesChrom);
		pProcParams->pSpeciesAlign[NumSpecies].szSpeciesName[cMaxDatasetSpeciesChrom-1] = '\0';
		}
	pszSpeciesList[-1] = Chr;
	NumSpecies++;
	}
return(NumSpecies);
}

int										// number of replacements made or -1 if MaxLen too small to hold replacements
ReplacePlaceHolders(int MaxLen,			// size of buffer ptd to by pszDest
				   char *pszDest,		// source is copied into this buffer with any replacements
				   char *pszSource,		// source with place holder chars
				   char *pszReplace,	// replacement text for any place holders in pszSource
				   char PlaceHolder)	// replacement place holder char - normally '#'
{
int NumReplacements = 0;			// number of replacements made
int ReplaceLen = (int)strlen(pszReplace);
char Chr;
if(MaxLen < 2)
	return(-1);
while(Chr = *pszSource++)
	{
	if(Chr == PlaceHolder)
		{
		if(MaxLen <= ReplaceLen)
			break;
		MaxLen -= ReplaceLen;
		strcpy(pszDest,pszReplace);
		pszDest += ReplaceLen;
		}
	else
		{
		if(MaxLen <= 1)
			break;
		MaxLen -= 1;
		*pszDest++ = Chr;
		}
	}
*pszDest = '\0';
return(MaxLen >= 1 ? NumReplacements : -1);
}


int 
Process(etProcMode ProcMode,		// processing mode 0: default, 1:filter by ref species in loci file
		char *pszRefLociFile,		// CSV file containing reference loci to be mapped
		char *pszFilterRefIDFile, // exclude any RefIDs in this filter file
		char *pszAlignFile,			// bio multialignment (.algn) file to process
			char *pszMappedLociSpec,// where to write out mapped loci (if contains place holder # then species name inserted)
			char *pszOverlapFile,	// where to write overlap detail
			char *pszSpeciesList,	// list of species to be processed
			int	MinLen,		// elements must be of at least this length
			int	MaxLen)		// elements must be no longer than this length
{
int Rslt;
tsProcParams ProcParams;
int CurSpeciesID;
char *pszCurSpeciesName;
char szCSVSpecies[2048];				// to hold comma separated species list
int SpeciesIdx;
tsSpeciesAlign *pSpeciesAlign;

bool bSubFileSpec;					// set true if pszOverlapsFileSpec contains placeholders for species name 

memset(&ProcParams,0,sizeof(ProcParams));
ProcParams.iProcMode = ProcMode;
memset(ProcParams.hMappedLociFiles,-1,sizeof(ProcParams.hMappedLociFiles));
ProcParams.pszRefLociFile = pszRefLociFile;
ProcParams.pszAlignFile = pszAlignFile;
ProcParams.pszMappedLociSpec = pszMappedLociSpec;
ProcParams.MaxLen = MaxLen;
ProcParams.MinLen = MinLen;
ProcParams.bFilterRefSpecies = ProcMode == eProcModeFilterRef ? true : false;
ProcParams.bOverlapDetail = (pszOverlapFile == NULL || pszOverlapFile[0]=='\0') ? false : true;


if(pszFilterRefIDFile != NULL && pszFilterRefIDFile[0] != '\0')
	{
	ProcParams.pFilterRefIDs = new CFilterRefIDs;
	if((Rslt=ProcParams.pFilterRefIDs->Open(pszFilterRefIDFile)) < 0)
		{
		while(ProcParams.pFilterRefIDs->NumErrMsgs())
			gDiagnostics.DiagOut(eDLFatal,gszProcName,ProcParams.pFilterRefIDs->GetErrMsg());
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open/parse input filter RefID file %s",pszFilterRefIDFile);
		return(Rslt);
		}
	}

if((ProcParams.pCSVLoci = new CCSVFile())==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create new instance of CCSVFile");
	if(ProcParams.pFilterRefIDs != NULL)
		delete ProcParams.pFilterRefIDs;
	return(eBSFerrObj);
	}
if((Rslt=ProcParams.pCSVLoci->Open(pszRefLociFile))!=eBSFSuccess)
	{
	while(ProcParams.pCSVLoci->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,ProcParams.pCSVLoci->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open input ref loci file %s",pszRefLociFile);
	delete ProcParams.pCSVLoci;
	if(ProcParams.pFilterRefIDs != NULL)
		delete ProcParams.pFilterRefIDs;
	return(Rslt);
	}

if((ProcParams.pAlignments = new CMAlignFile())==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create new instance of CMAlignFile");
	delete ProcParams.pCSVLoci;
	if(ProcParams.pFilterRefIDs != NULL)
		delete ProcParams.pFilterRefIDs;
	return(eBSFerrObj);
	}

if((Rslt=ProcParams.pAlignments->Open(pszAlignFile))!=eBSFSuccess)
	{
	while(ProcParams.pAlignments->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,ProcParams.pAlignments->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open alignment file %s",pszAlignFile);
	delete ProcParams.pCSVLoci;
	delete ProcParams.pAlignments;
	if(ProcParams.pFilterRefIDs != NULL)
		delete ProcParams.pFilterRefIDs;
	return(Rslt);
	}

ProcParams.RefSpeciesID = ProcParams.pAlignments->GetRefSpeciesID();
ProcParams.RefSpeciesIdx = -1;

szCSVSpecies[0] = '\0';
if(pszSpeciesList[0] != '\0')	// if species list was specified by user
	ProcParams.NumSpecies = ParseNumSpecies(pszSpeciesList,NULL);
else
	ProcParams.NumSpecies = ProcParams.pAlignments->GetNumSpecies(0);

if((ProcParams.pSpeciesAlign = new tsSpeciesAlign [ProcParams.NumSpecies])==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName," Unable to allocate memory for holding species");
	delete ProcParams.pCSVLoci;
	delete ProcParams.pAlignments;
	if(ProcParams.pFilterRefIDs != NULL)
		delete ProcParams.pFilterRefIDs;
	return(eBSFerrMem);
	}
memset(ProcParams.pSpeciesAlign,0,sizeof(tsSpeciesAlign) * ProcParams.NumSpecies);


if(pszSpeciesList[0] != '\0')
	{
	ParseNumSpecies(pszSpeciesList,&ProcParams);
	szCSVSpecies[0]='\0';
	pSpeciesAlign = ProcParams.pSpeciesAlign;
	for(SpeciesIdx = 0; SpeciesIdx < ProcParams.NumSpecies; SpeciesIdx++,pSpeciesAlign++)
		{
		if((Rslt = pSpeciesAlign->SpeciesID = ProcParams.pAlignments->LocateSpeciesID(pSpeciesAlign->szSpeciesName))<1)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Species '%s' not represented in %s",pSpeciesAlign->szSpeciesName,pszAlignFile);
			Rslt = ProcParams.pAlignments->GetNumSpecies();
			for(CurSpeciesID=1;CurSpeciesID<=Rslt;CurSpeciesID++)
				gDiagnostics.DiagOut(eDLFatal,gszProcName," Represented species: %s",ProcParams.pAlignments->GetSpeciesName(CurSpeciesID));
			delete ProcParams.pCSVLoci;
			delete ProcParams.pAlignments;
			if(ProcParams.pFilterRefIDs != NULL)
				delete ProcParams.pFilterRefIDs;
			delete ProcParams.pSpeciesAlign;
			return(eBSFerrEntry);
			}
		if(SpeciesIdx > 0)
			strcat(szCSVSpecies,",");
		strcat(szCSVSpecies,pSpeciesAlign->szSpeciesName);
		if(pSpeciesAlign->SpeciesID == ProcParams.RefSpeciesID)
			ProcParams.RefSpeciesIdx = SpeciesIdx;
		}
	}
else
	{
	pSpeciesAlign = ProcParams.pSpeciesAlign;
	for(CurSpeciesID=1;CurSpeciesID <= ProcParams.NumSpecies;CurSpeciesID++,pSpeciesAlign++)
		{
		pSpeciesAlign->SpeciesID = CurSpeciesID;
		if((pszCurSpeciesName = ProcParams.pAlignments->GetSpeciesName(CurSpeciesID))==NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"SpeciesID '%d' not in alignment file '%s'",CurSpeciesID,ProcParams.pszAlignFile);
			delete ProcParams.pCSVLoci;
			delete ProcParams.pAlignments;
			if(ProcParams.pFilterRefIDs != NULL)
				delete ProcParams.pFilterRefIDs;
			delete ProcParams.pSpeciesAlign;
			return(eBSFerrEntry);
			}
		strcpy(pSpeciesAlign->szSpeciesName,pszCurSpeciesName);
		if(CurSpeciesID > 1)
			strcat(szCSVSpecies,",");
		strcat(szCSVSpecies,pszCurSpeciesName);
		if(CurSpeciesID == ProcParams.RefSpeciesID)
			ProcParams.RefSpeciesIdx = CurSpeciesID-1;
		}
	}
if(ProcParams.RefSpeciesIdx == -1)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Reference species not specified in species list with '-s<specieslist>'");
	delete ProcParams.pCSVLoci;
	delete ProcParams.pAlignments;
	if(ProcParams.pFilterRefIDs != NULL)
		delete ProcParams.pFilterRefIDs;
	delete ProcParams.pSpeciesAlign;
	return(eBSFerrEntry);
	}

ProcParams.pszSpeciesList = szCSVSpecies;

bSubFileSpec = strchr(pszMappedLociSpec,'#')==NULL ? false : true;
if(!bSubFileSpec)
	{
	strcpy((char *)&ProcParams.szMappedLociFiles[0],pszMappedLociSpec);
	ProcParams.NumMappedLociFiles = 1;
#ifdef _WIN32
	if((ProcParams.hMappedLociFiles[0] = open(ProcParams.szMappedLociFiles[0], _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
	if((ProcParams.hMappedLociFiles[0] = open(ProcParams.szMappedLociFiles[SpeciesIdx],O_RDWR | O_CREAT |O_TRUNC, S_IREAD | S_IWRITE ))==-1)
#endif
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create %s - %s",ProcParams.szMappedLociFiles[0],strerror(errno));
		delete ProcParams.pCSVLoci;
		delete ProcParams.pAlignments;
		if(ProcParams.pFilterRefIDs != NULL)
			delete ProcParams.pFilterRefIDs;
		delete ProcParams.pSpeciesAlign;
		return(eBSFerrCreateFile);
		}
	}
else
	{
	bSubFileSpec = false;
	ProcParams.NumMappedLociFiles = ProcParams.NumSpecies;
	pSpeciesAlign = ProcParams.pSpeciesAlign;
	for(SpeciesIdx = 0; SpeciesIdx < ProcParams.NumSpecies; SpeciesIdx++,pSpeciesAlign++)	
		ReplacePlaceHolders(_MAX_PATH,ProcParams.szMappedLociFiles[SpeciesIdx],pszMappedLociSpec,pSpeciesAlign->szSpeciesName,'#');

	for(SpeciesIdx = 0; SpeciesIdx < ProcParams.NumSpecies; SpeciesIdx++)	
		{
#ifdef _WIN32
		if((ProcParams.hMappedLociFiles[SpeciesIdx] = open(ProcParams.szMappedLociFiles[SpeciesIdx], _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
		if((ProcParams.hMappedLociFiles[SpeciesIdx] = open(ProcParams.szMappedLociFiles[SpeciesIdx],O_RDWR | O_CREAT |O_TRUNC, S_IREAD | S_IWRITE ))==-1)
#endif
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create %s - %s",ProcParams.szMappedLociFiles[SpeciesIdx],strerror(errno));
			while(--SpeciesIdx >= 0)
				{
				if(ProcParams.hMappedLociFiles[SpeciesIdx] != -1)
					close(ProcParams.hMappedLociFiles[SpeciesIdx]);
				ProcParams.hMappedLociFiles[SpeciesIdx] = -1;
				}
			delete ProcParams.pCSVLoci;
			delete ProcParams.pAlignments;
			if(ProcParams.pFilterRefIDs != NULL)
				delete ProcParams.pFilterRefIDs;
			delete ProcParams.pSpeciesAlign;
			return(eBSFerrCreateFile);
			}
		}
	}

Rslt = MapLoci(&ProcParams);
if(Rslt >= 0 && ProcParams.bOverlapDetail && ProcParams.NumCoreLocs)
	Rslt = OutputOverlaps(pszOverlapFile,&ProcParams);


for(SpeciesIdx = 0; SpeciesIdx < ProcParams.NumSpecies; SpeciesIdx++)
	{
	if(ProcParams.hMappedLociFiles[SpeciesIdx] != -1)
		close(ProcParams.hMappedLociFiles[SpeciesIdx]);
	ProcParams.hMappedLociFiles[SpeciesIdx] = -1;
	}

if(ProcParams.pCSVLoci != NULL)
	delete ProcParams.pCSVLoci;
if(ProcParams.pFilterRefIDs != NULL)
	delete ProcParams.pFilterRefIDs;
if(ProcParams.pAlignments != NULL)
	delete ProcParams.pAlignments;
if(ProcParams.pSpeciesAlign != NULL)
	delete ProcParams.pSpeciesAlign;
if(ProcParams.pCoreLocs != NULL)
	delete ProcParams.pCoreLocs;
if(ProcParams.pCoreChroms != NULL)
	delete ProcParams.pCoreChroms;
return(eBSFSuccess);
}

int
OutputOverlaps(char *pszOverlapFile,tsProcParams *pProcParams)
{
int CurIdx;
int ProbeIdx;
int MaxSizeCore;
int CurSizeCore;
tsCoreLoci *pCurCore;
tsCoreLoci *pProbe;
int NumOverlaps;
int TotOverlaps;
int OverlapRefIDs[cMaxOverlapIDs];
int Len;
char szLineBuff[32000];
int hOverFile = -1;

#ifdef _WIN32
if((hOverFile = open(pszOverlapFile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
	if((hOverFile = open(pszOverlapFile,O_RDWR | O_CREAT |O_TRUNC, S_IREAD | S_IWRITE ))==-1)
#endif
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create %s - %s",pszOverlapFile,strerror(errno));
	return(eBSFerrCreateFile);
	}

pCurCore = pProcParams->pCoreLocs;
MaxSizeCore = 0;
for(CurIdx = 0; CurIdx < pProcParams->NumCoreLocs; CurIdx++,pCurCore++)
	{
	CurSizeCore = 1 + pCurCore->EndLoci - pCurCore->StartLoci;
	if(CurSizeCore > MaxSizeCore)
		MaxSizeCore = CurSizeCore;
	}
pCurCore = pProcParams->pCoreLocs;
for(CurIdx = 0; CurIdx < pProcParams->NumCoreLocs; CurIdx++,pCurCore++)
	{
	NumOverlaps = 0;
	TotOverlaps = 0;
	if(CurIdx > 0)
		{
		pProbe = pCurCore - 1;
		for(ProbeIdx = CurIdx-1; ProbeIdx >= 0; ProbeIdx--,pProbe--)
			{
			if(pProbe->ChromID != pCurCore->ChromID)
				break;
			if(pProbe->EndLoci >= pCurCore->StartLoci)
				{
				if(NumOverlaps < cMaxOverlapIDs)
					OverlapRefIDs[NumOverlaps++] = pProbe->RefID;
				TotOverlaps++;
				}
			if((pCurCore->StartLoci - pProbe->StartLoci) > MaxSizeCore)
				break;
			}
		}
	if(CurIdx < pProcParams->NumCoreLocs-1)
		{
		pProbe = pCurCore + 1;
		for(ProbeIdx = CurIdx + 1; ProbeIdx < pProcParams->NumCoreLocs; ProbeIdx++,pProbe++)
			{
			if(pProbe->ChromID != pCurCore->ChromID)
				break;
			if(pProbe->StartLoci > pCurCore->EndLoci)
				break;
			if(NumOverlaps < cMaxOverlapIDs)
				OverlapRefIDs[NumOverlaps++] = pProbe->RefID;
			TotOverlaps++;
			}
		}

	Len = sprintf(szLineBuff,"%d,%d,%d,",pCurCore->RefID,TotOverlaps,NumOverlaps);
	if(!NumOverlaps)
		Len += sprintf(&szLineBuff[Len],"\"\"\n");
	else
		{
		Len += sprintf(&szLineBuff[Len],"\"%d",OverlapRefIDs[0]);
		for(ProbeIdx = 1; ProbeIdx < NumOverlaps; ProbeIdx++)
			Len += sprintf(&szLineBuff[Len],",%d",OverlapRefIDs[ProbeIdx]);
		Len += sprintf(&szLineBuff[Len],"\"\n");
		}
	CUtility::SafeWrite(hOverFile,szLineBuff,Len);
	}
if(hOverFile != -1)
	close(hOverFile);
return(eBSFSuccess);
}

int
MapLoci(tsProcParams *pProcParams)
{
int Rslt;
CCSVFile *pCSV = pProcParams->pCSVLoci;
CMAlignFile *pAlignments = pProcParams->pAlignments;
int NumFields;
int SrcID;
int Features;
char *pszChrom;
char *pszElType;
char *pszRefSpecies;
int StartLoci;
int EndLoci;
int Len;
int RefChromID;
int StartBlockID;
int NextBlockID;
int EndBlockID;
int NumElsRead =0;		// number of elements before filtering
int NumElsAccepted =0;	// number of elements accepted after filtering
int NumFiltRefIDs =0;	// number of elements filtered out because of matching RefID
int NumFiltLen=0;		// number of elements filtered out because of length
int NumFiltSpecies=0;	// number of elements filtered out because of species
int NumFiltRefLoci=0;	// number of elements filtered out because of ref chrom.start or end
int NumFiltRelLoci=0;	// number of elements filtered out because of no mapped rel start or end



while((Rslt=pCSV->NextLine()) > 0)	// onto next line containing fields
	{
	NumFields = pCSV->GetCurFields();
	if(NumFields < 9)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expected at least 9 fields in '%s', GetCurFields() returned '%d'",pProcParams->pszRefLociFile,NumFields);
		return(eBSFerrFieldCnt);
		}
	NumElsRead += 1;

	// apply filtering - don't process those loci which don't meet filtering criteria
	pCSV->GetInt(1,&SrcID);
	if(pProcParams->pFilterRefIDs != NULL && pProcParams->pFilterRefIDs->Locate(SrcID))
		{	
		NumFiltRefIDs += 1;
		continue;
		}

	pProcParams->pCSVLoci->GetInt(7,&Len);
	if(Len < pProcParams->MinLen || Len > pProcParams->MaxLen)
		{
		NumFiltLen += 1;
		continue;
		}

	pCSV->GetText(3,&pszRefSpecies);

	if(pProcParams->bFilterRefSpecies && stricmp(pszRefSpecies,pProcParams->pSpeciesAlign[pProcParams->RefSpeciesIdx].szSpeciesName))
		{
		NumFiltSpecies += 1;
		continue;
		}

	pCSV->GetText(2,&pszElType);
	pCSV->GetText(4,&pszChrom);
	pCSV->GetInt(5,&StartLoci);
	pCSV->GetInt(6,&EndLoci);
	pCSV->GetInt(9,&Features);

	pProcParams->RefSrcID = SrcID;			// gets copied into output loci file
	pProcParams->pszElType = pszElType;		// gets copied into output loci file
	pProcParams->RefFeatures = Features;		// gets copied into output loci file
	pProcParams->pszRefSpecies = pszRefSpecies; // gets copied into output loci file

	RefChromID = pAlignments->LocateChromID(pProcParams->RefSpeciesID,pszChrom);
	if(RefChromID == eBSFerrChrom)			// if not located then can't be any alignment onto this chromosome
		{
		NumFiltRefLoci += 1;
		continue;
		}
	else
		if(RefChromID <= 0)
			break;

	// determine starting block containing reference start chrom.loci
	if((Rslt = StartBlockID = pAlignments->GetBlockID(RefChromID,StartLoci)) < 1)
		{
		NumFiltRefLoci += 1;
		continue;								// no block with the chrom and start loci
		}

	// process starting block looking for at least one relative species start
	if((Rslt = ProcessStartBlock(StartBlockID,StartLoci,pProcParams)) < 0)
		break;
	if(!Rslt)	// if no relative species start then process next ref start chrom.loci
		{
		NumFiltRelLoci += 1;
		continue;								// no block with the chrom and start loci
		}

	// determine end block containing reference end chrom.loci
	if((Rslt = EndBlockID = pAlignments->GetBlockID(RefChromID,EndLoci)) < 1)
		{
		NumFiltRefLoci += 1;
		continue;								// no block with the chrom and end loci
		}

	// if intermediate blocks between EndBlockID and StartBlockID then ensure rel species are contiguous 
	if(EndBlockID != StartBlockID)
		{
		NextBlockID = StartBlockID;
		while((Rslt = NextBlockID = pAlignments->NxtBlock(NextBlockID))!= EndBlockID)
			{
			if(Rslt <= 0)
				break;
			if((Rslt = ProcessNextBlock(NextBlockID,pProcParams)) <= 0)
				break;
			}
		}
	if(Rslt < 0)
		break;
	if(!Rslt)	// if no contiguous relative species then process next ref start chrom.loci
		{
		NumFiltRelLoci += 1;
		continue;							
		}

 	// process ending block looking for at least one relative species
	if((Rslt = ProcessEndBlock(EndBlockID,StartBlockID,EndLoci,pProcParams)) < 0)
		break;
	if(!Rslt)
		{
		NumFiltRelLoci += 1;
		continue;
		}

	if((Rslt=OutputLoci(pProcParams)) < 0)
		break;
	NumElsAccepted += 1;
	}

if(Rslt >= 0)
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"%s: Read: %d, Accepted: %d, FiltRefIDs: %d, FiltLen: %d, FiltSpecies: %d, FiltRefLoci: %d, FiltRelLoci: %d",
		pProcParams->pszRefLociFile,NumElsRead,NumElsAccepted,NumFiltRefIDs,NumFiltLen,NumFiltSpecies,NumFiltRefLoci,NumFiltRelLoci);

// sort cores by chrom.start.end
if(Rslt >= 0 && pProcParams->bOverlapDetail)
	{
	if(pProcParams->NumCoreLocs > 1)
		qsort(pProcParams->pCoreLocs,pProcParams->NumCoreLocs,sizeof(tsCoreLoci),CompareCoreEls);
	// add a sentenil element as last element (allocation ensures at least 1 free element!) 
	if(pProcParams->NumCoreLocs)
		pProcParams->pCoreLocs[pProcParams->NumCoreLocs].ChromID = -1;
	}

return(Rslt);
}

int
OutputLoci(tsProcParams *pProcParams)
{
char szBuff[cMaxReadLen+1];
int BuffLen;
tsSpeciesAlign *pRelSpecies;
int SpeciesIdx;
int StartLoci;
int EndLoci;
char *pszChrom;
int hFile;
int Rslt;
int Cnt = 0;
pRelSpecies = &pProcParams->pSpeciesAlign[0];
for(SpeciesIdx = 0; SpeciesIdx < pProcParams->NumSpecies; SpeciesIdx++,pRelSpecies++)
	{
	if(pRelSpecies->bIsAligned == false)
		continue;
	if(pProcParams->NumMappedLociFiles == 1)
		hFile = pProcParams->hMappedLociFiles[0];
	else
		hFile = pProcParams->hMappedLociFiles[SpeciesIdx];

	if(pRelSpecies->Strand == '-')
		{
		EndLoci = pProcParams->pAlignments->MapChromOfsToOtherStrand(pRelSpecies->ChromID,pRelSpecies->StartLoci);
		StartLoci = pProcParams->pAlignments->MapChromOfsToOtherStrand(pRelSpecies->ChromID,pRelSpecies->EndLoci);
		}
	else
		{
		EndLoci = pRelSpecies->EndLoci;
		StartLoci = pRelSpecies->StartLoci;
		}

	pszChrom = pProcParams->pAlignments->GetChromName(pRelSpecies->ChromID);
	if(hFile != -1)
		{
		BuffLen = sprintf(szBuff,"%d,\"%s\",\"%s\",\"%s\",%d,%d,%d,\"%s\",%d,\"%c\"\n",
			pProcParams->RefSrcID,pProcParams->pszElType,pRelSpecies->szSpeciesName,
			pszChrom,
			StartLoci,EndLoci,1 + EndLoci - StartLoci,
			pProcParams->pszRefSpecies,pProcParams->RefFeatures,
			pRelSpecies->Strand);
		CUtility::SafeWrite(hFile,szBuff,BuffLen);
		Cnt+=1;
		}

	if(pProcParams->bOverlapDetail && SpeciesIdx != pProcParams->RefSpeciesIdx)
		{
		if((Rslt=AddLoci(pProcParams->RefSrcID,pszChrom,StartLoci,1 + EndLoci - StartLoci,pProcParams)) < 1)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddLoci() returned error '%d'", Rslt);
			return(Rslt);
			}
		}
	}
return(Cnt);
}


int			// returned number of relative species loci mapped (ref species not counted!)					
ProcessStartBlock(int CurBlockID,int StartRefLoci2Map,tsProcParams *pProcParams)
{
etSeqBase *pSeq;
int CurRefAlignLen;
int SpeciesIdx;
tsSpeciesAlign *pRefSpecies;
tsSpeciesAlign *pRelSpecies;
int Psn;
int NumMapped;
etSeqBase *pRefSeq;
etSeqBase *pRelSeq;
etSeqBase RefBase;
etSeqBase RelBase;

CMAlignFile *pAlignments = pProcParams->pAlignments;
CurRefAlignLen = pAlignments->GetAlignLen(CurBlockID,pProcParams->RefSpeciesID);

NumMapped = 0;
pRefSpecies = NULL;
pRelSpecies = &pProcParams->pSpeciesAlign[0];
for(SpeciesIdx = 0; SpeciesIdx < pProcParams->NumSpecies; SpeciesIdx++,pRelSpecies++)
	{
	if((pRelSpecies->ChromID = pAlignments->GetRelChromID(CurBlockID,pRelSpecies->SpeciesID))<=0)
		{	
		pRelSpecies->bIsAligned = false;
		continue;
		}
	pRelSpecies->bIsAligned = true;		// assume will get full alignment
	pRelSpecies->Strand = pAlignments->GetStrand(CurBlockID,pRelSpecies->SpeciesID);
	pRelSpecies->BlockStartLoci = pAlignments->GetRelChromOfs(CurBlockID,pRelSpecies->SpeciesID);
	pRelSpecies->NxtStartLoci = pAlignments->GetRelChromEndOfs(CurBlockID,pRelSpecies->SpeciesID) + 1;
	pRelSpecies->AlignLen = CurRefAlignLen;
	pSeq = pAlignments->GetSeq(CurBlockID,pRelSpecies->SpeciesID);
	memcpy(pRelSpecies->Seq,pSeq,CurRefAlignLen);
	if(pRelSpecies->SpeciesID == pProcParams->RefSpeciesID)
		pRefSpecies = pRelSpecies;
	else
		NumMapped += 1;
	}
if(!NumMapped || pRefSpecies == NULL)
	return(0);

// now to determine the start mapped loci
pRelSpecies = &pProcParams->pSpeciesAlign[0];
for(SpeciesIdx = 0; SpeciesIdx < pProcParams->NumSpecies; SpeciesIdx++,pRelSpecies++)
	{
	if(pRelSpecies->SpeciesID == pProcParams->RefSpeciesID ||
		pRelSpecies->bIsAligned == false)
		continue;
	pRefSpecies->StartLoci = pRefSpecies->BlockStartLoci;
	pRelSpecies->StartLoci = pRelSpecies->BlockStartLoci;
	pRefSeq = pRefSpecies->Seq;
	pRelSeq = pRelSpecies->Seq;

	for(Psn = 0; Psn < CurRefAlignLen; Psn++, pRefSeq++, pRelSeq++)
		{
		RefBase = (*pRefSeq & ~cRptMskFlg);
		RelBase = (*pRelSeq & ~cRptMskFlg);

		if(RefBase <= eBaseN)
			{
			if(pRefSpecies->StartLoci == StartRefLoci2Map)
				{
				if(RelBase == eBaseUndef || RelBase >= eBaseEOS)
					{
					NumMapped -= 1;
					pRelSpecies->bIsAligned = false;
					}
				break;
				}
			pRefSpecies->StartLoci += 1;
			}
		if(RelBase <= eBaseN)
			pRelSpecies->StartLoci += 1;
		}
	}
return(NumMapped);
}

int						// returned number of relative species loci mapped (ref species not counted!)				
ProcessNextBlock(int CurBlockID,tsProcParams *pProcParams)
{
int CurRefAlignLen;
int SpeciesIdx;
tsSpeciesAlign *pRelSpecies;
int NumMapped;
bool bRefOk;			// set if reference alignment still contiguous
CMAlignFile *pAlignments = pProcParams->pAlignments;
CurRefAlignLen = pAlignments->GetAlignLen(CurBlockID,pProcParams->RefSpeciesID);

bRefOk = false;
NumMapped = 0;
pRelSpecies = &pProcParams->pSpeciesAlign[0];
for(SpeciesIdx = 0; SpeciesIdx < pProcParams->NumSpecies; SpeciesIdx++,pRelSpecies++)
	{
	if(!pRelSpecies->bIsAligned)
		continue;

	// intemediate must be colinear and on same strand and chromosome to be contiguous
	if(pRelSpecies->ChromID != pAlignments->GetRelChromID(CurBlockID,pRelSpecies->SpeciesID))
		{	
		pRelSpecies->bIsAligned = false;
		continue;
		}
	if(pRelSpecies->Strand != pAlignments->GetStrand(CurBlockID,pRelSpecies->SpeciesID))
		{	
		pRelSpecies->bIsAligned = false;
		continue;
		}

	if(pRelSpecies->NxtStartLoci != pAlignments->GetRelChromOfs(CurBlockID,pRelSpecies->SpeciesID))
		{	
		pRelSpecies->bIsAligned = false;
		continue;
		}
	pRelSpecies->NxtStartLoci = pAlignments->GetRelChromEndOfs(CurBlockID,pRelSpecies->SpeciesID) + 1;
 	if(pRelSpecies->SpeciesID == pProcParams->RefSpeciesID)
		bRefOk = true;
	else
		NumMapped += 1;
	}
return(bRefOk ? NumMapped : 0);
}

int					// returned number of relative species loci mapped (ref species not counted!)
ProcessEndBlock(int CurBlockID,int StartBlockID,int RefLoci2Map,tsProcParams *pProcParams)
{
etSeqBase *pSeq;
int CurRefAlignLen;
int SpeciesIdx;
tsSpeciesAlign *pRefSpecies;
tsSpeciesAlign *pRelSpecies;
int Psn;
int NumMapped;

etSeqBase *pRefSeq;
etSeqBase *pRelSeq;
etSeqBase RefBase;
etSeqBase RelBase;

CMAlignFile *pAlignments = pProcParams->pAlignments;
CurRefAlignLen = pAlignments->GetAlignLen(CurBlockID,pProcParams->RefSpeciesID);


NumMapped = 0;
pRefSpecies = NULL;
pRelSpecies = &pProcParams->pSpeciesAlign[0];
for(SpeciesIdx = 0; SpeciesIdx < pProcParams->NumSpecies; SpeciesIdx++,pRelSpecies++)
	{
	if(!pRelSpecies->bIsAligned)
		continue;

	if(pRelSpecies->ChromID != pAlignments->GetRelChromID(CurBlockID,pRelSpecies->SpeciesID))
		{	
		pRelSpecies->bIsAligned = false;
		continue;
		}
	if(pRelSpecies->Strand != pAlignments->GetStrand(CurBlockID,pRelSpecies->SpeciesID))
		{	
		pRelSpecies->bIsAligned = false;
		continue;
		}

	pRelSpecies->BlockEndLoci = pAlignments->GetRelChromOfs(CurBlockID,pRelSpecies->SpeciesID);
	if(StartBlockID != CurBlockID && pRelSpecies->NxtStartLoci != pRelSpecies->BlockEndLoci)
		{	
		pRelSpecies->bIsAligned = false;
		continue;
		}
	
	pRelSpecies->AlignLen = CurRefAlignLen;
	pSeq = pAlignments->GetSeq(CurBlockID,pRelSpecies->SpeciesID);
	memcpy(pRelSpecies->Seq,pSeq,CurRefAlignLen);
	if(pRelSpecies->SpeciesID == pProcParams->RefSpeciesID)
		pRefSpecies = pRelSpecies;
	else
		NumMapped += 1;
	}
if(!NumMapped || pRefSpecies == NULL)
	return(0);

// now to determine the mapped loci!!!!!
pRelSpecies = &pProcParams->pSpeciesAlign[0];
for(SpeciesIdx = 0; SpeciesIdx < pProcParams->NumSpecies; SpeciesIdx++,pRelSpecies++)
	{
	if(pRelSpecies->SpeciesID == pProcParams->RefSpeciesID ||
		pRelSpecies->bIsAligned == false)
		continue;
	pRefSpecies->EndLoci = pRefSpecies->BlockEndLoci;
	pRelSpecies->EndLoci = pRelSpecies->BlockEndLoci;
	pRefSeq = pRefSpecies->Seq;
	pRelSeq = pRelSpecies->Seq;

	for(Psn = 0; Psn < CurRefAlignLen; Psn++, pRefSeq++, pRelSeq++)
		{
		RefBase = (*pRefSeq & ~cRptMskFlg);
		RelBase = (*pRelSeq & ~cRptMskFlg);

		if(RefBase <= eBaseN)
			{
			if(pRefSpecies->EndLoci == RefLoci2Map)
				{
				if(RelBase == eBaseUndef || RelBase >= eBaseEOS)
					{
					NumMapped -= 1;
					pRelSpecies->bIsAligned = false;
					}
				break;
				}
			pRefSpecies->EndLoci += 1;
			}
		if(RelBase <= eBaseN)
			pRelSpecies->EndLoci += 1;
		}
	}

return(NumMapped);
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


tsCoreChrom *
LocateCoreChrom(char *pszChrom,
				tsProcParams *pProcParams) // global processing parameters
{
tsCoreChrom *pChrom;
UINT16 Hash;
int Idx;
if(pProcParams->pCoreChroms == NULL || pProcParams->NumCoreChroms < 1)
	return(NULL);
Hash = GenNameHash(pszChrom);

// quick check against cached chrom name
if(pProcParams->CachCoreChromID > 0)
	{
	pChrom = &pProcParams->pCoreChroms[pProcParams->CachCoreChromID-1];
	if(Hash == pChrom->Hash && !stricmp(pChrom->szChrom,pszChrom))
		return(pChrom);
	}
pChrom = pProcParams->pCoreChroms;
for(Idx = 0; Idx < pProcParams->NumCoreChroms; Idx++,pChrom++)
	{
	if(Hash == pChrom->Hash && !stricmp(pChrom->szChrom,pszChrom))
		{
		pProcParams->CachCoreChromID = pChrom->ChromID;
		return(pChrom);
		}
	}
return(NULL);
}


int
AddChrom(char *pszChrom,tsProcParams *pParams)
{
tsCoreChrom *pChrom;
int Idx;
UINT16 Hash;
Hash = GenNameHash(pszChrom);

if(pParams->pCoreChroms != NULL && pParams->NumCoreChroms)
	{
	if(pParams->CachCoreChromID > 0)
		{
		pChrom = &pParams->pCoreChroms[pParams->CachCoreChromID-1];
		if((pChrom->Hash == Hash) && !stricmp(pszChrom,pChrom->szChrom))
			return(pChrom->ChromID);
		}

	pChrom = pParams->pCoreChroms;
	for(Idx = 0; Idx < pParams->NumCoreChroms; Idx++,pChrom++)
		{
		if((pChrom->Hash == Hash) && !stricmp(pszChrom,pChrom->szChrom))
			return(pParams->CachCoreChromID = pChrom->ChromID);
		}
	}

if(pParams->pCoreChroms == NULL || pParams->NumCoreChroms == pParams->MaxCoreChroms)
	{
	if((pChrom = new tsCoreChrom[pParams->MaxCoreChroms + cCoreChromAllocChunk])==NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to alloc memory for chrom names");
		return(eBSFerrMem);
		}
	if(pParams->pCoreChroms != NULL)
		{
		memcpy(pChrom,pParams->pCoreChroms,sizeof(tsCoreChrom) * pParams->NumCoreChroms);
		delete pParams->pCoreChroms;
		}
	else
		{
		pParams->MaxCoreChroms = 0;
		pParams->NumCoreLocs = 0;
		}
	pParams->pCoreChroms = pChrom;
	pParams->MaxCoreChroms += cCoreChromAllocChunk;
	}
pChrom = &pParams->pCoreChroms[pParams->NumCoreChroms++];
pChrom->ChromID = pParams->NumCoreChroms;
pChrom->Hash = Hash;
pChrom->pCoreLoci = NULL;
strcpy(pChrom->szChrom,pszChrom);
return(pParams->CachCoreChromID = pChrom->ChromID);
}


int
AddLoci(int RefID,char *pszChrom,int StartLoci,int Len,tsProcParams *pParams)
{
tsCoreLoci *pNxtCore;
int ChromID;

if((ChromID = AddChrom(pszChrom,pParams)) < 1)
	return(ChromID);

if(pParams->pCoreLocs == NULL || (pParams->NumCoreLocs + 1) == pParams->MaxCoreLocs)
	{
	if((pNxtCore = new tsCoreLoci[pParams->MaxCoreLocs + cCoreLociAllocChunk])==NULL)
		{
		return(eBSFerrMem);
		}
	if(pParams->pCoreLocs != NULL)
		{
		memcpy(pNxtCore,pParams->pCoreLocs,sizeof(tsCoreLoci) * pParams->NumCoreLocs);
		pParams->pCoreLocs = pNxtCore;
		pParams->MaxCoreLocs += cCoreLociAllocChunk;
		}
	else
		{
		pParams->MaxCoreLocs = cCoreLociAllocChunk;
		pParams->NumCoreLocs = 0;
		}
	pParams->pCoreLocs = pNxtCore;
	}
pNxtCore = &pParams->pCoreLocs[pParams->NumCoreLocs++];
pNxtCore->RefID = RefID;
pNxtCore->ChromID = ChromID;
pNxtCore->StartLoci = StartLoci;
pNxtCore->EndLoci = StartLoci + Len - 1;
return(pParams->NumCoreLocs);
}



tsCoreLoci *
GetLociOverlaps(char *pszChrom,			// chromosome
				int ChromOfs,			// start loci
				int ChromEndOfs,		// end loci
				tsCoreLoci *pPrevCoreLoci, //
				tsProcParams *pProcParams) // global processing parameters
{
tsCoreLoci *pLoci;

if(pPrevCoreLoci == NULL)
	return(GetFirstLociOverlaps(pszChrom,ChromOfs,ChromEndOfs,pProcParams));

pLoci = pPrevCoreLoci + 1;
while(pLoci->ChromID > 0 && pLoci->ChromID == pPrevCoreLoci->ChromID)
	{
	if(pLoci->StartLoci > ChromEndOfs)
		break;

	if(pLoci->StartLoci <= ChromEndOfs && pLoci->EndLoci >= ChromOfs)
		return(pLoci);
	pLoci += 1;
	}
return(NULL);
}


tsCoreLoci *
GetFirstLociOverlaps(char *pszChrom,	// species chromosome
				int ChromOfs,		// start loci
				int ChromEndOfs,		// end loci
				tsProcParams *pProcParams) // global processing parameters
{
tsCoreChrom *pChrom;
tsCoreLoci *pLoci;

if((pChrom = LocateCoreChrom(pszChrom,pProcParams))==NULL)
   return(NULL);
if((pLoci = pChrom->pCoreLoci)==NULL)
	return(NULL);
while(pLoci->ChromID > 0 && pLoci->ChromID == pChrom->ChromID)
	{
	if(pLoci->StartLoci > ChromEndOfs)
		break;
	if(pLoci->StartLoci <= ChromEndOfs && pLoci->EndLoci >= ChromOfs)
		return(pLoci);
	pLoci += 1;
	}
return(NULL);
}



// CompareCoreEls
// Used to sort core loci elements 
int 
CompareCoreEls( const void *arg1, const void *arg2)
{
tsCoreLoci *pEl1 = (tsCoreLoci *)arg1;
tsCoreLoci *pEl2 = (tsCoreLoci *)arg2;
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

