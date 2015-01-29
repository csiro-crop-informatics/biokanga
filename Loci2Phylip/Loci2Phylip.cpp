// Loci2Phylip.cpp : Defines the entry point for the console application.


#include "stdafx.h"
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#if _WIN32

#include "../libbiokanga/commhdrs.h"
#else
#include "../libbiokanga/commhdrs.h"
#endif

const unsigned int cProgVer = 201;			// increment with each release

const int cMAFminSpecies = 2;				// minimum number of species that can be handled in alignment
const int cMAFmaxSpecies = 50;				// maximum number of species that can be handled in alignment
const int cDNADistRsltLen = cMAFmaxSpecies * cMAFmaxSpecies * 10;	// max length of dnadist output file containing distances

const int cMAtotMaxSeqAlignLen = 0x0ffffff; // total (over all aligned species) max seq length that can be buffered in concatenated seqs

const int cMinCoreLen = 4;				// allow core lengths to be specified down to cMinCoreLen
const int cDfltCoreLen= 50;				// if core lengths not specified then default to cDfltMinCoreLen
const int cMaxCoreLen = 10000;			// minimum core lengths can be specified upto this length

const int cMaxFlankLen = 10000;			// allow flanks of upto this many bases to be requested

const int cMaxIncludeFiles = 10;		// maximun number of include region filter files
const int cMaxExcludeFiles = 10;		// maximun number of exclude region filter files

const int cMaxIncludeChroms = 20;		// max number of include chromosomes regular expressions
const int cMaxExcludeChroms = 20;		// max number of exclude chromosomes regular expressions

const int cCoreLociAllocChunk = 10000;	// grow core loci array by this many elements
const int cCoreChromAllocChunk = 1000;	// grow core chrom array by this many elements

const int cMaxTotPhylipSeqLen = 0x0fffffff; // allow total of all Phylip concatenated sequences to be at most this many bases

// processing mode
typedef enum eProcMode {
	eProcModePhylip,					// output concatenated species sequences as single Phylip
	eProcModeMultiPhylip				// output each aligned block into separate Phylip file
	} etProcMode;

// expected input loci file type
typedef enum eLociFileType {
	eLociFileTypeCSV = 0,				// default processing is from CSV file
	eLociFileTypeBED,					// or from biobed feature file
	eLociFileTypeNone					// loci specific filtering not used
} etLociFileType;

typedef enum eBEDRegion {
	ePhylipRAny = 0,				// process any (all) region
	ePhylipRIntergenic,				// only process intergenic
	ePhylipRExons,					// only process exons
	ePhylipRIntrons,				// only process introns
	ePhylipRCDS,					// only process CDSs
	ePhylipUTR,						// only process UTRs
	ePhylip5UTR,					// only process 5'UTRs
	ePhylip3UTR						// only process 3'UTRs
} etBEDRegion;


typedef struct TAG_sExcludeSpeciesChrom {
	int SpeciesID;			// which species
	int ChromID;			// chromosome is not to be processed
	} tsExcludeSpeciesChrom;

typedef struct TAG_sCoreLoci {
	int RefID;
	int ChromID;
	char Strand;
	int StartLoci;
	int EndLoci;
} tsCoreLoci;

typedef struct TAG_sCoreChrom {
	int ChromID;							// unique identifier for this chromosome
	tsCoreLoci *pCoreLoci;					// pts to first core on this chromosome
	UINT16 Hash;								// hash over chromosome name
	char szChrom[cMaxDatasetSpeciesChrom];			// core chromosome
} tsCoreChrom;

typedef struct TAG_sProcParams 
	{
	etProcMode ProcMode;			// processing mode

	char *pszExeFile;				// process with this phylip program
	char *pszExeParFile;			// using this parameter file
	char *pszExeInFile;			// Phylip processing input file
	bool bWCExeInFile;			// true if processing each aligned sequence individually and pszExeInFileTpl contains sprintf template to use
	char *pszExeInFileTpl;		// sprintf Phylip file template if multiple Phylip input files
	char *pszExeRsltFile;		// Phylip processing to capture phylip stdout into
	bool bWCExeRsltFile;		// true if processing each aligned sequence individually and pszExeRsltFileTpl contains sprintf template to use
	char *pszExeRsltFileTpl;	// sprintf Phylip file template if multiple Phylip stdout files

	int MinNumSpecies;				// minimum number of species required in alignment block
	int MaxNumSpecies;				// maximum number of species required in alignment block
	int NumBlockSpecies;			// number of species in currently processed block(s)
	int AlignBlockID;				// current alignment block identifier
	char *pszSpeciesList;			// comma separated species list starting with reference species 
	int RefSpeciesIdx;				// current index into pSeqs for the reference species
	char szSpecies[cMaxAlignedSpecies][cMaxDatasetSpeciesChrom];	// species names of interest - other species are sloughed
									// only alignments with species included will be processed
									// first species is the reference species
	int MaxAlignIdxSpecies;			// pSeqs[MaxAlignIdxSpecies-1] is last actual alignment sequence
	char szRefChrom[cMaxDatasetSpeciesChrom]; // current reference chromosome
	int RefChromID;					// current reference chromosome identifier
	int NxtOutputOffset;			// when to next output results - i.e end of current window
	int SeqOfs;						// offset into sequences
	int SeqLen;						// sequence length
	etSeqBase *pSeqs[cMaxAlignedSpecies];  // ptrs to each alignment sequence

	char *pszOutputFile;			// output file to create/write
	bool bWrtHdr;					// true if header to be generated in output file
	bool bWCOutputFile;				// true if processing each aligned sequence individually and pszOutputFileTpl contains sprintf template to use
	char *pszOutputFileTpl;			// sprintf Phylip file template if multiple Phylip files to generate
	
	int hOutFile;					// opened file handle for output Phylip file or results after phylip processing
	int PhylipLen;					// current length of each Phylip sequence
	int MaxPhylipSeqLen;			// Phylip sequences can contain at most this many bases before being truncated
	int PhylipBlockStart;			// length when block was started - used for checkpointing in case rollback required
	unsigned char *pPhylipSeqs[cMaxAlignedSpecies];  // ptrs to each Phylip sequence
	bool bPhylipLenTrunc;			// true if concatenated Phylip sequences are being truncated because overlength 

	int MaxSeqAlignLen;				// max length alignment which can be processed (how much mem was alloc'd to pSeq[n])			
	int NumIncludes;				// number of biobed files containing regions to include
	int NumExcludes;				// number of biobed files containing regions to exclude
	CBEDfile *pIncludes[cMaxIncludeFiles];	// if opened biobed files for regions to include - all other regions are to be excluded
	CBEDfile *pExcludes[cMaxExcludeFiles];	// if opened biobed files for regions to exclude 
	int hRsltsFile;					// write results into this Phylip file
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
	etLociFileType LociFileType;	// expected loci file type
	char *pszLociFile;				// name of file containing core loci
	etBEDRegion FiltRegion;			// regions 0=ANY,1=Exons,2=Introns,3=CDS,4=UTRs,5=5'UTR,6=3'UTR
	char *pszRegionFile;			// region charaterisation biobed file
	CBEDfile *pBEDregions;			// if not NULL then opened biobed file with regional characterisations
	CBEDfile *pBEDloci;				// if not NULL then opened biobed file with core loci
	CCSVFile *pCSVloci;				// if not NULL then opened CSV file with core loci

	tsCoreLoci *pCoreLocs;			// pts to mem allocd to hold array of core loci
	int MaxCoreLocs;				// max allocd core locii
	int NumCoreLocs;				// number of core loci ptd at by pCoreLocs

	int CachCoreChromID;			// identifier of last retrieved core chromosome from pCoreChroms
	tsCoreChrom *pCoreChroms;		// pts to mem allocd to hold array of core chromosome names
	int MaxCoreChroms;				// max allocd core chroms
	int NumCoreChroms;				// number of core chroms ptd at by pCoreLocs

	int	MinCoreLen;					// minimum core length required
	int	MaxCoreLen;					// maximum core length required
	int FlankLen;					// flank loci left+right by this many bases

	tsCoreLoci AcceptAllLoci;		// used when locating first and next loci overlaps
	
	CFilterRefIDs *pFilterRefIDs;	// used when filtering by RefID
} tsProcParams; 

const int cMaxExcludeHistory = 100;
typedef struct TAG_sExcludeEl {
	struct TAG_sExcludeEl *pNext;
	struct TAG_sExcludeEl *pPrev;
	int SpeciesID;					// identifies species
	int ChromID;					// identifies chromosome
	bool bExclude;					// true if to be excluded, false if not
	} tsExcludeEl;

int gNumExcludeEls = 0;				// current number of elements in gExcludeChroms
tsExcludeEl *gpMRA = NULL;			// pts to most recently accessed or added
tsExcludeEl *gpLRA = NULL;			// pts to least recently accessed
tsExcludeEl gExcludeChroms[cMaxExcludeHistory];


#ifdef _WIN32
// required by str library
#if !defined(__AFX_H__)  ||  defined(STR_NO_WINSTUFF)
HANDLE STR_get_stringres()
{
	return NULL;					//Works for EXEs; in a DLL, return the instance handle
}
#endif

const STRCHAR* STR_get_debugname()
{
	return _T("Loci2Phylip");
}
// end of str library required code
#endif

int
Process(etProcMode ProcMode,					// processing mode

				char *pszExeFile,				// process with this phylip program
				char *pszExeParFile,			// using this parameter file
				char *pzExeInFile,				// use this file to pass alignments into phylip program
				char *pszExeRsltFile,		// capture phylip stdout into this file

				etLociFileType LociFileType,	// expected loci file type
				 char *pszInputFile,			// core loci (CSV or Biobed .bsb) file to process
		 		char *pszFilterRefIDFile,		// exclude any RefIDs in this filter file
 				 char *pszMAFFile,				// multiple alignment file
					char *pszOutputFile,		// where to write out Phylip file or final results file
					char *pszRegionFile,		// region charaterisation biobed file
					int NumIncludeFiles,		// number of include region files
					char **ppszIncludeFiles,	// biobed files containing regions to include - default is to exclude none
					int NumExcludeFiles,		// number of exclude region files
					char **ppszExcludeFiles,	// biobed file containing regions to exclude - default is to include all
					int MinNumSpecies,			// minimum number of of species to be in alignment
					char *pszSpeciesList,		// space or comma separated list of species, priority ordered
					int FlankLen,				// extend the core loci left + right by this many bases
					int	MinCoreLen,				// minimum core length required
					int	MaxCoreLen,				// maximum core length required
					etBEDRegion Region,			// if biobed file then which regions 0=ANY,1=Exons,2=Introns,3=CDS,4=UTRs,5=5'UTR,6=3'UTR
					int NumIncludeChroms,		// number of chromosomes explicitly defined to be included
					char **ppszIncludeChroms,	// ptr to array of reg expressions defining chroms to include - overides exclude
					int NumExcludeChroms,		// number of chromosomes explicitly defined to be excluded
					char **ppszExcludeChroms);	// ptr to array of reg expressions defining chroms to include

char *Region2Txt(etBEDRegion Region);
char *ProcMode2Txt(etProcMode ProcMode);
int ParseNumSpecies(char *pszSpeciesList,tsProcParams *pProcParams);
int TrimQuotes(char *pszTxt);
bool IncludeFilter(int RefChromID,int SubRefOfs,int SubRefEndOfs,tsProcParams *pProcParams);
bool CleanupResources(tsProcParams *pProcParams);
CBEDfile *OpenBedfile(char *pToOpen,bool bReportErrs);
bool ExcludeThisChrom(CMAlignFile *pAlignments,int SpeciesID,int ChromID,tsProcParams *pProcParams);
tsExcludeEl *LocateExclude(int SpeciesID,int ChromID);
bool AddExcludeHistory(int SpeciesID,int ChromID,bool bExclude);
int AddLoci(int RefID,char *pszChrom,char Strand,int StartLoci,int Len,tsProcParams *pParams);
int AddChrom(char *pszChrom,tsProcParams *pParams);
int LoadCSVloci(tsProcParams *pParams);
int LoadBEDloci(tsProcParams *pParams);
int ProcessLoci(tsProcParams *pParams);

int StartAlignBlock(tsProcParams *pProcParams);
int ResetAlignBlock(tsProcParams *pProcParams);
bool OutputAlignColumn(int SeqIdx,tsProcParams *pProcParams);
int EndAlignBlock(int RefID,char *pszRefChrom,int StartLoci,int EndLoci,tsProcParams *pProcParams);
int ProcessAlignedBlock(int RefID,int CoreLen,tsProcParams *pProcParams);

int ParsePhylipRslts(int RefID,int CoreLen,tsProcParams *pProcParams);

int 
ProcessAlignments(char *pszMAF,			 // source bioseq multialignment file
				  tsProcParams *pProcParams); // processing parameters
int NormaliseInDelColumns(tsProcParams *pProcParams,int AlignLen);
int										// returned blockid to next start loading from
LoadContiguousBlocks(int RefSpeciesID,	// reference species identifier
 			   int  BlockID,			// which block to initially start loading from
			   bool *pbLoaded,			// returned indicator as to if any loaded blocks meet processing requirements
			   int *pRefChromID,		// returned reference chromosome identifier 
			   char *pRefStrand,		// returned reference strand
			   int *pRefAlignLen,		// returned alignment (incl InDels) length
   			   int *pRefChromOfs,		// returned alignment start offset
   			   int *pSpeciesIDs,		// input - species of interest identifier array
			   CMAlignFile *pAlignments,
			   tsProcParams *pProcParams);
bool 
ProcAlignBlock(int RefChromID,	    // reference chromosome
			   int RefChromOfs,		// offset along reference chromosome 
			   int AlignLen,		// alignment length incl any InDels
			   tsProcParams *pProcParams); // global processing parameters

tsCoreLoci *GetFirstLociOverlaps(char *pszChrom,	// reference species chromosome
   				int RefChromID,			// chrom identifier
				char RefStrand,			// strand
				int RefChromOfs,		// start loci
				int RefChromEndOfs,		// end loci
				tsProcParams *pProcParams); // global processing parameters

tsCoreLoci *GetLociOverlaps(char *pszChrom,	// reference species chromosome
				int RefChromID,			// chrom identifier
				char RefStrand,			// strand
				int RefChromOfs,		// start loci
				int RefChromEndOfs,		// end loci
				tsCoreLoci *pPrevCoreLoci, // previous loci or NULL if none
				tsProcParams *pProcParams); // global processing parameters

UINT16 GenNameHash(char *pszName);

int WritePhylipFile(char *pszOutputFile,tsProcParams *pProcParams);

int RunPhylip(char *pszInFile,char *pszOutFile,tsProcParams *pProcParams);

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

int iProcMode;
int iLociFileType;
int Idx;
char szOutputFile[_MAX_PATH];
char szInputFile[_MAX_PATH];
char szFilterRefIDFile[_MAX_PATH];  // exclude any RefIDs in this filter file

char szMAFFile[_MAX_PATH];
char szRegionFile[_MAX_PATH];

char szExeFile[_MAX_PATH];	
char szExeParFile[_MAX_PATH];	
char szExeInFile[_MAX_PATH];	
char szExeRsltsFile[_MAX_PATH];	

int NumIncludeChroms;
char *pszIncludeChroms[cMaxIncludeChroms];
int NumExcludeChroms;
char *pszExcludeChroms[cMaxExcludeChroms];
int NumIncludeFiles;
char *pszIncludeFiles[cMaxExcludeFiles];
int NumExcludeFiles;
char *pszExcludeFiles[cMaxIncludeFiles];
char szSpeciesList[512];
int iNumSpecies;
int iFlankLen;			// number of bases to extend flanks of cores left+right
int	iMinCoreLen;		// minimum core length required
int iMaxCoreLen;		// maximum core length required
int iMinNumSpecies;		// minimum number of species required before alignment block is processed
int iRegion;			// genomic regions to include: 0:ANY,1:Exons,2:Introns,3:CDS,4:UTRs,5:5'UTR,6:3'UTR (default = ANY)
int LenFileList;

// command line args
struct arg_lit  *help    = arg_lit0("hH","help",                "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_int *ScreenLogLevel=arg_int0("S", "ScreenLogLevel",	"<int>","Level of diagnostics written to screen 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"diagnostics log file");

struct arg_file *InFile = arg_file1("i",NULL,"<file>",				"input loci of interest from .csv or BED file");
struct arg_file *FilterRefIDFile = arg_file0("X",NULL,"<file>",		"filter out any loci (must be CSV file type) with RefIDs in this filter file");
struct arg_int  *LociFileType = arg_int0("l","locitype","<int>",	"loci file type 0:CSV loci, 1:Biobed BED (default: CSV)");
struct arg_file *MAFFile = arg_file1("I",NULL,"<file>",				"input from bioseq multiple alignment file");
struct arg_int  *ProcMode = arg_int0("x","procmode","<int>",		"processing output to file as 0: Concatenated Phylip, 1: Multiple Phylip files (embed # in filename where RefID is to be substituted)");
struct arg_file *OutFile = arg_file1("o",NULL,"<file>",				"output to file");
struct arg_str  *SpeciesList = arg_str1("s","species","<string>",	"species list, ordered by processing priority with Reference species the first");
struct arg_int  *MinNumSpecies = arg_int0("n","minspecies","<int>",	"minimum number of species required in alignment blocks (default = number of species in list)");

struct arg_int  *FlankLen = arg_int0("e","extendflanks","<int>",	"extend cores left+right by this many bases (default = 0)");
struct arg_int  *MinCoreLen = arg_int0("m","mincorelen","<int>",	"only process elements >= this length (default= 50)");
struct arg_int  *MaxCoreLen = arg_int0("M","maxcorelen","<int>",	"only process elements <= this length (default = 100000)");

struct arg_file *RegionFile = arg_file0("b",NULL,"<file>",			"Process regions from this biobed file");
struct arg_int *Region = arg_int0("g","genomicregion","<int>",		"Process regions 0:ALL,1:Intergenic,2:Exons,3:Introns,4:CDS,5:UTRs,6:5'UTR,7:3'UTR (default = ALL)");

struct arg_file *ExcludeFile = arg_filen("E","exclude","<file>",0,cMaxExcludeFiles,	"exclude all regions in biobed file from processing ");
struct arg_file *IncludeFile = arg_filen("I","include","<file>",0,cMaxExcludeFiles,	"include all regions (unless specific regions excluded) in biobed file");
struct arg_str  *IncludeChroms = arg_strn("z","chromeinclude","<string>",0,cMaxIncludeChroms,"low priority - regular expressions defining species.chromosomes to include for processing");
struct arg_str  *ExcludeChroms = arg_strn("Z","chromexclude","<string>",0,cMaxExcludeChroms,"high priority - regular expressions defining species.chromosomes to exclude from processing");

struct arg_file *ExeFile    = arg_file0("y","exefile","<file>",		"pass blocks to this Phylip program for processing");
struct arg_file *ExeParFile = arg_file0("Y","exepar","<file>",		"file containing phylip processing parameters");
struct arg_file *ExeInFile = arg_file0("r","exein","<file>",		"file containing phylip alignment to pass to Phylip program");
struct arg_file *ExeRsltsFile = arg_file0("R","exeout","<file>",	"file to capture phylip stdout into");

struct arg_end *end = arg_end(20);

void *argtable[] = {help,version,FileLogLevel,ScreenLogLevel,LogFile,
					InFile,FilterRefIDFile,MAFFile,OutFile,
					SpeciesList,MinNumSpecies,
					FlankLen,MinCoreLen,MaxCoreLen,
					ProcMode,LociFileType,Region,RegionFile,ExcludeFile,IncludeFile,IncludeChroms,ExcludeChroms,ExeFile,ExeParFile,ExeInFile,ExeRsltsFile,
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
		printf("\n%s: Version: %d.%2.2d\n",gszProcName,cProgVer/100,cProgVer%100);
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

	iProcMode = ProcMode->count ? ProcMode->ival[0] : eProcModePhylip;
	if(iProcMode < eProcModePhylip || iProcMode > eProcModeMultiPhylip)
		{
		printf("\nError: Requested processing mode '-x%d' not supported",iProcMode);
		exit(1);
		}

	if(ExeFile->count)
		{
		strncpy(szExeFile,ExeFile->filename[0],sizeof(szExeFile));
		szExeFile[sizeof(szExeFile)-1] = '\0';
		TrimQuotes(szExeFile);

		if(!ExeParFile->count)
			{
			printf("\nError: Requested processing by Phylip with '-y%s' but no parameter file specified with '-Y<parfile>'",szExeFile);
			exit(1);
			}
		strncpy(szExeParFile,ExeParFile->filename[0],sizeof(szExeParFile));
		szExeParFile[sizeof(szExeParFile)-1] = '\0';
		TrimQuotes(szExeParFile);

		if(!ExeInFile->count)
			{
			printf("\nError: Requested processing by Phylip with '%s' but no Pylip input file specified with '-r<infile>'",szExeFile);
			exit(1);
			}
		strncpy(szExeInFile,ExeInFile->filename[0],sizeof(szExeInFile));
		szExeInFile[sizeof(szExeInFile)-1] = '\0';
		TrimQuotes(szExeInFile);

		if(!ExeRsltsFile->count)
			{
			printf("\nError: Requested processing by Phylip with '%s' but no stdout capture results file specified with '-R<rsltsfile>'",szExeFile);
			exit(1);
			}
		strncpy(szExeRsltsFile,ExeRsltsFile->filename[0],sizeof(szExeRsltsFile));
		szExeRsltsFile[sizeof(szExeRsltsFile)-1] = '\0';
		TrimQuotes(szExeRsltsFile);
		}
	else
		{
		szExeFile[0] = '\0';
		szExeParFile[0] = '\0';
		szExeInFile[0] = '\0';
		szExeRsltsFile[0] = '\0';
		}

	iLociFileType = LociFileType->count ? LociFileType->ival[0] : eLociFileTypeCSV;
	if(iLociFileType < eLociFileTypeCSV || iLociFileType > eLociFileTypeNone)
		{
		printf("\nError: Requested loci file type '-l%d' not supported",iLociFileType);
		exit(1);
		}
	else
		if(iLociFileType != eLociFileTypeNone && InFile->count == 0)
			{
			printf("\nError: Input loci file '-i<file>' must be specified");
			exit(1);
			}

	if(FilterRefIDFile->count)
		{
		if(iLociFileType != eLociFileTypeCSV)
			{
			printf("\nError: Filtering out RefID permitted only if loci file type is CSV");
			exit(1);
			}
		strncpy(szFilterRefIDFile,FilterRefIDFile->filename[0],sizeof(szFilterRefIDFile));
		szFilterRefIDFile[sizeof(szFilterRefIDFile)-1] = '\0';
		TrimQuotes(szFilterRefIDFile);
		}
	else
		szFilterRefIDFile[0] = '\0';

	iFlankLen = FlankLen->count ? FlankLen->ival[0] : 0;
	if(iFlankLen < 0)
		{
		printf("\nSpecified flanks '-x%d' is less than minimum, assuming you meant 0",iFlankLen);
		iFlankLen = 0;
		}
	else
		if(iFlankLen > cMaxFlankLen)
		{
		printf("\nSpecified flanks '-x%d' is more than %d, assuming you meant %d",iFlankLen,cMaxFlankLen,cMaxFlankLen);
		iFlankLen = cMaxFlankLen;
		}

	iMinCoreLen = MinCoreLen->count ? MinCoreLen->ival[0] : cDfltCoreLen;
	if(iMinCoreLen < cMinCoreLen)
		{
		printf("\nSpecified minimum core length '-n%d' < %d, assuming you meant %d",iMinCoreLen,cMinCoreLen,cMinCoreLen);
		iMinCoreLen = cMinCoreLen;
		}
	else
		{
		if(iMinCoreLen > cMaxCoreLen)
			{
			printf("\nSpecified minimum core length '-n%d' > %d, assuming you meant %d",iMinCoreLen,cMaxCoreLen,cMaxCoreLen);
			iMinCoreLen = cMaxCoreLen;
			}
		}


	iMaxCoreLen = MaxCoreLen->count ? MaxCoreLen->ival[0] : cMaxCoreLen;
	if(iMaxCoreLen < iMinCoreLen)
		{
		printf("\nSpecified maximum hyper length '-n%d' < %d, assuming you meant %d",iMaxCoreLen,iMinCoreLen,iMinCoreLen);
		iMaxCoreLen = iMinCoreLen;
		}
	else
		{
		if(iMaxCoreLen > cMaxCoreLen)
			{
			printf("\nSpecified maximum core length was '-n%d' > %d, assuming you meant '-n%d'",iMaxCoreLen,cMaxCoreLen,cMaxCoreLen);
			iMaxCoreLen = cMaxCoreLen;
			}
		}

	iRegion = Region->count ? Region->ival[0] : ePhylipRAny;	// default as being any region
	if(iRegion < ePhylipRAny)
		{
		printf("\nSpecified region '-g%d' < 0, assuming you meant ANY region",iRegion);
		iRegion = ePhylipRAny;
		}
	else
		{
		if(iRegion > ePhylip3UTR)
			{
			printf("\nSpecified region '-n%d' > %d, assuming you meant 3'DS",iRegion,ePhylip3UTR);
			iRegion = ePhylip3UTR;
			}
		}
	if(iRegion != ePhylipRAny)
		{
		if(RegionFile->count)
			{
			strncpy(szRegionFile,RegionFile->filename[0],sizeof(szRegionFile));
			szRegionFile[sizeof(szRegionFile)-1] = '\0';
			TrimQuotes(szRegionFile);
			}
		else
			{
			printf("\nError: Specified region '-n%d' (%s) but no region file specified with '-b<file>'",iRegion,Region2Txt((etBEDRegion)iRegion));
			exit(1);
			}
		}
	else
		szRegionFile[0] = '\0';

	if(iLociFileType != eLociFileTypeNone)
		{
		strncpy(szInputFile,InFile->filename[0],sizeof(szInputFile));
		szInputFile[sizeof(szInputFile)-1] = '\0';
		TrimQuotes(szInputFile);
		}
	strncpy(szOutputFile,OutFile->filename[0],sizeof(szOutputFile));
	szOutputFile[sizeof(szOutputFile)-1] = '\0';
	TrimQuotes(szOutputFile);
	strncpy(szMAFFile,MAFFile->filename[0],sizeof(szMAFFile));
	szMAFFile[sizeof(szMAFFile)-1] = '\0';
	TrimQuotes(szMAFFile);
	NumIncludeFiles = IncludeFile->count;
	for(Idx=0;Idx < IncludeFile->count; Idx++)
		{
		LenFileList = (int)strlen(IncludeFile->filename[Idx]);
		pszIncludeFiles[Idx] = new char [LenFileList+1];
		strcpy(pszIncludeFiles[Idx],IncludeFile->filename[Idx]);
		TrimQuotes(pszIncludeFiles[Idx]);
		}
	NumExcludeFiles = ExcludeFile->count;
	for(Idx=0;Idx < ExcludeFile->count; Idx++)
		{
		LenFileList = (int)strlen(ExcludeFile->filename[Idx]);
		pszExcludeFiles[Idx] = new char [LenFileList+1];
		strcpy(pszExcludeFiles[Idx],ExcludeFile->filename[Idx]);
		TrimQuotes(pszExcludeFiles[Idx]);
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

	if(!SpeciesList->count)
		{
		printf("\nError: Species list '-s<specieslist>' is empty\n");
		exit(1);
		}
	strcpy(szSpeciesList,SpeciesList->sval[0]);
	TrimQuotes(szSpeciesList);

	iNumSpecies = ParseNumSpecies(szSpeciesList,NULL);
	if(iNumSpecies < cMAFminSpecies)
		{
		printf("Error: Species list %s",iNumSpecies < 0 ? "is unable to be parsed" : "must contain at least 2 species");
		exit(1);
		}
	else
		if(iNumSpecies > cMaxAlignedSpecies)
			{
			printf("Error: Species list contains more than %d species",cMaxAlignedSpecies);
			exit(1);
			}

	iMinNumSpecies = MinNumSpecies->count ? MinNumSpecies->ival[0] : iNumSpecies;
	if(iMinNumSpecies < 2)
		{
		printf("Minimum number of species '-n%d' < 2, assuming you meant 2",iMinNumSpecies);
		iMinNumSpecies = 2;
		}
	else
		if(iMinNumSpecies > iNumSpecies)
			{
			printf("Minimum number of species '-n%d' > species in list, assuming you meant %d",iMinNumSpecies,iNumSpecies);
			iMinNumSpecies = iNumSpecies;
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
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Processing mode: %s",ProcMode2Txt((etProcMode)iProcMode));
	switch(iLociFileType) {
		case eLociFileTypeCSV:
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Process CSV input file %s",szInputFile);
			break;

		case eLociFileTypeBED:
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Process Biobed BED input file %s",szInputFile);
			break;

		case eLociFileTypeNone:
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"No input loci processing");
			break;

		default:
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Unrecognised input file type",iLociFileType);
			break;
		}

	if(szFilterRefIDFile[0])
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Exclude any RefIDs in this filter file: '%s'",szFilterRefIDFile);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Output file: '%s'",szOutputFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Input (.algn) multialignment file to process: '%s'",szInputFile);

	for(Idx = 0; Idx < NumIncludeFiles; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"biobed file containing regions to include: '%s'",pszIncludeFiles[Idx]); 
	for(Idx = 0; Idx < NumExcludeFiles; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"biobed file containing regions to exclude: '%s'",pszExcludeFiles[Idx]); 
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Number of species: %d",iNumSpecies);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Minimum number species required: %d",iMinNumSpecies);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"species list: '%s'",	szSpeciesList);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"minimum core length required: %d",iMinCoreLen);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"maximum core length required: %d",iMaxCoreLen);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"cores will be extended left+right by %d bases",iFlankLen);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Include Region: %s",Region2Txt((etBEDRegion)iRegion));
	if(iRegion != ePhylipRAny)
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Characterise regions from: '%s'",szRegionFile);

	for(Idx = 0; Idx < NumIncludeChroms; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"reg expressions defining chroms to include: '%s'",pszIncludeChroms[Idx]);
	for(Idx = 0; Idx < NumExcludeChroms; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"reg expressions defining chroms to exclude: '%s'",pszExcludeChroms[Idx]); 

	if(szExeFile[0])
		{
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Phylip processing with: '%s'",szExeFile);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Phylip parameter file: '%s'",szExeParFile);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Phylip input file: '%s'",szExeInFile);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Phylip stdout capture file: '%s'",szExeRsltsFile);
		}

	gStopWatch.Start();
#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	Rslt = Process((etProcMode)iProcMode,	// processing mode
					szExeFile,			// process with this phylip program
					szExeParFile,		// using this parameter file
					szExeInFile,		// use this file to pass alignments into phylip program
					szExeRsltsFile,		// capture into this file
					(etLociFileType)iLociFileType,		// expected input file type
					szInputFile,		// core loci (CSV or Biobed .bsb) file to process
					szFilterRefIDFile, // exclude any RefIDs in this filter file
					szMAFFile,			// multiple alignment file
					szOutputFile,		// where to write out sequences
					szRegionFile,		// region charaterisation biobed file
					NumIncludeFiles,	// number of include region files
					pszIncludeFiles,	// biobed files containing regions to include - default is to exclude none
					NumExcludeFiles,	// number of exclude region files
					pszExcludeFiles,	// biobed file containing regions to exclude - default is to include all
					iMinNumSpecies,		// minimum number of of species to be in alignment
					szSpeciesList,		// space or comma separated list of species, priority ordered
					iFlankLen,			// extend the core loci left + right by this many bases
					iMinCoreLen,		// minimum core length required
					iMaxCoreLen,		// maximum core length required
					(etBEDRegion)iRegion,			// if biobed file then which regions 0=ANY,1=Exons,2=Introns,3=CDS,4=UTRs,5=5'UTR,6=3'UTR
					NumIncludeChroms,	// number of chromosomes explicitly defined to be included
					pszIncludeChroms,	// ptr to array of reg expressions defining chroms to include - overides exclude
					NumExcludeChroms,	// number of chromosomes explicitly defined to be excluded
					pszExcludeChroms);	// ptr to array of reg expressions defining chroms to include

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
ProcMode2Txt(etProcMode ProcMode)
{
switch(ProcMode) {
	case eProcModePhylip:
		return((char *)"Single Concatenated Phylip File");
	case eProcModeMultiPhylip:
		return((char *)"Multiple Phylip Files");
	default:
		break;
	}
return((char *)"Unsupported");
}

char *
Region2Txt(etBEDRegion Region)
{
switch(Region) {
	case ePhylipRAny:		// process any region
		return((char *)"All");

	case ePhylipRIntergenic:	// only process intergenic
		return((char *)"Intergenic");

	case ePhylipRExons:	// only process exons
		return((char *)"EXONS");

	case ePhylipRIntrons:	// only process introns
		return((char *)"INTRONS");

	case ePhylipRCDS:		// only process CDSs
		return((char *)"CDS");

	case ePhylipUTR:		// only process UTRs
		return((char *)"UTR");

	case ePhylip5UTR:		// only process 5'UTRs
		return((char *)"5'UTR");

	case ePhylip3UTR:		// only process 3'UTRs
		return((char *)"3'UTR");

	default:
		break;
	}
return((char *)"Unsupported");
}

int
Process(etProcMode ProcMode,				// processing mode
				char *pszExeFile,			// process with this phylip program
				char *pszExeParFile,		// using this parameter file
				char *pszExeInFile,			// use this file to pass alignments into phylip program
				char *pszExeRsltFile,		// capture into this file
				etLociFileType LociFileType,// expected loci file type
				char *pszInputFile,			// core loci (CSV or Biobed .bsb) file to process
		 		char *pszFilterRefIDFile,	// exclude any RefIDs in this filter file
				char *pszMAFFile,			// multiple alignment file
				char *pszOutputFile,	// where to write out Phylip file
				char *pszRegionFile,	// region characterisation biobed file
				int NumIncludeFiles,	// number of include region files
				char **ppszIncludeFiles,// biobed files containing regions to include - default is to exclude none
				int NumExcludeFiles,	// number of exclude region files
				char **ppszExcludeFiles,// biobed file containing regions to exclude - default is to include all
				int MinNumSpecies,		// minimum number of of species to be in alignment
				char *pszSpeciesList,	// space or comma separated list of species, priority ordered
				int FlankLen,		// extend the core loci left + right by this many bases
				int	MinCoreLen,		// minimum core length required
				int	MaxCoreLen,		// maximum core length required
				etBEDRegion Region,			// which regions 0=ANY,1=Exons,2=Introns,3=CDS,4=UTRs,5=5'UTR,6=3'UTR
				int NumIncludeChroms,	// number of chromosomes explicitly defined to be included
				char **ppszIncludeChroms,	// ptr to array of reg expressions defining chroms to include - overides exclude
				int NumExcludeChroms,	// number of chromosomes explicitly defined to be excluded
				char **ppszExcludeChroms)	// ptr to array of reg expressions defining chroms to include
{
int Rslt;
int Idx;
char szOutputFileTpl[_MAX_PATH];
char szExeRsltFileTpl[_MAX_PATH];
char szExeInFileTpl[_MAX_PATH];

char *pDst;
char Chr;

tsProcParams ProcParams;
char szCSVSpecies[1024];				// to hold comma separated species list

#ifndef _WIN32
int RegErr;				// regular expression parsing error
char szRegErr[128];			// to hold RegErr as textual representation ==> regerror();
#endif

memset(&ProcParams,0,sizeof(tsProcParams));

ProcParams.ProcMode = ProcMode;
ProcParams.LociFileType = LociFileType;
ProcParams.MinCoreLen = MinCoreLen;
ProcParams.MaxCoreLen = MaxCoreLen;
ProcParams.FlankLen = FlankLen;
ProcParams.MinNumSpecies = MinNumSpecies;
ProcParams.FiltRegion = Region;
ProcParams.pszExeFile = pszExeFile;
ProcParams.pszExeInFile = pszExeInFile;
ProcParams.pszExeParFile = pszExeParFile;
ProcParams.pszExeRsltFile = pszExeRsltFile;
ProcParams.pszOutputFile = pszOutputFile;
ProcParams.bWrtHdr = true;		// always generate header in output file
if(LociFileType != eLociFileTypeNone)
	ProcParams.pszLociFile = pszInputFile;
else
	ProcParams.pszLociFile = NULL;

if(ProcMode == eProcModePhylip)
	{
	ProcParams.pszOutputFileTpl = pszOutputFile;
	ProcParams.bWCOutputFile = false;

	ProcParams.pszExeInFileTpl = pszExeInFile;
	ProcParams.bWCExeInFile = false;

	ProcParams.pszExeRsltFileTpl = pszExeRsltFile;
	ProcParams.bWCExeRsltFile = false;
	}
else	// multiple output files or each aligned sequence to be separately processed
	{
	pDst = szOutputFileTpl;
	while(Chr = *pszOutputFile++)
		{
		if(Chr == '#')
			{
			*pDst++ = '%';
			*pDst++ = 'd';
			ProcParams.bWCOutputFile = true;
			}
		else
			*pDst++ = Chr;
		}
	*pDst = '\0';
	ProcParams.pszOutputFileTpl = szOutputFileTpl;

	pDst = szExeRsltFileTpl;
	while(Chr = *pszExeRsltFile++)
		{
		if(Chr == '#')
			{
			*pDst++ = '%';
			*pDst++ = 'd';
			ProcParams.bWCExeRsltFile = true;
			}
		else
			*pDst++ = Chr;
		}
	*pDst = '\0';
	ProcParams.pszExeRsltFileTpl = szExeRsltFileTpl;


	pDst = szExeInFileTpl;
	while(Chr = *pszExeInFile++)
		{
		if(Chr == '#')
			{
			*pDst++ = '%';
			*pDst++ = 'd';
			ProcParams.bWCExeInFile = true;
			}
		else
			*pDst++ = Chr;
		}
	*pDst = '\0';
	ProcParams.pszExeInFileTpl = szExeInFileTpl;
	}

if(Region != ePhylipRAny)
	ProcParams.pszRegionFile = pszRegionFile;
else
	ProcParams.pszRegionFile = NULL;

ProcParams.NumIncludeChroms = NumIncludeChroms;		// number of chromosomes explicitly defined to be included
ProcParams.ppszIncludeChroms = ppszIncludeChroms;	// ptr to array of reg expressions defining chroms to include - overides exclude
ProcParams.NumExcludeChroms = NumExcludeChroms;		// number of chromosomes explicitly defined to be excluded
ProcParams.ppszExcludeChroms = ppszExcludeChroms;	// ptr to array of reg expressions defining chroms to include

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
	CleanupResources(&ProcParams);
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to process include regexpr chrom '%s'",ppszIncludeChroms[Idx]);
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
	CleanupResources(&ProcParams);
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to process exclude regexpr chrom '%s'",ppszExcludeChroms[Idx]);
	return(eBSFerrMem);
	}

#else
for(Idx=0;Idx < NumIncludeChroms;Idx++)
	{

	RegErr=regcomp(&ProcParams.IncludeChromsRE[Idx],ppszIncludeChroms[Idx],REG_EXTENDED | REG_ICASE);	// note case insensitive
	if(RegErr)
		{
		regerror(RegErr,&ProcParams.IncludeChromsRE[Idx],szRegErr,sizeof(szRegErr));
		CleanupResources(&ProcParams);
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to process include chrom '%s' error: %s",ppszIncludeChroms[Idx],szRegErr);
		return(eBSFerrMem);
		}
	}
for(Idx=0;Idx < NumExcludeChroms;Idx++)
	{
	RegErr = regcomp(&ProcParams.ExcludeChromsRE[Idx],ppszExcludeChroms[Idx],REG_EXTENDED | REG_ICASE);	// note case insensitive
	if(RegErr)
		{
		regerror(RegErr,&ProcParams.ExcludeChromsRE[Idx],szRegErr,sizeof(szRegErr));
		CleanupResources(&ProcParams);
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to process exclude chrom '%s' error: %s",ppszExcludeChroms[Idx],szRegErr);
		return(eBSFerrMem);
		}
	}
#endif


// parse out species list
ProcParams.MaxNumSpecies = ParseNumSpecies(pszSpeciesList,&ProcParams);
szCSVSpecies[0]='\0';
for(Idx = 0; Idx < ProcParams.MaxNumSpecies; Idx++)
	{
	if(Idx > 0)
		strcat(szCSVSpecies,",");
	strcat(szCSVSpecies,ProcParams.szSpecies[Idx]);
	}
ProcParams.pszSpeciesList = szCSVSpecies;

if(Region != ePhylipRAny)
	{
	if((ProcParams.pBEDregions = OpenBedfile(pszRegionFile,true))==NULL)
		{
		CleanupResources(&ProcParams);
		return(eBSFerrObj);
		}
	if(!ProcParams.pBEDregions->ContainsGeneDetail())
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Genomic region %s specified but file '%s' does not contain regions",Region2Txt(Region),pszRegionFile);
		CleanupResources(&ProcParams);
		return(eBSFerrFileType);
		}
	}

for(Idx=0;Idx<NumIncludeFiles; Idx++)
	{
	if((ProcParams.pIncludes[Idx] = OpenBedfile(ppszIncludeFiles[Idx],true))==NULL)
		{
		CleanupResources(&ProcParams);
		return(eBSFerrObj);
		}
	ProcParams.NumIncludes++;
	}

for(Idx=0;Idx<NumExcludeFiles; Idx++)
	{
	if((ProcParams.pExcludes[Idx] = OpenBedfile(ppszExcludeFiles[Idx],true))==NULL)
		{
		CleanupResources(&ProcParams);
		return(eBSFerrObj);
		}
	ProcParams.NumExcludes++;
	}


switch(LociFileType) {
	case eLociFileTypeCSV:
		if((ProcParams.pCSVloci = new CCSVFile())==NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"\nUnable to instantiate CCSVFile");
			CleanupResources(&ProcParams);
			return(eBSFerrObj);
			}

		if((Rslt = ProcParams.pCSVloci->Open(pszInputFile))!= eBSFSuccess)
			{
			while(ProcParams.pCSVloci->NumErrMsgs())
				gDiagnostics.DiagOut(eDLFatal,gszProcName,ProcParams.pCSVloci->GetErrMsg());
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"\nUnable to open loci file '%s'",pszInputFile);
			CleanupResources(&ProcParams);
			return(Rslt);
			}

		if((Rslt = LoadCSVloci(&ProcParams)) < 1)
			{
			CleanupResources(&ProcParams);
			return(Rslt);
			}
		delete ProcParams.pCSVloci;
		ProcParams.pCSVloci = NULL;
		break;

	case eLociFileTypeBED:
		if((ProcParams.pBEDloci = OpenBedfile(pszInputFile,true))==NULL)
			{
			CleanupResources(&ProcParams);
			return(eBSFerrObj);
			}

		if((Rslt = LoadBEDloci(&ProcParams)) < 1)
			{
			CleanupResources(&ProcParams);
			return(Rslt);
			}
		delete ProcParams.pBEDloci;
		ProcParams.pBEDloci = NULL;
		break;

	case eLociFileTypeNone:		// simple, no loci to filter on
		break;
	}

// determine max aligned sequence length for any single species which can be handled
ProcParams.MaxSeqAlignLen = cMAtotMaxSeqAlignLen/ProcParams.MaxNumSpecies;
for(Idx = 0; Idx < ProcParams.MaxNumSpecies; Idx++)
	{
	if((ProcParams.pSeqs[Idx] = new unsigned char [ProcParams.MaxSeqAlignLen])==NULL)
		{
	    CleanupResources(&ProcParams);
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate memory (%d bytes) for holding species sequences",ProcParams.MaxSeqAlignLen);
		return(eBSFerrMem);
		}
	}

// determine max Phylip sequence length for any single species which can be handled
ProcParams.MaxPhylipSeqLen = cMaxTotPhylipSeqLen/ProcParams.MaxNumSpecies;
for(Idx = 0; Idx < ProcParams.MaxNumSpecies; Idx++)
	{
	if((ProcParams.pPhylipSeqs[Idx] = new unsigned char [ProcParams.MaxPhylipSeqLen])==NULL)
		{
	    CleanupResources(&ProcParams);
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate memory (%d bytes) for holding Phylip sequences",ProcParams.MaxPhylipSeqLen);
		return(eBSFerrMem);
		}
	}

if(ProcParams.pszExeFile != NULL && ProcParams.pszExeFile[0] != '\0')
	{
#ifdef _WIN32
	ProcParams.hOutFile = open(ProcParams.pszOutputFile, O_CREATETRUNC);
#else
     if((ProcParams.hOutFile = open64(ProcParams.pszOutputFile,O_RDWR | O_CREAT,S_IREAD | S_IWRITE))!=-1)
		 if(ftruncate(ProcParams.hOutFile,0)){};
#endif

if(ProcParams.hOutFile == -1)					// check if file open succeeded
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create/open output file '%s' - %s",ProcParams.pszOutputFile,strerror(errno));
	return(eBSFerrCreateFile);
	}
	}

Rslt = ProcessAlignments(pszMAFFile,	// source bioseq multialignment file
				  &ProcParams); // processing parameters

if(Rslt >= 0)
	Rslt = ProcessAlignedBlock(0,0,&ProcParams);

CleanupResources(&ProcParams);
return(Rslt);
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
		return(eBSFerrMem);
		}
	if(pParams->pCoreChroms != NULL)
		{
		memcpy(pChrom,pParams->pCoreChroms,sizeof(tsCoreChrom) * pParams->NumCoreChroms);
		pParams->pCoreChroms = pChrom;
		pParams->MaxCoreChroms += cCoreChromAllocChunk;
		}
	else
		{
		pParams->MaxCoreChroms = cCoreChromAllocChunk;
		pParams->NumCoreLocs = 0;
		}
	pParams->pCoreChroms = pChrom;
	}
pChrom = &pParams->pCoreChroms[pParams->NumCoreChroms++];
pChrom->ChromID = pParams->NumCoreChroms;
pChrom->Hash = Hash;
pChrom->pCoreLoci = NULL;
strcpy(pChrom->szChrom,pszChrom);
return(pParams->CachCoreChromID = pChrom->ChromID);
}


int					// returns > 0 if loci added, 0 if not added due to filtering, < 0 if error
AddLoci(int RefID,char *pszChrom,char Strand,int StartLoci,int Len,tsProcParams *pParams)
{
tsCoreLoci *pNxtCore;
int ChromID;

if(Len < pParams->MinCoreLen || Len > pParams->MaxCoreLen)
	return(0);

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
		if(pParams->NumCoreLocs)
			memcpy(pNxtCore,pParams->pCoreLocs,sizeof(tsCoreLoci) * pParams->NumCoreLocs);
		delete pParams->pCoreLocs;
		pParams->pCoreLocs = pNxtCore;
		}
	else
		{
		pParams->MaxCoreLocs = 0;
		pParams->NumCoreLocs = 0;
		}
	pParams->MaxCoreLocs += cCoreLociAllocChunk;
	pParams->pCoreLocs = pNxtCore;
	}
pNxtCore = &pParams->pCoreLocs[pParams->NumCoreLocs++];
pNxtCore->RefID = RefID;
pNxtCore->ChromID = ChromID;
pNxtCore->Strand = Strand;
pNxtCore->StartLoci = StartLoci;
pNxtCore->EndLoci = StartLoci + Len - 1;
return(pParams->NumCoreLocs);
}

int				// returns > 0 if region accepted, 0 if filtered out, < 0 if error
FilterRegionIn(int RefID,char *pszChrom,char Strand,int StartLoci,int EndLoci,tsProcParams *pParams)
{
int Rslt;
int FeatureLen;
int CurFeatureID;
int ChromID;
int NumEls;

int LociLen;

int NumExons;
int CurExon;
int ExonStartOfs;
int ExonEndOfs;

int NumIntrons;
int CurIntron;
int IntronStartOfs;
int IntronEndOfs;

int CDSStart;
int CDSEnd;
NumEls = 0;
CBEDfile *pBED;

FeatureLen = (EndLoci - StartLoci) + 1;
if(pParams->FiltRegion == ePhylipRAny)
	return(AddLoci(RefID,pszChrom,'+',StartLoci,FeatureLen,pParams));

pBED = pParams->pBEDregions;
ChromID = pBED->LocateChromIDbyName(pszChrom);

if(ChromID < 1 || (CurFeatureID = pBED->LocateFeatureIDinRangeOnChrom(ChromID,StartLoci,EndLoci,1)) < 1)
	{
	if(pParams->FiltRegion == ePhylipRIntergenic)
		return(AddLoci(RefID,pszChrom,'+',StartLoci,FeatureLen,pParams) > 0 ? 1 : 0);
	return(0);
	}

NumEls = 0;
switch(pParams->FiltRegion) {
	case ePhylipRIntrons:		// only interested in INTRONS
		NumIntrons = pBED->GetNumIntrons(CurFeatureID);
		for(CurIntron = 1; CurIntron <= NumIntrons; CurIntron++)
			{
			IntronStartOfs = pBED->GetIntronStart(CurFeatureID,CurIntron);
			IntronEndOfs = pBED->GetIntronEnd(CurFeatureID,CurIntron);
			if((Rslt = AddLoci(RefID,pszChrom,'+',StartLoci+IntronStartOfs,(IntronEndOfs - IntronStartOfs) + 1,pParams)) < 0)
				return(Rslt);
			NumEls++;
			}
		break;

	default:
		NumExons = pBED->GetNumExons(CurFeatureID);
		if(pParams->FiltRegion >= ePhylipRCDS)
			{
			CDSStart = StartLoci + pBED->GetCDSStart(CurFeatureID);
			CDSEnd = StartLoci + pBED->GetCDSEnd(CurFeatureID);
			}
	
		for(CurExon = 1; CurExon <= NumExons; CurExon++)
			{
			ExonStartOfs = pBED->GetExonStart(CurFeatureID,CurExon);
			ExonEndOfs = pBED->GetExonEnd(CurFeatureID,CurExon);
			switch(pParams->FiltRegion) {
				case ePhylipRExons:	// only process exons
					break;

				case ePhylipRCDS:		// only process CDSs
					if(ExonStartOfs > CDSEnd || ExonEndOfs < CDSStart)
						continue;
					if(ExonStartOfs < CDSStart)
						ExonStartOfs = CDSStart;
					else
						if(ExonEndOfs > CDSEnd)
							ExonEndOfs = CDSEnd;
					break;


				case ePhylipUTR:		// only process UTRs
					if(ExonStartOfs >= CDSStart && ExonEndOfs <= CDSEnd)
						continue;
					// special handling for instances where 5'UTR + CDS + 3'UTR in a single exon 
					if(ExonStartOfs < CDSStart && ExonEndOfs > CDSEnd)
						{
						if((Rslt = AddLoci(RefID,pszChrom,'+',ExonStartOfs,CDSStart - ExonStartOfs ,pParams)) < 0)
							return(Rslt);
						NumEls++;
						ExonStartOfs = CDSEnd+1;
						break;
						}

					if(ExonStartOfs < CDSStart && ExonEndOfs >= CDSStart)
						ExonEndOfs = CDSStart - 1;
					else
						if(ExonStartOfs <= CDSEnd)
							ExonStartOfs = CDSEnd+1;
					break;


				case ePhylip5UTR:		// only process 5'UTRs
					if(Strand == '+')
						{
						if(ExonStartOfs >= CDSStart)
							continue;
						if(ExonEndOfs >= CDSStart)
							ExonEndOfs = CDSStart - 1;
						}
					else
						{
						if(ExonEndOfs <= CDSEnd)
							continue;
						if(ExonStartOfs <= CDSEnd)
							ExonStartOfs = CDSEnd + 1;
						}

					break;
				
				
				case ePhylip3UTR:		// only process 3'UTRs
					if(Strand == '+')
						{
						if(ExonEndOfs <= CDSEnd)
							continue;
						if(ExonStartOfs <= CDSEnd)
							ExonStartOfs = CDSEnd+1;
						}
					else
						{
						if(ExonStartOfs >= CDSStart)
							continue;
						if(ExonEndOfs >= CDSStart)
							ExonEndOfs = CDSStart-1;
						}
					break;
				}

			if(ExonStartOfs >= 0 &&  ExonEndOfs >= 0)
				{
				LociLen = (ExonEndOfs - ExonStartOfs) + 1;
				if((Rslt = AddLoci(RefID,pszChrom,'+',ExonStartOfs,LociLen,pParams)) < 0)
					return(Rslt);
				NumEls++;
				}
			}
		break;
	}
return(NumEls);
}

int
LoadCSVloci(tsProcParams *pParams)
{
int NumFields;
int Rslt;
int SrcID;
char *pszChrom;
int StartLoci;
int EndLoci;
int Len;
int NumElsRead;		// number of elements before filtering
int NumElsAccepted;	// number of elements accepted after filtering
int NumFiltRefIDs;	// number of elements filtered out because of matching RefID
int NumFiltLen;		// number of elements filtered out because of length
int NumFiltRegion;	// number of elements filtered out because of region

CCSVFile *pCSV = pParams->pCSVloci;

NumElsRead =0;		// number of elements before filtering
NumElsAccepted =0;	// number of elements accepted after filtering
NumFiltRefIDs =0;	// number of elements filtered out because of matching RefID
NumFiltLen=0;		// number of elements filtered out because of length
NumFiltRegion=0;	// number of elements filtered out because of region

while((Rslt=pCSV->NextLine()) > 0)	// onto next line containing fields
	{
	NumFields = pCSV->GetCurFields();
	if(NumFields < 7)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expected at least 7 fields in '%s', GetCurFields() returned '%d'",pParams->pszLociFile,NumFields);
		return(eBSFerrFieldCnt);
		}
	NumElsRead += 1;

	pCSV->GetInt(1,&SrcID);
	if(pParams->pFilterRefIDs != NULL && pParams->pFilterRefIDs->Locate(SrcID))
		{	
		NumFiltRefIDs += 1;
		continue;
		}

	pCSV->GetInt(7,&Len);
	if(Len < pParams->MinCoreLen || Len > pParams->MaxCoreLen)
		{
		NumFiltLen += 1;
		continue;
		}

	pCSV->GetText(4,&pszChrom);
	pCSV->GetInt(5,&StartLoci);
	pCSV->GetInt(6,&EndLoci);

	if(StartLoci > pParams->FlankLen)
		StartLoci -= pParams->FlankLen;
	else
		StartLoci = 0;
	EndLoci += pParams->FlankLen;

	if((Rslt = FilterRegionIn(SrcID,pszChrom,'+',StartLoci,EndLoci,pParams)) < 0)
		return(Rslt);
	if(Rslt == 0)
		NumFiltRegion += 1;
	else
		NumElsAccepted += 1;
	}

if(Rslt >= 0)
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"%s: Read: %d, Accepted: %d, FiltRefIDs: %d, FiltLen: %d, FiltRegion: %d",
		pParams->pszLociFile,NumElsRead,NumElsAccepted,NumFiltRefIDs,NumFiltLen,NumFiltRegion);

return(ProcessLoci(pParams));
}

int
LoadBEDloci(tsProcParams *pParams)
{
int Rslt;
CBEDfile *pBED = pParams->pBEDloci;
int CurFeatureID;
int NumFilteredLoci;
int NumUnfilteredLoci;
char szChrom[cMaxDatasetSpeciesChrom];
int StartLoci;
int EndLoci;
int FeatureLen;
char Strand;

CurFeatureID = 0;
NumFilteredLoci = 0;
NumUnfilteredLoci = 0;
while((CurFeatureID = pBED->GetNextFeatureID(CurFeatureID)) > 0)
	{
	NumUnfilteredLoci += 1;
	pBED->GetFeature(CurFeatureID,NULL,szChrom,&StartLoci,&EndLoci,0,&Strand);
	FeatureLen = (EndLoci - StartLoci) + 1;
	if(FeatureLen < pParams->MinCoreLen || FeatureLen > pParams->MaxCoreLen)
		continue;

	if(StartLoci > pParams->FlankLen)
		StartLoci -= pParams->FlankLen;
	else
		StartLoci = 0;
	EndLoci += pParams->FlankLen;

	if((Rslt = FilterRegionIn(CurFeatureID,szChrom,Strand,StartLoci,EndLoci,pParams)) < 0)
		return(Rslt);
	NumFilteredLoci += 1;
	}

return(ProcessLoci(pParams));
}

int
ProcessLoci(tsProcParams *pParams)
{
int CurChromID;
int ChromIdx;
int LociIdx;
tsCoreLoci *pLoci;
tsCoreChrom *pChrom;

if(pParams->pCoreLocs == NULL || pParams->NumCoreLocs < 1)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"No accepted loci to process in '%s'",pParams->pszLociFile);
	return(-1);
	}

// sort cores by chrom.start.end
if(pParams->NumCoreLocs > 1)
	qsort(pParams->pCoreLocs,pParams->NumCoreLocs,sizeof(tsCoreLoci),CompareCoreEls);

// add a sentenil element as last element (allocation ensures at least 1 free element!) 
pParams->pCoreLocs[pParams->NumCoreLocs].ChromID = -1;

// process core chroms setting ptrs to first core in each chrom
pLoci = pParams->pCoreLocs;
pChrom = pParams->pCoreChroms;
CurChromID = -1;
for(ChromIdx = LociIdx = 0; LociIdx < pParams->NumCoreLocs; LociIdx++,pLoci++)
	{
	if(pLoci->ChromID != CurChromID)
		{
		pChrom->pCoreLoci = pLoci;
		CurChromID = pLoci->ChromID;
		pChrom += 1;
		}
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Parsed and accepted %d loci for processing from '%s'",LociIdx,pParams->pszLociFile);
return(LociIdx);
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





// AddExcludeHistory
// Adds a tsExcludeEl to the cached history
// The newly added element will be the MRA
bool
AddExcludeHistory(int SpeciesID,int ChromID,bool bExclude)
{
tsExcludeEl *pEl;
if(gNumExcludeEls < cMaxExcludeHistory)
	pEl = &gExcludeChroms[gNumExcludeEls++];
else
	{
	pEl = gpLRA;		// reuse the least recently accessed element
	gpLRA = pEl->pPrev;	// 2nd LRA now becomes the LRA
	gpLRA->pNext = NULL;
	}
if(gpMRA != NULL)
	gpMRA->pPrev = pEl;
pEl->pNext = gpMRA;
pEl->pPrev = NULL;
gpMRA = pEl;
if(gpLRA == NULL)
	gpLRA = pEl;
pEl->bExclude = bExclude;
pEl->SpeciesID = SpeciesID;
pEl->ChromID = ChromID;
return(bExclude);
}

// LocateExclude
// Locates - starting from the MRA - a tsExcludeEl which matches on SpeciesID and ChromID
// If matches then this tsExcludeEl is made the MRA
// Returns ptr to matching tsExcludeEl if located or NULL
tsExcludeEl *
LocateExclude(int SpeciesID,int ChromID)
{
tsExcludeEl *pEl;
pEl = gpMRA;
while(pEl != NULL)
	{
	if(pEl->SpeciesID == SpeciesID && pEl->ChromID == ChromID)
		{
		if(gNumExcludeEls ==1 || gpMRA == pEl)	// if only, or already the MRA then no need for any relinking 
			return(pEl);

		if(gpLRA == pEl)						// if was the LRA then the 2nd LRA becomes the LRA
			{
			gpLRA = pEl->pPrev;
			gpLRA->pNext = NULL;
			}
		else									// not the LRA, and not the MRA
			{
			pEl->pPrev->pNext = pEl->pNext;
			pEl->pNext->pPrev = pEl->pPrev;
			}
		gpMRA->pPrev = pEl;
		pEl->pNext = gpMRA;
		pEl->pPrev = NULL;
		gpMRA = pEl;
		return(pEl);
		}
	pEl = pEl->pNext;
	}
return(NULL);
}


// ExcludeThisChrom
// Returns true if SpeciesID.ChromosomeID is to be excluded from processing
// ExcludeThisChrom
bool
ExcludeThisChrom(CMAlignFile *pAlignments,int SpeciesID,int ChromID,tsProcParams *pProcParams)
{
#ifdef _WIN32
RegexpMatch mc;
#else
regmatch_t mc;
#endif

char *pszSpecies;
char *pszChrom;
char szSpeciesChrom[200];
int Idx;
if(!pProcParams->NumExcludeChroms && !pProcParams->NumIncludeChroms)
	return(false);
// check if this species and chromosome are already known to be included/excluded
tsExcludeEl *pEl;
if((pEl = LocateExclude(SpeciesID,ChromID))!=NULL)
	return(pEl->bExclude);
// haven't seen this species or chromosome before - or else they have been discarded from history...
pszSpecies = pAlignments->GetSpeciesName(SpeciesID);
pszChrom = pAlignments->GetChromName(ChromID);
sprintf(szSpeciesChrom,"%s.%s",pszSpecies,pszChrom);

// to be excluded?
for(Idx = 0; Idx < pProcParams->NumExcludeChroms; Idx++)
#ifdef _WIN32	
	if(pProcParams->ExcludeChromsRE[Idx]->Match(szSpeciesChrom,&mc))
#else
	if(!regexec(&pProcParams->ExcludeChromsRE[Idx],szSpeciesChrom,1,&mc,0))
#endif
		return(AddExcludeHistory(SpeciesID,ChromID,true));

// to be included?
for(Idx = 0; Idx < pProcParams->NumIncludeChroms; Idx++)
	{
#ifdef _WIN32
	if(pProcParams->IncludeChromsRE[Idx]->Match(szSpeciesChrom,&mc))
#else
	if(!regexec(&pProcParams->IncludeChromsRE[Idx],szSpeciesChrom,1,&mc,0))
#endif
		return(AddExcludeHistory(SpeciesID,ChromID,false));
	}


return(AddExcludeHistory(SpeciesID,ChromID,pProcParams->NumIncludeChroms > 0 ? true : false));
}

// OpenBedfile
// Attempts to open specified bedfile
// Returns ptr to opened bedfile or NULL
CBEDfile *
OpenBedfile(char *pToOpen,bool bReportErrs)
{
int Rslt;
CBEDfile *pBed;
if(pToOpen != NULL && pToOpen[0] != '\0')
	{
	if((pBed = (CBEDfile *)new CBEDfile())==NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"\nUnable to instantiate CBEDfile '%s'",pToOpen);
		return(NULL);
		}

	if((Rslt = pBed->Open(pToOpen,eBTAnyBed))!=eBSFSuccess)
		{
		if(bReportErrs)
			{
			while(pBed->NumErrMsgs())
				gDiagnostics.DiagOut(eDLFatal,gszProcName,pBed->GetErrMsg());
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open biobed file %s",pToOpen);
			}
		delete pBed;
		return(NULL);
		}
	return(pBed);
	}
return(NULL);
}

//CleanupResources
//Closes and deletes all opened resources including files
bool
CleanupResources(tsProcParams *pProcParams)
{
int Idx;

if(pProcParams->hOutFile > 0)
	{
	close(pProcParams->hOutFile);
	pProcParams->hOutFile = 0;
	}

for(Idx=0; Idx < pProcParams->NumIncludes; Idx++)
	{
	if(pProcParams->pIncludes[Idx] != NULL)
		delete pProcParams->pIncludes[Idx];
	pProcParams->pIncludes[Idx] = NULL;
	}
for(Idx=0; Idx < pProcParams->NumExcludes; Idx++)
	{
	if(pProcParams->pExcludes[Idx] != NULL)
		delete pProcParams->pExcludes[Idx];
	pProcParams->pExcludes[Idx] = NULL;
	}


#ifdef _WIN32
for(Idx=0;Idx < pProcParams->NumIncludeChroms;Idx++)
	if(pProcParams->IncludeChromsRE[Idx] != NULL)
		{
		delete pProcParams->IncludeChromsRE[Idx];
		pProcParams->IncludeChromsRE[Idx] = NULL;
		}
for(Idx=0;Idx < pProcParams->NumExcludeChroms;Idx++)
	if(pProcParams->ExcludeChromsRE[Idx] != NULL)
		{
		delete pProcParams->ExcludeChromsRE[Idx];
		pProcParams->ExcludeChromsRE[Idx] = NULL;
		}
pProcParams->NumIncludeChroms = 0;
pProcParams->NumExcludeChroms = 0;
#endif

if(pProcParams->pFilterRefIDs != NULL)
	delete pProcParams->pFilterRefIDs;
pProcParams->pFilterRefIDs = NULL;

if(pProcParams->pBEDregions != NULL)
	delete pProcParams->pBEDregions;
pProcParams->pBEDregions = NULL;

if(pProcParams->pBEDloci != NULL)
	delete pProcParams->pBEDloci;
pProcParams->pBEDloci = NULL;
if(pProcParams->pCSVloci != NULL)
	delete pProcParams->pCSVloci;
pProcParams->pCSVloci = NULL;
if(pProcParams->pCoreLocs != NULL)
	delete pProcParams->pCoreLocs;
pProcParams->pCoreLocs = NULL;
pProcParams->MaxCoreLocs = 0;
pProcParams->NumCoreLocs = 0;
if(pProcParams->pCoreChroms != NULL)
	delete pProcParams->pCoreChroms;
pProcParams->pCoreChroms = NULL;
pProcParams->MaxCoreChroms = 0;
pProcParams->NumCoreChroms = 0;
pProcParams->CachCoreChromID = 0;

for(Idx = 0; Idx < pProcParams->MaxNumSpecies; Idx++)
	{
	if(pProcParams->pSeqs[Idx] != NULL)
		{
		delete pProcParams->pSeqs[Idx];
		pProcParams->pSeqs[Idx] = NULL;
		}
	if(pProcParams->pPhylipSeqs[Idx] != NULL)
		{
		delete pProcParams->pPhylipSeqs[Idx];
		pProcParams->pPhylipSeqs[Idx] = NULL;
		}
	}
pProcParams->PhylipLen = 0;
pProcParams->SeqLen = 0;

return(true);
}

// IncludeFilter
// Returns true if subsequence can be processed or false if subsequence is not to be processed
// To be processed a subsequence -
// A) If no Include BED files specified then subsequence is assumed included unless excluded by following rule B).
//    If at least one Include BED file was specified then if subsequence not in any of the Include files then false is returned
// B) If no Exclude files specified then subsequence is returned as Ok (true) to process. If subsequence not in any of the Exclude files
//    then subsequence is returned as Ok (true) to process, otherwise false is returned
bool
IncludeFilter(int RefChromID,int SubRefOfs,int SubRefEndOfs,tsProcParams *pProcParams)
{
int Idx;
int BEDChromID;
if(pProcParams->NumIncludes)
	{
	for(Idx = 0; Idx < pProcParams->NumIncludes; Idx++)
		{
		if(pProcParams->pIncludes[Idx] == NULL) // should'nt ever be NULL but...
			continue;
		if((BEDChromID = pProcParams->pIncludes[Idx]->LocateChromIDbyName(pProcParams->szRefChrom))<1)
			continue;
		if(pProcParams->pIncludes[Idx]->InAnyFeature(BEDChromID,SubRefOfs,SubRefEndOfs))
			break;
		}
	if(Idx == pProcParams->NumIncludes)
		return(false);
	}

if(pProcParams->NumExcludes)
	{
	for(Idx = 0; Idx < pProcParams->NumExcludes; Idx++)
		{
		if(pProcParams->pExcludes[Idx] == NULL) // should'nt ever be NULL but...
			continue;
		if((BEDChromID = pProcParams->pExcludes[Idx]->LocateChromIDbyName(pProcParams->szRefChrom))<1)
			continue;
		if(pProcParams->pExcludes[Idx]->InAnyFeature(BEDChromID,SubRefOfs,SubRefEndOfs))
			return(false);
		}
	}
return(true);
}


int 
ProcessAlignments(char *pszMAF,			 // source bioseq multialignment file
				  tsProcParams *pProcParams) // processing parameters
{
int RefSpeciesID;
int RefChromID;
int PrevRefChromID;
int PrevDispRefChromID;
char *pszRefChrom;
int RefChromOfs;
int RefChromEndOfs;
char RefStrand;
int SpeciesIDs[cMaxAlignedSpecies];
int Rslt;
int RefAlignLen;
CMAlignFile *pAlignments;
int Idx;
int CurBlockID;
int PrevBlockID;
bool bLoaded;

tsCoreLoci *pLoci;

if((pAlignments = new CMAlignFile())==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create new instance of CMAlignFile");
	return(eBSFerrObj);
	}

if((Rslt=pAlignments->Open(pszMAF))!=eBSFSuccess)
	{
	while(pAlignments->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,pAlignments->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open MAF file '%s'",pszMAF);
	return(Rslt);
	}

// ensure all species are represented in multispecies alignment file plus get their species identifiers
for(Idx = 0; Idx < pProcParams->MaxNumSpecies; Idx++)
	{
	if((Rslt = SpeciesIDs[Idx] = pAlignments->LocateSpeciesID(pProcParams->szSpecies[Idx]))<1)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Species '%s' not represented in %s",pProcParams->szSpecies[Idx],pszMAF);
		Rslt = pAlignments->GetNumSpecies();
		for(Idx=1;Idx<=Rslt;Idx++)
			gDiagnostics.DiagOut(eDLFatal,gszProcName," Represented species: %s",pAlignments->GetSpeciesName(Idx));
		delete pAlignments;
		return(eBSFerrEntry);
		}
	if(!Idx)	// reference species is always the first species in the species name list
		RefSpeciesID = SpeciesIDs[0];
	}


// iterate over reference blocks which are sorted by chrom then offset
CurBlockID = 0;
PrevRefChromID = 0;
PrevDispRefChromID = 0;
pProcParams->RefSpeciesIdx = 0;		// reference sequence will always be 1st
pProcParams->AlignBlockID = 0;
// iterate over all alignment blocks
// if a block could contain a loci of interest then start loading as a contiguous block
while((PrevBlockID = CurBlockID) >= 0 && (CurBlockID =	pAlignments->NxtBlock(CurBlockID)))	// returned blockid to next start loading from
	{
	// block must contain the reference species
	RefChromID  = pAlignments->GetRelChromID(CurBlockID,RefSpeciesID);
	if(RefChromID < 1)	// if no chrom for reference species in block then not interested!
		continue;
		
	if(RefChromID != PrevRefChromID)
		{
		PrevRefChromID = RefChromID;
		pszRefChrom = pAlignments->GetChromName(RefChromID);
		strcpy(pProcParams->szRefChrom,pszRefChrom);
		pProcParams->RefChromID = RefChromID;
		}

	RefStrand   = pAlignments->GetStrand(CurBlockID,RefSpeciesID);
	RefChromOfs = pAlignments->GetRelChromOfs(CurBlockID,RefSpeciesID);
	RefChromEndOfs = pAlignments->GetRelChromEndOfs(CurBlockID,RefSpeciesID);

	// if no loci overlaps block then keep iterating over the blocks
	if((pLoci = GetFirstLociOverlaps(pszRefChrom,RefChromID,RefStrand,RefChromOfs,RefChromEndOfs,pProcParams))==NULL)
		continue;

	// block is of interest, load any contiguous blocks
	if((CurBlockID =						// returned blockid to next start loading from
		LoadContiguousBlocks(RefSpeciesID,	// reference species identifier
 			   PrevBlockID,			// which block to initially start loading from
			   &bLoaded,			// returned indicator as to if any loaded blocks meet processing requirements
			   &RefChromID,			// returned reference chromosome identifier 
			   &RefStrand,			// returned reference strand
			   &RefAlignLen,		// returned alignment (incl InDels) length
   			   &RefChromOfs,		// returned alignment start offset
   			   SpeciesIDs,			// input - species of interest identifier array
			   pAlignments,
			   pProcParams)) > 0 || (CurBlockID == eBSFerrAlignBlk && RefAlignLen > 0))
		{
		if(!bLoaded || RefChromID < 1 || RefAlignLen < pProcParams->MinCoreLen)
			continue;

		pProcParams->RefChromID = RefChromID;

		// here we need to normalise the alignments so that there will be no case of all InDels in
		// any column which can occur if none of the processed species is not the reference species
		if((RefAlignLen = NormaliseInDelColumns(pProcParams,RefAlignLen))< pProcParams->MinCoreLen)
			continue;

		if(RefChromID != PrevRefChromID)
			{
			pszRefChrom = pAlignments->GetChromName(RefChromID);
			strcpy(pProcParams->szRefChrom,pszRefChrom);
			}
		if(PrevDispRefChromID != RefChromID)
			{
			PrevDispRefChromID = RefChromID;
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing chromosome %s",pszRefChrom);
			}

		// now have contiguous block that can be processed
		if(!ProcAlignBlock(RefChromID,RefChromOfs,RefAlignLen,pProcParams))
			break;
		}
	}

delete pAlignments;
return(eBSFSuccess);
}

int										// returned blockid to next start loading from
LoadContiguousBlocks(int RefSpeciesID,	// reference species identifier
 			   int  BlockID,			// which block to initially start loading from
			   bool *pbLoaded,			// returned indicator as to if any loaded blocks meet processing requirements
			   int *pRefChromID,		// returned reference chromosome identifier 
			   char *pRefStrand,		// returned reference strand
			   int *pRefAlignLen,		// returned alignment (incl InDels) length
   			   int *pRefChromOfs,		// returned alignment start offset
   			   int *pSpeciesIDs,		// input - species of interest identifier array
			   CMAlignFile *pAlignments,
			   tsProcParams *pProcParams)
{
int CurBlockID;
int PrvBlockID;
int bFirst;
int CurNumAlignments;
int NumBlockSpecies;
int MaxAlignIdxSpecies;
int Idx;
etSeqBase *pSeq;


int	 CurRefChromID;
char CurRefStrand;
int	 CurRefChromOfs;
int	 CurRefChromEndOfs;
int  CurRefAlignLen;

int	 RefChromID;
char RefStrand;
int	 RefChromOfs;
int	 RefChromEndOfs;
int  RefAlignLen;

int RelSpeciesID;
int RelChromID;

bFirst = true;
CurBlockID = BlockID;
PrvBlockID = CurBlockID;
pProcParams->NumBlockSpecies = 0;
pProcParams->MaxAlignIdxSpecies = 0;
RefAlignLen=   0;
CurRefChromID = 0;
CurRefStrand = '+';
CurRefChromOfs = -1;
CurRefChromEndOfs = -1;
CurRefAlignLen = 0;

while(CurBlockID >= 0 && ((CurBlockID = pAlignments->NxtBlock(CurBlockID)) > 0))
	{
	CurRefChromID  = pAlignments->GetRelChromID(CurBlockID,RefSpeciesID);
	CurRefStrand   = pAlignments->GetStrand(CurBlockID,RefSpeciesID);
	CurRefChromOfs = pAlignments->GetRelChromOfs(CurBlockID,RefSpeciesID);
	CurRefChromEndOfs = pAlignments->GetRelChromEndOfs(CurBlockID,RefSpeciesID);
	CurRefAlignLen = pAlignments->GetAlignLen(CurBlockID,RefSpeciesID);
	
			// terminate with blocks in which there are less species than the minimum number we need!
	if((CurNumAlignments = pAlignments->GetNumSpecies(CurBlockID)) < pProcParams->MinNumSpecies)
		break;
		
	// check if this aligned block would overflow alloc'd memory for holding concatenated alignments
	if((CurRefAlignLen + RefAlignLen) > pProcParams->MaxSeqAlignLen)
		{
		CurBlockID = PrvBlockID; // force re-read of block next time
		break;
		}

	if(bFirst) 
		{
		RefChromID =   CurRefChromID;
		RefStrand  =   CurRefStrand;
		RefChromOfs=   CurRefChromOfs;
		RefChromEndOfs=CurRefChromOfs;
		}
	else
		{
		if(CurRefChromID != RefChromID || 
		   CurRefStrand != RefStrand || 
		   (CurRefChromOfs != (RefChromEndOfs + 1)))
			{
			CurBlockID = PrvBlockID; // force re-read of block next time
			break;
			}
		}

		// iterate over species aligned in current block, species are assumed to be in priority order
	for(NumBlockSpecies = Idx = 0; Idx < pProcParams->MaxNumSpecies; Idx++) 
		{
		RelSpeciesID = pSpeciesIDs[Idx];
		RelChromID = pAlignments->GetRelChromID(CurBlockID,RelSpeciesID);
		if(RelChromID < 1)			// assume < 1 is because species not in alignment
			{
			if(Idx < pProcParams->MinNumSpecies)	// first pParams->MinNumSpecies must always be present
				break;
			memset(&pProcParams->pSeqs[Idx][RefAlignLen],eBaseUndef,CurRefAlignLen);
			continue;		
			}

		if(ExcludeThisChrom(pAlignments,RelSpeciesID,RelChromID,pProcParams))	// should this chromosome be accepted for processing or sloughed?
			{
			if(Idx < pProcParams->MinNumSpecies)	// first pParams->MinNumSpecies must always be present
				break;
			memset(&pProcParams->pSeqs[Idx][RefAlignLen],eBaseUndef,CurRefAlignLen);
			continue;	
			}

			// get each sequence
		if((pSeq = pAlignments->GetSeq(CurBlockID,RelSpeciesID))==NULL) 
			{
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"LoadContiguousBlocks: Unexpected missing sequence in alignments");
			break;
			}
		
		memcpy(&pProcParams->pSeqs[Idx][RefAlignLen],pSeq,CurRefAlignLen);
		NumBlockSpecies++;
		MaxAlignIdxSpecies = Idx + 1;
		}
		
	// check if required minimum number of species required were in multialignment
	if(NumBlockSpecies < pProcParams->MinNumSpecies)
		break;
		
	if(NumBlockSpecies > pProcParams->NumBlockSpecies)
		pProcParams->NumBlockSpecies = NumBlockSpecies;
	if(MaxAlignIdxSpecies > pProcParams->MaxAlignIdxSpecies)
		pProcParams->MaxAlignIdxSpecies = MaxAlignIdxSpecies;
	RefAlignLen += CurRefAlignLen;
	RefChromEndOfs=CurRefChromEndOfs;
	PrvBlockID = CurBlockID; 
	bFirst = false;
	}

if(bFirst) 
	{
	*pRefChromID  =   CurRefChromID;
	*pRefStrand   =   CurRefStrand;
	*pRefChromOfs =   CurRefChromOfs;
	*pRefAlignLen =   CurRefAlignLen;
	*pbLoaded = false;
	}
else
	{
	*pRefChromID    =RefChromID;
	*pRefStrand     =RefStrand;
	*pRefChromOfs   =RefChromOfs;
	*pRefAlignLen   =RefAlignLen;
	*pbLoaded = true;
	}
return(CurBlockID);
}



// NormaliseInDelColumns
// Because multialignments may be the result of merged alignments resulting in InDels being generated back into the 
// existing reference sequence then if only a subset of species are being processed there could be all InDels in any column
// of the subset sequences. In addition, where there are optional species then there could be eBaseUndefs in InDel columns
// This function will delete all columns which only contain InDels and/or eBaseUndef
// Returns the subset sequence length after eBaseInDel/eBaseUndef columns have been deleted
int
NormaliseInDelColumns(tsProcParams *pProcParams,int AlignLen)
{
etSeqBase *pSeq;
etSeqBase *pSrcSeq;
etSeqBase SeqBase;
int NumVInDels;
int SeqIdx;
int VIdx;
int FirstInDel=-1;
int NormAlignLen=0;
for(SeqIdx = 0; SeqIdx < AlignLen; SeqIdx++)
	{
	NumVInDels = 0;					// assume no InDels in column
	for(VIdx = 0; VIdx < pProcParams->MaxAlignIdxSpecies; VIdx++)
		{
		pSeq = pProcParams->pSeqs[VIdx];
		SeqBase = pSeq[SeqIdx] & ~cRptMskFlg;
		if(SeqBase == eBaseInDel || SeqBase == eBaseUndef)
			NumVInDels++;
		else
			break;				
		}

	if(NumVInDels == pProcParams->MaxAlignIdxSpecies)	// all were InDels or undefined?
		{
		// mark col for deletion
		for(VIdx = 0; VIdx < pProcParams->MaxAlignIdxSpecies; VIdx++)
			{
			pSeq = pProcParams->pSeqs[VIdx];
			pSeq[SeqIdx] = (etSeqBase)0x0ff;
			}
		if(FirstInDel == -1)		// note idx of first InDel
			FirstInDel = SeqIdx;
		}
	else
		NormAlignLen++;				// accept this column
	}
if(NormAlignLen == AlignLen)	// if no columns to delete
	return(AlignLen);

// have at least one column which is all InDels and is to be deleted
for(VIdx = 0; VIdx < pProcParams->MaxAlignIdxSpecies; VIdx++)
	{
	pSeq = pProcParams->pSeqs[VIdx];
	pSeq = &pSeq[FirstInDel];
	pSrcSeq = pSeq;
	for(SeqIdx = FirstInDel; SeqIdx < AlignLen; SeqIdx++,pSrcSeq++)
		{
		if(*pSrcSeq != (etSeqBase)0x0ff)
			*pSeq++ = *pSrcSeq;
		}
	}
return(NormAlignLen);
}



// ProcAlignBlock
// Process an alignment block which may contain aligned subsequences meeting processing requirements
bool 
ProcAlignBlock(int RefChromID,	    // reference chromosome
			   int RefChromOfs,		// offset along reference chromosome 
			   int AlignLen,		// alignment length incl any InDels
			   tsProcParams *pProcParams) // global processing parameters
{
int SeqIdx;
int AlignLenNoInDels;
etSeqBase *pRefSeq;
etSeqBase RefBase;
int StartOfs;
int EndOfs;
int CurLoci;
int StartLoci;
int EndLoci;
bool bInLoci;
tsCoreLoci *pCoreLoci;

// firstly determine loci range on reference sequence
AlignLenNoInDels = 0;
pRefSeq = pProcParams->pSeqs[0];
StartOfs = -1;
for(SeqIdx = 0; SeqIdx < AlignLen; SeqIdx++,pRefSeq++)
	{
	RefBase = *pRefSeq & ~cRptMskFlg;
	if(RefBase <= eBaseN)
		{
		if(StartOfs < 0)
			StartOfs = SeqIdx;
		AlignLenNoInDels += 1;
		EndOfs = SeqIdx;
		}
	}

pCoreLoci = NULL;
EndLoci = -1;
StartLoci = -1;
while((pCoreLoci = GetLociOverlaps(pProcParams->szRefChrom,RefChromID,'+',RefChromOfs,RefChromOfs + AlignLenNoInDels -1, pCoreLoci, pProcParams))!=NULL)
	{
	pRefSeq = pProcParams->pSeqs[0];
	pRefSeq += StartOfs;
	CurLoci = RefChromOfs;
	bInLoci = false;
	for(SeqIdx = StartOfs; SeqIdx <= EndOfs; SeqIdx++,pRefSeq++)
		{
		if(CurLoci > pCoreLoci->EndLoci)
			break;
		RefBase = *pRefSeq & ~cRptMskFlg;
		if(RefBase <= eBaseN)
			{
			if(CurLoci >= pCoreLoci->StartLoci)
				{
				if(!bInLoci)
					{
					StartLoci = CurLoci;
					EndLoci = CurLoci;
					if(StartAlignBlock(pProcParams) < 0)
						return(false);
					}
				else
					EndLoci += 1;
				bInLoci = true;
				}
			CurLoci += 1;
			}
		else
			{
			if(RefBase != eBaseInDel)
				break;					// can't handle unaligned or missing reference bases (should never occur on ref but who knows)
			}
		if(bInLoci && !OutputAlignColumn(SeqIdx,pProcParams))
			return(false);
		}

	if(bInLoci)
		{
		// check if this block should be filtered out
		if(!IncludeFilter(RefChromID,StartLoci,EndLoci,pProcParams))
			{
			ResetAlignBlock(pProcParams);
			return(true);
			}
		pProcParams->AlignBlockID += 1;
		if(EndAlignBlock(pCoreLoci->RefID,pProcParams->szRefChrom,StartLoci,EndLoci,pProcParams) < 0)
			return(false);
		}
	}

return(true);
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


tsCoreLoci *
GetLociOverlaps(char *pszChrom,			// reference species chromosome
				int RefChromID,			// chrom identifier
				char RefStrand,			// strand
				int RefChromOfs,		// start loci
				int RefChromEndOfs,		// end loci
				tsCoreLoci *pPrevCoreLoci, //
				tsProcParams *pProcParams) // global processing parameters
{
tsCoreLoci *pLoci;
if(pProcParams->LociFileType == eLociFileTypeNone)
	{
	if(pPrevCoreLoci == NULL ||
		pProcParams->AcceptAllLoci.ChromID != RefChromID ||
		pProcParams->AcceptAllLoci.EndLoci != RefChromEndOfs ||
		pProcParams->AcceptAllLoci.StartLoci != RefChromOfs ||
		pProcParams->AcceptAllLoci.Strand != RefStrand)
		{
		pProcParams->AcceptAllLoci.ChromID = RefChromID;
		pProcParams->AcceptAllLoci.EndLoci = RefChromEndOfs;
		pProcParams->AcceptAllLoci.StartLoci = RefChromOfs;
		pProcParams->AcceptAllLoci.Strand = RefStrand;
		return(&pProcParams->AcceptAllLoci);
		}
	return(NULL);
	}

if(pPrevCoreLoci == NULL)
	return(GetFirstLociOverlaps(pszChrom,RefChromID,RefStrand,RefChromOfs,RefChromEndOfs,pProcParams));

pLoci = pPrevCoreLoci + 1;
while(pLoci->ChromID > 0 && pLoci->ChromID == pPrevCoreLoci->ChromID)
	{
	if(pLoci->StartLoci > RefChromEndOfs)
		break;

	if(pLoci->Strand == RefStrand && pLoci->StartLoci <= RefChromEndOfs && pLoci->EndLoci >= RefChromOfs)
		return(pLoci);
	pLoci += 1;
	}
return(NULL);
}


tsCoreLoci *
GetFirstLociOverlaps(char *pszChrom,	// reference species chromosome
   				int RefChromID,			// chrom identifier
				char RefStrand,			// strand
				int RefChromOfs,		// start loci
				int RefChromEndOfs,		// end loci
				tsProcParams *pProcParams) // global processing parameters
{
tsCoreChrom *pChrom;
tsCoreLoci *pLoci;

if(pProcParams->LociFileType == eLociFileTypeNone)
	{
	if(	pProcParams->AcceptAllLoci.ChromID != RefChromID ||
		pProcParams->AcceptAllLoci.EndLoci != RefChromEndOfs ||
		pProcParams->AcceptAllLoci.StartLoci != RefChromOfs ||
		pProcParams->AcceptAllLoci.Strand != RefStrand)
		{
		pProcParams->AcceptAllLoci.ChromID = RefChromID;
		pProcParams->AcceptAllLoci.EndLoci = RefChromEndOfs;
		pProcParams->AcceptAllLoci.StartLoci = RefChromOfs;
		pProcParams->AcceptAllLoci.Strand = RefStrand;
		return(&pProcParams->AcceptAllLoci);
		}
	return(NULL);
	}

if((pChrom = LocateCoreChrom(pszChrom,pProcParams))==NULL)
   return(NULL);
if((pLoci = pChrom->pCoreLoci)==NULL)
	return(NULL);
while(pLoci->ChromID > 0 && pLoci->ChromID == pChrom->ChromID)
	{
	if(pLoci->StartLoci > RefChromEndOfs)
		break;
	if(pLoci->Strand == RefStrand && pLoci->StartLoci <= RefChromEndOfs && pLoci->EndLoci >= RefChromOfs)
		return(pLoci);
	pLoci += 1;
	}
return(NULL);
}


int 
StartAlignBlock(tsProcParams *pProcParams)
{
pProcParams->PhylipBlockStart = pProcParams->PhylipLen;
return(0);
}

int 
ResetAlignBlock(tsProcParams *pProcParams)
{
pProcParams->PhylipLen = pProcParams->PhylipBlockStart;
return(0);
}


typedef struct TAG_sPhylipRslt {
	int	SpeciesID;
	int SpeciesIdx;
	int NumValues;
	float Values[cMAFmaxSpecies];
	char szSpeciesName[80];
} tsPhylipRslt;

int
ParsePhylipRslts(int RefID,int CoreLen,tsProcParams *pProcParams)
{
int hPhylipRslt;
int NumSpecies;
int SpeciesIdx;
int ValueIdx;
tsPhylipRslt Distances[cMAFmaxSpecies];
char szRsltBuffer[cDNADistRsltLen];
char *pBuff;
bool bExpectNL;
int NumChrsProcessed;
int ParseState;
int NumRead;
char Chr;
float Value;
#ifdef _WIN32
hPhylipRslt = open("outfile",_O_RDONLY);
#else
hPhylipRslt = open("outfile",O_RDONLY);
#endif
NumRead = read(hPhylipRslt,szRsltBuffer,sizeof(szRsltBuffer));
szRsltBuffer[NumRead] = '\0';
close(hPhylipRslt);
pBuff = szRsltBuffer;
ParseState = 0;
bExpectNL = false;
ValueIdx = 0;
SpeciesIdx = 0;

while(*pBuff != '\0') {
	if(*pBuff == '\n' || *pBuff == '\r' || *pBuff == '\t' || *pBuff == ' ')
		{
		pBuff += 1;
		continue;
		}
	if(ParseState == 0)
		{
		sscanf(pBuff," %d%n",&NumSpecies,&NumChrsProcessed);
		SpeciesIdx = 0;
		pBuff += NumChrsProcessed;
		ParseState = 1;		// will look for species name
		while((Chr = *pBuff++) != '\0' && Chr != '\n');
		if(Chr == '\0')
			break;
		}

	if(ParseState == 1)
		{
		Distances[SpeciesIdx].SpeciesIdx = SpeciesIdx;
		Distances[SpeciesIdx].SpeciesID = 0;
		Distances[SpeciesIdx].szSpeciesName[0] = '\0';
		Distances[SpeciesIdx].NumValues = 0;
		sscanf(pBuff," %s%n",Distances[SpeciesIdx].szSpeciesName,&NumChrsProcessed);
		pBuff += NumChrsProcessed;
		ValueIdx = 0;
		ParseState = 2;
		continue;
		}
	if(ParseState == 2)
		{
		sscanf(pBuff," %f%n",&Value,&NumChrsProcessed);
		Distances[SpeciesIdx].Values[ValueIdx] = Value;
		Distances[SpeciesIdx].NumValues = ValueIdx;
		pBuff += NumChrsProcessed;
		ValueIdx += 1;
		if(ValueIdx == NumSpecies)
			{
			SpeciesIdx += 1;
			ParseState = 1;
			}
		}
	}

// write out to file here
if(pProcParams->hOutFile != -1)
	{
	if(pProcParams->bWrtHdr)
		{
		NumChrsProcessed = sprintf(szRsltBuffer,"\"RefID\",\"CoreLen\",\"Centric\"");
		for(SpeciesIdx = 0; SpeciesIdx < NumSpecies; SpeciesIdx++)
			NumChrsProcessed += sprintf(&szRsltBuffer[NumChrsProcessed],",\"%s\"",Distances[SpeciesIdx].szSpeciesName);
		NumChrsProcessed += sprintf(&szRsltBuffer[NumChrsProcessed],"\n");
		pProcParams->bWrtHdr = false;
		CUtility::SafeWrite(pProcParams->hOutFile,szRsltBuffer,NumChrsProcessed);
		}

	NumChrsProcessed = 0;
	for(SpeciesIdx = 0; SpeciesIdx < NumSpecies; SpeciesIdx++)
		{
		NumChrsProcessed += sprintf(&szRsltBuffer[NumChrsProcessed],"%d,%d",RefID,CoreLen);
		NumChrsProcessed += sprintf(&szRsltBuffer[NumChrsProcessed],",\"%s\"",Distances[SpeciesIdx].szSpeciesName);
		for(ValueIdx = 0; ValueIdx < NumSpecies; ValueIdx++)
			NumChrsProcessed += sprintf(&szRsltBuffer[NumChrsProcessed],",%.6f",Distances[SpeciesIdx].Values[ValueIdx]);
		NumChrsProcessed += sprintf(&szRsltBuffer[NumChrsProcessed],"\n");
		}
	CUtility::SafeWrite(pProcParams->hOutFile,szRsltBuffer,NumChrsProcessed);
	}
return(0);
}

int
ProcessAlignedBlock(int RefID,int CoreLen, tsProcParams *pProcParams)
{
int Rslt;
char szExeInFile[1024];
char szExeRsltFile[1024];
char *pszExeInFile;
char *pszExeRsltFile;

if(!pProcParams->PhylipLen || 
   pProcParams->PhylipBlockStart == pProcParams->PhylipLen)
   return(0);

if(pProcParams->pszExeRsltFile != NULL && pProcParams->pszExeRsltFile[0] != '\0')
	{
	if(pProcParams->bWCExeInFile)
		{
		sprintf(szExeInFile,pProcParams->pszExeInFileTpl,RefID);
		pszExeInFile = szExeInFile;
		}
	else
		pszExeInFile = pProcParams->pszExeInFile;

	if(pProcParams->bWCExeRsltFile)
		{
		sprintf(szExeRsltFile,pProcParams->pszExeRsltFileTpl,RefID);
		pszExeRsltFile = szExeRsltFile;
		}
	else
		pszExeRsltFile = pProcParams->pszExeRsltFile;

	Rslt = WritePhylipFile(pszExeInFile,pProcParams);
	if(Rslt >= 0)
		Rslt = RunPhylip(pszExeInFile,pszExeRsltFile,pProcParams);
	if(Rslt >= 0)
		{
		Rslt = ParsePhylipRslts(RefID,CoreLen,pProcParams);
		}
	}
else
	{
	if(pProcParams->bWCOutputFile)
		{
		sprintf(szExeRsltFile,pProcParams->pszOutputFileTpl,RefID);
		pszExeRsltFile = szExeRsltFile;
		}
	else
		pszExeRsltFile  = pProcParams->pszOutputFile;

	Rslt = WritePhylipFile(pszExeRsltFile,pProcParams);
	}

pProcParams->PhylipLen = 0;
pProcParams->PhylipBlockStart = 0;
return(Rslt);
}

int 
EndAlignBlock(int RefID,char *pszRefChrom,int StartLoci,int EndLoci,tsProcParams *pProcParams)
{
if(pProcParams->ProcMode == eProcModePhylip || 
   !pProcParams->PhylipLen || 
   pProcParams->PhylipBlockStart == pProcParams->PhylipLen)
   return(0);
return(ProcessAlignedBlock(RefID,1 + EndLoci - StartLoci,pProcParams));
}

bool
OutputAlignColumn(int SeqIdx,tsProcParams *pProcParams)
{
int Idx;
etSeqBase SeqBase;
etSeqBase *pSeq;
unsigned char Base;
unsigned char *pPhylip;

if(pProcParams->PhylipLen >= (pProcParams->MaxPhylipSeqLen-1))
	{
	if(!pProcParams->bPhylipLenTrunc)
		{
		gDiagnostics.DiagOut(eDLWarn,gszProcName,"Truncated concatenated Phylip sequences at %d bases",pProcParams->MaxPhylipSeqLen-1);
		pProcParams->bPhylipLenTrunc = true;
		}
	return(false);
	}

for(Idx = 0; Idx < pProcParams->MaxNumSpecies; Idx++)
	{
	pSeq = pProcParams->pSeqs[Idx];
	SeqBase = pSeq[SeqIdx];
	switch(SeqBase & ~cRptMskFlg) {
		case eBaseA:
			Base = 'a';
			break;
		case eBaseC:
			Base = 'c';
			break;
		case eBaseG:
			Base = 'g';
			break;
		case eBaseT:
			Base = 't';
			break;
		case eBaseInDel:
			Base = '-';
			break;
		default:
			Base = '?';
			break;
		}
	pPhylip = pProcParams->pPhylipSeqs[Idx];
	pPhylip[pProcParams->PhylipLen] = Base;
	pPhylip[pProcParams->PhylipLen+1] = '\0';
	}
pProcParams->PhylipLen += 1;
return(true);
}

int
WritePhylipFile(char *pszPhylipFile,tsProcParams *pProcParams)
{
int LineLen;
int Idx;
char Base;
char *pSeq;
char *pDst;
int LenRemaining;
int SeqLen;
char szLineBuff[1024];
bool bFirst;
int BasesEmitted;
int hOutFile;

if(!pProcParams->PhylipLen)
	return(0);

#ifdef _WIN32
	hOutFile = open(pszPhylipFile, O_CREATETRUNC);
#else
     if((hOutFile = open64(pszPhylipFile,O_RDWR | O_CREAT,S_IREAD | S_IWRITE))!=-1)
		 if(ftruncate(hOutFile,0)){};
#endif

if(hOutFile == -1)					// check if file open succeeded
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create/open output file '%s' - %s",pszPhylipFile,strerror(errno));
	return(eBSFerrCreateFile);
	}

LineLen = sprintf(szLineBuff,"  %d   %d\n",pProcParams->MaxNumSpecies,pProcParams->PhylipLen);
CUtility::SafeWrite(hOutFile,szLineBuff,LineLen);
bFirst = true;
BasesEmitted = 0;
while(BasesEmitted < pProcParams->PhylipLen) 
	{
	for(Idx = 0; Idx < pProcParams->MaxNumSpecies; Idx++)
		{
		if(bFirst)
			LineLen = sprintf(szLineBuff,"%-10.10s",pProcParams->szSpecies[Idx]);
		else
			LineLen = 0;
		pSeq = (char *)pProcParams->pPhylipSeqs[Idx];	// sequence for this species
		pSeq += BasesEmitted;							// where to start from this iteration

		pDst = &szLineBuff[LineLen];
		SeqLen = 0;
		LenRemaining = pProcParams->PhylipLen - BasesEmitted;
		if(bFirst)
			{
			if(LenRemaining > 60)
				LenRemaining = 60;
			}
		else
			if(LenRemaining > 70)
				LenRemaining = 70;
		while(LenRemaining-- && (Base = *pSeq++))
			{
			if(!SeqLen && !((SeqLen) % 10))
				{
				*pDst++ = ' ';
				LineLen++;
				}
			*pDst++ = Base;
			SeqLen++;
			LineLen++;
			}
		LineLen += sprintf(&szLineBuff[LineLen],"\n");
		CUtility::SafeWrite(hOutFile,szLineBuff,LineLen);
		}
	LineLen = sprintf(szLineBuff,"\n");
	CUtility::SafeWrite(hOutFile,szLineBuff,LineLen);
	BasesEmitted += SeqLen;
	bFirst = false;
	}
close(hOutFile);
pProcParams->PhylipLen = 0;
return(1);
}

int RunPhylip(char *pszInFile,char *pszOutFile,tsProcParams *pProcParams)
{
char szExecFileCmd[2000];
int SysRslt;
remove("outfile");
sprintf(szExecFileCmd,"%s < %s >> %s", pProcParams->pszExeFile,pProcParams->pszExeParFile,pszOutFile);
if((SysRslt = system(szExecFileCmd)) < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"RunPhylip Unable to system('%s'), returned %d",szExecFileCmd,SysRslt);
		return(eBSFerrExecFile);
		}
return(eBSFSuccess);
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



