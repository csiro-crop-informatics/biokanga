// filterreads.cpp : Defines the entry point for the console application.
// Purpose -
// Loads aligned short reads generated alignreads (elements) and processes these after user specified filtering into an output file
// for subsequent efficent pipeline processing
//
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

const char *cpszProgVer = "1.1.0";		// increment with each release

const int cMaxNumBEDs = 20;			// allow at most this number of input Bed filter files

const int cMaxAlignChroms = 10000;		// can handle at most this many chromosomes in CSV alignment file

// processing modes
typedef enum TAG_ePMode {		
	ePMdefault,					// default processing mode
	ePMplaceholder				// used to set the enumeration range
	} etPMode;

// output format modes
typedef enum TAG_eFRMode {		
	eFRMdefault,					// default
	eFRMplaceholder					// used to set the enumeration range
	} etFRMode;

typedef enum TAG_eReadsSortMode {
		eRSMReadID,				// index by ascending ReadID
		eRSMHitMatch,			// index by ascending chrom, loci, strand, level
		eRSMplaceholder			// used to limit the enumeration range
} etReadsSortMode;


#pragma pack(1)
typedef struct TAG_sAlignChrom {
	UINT32 ChromID;				// uniquely identifies this chromosome
	char szChromName[cMaxDatasetSpeciesChrom];	// chromosome name
	} tsAlignChrom;


typedef struct TAG_sAlignHit {
	UINT32 AlignHitIdx;			// current read hit index + 1 for this read
	UINT8 State;				// alignment state - 1 if this read is to be filtered out
	UINT32 ReadID;				// read identifier
	UINT32 ChromID;				// identifies chromosome
	UINT8 Strand;				// hit strand - '+' or '-'
	UINT32 Loci;				// offset on chromosome of hit
	UINT16 LociLen;				// length
	int RawLineLen;				// length of original raw line from  CSV or BED file
	size_t RawLineOfs;			// offset into m_pRawLines of the original CSV or BED raw line 
} tsAlignHit;
#pragma pack()


int
Process(etPMode PMode,					// processing mode
		char Strand,					// retain regions on this strand
		etFRMode FMode,					// output format mode
		char *pszInFile,				// file containing alignments
		char *pszFiltInFile,			// where to write filtered in, retained, aligned reads or regions
		char *pszFiltOutFile,			// where to write filtered out aligned reads or regions
		int RegLen,						// 5'US and 3'DS regulatory region length
		int RegionsIncl,				// accept only these regions
		int NumBedFiles,				// number of input BED filtering files containing regions
		char *pszInBedFileSpecs[]);		// input (wildcards allowed) bed region filtering files

int TrimQuotes(char *pszTxt);			// inplace trim leading/trailing single/double quotes
int ParseFeatures(char *pszFeaturesList); // parse features
char *Features2Txt(int Feature);

const int cInitialNumReads = 1000;		// initially alloc to hold this many reads
const size_t cReadsHitAlloc   = (sizeof(tsAlignHit) * cInitialNumReads);// initial allocation 
const size_t cReadsHitReAlloc = (cReadsHitAlloc/2);	// realloc allocation to hold this many reads

const size_t cRawLineAlloc = cInitialNumReads * 40;	// initially alloc mem for holding raw CSV or BED (40 is about the average length line)
const size_t cRawLineReAlloc = (cRawLineAlloc/2);	// then realloc by this size

size_t m_AllocdRawLineMem;					// how much memory is currently allocated for holding raw CSV or BED  lines
size_t m_UsedRawLineMem;					// how much memory is currently used for holding raw CSV  or BED lines
char *m_pRawLines = NULL;					// mem allocated for holding raw CSV or BED lines

CBEDfile *m_pBEDFile = NULL;
CCSVFile *m_pCSVAligns = NULL;
tsAlignHit *m_pAlignHits = NULL;	// memory allocated to hold reads, reads are written contiguously into this memory
UINT32 m_AllocdAlignHitsMem = 0;	// how many bytes of memory  for reads have been allocated
UINT32 m_UsedAlignHitsMem = 0;	// how many bytes of allocated reads memory is currently used 
UINT32 m_NumAlignHits = 0;		// m_pAlignHits contains this many reads
UINT32 m_FinalReadID = 0;		// final read identifier loaded as a preprocessed read (tsProcRead)

tsAlignHit **m_ppAlignHitsIdx = NULL;	// memory allocated to hold array of ptrs to read hits in m_pAlignHits - usually sorted by some critera
UINT32 m_AllocdAlignHitsIdx = 0;		// how many elements for m_pAlignHitsIdx have been allocated
etReadsSortMode	m_CurReadsSortMode;	// sort mode last used on m_ppAlignHitsIdx

int m_NumAlignChroms;							// number of chroms in m_AlignChroms
tsAlignChrom m_AlignChroms[cMaxAlignChroms];	// to hold all chroms from CSV alignment file

tsAlignHit *LocateRead(UINT32 ReadID);			// locate read identified by ReadID

tsAlignHit *IterReads(tsAlignHit *pCurAlignHit);	// iterate over unsorted reads
tsAlignHit *IterSortedReads(tsAlignHit *pCurAlignHit);	// iterate over unsorted reads

int SortAlignHits(etReadsSortMode SortMode);		// index read hits according to specified SortMode

static int SortHitReadIDs(const void *arg1, const void *arg2);
static int SortHitMatch(const void *arg1, const void *arg2);

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
	return _T("kangaf");
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
int Rslt = 0;   			// function result code >= 0 represents success, < 0 on failure

int Idx;

etPMode PMode;				// processing mode
etFRMode FRMode;			// format output mode
int Strand;					// retain reads on this strand
char szFiltInFile[_MAX_PATH];	// write filtered in elements to this file
char szFiltOutFile[_MAX_PATH];  // write filtered out elements to this file
char szAlignsFile[_MAX_PATH];

int FeaturesIn;
char szFeaturesIn[128];		// features to retain
int RegLen;

int NumBedFiles;	// number of input BED filtering files
char *pszInBedFileSpecs[cMaxNumBEDs];  // input (wildcards allowed) bed filtering files


// command line args
struct arg_lit  *help    = arg_lit0("hH","help",                "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to screen and logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"diagnostics log file");

struct arg_int *pmode = arg_int0("m","mode","<int>",		    "processing mode: 0 - default");
struct arg_int *format = arg_int0("M","format","<int>",		    "output format: 0 - (default: 0)");
struct arg_int  *strand = arg_int0("s","strand","<int>",		"Retain reads on this strand: 0 - any, 1 - Watson '+', 2 - Crick '-' (default is any)");
struct arg_file *infile = arg_file1("i","in","<file>",			"input from this aligned reads or region file");
struct arg_file *filtinfile = arg_file0("o","filtinfile","<file>",		"output filtered in, or retained, reads to this file (optional if '-O<filtoutfile>' specified)");
struct arg_file *filtoutfile = arg_file0("O","filtoutfile","<file>",	"output filtered out reads to this file (optional if '-o<filtintfile>' specified)");

struct arg_int  *reglen = arg_int0("L","updnstream","<int>",	"length of 5'up or 3'down  stream regulatory region length (default = 2000) 0..1000000");
struct arg_str  *featuresin = arg_str0("r","regionsin","<string>","Retain any overlaps into any of these annotated features (space or comma delimit), 1: Intergenic, 2: US, 3: 5'UTR, 4: CDS, 5: Intron, 6: 3'UTR, 7: DS");

struct arg_file *bedfiles = arg_filen("I","bedfiles","<file>",0,cMaxNumBEDs, "BED files containing annotated features");

struct arg_end *end = arg_end(20);

void *argtable[] = {help,version,FileLogLevel,LogFile,
					pmode,format,strand,reglen,featuresin,bedfiles,infile,filtinfile,filtoutfile,
					end};

char **pAllArgs;
int argerrors;
argerrors = CUtility::arg_parsefromfile(argc,(char **)argv,&pAllArgs);
if(argerrors >= 0)
	argerrors = arg_parse(argerrors,pAllArgs,argtable);

/* special case: '--help' takes precedence over error reporting */
if (help->count > 0)
        {
		printf("\n%s the K-mer Adaptive Next Generation Aligner Filter, Version %s\nOptions ---\n", gszProcName,cpszProgVer);
        arg_print_syntax(stdout,argtable,"\n");
        arg_print_glossary(stdout,argtable,"  %-25s %s\n");
		printf("\nNote: Parameters can be entered into a parameter file, one parameter per line.");
		printf("\n      To invoke this parameter file then precede it's name with '@'");
		printf("\n      e.g. %s @myparams.txt\n",gszProcName);
		printf("\nPlease report any issues regarding usage of %s at https://github.com/csiro-crop-informatics/biokanga/issues\n\n",gszProcName);
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

	if(filtinfile->count == 0 && filtoutfile->count == 0)
		{
		printf("\nError: Either one or both output files must be specified with '-o<filtinfile>' or '-O<filtoutfile>'");
		exit(1);
		}


	PMode = (etPMode)(pmode->count ? pmode->ival[0] : 0);
	if(PMode < 0 || PMode >= ePMplaceholder)
		{
		printf("\nError: Processing mode '-m%d' specified outside of range %d..%d",PMode,0,(int)ePMplaceholder-1);
		exit(1);
		}

	FRMode = (etFRMode)(format->count ? format->ival[0] : 0);
	if(FRMode < 0 || FRMode >= eFRMplaceholder)
		{
		printf("\nError: Output format mode '-m%d' specified outside of range %d..%d",FRMode,0,(int)eFRMplaceholder-1);
		exit(1);
		}

	Strand = strand->count ? strand->ival[0] : 0;
	if(Strand < 0 || Strand > 2)
		{
		printf("\nError: Strand '-s%d' specified outside of range 0..2",Strand);
		exit(1);
		}
	switch(Strand) {
		case 1: Strand = (int)'+'; break;
		case 2: Strand = (int)'-'; break;
		case 0: Strand = (int)'*'; break;
		}
	
	if(featuresin->count)
		{
		strcpy(szFeaturesIn,featuresin->sval[0]);
		TrimQuotes(szFeaturesIn);
		if((FeaturesIn = ParseFeatures(szFeaturesIn)) < 0)
			{
			printf("Error: unable to parse '-r%s' into features to retain",szFeaturesIn);
			exit(1);
			}
		}
	else
		{
		printf("Error: No features specified to be retained with '-f<features>'");
		exit(1);
		}

	RegLen = reglen->count ? reglen->ival[0] : cDfltRegLen;
	if(RegLen < cMinRegLen)
		{
		printf("\nRegulatory region length '-L%d' less than minimum %d, assuming you meant to use '-L%d'",RegLen,cMinRegLen,cMinRegLen);
		RegLen = cMinRegLen;
		}
	else
		{
		if(RegLen > cMaxRegLen)
			{
			printf("\nRegulatory region length '-L%d' more than maximum %d, assuming you meant to use '-L%d'",RegLen,cMaxRegLen,cMaxRegLen);
			RegLen = cMaxRegLen;
			}
		}

	NumBedFiles = 0;
	if(bedfiles->count)
		{
		for(NumBedFiles=Idx=0;NumBedFiles < cMaxNumBEDs && Idx < bedfiles->count; Idx++)
			{
			pszInBedFileSpecs[Idx] = NULL;
			if(pszInBedFileSpecs[NumBedFiles] == NULL)
				pszInBedFileSpecs[NumBedFiles] = new char [_MAX_PATH];
			strncpy(pszInBedFileSpecs[NumBedFiles],bedfiles->filename[Idx],_MAX_PATH);
			pszInBedFileSpecs[NumBedFiles][_MAX_PATH-1] = '\0';
			CUtility::TrimQuotedWhitespcExtd(pszInBedFileSpecs[NumBedFiles]);
			if(pszInBedFileSpecs[NumBedFiles][0] != '\0')
				NumBedFiles++;
			}
		}
	if(!NumBedFiles)
		{
		printf("\nError: After removal of whitespace, no input BED feature filter file(s) specified with '-I<filespec>' option)");
		exit(1);
		}

	szFiltInFile[0] = '\0';
	szFiltOutFile[0] = '\0';
	if(filtinfile->count)
		strcpy(szFiltInFile,filtinfile->filename[0]);
	if(filtoutfile->count)
		strcpy(szFiltOutFile,filtoutfile->filename[0]);
	strcpy(szAlignsFile,infile->filename[0]);

			// now that command parameters have been parsed then initialise diagnostics log system
	if(!gDiagnostics.Open(szLogFile,(etDiagLevel)iScreenLogLevel,(etDiagLevel)iFileLogLevel,true))
		{
		printf("\nError: Unable to start diagnostics subsystem.");
		if(szLogFile[0] != '\0')
			printf(" Most likely cause is that logfile '%s' can't be opened/created",szLogFile);
		exit(1);
		}

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Version: %s Processing parameters:",cpszProgVer);

	const char *pszDescr;
	switch(PMode) {
		case ePMdefault:				
			pszDescr = "Default";
			break;
		}
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"processing mode is : '%s'",pszDescr);

	switch(FRMode) {
		case eFRMdefault:				
			pszDescr = "Default";
			break;
		}

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"output format is : '%s'",pszDescr);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"process for this strand only: '%c'",(char)Strand);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"input reads or regions alignments file: '%s'",szAlignsFile);

	if(szFiltInFile[0] != '\0')
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"output filtered in, retained, reads/regions to file: '%s'",szFiltInFile);
	if(szFiltOutFile[0] != '\0')
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"output filtered out reads/regions to file: '%s'",szFiltOutFile);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"5'US and 3'DS regulatory region length: %d",RegLen);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Filter in, retain, alignment overlaps onto these features: '%s'",Features2Txt(FeaturesIn));
	for(Idx=0; Idx < NumBedFiles; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Filter on features in these BED files (%d): '%s'",Idx+1,pszInBedFileSpecs[Idx]);

#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	gStopWatch.Start();
	Rslt = Process(PMode,(char)Strand,FRMode,szAlignsFile,szFiltInFile,szFiltOutFile,RegLen,FeaturesIn,NumBedFiles,pszInBedFileSpecs);
	gStopWatch.Stop();
	Rslt = Rslt >=0 ? 0 : 1;
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Exit code: %d Total processing time: %s",Rslt,gStopWatch.Read());
	exit(Rslt);
	}
else
	{
	printf("\n%s the K-mer Adaptive Next Generation Aligner Filter, Version %s\n",gszProcName,cpszProgVer);
	arg_print_errors(stdout,end,gszProcName);
	arg_print_syntax(stdout,argtable,"\nUse '-h' to view option and parameter usage\n");
	exit(1);
	}
return 0;
}



void
Reset(void)
{
if(m_pBEDFile != NULL)
	{
	delete m_pBEDFile;
	m_pBEDFile = NULL;
	}
if(m_pCSVAligns != NULL)
	{
	delete m_pCSVAligns;
	m_pCSVAligns = NULL;
	}

if(m_ppAlignHitsIdx != NULL)
	{
	delete m_ppAlignHitsIdx;
	m_ppAlignHitsIdx = NULL;
	}

if(m_pAlignHits != NULL)
	{
#ifdef _WIN32
	free(m_pAlignHits);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pAlignHits != MAP_FAILED)
		munmap(m_pAlignHits,m_AllocdAlignHitsMem);
#endif
	m_pAlignHits = NULL;
	}

if(m_pRawLines != NULL)
	{
#ifdef _WIN32
	free(m_pRawLines);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pRawLines != MAP_FAILED)
		munmap(m_pAlignHits,m_AllocdAlignHitsMem);
#endif
	m_pRawLines = NULL;
	}

m_AllocdRawLineMem = 0;					// how much memory is currently allocated for holding raw CSV or BED  lines
m_UsedRawLineMem = 0;					// how much memory is currently used for holding raw CSV  or BED lines

m_NumAlignChroms = 0;
m_AllocdAlignHitsMem = 0;
m_UsedAlignHitsMem = 0;
m_NumAlignHits = 0;
m_FinalReadID = 0;
m_AllocdAlignHitsIdx = 0;
}

void
Init(void)
{
m_pBEDFile = NULL;
m_pCSVAligns = NULL;
m_pAlignHits = NULL;
m_ppAlignHitsIdx = NULL;
m_pRawLines = NULL;
Reset();
}

char *
LocateAlignChromName(int ChromID)
{
if(!m_NumAlignChroms || ChromID <= 0 || ChromID > m_NumAlignChroms)
	return(NULL);
return(m_AlignChroms[ChromID-1].szChromName);
}

int
LocateAlignChrom(char *pszChromName)
{
int Idx;
tsAlignChrom *pChrom;
if(!m_NumAlignChroms)
	return(0);
pChrom = &m_AlignChroms[m_NumAlignChroms-1];
for(Idx = 0; Idx < m_NumAlignChroms; Idx++,pChrom--)
	if(!stricmp(pChrom->szChromName,pszChromName))
		return(pChrom->ChromID);
return(-1);
}

int // returned chrom identifier
AddAlignChrom(char *pszChromName)
{
int ChromID;
tsAlignChrom *pChrom;
if((ChromID = LocateAlignChrom(pszChromName))>0)
	return(ChromID);

if(m_NumAlignChroms == cMaxAlignChroms)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Too many different chromosomes in alignment file - can only process at most %d",cMaxAlignChroms);
	Reset();
	return(-1);
	}
pChrom = &m_AlignChroms[m_NumAlignChroms++];
pChrom->ChromID = m_NumAlignChroms;
strcpy(pChrom->szChromName,pszChromName);
return(m_NumAlignChroms);
}


// ParseFeatures
// Parses space or comma delimited list of Features in which
// 1 == Intergenic, 2 == US, 3 == 5'UTR, 4 == CDS, 5 == Intron, 6 == 3'UTR, 7 == DS
//
// Returns bitmap of Features or -1 if parse errors
// If no Features specified then assumes all Features are selected
int
ParseFeatures(char *pszFeaturesList)
{
// parse out Features list
char Chr;
int Features = 0;
if(pszFeaturesList == NULL || *pszFeaturesList == '\0')
	return(cFeatBitIG | cFeatBitUpstream | cFeatBit5UTR | cFeatBitCDS | cFeatBitIntrons | cFeatBit3UTR | cFeatBitDnstream);

while(Chr = *pszFeaturesList++) {
	if(isspace(Chr) || Chr == ',')		// accept spaces and commas as separators
		continue;
	switch(Chr) {
		case '1':						// intergenic to be filtered
			Features |= cFeatBitIG;
			break;
		case '2':						// 5'US to be filtered
			Features |= cFeatBitUpstream;
			break;
		case '3':						// 5'UTR to be filtered
			Features |= cFeatBit5UTR;
			break;
		case '4':
			Features |= cFeatBitCDS;		// CDS to be filtered
			break;
		case '5':
			Features |=  cFeatBitIntrons;	// any intronic to be filtered
			break;
		case '6':
			Features |=  cFeatBit3UTR;	// any 3'UTR to be filtered
			break;
		case '7':
			Features |=  cFeatBitDnstream;	// any 3'DS to be filtered 	
			break;
		default:
			return(-1);
		}
	}
return(Features);
}

// Features2Txt
// Returns textual representation of regions
char *
Features2Txt(int Features)
{
static char szFeatures[200];
if(!Features)
	return((char *)"None specified");
if(Features & cFeatBitIG || Features == 0)
	strcpy(szFeatures,"Intergenic");
else
	szFeatures[0] = '\0';
if(Features & cFeatBitUpstream)
	{
	if(szFeatures[0] != '\0')
		strcat(szFeatures,", ");
	strcat(szFeatures,"5'US");
	}
if(Features & cFeatBit5UTR)
	{
	if(szFeatures[0] != '\0')
		strcat(szFeatures,", ");
	strcat(szFeatures,"5'UTR");
	}
if(Features & cFeatBitCDS)
	{
	if(szFeatures[0] != '\0')
		strcat(szFeatures,", ");
	strcat(szFeatures,"CDS");
	}
if(Features & cFeatBitIntrons)
	{
	if(szFeatures[0] != '\0')
		strcat(szFeatures,", ");
	strcat(szFeatures,"Introns");
	}
if(Features & cFeatBit3UTR)
	{
	if(szFeatures[0] != '\0')
		strcat(szFeatures,", ");
	strcat(szFeatures,"3'UTR");
	}
if(Features & cFeatBitDnstream)
	{
	if(szFeatures[0] != '\0')
		strcat(szFeatures,", ");
	strcat(szFeatures,"3'DS");
	}
return(szFeatures);
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

int
InitialMemAlloc(void)
{
size_t memreq;

if(m_pAlignHits != NULL)
	{
#ifdef _WIN32
	free(m_pAlignHits);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pAlignHits != MAP_FAILED)
		munmap(m_pAlignHits,m_AllocdAlignHitsMem);
#endif
	m_pAlignHits = NULL;
	}

if(m_pRawLines != NULL)
	{
#ifdef _WIN32
	free(m_pRawLines);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pRawLines != MAP_FAILED)
		munmap(m_pAlignHits,m_AllocdAlignHitsMem);
#endif
	m_pRawLines = NULL;
	}

memreq = cReadsHitAlloc;
#ifdef _WIN32
m_pAlignHits = (tsAlignHit *) malloc(memreq);	// initial and perhaps the only allocation
if(m_pAlignHits == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"InitialMemAlloc: Memory allocation of %lld bytes - %s",(INT64)memreq,strerror(errno));
	Reset();
	return(eBSFerrMem);
	}
#else
	// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
m_pAlignHits = (tsAlignHit *)mmap(NULL,memreq, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
if(m_pAlignHits == MAP_FAILED)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"InitialMemAlloc: Memory allocation of %lld bytes through mmap()  failed - %s",(INT64)memreq,strerror(errno));
	m_pAlignHits = NULL;
	Reset();
	return(eBSFerrMem);
	}
#endif
m_AllocdAlignHitsMem = cReadsHitAlloc;
m_UsedAlignHitsMem = 0;
m_NumAlignHits = 0;


memreq = cRawLineAlloc;
#ifdef _WIN32
m_pRawLines = (char *) malloc(memreq);	// initial and perhaps the only allocation
if(m_pRawLines == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"InitialMemAlloc: Memory allocation of %lld bytes - %s",(INT64)memreq,strerror(errno));
	Reset();
	return(eBSFerrMem);
	}
#else
	// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
m_pRawLines = (char *)mmap(NULL,memreq, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
if(m_pRawLines == MAP_FAILED)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"InitialMemAlloc: Memory allocation of %lld bytes through mmap()  failed - %s",(INT64)memreq,strerror(errno));
	m_pRawLines = NULL;
	Reset();
	return(eBSFerrMem);
	}
#endif
m_AllocdRawLineMem = cRawLineAlloc;
m_UsedRawLineMem = 0;
return(eBSFSuccess);
}

int
LoadReads(char Strand,		// retain reads/regions on this strand
		  char *pszInFile)	// load reads from this file
{
int ChromID;
int Rslt;
int NumProcessed;
int ReadID;
int LociLen;
char *pszTargSpecies;
char *pszChromName;
int Loci;
int LociEnd;
char *pszStrand;
int FiltByStrand;
int NumFields;

tsAlignHit *pAlignHit;					// current read hit
int BuffLen = 0;
int BuffOfs = 0;

char *pRawLine;

m_NumAlignHits = 0;

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Loading read alignments from file: %s",pszInFile);




// load into memory
if((m_pCSVAligns = new CCSVFile) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CCSVfile");
	Reset();
	return(eBSFerrObj);
	}


m_pCSVAligns->SetMaxFields(8);	// only interested in first 8 fields in reads alignment CSV file
if((Rslt=m_pCSVAligns->Open(pszInFile))!=eBSFSuccess)
	{
	while(m_pCSVAligns->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,m_pCSVAligns->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open file: %s",pszInFile);
	Reset();
	return(Rslt);
	}
NumProcessed = 0;
FiltByStrand = 0;
while((Rslt=m_pCSVAligns->NextLine()) > 0)	// onto next line containing fields
	{
	NumFields = m_pCSVAligns->GetCurFields();
	if(NumFields < 8)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"file: %s contains % fields, expected at least %d",pszInFile,NumFields,8);
		Reset();
		return(eBSFerrParams);
		}

	if(!NumProcessed && m_pCSVAligns->IsLikelyHeaderLine())
		continue;
	if(!(m_NumAlignHits % 100000))
		{
		if(!NumProcessed)
			printf("\n     processing alignment %8.8d",m_NumAlignHits);
		else
			printf("\b\b\b\b\b\b\b\b%8.8d",m_NumAlignHits);
		}
	NumProcessed += 1;

	m_pCSVAligns->GetInt(1,&ReadID);
	m_pCSVAligns->GetText(3,&pszTargSpecies);
	m_pCSVAligns->GetText(4,&pszChromName);
	m_pCSVAligns->GetInt(5,&Loci);
	m_pCSVAligns->GetInt(6,&LociEnd);
	m_pCSVAligns->GetInt(7,&LociLen);
	m_pCSVAligns->GetText(8,&pszStrand);
	if(Strand != '*' && *pszStrand != Strand)
		{
		FiltByStrand += 1;
		continue;
		}

	if(m_UsedAlignHitsMem + (10 * sizeof(tsAlignHit)) >= m_AllocdAlignHitsMem)
		{
#ifdef _WIN32
		pAlignHit = (tsAlignHit *) realloc(m_pAlignHits,m_AllocdAlignHitsMem + cReadsHitReAlloc);
#else
		pAlignHit = (tsAlignHit *) mremap(m_pAlignHits,m_AllocdAlignHitsMem,m_AllocdAlignHitsMem + cReadsHitReAlloc,MREMAP_MAYMOVE);
		if(pAlignHit == MAP_FAILED)
			pAlignHit = NULL;
#endif
		if(pAlignHit == NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadReads: Memory re-allocation to %d bytes - %s",m_AllocdAlignHitsMem + cReadsHitReAlloc,strerror(errno));
			return(eBSFerrMem);
			}
		m_pAlignHits = pAlignHit;
		m_AllocdAlignHitsMem += cReadsHitReAlloc;
		}


	if(m_UsedRawLineMem + (cReadsHitReAlloc/100) >= m_AllocdRawLineMem)
		{
#ifdef _WIN32
		pRawLine = (char *) realloc(m_pRawLines,m_AllocdRawLineMem + cRawLineReAlloc);
#else
		pRawLine = (char *) mremap(m_pRawLines,m_AllocdRawLineMem,m_AllocdRawLineMem + cRawLineReAlloc,MREMAP_MAYMOVE);
		if(pAlignHit == MAP_FAILED)
			pAlignHit = NULL;
#endif
		if(pRawLine == NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadReads: Memory re-allocation to %d bytes - %s",m_AllocdRawLineMem + cRawLineReAlloc,strerror(errno));
			return(eBSFerrMem);
			}
		m_pRawLines = pRawLine;
		m_AllocdRawLineMem += cRawLineReAlloc;
		}


	
	pAlignHit = (tsAlignHit *)((UINT8 *)m_pAlignHits + m_UsedAlignHitsMem);
	m_UsedAlignHitsMem += sizeof(tsAlignHit); 
	memset(pAlignHit,0,sizeof(tsAlignHit));
	pAlignHit->ReadID = ReadID;
	pAlignHit->LociLen = LociLen;
	pAlignHit->RawLineLen = m_pCSVAligns->GetLine((int)(m_AllocdRawLineMem - m_UsedRawLineMem), &m_pRawLines[m_UsedRawLineMem]);
	pAlignHit->RawLineOfs = m_UsedRawLineMem;
	m_UsedRawLineMem += pAlignHit->RawLineLen + 1; 
	ChromID = AddAlignChrom(pszChromName);
	if(ChromID < 1)
		return(ChromID);
	pAlignHit->ChromID = ChromID;
	pAlignHit->Loci = Loci;
	pAlignHit->Strand = *pszStrand;
	m_FinalReadID = ReadID;
	m_NumAlignHits += 1;
	}
delete m_pCSVAligns;
m_pCSVAligns = NULL;

printf("\b\b\b\b\b\b\b\b%8.8d",m_NumAlignHits);

gDiagnostics.DiagOut(eDLInfo,gszProcName,"File: %s contained %d elements of which %d were filtered out because of strand",pszInFile,m_NumAlignHits+FiltByStrand,FiltByStrand);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"There are %d elements remaining to be filtered by genomic feature",m_NumAlignHits);

// finally, create sorted index by chrom, loci, strand over loaded reads
Rslt = SortAlignHits(eRSMHitMatch);
return(Rslt);
}

int BEDFilter(char ROIStrand,				// filter for this strand of interest
	  		int RegLen,						// 5'US and 3'DS regulatory region length
			 int FeatureIncl,				// reads overlaying these features are to be retained
			 char *pszBedFile)				// UCSC BED input file containing features of interest
{
int Rslt;
UINT32 NumEls;
UINT32 NumElsOut;
int OverlapFeats;
UINT32 BEDChromID;
UINT32 ReadsChromID;
char *pszChrom;
tsAlignHit *pCurAlignHit;

if((m_pBEDFile = new CBEDfile)==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CBEDfile");
	Reset();
	return(eBSFerrObj);
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Starting to load %s",pszBedFile);
if((Rslt=m_pBEDFile->Open(pszBedFile,eBTAnyBed)) !=eBSFSuccess)
	{
	while(m_pBEDFile->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,m_pBEDFile->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open '%s' for processing",pszBedFile);
	Reset();
	return(eBSFerrOpnFile);
	}
if(!m_pBEDFile->ContainsGeneDetail())			// returns true if file contains gene detail (utr/cds/intron etc)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"bed file %s does not contain gene features",pszBedFile);
	Reset();
	return(eBSFerrFeature);
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"BED file '%s' Loaded, iterating elements...",pszBedFile);

m_pBEDFile->SetStrand(ROIStrand);

// now iterate over the reads
BEDChromID = 0;
ReadsChromID = 0;
pCurAlignHit = NULL;
NumEls = 0;
NumElsOut = 0;
while((pCurAlignHit = IterSortedReads(pCurAlignHit))!=NULL)
	{
	if(pCurAlignHit->State & 0x01)			// already filtered out?
		continue;
	NumEls += 1;
	if(pCurAlignHit->ChromID != ReadsChromID)
		{
		pszChrom = LocateAlignChromName(pCurAlignHit->ChromID);
		if((BEDChromID = m_pBEDFile->LocateChromIDbyName(pszChrom)) < 1)	// if chrom not known in BED then filter out this read unless intergenic features requested
			{
			if(!(FeatureIncl & cFeatBitIG))
				{
				pCurAlignHit->State |= 0x01;
				NumElsOut += 1;
				}
			continue;
			}
		ReadsChromID = pCurAlignHit->ChromID;
		}
	OverlapFeats = m_pBEDFile->GetFeatureBits(BEDChromID,pCurAlignHit->Loci,pCurAlignHit->Loci+pCurAlignHit->LociLen-1,cRegionFeatBits,RegLen);
	if(OverlapFeats == 0)
		OverlapFeats |= cFeatBitIG;
	if(!(FeatureIncl & OverlapFeats))
		{
		pCurAlignHit->State |= 0x01;
		NumElsOut += 1;
		}
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Filtering of remaining %d reads/regions resulted in %d filtered out",NumEls,NumElsOut);
m_pBEDFile->Close();
return(Rslt < 0 ? Rslt : NumEls - NumElsOut);
}


int
WriteResults(char *pszFiltInFile,			// where to write filtered in, retained, aligned reads or regions
			char *pszFiltOutFile)			// where to write filtered out aligned reads or regions
{
int Rslt;
int hFiltOutFile;
int hFiltInFile;
tsAlignHit *pCurAlignHit;
char szLineBufIn[16000];
int BufIdxIn;
char szLineBufOut[16000];
int BufIdxOut;

int NumFilteredIn;
int NumFilteredOut;

hFiltOutFile = -1;
hFiltInFile = -1;
if(pszFiltInFile != NULL && pszFiltInFile[0] != '\0')
	{
#ifdef _WIN32
	if((hFiltInFile = open(pszFiltInFile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
	if((hFiltInFile = open(pszFiltInFile, O_RDWR | O_CREAT | O_TRUNC, S_IREAD | S_IWRITE ))==-1)
#endif
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create/truncate output filtered in, retained, file: %s - %s",pszFiltInFile,strerror(errno));
		return(eBSFerrCreateFile);
		}
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Writing filtered in, retained, reads/regions to results file '%s'",pszFiltInFile);
	}

if(pszFiltOutFile != NULL && pszFiltOutFile[0] != '\0')
	{
#ifdef _WIN32
	if((hFiltOutFile = open(pszFiltOutFile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
	if((hFiltOutFile = open(pszFiltOutFile, O_RDWR | O_CREAT | O_TRUNC, S_IREAD | S_IWRITE ))==-1)
#endif
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create/truncate output filtered in, retained, file: %s - %s",pszFiltOutFile,strerror(errno));
		if(hFiltInFile != -1)
			close(hFiltInFile);
		return(eBSFerrCreateFile);
		}
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Writing filtered out reads/regions to results file '%s'",pszFiltOutFile);
	}

Rslt = eBSFSuccess;
BufIdxIn = 0;
BufIdxOut = 0;
NumFilteredOut = 0;
NumFilteredIn = 0;
pCurAlignHit = NULL;
while((pCurAlignHit = IterSortedReads(pCurAlignHit))!=NULL)
	{
	if(pCurAlignHit->State & 0x01)			// filtered out?
		{
		NumFilteredOut += 1;
		if(hFiltOutFile == -1)
			continue;
		BufIdxOut += sprintf(&szLineBufOut[BufIdxOut],"%s\n",&m_pRawLines[pCurAlignHit->RawLineOfs]);
		if((BufIdxOut+500) > sizeof(szLineBufOut))
			{
			CUtility::SafeWrite(hFiltOutFile,szLineBufOut,BufIdxOut);
			BufIdxOut = 0;
			}
		}
	else
		{
		NumFilteredIn += 1;
		if(hFiltInFile == -1)
			continue;
		BufIdxIn += sprintf(&szLineBufIn[BufIdxIn],"%s\n",&m_pRawLines[pCurAlignHit->RawLineOfs]);
		if((BufIdxIn+500) > sizeof(szLineBufIn))
			{
			CUtility::SafeWrite(hFiltInFile,szLineBufIn,BufIdxIn);
			BufIdxIn = 0;
			}
		}
	}

if(hFiltInFile != -1)
	{
	if(BufIdxIn > 0)
		CUtility::SafeWrite(hFiltInFile,szLineBufIn,BufIdxIn);
	close(hFiltInFile);
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Writing %d filtered in, retained, reads or regions to '%s' completed",NumFilteredIn,pszFiltInFile);
	}
if(hFiltOutFile != -1)
	{
	if(BufIdxOut > 0)
		CUtility::SafeWrite(hFiltOutFile,szLineBufOut,BufIdxOut);
	close(hFiltOutFile);
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Writing %d filtered out reads or regions to file '%s' completed",NumFilteredOut,pszFiltOutFile);
	}
return(Rslt);
}

int
Process(etPMode PMode,					// processing mode
		char Strand,					// retain regions on this strand
		etFRMode FMode,					// output format mode
		char *pszInFile,				// file containing alignments
		char *pszFiltInFile,			// where to write filtered in, retained, aligned reads or regions
		char *pszFiltOutFile,			// where to write filtered out aligned reads or regions
		int RegLen,						// 5'US and 3'DS regulatory region length
		int FeaturesIncl,				// retain reads which are overlapping these features
		int NumBedFiles,				// number of input BED filtering files containing features
		char *pszInBedFileSpecs[])		// input (wildcards allowed) bed feature filtering files
{
int Rslt;
int Idx;
char *pszBedFile;
Init();

if((Rslt=InitialMemAlloc())!=eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate required memory");
	return(Rslt);
	}

Rslt = LoadReads(Strand,pszInFile);
if(Rslt >= 0)
	{
	CSimpleGlob glob(SG_GLOB_FULLSORT);
	for(Idx = 0; Idx < NumBedFiles; Idx++)
		{
		glob.Init();
		if(glob.Add(pszInBedFileSpecs[Idx]) < SG_SUCCESS)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to glob '%s",pszInBedFileSpecs[Idx]);
			Reset();
			return(eBSFerrOpnFile);	// treat as though unable to open file
			}

		if(glob.FileCount() <= 0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to locate any source raw reads file matching '%s",pszInBedFileSpecs[Idx]);
			continue;
			}

		Rslt = eBSFSuccess;
		for (int FileID = 0; Rslt >= eBSFSuccess &&  FileID < glob.FileCount(); ++FileID)
			{
			pszBedFile = glob.File(FileID);
			Rslt = BEDFilter(Strand,			// filter for this strand of interest	
						RegLen,					// 5'US and 3'DS regulatory region length
						FeaturesIncl,			// which features are to be included
						pszBedFile);			// UCSC BED input file containing features of interest
			if(Rslt == eBSFSuccess)
				{
				gDiagnostics.DiagOut(eDLInfo,gszProcName,"All reads/regions have been filtered out");
				continue;
				}
			}	
		}
	}
if(Rslt >= eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"After filtering there are %d reads/regions remaining",Rslt);
	Rslt = WriteResults(pszFiltInFile,pszFiltOutFile);
	}
Reset();
return(Rslt);
}




// LocateRead
// Locate read with requested ReadID
tsAlignHit *
LocateRead(UINT32 ReadID)
{
int Rslt;
tsAlignHit *pProbe;
int Lo,Mid,Hi;	// search limits

if(m_ppAlignHitsIdx == NULL || m_AllocdAlignHitsIdx < m_NumAlignHits)
	return(NULL);

Lo = 0; 
Hi = m_NumAlignHits-1;
while(Hi >= Lo) {
	Mid = (Hi + Lo)/2;
	pProbe = m_ppAlignHitsIdx[Mid];
	Rslt = pProbe->ReadID - ReadID;
	if(Rslt > 0)	
		{
		Hi = Mid - 1;
		continue;
		}
	if(Rslt < 0)	
		{
		Lo = Mid + 1;
		continue;
		}
	return(pProbe);
	}
return(NULL);
}



// IterReads
// use to iterate over reads returning ptr to next read following the current read
// To start from first read then pass in NULL as pCurAlignHit
// Returns NULL if all read hits have been iterate
tsAlignHit *
IterReads(tsAlignHit *pCurAlignHit)
{
tsAlignHit *pNxtAlignHit = NULL;
if(pCurAlignHit == NULL)
	pNxtAlignHit = m_pAlignHits;
else
	if(pCurAlignHit->ReadID != m_FinalReadID)
		pNxtAlignHit = (tsAlignHit *)((UINT8 *)pCurAlignHit + sizeof(tsAlignHit));
return(pNxtAlignHit);
}

// IterSortedReads
// use to iterate over sorted reads returning ptr to next read following the current read
// To start from first read then pass in NULL as pCurAlignHit
// Returns NULL if all read hits have been iterated
tsAlignHit *
IterSortedReads(tsAlignHit *pCurAlignHit)
{
tsAlignHit *pNxtAlignHit = NULL;
if(pCurAlignHit == NULL)
	pNxtAlignHit = m_ppAlignHitsIdx[0];
else
	if(pCurAlignHit->AlignHitIdx < m_NumAlignHits)
		pNxtAlignHit = m_ppAlignHitsIdx[pCurAlignHit->AlignHitIdx];
return(pNxtAlignHit);
}

int
SortAlignHits(etReadsSortMode SortMode)
{
tsAlignHit *pAlignHit;
int Idx;

if(SortMode == m_CurReadsSortMode && m_ppAlignHitsIdx != NULL && m_AllocdAlignHitsIdx >= m_NumAlignHits)
	return(eBSFSuccess);

if(m_ppAlignHitsIdx == NULL || m_AllocdAlignHitsIdx < m_NumAlignHits)
	{
	if(m_ppAlignHitsIdx != NULL)
		{
		delete m_ppAlignHitsIdx;
		m_ppAlignHitsIdx = NULL;
		m_AllocdAlignHitsIdx = 0;
		}
	if((m_ppAlignHitsIdx = (tsAlignHit **) new tsAlignHit * [m_NumAlignHits])==NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"SortAlignHits: Memory reads index allocation for %d ptrs - %s",m_NumAlignHits,strerror(errno));
		return(eBSFerrMem);
		}
	m_AllocdAlignHitsIdx = m_NumAlignHits;
	}

pAlignHit = NULL;
tsAlignHit **pIdx = m_ppAlignHitsIdx;
while((pAlignHit = IterReads(pAlignHit))!=NULL)
	*pIdx++ = pAlignHit;

switch(SortMode) {
	case eRSMReadID: 
		qsort(m_ppAlignHitsIdx,m_NumAlignHits,sizeof(tsAlignHit *),SortHitReadIDs);
		break;
	case eRSMHitMatch:
		qsort(m_ppAlignHitsIdx,m_NumAlignHits,sizeof(tsAlignHit *),SortHitMatch);
		break;
	default:
		break;
	}

for(Idx = 0; Idx < (int)m_NumAlignHits; Idx++)
	m_ppAlignHitsIdx[Idx]->AlignHitIdx = Idx + 1;

m_CurReadsSortMode = SortMode;
return(eBSFSuccess);
}


// SortHitReadIDs
// Sort reads by ascending read identifiers
static int
SortHitReadIDs(const void *arg1, const void *arg2)
{
tsAlignHit *pEl1 = *(tsAlignHit **)arg1;
tsAlignHit *pEl2 = *(tsAlignHit **)arg2;

if(pEl1->ReadID < pEl2->ReadID )
		return(-1);
if(pEl1->ReadID > pEl2->ReadID )
	return(1);
return(0);
}

// SortHitmatch
// Sort by ascending chromid, loci, strand
static int
SortHitMatch(const void *arg1, const void *arg2)
{
tsAlignHit *pEl1 = *(tsAlignHit **)arg1;
tsAlignHit *pEl2 = *(tsAlignHit **)arg2;

if(pEl1->ChromID < pEl2->ChromID )
		return(-1);
if(pEl1->ChromID > pEl2->ChromID )
	return(1);
if(pEl1->Loci < pEl2->Loci )
		return(-1);
if(pEl1->Loci > pEl2->Loci )
	return(1);
if(pEl1->Strand < pEl2->Strand )
		return(-1);
if(pEl1->Strand > pEl2->Strand )
	return(1);

return(0);
}


