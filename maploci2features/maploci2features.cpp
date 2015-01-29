// maploci2features.cpp : Defines the entry point for the console application.
// Update: 4/08/2010 feature isoform processing added
// 2.0.2 bugfix for isoform determination (names could be reduced in length)
// 2.0.7 bugfix for BED files with large numbers ( > 60K chromosomes ) and missing BED title/descriptions
// 2.0.8 additional bugfix for BED files with large numbers ( > 60K chromosomes ) and missing BED title/descriptions
// 2.1.0 increased accepted size of the feature name when processing SAM files
// 2.1.1 increased diagnostics output
// 2.1.2 filtered out all unmapped reads
// 2.1.3 BED feature file must be specified (not optional)

#include "stdafx.h"
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "../libbiokanga/commhdrs.h"

const char *cpszProgVer = "2.1.4";		// increment with each release

const int cDfltJoinOverlap = 0;			// join cores which overlap and have starts which at most differ by cMaxJoinOverlap
const int cMaxJoinOverlap  = 100;		// max allowed join overlap

const int cMaxLengthRange = 1000000;	// maximal element length

// processing modes
typedef enum TAG_ePMode {
	ePMdefault,					// default is to associate features along full length of cores
	ePMstarts,					// associate features with starts only
	ePMdyad,					// associate features with assumed dyad positioned 73nt downstream of 5' core start
	ePMplaceholder				// used to set the enumeration range
	} etPMode;

// strand processing modes
typedef enum TAG_eStrandProc {
		eStrandDflt,			// default is to ignore the strand
		eStrandSense,			// process for sense
		eStrandAnti,			// process for antisense
		eStrandPlaceholder
} etStrandProc;

// feature isoform processing
typedef enum TAG_eISOFProc {
	eISOFPRPKM,					// default is to report the feature isoform with maximal RPKM
	eISOFReads,					// report the feature isoform with maximal total reads
	eISOFPall,				    // report all feature isoforms
	eISOFPlaceholder
} etISOFProc;

int Process(etPMode PMode,				// processing mode
			bool bDedupe,				// true if input elements are to be deduped
			bool bFeatinsts,			// true if input elements are to be associated to individual features, false if to all features at that locus
			bool bOneCntRead,			// true if one count per read rule to be applied (functional regions are prioritised with CDS as the highest) 
			etISOFProc IsoformRprt,		// feature isoform reporting mode
			etStrandProc StrandProc,	// how to process read + element strand
			int Ftype,					// input element file format: 0 - auto, 1 - CSV, 2 - BED, 3 - SAM
			teCSVFormat CSVFormat,		// expected input CSV loci file format
			char *pszInLociFile,		// input CSV loci file
			char *pszInBEDFile,			// input BED file
			char *pszRsltsFile,			// optional output loci mapping file
			char *pszFeatRsltsFile,		// optional output feature mapping results file
			char *pszSummRsltsFile,		// optional output chrom summary results file
			int RegRegionLen,			// regulatory region length
			int MinLength,				// minimum element length
			int MaxLength,				// maximum element length
			int JoinOverlap);			// deduping join overlap



char *CSVFormat2Text(teCSVFormat Format);


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

etPMode PMode;				// processing mode
int iCSVFormat;
int FType;					// expected input element file type - auto, CSV, BED or SAM
int iRegLen;				// up/down stream regulatory region length
int iMinLength;
int iMaxLength;
int iJoinOverlap;
bool bDedupe;				// if true then dedupe input elements
bool bFeatinsts;			// if true then associate to individual features
bool bOneCntRead;			// true if one count per read rule to be applied (functional regions are prioritised with CDS as the highest) 
etISOFProc IsoformRprt;			// feature isoform report processing
etStrandProc StrandProc;	// how to process read + element strand

char szInLociFile[_MAX_PATH];	// input element loci from this file
char szInBEDFile[_MAX_PATH];	// input bed file containing gene features
char szRsltsFile[_MAX_PATH];	// output loci + features to this file
char szFeatRsltsFile[_MAX_PATH]; // optional features output file
char szSummRsltsFile[_MAX_PATH]; // optional summary output file

// command line args
struct arg_lit  *help    = arg_lit0("hH","help",                "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"diagnostics log file");

struct arg_int *pmode = arg_int0("m","mode","<int>",		    "processing mode: 0 - full length attributes, 1 - 5' attributes associations, 2 - 5' dyad attribute associations (default = 0)");

struct arg_int *ftype = arg_int0("t","filetype","<int>",		"input element file format: 0 - auto, 1 - CSV, 2 - BED, 3 - SAM (default = 0)");
struct arg_lit  *dedupe = arg_lit0("d","dedupe",                "dedupe the input elements (default is to process all input elements)");
struct arg_int  *isoformrprt = arg_int0("x","isoform","<int>",  "feature isoform processing: 0 - report highest RPKM, 1 - report highest total reads, 2 - report all)");
struct arg_int  *strandproc = arg_int0("s","strandproc","<int>","strand processing: 0 - independent, 1 - sense, 2 - antisense (default is independent)");
struct arg_int  *CSVFormat = arg_int0("c","informat","<int>",	"if input CSV file processing 0:Loci 1:locsfx probe 2:locsfx target");
struct arg_file *InLociFile = arg_file1("i","inloci","<file>",	"input element (CSV, BED, SAM) file");
struct arg_file *InBEDFile = arg_file1("I","inbed","<file>",	"input gene or feature biobed BED file");
struct arg_file *RsltsFile = arg_file0("o","output","<file>",	 "optional mapped element output file");
struct arg_file *featrsltsfile = arg_file0("O","feats","<file>", "optional mapped features output file");
struct arg_file *summrsltsfile = arg_file0("C","chromsumm","<file>","optional chrom mapping summary output file");
struct arg_int  *RegLen = arg_int0("r","updnstream","<int>",	"length of 5'up or 3'down  stream regulatory region length (default = 2000) 0..1000000");
struct arg_int  *MinLength = arg_int0("l","minlength","<int>",   "minimum element length (default 0)");
struct arg_int  *MaxLength = arg_int0("L","maxlength","<int>",   "maximum element length (default 1000000)");
struct arg_int  *JoinOverlap = arg_int0("j","joinoverlap","<int>","joins cores which only differ by this many bases in start loci (default 0)");
struct arg_lit  *featinsts    = arg_lit0("z","featinsts",        "associate to individual features");
struct arg_lit  *onecntread   = arg_lit0("y","onecntread",       "prioritise functional regions so that there is only one count per read");

struct arg_end *end = arg_end(20);

void *argtable[] = {help,version,FileLogLevel,LogFile,
					pmode,dedupe,ftype,featinsts,isoformrprt,onecntread,strandproc,CSVFormat,InLociFile,InBEDFile,RsltsFile,featrsltsfile,summrsltsfile,RegLen,
					MinLength,MaxLength,JoinOverlap,
					end};

char **pAllArgs;
int argerrors;
argerrors = CUtility::arg_parsefromfile(argc,(char **)argv,&pAllArgs);
if(argerrors >= 0)
	argerrors = arg_parse(argerrors,pAllArgs,argtable);

    /* special case: '--help' takes precedence over error reporting */
if (help->count > 0)
        {
		printf("\n%s map loci to features, Version %s\nOptions ---\n", gszProcName,cpszProgVer);
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

	PMode = (etPMode)(pmode->count ? pmode->ival[0] : ePMdefault);
	if(PMode < ePMdefault || PMode >= ePMplaceholder)
		{
		printf("\nError: Processing mode '-m%d' specified outside of range %d..%d",PMode,0,(int)ePMplaceholder-1);
		exit(1);
		}
	bDedupe = dedupe->count ? true : false;
	bFeatinsts = featinsts->count ? true : false;
	bOneCntRead = onecntread->count ? true : false;

	FType = ftype->count ? ftype->ival[0] : 0;
	if(FType < 0 || FType >= 4)
		{
		printf("\nError: Expected input element file format '-t%d' specified outside of range %d..%d",FType,0,3);
		exit(1);
		}

	IsoformRprt = (etISOFProc)(isoformrprt->count ? isoformrprt->ival[0] : eISOFPRPKM);
	if(IsoformRprt < eISOFPRPKM || IsoformRprt >= eISOFPlaceholder)
		{
		printf("\nError: Isoform report mode '-x%d' must be in range %d..%d",IsoformRprt,eISOFPRPKM,eISOFPlaceholder-1);
		exit(1);
		}

	StrandProc = (etStrandProc)(strandproc->count ? strandproc->ival[0] : eStrandDflt);
	if(StrandProc < eStrandDflt || StrandProc >= eStrandPlaceholder)
		{
		printf("\nError: Strand processing mode '-s%d' must be in range %d..%d",StrandProc,eStrandDflt,eStrandAnti);
		exit(1);
		}
	
	if(FType <= 1)
		{
		iCSVFormat = CSVFormat->count ? CSVFormat->ival[0] : 0;
		if(iCSVFormat < eCSVFdefault || iCSVFormat > eCSVFtarget)
			{
			printf("\nError: expected input CSV format specified '-c%d' must be in range 0..2",iCSVFormat);
			exit(1);
			}
		}
	else
		iCSVFormat = 0;

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

	iJoinOverlap = JoinOverlap->count ? JoinOverlap->ival[0] : cDfltJoinOverlap;
	if(iJoinOverlap < 0 || iJoinOverlap > cMaxJoinOverlap)
		{
		printf("Error: join overlap length '-j%d' is not in range 0..%d",iJoinOverlap,cMaxJoinOverlap);
		exit(1);
		}

	if(InBEDFile->count)
		{
		strncpy(szInBEDFile,InBEDFile->filename[0],_MAX_PATH);
		szInBEDFile[_MAX_PATH-1] = '\0';
		iRegLen = RegLen->count ? RegLen->ival[0] : cDfltRegLen;
		if(iRegLen < cMinRegLen)
			{
			printf("\nRegulatory region length '-L%d' less than minimum %d, assuming you meant to use '-L%d'",iRegLen,cMinRegLen,cMinRegLen);
			iRegLen = cMinRegLen;
			}
		else
			{
			if(iRegLen > cMaxRegLen)
				{
				printf("\nRegulatory region length '-L%d' more than maximum %d, assuming you meant to use '-L%d'",iRegLen,cMaxRegLen,cMaxRegLen);
				iRegLen = cMaxRegLen;
				}
			}
		}
	else
		{
		printf("Error: No BED gene/feature file specified");
		exit(1);
		}

	strncpy(szInLociFile,InLociFile->filename[0],_MAX_PATH);
	szInLociFile[_MAX_PATH-1] = '\0';

	if(RsltsFile->count)
		{
		strncpy(szRsltsFile,RsltsFile->filename[0],_MAX_PATH);
		szRsltsFile[_MAX_PATH-1] = '\0';
		}
	else
		szRsltsFile[0] = '\0';

	if(featrsltsfile->count)
		{
		strncpy(szFeatRsltsFile,featrsltsfile->filename[0],_MAX_PATH);
		szFeatRsltsFile[_MAX_PATH-1] = '\0';
		}
	else
		szFeatRsltsFile[0] = '\0';

	if(szFeatRsltsFile[0] == '\0' && szRsltsFile[0] == '\0')
		{
		printf("Error: at least one of '-o<rsltsfile>' or '-O<featrsltsfile>' must be specified");
		exit(1);
		}

	if(summrsltsfile->count)
		{
		strncpy(szSummRsltsFile,summrsltsfile->filename[0],_MAX_PATH);
		szSummRsltsFile[_MAX_PATH-1] = '\0';
		}
	else
		szSummRsltsFile[0] = '\0';

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
			pszDescr = "associate features along full length of cores";
			break;
		case ePMstarts:
			pszDescr = "associate features with 5' core starts only";
			break;
		case ePMdyad:
			pszDescr = "associate features with assumed nucleosome dyad positioned 73nt downstream of 5' core start";
			break;

		}
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Processing mode is : '%s'",pszDescr);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Dedupe input elements : '%s'",bDedupe ? "yes" : "no");

	switch(IsoformRprt) {
		case eISOFPRPKM:
			pszDescr = "report the feature isoform with maximal RPKM";
			break;
		case eISOFReads:
			pszDescr = "report the feature isoform with maximal total reads";
			break;
		case eISOFPall:
			pszDescr = "report all feature isoforms";
			break;
		}
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Feature isoform report mode is : '%s'",pszDescr);

	switch(StrandProc) {
		case eStrandDflt:
			pszDescr = "ignor strand";
			break;
		case eStrandSense:
			pszDescr = "process for sense";
			break;
		case eStrandAnti:			
			pszDescr = "process for antisense";
			break;
		}

	
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Feature association is : '%s'",bFeatinsts ? "individual" : "all");
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"One count per read : '%s'",bOneCntRead ? "yes" : "no");
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Strand processing : '%s'",pszDescr);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Input element loci file: '%s'",szInLociFile);
	switch(FType) {
		case 0:
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Auto-classify input element file as either CSV, BED or SAM");
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"If CSV then as: %s",CSVFormat2Text((teCSVFormat)iCSVFormat));
			break;

		case 1:
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Expecting input element file to be CSV format as: %s",CSVFormat2Text((teCSVFormat)iCSVFormat));
			break;

		case 2:
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Expecting input element file to be BED format");
			break;

		case 3:
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Expecting input element file to be SAM format");
			break;
		}
	
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Element mapping results to file: '%s'",szRsltsFile[0]=='\0'?"None Specified":szRsltsFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Feature mapping results to file: '%s'",szFeatRsltsFile[0]=='\0'?"None Specified":szFeatRsltsFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Chrom summary results to file: '%s'",szSummRsltsFile[0]=='\0'?"None Specified":szSummRsltsFile);


	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Input biobed BED gene feature file: '%s'",szInBEDFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Feature mappings with regulatory region length: %d",iRegLen);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Minimum element length: %d",iMinLength);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Maximum element length: %d",iMaxLength);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Join overlap: %d",iJoinOverlap);

#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	gStopWatch.Start();
	Rslt = Process(PMode,bDedupe,bFeatinsts,bOneCntRead,IsoformRprt,StrandProc,FType,(teCSVFormat)iCSVFormat,szInLociFile,szInBEDFile,szRsltsFile,szFeatRsltsFile,szSummRsltsFile,iRegLen,iMinLength,iMaxLength,iJoinOverlap);

	gStopWatch.Stop();
	Rslt = Rslt >=0 ? 0 : 1;
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Exit code: %d Total processing time: %s",Rslt,gStopWatch.Read());
	exit(Rslt);
	}
else
	{
	printf("\n%s map loci to features, Version %s\n",gszProcName,cpszProgVer);
	arg_print_errors(stdout,end,gszProcName);
	arg_print_syntax(stdout,argtable,"\nUse '-h' to view option and parameter usage\n");
	exit(1);
	}
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
	default:
		break;
	}
return((char *)"Unsupported");
}


typedef struct TAG_sFeatCntDist {
	int FeatID;							// identifies this feature
	int GeneLen;						// gene length
	int TranscribedLen;					// transcribed length
	char szName[80];					// gene name
	int RegionCnts[8];					// region counts
	double RPKM;						// RPKM for this feature
	double RelAbundance;				// relative abundance (sum of all the reciprocals for each read's RelScale)
	double LenNormRelAbundance;			// RelAbundance normalised to transcript length
	double TransRelAbundance;			// RelAbundance as a proportion of the sum of all LenNormRelAbundance
	int UniqueReadLociHits;				// number of unique loci in this transcript to which at least one read mapped
	double UniqueHitsRelAbundance;		// RelAbundance normalised to number of UniqueReadLociHits in this transcript
	double TransUniqueHitsRelAbundance; // RelAbundance as a proportion of the sum of all UniqueHitsRelAbundance
	int NumExonReads;					// number of reads in UTR's plus CDS
	bool bIsIsoform;					// set true if feature is an isoform, assumed if name at least 9 long and suffixed by ".[0-99]"
	bool bMaxRPKM;						// set true if feature is an isoform with the maximal RPKM
	bool bMaxExonReads;				    // set true if feature is an isoform with the maximal reads in UTR's+CDS
	
} tsFeatCntDist;


int	m_FeatureDist[cMaxNumChroms][9];	

CHyperEls *m_pHypers;
tsFeatCntDist *m_pFeatCntDists;		// to hold feature count distributions
CBEDfile *m_pBiobed;

etPMode m_PMode;					// processing mode
etStrandProc m_StrandProc;			// strand processing
etISOFProc m_IsoformRprt;    		// feature isoform reporting mode
int m_RegRegionLen;					// regulatory region length
UINT32 m_NumEls;					// total number of elements to be processed

int m_NumSplitEls;					// number of elements identified as being split elements
bool m_bFeatinsts;					// true if counts are to be accumulated on a feature unique basis instead of any features covering the genomic locus
bool m_bOneCntRead;					// true if only single counts to be accumulated per read
int m_hRsltFile;

int m_hFeatRsltFile;

void
Reset(void)
{
if(m_hRsltFile != -1)
	{
	close(m_hRsltFile);
	m_hRsltFile = -1;
	}
if(m_hFeatRsltFile != -1)
	{
	close(m_hFeatRsltFile);
	m_hFeatRsltFile = -1;
	}
if(m_pBiobed != NULL)
	{
	delete m_pBiobed;
	m_pBiobed = NULL;
	}
if(m_pHypers != NULL)
	{
	delete m_pHypers;
	m_pHypers = NULL;
	}
if(m_pFeatCntDists != NULL)
	{
	delete m_pFeatCntDists;
	m_pFeatCntDists = NULL;
	}
m_PMode = ePMdefault; 
m_StrandProc = eStrandDflt;
m_IsoformRprt = eISOFPRPKM;
m_RegRegionLen = cDfltRegLen;
m_NumEls = 0;
m_NumSplitEls = 0;
m_bFeatinsts = false;
m_bOneCntRead = false;
}

void
Init(void)
{
m_pBiobed = NULL;
m_pHypers = NULL;
m_pFeatCntDists = NULL;
m_hRsltFile = -1;
m_hFeatRsltFile = -1;
Reset();
}

// MapLoci2Features
// Maps element loci to features and writes mapping for each element loci to pszRsltsFile
// Assumes that the bed file containing features (m_pBiobed) has been opened and that all element loci have been parsed into m_pHypers
int
MapLoci2Features(char *pszRsltsFile)
{
int Rslt;
char szLineBuff[0x03fff];
int BuffIdx;
int SrcID;
char Strand;
char *pszChrom;
int ChromID;
char *pszElType;
char *pszRefSpecies;
char *pszRelSpecies;
int PrevChromID;
int PrevElTypeID;
int PrevRefSpeciesID;
int PrevRelSpeciesID;
int StartLoci;
int EndLoci;
int Len;
int Features;
int AccumFeatures;
int NumFeatsOverlap;
int FeatMsk;
int FeatIdx;
UINT32 ElID;
int FeatID;
int NxtFeatID;
int NxtFeatStart;
int NxtFeatEnd;
int PrvFeatID;
int PrvFeatStart;
int PrvFeatEnd;
int RelScale;


int PrevStartChromID = -1;
int	PrevStartLoci = -1;
char PrevStrand = '*';
bool bStartUniqLoci = true;

int MaxChromID;
int CoreStartLoci;
int CoreEndLoci;
int TotNumFeatures;
tsFeatCntDist *pCurFeatCntDist;		// to hold currently being processed feature count distribution
tsHyperElement *pEl;

if(m_hRsltFile != -1)				// ensure closed
	{
	close(m_hRsltFile);
	m_hRsltFile = -1;
	}
if(pszRsltsFile != NULL && pszRsltsFile[0] != '\0')
	{
#ifdef _WIN32
	if((m_hRsltFile = open(pszRsltsFile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
	if((m_hRsltFile = open(pszRsltsFile,O_RDWR | O_CREAT | O_TRUNC,S_IREAD | S_IWRITE))==-1)
#endif
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create %s - %s",pszRsltsFile,strerror(errno));
		Reset();
		return(eBSFerrCreateFile);
		}
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Output file created/truncated: '%s'",pszRsltsFile);
	}

pszChrom = NULL;
PrevChromID = -1;
pszElType = NULL;
PrevElTypeID = -1;
pszRefSpecies = NULL;
PrevRefSpeciesID = -1;
pszRelSpecies = NULL;
PrevRelSpeciesID = -1;
BuffIdx = 0;
memset(m_FeatureDist,0,sizeof(m_FeatureDist));
MaxChromID = 0;
AccumFeatures = 0;

	// determine how many features
TotNumFeatures = m_pBiobed->GetNumFeatures();
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Feature file contains %1.9d features",TotNumFeatures);
if(m_pFeatCntDists == NULL)
	{
	if((m_pFeatCntDists = new tsFeatCntDist [TotNumFeatures])==NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate for %d tsFeatCntDist instances",TotNumFeatures);
		Reset();
		return(eBSFerrMem);
		}
	memset(m_pFeatCntDists,0,sizeof(tsFeatCntDist) * TotNumFeatures);
	pCurFeatCntDist = m_pFeatCntDists;
	for(FeatID=1;FeatID<=TotNumFeatures;FeatID++,pCurFeatCntDist++)
		{
		pCurFeatCntDist->FeatID = FeatID;
		pCurFeatCntDist->TranscribedLen = m_pBiobed->GetTranscribedLen(FeatID);
		pCurFeatCntDist->GeneLen = m_pBiobed->GetFeatLen(FeatID);
		m_pBiobed->GetFeature(FeatID,pCurFeatCntDist->szName);
		if(pCurFeatCntDist->szName[0] == '\0' || pCurFeatCntDist->GeneLen < 1)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Inconsistency in feature: %d",FeatID);
			Reset();
			return(eBSFerrInternal);
			}


		}
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Associating (from %d) element %9.9d",m_NumEls,1);
// Note: alignments will have been sorted ascending by chrom, start loci
for(ElID = 1; ElID <= m_NumEls; ElID++)
	{
	AccumFeatures = 0;
	if(!(ElID % 100000))
		printf("\b\b\b\b\b\b\b\b\b%9.9d",ElID);
	pEl = m_pHypers->GetElement(ElID);
	if(pEl == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to get details for element: %d",ElID);
		Reset();
		return(eBSFerrInternal);
		}

	// if only associating the actual start loci, and if a split element then need to process first if element on '+' and last if element on '-' strand
	if(m_PMode == ePMstarts && pEl->SplitElement == 1)
		{
		if(pEl->PlusStrand == 1 && pEl->SplitFirst != 1)
			continue;
		if(pEl->PlusStrand == 0 && pEl->SplitLast != 1)
			continue;
		}

	if(pszChrom == NULL || PrevChromID != pEl->ChromID)
		{
		pszChrom = m_pHypers->GetChrom(pEl->ChromID);
		if(pszChrom == NULL || pszChrom[0] == '\0')
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to get chrom text for element: %d ChromID: %d",ElID,pEl->ChromID);
			Reset();
			return(eBSFerrInternal);
			}

		// some old datasets may be referencing ChrM as mitochondria, or ChrC as chloroplast
		// so need to check for these
		if(!stricmp(pszChrom,"chloroplast"))
			pszChrom = (char *)"ChrC";
		else
			if(!stricmp(pszChrom,"mitochondria"))
				pszChrom = (char *)"ChrM";
		PrevChromID = pEl->ChromID;
		PrevStartChromID = -1;
		PrevStartLoci = -1;
		PrevStrand = '*';
		bStartUniqLoci = true;
		}
	
	if(pszElType == NULL || PrevElTypeID != pEl->ElTypeID)
		{
		pszElType = m_pHypers->GetType(pEl->ElTypeID);
		PrevElTypeID = pEl->ElTypeID;
		}

	if(pszRefSpecies == NULL || PrevRefSpeciesID != pEl->RefSpeciesID)
		{
		pszRefSpecies = m_pHypers->GetRefSpecies(pEl->RefSpeciesID);
		PrevRefSpeciesID = pEl->RefSpeciesID;
		}

	if(pszRelSpecies == NULL || PrevRelSpeciesID != pEl->RelSpeciesID)
		{
		pszRelSpecies = m_pHypers->GetRelSpecies(pEl->RelSpeciesID);
		PrevRelSpeciesID = pEl->RelSpeciesID;
		}

	
	StartLoci = pEl->StartLoci;
	EndLoci = pEl->StartLoci + pEl->Len - 1;

	Len = pEl->Len;
	Strand = pEl->PlusStrand ? '+' : '-';

	if(PrevStartChromID != pEl->ChromID || StartLoci != PrevStartLoci || Strand != PrevStrand) // not processed this loci previously?
		{
		PrevStartChromID = pEl->ChromID;
		PrevStartLoci = StartLoci;
		PrevStrand = Strand;
		bStartUniqLoci = true;
		}
	else
		bStartUniqLoci = false;

	SrcID = pEl->SrcID;
	Features = pEl->Features;
	RelScale = pEl->RelScale;
	switch(m_StrandProc) {
		case eStrandDflt:
			break;
		case eStrandSense:
			m_pBiobed->SetStrand(Strand);
			break;
		case eStrandAnti:
			m_pBiobed->SetStrand(Strand == '+' ? '-' : '+');
			break;
		}
	FeatID = 0;

	if(m_pBiobed != NULL)
		{
		Rslt= ChromID = m_pBiobed->LocateChromIDbyName(pszChrom);
		if(Rslt == eBSFerrChrom)
			{
			// some old datasets may be referencing ChrM as mitochondria, or ChrC as chloroplast
			// so need to check for these
			if(!stricmp(pszChrom,"ChrM"))
				Rslt= ChromID = m_pBiobed->LocateChromIDbyName((char *)"mitochondria");
			else
				if(!stricmp(pszChrom,"ChrC"))
					Rslt= ChromID = m_pBiobed->LocateChromIDbyName((char *)"chloroplast");
			}

		if(Rslt == eBSFerrChrom)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to locate chromosome %s in BED file",pszChrom);
			Reset();
			return(eBSFerrInternal);
			}

		if(Rslt > 0)
			{
			if(MaxChromID < ChromID)
				MaxChromID = ChromID;

			CoreStartLoci = StartLoci;
			CoreEndLoci = EndLoci;

			switch(m_PMode) {
				case ePMdefault:
					break;

				case ePMstarts:
					if(Strand == '+')
						CoreEndLoci = CoreStartLoci;
					else
						CoreStartLoci = CoreEndLoci;
					break;

				case ePMdyad:
					if(Strand == '+')
						{
						CoreStartLoci += 73;
						CoreEndLoci = CoreStartLoci;
						}
					else
						{
						CoreEndLoci -= 73;
						if(CoreEndLoci < 0)
							CoreEndLoci = 0;
						CoreStartLoci = CoreEndLoci;
						}
					break;
				}

			// see if overlapping any features
			AccumFeatures = 0;
			NumFeatsOverlap = 0;
			do
				{
				FeatID=m_pBiobed->LocateFeatureIDinRangeOnChrom(ChromID,	// feature is on which chromsome
										 CoreStartLoci,						// feature must end on or after Start
										 CoreEndLoci,						// and start on or before End 
										 NumFeatsOverlap+1);				// Ith instance to return (1..n)

				if(FeatID > 0 && m_pFeatCntDists != NULL)
					{
					if(m_bFeatinsts)
						{
						Features = m_pBiobed->GetFeatureOverlaps(cRegionFeatBits,FeatID,CoreStartLoci,CoreEndLoci,m_RegRegionLen);
						Features |= m_pBiobed->GetFeatureBitsSpliceOverlaps(FeatID,CoreStartLoci,CoreEndLoci,cMinSpliceOverlap);
						}
					else
						{
						Features = m_pBiobed->GetFeatureBits(ChromID,			// feature is on which chromsome
									 CoreStartLoci,							// feature must end on or after Start
									 CoreEndLoci,							// and start on or before End
									cRegionFeatBits,
									m_RegRegionLen);

						Features |= m_pBiobed->GetSpliceSiteBits(ChromID,CoreStartLoci,CoreEndLoci,cMinSpliceOverlap);
						}

					if(m_bOneCntRead)
						{
						for(FeatMsk = 0x01, FeatIdx = 0; FeatIdx < 8; FeatIdx++, FeatMsk <<= 1)
							if(Features & FeatMsk)
								{
								Features = FeatMsk;
								break;
								}
						}
					AccumFeatures |= Features;
					pCurFeatCntDist = &m_pFeatCntDists[FeatID-1];
					for(FeatMsk = 0x01, FeatIdx = 0; FeatIdx < 8; FeatIdx++, FeatMsk <<= 1)
						if(Features & FeatMsk)
							pCurFeatCntDist->RegionCnts[FeatIdx] += 1;

					// only accumulate relative abundance for reads in exons
					if(Features & (cFeatBitCDS | cFeatBit5UTR | cFeatBit3UTR))
						{
						pCurFeatCntDist->RelAbundance += 999.0 / (double)max(1,RelScale);
						if(bStartUniqLoci)
							pCurFeatCntDist->UniqueReadLociHits += 1;
						}
					}
				if(FeatID > 0)
					NumFeatsOverlap += 1;
				}
			while(FeatID > 0);
		 
			if(!NumFeatsOverlap) // if not overlapping or not contained in any feature then locate nearest feature up/dnstream
				{
				// find feature starting after core end loci
				NxtFeatID = m_pBiobed->LocateFeatureAfter(ChromID,	// feature is on this chromosome
							 CoreEndLoci);					         // feature starts on or immediately after this offset
				if(NxtFeatID > 0)
					m_pBiobed->GetFeature(NxtFeatID,		// feature instance identifier
							 NULL,							// where to return feature name
							 NULL,							// where to return chromosome name
							 &NxtFeatStart,					// where to return feature start on chromosome (0..n) 
							 &NxtFeatEnd);					// where to return feature end on chromosome
							
				// find feature ending before or at core start loci
				PrvFeatID = m_pBiobed->LocateFeatureBefore(ChromID,	// feature is on this chromosome
							 CoreStartLoci);			// feature ends on or immediately before this offset
				if(PrvFeatID > 0)
					m_pBiobed->GetFeature(PrvFeatID,		// feature instance identifier
							 NULL,	// where to return feature name
							 NULL,	// where to return chromosome name
							 &PrvFeatStart,		// where to return feature start on chromosome (0..n) 
							 &PrvFeatEnd);		// where to return feature end on chromosome
							
				if(NxtFeatID < 1)
					FeatID = PrvFeatID;
				else
					{
					if(PrvFeatID < 1)
						FeatID = NxtFeatID;
					else
						{
						if((NxtFeatStart - CoreEndLoci) < (CoreStartLoci - PrvFeatEnd))
							FeatID = NxtFeatID;
						else
							FeatID = PrvFeatID;
						}
					}


				if(m_bFeatinsts)
					AccumFeatures = m_pBiobed->GetFeatureOverlaps(cRegionFeatBits,FeatID,CoreStartLoci,CoreEndLoci,m_RegRegionLen);
				else
					AccumFeatures = m_pBiobed->GetFeatureBits(ChromID,			// feature is on which chromsome
									 CoreStartLoci,							// feature must end on or after Start
									 CoreEndLoci,							// and start on or before End
									cRegionFeatBits,
									m_RegRegionLen);
	
				if(m_bOneCntRead)
					{
					for(FeatMsk = 0x01, FeatIdx = 0; FeatIdx < 8; FeatIdx++, FeatMsk <<= 1)
						if(AccumFeatures & FeatMsk)
							{
							AccumFeatures = FeatMsk;
							break;
							}
					}
				if(AccumFeatures && FeatID > 0 && m_pFeatCntDists != NULL)
					{
					pCurFeatCntDist = &m_pFeatCntDists[FeatID-1];

					for(FeatMsk = 0x01, FeatIdx = 0; FeatIdx < 8; FeatIdx++, FeatMsk <<= 1)
						if(AccumFeatures & FeatMsk)
							pCurFeatCntDist->RegionCnts[FeatIdx] += 1;
					}
				}
			}
		}

	if(m_hRsltFile != -1)
		BuffIdx += sprintf(&szLineBuff[BuffIdx],"%d,\"%s\",\"%s\",\"%s\",%d,%d,%d,\"%c\",%d\n",
					SrcID,pszElType,pszRefSpecies,pszChrom,StartLoci,EndLoci,Len,Strand,AccumFeatures,RelScale);
	
	if(m_bOneCntRead)
		{
		for(FeatMsk = 0x01, FeatIdx = 0; FeatIdx < 8; FeatIdx++, FeatMsk <<= 1)
			if(AccumFeatures & FeatMsk)
				{
				AccumFeatures = FeatMsk;
				break;
				}
			}

	if(!AccumFeatures)
		m_FeatureDist[ChromID-1][8] += 1;
	else
		for(FeatMsk = 0x01, FeatIdx = 0; FeatIdx < 8; FeatIdx++, FeatMsk <<= 1)
			if(AccumFeatures & FeatMsk)
				m_FeatureDist[ChromID-1][FeatIdx] += 1;

	if(m_hRsltFile != -1 && ((BuffIdx + 1000) > sizeof(szLineBuff)))
		{
		CUtility::SafeWrite(m_hRsltFile,szLineBuff,BuffIdx);
		BuffIdx = 0;
		}
	}
printf("\b\b\b\b\b\b\b\b\b%9.9d",ElID-1);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"All elements (%d) now associated",ElID-1);
if(m_hRsltFile != -1)
	{
	if(BuffIdx > 0)
		CUtility::SafeWrite(m_hRsltFile,szLineBuff,BuffIdx);
	close(m_hRsltFile);
	m_hRsltFile = -1;
	}
return(eBSFSuccess);
}

// CompareFeatName
// Used to sort feature names ascending
static int 
CompareFeatName( const void *arg1, const void *arg2)
{
tsFeatCntDist *pEl1 = (tsFeatCntDist *)arg1;
tsFeatCntDist *pEl2 = (tsFeatCntDist *)arg2;
return(stricmp(pEl1->szName,pEl2->szName));
}


//
// IsSameFeature
// Feature names must be at least cMinNameRootLen chars long
// To be an isoform they must have suffixes of ".[0-99]"
// feature names first have any suffix trimmed off and then the remaining root names are compared for equality
// 
bool
IsSameFeature(char *pszFeatA,char *pszFeatB)
{
char ChrA;
char ChrB;
char *pSfxA;
char *pSfxB;
int NameLen;
bool bIsIsoform;
if(pszFeatA == NULL || *pszFeatA == '\0')
	return(false);

if(pszFeatB == NULL  || *pszFeatB == '\0')
	return(false);

NameLen = (int)strlen(pszFeatA);
if(NameLen < cMinIsonameLen)
	return(false);
pSfxA = &pszFeatA[NameLen-1];
if(*pSfxA >= '0' && *pSfxA <= '9')
	{
	pSfxA -= 1;
	if(*pSfxA >= '0' && *pSfxA <= '9')
		pSfxA -= 1;
	if(*pSfxA == '.')
		{
		*pSfxA = '\0';
		ChrA = '.';
		}
	else
		{
		ChrA = '\0';
		pSfxA = &pszFeatA[NameLen];
		}
	}
else
	{
	pSfxA = &pszFeatA[NameLen];
	ChrA = '\0';
	}

NameLen = (int)strlen(pszFeatB);
if(NameLen < cMinIsonameLen)
	return(false);
pSfxB = &pszFeatB[NameLen-1];
if(*pSfxB >= '0' && *pSfxB <= '9')
	{
	pSfxB -= 1;
	if(*pSfxB >= '0' && *pSfxB <= '9')
		pSfxB -= 1;
	if(*pSfxB == '.')
		{
		*pSfxB = '\0';
		ChrB = '.';
		}
	else
		{
		ChrB = '\0';
		pSfxB = &pszFeatB[NameLen];
		}
	}
else
	{
	ChrB = '\0';
	pSfxB = &pszFeatB[NameLen];
	}

bIsIsoform = stricmp(pszFeatA,pszFeatB) == 0 ? true : false;
*pSfxA = ChrA;
*pSfxB = ChrB;
return(bIsIsoform);
}


// TrimNameIso
// Inplace remove any name isoform suffix of the form '.[0-99]'
bool			// true if suffix was trimmed
TrimNameIso(char *pszName)
{
char *pSfx;
int NameLen;
if(pszName == NULL || pszName[0] == '\0')
	return(false);
NameLen = (int)strlen(pszName);
if(NameLen < cMinIsonameLen)
	return(false);
pSfx = &pszName[NameLen-1];
if(*pSfx >= '0' && *pSfx <= '9')
	{
	pSfx -= 1;
	if(*pSfx >= '0' && *pSfx <= '9')
		pSfx -= 1;
	if(*pSfx == '.')
		{
		*pSfx = '\0';
		return(true);
		}
	}
return(false);
}

// IsIsoform
// assumes that if the feature name is at least cMinNameRootLen long and suffixed by '.[0-99]' then thats an isoform
bool IsIsoform(char *pszName)
{
int NameLen;
char *pSfx;
NameLen = (int)strlen(pszName);
if(NameLen < cMinIsonameLen)
	return(false);
pSfx = &pszName[NameLen-1];
if(*pSfx >= '0' && *pSfx <= '9')
	{
	pSfx -= 1;
	if(*pSfx >= '0' && *pSfx <= '9')
		pSfx -= 1;
	if(*pSfx == '.')
		return(true);
	}
return(false);

}

int
MapFeatures2Loci(char *pszFeatRsltsFile)
{
char szLineBuff[0x03fff];
int BuffIdx;
int TotNumFeatures;
int FeatID;
int FeatIdx;
double RPKM;
double LenRelAbundance;
double SumTransLenRelAbundance;
double UniqueHitsRelAbundance;
double SumTransUniqueHitsRelAbundance;

tsFeatCntDist *pCurFeatCntDist;		// to hold currently being processed feature count distribution
tsFeatCntDist *pMaxRPKM;		// best RPKM isoform instance
tsFeatCntDist *pMaxExonReads;		// best reads isoform instance
bool bSameMaxRPKMFeature;		// true if processing same feature, different isoforms, for maximal RPKMs
bool bSameMaxExonReadsFeature;  // true if processing same feature,  different isoforms, for maximal exonic reads

if(m_hFeatRsltFile != -1)		// shouldn't be open but let's make sure...
	{
	close(m_hFeatRsltFile);
	m_hFeatRsltFile = -1;
	}
if(m_pBiobed != NULL && pszFeatRsltsFile != NULL && pszFeatRsltsFile[0] != '\0')
	{
	// determine how many features
	TotNumFeatures = m_pBiobed->GetNumFeatures();
	if(m_pFeatCntDists == NULL)
		{
		if((m_pFeatCntDists = new tsFeatCntDist [TotNumFeatures])==NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate for %d tsFeatCntDist instances",TotNumFeatures);
			Reset();
			return(eBSFerrMem);
			}
		memset(m_pFeatCntDists,0,sizeof(tsFeatCntDist) * TotNumFeatures);
		pCurFeatCntDist = m_pFeatCntDists;
		for(FeatID=1;FeatID<=TotNumFeatures;FeatID++,pCurFeatCntDist++)
			{
			pCurFeatCntDist->FeatID = FeatID;
			pCurFeatCntDist->TranscribedLen = m_pBiobed->GetTranscribedLen(FeatID);
			pCurFeatCntDist->GeneLen = m_pBiobed->GetFeatLen(FeatID);
			m_pBiobed->GetFeature(FeatID,pCurFeatCntDist->szName);
			}
		}

#ifdef _WIN32
	if((m_hFeatRsltFile = open(pszFeatRsltsFile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
	if((m_hFeatRsltFile = open(pszFeatRsltsFile, O_RDWR | O_CREAT | O_TRUNC,S_IREAD | S_IWRITE))==-1)
#endif
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create %s - %s",pszFeatRsltsFile,strerror(errno));
		Reset();
		return(eBSFerrCreateFile);
		}
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Features output file created/truncated: '%s'",pszFeatRsltsFile);
	}
else
	TotNumFeatures = 0;



if(m_hFeatRsltFile != -1 && m_pFeatCntDists != NULL)
	{
	// sort by feature name 
	qsort(m_pFeatCntDists,TotNumFeatures,sizeof(tsFeatCntDist),CompareFeatName);
	
	// generate RPKMs and ExonReads over all features and flag if feature is an isoform
	SumTransLenRelAbundance = 0.0;
	SumTransUniqueHitsRelAbundance = 0.0;
	pCurFeatCntDist = m_pFeatCntDists;
	pMaxRPKM = NULL;
	pMaxExonReads = NULL;
	for(FeatID=1;FeatID<=TotNumFeatures;FeatID++,pCurFeatCntDist++)
		{
		pCurFeatCntDist = &m_pFeatCntDists[FeatID-1];
		pCurFeatCntDist->NumExonReads = pCurFeatCntDist->RegionCnts[0] + pCurFeatCntDist->RegionCnts[1] + pCurFeatCntDist->RegionCnts[2];
		if(pCurFeatCntDist->NumExonReads)	// if at least 1 read mapping to feature...
			{
			RPKM = ((double)pCurFeatCntDist->NumExonReads * 1000.0f);
			RPKM /= (double)pCurFeatCntDist->TranscribedLen;
			RPKM *= 1000000.0f/(double)m_NumEls;
			LenRelAbundance = pCurFeatCntDist->RelAbundance / (double)pCurFeatCntDist->TranscribedLen;
			UniqueHitsRelAbundance = pCurFeatCntDist->RelAbundance / (double)pCurFeatCntDist->UniqueReadLociHits;
			}
		else
			{
			RPKM = 0.0f;
			LenRelAbundance = 0;
			UniqueHitsRelAbundance = 0;
			}
		pCurFeatCntDist->RPKM = RPKM;
		pCurFeatCntDist->LenNormRelAbundance = LenRelAbundance;
		pCurFeatCntDist->UniqueHitsRelAbundance = UniqueHitsRelAbundance;
		SumTransLenRelAbundance += LenRelAbundance;
		SumTransUniqueHitsRelAbundance += UniqueHitsRelAbundance;
		pCurFeatCntDist->bIsIsoform = IsIsoform(pCurFeatCntDist->szName);
		if(m_IsoformRprt == eISOFPall)		// if reporting all isoforms then that's easy...
			{
			pCurFeatCntDist->bMaxRPKM = true;
			pCurFeatCntDist->bMaxExonReads= true;
			continue;
			}	

		
		bSameMaxRPKMFeature = pMaxRPKM == NULL ? false : IsSameFeature(pCurFeatCntDist->szName,pMaxRPKM->szName);
		bSameMaxExonReadsFeature = pMaxExonReads == NULL ? false : IsSameFeature(pCurFeatCntDist->szName,pMaxExonReads->szName);

		// will report on 'maximal' isoforms
		// if first iteration (pMaxRPKM will be NULL) or not an isoform of current maximal RPKM ...
		if(pMaxRPKM==NULL || !bSameMaxRPKMFeature)
			{
			pCurFeatCntDist->bMaxRPKM = true;
			pMaxRPKM = pCurFeatCntDist;
			}

		// if first iteration (pMaxExonReads will be NULL) or not an isoform of current maximal MaxExonReads ...
		if(pMaxExonReads == NULL || !bSameMaxExonReadsFeature)
			{
			pCurFeatCntDist->bMaxExonReads= true;
			pMaxExonReads = pCurFeatCntDist;
			}

		// same feature - but different isoforms
		if(bSameMaxRPKMFeature && pCurFeatCntDist->RPKM > pMaxRPKM->RPKM)
			{
			pMaxRPKM->bMaxRPKM = false;
			pCurFeatCntDist->bMaxRPKM = true;
			pMaxRPKM = pCurFeatCntDist;
			}

		if(bSameMaxExonReadsFeature && pCurFeatCntDist->NumExonReads > pMaxExonReads->NumExonReads)
			{
			pMaxExonReads->bMaxExonReads = false;
			pCurFeatCntDist->bMaxExonReads = true;
			pMaxExonReads = pCurFeatCntDist;
			}
		}

	// now total is known, can determine abundance proportions for each transcript 
	if(SumTransLenRelAbundance > 0.0)
		{
		for(FeatID=1;FeatID<=TotNumFeatures;FeatID++,pCurFeatCntDist++)
			{
			pCurFeatCntDist = &m_pFeatCntDists[FeatID-1];
			pCurFeatCntDist->TransRelAbundance = pCurFeatCntDist->LenNormRelAbundance / SumTransLenRelAbundance;
			pCurFeatCntDist->TransUniqueHitsRelAbundance = pCurFeatCntDist->UniqueHitsRelAbundance / SumTransUniqueHitsRelAbundance;
			}
		}

	BuffIdx = sprintf(szLineBuff,"\"FeatID\",\"Feature\",\"GeneLen\",\"TransLen\",\"CDS\",\"5'UTR\",\"3'UTR\",\"Introns\",\"5'upstream\",\"3'downstream\",\"intron3'/5'exon\",\"exon3'/5'intron\",\"RPKM\",\"ExonReads\",\"UniqueLociHits\",\"RelAbundance\",\"UniqueRelAbundance\",\"TransUniqueRelAbundance\",\"LenRelAbundance\",\"TransLenRelAbundance\"");
	pCurFeatCntDist = m_pFeatCntDists;

	for(FeatID=1;FeatID<=TotNumFeatures;FeatID++,pCurFeatCntDist++)
		{
		switch(m_IsoformRprt) {
			case eISOFPRPKM:
				if(!pCurFeatCntDist->bMaxRPKM)
					continue;
				break;
			case eISOFReads:
				if(!pCurFeatCntDist->bMaxExonReads)
					continue;
				break;
			default:
				break;
			}
	  	BuffIdx += sprintf(&szLineBuff[BuffIdx],"\n%d,\"%s\",%d,%d",
				FeatID,pCurFeatCntDist->szName,pCurFeatCntDist->GeneLen,pCurFeatCntDist->TranscribedLen);
		pCurFeatCntDist = &m_pFeatCntDists[FeatID-1];
		for(FeatIdx = 0; FeatIdx < 8; FeatIdx++)
			BuffIdx += sprintf(&szLineBuff[BuffIdx],",%d",pCurFeatCntDist->RegionCnts[FeatIdx]);
		BuffIdx += sprintf(&szLineBuff[BuffIdx],",%1.6f,%d,%d,%1.9f,%1.9f,%1.9f,%1.9f,%1.9f",
						pCurFeatCntDist->RPKM,pCurFeatCntDist->NumExonReads,
						pCurFeatCntDist->UniqueReadLociHits,
						pCurFeatCntDist->RelAbundance,
						pCurFeatCntDist->UniqueHitsRelAbundance, 
						pCurFeatCntDist->TransUniqueHitsRelAbundance,
						pCurFeatCntDist->LenNormRelAbundance,
						pCurFeatCntDist->TransRelAbundance);

		// something really strange on Ubuntu - sometimes a number of features are replicated at the end of the generated file
		// almost as though the final write was occuring 2x but I can't determine why so am playing safe now and manually writing out the
		// buffer on the final iteration and also doing a fsync immediately before closing the handle
		if(FeatID == TotNumFeatures || ((BuffIdx + 1000) > sizeof(szLineBuff)))
			{
			if(!CUtility::SafeWrite(m_hFeatRsltFile,szLineBuff,BuffIdx))
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to write %s",pszFeatRsltsFile);
				Reset();
				return(eBSFerrCreateFile);
				}
			BuffIdx = 0;
			}
		}
	if(BuffIdx)
		CUtility::SafeWrite(m_hFeatRsltFile,szLineBuff,BuffIdx);
#ifdef _WIN32
	_commit(m_hFeatRsltFile);
#else
	fsync(m_hFeatRsltFile);
#endif
	BuffIdx = 0;
	close(m_hFeatRsltFile);
	m_hFeatRsltFile = -1;
	}
if(m_pFeatCntDists != NULL)
	{
	delete m_pFeatCntDists;
	m_pFeatCntDists = NULL;
	}
return(eBSFSuccess);
}


int Process(etPMode PMode,				// processing mode
			bool bDedupe,				// true if input elements are to be deduped
			bool bFeatinsts,			// true if input elements are to be associated to individual features, false if to all features at that locus
			bool bOneCntRead,			// true if one count per read rule to be applied (functional regions are prioritised with CDS as the highest) 
			etISOFProc IsoformRprt,		// feature isoform reporting mode
			etStrandProc StrandProc,	// how to process read + element strand
			int FType,					// input element file format: 0 - auto, 1 - CSV, 2 - BED, 3 - SAM
			teCSVFormat CSVFormat,		// if CSV input then expected file format
			char *pszInLociFile,		// input CSV, BED or SAM loci file
			char *pszInBEDFile,			// input BED file
			char *pszRsltsFile,			// output loci to feature mapping file
			char *pszFeatRsltsFile,		// optional feature mapping results file
			char *pszSummRsltsFile,		// optional output chrom summary results file
			int RegRegionLen,			// regulatory region length
			int MinLength,				// minimum element length
			int MaxLength,				// maximum element length
			int JoinOverlap)			// deduping join overlap
{
int Rslt;
char *pszChrom;
int ChromID;
char szChrom[cMaxDatasetSpeciesChrom];
int TmpNumEls;

Init();
m_PMode = PMode;
m_StrandProc = StrandProc;
m_IsoformRprt = IsoformRprt;
m_RegRegionLen = RegRegionLen;
m_bFeatinsts = bFeatinsts;
m_bOneCntRead = bOneCntRead;

m_pHypers = new CHyperEls;

etClassifyFileType FileType;

if(pszInBEDFile != NULL && pszInBEDFile[0] != '\0')
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Loading features file '%s'",pszInBEDFile);
	if((m_pBiobed = new CBEDfile) == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CBEDfile gene/exon file '%s'",pszInBEDFile);
		Reset();
		return(eBSFerrObj);
		}

	if((Rslt=m_pBiobed->Open(pszInBEDFile,eBTGeneExons))!=eBSFSuccess)
		{
		while(m_pBiobed->NumErrMsgs())
			gDiagnostics.DiagOut(eDLFatal,gszProcName,m_pBiobed->GetErrMsg());
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open features file '%s'",pszInBEDFile);
		Reset();
		return(eBSFerrObj);
		}
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Opened features file '%s'",pszInBEDFile);
	}
else
	{
	m_pBiobed = NULL;
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"No gene/exon BED file specified");
	Reset();
	return(eBSFerrObj);
	}

if(FType == 0)
	FileType = CUtility::ClassifyFileType(pszInLociFile);
else
	FileType = (etClassifyFileType)(FType - 1);

switch(FileType) {
	case eCFTopenerr:		// unable to open file for reading
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open file: '%s'",pszInLociFile);
		return(eBSFerrOpnFile);

	case eCFTlenerr:		// file length is insufficent to classify type
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Unable to classify file type (insuffient data points): '%s'",pszInLociFile);
		return(eBSFerrFileAccess);

	case eCFTunknown:		// unable to reliably classify
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Unable to reliably classify file type: '%s'",pszInLociFile);
		return(eBSFerrFileType);

	case eCFTCSV:			// file has been classified as being CSV
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Parsing CSV file (%d..%d): '%s'",MinLength,MaxLength,pszInLociFile);
		if((Rslt = m_pHypers->ParseCSVFileElements(pszInLociFile,MinLength,MaxLength,CSVFormat)) < 0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Parse errors in CSV file (%d..%d): '%s'",MinLength,MaxLength,pszInLociFile);
			Reset();
			return(Rslt);
			}
		break;

	case eCFTBED:			// file has been classified as being BED
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Parsing BED file (%d..%d): '%s'",MinLength,MaxLength,pszInLociFile);
		if((Rslt = m_pHypers->ParseBEDFileElements(pszInLociFile,MinLength,MaxLength)) < 0)
			{	
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Parse errors in BED file (%d..%d): '%s'",MinLength,MaxLength,pszInLociFile);
			Reset();
			return(Rslt);
			}
		break;

	case eCFTSAM:			// file has been classified as being SAM
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Parsing SAM file (%d..%d): '%s'",MinLength,MaxLength,pszInLociFile);
		if((Rslt = m_pHypers->ParseSAMFileElements(pszInLociFile,MinLength,MaxLength)) < 0)
			{	
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Parse errors in SAM file (%d..%d): '%s'",MinLength,MaxLength,pszInLociFile);
			Reset();
			return(Rslt);
			}
		break;
	}


m_NumEls = m_pHypers->NumEls();
if(m_NumEls == 0)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"No elements with length range %d..%d in file: '%s'",MinLength,MaxLength,pszInLociFile);
	Reset();
	return(Rslt);
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Loaded and parsed %d elements",m_NumEls);

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Now identifying split (microInDels or splice junction spanning?) elements...",m_NumEls);
m_NumSplitEls = m_pHypers->IdentifySplitElements();					// identify any split elements which may be present
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Identified %d split (microInDels or splice junction spanning?) elements",m_NumSplitEls);

if(JoinOverlap > 0 || bDedupe)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Now joining/deduping elements...",m_NumEls);
	TmpNumEls = m_NumEls;
	m_NumEls = m_pHypers->DedupeSort(JoinOverlap,bDedupe);				// dedupe, join and sort elements ascending
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"After join/dedupe, %d elements removed, %d elements remaining",
					 TmpNumEls - m_NumEls, m_NumEls);
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Mapping elements to features...");
if((Rslt=MapLoci2Features(pszRsltsFile)) < eBSFSuccess)
	{
	Reset();
	return(Rslt);
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Mapping elements to features completed");

if(pszFeatRsltsFile != NULL && pszFeatRsltsFile[0] != '\0')
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Mapping features to elements...");
	if((Rslt=MapFeatures2Loci(pszFeatRsltsFile)) < eBSFSuccess)
		{
		Reset();
		return(Rslt);
		}
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Mapping features to elements completed");
	}

if(pszSummRsltsFile != NULL && pszSummRsltsFile[0] != '\0')
	{
	char szOutBuff[8000];
	int BufIdx;
	int hOutFile;
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Generating summary chromosome results file...");

#ifdef _WIN32
	if((hOutFile = open(pszSummRsltsFile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
	if((hOutFile = open(pszSummRsltsFile, O_RDWR | O_CREAT | O_TRUNC, S_IREAD | S_IWRITE ))==-1)
#endif
		{
		gDiagnostics.DiagOutMsgOnly(eDLFatal,"Unable to create or truncate  output file %s error: %s",pszSummRsltsFile,strerror(errno));
		return(-1);
		}
	BufIdx = sprintf(szOutBuff,"Chrom, CDS, 5'UTR, 3'UTR, Introns, 5'upstream, 3'downstream, intron3'/5'exon, exon3'/5'intron, Intergenic\n");
	for(ChromID = 0; ChromID < m_pBiobed->GetNumChromosomes(); ChromID++)
		{
		m_pBiobed->GetChromosome(ChromID+1,szChrom);
		pszChrom = szChrom;
		if(!stricmp(pszChrom,"chloroplast"))
			pszChrom = (char *)"ChrC";
		else
			if(!stricmp(pszChrom,"mitochondria"))
				pszChrom = (char *)"ChrM";
		BufIdx += sprintf(&szOutBuff[BufIdx],"\"%s\",%d,%d,%d,%d,%d,%d,%d,%d,%d\n",
						 pszChrom,m_FeatureDist[ChromID][0],m_FeatureDist[ChromID][1],m_FeatureDist[ChromID][2],m_FeatureDist[ChromID][3],m_FeatureDist[ChromID][4],m_FeatureDist[ChromID][5],m_FeatureDist[ChromID][6],m_FeatureDist[ChromID][7],m_FeatureDist[ChromID][8]);
		if((BufIdx + 1000) > sizeof(szOutBuff))
			{
			CUtility::SafeWrite(hOutFile,szOutBuff,BufIdx);
			BufIdx = 0;
			}
		}
	if(BufIdx)
		CUtility::SafeWrite(hOutFile,szOutBuff,BufIdx);
	close(hOutFile);
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Generating summary chromosome results file completed");
	}

Reset();
return(Rslt < 0 ? m_NumEls : Rslt);
}


