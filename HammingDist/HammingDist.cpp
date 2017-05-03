// HammingDist.cpp : Defines the entry point for the console application.
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

const char *cpszProgVer = "0.1.1";		// increment with each release
const int cMaxInFiles = 50;				// allow for at most 50 input files
const int cMaxInFileSpecs = 10;			// allow for at most 10 input file specs which may be wildcarded

// input loci file processing modes
typedef enum TAG_ePMode {		
	ePMdefault,					// CSV hamming loci
	ePMplaceholder				// used to set the enumeration range
	} etPMode;


// strand processing modes
typedef enum TAG_eStrandProc {
		eStrandDflt,			// default is to ignore the strand
		eStrandWatson,			// process for Watson
		eStrandCrick,			// process for Crick
		eStrandPlaceholder
} etStrandProc;


int Process(etPMode PMode,				// processing mode
			char Strand,				// which element strand to filter (retain) on
			int RegRegionLen,			// regulatory region length
			int OfsLoci,				// offset start loci by this relative offset
			int NumInputFiles,			// number of input files
			char *pszInfileSpecs[],		// input element hamming loci from these files
			char *pszInBEDFile,			// input BED feature file
			char *pszRsltsFile); 		// output loci file

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
	return _T("kanga");
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
int Idx;

etPMode PMode;				// processing mode
int iRegLen;				// up/down stream regulatory region length
int RelOfs;					// relative start offset 
int Strand;					// which element strand to filter (retain) on

int NumInputFiles;			// number of input files
char *pszInfileSpecs[cRRMaxInFileSpecs];  // input element hamming loci from these files
char szInBEDFile[_MAX_PATH];	// input bed file containing gene features
char szRsltsFile[_MAX_PATH];	// output result distributions to this file

// command line args
struct arg_lit  *help    = arg_lit0("hH","help",                "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"diagnostics log file");

struct arg_int *pmode = arg_int0("m","mode","<int>",		    "input reads loci file: 0 - CSV (default = 0)");

struct arg_int  *strandproc = arg_int0("s","strandproc","<int>","strand processing: 0 - independent, 1 - Watson, 2 - Crick (default is independent)");

struct arg_file *InLociFiles = arg_filen("i","incsv","<file>",0,cRRMaxInFileSpecs,"input element CSV or BED files");
struct arg_file *InBEDFile = arg_file0("I","infeats","<file>",	"input gene or feature biobed BED file");
struct arg_file *RsltsFile = arg_file1("o","output","<file>",	 "length distributions output file");
struct arg_int  *RegLen = arg_int0("r","updnstream","<int>",	 "length of 5'up or 3'down  stream regulatory region length (default = 2000, range 0..1000000)");
struct arg_int  *relofs = arg_int0("R","relofs","<int>",	     "relative loci offset (default = 0, range +/-200)");

struct arg_end *end = arg_end(20);

void *argtable[] = {help,version,FileLogLevel,LogFile,
					pmode,strandproc,RegLen,relofs,InLociFiles,InBEDFile,RsltsFile,
					end};

char **pAllArgs;
int argerrors;
argerrors = CUtility::arg_parsefromfile(argc,(char **)argv,&pAllArgs);
if(argerrors >= 0)
	argerrors = arg_parse(argerrors,pAllArgs,argtable);

/* special case: '--help' takes precedence over error reporting */
if (help->count > 0)
        {
		printf("\n%s HammingDist - Hamming Distance Distributions , Version %s\nOptions ---\n", gszProcName,cpszProgVer);
        arg_print_syntax(stdout,argtable,"\n");
        arg_print_glossary(stdout,argtable,"  %-25s %s\n");
		printf("\nNote: Parameters can be entered into a parameter file, one parameter per line.");
		printf("\n      To invoke this parameter file then precede its name with '@'");
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
		printf("\nError: FileLogLevel '-f%d' specified but no logfile '-F<logfile>\n'",FileLogLevel->ival[0]);
		exit(1);
		}

	iScreenLogLevel = iFileLogLevel = FileLogLevel->count ? FileLogLevel->ival[0] : eDLInfo;
	if(iFileLogLevel < eDLNone || iFileLogLevel > eDLDebug)
		{
		printf("\nError: FileLogLevel '-l%d' specified outside of range %d..%d\n",iFileLogLevel,eDLNone,eDLDebug);
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

	// now that log parameters have been parsed then initialise diagnostics log system
	if(!gDiagnostics.Open(szLogFile,(etDiagLevel)iScreenLogLevel,(etDiagLevel)iFileLogLevel,true))
		{
		printf("\nError: Unable to start diagnostics subsystem\n");
		if(szLogFile[0] != '\0')
			printf(" Most likely cause is that logfile '%s' can't be opened/created\n",szLogFile);
		exit(1);
		}

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Version: %s",cpszProgVer);

	PMode = (etPMode)(pmode->count ? pmode->ival[0] : ePMdefault);
	if(PMode < ePMdefault || PMode >= ePMplaceholder)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Processing mode '-m%d' specified outside of range %d..%d",PMode,0,(int)ePMplaceholder-1);
		exit(1);
		}
	
	Strand = strandproc->count ? strandproc->ival[0] : eStrandDflt;
	if(Strand < eStrandDflt || Strand >= eStrandPlaceholder)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Strand processing mode '-s%d' must be in range %d..%d",Strand,eStrandDflt,eStrandCrick);
		exit(1);
		}
	switch(Strand) {
		case 1: Strand = (int)'+'; break;
		case 2: Strand = (int)'-'; break;
		case 0: Strand = (int)'*'; break;
		}	

	if(InBEDFile->count)
		{
		strncpy(szInBEDFile,InBEDFile->filename[0],_MAX_PATH);
		szInBEDFile[_MAX_PATH-1] = '\0';
		iRegLen = RegLen->count ? RegLen->ival[0] : cDfltRegLen;
		if(iRegLen < cMinRegLen || iRegLen > cMaxRegLen)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Regulatory region length '-r%d' must be in range %d..%d",iRegLen,cMinRegLen,cMinRegLen);
			exit(1);
			}

		RelOfs = relofs->count ? relofs->ival[0] : 0;
		if(RelOfs < -200 || RelOfs > 200)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Relative offset '-R%d' must be in range %d..%d",RelOfs,-200,200);
			exit(1);
			}
		}
	else
		{
		szInBEDFile[0] = '\0';
		iRegLen = 0;
		RelOfs = 0;
		}


	for(NumInputFiles=Idx=0;NumInputFiles < cMaxInFileSpecs && Idx < InLociFiles->count; Idx++)
		{
		pszInfileSpecs[Idx] = NULL;
		if(pszInfileSpecs[NumInputFiles] == NULL)
			pszInfileSpecs[NumInputFiles] = new char [_MAX_PATH];
		strncpy(pszInfileSpecs[NumInputFiles],InLociFiles->filename[Idx],_MAX_PATH);
		pszInfileSpecs[NumInputFiles][_MAX_PATH-1] = '\0';
		CUtility::TrimQuotedWhitespcExtd(pszInfileSpecs[NumInputFiles]);
		if(pszInfileSpecs[NumInputFiles][0] != '\0')
			NumInputFiles++;
		}

	if(!NumInputFiles)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: After removal of whitespace, no input file(s) specified with '-i<filespec>' option)");
		exit(1);
		}

	strncpy(szRsltsFile,RsltsFile->filename[0],_MAX_PATH);
	szRsltsFile[_MAX_PATH-1] = '\0';

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing parameters:");

	const char *pszDescr;
	switch(PMode) {
		case ePMdefault:
			pszDescr = "CSV loci";
			break;
		}
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"processing mode for input loci file is : '%s'",pszDescr);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Strand processing : '%c'",(char)Strand);

	for(Idx=0; Idx < NumInputFiles; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"input raw reads files (%d): '%s'",Idx+1,pszInfileSpecs[Idx]);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Length distributions to file: '%s'",szRsltsFile);

	if(szInBEDFile[0] != '\0')
		{
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Input biobed BED gene feature file: '%s'",szInBEDFile);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Feature mappings with regulatory region length: %d",iRegLen);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Feature mappings with relative loci offset: %d",RelOfs);
		}

#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	gStopWatch.Start();
	Rslt = Process(PMode,Strand,iRegLen,RelOfs,NumInputFiles,pszInfileSpecs,szInBEDFile,szRsltsFile);

	gStopWatch.Stop();
	Rslt = Rslt >=0 ? 0 : 1;
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Exit code: %d Total processing time: %s",Rslt,gStopWatch.Read());
	exit(Rslt);
	}
else
	{
	printf("\n%s HammingDist, Version %s\n",gszProcName,cpszProgVer);
	arg_print_errors(stdout,end,gszProcName);
	arg_print_syntax(stdout,argtable,"\nUse '-h' to view option and parameter usage\n");
	exit(1);
	}
return 0;
}


const int cNumRegions = 7;		// regions: cFeatBitUpstream, cFeatBit5UTR, cFeatBitCDS, cFeatBitIntrons, cFeatBit3UTR, cFeatBitDnstream and intergenic
const int cMaxHamming = 200;   // allow for hammings of up to this distance 

typedef struct TAG_sRegHammDist {
	int Region;
	char szRegion[50];
	int HighestHamm;			// highest hamming in this region
	UINT32 Cnts[cMaxHamming+1];
    } tsRegHammDist;

tsRegHammDist m_RegHammDists[cNumRegions];	// counts for each region

CBEDfile *m_pFeatFile;		// features file
CCSVFile *m_pHammings;		// hammings file
int m_hRsltsFile;;			// results file handle

bool m_bRegions;			// set true if regional distributions being processed
int m_MaxHamming;			// maximum Hamming at any loci
int m_TotNumProcessed;		// total number of loci over all files 

void
Init(void)
{
int Region;
int RegionMsk;
tsRegHammDist *pRegHammDist;

m_pFeatFile = NULL;
m_pHammings = NULL;
m_hRsltsFile = -1;

m_MaxHamming = -1;
m_TotNumProcessed = 0;
m_bRegions = false;
memset(m_RegHammDists,0,sizeof(tsRegHammDist));
pRegHammDist = m_RegHammDists; 
for(RegionMsk = Region = 1; Region <= cNumRegions; Region++, pRegHammDist++, RegionMsk <<= 1)
	{
	pRegHammDist->Region = Region;
	switch(RegionMsk) {
		case cFeatBitCDS:			// part of feature overlaps CDS
			strcpy(pRegHammDist->szRegion,"CDS"); 
			break;
		case cFeatBit5UTR:				// part of feature overlaps 5'UTR
			strcpy(pRegHammDist->szRegion,"UTR5"); 
			break;
		case cFeatBit3UTR:			// part of feature overlaps 3'UTR
			strcpy(pRegHammDist->szRegion,"UTR3"); 
			break;
		case cFeatBitIntrons:		// part of feature overlaps Intron
			strcpy(pRegHammDist->szRegion,"Intron"); 
			break;
		case cFeatBitUpstream:		// part of feature overlaps 5'upstream regulatory
			strcpy(pRegHammDist->szRegion,"UP5"); 
			break;
		case cFeatBitDnstream:		// part of feature overlaps 3'downstream regulatory
			strcpy(pRegHammDist->szRegion,"DN3"); 
			break;
		default:					// any other is intergenic
			strcpy(pRegHammDist->szRegion,"Intergenic"); 
			break;
		}
	}
}

void
Reset(bool bSync)
{
if(m_hRsltsFile != -1)
	{
	if(bSync)
#ifdef _WIN32
		_commit(m_hRsltsFile);
#else
		fsync(m_hRsltsFile);
#endif
	close(m_hRsltsFile);
    m_hRsltsFile = -1;
	}
if(m_pFeatFile != NULL)
	{
	delete m_pFeatFile;
	m_pFeatFile = NULL;
	}

if(m_pHammings != NULL)
	{
	delete m_pHammings;
	m_pHammings = NULL;
	}
}



int
LoadHammings(int OfsLoci,int RegRegionLen,char *pszInLociFile)
{
int Rslt;
int ChromID;
int NumProcessed;
int MaxHamming;
int NumFields;
char *pszChrom;
int StartLoci;
int Hamming;
int FeatureBits;
int Region;
char szPrevChrom[128];
tsRegHammDist *pRegHammDist;

if(m_pHammings != NULL)
	delete m_pHammings;

m_pHammings= new CCSVFile;
if(m_pHammings == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CCSVfile");
	Reset(false);
	return(eBSFerrObj);
	}

if((Rslt=m_pHammings->Open(pszInLociFile))!=eBSFSuccess)
	{
	while(m_pHammings->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,m_pHammings->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open file: %s",pszInLociFile);
	Reset(false);
	return(Rslt);
	}

// iterate over hamming loci
NumProcessed = 0;
MaxHamming = -1;
ChromID = 0;
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processed %d",m_TotNumProcessed);
while((Rslt=m_pHammings->NextLine()) > 0)	// onto next line containing fields
	{
	NumFields = m_pHammings->GetCurFields();
	if(NumFields < 3)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expected 3+ fields in '%s', GetCurFields() returned '%d'",pszInLociFile,NumFields);
		Rslt = eBSFerrFieldCnt;
		break;
		}

	if(!NumProcessed && m_pHammings->IsLikelyHeaderLine())
		continue;

	NumProcessed += 1;

	if((Rslt = m_pHammings->GetText(1,&pszChrom))!=eBSFSuccess)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to parse chrom name at element %d",NumProcessed);
		break;
		}
	if((Rslt = m_pHammings->GetInt(2,&StartLoci))!=eBSFSuccess)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to parse Hamming start loci element %d",NumProcessed);
		break;
		}

	if(m_bRegions)
		{
		if((Rslt = m_pHammings->GetInt(3,&Hamming))!=eBSFSuccess)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to parse chrom name at element %d",NumProcessed);
			break;
			}

		// now try and characterise into region
		if(ChromID == 0 || pszChrom[0] != szPrevChrom[0] || stricmp(pszChrom,szPrevChrom))
			{
			if((ChromID = m_pFeatFile->LocateChromIDbyName(pszChrom)) < 1)
				{
				Rslt = ChromID;
				break;
				}
			strcpy(szPrevChrom,pszChrom);
			}

		StartLoci += OfsLoci;
		if(StartLoci < 0)
			StartLoci = 0;
		FeatureBits = m_pFeatFile->GetFeatureBits(ChromID,StartLoci,StartLoci,cRegionFeatBits,RegRegionLen);

		// now allocate counts...
		for(Region = 0; Region < (cNumRegions-1); Region++)
			{
			if(FeatureBits & 0x01)
				break;
			FeatureBits >>= 1;
			}
		}
	else
		Region = 0;

	pRegHammDist = &m_RegHammDists[Region];
	if(pRegHammDist->HighestHamm < Hamming)
		pRegHammDist->HighestHamm = Hamming;
	if(MaxHamming < Hamming)
		MaxHamming = Hamming;
	pRegHammDist->Cnts[Hamming] += 1;

	if(!(NumProcessed % 10000000))
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processed %d",NumProcessed);
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Rslt: %d, Processed %d",Rslt,NumProcessed);
if(Rslt >= 0)
	{
	m_TotNumProcessed += NumProcessed;
	if(m_MaxHamming < MaxHamming)
		m_MaxHamming = MaxHamming;
	}
return(NumProcessed);
}


int Process(etPMode PMode,				// processing mode
			char Strand,				// which element strand to filter (retain) on
			int RegRegionLen,			// regulatory region length
			int OfsLoci,				// offset start loci by this relative offset
			int NumInputFiles,			// number of input files
			char *pszInfileSpecs[],		// input element hamming loci from these files
			char *pszInBEDFile,			// input BED feature file
			char *pszRsltsFile) 		// output loci file
{
int Rslt;
int Idx;

int NumInputFilesProcessed;
char *pszInfile;
UINT32 RegionTotCnts[cNumRegions];
double RegionCummDist[cNumRegions];

tsRegHammDist *pRegHammDist;
int Region;

int Hamming;

Init();

if(pszInBEDFile != NULL && pszInBEDFile[0] != '\0')
	{
	if((m_pFeatFile = new CBEDfile) == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CBEDfile object");
		Reset(false);
		return(eBSFerrObj);
		}

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Loading: %s",pszInBEDFile);
	if((Rslt = m_pFeatFile->Open(pszInBEDFile))!=eBSFSuccess)
		{
		while(m_pFeatFile->NumErrMsgs())
			gDiagnostics.DiagOut(eDLFatal,gszProcName,m_pFeatFile->GetErrMsg());
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open bioseq feature file '%s'",pszInBEDFile);
		Reset(false);
		return(Rslt);
		}
	m_bRegions = true;
	}
else
	{
	m_bRegions = false;
	strcpy(m_RegHammDists[0].szRegion,"All");
	}

#ifdef _WIN32
if((m_hRsltsFile = open(pszRsltsFile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
if((m_hRsltsFile = open(pszRsltsFile, O_RDWR | O_CREAT |O_TRUNC, S_IREAD | S_IWRITE))==-1)
#endif
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create %s - %s",pszRsltsFile,strerror(errno));
	Reset(false);
	return(eBSFerrCreateFile);
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Output file created/truncated: '%s'",pszRsltsFile);


m_MaxHamming = -1;
m_TotNumProcessed = 0;
NumInputFilesProcessed = 0;
pszInfile = NULL;
CSimpleGlob glob(SG_GLOB_FULLSORT);
for(Idx = 0; Idx < NumInputFiles; Idx++)
	{
	glob.Init();
	if(glob.Add(pszInfileSpecs[Idx]) < SG_SUCCESS)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to glob '%s",pszInfileSpecs[Idx]);
		Reset(false);
		return(eBSFerrOpnFile);	// treat as though unable to open file
		}

	if(glob.FileCount() <= 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to locate any input loci Hamming file matching '%s",pszInfileSpecs[Idx]);
		continue;
		}

	Rslt = eBSFSuccess;
	for (int FileID = 0; Rslt >= eBSFSuccess &&  FileID < glob.FileCount(); ++FileID)
		{
		pszInfile = glob.File(FileID);
		NumInputFilesProcessed += 1;
		if(NumInputFilesProcessed > cMaxInFiles)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Too many input files (max allowed is %d)",cRRMaxInFiles);
			Reset(false);
			return(eBSFerrNumSrcFiles);
			}
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing Hammings from loci file '%s'\n",pszInfile);
		Rslt = LoadHammings(OfsLoci,RegRegionLen,pszInfile);
		if(Rslt < eBSFSuccess)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: Load failed for input Hammings loci file '%s'\n",pszInfile);
			Reset(false);
			return((teBSFrsltCodes)Rslt);
			}
		}
	}


memset(RegionTotCnts,0,sizeof(RegionTotCnts));
memset(RegionCummDist,0,sizeof(RegionCummDist));
if(m_MaxHamming >= 0 && m_TotNumProcessed > 0)		// any to report?
	{
	int BuffIdx;
	int NumRegions;
	char szBuff[4096];
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Write distributions...");
	BuffIdx = 0;
	NumRegions = m_bRegions ? cNumRegions : 1;
	pRegHammDist = m_RegHammDists;
	for(Region = 0; Region < NumRegions; Region++, pRegHammDist++)
		BuffIdx += sprintf(&szBuff[BuffIdx],",\"%s\"",pRegHammDist->szRegion);
	pRegHammDist = m_RegHammDists;
	for(Region = 0; Region < NumRegions; Region++, pRegHammDist++)
		BuffIdx += sprintf(&szBuff[BuffIdx],",\"Proportion %s\"",pRegHammDist->szRegion);
	pRegHammDist = m_RegHammDists;
	for(Region = 0; Region < NumRegions; Region++, pRegHammDist++)
		BuffIdx += sprintf(&szBuff[BuffIdx],",\"Cumulative %s\"",pRegHammDist->szRegion);

	if(write(m_hRsltsFile,szBuff,BuffIdx)!=BuffIdx)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error write to file '%s",pszRsltsFile);
		Reset(false);
		return(eBSFerrFileAccess);
		}

	for(Hamming = 0; Hamming < m_MaxHamming; Hamming++)
		{
		pRegHammDist = m_RegHammDists;
		for(Region = 0; Region < NumRegions; Region++, pRegHammDist++)
			RegionTotCnts[Region] += pRegHammDist->Cnts[Hamming]; 
		}


	BuffIdx = 0;
	for(Hamming = 0; Hamming < m_MaxHamming; Hamming++)
		{
		pRegHammDist = m_RegHammDists;
		BuffIdx += sprintf(&szBuff[BuffIdx],"\n%d",Hamming);
		for(Region = 0; Region < NumRegions; Region++, pRegHammDist++)
			{
			BuffIdx += sprintf(&szBuff[BuffIdx],",%d",pRegHammDist->Cnts[Hamming]);
			if((BuffIdx + 200) > sizeof(szBuff))
				{
				if(write(m_hRsltsFile,szBuff,BuffIdx)!=BuffIdx)
					{
					gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error on write to file '%s",pszRsltsFile);
					Reset(false);
					return(eBSFerrFileAccess);
					}
				BuffIdx = 0;
				}
			}

		pRegHammDist = m_RegHammDists;
		for(Region = 0; Region < NumRegions; Region++, pRegHammDist++)
			{
			BuffIdx += sprintf(&szBuff[BuffIdx],",%f",RegionTotCnts[Region] > 0 ? pRegHammDist->Cnts[Hamming]/(double)RegionTotCnts[Region] : 0.0);
			if((BuffIdx + 200) > sizeof(szBuff))
				{
				if(write(m_hRsltsFile,szBuff,BuffIdx)!=BuffIdx)
					{
					gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error on write to file '%s",pszRsltsFile);
					Reset(false);
					return(eBSFerrFileAccess);
					}
				BuffIdx = 0;
				}
			}

		pRegHammDist = m_RegHammDists;
		for(Region = 0; Region < NumRegions; Region++, pRegHammDist++)
			{
			RegionCummDist[Region] += RegionTotCnts[Region] > 0 ? (double)pRegHammDist->Cnts[Hamming]/(double)RegionTotCnts[Region] : 0.0;
			BuffIdx += sprintf(&szBuff[BuffIdx],",%f",RegionCummDist[Region] );
			if((BuffIdx + 200) > sizeof(szBuff))
				{
				if(write(m_hRsltsFile,szBuff,BuffIdx)!=BuffIdx)
					{
					gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error on write to file '%s",pszRsltsFile);
					Reset(false);
					return(eBSFerrFileAccess);
					}
				BuffIdx = 0;
				}
			}
		}

	if(BuffIdx)
		{
		if(write(m_hRsltsFile,szBuff,BuffIdx)!=BuffIdx)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error on write to file '%s",pszRsltsFile);
			Reset(false);
			return(eBSFerrFileAccess);
			}
		BuffIdx = 0;
		}
	}
Reset(true);
return(eBSFSuccess);
}

