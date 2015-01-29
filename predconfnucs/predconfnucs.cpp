// prednucleosomes.cpp : Defines the entry point for the console application.
// predicts nucleosome postitions from conformational characteristics - minor groove + groove orientation (from accumulated twist)
// predict nucleosome dyad potential from -
// a) dyad (ChkGroove[6]) above background by default 1.5%
// b) average flanking ChkGroove[5 and 7] above background by default 1.25%
// c) average all other ChkGrooves above background by default 1.0%
// putative dyads are crudely scored by their percentages above background and if multiple dyads are separated by a single nucleotide then
// that dyad with a local maxima is the dyad reported 


#include "stdafx.h"
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#if _WIN32
#include "../libbiokanga/commhdrs.h"
#else
#include "../libbiokanga/commhdrs.h"
#endif

const unsigned int cProgVer = 120;		// increment with each release

const double cDfltDyadratio = 1.020f;	// default dyad grooves must be at least this ratio to background
const double cDfltDyad2ratio = 1.015f;	// default immediately flanking grooves must be at least this ratio to background
const double cDfltDyad3ratio = 1.010f;	// default remainder of flanking grooves must be at least this ration to background

const int cAllocOutBuff = 1024000;		// results output buffer allocation size

// processing modes
typedef enum TAG_ePMode {		
	ePMdefault,					// default processing mode
	ePMplaceholder				// used to set the enumeration range
	} etPMode;

// output format modes
typedef enum TAG_eFMode {		
	eFMbedGraphDyads,			// default is for ucsc bedGraph dyads
	eFMbedDyads,				// ucsc bed format dyads
	eFMcsvDyads,				// CSV format dyads
	eFMbedGraphNucs,			// ucsc bedGraph nucleosomes
	eFMbedNucs,					// ucsc bed nucleosomes
	eFMcsvNucs,					// CSV format nucleosomes
	eFMMcsvScores,				// CSV all scores along genome
	eFMplaceholder				// used to set the enumeration range
	} etFMode;

typedef struct TAG_sRegion {
	int RegionID;					// uniquely identifies this region
	int NumHits;					// number of times this region contained a putative nucleosome
	int ChromID;					// region is on this chromosome
	int StartLoci;					// loci at which region starts
	int EndLoci;					// loci at which region ends
} tsRegion;

int
Process(etPMode PMode,				// processing mode
		etFMode FMode,				// output format mode
		int MovAvgFilter,				// apply this moving average window size filter
		int BaselineFilter,				// baseline normalisation window size
		char *pszTrackName,			// UCSC track name
		double DyadratioThres,		// dyad grooves must be at least this ratio to background
		double Dyad2ratioThres,		// immediately flanking grooves must be at least this ratio to background
		double Dyad3ratioThres,		// remainder of flanking grooves must be at least this ration to background
		char *pszInGenomefile,		// bioseq genome file
		char *pszInConfFile,		// file containing conformation characteristics
		char *pszOutFile,			// where to write nucleosome predictions
		char *pszInclRegionFile,	//  only report predicted nucleosomes if intersecting with regions in this file
		int OfsLoci,				// offset region start loci by this many nt
		int DeltaLen,				// change region length by this many nt
		int TruncLength);			// truncate regions to be a maximum of this length


void Init(void);
void Reset(void);

int SortRegions( const void *arg1, const void *arg2);

CStopWatch gStopWatch;
CDiagnostics gDiagnostics;				// for writing diagnostics messages to log file
char gszProcName[_MAX_FNAME];			// process name

etPMode m_PMode;						// processing mode
etFMode m_FMode;						// output format mode

CCSVFile *m_pCSV;					// used if loading include/exclude regions from CSV file
CBioSeqFile *m_pBioSeqFile;			// holds instantiated genome assembly sequences, used if MNase scoring
CTwister *m_pTwister;				// contains hunter group dsDNA conformational characteristic values

int m_hOutFile;						// file handle for output results
char *m_pszOutBuff;					// used for buffering output results
int m_AllocdOutBuff;				// how many chars have been allocated to m_pszOutBuff
int m_UsedOutBuff;					// how many chars of results currently in m_pszOutBuff

etSeqBase *m_pChromSeq;				// holds current assembly chromosome sequence
int *m_pConfGroove;					// conformational values corresponding to sequence in m_pChromSeq for minor groove
int *m_pConfTwist;					// conformational values corresponding to sequence in m_pChromSeq for rotational twist
int *m_pScores;						// scores along m_pChromSeq


int m_AllocdChromSeq;				// allocated assembly chromosome sequence length
int m_ChromSeqLen;					// current assembly chromosome sequence length
char m_szCurChrom[cMaxDatasetSpeciesChrom+1];		// chromosome currently being processed
int m_NumPutDyads;					// total number of putative dyads predicted
int m_ChromPutDyads;				// number of putative dyads on current chromosome
char m_szTrackName[cMaxDatasetSpeciesChrom+1];	// user specified track name

tsRegion *m_pRegions;				// pts to include/exclude regions
int m_AllocdRegions;				// memory allocated for this many regions
int m_NumRegions;					// this many regions have been loaded
int m_NumHitRegions;				// number of times regions were hit at least once

double m_DyadratioThres;			// dyad (central) grooves must be at least this ratio to background
double m_Dyad2ratioThres;			// immediately flanking grooves must be at least this ratio to background
double m_Dyad3ratioThres;			// remainder of flanking grooves must be at least this ratio to background

double m_DyadRatio;					// current dyad (central) groove ratio to background
double m_Dyad2Ratio;				// current, immediately flanking, groove ratio to background
double m_Dyad3Ratio;				// current, remainder of flanking grooves, ratio to background

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
	return _T("predconfnucs");
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

etPMode PMode;				// processing mode
etFMode FMode;				// format output mode

int TruncLength;			// truncate regions to be no longer than this length
int OfsLoci;				// offset region loci by this many bases
int DeltaLen;				// change region lengths by this many bases

int MovAvgFilter;			// apply this moving average window size filter
int BaselineFilter;			// use this window size when normalising for baseline
char szTrackTitle[cMaxDatasetSpeciesChrom];			// title to identify predicted nucleosome track
double Dyadratio;
double Dyad2ratio;
double Dyad3ratio;
char szRsltsFile[_MAX_PATH];
char szInConfFile[_MAX_PATH];
char szInGenomeFile[_MAX_PATH];

char szInclRegionFile[_MAX_PATH]; //  only report predicted nucleosomes if intersecting with regions in this file


// command line args
struct arg_lit  *help    = arg_lit0("hH","help",                "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"diagnostics log file");

struct arg_int *pmode = arg_int0("m","mode","<int>",		    "processing mode: 0 - predict from minor groove (default = 0)");
struct arg_file *inGenomefile = arg_file1("i","in","<file>",	"input from this bioseq genome file");
struct arg_file *inConfFile = arg_file1("I","conf","<file>",	"input conformation characteristics from this file");

struct arg_file *inclRegionFile = arg_file0("r","inclregions","<file>",	"only report predicted nucleosomes if intersecting with regions in this file");
struct arg_int  *trunclength = arg_int0("T","trunclength","<int>","truncate regions to be no longer than this length (default 147)");
struct arg_int  *ofsloci   = arg_int0("u","offset","<int>",	    "offset region loci by this many bases, -2048 to +2048 (default 0)");
struct arg_int  *deltalen  = arg_int0("U","deltalen","<int>",	"delta region lengths by this many bases, -2048 to +2048 (default 0)");

struct arg_dbl *dyadratio=arg_dbl0("d", "dyadratio",	"<double>","dyad minor grooves must be at least this ratio to background (default 1.020");
struct arg_dbl *dyad2ratio=arg_dbl0("D", "dyad2ratio",	"<double>","immediately flanking minor grooves must be at least this ratio to background (default 1.015");
struct arg_dbl *dyad3ratio=arg_dbl0("e", "dyad3ratio",	"<double>","remainder flanking minor grooves must be at least this ratio to background (default 1.010");

struct arg_int *movavgfilter = arg_int0("a","avgwindow","<int>","apply lowpass filter as a moving average window of this size, 0 if none else 5..100 (default: 10)");
struct arg_int *baselinefilter = arg_int0("A","basewindow","<int>","baseline normalisation window size, 0 if none else 25..5000 (default: 250)");

struct arg_file *outfile = arg_file1("o","out","<file>",		"output nucleosome predictions to this file");
struct arg_int *format = arg_int0("M","format","<int>",		    "output format: 0 - UCSC bedGraph dyads, 1 - UCSC BED dyads, 2 - CSV dyads,  3 - UCSC bedGraph nucleosomes, 4 - UCSC BED nucleosomes, 5 - CSV nucleosomes, 6 - CSV scores (default: 0)");
struct arg_str  *title = arg_str0("t","title","<string>",       "track title");

struct arg_end *end = arg_end(20);

void *argtable[] = {help,version,FileLogLevel,LogFile,
					pmode,format,movavgfilter,baselinefilter,title,dyadratio,dyad2ratio,dyad3ratio,inGenomefile,inConfFile,outfile,
					inclRegionFile,trunclength,ofsloci,deltalen,
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

	PMode = (etPMode)(pmode->count ? pmode->ival[0] : 0);
	if(PMode < 0 || PMode >= ePMplaceholder)
		{
		printf("\nError: Processing mode '-m%d' specified outside of range %d..%d",PMode,0,(int)ePMplaceholder-1);
		exit(1);
		}

	FMode = (etFMode)(format->count ? format->ival[0] : eFMbedGraphDyads);
	if(FMode < eFMbedGraphDyads || FMode >= eFMplaceholder)
		{
		printf("\nError: Output format mode '-m%d' specified outside of range %d..%d",FMode,eFMbedGraphDyads,(int)eFMplaceholder-1);
		exit(1);
		}

		MovAvgFilter = movavgfilter->count ? movavgfilter->ival[0] : 10;
	if(MovAvgFilter != 0 && (MovAvgFilter < 5 || MovAvgFilter > 100))
		{
		printf("\nError: Moving average filter width '-a%d' specified outside of range 0 or 5..100",MovAvgFilter);
		exit(1);
		}


	BaselineFilter = baselinefilter->count ? baselinefilter->ival[0] : 250;
	if(BaselineFilter != 0 && (BaselineFilter < 25 || MovAvgFilter > 5000))
		{
		printf("\nError: Baseline normalisation window width '-A%d' specified outside of range 0 or 25..5000",BaselineFilter);
		exit(1);
		}

	Dyadratio = dyadratio->count ? dyadratio->dval[0] : cDfltDyadratio;
	Dyad2ratio = dyad2ratio->count ? dyad2ratio->dval[0] : cDfltDyad2ratio;
	Dyad3ratio = dyad3ratio->count ? dyad3ratio->dval[0] : cDfltDyad3ratio;
	if(Dyadratio < 1.0f || Dyad2ratio < 1.0f || Dyad3ratio < 1.0f)
		{
		printf("\nError: All dyad and flanking threshold ratios must be at least 1.0");
		exit(1);
		}


	if(!title->count)
		{
		printf("\nError: no track title has been specified with '-t<title>' option");
		exit(1);
		}
	if(strlen(title->sval[0]) >= cMaxDatasetSpeciesChrom)
		{
		printf("\nError: track title '%s' is too long, must be shorter than %d chars",title->sval[0],cMaxDatasetSpeciesChrom-1);
		exit(1);
		}
	strncpy(szTrackTitle,title->sval[0],sizeof(szTrackTitle));
	szTrackTitle[sizeof(szTrackTitle)-1] = '\0';
	CUtility::TrimQuotedWhitespcExtd(szTrackTitle);
	CUtility::ReduceWhitespace(szTrackTitle);
	if(szTrackTitle[0] == '\0')
		{
		printf("\nError: specified track title is empty or whitespace only");
		exit(1);
		}

	strcpy(szInGenomeFile,inGenomefile->filename[0]);
	strcpy(szInConfFile,inConfFile->filename[0]);
	strcpy(szRsltsFile,outfile->filename[0]);

	if(inclRegionFile->count)
		{
		strncpy(szInclRegionFile,inclRegionFile->filename[0],_MAX_PATH);
		szInclRegionFile[_MAX_PATH-1] = '\0';

		OfsLoci = ofsloci->count ? ofsloci->ival[0] : 0;
		if(abs(OfsLoci) > 1024)
			{
			printf("Error: loci offset '-u%d' must be in the range -2048 to +2048",OfsLoci);
			exit(1);
			}

		DeltaLen = deltalen->count ? deltalen->ival[0] : 0;
		if(abs(DeltaLen) > 1024)
			{
			printf("Error: delta length '-U%d' must be in the range -2048 to +2048",DeltaLen);
			exit(1);
			}

		TruncLength = trunclength->count ? trunclength->ival[0] : 147;
		if(TruncLength < 147)
			{
			printf("Error: regions truncation length '-T%d' must be at least 147",TruncLength);
			exit(1);
			}
		}
	else
		{
		szInclRegionFile[0] = '\0';
		OfsLoci = 0;
		DeltaLen = 0;
		TruncLength = 0;
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
	const char *pszDescr;
	switch(PMode) {
		case ePMdefault:				
			pszDescr = "Default - predict from minor groove";
			break;
		}

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"processing mode is : '%s'",pszDescr);

	switch(FMode) {
		case eFMbedGraphDyads:				
			pszDescr = "UCSC bedGraph dyads";
			break;
		case eFMbedDyads:					
			pszDescr = "UCSC BED dyads";
			break;
		case eFMcsvDyads:					
			pszDescr = "CSV dyads";
			break;
		case eFMbedGraphNucs:					
			pszDescr = "UCSC bedGraph nucleosomes";
			break;
		case eFMbedNucs:					
			pszDescr = "UCSC BED nucleosomes";
			break;
		case eFMcsvNucs:					
			pszDescr = "CSV nucleosomes";
			break;
		case eFMMcsvScores:					
			pszDescr = "CSV dyad scores along genome";
			break;

		}
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"output format is : '%s'",pszDescr);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"moving average filter width : %d",MovAvgFilter);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"baseline noramlisation width : %d",BaselineFilter);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"dyad name: '%s'",szTrackTitle);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"dyad minor groove threshold: %1.4f",Dyadratio);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"dyad immediately flanking minor groove threshold: %1.4f",Dyad2ratio);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"dyad remainder flanking minor groove threshold: %1.4f",Dyad3ratio);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"input genome file: '%s'",szInGenomeFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"input conformation file: '%s'",szInConfFile);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"output predicted nucleosomes file: '%s'",szRsltsFile);
	
	if(szInclRegionFile[0] != '0')
		{
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Offset region start loci by : %d",OfsLoci);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Delta region length by : %d",DeltaLen);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Truncate region length: %d",TruncLength);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"include regions file: '%s'",szInclRegionFile);
		}
#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	gStopWatch.Start();
	Rslt = Process(PMode,FMode,MovAvgFilter,BaselineFilter,szTrackTitle,Dyadratio,Dyad2ratio,Dyad3ratio,szInGenomeFile,szInConfFile,szRsltsFile,
									szInclRegionFile,OfsLoci,DeltaLen,TruncLength);
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
return 0;
}



void
Init(void)
{
m_pCSV = NULL;
m_pBioSeqFile = NULL;
m_pTwister = NULL;
m_pChromSeq = NULL;
m_pConfGroove = NULL;
m_pConfTwist = NULL;
m_pScores = NULL;
m_pszOutBuff = NULL;
m_pRegions = NULL;
m_hOutFile = -1;
Reset();
}

void
Reset(void)
{
if(m_hOutFile != -1)
	{
	close(m_hOutFile);
	m_hOutFile = -1;
	}
if(m_pCSV != NULL)
	{
	delete m_pCSV;
	m_pCSV = NULL;
	}
if(m_pszOutBuff != NULL)
	{
	delete m_pszOutBuff;
	m_pszOutBuff = NULL;
	}
if(m_pBioSeqFile != NULL)
	{
	delete m_pBioSeqFile;
	m_pBioSeqFile = NULL;
	}
if(m_pTwister != NULL)
	{
	delete m_pTwister;
	m_pTwister = NULL;
	}
if(m_pChromSeq != NULL)
	{
	delete m_pChromSeq;
	m_pChromSeq = NULL;
	}
if(m_pConfGroove != NULL)
	{
	delete m_pConfGroove;
	m_pConfGroove = NULL;
	}
if(m_pConfTwist != NULL)
	{
	delete m_pConfTwist;
	m_pConfTwist = NULL;
	}
if(m_pScores != NULL)
	{
	delete m_pScores;
	m_pScores = NULL;
	}
if(m_pRegions != NULL)
	{
	free(m_pRegions);				// NOTE - use free() as memory malloc'd/realloc'd
	m_pRegions = NULL;
	}
m_AllocdOutBuff = 0;
m_UsedOutBuff = 0;
m_AllocdRegions = 0;				// memory allocated for this many regions
m_NumRegions = 0;					// this many regions have been loaded
m_NumHitRegions = 0;				// no regions hit as yet
m_AllocdChromSeq = 0;				// allocated assembly chromosome sequence length
m_ChromSeqLen = 0;					// current assembly chromosome sequence length
m_szCurChrom[0] = '\0';				// no chromosome currently loaded
m_NumPutDyads = 0;					// no dyads yet predicted
m_ChromPutDyads = 0;				// none yet 
m_szTrackName[0] = '\0';			// no track name
}

const int cRegionAlloc = 10000000;	// allocate for inclusive regions in this many sized region blocks
int
AddRegion(int RegionID,char *pszChrom,int StartLoci,int EndLoci)
{
int ChromID;
tsRegion *pRegion;
if((ChromID = m_pBioSeqFile->LocateEntryIDbyName(pszChrom)) < 1)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to locate region chromosome '%s'",pszChrom);
	return(eBSFerrChrom);
	}

if(m_pRegions == NULL || m_NumRegions == m_AllocdRegions)
	{
	if(m_pRegions == NULL)
		{
		if((pRegion = (tsRegion *)malloc(sizeof(tsRegion) * cRegionAlloc))==NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to malloc memory (%d bytes requested) for region filtering",sizeof(tsRegion) * cRegionAlloc);
			Reset();
			return(eBSFerrMem);
			}
		m_AllocdRegions = 0;
		m_NumRegions = 0;
		}
	else
		{
		if((pRegion = (tsRegion *)realloc(m_pRegions,sizeof(tsRegion) * (m_AllocdRegions+cRegionAlloc)))==NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to realloc memory (%d bytes requested) for region filtering",sizeof(tsRegion) * (m_AllocdRegions+cRegionAlloc));
			Reset();
			return(eBSFerrMem);
			}
		}
	m_pRegions = pRegion;
	m_AllocdRegions += cRegionAlloc;
	}
pRegion = &m_pRegions[m_NumRegions++];
pRegion->RegionID = RegionID;
pRegion->ChromID = ChromID;
pRegion->StartLoci = StartLoci;
pRegion->EndLoci = EndLoci;
pRegion->NumHits = 0;
return(m_NumRegions);
}


int			// actual length, or if under MinLen then -1, or if start would be < 0 or end < start then returns -2
AdjLoci(bool bOnMinStrand,		// true if element is on '-' strand
		int *pStartLoci,		// starting loci on '+' strand
		int *pEndLoci,			// ending loci on '+' strand
		int OfsLoci,			// offset loci by
		int DeltaLen,			// adjust length by
		int MinLen,				// must be at least this length
		int TruncLen)			// truncate to this length
{
int Ofs;
int DLen;

if(*pStartLoci < 0 || *pStartLoci > *pEndLoci)
	return(-2);

if(OfsLoci == 0 && DeltaLen == 0 && TruncLen == 0)
	return(1 + *pEndLoci -  *pStartLoci);

if(OfsLoci != 0 || DeltaLen != 0)
	{
	if(bOnMinStrand)
		{
		Ofs = OfsLoci * -1;
		DLen = DeltaLen * -1;
		}
	else
		{
		Ofs = OfsLoci;
		DLen = DeltaLen;
		}

	*pStartLoci += Ofs;
	*pEndLoci += Ofs;
	if(bOnMinStrand)
		*pStartLoci += DLen;
	else
		*pEndLoci += DLen;
	}

if(*pStartLoci > *pEndLoci || (!bOnMinStrand && *pStartLoci < 0))
	return(-2);
DLen = 1 + *pEndLoci - *pStartLoci;
if(DLen < MinLen)
	return(-1);
if(DLen > TruncLen)
	{
	if(bOnMinStrand)
		{
		*pStartLoci = 1 + *pEndLoci - TruncLen;
		if(*pStartLoci < 0)
			return(-2);
		}
	else
		*pEndLoci = *pStartLoci + TruncLen - 1;
	DLen= TruncLen;
	}
return(DLen);
}

void
DumpRegionHitStats(void)
{
int Idx;
int RegionHits[4];
int TotRegionHits;
int HitRegionInsts;
int NumRegions;
int ChromID = -1;
tsRegion *pRegion = m_pRegions;
if(!m_NumRegions)
	return;
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Of %d filtering regions there were %d regions containing nucleosome dyads",m_NumRegions,m_NumHitRegions);

ChromID = -1;
for(Idx = 0; Idx <= m_NumRegions; Idx++,pRegion++)
	{
	if(Idx == m_NumRegions || pRegion->ChromID != ChromID)
		{
		if(Idx == m_NumRegions || ChromID != -1)
			{
			m_pBioSeqFile->GetName(ChromID,sizeof(m_szCurChrom),m_szCurChrom);
			gDiagnostics.DiagOutMsgOnly(eDLInfo,
							"Chrom '%s' with %d filtering regions had hits to %d regions totaling %d - hits==1 %d, hits==2 %d, hits==3 %d and hits==4+ %d regions",
							m_szCurChrom,  NumRegions,  HitRegionInsts,    TotRegionHits,RegionHits[0],RegionHits[1],RegionHits[2],RegionHits[3]);
			}
		memset(RegionHits,0,sizeof(RegionHits));
		NumRegions = 0;
		TotRegionHits = 0;
		HitRegionInsts = 0;
		ChromID = pRegion->ChromID;
		}
	NumRegions += 1;
	if(pRegion->NumHits > 0)
		{
		switch(pRegion->NumHits) {
			case 1:
				RegionHits[0] += 1;
				break;
			case 2:
				RegionHits[1] += 1;
				break;
			case 3:
				RegionHits[2] += 1;
				break;
			default:
				RegionHits[3] += 1;
				break;
			}
		HitRegionInsts += 1;
		TotRegionHits += pRegion->NumHits;
		}
	}
}


int
LoadRegions(char *pszRegionFile,	// csv file containing region loci
		int OfsLoci,				// offset region start loci by this many nt
		int DeltaLen,				// change region length by this many nt
		int TruncLength)			// truncate regions to this length
{
int Rslt;
int NumFields;
int NumProcessed;

int Len;
int SrcID;
char *pszElType;
char *pszRefSpecies;
char *pszChrom;
char *pszStrand;
int StartLoci;
int EndLoci;
bool bOnMinStrand;
int NumUnderlen;
int NumStartLess0;

if(pszRegionFile == NULL || pszRegionFile[0] == '\0')
	return(eBSFSuccess);

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Loading regions from CSV loci file '%s'",pszRegionFile);

m_pCSV = new CCSVFile;
if(m_pCSV == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CCSVfile");
	return(eBSFerrObj);
	}

if((Rslt=m_pCSV->Open(pszRegionFile))!=eBSFSuccess)
	{
	while(m_pCSV->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,m_pCSV->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open file: %s",pszRegionFile);
	Reset();
	return(Rslt);
	}

NumUnderlen = 0;
NumStartLess0 = 0;
NumProcessed = 0;
while((Rslt=m_pCSV->NextLine()) > 0)	// onto next line containing fields
	{
	NumFields = m_pCSV->GetCurFields();
	if(NumFields < 7)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expected 7+ fields in '%s', GetCurFields() returned '%d'",pszRegionFile,NumFields);
		Rslt = eBSFerrFieldCnt;
		break;
		}
	if(!NumProcessed && m_pCSV->IsLikelyHeaderLine())
		continue;
	NumProcessed += 1;

	m_pCSV->GetInt(7,&Len);
	m_pCSV->GetInt(1,&SrcID);
	m_pCSV->GetText(2,&pszElType);
	m_pCSV->GetText(3,&pszRefSpecies);
	m_pCSV->GetText(4,&pszChrom);
	m_pCSV->GetInt(5,&StartLoci);
	m_pCSV->GetInt(6,&EndLoci);

	bOnMinStrand = false;
	if(NumFields >= 8)					// check if strand has been specified
		{
		m_pCSV->GetText(8,&pszStrand);
		if(pszStrand[0] == '-')		// assume anything other than '-' is on the plus strand
			bOnMinStrand = true;
		}

	if((Rslt = AdjLoci(bOnMinStrand,		// true if element is on '-' strand
		&StartLoci,		// starting loci on '+' strand
		&EndLoci,			// ending loci on '+' strand
		OfsLoci,			// offset loci by
		DeltaLen,			// adjust length by
		147,				// must be at least this length
		TruncLength)) < 1)	// truncate to this length
		{
		if(Rslt == -1)
			NumUnderlen += 1;
		else
			if(Rslt == -2)
				NumStartLess0 += 1;
		continue;
		}

	if((Rslt=AddRegion(SrcID,pszChrom,StartLoci,EndLoci))<0)
		{
		Reset();
		return(Rslt);
		}
	}
delete m_pCSV;
m_pCSV = NULL;
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Accepted %d regions from %d loaded - %d too short, %d adjusted loci before chrom starts",m_NumRegions,NumProcessed,NumUnderlen,NumStartLess0);
if(m_NumRegions > 1)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Sorting regions...");
	qsort(m_pRegions,m_NumRegions,sizeof(tsRegion),SortRegions);
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Finished sorting regions...");
	}
return(Rslt);
}

// InAnyRegion
// Returns true if sequence at ChromID.StartLoci to EndLoci is completely contained within any loaded filtering region
// Assumes that will be called with ascending start chrom.loci with chrom order the same as that in the loaded sequence file 
bool 
InAnyRegion(int ChromID,int StartLoci,int EndLoci,int NumHits)
{
tsRegion *pProbe;
static int ProbeIdx = 0;
int CurProbeIdx;

if(m_pRegions == NULL || !m_NumRegions)		// if no filtering regions then always return true
	return(true);
if(ProbeIdx == m_NumRegions)
	return(false);

pProbe = &m_pRegions[ProbeIdx];
if(pProbe->ChromID < ChromID)				// if new chrom then iterate forward until matching chrom
	{
	while(++ProbeIdx <= m_NumRegions && pProbe++)
		{
		if(ProbeIdx == m_NumRegions)
			return(false);					

		if(pProbe->ChromID == ChromID)
			break;
		}
	}
// same chromosome, iterate forward until StartLoci >= region StartLoci whilst still on same chromosome
while(ChromID == pProbe->ChromID && StartLoci > pProbe->EndLoci)
	{
	ProbeIdx+=1;
	pProbe+=1;
	}
if(ChromID != pProbe->ChromID)
	return(false);
	
for(CurProbeIdx = ProbeIdx; CurProbeIdx < m_NumRegions; CurProbeIdx++,pProbe++)
	{
	if(pProbe->ChromID != ChromID)
		return(false);
	if(pProbe->StartLoci <= StartLoci && pProbe->EndLoci >= EndLoci)
		{
		if(NumHits)
			{
			if(!pProbe->NumHits)
				m_NumHitRegions += 1;
			pProbe->NumHits += 1;
			}
		return(true);
		}
	if(pProbe->StartLoci > EndLoci)
		return(false);
	}
return(false);
}

double 
Log2(double a)
{
return log(a) / log(2.0);
}

// Simple inplace moving average filter for smoothing
int
MoveAvgFilter(int WinSize,int ValueLen,int *pValues)
{
int WinValues[8196];
int WinIdx;
int Idx;
int Sum;
int *pAvg;

if(ValueLen < 10)
	return(0);

if(WinSize > 8196)
	WinSize = 8196;
else
	if(WinSize < 5)
		WinSize = 5;

if(ValueLen < WinSize)
	WinSize = ValueLen;

Sum = 0;

// leadin..
WinIdx = 0;
pAvg = pValues;
for(Idx = 0; Idx < WinSize/2; Idx++)
	{
	WinValues[WinIdx++] = *pValues;
	Sum += *pValues++;
	}
for(; Idx < WinSize; Idx++)
	{
	*pAvg++ = Sum/Idx;
	WinValues[WinIdx++] = *pValues;
	Sum += *pValues++;
	}
WinIdx %= WinSize;
for(; Idx < ValueLen; Idx++)
	{
	*pAvg++ = Sum/WinSize;
	Sum -= WinValues[WinIdx];
	WinValues[WinIdx++] = *pValues;
	Sum += *pValues++;
	WinIdx %= WinSize;
	}

// lead out
for(Idx = 0; Idx < (WinSize - WinSize/2); Idx++)
	{
	*pAvg++ = Sum/(WinSize - Idx);
	Sum -= WinValues[WinIdx++];
	WinIdx %= WinSize;
	}
return(1);
}

// Simple inplace base line correction
// determines baseline through a wide moving average filter and then subtracts baseline from values
int
BaselineCorrection(int WinSize,int ValueLen,int *pValues)
{
int WinValues[8196];
int WinIdx;
int Idx;
int Sum;
int *pAvg;
int Baseline;

if(ValueLen < 10)
	return(0);

if(WinSize > 8196)
	WinSize = 8196;
else
	if(WinSize < 5)
		WinSize = 5;

if(ValueLen < WinSize)
	WinSize = ValueLen;

Sum = 0;

// leadin..
WinIdx = 0;
pAvg = pValues;
for(Idx = 0; Idx < WinSize/2; Idx++)
	{
	WinValues[WinIdx++] = *pValues;
	Sum += *pValues++;
	}
for(; Idx < WinSize; Idx++)
	{
	Baseline = Sum/Idx;
	if(*pAvg <= Baseline)
		*pAvg++ = 0;
	else
		*pAvg++ -= Baseline;
	WinValues[WinIdx++] = *pValues;
	Sum += *pValues++;
	}
WinIdx %= WinSize;
for(; Idx < ValueLen; Idx++)
	{
	Baseline = Sum/WinSize;
	if(*pAvg <= Baseline)
		*pAvg++ = 0;
	else
		*pAvg++ -= Baseline;
	Sum -= WinValues[WinIdx];
	WinValues[WinIdx++] = *pValues;
	Sum += *pValues++;
	WinIdx %= WinSize;
	}

// lead out
for(Idx = 0; Idx < (WinSize - WinSize/2); Idx++)
	{
	Baseline = Sum/(WinSize - Idx);
	if(*pAvg <= Baseline)
		*pAvg++ = 0;
	else
		*pAvg++ -= Baseline;
	Sum -= WinValues[WinIdx++];
	WinIdx %= WinSize;
	}
return(1);
}

int
OutputDyads(int hFile,					// results file
			etFMode FMode,				// output format mode
			int MovAvgFilt,				// apply this width moving average filter
			int BaselineFilter,			// use this window size when normalising for baseline
			char *pszTrackName,			// track name
			char *pszCurChromName,		// dyads on this chromosome
			int MaxDyadLoci,			// last dyad loci
			int *pDyadScores)			// dyad scores
{
int Rslt;
char szBuff[4096];
int BuffIdx;
int CurDyadLoci;
int StartDyadLoci;
int CurDyadScore;
int PeakScore;
BuffIdx = 0;
int InitialPutDyads =  m_NumPutDyads;

if(MovAvgFilt > 0)
	MoveAvgFilter(MovAvgFilt,MaxDyadLoci,pDyadScores);

if(BaselineFilter > 0)
	BaselineCorrection(BaselineFilter,MaxDyadLoci,pDyadScores);
	
if(FMode == eFMMcsvScores)
	{
	for(CurDyadLoci = 0; CurDyadLoci < MaxDyadLoci; CurDyadLoci += 1)
		{
		CurDyadScore = pDyadScores[CurDyadLoci];

		BuffIdx += sprintf(&szBuff[BuffIdx],"%s,%d,%d\n",pszCurChromName,CurDyadLoci,CurDyadScore);
		if((BuffIdx + 500) > sizeof(szBuff))
			{
			if((Rslt=write(hFile,szBuff,BuffIdx))!=BuffIdx)
				return(eBSFerrWrite);
			BuffIdx = 0;
			}
		}
	if((BuffIdx + 500) > sizeof(szBuff))
		{
		if((Rslt=write(hFile,szBuff,BuffIdx))!=BuffIdx)
			return(eBSFerrWrite);
		BuffIdx = 0;
		}
	return(eBSFSuccess);
	}

StartDyadLoci = -1;
PeakScore = 0;
for(CurDyadLoci = 0; CurDyadLoci < MaxDyadLoci; CurDyadLoci += 1)
	{
	CurDyadScore = pDyadScores[CurDyadLoci];
	if(CurDyadScore <= 0 && StartDyadLoci == -1)
		continue;
	if(CurDyadScore >= PeakScore) 
		{
		if(CurDyadScore > PeakScore)
			{
			StartDyadLoci = CurDyadLoci;
			PeakScore = CurDyadScore;
			}
		continue;
		}
	// past local peak?
	if(StartDyadLoci != -1)
		{
		PeakScore = (int)(10.0 * Log2((double)PeakScore));
		if(StartDyadLoci > 73)
			{
			m_NumPutDyads += 1;
			switch(FMode) {
				case eFMbedDyads:
					BuffIdx += sprintf(&szBuff[BuffIdx],"%s\t%d\t%d\t%s\t%d\n",pszCurChromName,StartDyadLoci,CurDyadLoci,pszTrackName,PeakScore);
					break;	
				case eFMbedGraphDyads:
					BuffIdx += sprintf(&szBuff[BuffIdx],"%s\t%d\t%d\t%d\n",pszCurChromName,StartDyadLoci,CurDyadLoci,PeakScore);
					break;
				case eFMcsvDyads:
					BuffIdx += sprintf(&szBuff[BuffIdx],"%d,Dyad,%s,%s,%d,%d,%d,%d\n",m_NumPutDyads,pszTrackName,pszCurChromName,StartDyadLoci,CurDyadLoci-1,CurDyadLoci-StartDyadLoci,PeakScore);
					break;

				case eFMbedNucs:
					BuffIdx += sprintf(&szBuff[BuffIdx],"%s\t%d\t%d\t%s\t%d\n",pszCurChromName,StartDyadLoci-73,CurDyadLoci+73,pszTrackName,PeakScore);
					break;
				case eFMbedGraphNucs:
					BuffIdx += sprintf(&szBuff[BuffIdx],"%s\t%d\t%d\t%d\n",pszCurChromName,StartDyadLoci-73,CurDyadLoci+73,PeakScore);
					break;
				case eFMcsvNucs:
					BuffIdx += sprintf(&szBuff[BuffIdx],"%d,Nucleosome,%s,%s,%d,%d,%d,%d\n",m_NumPutDyads,pszTrackName,pszCurChromName,StartDyadLoci-73,CurDyadLoci+72,146+CurDyadLoci-StartDyadLoci,PeakScore);
					break;
				}
			if((BuffIdx + 500) > sizeof(szBuff))
				{
				if((Rslt=write(hFile,szBuff,BuffIdx))!=BuffIdx)
					return(eBSFerrWrite);
				BuffIdx = 0;
				}
			}
		StartDyadLoci = -1;
		}
	PeakScore = (CurDyadScore * 110) / 100;			// allow 10% so as to provide some hysteresis
	}

if(StartDyadLoci != -1)
	{
	PeakScore = (int)(10.0 * Log2((double)PeakScore));
	m_NumPutDyads += 1;
	switch(FMode) {
		case eFMbedDyads:
			BuffIdx += sprintf(&szBuff[BuffIdx],"%s\t%d\t%d\t%s\t%d\n",pszCurChromName,StartDyadLoci,CurDyadLoci,pszTrackName,PeakScore);
			break;	
		case eFMbedGraphDyads:
			BuffIdx += sprintf(&szBuff[BuffIdx],"%s\t%d\t%d\t%d\n",pszCurChromName,StartDyadLoci,CurDyadLoci,PeakScore);
			break;
		case eFMcsvDyads:
			BuffIdx += sprintf(&szBuff[BuffIdx],"%d,Dyad,%s,%s,%d,%d,%d,%d\n",m_NumPutDyads,pszTrackName,pszCurChromName,StartDyadLoci,CurDyadLoci-1,CurDyadLoci-StartDyadLoci,PeakScore);
			break;

		case eFMbedNucs:
			BuffIdx += sprintf(&szBuff[BuffIdx],"%s\t%d\t%d\t%s\t%d\n",pszCurChromName,StartDyadLoci-73,CurDyadLoci+73,pszTrackName,PeakScore);
			break;
		case eFMbedGraphNucs:
			BuffIdx += sprintf(&szBuff[BuffIdx],"%s\t%d\t%d\t%d\n",pszCurChromName,StartDyadLoci-73,CurDyadLoci+73,PeakScore);
			break;
		case eFMcsvNucs:
			BuffIdx += sprintf(&szBuff[BuffIdx],"%d,Nucleosome,%s,%s,%d,%d,%d,%d\n",m_NumPutDyads,pszTrackName,pszCurChromName,StartDyadLoci-73,CurDyadLoci+72,146+CurDyadLoci-StartDyadLoci,PeakScore);
			break;
		}		
	}

if(BuffIdx)
	if((Rslt=write(hFile,szBuff,BuffIdx))!=BuffIdx)
		return(eBSFerrWrite);
m_ChromPutDyads = m_NumPutDyads - InitialPutDyads; 
return(eBSFSuccess);
}

int
Process(etPMode PMode,					// processing mode
		etFMode FMode,					// output format mode
		int MovAvgFilter,				// apply this moving average window size filter
		int BaselineFilter,				// baseline normalisation window size
		char *pszTrackName,				// UCSC track name
		double DyadratioThres,			// dyad grooves must be at least this ratio to background
		double Dyad2ratioThres,			// immediately flanking grooves must be at least this ratio to background
		double Dyad3ratioThres,			// remainder of flanking grooves must be at least this ration to background
		char *pszInGenomefile,			// bioseq genome file
		char *pszInConfFile,			// file containing conformation characteristics
		char *pszOutFile,				// where to write nucleosome predictions
		char *pszInclRegionFile,		// only report predicted nucleosomes if intersecting with regions in this file
		int OfsLoci,				// offset region start loci by this many nt
		int DeltaLen,				// change region length by this many nt
		int TruncLength)			// truncate regions to be a maximum of this length
{
int Rslt = 0;
bool bRegionFilter;		// set true if region filtering
int SeqIdx;
int WindLen = 147;		// nucleosomes are this length - agreed :-)
int ChkGroove[13];		// to hold dyad (ChkGroove[6]) and +/- 6 at approx decimer (depends on twist) offsets 
int DecIdx;				// index into ChkGroove, incr/decr every 360 degree twist
int GrooveCnt;			// to hold number of groove values contributing to current ChkGroove[DecIdx] so average can be calculated
int AccumTwist;			// to hold accumulated twist relative to dyad
int ChkTwist;			// AccumTwist % 360 used to determine if minor groove back on same plane as at dyad
int BaseLineValueSum;
double BaseLineAv;
double	ChromBaseLineSum;
int NumBaseLines;



Init();
m_PMode = PMode;	
m_FMode = FMode;	
strcpy(m_szTrackName,pszTrackName);
m_DyadratioThres = DyadratioThres;
m_Dyad2ratioThres = Dyad2ratioThres;
m_Dyad3ratioThres = Dyad3ratioThres;

if((m_pszOutBuff = new char [cAllocOutBuff])==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to alloc memory (%d bytes requested) for bufering output results",cAllocOutBuff);
	Reset();
	return(eBSFerrMem);
	}
m_AllocdOutBuff = cAllocOutBuff;

if((m_pBioSeqFile = new CBioSeqFile) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CBioSeqFile object");
	Reset();
	return(eBSFerrObj);
	}

if((Rslt = m_pBioSeqFile->Open(pszInGenomefile))!=eBSFSuccess)
	{
	while(m_pBioSeqFile->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,m_pBioSeqFile->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open genome assembly sequence file '%s'",pszInGenomefile);
	Reset();
	return(Rslt);
	}

if(pszInclRegionFile != NULL && pszInclRegionFile[0] != '\0')
	{
	if((Rslt = LoadRegions(pszInclRegionFile,OfsLoci,DeltaLen,TruncLength))!=eBSFSuccess)
		{
		Reset();
		return(Rslt);
		}
	bRegionFilter = true;
	}
else
	bRegionFilter = false;

if((m_pTwister = new CTwister)==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,"ProcessFastaStruct","Unable to create CTwister object");
	Reset();
	return(eBSFerrObj);
	}

if((Rslt = m_pTwister->LoadStructParams(pszInConfFile))  < eBSFSuccess)
	{
	while(m_pTwister->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,m_pTwister->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,"ProcessFastaStruct","LoadStructParams(%s) failed",pszInConfFile);
	Reset();
	return(Rslt);
	}

#ifdef _WIN32
if((m_hOutFile = open(pszOutFile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
if((m_hOutFile = open(pszOutFile, O_RDWR | O_CREAT | O_TRUNC, S_IREAD | S_IWRITE ))==-1)
#endif
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create predicted nucleosome output file: %s - %s",pszOutFile,strerror(errno));
	Reset();
	return(eBSFerrCreateFile);
	}


if(FMode == eFMbedGraphDyads || FMode == eFMbedGraphNucs)
	{
	m_UsedOutBuff = sprintf(m_pszOutBuff,
	 "track type=bedGraph name=\"%s\" description=\"%s\" visibility=full color=200,100,0 altColor=0,100,200 priority=20 autoScale=on alwaysZero=on graphType=bar smoothingWindow=4\n",
		pszTrackName,pszTrackName);
	CUtility::SafeWrite(m_hOutFile,m_pszOutBuff,m_UsedOutBuff);
	m_UsedOutBuff = 0;
	}

// iterate over chromosomes
tBSFEntryID ChromID = 0;
while((ChromID = m_pBioSeqFile->Next(ChromID))>0)
	{
	m_pBioSeqFile->GetName(ChromID,sizeof(m_szCurChrom),m_szCurChrom);
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing %s...",m_szCurChrom);
	m_ChromSeqLen = m_pBioSeqFile->GetDataLen(ChromID);
	if(m_pChromSeq == NULL || m_ChromSeqLen > m_AllocdChromSeq)
		{
		if(m_pChromSeq != NULL)
			{
			delete m_pChromSeq;
			m_pChromSeq = NULL;
			}
		int AllocLen = m_ChromSeqLen + m_ChromSeqLen/10; 
		if((m_pChromSeq = new unsigned char [AllocLen])==NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to alloc memory (%d requested) for holding raw sequence data",AllocLen);
			Reset();
			return(eBSFerrMem);
			}
		m_AllocdChromSeq = AllocLen;

		if(m_pScores != NULL)
			{
			delete m_pScores;
			m_pScores = NULL;
			}
		if((m_pScores = new int [AllocLen])==NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate memory (%d bytes) for holding scores",AllocLen);
			Reset();
			return(eBSFerrMem);
			}

		if(m_pConfGroove != NULL)
			{
			delete m_pConfGroove;
			m_pConfGroove = NULL;
			}
		if((m_pConfGroove = new int [AllocLen])==NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate memory (%d bytes) for holding minor groove conformation values",AllocLen);
			Reset();
			return(eBSFerrMem);
			}
		if(m_pConfTwist != NULL)
			{
			delete m_pConfTwist;
			m_pConfTwist = NULL;
			}
		if((m_pConfTwist = new int [AllocLen])==NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate memory (%d bytes) for holding rotational twist conformation values",AllocLen);
			Reset();
			return(eBSFerrMem);
			}
		}

	if((Rslt=m_pBioSeqFile->GetData(ChromID,eSeqBaseType,0,m_pChromSeq,m_ChromSeqLen)) != m_ChromSeqLen)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Loading sequence of length %d failed from chrom: %s file: '%s'",m_ChromSeqLen,m_szCurChrom,pszInGenomefile);
		Reset();
		return(Rslt);
		}

	// remove any repeat masking and randomly substitute bases for eBaseN's - not expecting too many of these say's he hopefully!
	etSeqBase *pSeq = m_pChromSeq;
	for(SeqIdx = 0; SeqIdx < m_ChromSeqLen; SeqIdx++,pSeq++)
		if((*pSeq &= ~cRptMskFlg) > eBaseT)
			*pSeq = rand() % 4;

	if((Rslt = m_pTwister->GetSequenceConformation(eSSminorgroove,	// process for this conformational parameter
				  0,						// initial starting offset (0..n) in pSeq
				  0,						// number of steps (0 for all) to process starting at pSeq[iStartPsn]|pSeq[iStartPsn+1]
				  m_ChromSeqLen,			// total length of sequence
				  m_pChromSeq,				// sequence to be processed
				  m_pConfGroove))!=eBSFSuccess) // where to return step conformational values
					{
					gDiagnostics.DiagOut(eDLFatal,"GetSequenceConformation","minor groove failed");
					Reset();
					return(Rslt);
					}

	if((Rslt = m_pTwister->GetSequenceConformation(eSStwist,	// process for this conformational parameter
				  0,						// initial starting offset (0..n) in pSeq
				  0,						// number of steps (0 for all) to process starting at pSeq[iStartPsn]|pSeq[iStartPsn+1]
				  m_ChromSeqLen,			// total length of sequence
				  m_pChromSeq,				// sequence to be processed
				  m_pConfTwist))!=eBSFSuccess) // where to return step conformational values
					{
					gDiagnostics.DiagOut(eDLFatal,"GetSequenceConformation","rotational twist failed");
					Reset();
					return(Rslt);
					}

	memset(m_pScores,0,m_ChromSeqLen * sizeof(int));	// reset all scores back to minimum

	// establish the baseline conformational characteristic value over initial window which is centered around the WindLen window
	// the baseline will be updated as the putative dyad is slid along the chromosome
	int *pConfGroove;
	int *pBaseLineWin5;
	int *pBaseLineWin3;
	int *pConfTwist;
	int BaseLineWin = min(5 * WindLen,m_ChromSeqLen);
	BaseLineValueSum = 0;
	pBaseLineWin3 = m_pConfGroove;
	for(SeqIdx = 0; SeqIdx < BaseLineWin; SeqIdx++)
		BaseLineValueSum += *pBaseLineWin3++;

	double BestMaxRatio = 0.0f;
	double BestMinRatio = 10000.0f;
	int DyadFirstOfs = 73;
	int DyadLastOfs = m_ChromSeqLen - 73;
	m_ChromPutDyads = 0;
	ChromBaseLineSum = 0.0;
	NumBaseLines = 0;
	pBaseLineWin5 = m_pConfGroove;
	for(SeqIdx = DyadFirstOfs; SeqIdx < DyadLastOfs; SeqIdx++)
		{
		if(SeqIdx > BaseLineWin/2 && SeqIdx < (m_ChromSeqLen - (BaseLineWin+1)/2))
			{
			BaseLineValueSum -= *pBaseLineWin5++;
			BaseLineValueSum += *pBaseLineWin3++;
			}
		BaseLineAv = (double)BaseLineValueSum/BaseLineWin;
		ChromBaseLineSum += BaseLineAv/10000.0f;
		NumBaseLines += 1;
		if(bRegionFilter && !InAnyRegion(ChromID,SeqIdx-74,SeqIdx+73,0))
			continue;

		DecIdx = 6;
		pConfGroove =  &m_pConfGroove[SeqIdx];
		ChkGroove[DecIdx++] = *pConfGroove++;
		m_DyadRatio = (double)ChkGroove[6]/BaseLineAv;
		if(m_DyadRatio < DyadratioThres)
			continue;
		pConfTwist =  &m_pConfTwist[SeqIdx+1];
		AccumTwist = *pConfTwist;
		ChkGroove[DecIdx] = 0;
		GrooveCnt = 0;
		
		// iterate over bases to right of putative dyad and every rotation of the dsDNA get the minor groove
		int Bases = 1;
		while(DecIdx <= 12)
			{
			Bases += 1;
			pConfTwist += 1;
			pConfGroove += 1;
			AccumTwist += *pConfTwist;
			ChkTwist = AccumTwist % 3600000;
			if(ChkTwist >= 3300000 || ChkTwist <= 300000)
				{
				ChkGroove[DecIdx] += *pConfGroove;
				GrooveCnt += 1;
				}
			else
				{
				if(GrooveCnt > 0)
					{
					ChkGroove[DecIdx] /= GrooveCnt;
					GrooveCnt = 0;
					if(DecIdx++ < 12)
						ChkGroove[DecIdx]= 0;
					}
				}
			}
		// now iterate over bases to left of putative dyad and every rotation of the dsDNA get the minor groove
		DecIdx = 5;
		pConfGroove =  &m_pConfGroove[SeqIdx-1];
		pConfTwist =  &m_pConfTwist[SeqIdx-1];
		AccumTwist = *pConfTwist;
		ChkGroove[DecIdx] = 0;
		GrooveCnt = 0;
		Bases = 1;
		while(DecIdx >= 0)
			{
			Bases += 1;
			pConfTwist -= 1;
			pConfGroove -= 1;
			AccumTwist += *pConfTwist;
			ChkTwist = AccumTwist % 3600000;
			if(ChkTwist >= 3300000 || ChkTwist <= 300000)
				{
				ChkGroove[DecIdx] += *pConfGroove;
				GrooveCnt += 1;
				}
			else
				{
				if(GrooveCnt > 0)
					{
					ChkGroove[DecIdx] /= GrooveCnt;
					GrooveCnt = 0;
					if(DecIdx-- > 0)
						ChkGroove[DecIdx] = 0;
					}
				}
			}

		m_Dyad2Ratio = (double)(ChkGroove[5] + ChkGroove[7])/(2*BaseLineAv);
		m_Dyad3Ratio = (double)(ChkGroove[0] + ChkGroove[1] + ChkGroove[2] + ChkGroove[3] + ChkGroove[4] +
					ChkGroove[8] + ChkGroove[9] + ChkGroove[10] + ChkGroove[11] + ChkGroove[12])/(10*BaseLineAv);

		
		if(m_Dyad2Ratio < Dyad2ratioThres || m_Dyad3Ratio < Dyad3ratioThres)
			continue;

		if(bRegionFilter)
			InAnyRegion(ChromID,SeqIdx-74,SeqIdx+73,1);

		int LocScore = (int)(1000 * ((m_DyadRatio - 1.0f) + ((m_Dyad2Ratio - 1.0f) * 0.85) + ((m_Dyad3Ratio - 1.0f) * 0.75)));
		if(LocScore < 0)
			printf("\nHave an issue!");
		m_pScores[SeqIdx] = LocScore;
		}

	OutputDyads(m_hOutFile,				// results file
			FMode,						// output format mode
			MovAvgFilter,				// apply this width moving average filter
			BaselineFilter,				// use this window size when normalising for baseline
			pszTrackName,				// track name
			m_szCurChrom,				// dyads on this chromosome
			m_ChromSeqLen,				// last dyad loci
			m_pScores);					// dyad scores
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"On chromosome '%s' (length %d), %d putative nucleosome dyads were identified - running total: %d",m_szCurChrom,m_ChromSeqLen,m_ChromPutDyads,m_NumPutDyads);
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"On chromosome '%s' (length %d), averaged minor groove baseline was %1.6f",m_szCurChrom,m_ChromSeqLen,ChromBaseLineSum/NumBaseLines);
	}
if(m_UsedOutBuff)
	{
	if((Rslt=write(m_hOutFile,m_pszOutBuff,m_UsedOutBuff))!=m_UsedOutBuff)
		Rslt = eBSFerrWrite;
	m_UsedOutBuff = 0;
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"After processing a total of %d putative nucleosome dyads have been identified",m_NumPutDyads);
DumpRegionHitStats();
Reset();
return(Rslt);
}


// SortRegions
// sort regions by ChromID-->StartLoci-->EndLoci ascending
int SortRegions( const void *arg1, const void *arg2)
{
tsRegion *pR1 = (tsRegion *)arg1;
tsRegion *pR2 = (tsRegion *)arg2;
if(pR1->ChromID < pR2->ChromID)
	return(-1);
if(pR1->ChromID > pR2->ChromID)
	return(1);
if(pR1->StartLoci < pR2->StartLoci)
	return(-1);
if(pR1->StartLoci > pR2->StartLoci)
	return(1);
return(pR1->EndLoci - pR2->EndLoci);
}