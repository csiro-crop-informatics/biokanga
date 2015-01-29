// proccentroids.cpp : Defines the entry point for the console application.
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

const unsigned int cProgVer = 300;		// increment with each release

const int cMaxOligoLen = 20;				// max length oligo supported
const int cMaxRegions = 7;					// IG/US5/UTR5/CDS/Intron/UTR3/DS3
const int cMaxChroms  = 500;				// supports at most this many chromosomes
const int cLineBuffSize = 32000;			// max source CVS file line length
const int cNumOligos2Alloc = 0x0ffff;		// block (number of oligos) size to prealloc memory for

typedef struct TAG_sRegionStats {
	int BaseA;								// centroid base counts
	int BaseC;
	int BaseG;
	int BaseT;
	double TransProbMatrix[4][4];			// transistional probability matrix
	double ChiSqrA;							// Chi-square for this oligo against all other oligos with same centroid irrespective of region 
	double ChiSqr;							// Chi-square for this oligo against all other oligos in same region
	double FMutA;							// fishers exact for mutation rates
	double FMut;							// fishers exact for mutation rates
	double FTransA;	     					// fishers exact for transistion vs transversion ratios
	double FTrans;	     					// fishers exact for transistion vs transversion ratios
	double ProbFixed;						// prob of this centroid not mutating and the mutation becoming fixed
	double ProbTransistion;					// when centroid is mutated/fixed what is prob of it being a transition		
	double ProbTransversion;				// when centroid is mutated/fixed what is prob of it being a transversion		
} tsRegionStats;

typedef struct TAG_sCentOligo {
	char szSeq[cMaxOligoLen];
	char szChrom[cMaxDatasetSpeciesChrom];
	int ChromID;
	char Centroid;							// centroid base in this oligo
	char szCentroid3[4];					// centroid plus each flanking base in this oligo
	int Offset;
	int SeqID;
	int NumInstances;
	tsRegionStats Region[cMaxRegions];
} tsCentOligo;

typedef struct TAG_sCChrom {
	char szChrom[cMaxDatasetSpeciesChrom];
	int ChromID;
	} tsCChrom;

typedef enum {
	eFTSpecies = 0,				// file contains genome specific oligo distribution counts
	eFTAlign,					// file contains alignment oligo distribution counts
	eFTTransMatrix,				// same as eFTAlign for input file, but transitional matrix probabilities will be generated
	eFTendmark					// used to mark end of enumeration
} teInterFileType;

bool ParseCSV(char *pszIntermediate,	// CSV file containing raw counts
		 teInterFileType FileType,	// intermediate CSV file type
		 bool bRegionalCnts);  	// true if file expected to contain regional countsbool ProcessRows(int NMer);
bool CalcAllStats(int Mode,int NMer,int NumOligos,tsCentOligo *pOligos);
int GetChromID(char *pszChrom); // returns unique identifier for this chromosome
bool ProcessRows(int Mode,int NMer);				// oligo length
bool OutputStats(int iMode,char *pszRsltsFile);
bool
StationaryProbs(int NumPeriods,		// number of time periods
				teSeqBases InitialBase,	// initially observed base
				double *pProbA,		// where to return probabilities in each time period for A
				double *pProbC,		// where to return probabilities in each time period for C
				double *pProbG,		// where to return probabilities in each time period for G
				double *pProbT,		// where to return probabilities in each time period for T
				double *pMatrix);	// 4*4 transitional probabilites


int m_CurNumOligos;		// current number of oligos used
int m_NumOligosAllocd;	// current number of centroids allocd
tsCentOligo *m_pAllocdOligos;	// memory allocd to hold centroids rows
int m_NumChroms;				// current number of chromsomes in m_Chroms[]
tsCChrom m_Chroms[cMaxChroms];	// to hold chromsome names

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


char szOutputFile[_MAX_PATH];
char szInputFile[_MAX_PATH];
int iNMer;
int iMode;
bool bRegionStats;

// command line args
struct arg_lit  *help    = arg_lit0("hH","help",                "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_int *ScreenLogLevel=arg_int0("S", "ScreenLogLevel",	"<int>","Level of diagnostics written to logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",	"diagnostics log file");

struct arg_int  *NMer = arg_int0("n","nmer","<int>",			"expected oligo length in CVS intermediate results file, must be 1,3,5,7,9,11");
struct arg_file *InFile = arg_file1("i",NULL,"<file>",			"CSV intermediate results file containing counts");
struct arg_file *OutFile = arg_file1("o",NULL,"<file>",			"output processed stats to results file as CSV");
struct arg_int  *Mode = arg_int0("m","mode","<int>",			"processing mode - 0 (default) == genome cnts, 1 == alignment cnts, 2 == transitional matrix probs, 3 == stationary probabilities");
struct arg_lit  *RegionStats = arg_lit0("r","region",			"generate regional stats in addition to global");
struct arg_end *end = arg_end(20);

void *argtable[] = {help,version,FileLogLevel,ScreenLogLevel,LogFile,NMer,InFile,OutFile,RegionStats,Mode,end};

char **pAllArgs;
int argerrors;
argerrors = CUtility::arg_parsefromfile(argc,(char **)argv,&pAllArgs);
if(argerrors >= 0)
	argerrors = arg_parse(argerrors,pAllArgs,argtable);

    /* special case: '--help' takes precedence over error reporting */
if (help->count > 0)
        {
		printf("Usage: %s ", gszProcName);
        arg_print_syntax(stdout,argtable,"\n");
        arg_print_glossary(stdout,argtable,"  %-25s %s\n");
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


	iMode = Mode->count ? Mode->ival[0] : 0;
	if(iMode > eFTendmark || iMode < 0)
		{
		printf("\nError: Processing mode '-m%d' specified, mode must be 0..2",iMode);
		exit(1);
		}
	bRegionStats = RegionStats->count ? true : false;

	iNMer = NMer->count ? NMer->ival[0] : 3;
	if(iNMer < 1)
		iNMer = 1;

	strcpy(szInputFile,InFile->filename[0]);
	strcpy(szOutputFile,OutFile->filename[0]);

	m_CurNumOligos = 0;
	m_NumOligosAllocd = 0;
	m_pAllocdOligos = NULL;

		// now that command parameters have been parsed then initialise diagnostics log system
	if(!gDiagnostics.Open(szLogFile,(etDiagLevel)iScreenLogLevel,(etDiagLevel)iFileLogLevel,true))
		{
		printf("\nError: Unable to start diagnostics subsystem.");
		if(szLogFile[0] != '\0')
			printf(" Most likely cause is that logfile '%s' can't be opened/created",szLogFile);
		exit(1);
		}

		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Version: %d.%2.2d Processing parameters:",cProgVer/100,cProgVer%100);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Mode: %d",iMode);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Regional Stats: %s",bRegionStats ? "yes" : "no");
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"CSV intermediate results file containing counts: '%s'",szInputFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"expected oligo length in CVS intermediate results file: %d",iNMer);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"where to write out stats: '%s'",szOutputFile);

	gStopWatch.Start();

#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif

	bool bRslt;
	if(bRslt = ParseCSV(szInputFile,(teInterFileType)iMode,bRegionStats))
		{
		if(bRslt = ProcessRows(iMode,iNMer))
			bRslt = OutputStats(iMode,szOutputFile);
		}

	if(m_pAllocdOligos)
		delete m_pAllocdOligos;
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Exit Code: %d Total processing time: %s",bRslt ? 0 : 1,gStopWatch.Read());
	exit(bRslt ? 0 : 1);
	}
else
	{
	arg_print_errors(stdout,end,gszProcName);
	arg_print_syntax(stdout,argtable,"\nend of help\n");
	exit(1);
	}
return 0;
}

int 
GetChromID(char *pszChrom)
{
int Cnt;
tsCChrom *pChrom;
char Chrom1 = tolower(*pszChrom);
pChrom = m_Chroms;
for(Cnt = 0; Cnt < m_NumChroms; Cnt++,pChrom++)
	{
	if(tolower(pChrom->szChrom[0]) == Chrom1 &&
		!stricmp(pChrom->szChrom,pszChrom))
			return(pChrom->ChromID);
	}
if(m_NumChroms == cMaxChroms)
	return(0);
strcpy(pChrom->szChrom,pszChrom);
pChrom->ChromID = ++m_NumChroms;
return(pChrom->ChromID);
}

// ProcessRows
// The intermediate output file containing the CSV data has been parsed
// Time to process this data..
bool
ProcessRows(int Mode,				// processing mode
			int NMer)				// oligo length
{
return(CalcAllStats(Mode,			// processing mode
					NMer,			// oligo length
			     m_CurNumOligos,	// number of oligos
		      m_pAllocdOligos));    // array of oligos

}

tsCentOligo *
AllocCentroids(int NumReq)
{
tsCentOligo *pNewCentroids;
pNewCentroids = new tsCentOligo [m_NumOligosAllocd + NumReq];
if(pNewCentroids == NULL)
	return(NULL);
if(m_pAllocdOligos != NULL)
	{
	memcpy(pNewCentroids,m_pAllocdOligos,m_NumOligosAllocd * sizeof(tsCentOligo));
	delete m_pAllocdOligos;
	}
m_pAllocdOligos = pNewCentroids;
memset(&m_pAllocdOligos[m_NumOligosAllocd],0,sizeof(tsCentOligo) * NumReq);
m_NumOligosAllocd += NumReq;
return(m_pAllocdOligos);
}

// ParseCSV
// Open and parse the intermediate results file into centroid counts structure which
// can then be post-processed to provide meaningful stats
// If first CSV line contains column titles then this is simply sloughed
//
bool
ParseCSV(char *pszIntermediate,		// CSV file containing raw counts
		 teInterFileType FileType,	// CSV file type
		 bool bRegionalCnts)		// true if file expected to contain regional counts
{
tsRegionStats *pRegion;
tsCentOligo *pCurOligo;
int NumChrs;
int LineLen;
FILE *pInterStream;
char szLineBuffer[cLineBuffSize];
int Cnt;
int Region;
int LineOfs;
char *pSrc;
char *pDst;
char Chr;
bool bHeader = true;		// expecting header
double dbltoload = 0.01;

if(m_pAllocdOligos == NULL)
	if(!AllocCentroids(cNumOligos2Alloc))
		return(false);

if((pInterStream = fopen(pszIntermediate,"r")) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"ParseCSV: Unable to open %s - %s",pszIntermediate,strerror(errno));
	delete m_pAllocdOligos;
	m_pAllocdOligos = NULL;
	return(false);
	}

m_NumChroms = 0;
while(fgets(szLineBuffer,sizeof(szLineBuffer),pInterStream)!=NULL)
	{
	if(strlen(szLineBuffer) < 5)	// simply slough lines which are too short to contain anything worth parsing
		continue;

	// check if more memory needs to be allocated to hold new row
	if(m_CurNumOligos == m_NumOligosAllocd)
		if(AllocCentroids(cNumOligos2Alloc)==NULL)
			{
			if(m_pAllocdOligos != NULL)
				delete m_pAllocdOligos;
			m_pAllocdOligos = NULL;
			fclose(pInterStream);
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"ParseCSV: Unable to allocate memory");
			return(false);
			}

	// strip out quotes, whitespace, makes it a lot easier on the parsing!
	pSrc = szLineBuffer;
	pDst = pSrc;
	while(Chr=*pSrc++)
		if(!(Chr == '\"' || Chr == '\'' || isspace(Chr)))
			*pDst++ = Chr;
	*pDst = '\0';
	if(!(LineLen = (int)strlen(szLineBuffer)))	// simply slough blank lines
		continue;

	pCurOligo = &m_pAllocdOligos[m_CurNumOligos];
	Cnt = sscanf(szLineBuffer," %25[^,],%d,%d,%c,%[^,],%[^,],%d %n",
			 &pCurOligo->szChrom[0],&pCurOligo->Offset,&pCurOligo->SeqID,&pCurOligo->Centroid,&pCurOligo->szCentroid3[0],&pCurOligo->szSeq[0],&pCurOligo->NumInstances,&NumChrs);
	if(Cnt < 7)		
		{
		if(bHeader)					// may be less than expected because this is a header row
			{
			bHeader = false;		// allow a single header instance
			continue;
			}
		if(m_pAllocdOligos != NULL)
			delete m_pAllocdOligos;
		m_pAllocdOligos = NULL;
		fclose(pInterStream);
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"ParseCSV: Unable to parse %s, unexpected layout",pszIntermediate);
		return(false);
		}
	if((pCurOligo->ChromID = GetChromID(pCurOligo->szChrom)) < 1)
		{
		if(m_pAllocdOligos != NULL)
			delete m_pAllocdOligos;
		m_pAllocdOligos = NULL;
		fclose(pInterStream);
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"ParseCSV: Too many chromosomes in %s",pszIntermediate);
		return(false);
		}

	if(bRegionalCnts)
		{
		LineOfs = NumChrs;
		pRegion = pCurOligo->Region;
		switch(FileType) {
			case eFTSpecies:						// species counts
				for(Region=0;Region<7;Region++,pRegion++)
					{
					Cnt = sscanf(&szLineBuffer[LineOfs],",%d %n",
						&pRegion->BaseA, &NumChrs);
					LineOfs += NumChrs;
					}
				break;
			case eFTAlign:						// aligned species counts
			case eFTTransMatrix:
				for(Region=0;Region<7;Region++,pRegion++)
					{
					Cnt = sscanf(&szLineBuffer[LineOfs],",%d,%d,%d,%d %n",
						&pRegion->BaseA,&pRegion->BaseC,&pRegion->BaseG,&pRegion->BaseT,&NumChrs);
					LineOfs += NumChrs;
					}
				break;
			}
		}
	bHeader = false;
	m_CurNumOligos++;
	}
fclose(pInterStream);
return(true);
}


// CalcFixProbs
// Probabilites that a given oligo will mutate globally and by region
// row generated has
// <seq> <global> <region>...<region>

// CalcTransProbs
// Probabilites that a given oligo will mutate transition vs transversion globally and by region
// <seq> <global> <region>...<region>


// CalcAllStats
// Calculate stats
bool	
CalcAllStats(int Mode,			// processing mode
			 int NMer,			// oligo length
			   int NumOligos,		// number of oligos
			   tsCentOligo *pOligos) // array of oligos
{
int CentIdx;
int SumACGT;
tsRegionStats *pRegion;
tsCentOligo *pOligo;
int Steps;
int Idx;
int CentMsk;

if(NMer < 3 || NMer > 7)
	return(false);

if(Mode == eFTTransMatrix)
	{
	int OligoOfs;
	tsCentOligo *pSrcOligo;
	tsRegionStats *pSrcRegion;
	switch(NMer) {
		case 3:
			OligoOfs = 0x0004;
			CentMsk =  0x000c;
			break;
		case 5:
			OligoOfs = 0x0010;
			CentMsk =  0x0030;
			break;
		case 7:
			OligoOfs = 0x0040;
			CentMsk =  0x00c0;
			break;
		}

	pOligo = pOligos;
	for(CentIdx = 0; CentIdx < NumOligos; CentIdx++, pOligo++)
		{
		pRegion = pOligo->Region;
		for(Steps = 0; Steps < cMaxRegions; Steps++,pRegion++)
			{
			pSrcOligo = &pOligos[CentIdx & ~CentMsk];
			for(Idx = 0; Idx < 4; Idx++, pSrcOligo += OligoOfs)
				{
				pSrcRegion = &pSrcOligo->Region[Steps];
				SumACGT = pSrcRegion->BaseA + pSrcRegion->BaseC + pSrcRegion->BaseG + pSrcRegion->BaseT;
				if(SumACGT > 0)
					{
					pRegion->TransProbMatrix[Idx][0] = (double)pSrcRegion->BaseA/SumACGT;
					pRegion->TransProbMatrix[Idx][1] = (double)pSrcRegion->BaseC/SumACGT;
					pRegion->TransProbMatrix[Idx][2] = (double)pSrcRegion->BaseG/SumACGT;
					pRegion->TransProbMatrix[Idx][3] = (double)pSrcRegion->BaseT/SumACGT;
					// fudge to ensure probabilities sum to 1.0
					// adjust the maximal probability so all do sum to 1.0
					double SumProbs = pRegion->TransProbMatrix[Idx][0] + pRegion->TransProbMatrix[Idx][1] + pRegion->TransProbMatrix[Idx][2] + pRegion->TransProbMatrix[Idx][3];
					if(SumProbs != 1.0)
						{
						if(pRegion->TransProbMatrix[Idx][0] >= pRegion->TransProbMatrix[Idx][1] && 
							pRegion->TransProbMatrix[Idx][0] >= pRegion->TransProbMatrix[Idx][2] && 
							pRegion->TransProbMatrix[Idx][0] >= pRegion->TransProbMatrix[Idx][3])
							{
							pRegion->TransProbMatrix[Idx][0] = 1.0 - pRegion->TransProbMatrix[Idx][1] - pRegion->TransProbMatrix[Idx][2] - pRegion->TransProbMatrix[Idx][3];
							}
						else
							if(pRegion->TransProbMatrix[Idx][1] >= pRegion->TransProbMatrix[Idx][2] && 
								pRegion->TransProbMatrix[Idx][1] >= pRegion->TransProbMatrix[Idx][3])
								{
								pRegion->TransProbMatrix[Idx][1] = 1.0 - pRegion->TransProbMatrix[Idx][0] - pRegion->TransProbMatrix[Idx][2] - pRegion->TransProbMatrix[Idx][3];
								}
							else
								if(pRegion->TransProbMatrix[Idx][2] >= pRegion->TransProbMatrix[Idx][3])
									pRegion->TransProbMatrix[Idx][2] = 1.0 - pRegion->TransProbMatrix[Idx][0] - pRegion->TransProbMatrix[Idx][1] - pRegion->TransProbMatrix[Idx][3];
								else
									pRegion->TransProbMatrix[Idx][3] = 1.0 - pRegion->TransProbMatrix[Idx][0] - pRegion->TransProbMatrix[Idx][1] - pRegion->TransProbMatrix[Idx][2];

						}
					}
				else
					{
					pRegion->TransProbMatrix[Idx][0] = 0.25;
					pRegion->TransProbMatrix[Idx][1] = 0.25;
					pRegion->TransProbMatrix[Idx][2] = 0.25;
					pRegion->TransProbMatrix[Idx][3] = 0.25;
					}
				}
			}
		}
	return(true);
	}

double ChiSqr;
int Instances;

int RegionCentDist[cMaxChroms][cMaxRegions][4][4];
int OligoCentDist[cMaxChroms][4][4];
int Cells[cMaxChiSqrRows*cMaxChiSqrCols];

int *pCell;
int *pCnt;
double Fishers;
int BkgndCntA,BkgndCntB,SmplCntA,SmplCntB;
int PartSum;

CStats Stats;

memset(RegionCentDist,0,sizeof(RegionCentDist));
memset(OligoCentDist,0,sizeof(OligoCentDist));

// pass one, calculate all centroid counts
pOligo = pOligos;
for(Idx = 0; Idx < NumOligos; Idx++, pOligo++)
	{
	pRegion = pOligo->Region;
	// first calc total number of instances over all regions for current oligo
	for(Instances = Steps = 0; Steps < cMaxRegions; Steps++,pRegion++)
		{
		Instances += pRegion->BaseA;
		Instances += pRegion->BaseC;
		Instances += pRegion->BaseG;
		Instances += pRegion->BaseT;
		}
	if(!Instances)
			continue;

	switch(NMer) {
		case 3:
			CentIdx = (pOligo->SeqID & 0x0c) >> 2;
			break;

		case 5:
			CentIdx = (pOligo->SeqID & 0x030) >> 4;
			break;

		case 7:
			CentIdx = (pOligo->SeqID & 0x0c0) >> 6;
			break;

		case 9:
			CentIdx = (pOligo->SeqID & 0x0300) >> 8;
			break;
		}


		// for each region calculate centroid stats - only if sums to at least 25 counts in region...
	pRegion = pOligo->Region;
	for(Steps = 0; Steps < cMaxRegions; Steps++,pRegion++)
		{
		SumACGT = pRegion->BaseA+pRegion->BaseC+pRegion->BaseG+pRegion->BaseT;
		if(SumACGT < 25)
			continue;
		pCnt = &RegionCentDist[pOligo->ChromID-1][Steps][CentIdx][0];
		*pCnt++ += pRegion->BaseA;
		*pCnt++  += pRegion->BaseC;
		*pCnt++  += pRegion->BaseG;
		*pCnt  += pRegion->BaseT;
		pCnt = &OligoCentDist[pOligo->ChromID-1][CentIdx][0];
		*pCnt++  += pRegion->BaseA;
		*pCnt++  += pRegion->BaseC;
		*pCnt++  += pRegion->BaseG;
		*pCnt  += pRegion->BaseT;
		}
	}

// pass two, now calculate the ChiSqr stats
pOligo = pOligos;
for(Idx = 0; Idx < NumOligos; Idx++, pOligo++)
	{
	pRegion = pOligo->Region;
	for(Steps=Instances = 0; Steps < cMaxRegions; Steps++,pRegion++)
		{
		Instances += pRegion->BaseA;
		Instances += pRegion->BaseC;
		Instances += pRegion->BaseG;
		Instances += pRegion->BaseT;
		}
	if(!Instances)
			continue;

	switch(NMer) {
		case 3:
			CentIdx = (pOligo->SeqID & 0x0c) >> 2;
			break;

		case 5:
			CentIdx = (pOligo->SeqID & 0x030) >> 4;
			break;

		case 7:
			CentIdx = (pOligo->SeqID & 0x0c0) >> 6;
			break;

		case 9:
			CentIdx = (pOligo->SeqID & 0x0300) >> 8;
			break;
		}


	pRegion = pOligo->Region;
	for(Steps = 0; Steps < cMaxRegions; Steps++,pRegion++)
		{
		SumACGT = pRegion->BaseA+pRegion->BaseC+pRegion->BaseG+pRegion->BaseT;
		if(!SumACGT) {
			pRegion->ProbFixed = -0.1;
			pRegion->ProbTransistion = -0.1;
			pRegion->ProbTransversion = -0.1;
			continue;
			}
		pRegion->ProbTransistion = -0.1;
		pRegion->ProbTransversion = -0.1;

		switch(CentIdx) {
			case 0:					// centroid == 'a'
				pRegion->ProbFixed = (pRegion->BaseA * 1.0) / (double)SumACGT;
				PartSum = pRegion->BaseC+pRegion->BaseG+pRegion->BaseT;
				if(!PartSum)
					break;
				pRegion->ProbTransistion = (pRegion->BaseG * 1.0) / (double)PartSum;
				pRegion->ProbTransversion = ((pRegion->BaseC + pRegion->BaseT) * 1.0) / (double)PartSum;
				break;

			case 1:					// centroid == 'c'
				pRegion->ProbFixed = (pRegion->BaseC * 1.0) / (double)SumACGT;
				PartSum = pRegion->BaseA+pRegion->BaseG+pRegion->BaseT;
				if(!PartSum)
					break;
				pRegion->ProbTransistion = (pRegion->BaseT * 1.0) / (double)PartSum;
				pRegion->ProbTransversion = ((pRegion->BaseA + pRegion->BaseG) * 1.0) / (double)PartSum;
				break;
			case 2:					// centroid == 'g'
				pRegion->ProbFixed = (pRegion->BaseG * 1.0) / (double)SumACGT;
				PartSum = pRegion->BaseA+pRegion->BaseC+pRegion->BaseT;
				if(!PartSum)
					break;
				pRegion->ProbTransistion = (pRegion->BaseA * 1.0) / (double)PartSum;
				pRegion->ProbTransversion = ((pRegion->BaseC + pRegion->BaseT) * 1.0) / (double)PartSum;
				break;

			case 3:					// centroid == 't'
			default:
				pRegion->ProbFixed = (pRegion->BaseT * 1.0) / (double)SumACGT;
				PartSum = pRegion->BaseA+pRegion->BaseC+pRegion->BaseG;
				if(!PartSum)
					break;
				pRegion->ProbTransistion = (pRegion->BaseC * 1.0) / (double)PartSum;
				pRegion->ProbTransversion = ((pRegion->BaseA + pRegion->BaseG) * 1.0) / (double)PartSum;
				break;
			}
		}

		// for each region calculate chi-square and fisher stats - but only if sums to at least 25 counts in region...
	pRegion = pOligo->Region;
	for(Steps = 0; Steps < cMaxRegions; Steps++,pRegion++)
		{
		SumACGT = pRegion->BaseA+pRegion->BaseC+pRegion->BaseG+pRegion->BaseT;
		
		if(SumACGT < 25)
			{
			pRegion->ChiSqr = 0.0;
			pRegion->ChiSqrA = 0.0;
			continue;
			}

		pCnt = &OligoCentDist[pOligo->ChromID-1][CentIdx][0];
		pCell = Cells;
		*pCell++=*pCnt++ - pRegion->BaseA;
		*pCell++=*pCnt++ - pRegion->BaseC;
		*pCell++=*pCnt++ - pRegion->BaseG;
		*pCell++=*pCnt++ - pRegion->BaseT;
		*pCell++=pRegion->BaseA;
		*pCell++=pRegion->BaseC;
		*pCell++=pRegion->BaseG;
		*pCell++=pRegion->BaseT;
		ChiSqr = Stats.CalcChiSqr(2,4,Cells);
		if(ChiSqr > 1000000.0)		// ensure no subsequent overflows
			ChiSqr = 1000000.0;
		pRegion->ChiSqrA = ChiSqr;
		pCnt = &RegionCentDist[pOligo->ChromID-1][Steps][CentIdx][0];
		pCell = Cells;
		*pCell++=*pCnt++ - pRegion->BaseA;
		*pCell++=*pCnt++ - pRegion->BaseC;
		*pCell++=*pCnt++ - pRegion->BaseG;
		*pCell++=*pCnt++ - pRegion->BaseT;
		ChiSqr = Stats.CalcChiSqr(2,4,Cells);
		if(ChiSqr > 1000000.0)		// ensure no subsequent overflows
			ChiSqr = 1000000.0;
		pRegion->ChiSqr = ChiSqr;
		}

	// calc Fishers exact for P(Transition vs Transversion)
	pRegion = pOligo->Region;
	for(Steps = 0; Steps < cMaxRegions; Steps++,pRegion++)
		{
		SumACGT = pRegion->BaseA+pRegion->BaseC+pRegion->BaseG+pRegion->BaseT;


		if(SumACGT < 5)
			{
			pRegion->FTransA = -1.0;
			pRegion->FTrans = -1.0;
			pRegion->FMutA = -1.0;
			pRegion->FMut = -1.0;
			continue;
			}
		
		pCnt = &OligoCentDist[pOligo->ChromID-1][CentIdx][0];
		switch(CentIdx) {
			case 0:						// A
				BkgndCntA = pCnt[1] + pCnt[3];		// A->(C|T) 
				BkgndCntB = pCnt[2];				// A->G
				SmplCntA = pRegion->BaseC + pRegion->BaseT;
				SmplCntB = pRegion->BaseG;
				break;
			case 1:						// C
				BkgndCntA = pCnt[0] + pCnt[2];		// C->(A|G)
				BkgndCntB = pCnt[3];				// C->T
				SmplCntA = pRegion->BaseA + pRegion->BaseG;
				SmplCntB = pRegion->BaseT;
				break;
			case 2:						// G
				BkgndCntA = pCnt[1] + pCnt[3];		// G->(C|T)
				BkgndCntB = pCnt[0];				// G->A
				SmplCntA = pRegion->BaseC + pRegion->BaseT;
				SmplCntB = pRegion->BaseA;
				break;		
			case 3:						// T
				BkgndCntA = pCnt[0] + pCnt[2];		// T->(A|G)
				BkgndCntB = pCnt[1];				// T->C
				SmplCntA = pRegion->BaseA + pRegion->BaseG;
				SmplCntB = pRegion->BaseC;
				break;
			}
		BkgndCntA -= SmplCntA;
		BkgndCntB -= SmplCntB;
		Fishers = Stats.FishersExactTest(BkgndCntA,BkgndCntB,SmplCntA,SmplCntB);

		pRegion->FTransA = Fishers;
		pCnt = &RegionCentDist[pOligo->ChromID-1][Steps][CentIdx][0];
		switch(CentIdx) {
			case 0:						// A
				BkgndCntA = pCnt[1] + pCnt[3];		// A->(C|T) 
				BkgndCntB = pCnt[2];				// A->G
				SmplCntA = pRegion->BaseC + pRegion->BaseT;
				SmplCntB = pRegion->BaseG;
				break;
			case 1:						// C
				BkgndCntA = pCnt[0] + pCnt[2];		// C->(A|G)
				BkgndCntB = pCnt[3];				// C->T
				SmplCntA = pRegion->BaseA + pRegion->BaseG;
				SmplCntB = pRegion->BaseT;
				break;
			case 2:						// G
				BkgndCntA = pCnt[1] + pCnt[3];		// G->(C|T)
				BkgndCntB = pCnt[0];				// G->A
				SmplCntA = pRegion->BaseC + pRegion->BaseT;
				SmplCntB = pRegion->BaseA;
				break;		
			case 3:						// T
				BkgndCntA = pCnt[0] + pCnt[2];		// T->(A|G)
				BkgndCntB = pCnt[1];				// T->C
				SmplCntA = pRegion->BaseA + pRegion->BaseG;
				SmplCntB = pRegion->BaseC;
				break;
			}
		BkgndCntA -= SmplCntA;
		BkgndCntB -= SmplCntB;

		Fishers = Stats.FishersExactTest(BkgndCntA,BkgndCntB,SmplCntA,SmplCntB);
		pRegion->FTrans = Fishers;

		// now for fishers exact on fixation rates
		pCnt = &OligoCentDist[pOligo->ChromID-1][CentIdx][0];
		switch(CentIdx) {
			case 0:						// A
				BkgndCntA = pCnt[1] + pCnt[2] + pCnt[3];  // C+G+T 
				BkgndCntB = pCnt[0];				      // A
				SmplCntA = pRegion->BaseC + pRegion->BaseG + pRegion->BaseT;
				SmplCntB = pRegion->BaseA;
				break;
			case 1:						// C
				BkgndCntA = pCnt[0] + pCnt[2] + pCnt[3];  // A+G+T 
				BkgndCntB = pCnt[1];				      // C
				SmplCntA = pRegion->BaseA + pRegion->BaseG + pRegion->BaseT;
				SmplCntB = pRegion->BaseC;
				break;
			case 2:						// G
				BkgndCntA = pCnt[0] + pCnt[1] + pCnt[3];  // A+C+T 
				BkgndCntB = pCnt[2];				      // G
				SmplCntA = pRegion->BaseA + pRegion->BaseC + pRegion->BaseT;
				SmplCntB = pRegion->BaseG;
				break;		
			case 3:						// T
				BkgndCntA = pCnt[0] + pCnt[1] + pCnt[2];  // C+G+T 
				BkgndCntB = pCnt[3];				      // T
				SmplCntA = pRegion->BaseA + pRegion->BaseC + pRegion->BaseG;
				SmplCntB = pRegion->BaseT;
				break;
			}
		BkgndCntA -= SmplCntA;
		BkgndCntB -= SmplCntB;
		Fishers = Stats.FishersExactTest(BkgndCntA,BkgndCntB,SmplCntA,SmplCntB);
		pRegion->FMutA = Fishers;

		pCnt = &RegionCentDist[pOligo->ChromID-1][Steps][CentIdx][0];
		switch(CentIdx) {
			case 0:						// A
				BkgndCntA = pCnt[1] + pCnt[2] + pCnt[3];  // C+G+T 
				BkgndCntB = pCnt[0];				      // A
				SmplCntA = pRegion->BaseC + pRegion->BaseG + pRegion->BaseT;
				SmplCntB = pRegion->BaseA;
				break;
			case 1:						// C
				BkgndCntA = pCnt[0] + pCnt[2] + pCnt[3];  // A+G+T 
				BkgndCntB = pCnt[1];				      // C
				SmplCntA = pRegion->BaseA + pRegion->BaseG + pRegion->BaseT;
				SmplCntB = pRegion->BaseC;
				break;
			case 2:						// G
				BkgndCntA = pCnt[0] + pCnt[1] + pCnt[3];  // A+C+T 
				BkgndCntB = pCnt[2];				      // G
				SmplCntA = pRegion->BaseA + pRegion->BaseC + pRegion->BaseT;
				SmplCntB = pRegion->BaseG;
				break;		
			case 3:						// T
				BkgndCntA = pCnt[0] + pCnt[1] + pCnt[2];  // C+G+T 
				BkgndCntB = pCnt[3];				      // T
				SmplCntA = pRegion->BaseA + pRegion->BaseC + pRegion->BaseG;
				SmplCntB = pRegion->BaseT;
				break;
			}
		BkgndCntA -= SmplCntA;
		BkgndCntB -= SmplCntB;
		Fishers = Stats.FishersExactTest(BkgndCntA,BkgndCntB,SmplCntA,SmplCntB);
		pRegion->FMut = Fishers;
		}
	}
return(true);
}


const int cMaxPeriods = 100;

bool
OutputStats(int Mode,				// processing mode
			char *pszRsltsFile)		// stats results file
{
int hRsltsFile;
char szLineBuff[cLineBuffSize];
int LineOfs;
tsRegionStats *pRegion;
tsCentOligo *pOligo;
int OligoIdx;
int RegionIdx;
int MatrixRow;


CStats Stats;

#ifdef _WIN32
if((hRsltsFile = open(pszRsltsFile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
if((hRsltsFile = open(pszRsltsFile,O_RDWR | O_CREAT |O_TRUNC, S_IREAD | S_IWRITE ))==-1)
#endif
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create %s - %s",pszRsltsFile,strerror(errno));
	return(false);
	}
// write header
LineOfs = sprintf(szLineBuff,"\"Chrom\",SeqID,\"Centroid\",\"Centroid3\",\"Seq\"");

if(Mode < eFTTransMatrix)
	{
	LineOfs+= sprintf(&szLineBuff[LineOfs],",\"IGPFixed\",\"IGPTrans\",\"IGPTransv\"");
	LineOfs+= sprintf(&szLineBuff[LineOfs],",\"US5PFixed\",\"US5PTrans\",\"US5PTransv\"");
	LineOfs+= sprintf(&szLineBuff[LineOfs],",\"UTR5PFixed\",\"UTR5PTrans\",\"UTR5PTransv\"");
	LineOfs+= sprintf(&szLineBuff[LineOfs],",\"CDSPFixed\",\"CDSPTrans\",\"CDSPTransv\"");
	LineOfs+= sprintf(&szLineBuff[LineOfs],",\"IntronPFixed\",\"IntronPTrans\",\"IntronPTransv\"");
	LineOfs+= sprintf(&szLineBuff[LineOfs],",\"UTR3PFixed\",\"UTR3PTrans\",\"UTR3PTransv\"");
	LineOfs+= sprintf(&szLineBuff[LineOfs],",\"DS3PFixed\",\"DS3PTrans\",\"DS3PTransv\"");

	LineOfs+= sprintf(&szLineBuff[LineOfs],",\"IGChiSqrA\",\"IGChiSqrAP\",\"IGChiSqr\",\"IGChiSqrP\",\"IGFMutA\",\"IGFMut\",\"IGFTransA\",\"IGFTrans\"");
	LineOfs+= sprintf(&szLineBuff[LineOfs],",\"US5ChiSqrA\",\"US5ChiSqrAP\",\"US5ChiSqr\",\"US5ChiSqrP\",\"US5FMutA\",\"US5FMut\",\"US5FTransA\",\"US5FTrans\"");
	LineOfs+= sprintf(&szLineBuff[LineOfs],",\"UTR5ChiSqrA\",\"UTR5ChiSqrAP\",\"UTR5ChiSqr\",\"UTRChiSqrP\",\"UTR5FMutA\",\"UTR5FMut\",\"UTR5FTransA\",\"UTR5FTrans\"");
	LineOfs+= sprintf(&szLineBuff[LineOfs],",\"CDSChiSqrA\",\"CDSChiSqrAP\",\"CDSChiSqr\",\"CDSChiSqrP\",\"CDSFMutA\",\"CDSFMut\",\"CDSFTransA\",\"CDSFTrans\"");
	LineOfs+= sprintf(&szLineBuff[LineOfs],",\"INTChiSqrA\",\"INTChiSqrAP\",\"INTChiSqr\",\"INTChiSqrP\",\"INTFMutA\",\"INTFMut\",\"INTFTransA\",\"INTFTrans\"");
	LineOfs+= sprintf(&szLineBuff[LineOfs],",\"UTR3ChiSqrA\",\"UTR3ChiSqrAP\",\"UTR3ChiSqr\",\"UTR3ChiSqrP\",\"UTR3FMutA\",\"UTR3FMut\",\"UTR3FTransA\",\"UTR3FTrans\"");
	LineOfs+= sprintf(&szLineBuff[LineOfs],",\"DS3ChiSqrA\",\"DS3ChiSqrAP\",\"DS3ChiSqr\",\"DS3ChiSqrP\",\"DS3FMutA\",\"DS3FMut\",\"DS3FTransA\",\"DS3FTrans\"");
	}
else
	{
	LineOfs+= sprintf(&szLineBuff[LineOfs],",\"IGPaa\",\"IGPac\",\"IGPag\",\"IGPat\"");
	LineOfs+= sprintf(&szLineBuff[LineOfs],",\"IGPca\",\"IGPcc\",\"IGPcg\",\"IGPct\"");
	LineOfs+= sprintf(&szLineBuff[LineOfs],",\"IGPga\",\"IGPgc\",\"IGPgg\",\"IGPgt\"");
	LineOfs+= sprintf(&szLineBuff[LineOfs],",\"IGPta\",\"IGPtc\",\"IGPtg\",\"IGPtt\"");
	LineOfs+= sprintf(&szLineBuff[LineOfs],",\"US5Paa\",\"US5Pac\",\"US5Pag\",\"US5Pat\"");
	LineOfs+= sprintf(&szLineBuff[LineOfs],",\"US5Pca\",\"US5Pcc\",\"US5Pcg\",\"US5Pct\"");
	LineOfs+= sprintf(&szLineBuff[LineOfs],",\"US5Pga\",\"US5Pgc\",\"US5Pgg\",\"US5Pgt\"");
	LineOfs+= sprintf(&szLineBuff[LineOfs],",\"US5Pta\",\"US5Ptc\",\"US5Ptg\",\"US5Ptt\"");
	LineOfs+= sprintf(&szLineBuff[LineOfs],",\"UTR5Paa\",\"UTR5Pac\",\"UTR5Pag\",\"UTR5Pat\"");
	LineOfs+= sprintf(&szLineBuff[LineOfs],",\"UTR5Pca\",\"UTR5Pcc\",\"UTR5Pcg\",\"UTR5Pct\"");
	LineOfs+= sprintf(&szLineBuff[LineOfs],",\"UTR5Pga\",\"UTR5Pgc\",\"UTR5Pgg\",\"UTR5Pgt\"");
	LineOfs+= sprintf(&szLineBuff[LineOfs],",\"UTR5Pta\",\"UTR5Ptc\",\"UTR5Ptg\",\"UTR5Ptt\"");
	LineOfs+= sprintf(&szLineBuff[LineOfs],",\"CDSPaa\",\"CDSPac\",\"CDSPag\",\"CDSPat\"");
	LineOfs+= sprintf(&szLineBuff[LineOfs],",\"CDSPca\",\"CDSPcc\",\"CDSPcg\",\"CDSPct\"");
	LineOfs+= sprintf(&szLineBuff[LineOfs],",\"CDSPga\",\"CDSPgc\",\"CDSPgg\",\"CDSPgt\"");
	LineOfs+= sprintf(&szLineBuff[LineOfs],",\"CDSPta\",\"CDSPtc\",\"CDSPtg\",\"CDSPtt\"");
	LineOfs+= sprintf(&szLineBuff[LineOfs],",\"INTPaa\",\"INTPac\",\"INTPag\",\"INTPat\"");
	LineOfs+= sprintf(&szLineBuff[LineOfs],",\"INTPca\",\"INTPcc\",\"INTPcg\",\"INTPct\"");
	LineOfs+= sprintf(&szLineBuff[LineOfs],",\"INTPga\",\"INTPgc\",\"INTPgg\",\"INTPgt\"");
	LineOfs+= sprintf(&szLineBuff[LineOfs],",\"INTPta\",\"INTPtc\",\"INTPtg\",\"INTPtt\"");
	LineOfs+= sprintf(&szLineBuff[LineOfs],",\"UTR3Paa\",\"UTR3Pac\",\"UTR3Pag\",\"UTR3Pat\"");
	LineOfs+= sprintf(&szLineBuff[LineOfs],",\"UTR3Pca\",\"UTR3Pcc\",\"UTR3Pcg\",\"UTR3Pct\"");
	LineOfs+= sprintf(&szLineBuff[LineOfs],",\"UTR3Pga\",\"UTR3Pgc\",\"UTR3Pgg\",\"UTR3Pgt\"");
	LineOfs+= sprintf(&szLineBuff[LineOfs],",\"UTR3Pta\",\"UTR3Ptc\",\"UTR3Ptg\",\"UTR3Ptt\"");
	LineOfs+= sprintf(&szLineBuff[LineOfs],",\"DS3Paa\",\"DS3Pac\",\"DS3Pag\",\"DS3Pat\"");
	LineOfs+= sprintf(&szLineBuff[LineOfs],",\"DS3Pca\",\"DS3Pcc\",\"DS3Pcg\",\"DS3Pct\"");
	LineOfs+= sprintf(&szLineBuff[LineOfs],",\"DS3Pga\",\"DS3Pgc\",\"DS3Pgg\",\"DS3Pgt\"");
	LineOfs+= sprintf(&szLineBuff[LineOfs],",\"DS3Pta\",\"DS3Ptc\",\"DS3Ptg\",\"DS3Ptt\"");
	}

LineOfs += sprintf(&szLineBuff[LineOfs],"\n");
CUtility::SafeWrite(hRsltsFile,szLineBuff,LineOfs);
LineOfs = 0;

// write stats
pOligo = m_pAllocdOligos;
for(OligoIdx = 0; OligoIdx < m_CurNumOligos; OligoIdx++,pOligo++)
	{
	LineOfs = sprintf(szLineBuff,"\"%s\",%d,\"%c\",\"%s\",\"%s\"",pOligo->szChrom,pOligo->SeqID,pOligo->Centroid,pOligo->szCentroid3,pOligo->szSeq);
	pRegion = pOligo->Region;
	if(Mode < eFTTransMatrix)
		{
		for(RegionIdx = 0; RegionIdx < cMaxRegions; RegionIdx++,pRegion++)
			{
			LineOfs += sprintf(&szLineBuff[LineOfs],",%f,%f,%f",
				pRegion->ProbFixed,pRegion->ProbTransistion,pRegion->ProbTransversion);
			}
		}

	pRegion = pOligo->Region;
	for(RegionIdx = 0; RegionIdx < cMaxRegions; RegionIdx++,pRegion++)
		{
		if(Mode != eFTTransMatrix)
			LineOfs += sprintf(&szLineBuff[LineOfs],",%f,%f,%f,%f,%f,%f,%f,%f",
				pRegion->ChiSqrA,Stats.ChiSqr2PVal(3,pRegion->ChiSqrA),
				pRegion->ChiSqr,Stats.ChiSqr2PVal(3,pRegion->ChiSqr),
				pRegion->FMutA,
				pRegion->FMut,
    			pRegion->FTransA,
				pRegion->FTrans);
		else
			{
			for(MatrixRow=0;MatrixRow<4;MatrixRow++)
					LineOfs += sprintf(&szLineBuff[LineOfs],",%f,%f,%f,%f",
						pRegion->TransProbMatrix[MatrixRow][0],
						pRegion->TransProbMatrix[MatrixRow][1],
						pRegion->TransProbMatrix[MatrixRow][2],
						pRegion->TransProbMatrix[MatrixRow][3]);
			}
		}
	LineOfs += sprintf(&szLineBuff[LineOfs],"\n");
	CUtility::SafeWrite(hRsltsFile,szLineBuff,LineOfs);
	LineOfs = 0;
	}

close(hRsltsFile);
return(true);
}



bool
StationaryProbs(int NumPeriods,		// number of time periods
				teSeqBases InitialBase,	// initially observed base
				double *pProbA,		// where to return probabilities in each time period for A
				double *pProbC,		// where to return probabilities in each time period for C
				double *pProbG,		// where to return probabilities in each time period for G
				double *pProbT,		// where to return probabilities in each time period for T
				double *pMatrix)	// 4*4 transitional probabilites
{
int Period;
double SumProbs;
if(NumPeriods < 1 ||
   InitialBase < eBaseA || InitialBase > eBaseT ||
   pProbA == NULL || pProbC == NULL || pProbG == NULL || pProbT == NULL ||
   pMatrix == NULL)
	return(false);

*pProbA++ = InitialBase == eBaseA ? 1.0 : 0.0;
*pProbC++ = InitialBase == eBaseC ? 1.0 : 0.0;
*pProbG++ = InitialBase == eBaseG ? 1.0 : 0.0;
*pProbT++ = InitialBase == eBaseT ? 1.0 : 0.0;

for(Period = 1; Period < NumPeriods; Period++,pProbA++,pProbC++,pProbG++,pProbT++)
	{
	*pProbA = pProbA[-1] * pMatrix[0] +			// a->a
			  pProbC[-1] * pMatrix[4] +			// c->a
			  pProbG[-1] * pMatrix[8] +			// g->a
			  pProbT[-1] * pMatrix[12];			// t->a
	
	*pProbC = pProbA[-1] * pMatrix[1] +			// a->c
			  pProbC[-1] * pMatrix[5] +			// c->c 
			  pProbG[-1] * pMatrix[9] +			// g->c
			  pProbT[-1] * pMatrix[13];			// t->c

	*pProbG = pProbA[-1] * pMatrix[2] +			// a->g
			  pProbC[-1] * pMatrix[6] +			// c->g
			  pProbG[-1] * pMatrix[10] +		// g->g	
			  pProbT[-1] * pMatrix[14];			// t->g

	*pProbT = pProbA[-1] * pMatrix[3] +			// a->t
			  pProbC[-1] * pMatrix[7] +			// c->t
			  pProbG[-1] * pMatrix[11] +		// g->t
			  pProbT[-1] * pMatrix[15];			// t->t

	// fudge to ensure probabilities sum to 1.0
	// adjust the maximal probability so all do sum to 1.0
	SumProbs = *pProbA + *pProbC + *pProbG + *pProbT;
	if(SumProbs != 1.0)
		{
		if(*pProbA >= *pProbC && *pProbA >= *pProbG && *pProbA >= *pProbT)
			{
			*pProbA = 1.0 - *pProbC - *pProbG - *pProbT;
			}
		else
			if(*pProbC >= *pProbG && *pProbC >= *pProbT)
				{
				*pProbC = 1.0 - *pProbA - *pProbG - *pProbT;
				}
			else
				if(*pProbG >= *pProbT)
					*pProbG = 1.0 - *pProbA - *pProbC - *pProbT;
				else
					*pProbT = 1.0 - *pProbA - *pProbC - *pProbG;
		}
	}

return(0);
}
