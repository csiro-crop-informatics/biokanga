// loci2dist.cpp : Defines the entry point for the console application.
// generates loci length distributions

#include "stdafx.h"
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "../libbiokanga/commhdrs.h"

const char *cpszProgVer = "0.0.1";		// increment with each release

const int cMaxLengthRange = 10000;		// maximal element length
const int cDfltLengthRange = 200;		// default element length

// input loci file processing modes
typedef enum TAG_ePMode {		
	ePMdefault,					// CSV loci
	ePMBed,						// UCSC BED format
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
			char *pszInLociFile,		// input CSV or BED loci file
			char *pszInBEDFile,			// input BED feature file
			char *pszRsltsFile,			// output loci file
			int RegRegionLen,			// regulatory region length
			int MinLength,				// minimum element length
			int MaxLength);				// maximum element length

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
int iRegLen;				// up/down stream regulatory region length
int iMinLength;
int iMaxLength;

int Strand;					// which element strand to filter (retain) on

char szInLociFile[_MAX_PATH];	// input element loci from this file
char szInBEDFile[_MAX_PATH];	// input bed file containing gene features
char szRsltsFile[_MAX_PATH];	// output result distributions to this file

// command line args
struct arg_lit  *help    = arg_lit0("hH","help",                "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"diagnostics log file");

struct arg_int *pmode = arg_int0("m","mode","<int>",		    "input reads loci file: 0 - CSV, 1 - BED (default = 0)");

struct arg_int  *strandproc = arg_int0("s","strandproc","<int>","strand processing: 0 - independent, 1 - Watson, 2 - Crick (default is independent)");

struct arg_file *InLociFile = arg_file1("i","incsv","<file>",	"input element CSV or BED file");
struct arg_file *InBEDFile = arg_file0("I","inbed","<file>",	"input gene or feature biobed BED file");
struct arg_file *RsltsFile = arg_file1("o","output","<file>",	 "length distributions output file");
struct arg_int  *RegLen = arg_int0("r","updnstream","<int>",	"length of 5'up or 3'down  stream regulatory region length (default = 2000) 0..1000000");
struct arg_int  *MinLength = arg_int0("l","minlength","<int>",   "minimum element length (default 1)");
struct arg_int  *MaxLength = arg_int0("L","maxlength","<int>",   "maximum element length (default 500)");

struct arg_end *end = arg_end(20);

void *argtable[] = {help,version,FileLogLevel,LogFile,
					pmode,strandproc,InLociFile,InBEDFile,RsltsFile,RegLen,
					MinLength,MaxLength,
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
		printf("\n      To invoke this parameter file then precede its name with '@'");
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
	
	Strand = strandproc->count ? strandproc->ival[0] : eStrandDflt;
	if(Strand < eStrandDflt || Strand >= eStrandPlaceholder)
		{
		printf("\nError: Strand processing mode '-s%d' must be in range %d..%d",Strand,eStrandDflt,eStrandCrick);
		exit(1);
		}
	switch(Strand) {
		case 1: Strand = (int)'+'; break;
		case 2: Strand = (int)'-'; break;
		case 0: Strand = (int)'*'; break;
		}	
	iMinLength = MinLength->count ? MinLength->ival[0] : 1;
	if(iMinLength < 1 || iMinLength > cMaxLengthRange)
		{
		printf("Error: Minimum element length '-m%d' is not in range 1..%d",iMinLength,cMaxLengthRange);
		exit(1);
		}

	iMaxLength = MaxLength->count ? MaxLength->ival[0] : cDfltLengthRange;
	if(iMaxLength < iMinLength || iMaxLength > cMaxLengthRange)
		{
		printf("Error: Maximum element length '-M%d' is not in range %d..%d",iMaxLength,iMinLength,cMaxLengthRange);
		exit(1);
		}

	if(InBEDFile->count)
		{
		printf("Sorry, regional processing from feature/genes not currently implemented");
		exit(1);

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
		szInBEDFile[0] = '\0';
		iRegLen = 0;
		}

	strncpy(szInLociFile,InLociFile->filename[0],_MAX_PATH);
	szInLociFile[_MAX_PATH-1] = '\0';

	strncpy(szRsltsFile,RsltsFile->filename[0],_MAX_PATH);
	szRsltsFile[_MAX_PATH-1] = '\0';

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
			pszDescr = "CSV loci";
			break;
		case ePMBed:			
			pszDescr = "UCSC BED loci";
			break;
		}
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"processing mode for input loci file is : '%s'",pszDescr);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Strand processing : '%c'",(char)Strand);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Input CSV or BED element loci file: '%s'",szInLociFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Length distributions to file: '%s'",szRsltsFile);

	if(szInBEDFile[0] != '\0')
		{
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Input biobed BED gene feature file: '%s'",szInBEDFile);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Feature mappings with regulatory region length: %d",iRegLen);
		}
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Minimum element length: %d",iMinLength);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Maximum element length: %d",iMaxLength);

#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	gStopWatch.Start();
	Rslt = Process(PMode,Strand,szInLociFile,szInBEDFile,szRsltsFile,iRegLen,iMinLength,iMaxLength);

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

typedef struct TAG_sFeatCntDist {
	int FeatID;							// identifies this feature
	int GeneLen;						// gene length
	int TranscribedLen;					// transcribed length
	char szName[80];					// gene name
	int RegionCnts[8];					// region counts
	double RPKM;						// RPKM for this feature
	int NumExonReads;					// number of reads in UTR's plus CDS
	bool bIsIsoform;					// set true if feature is an isoform, assumed if name at least 9 long and suffixed by ".[0-99]"
	bool bMaxRPKM;						// set true if feature is an isoform with the maximal RPKM
	bool bMaxExonReads;				// set true if feature is an isoform with the maximal reads in UTR's+CDS
} tsFeatCntDist;


typedef struct TAG_sChromDist {
	struct TAG_sChromDist *pNext;	// tsChromDist's are simply linked together in order of creation
	int ChromID;					// uniquely identifies this chromosome
	char szChrom[100];				// chromosome name
	int MinLen;						// minimum length encountered
	int MaxLen;						// maximum length encountered
	int TotCnts;					// total of all counts attributed to this chromosome
	int Dist[cMaxLengthRange+1];	// length distribution over this chromosome
	int m_FeatDist[9][cMaxLengthRange+1];	// to hold length distributions for each functional region
} tsChromDist;

int m_NumChromDists;				// current number of chromosome distribution histograms
tsChromDist *m_pChromDists;			// pts to linked list of chromosome distribution histograms
tsChromDist *m_pGenomeDists;		// pts to genome wide distribution histograms

CHyperEls *m_pHypers;
tsFeatCntDist *m_pFeatCntDists;		// to hold feature count distributions
CBEDfile *m_pBiobed;				// BED file containing features/genes

etPMode m_PMode;					// processing mode
char m_Strand;						// strand to filter (retain) on
int m_RegRegionLen;					// regulatory region length
UINT32 m_NumEls;					// total number of elements to be processed

int m_hRsltsFile;


CCSVFile *m_pCSVFile;					// used if processing input CSV file
CBEDfile *m_pBEDFile;					// used if processing input BED file


void
Reset(void)
{
if(m_hRsltsFile != -1)
	{
	close(m_hRsltsFile);
	m_hRsltsFile = -1;
	}
if(m_pBiobed != NULL)
	{
	delete m_pBiobed;
	m_pBiobed = NULL;
	}
if(m_pBEDFile != NULL)
	{
	delete m_pBEDFile;
	m_pBEDFile = NULL;
	}
if(m_pCSVFile != NULL)
	{
	delete m_pCSVFile;
	m_pCSVFile = NULL;
	}
if(m_pChromDists)
	{
	tsChromDist *pTmp;
	while(m_pChromDists != NULL)
		{
		pTmp = m_pChromDists->pNext;
		delete m_pChromDists;
		m_pChromDists = pTmp;
		}
	m_pChromDists = NULL;
	}

if(m_pGenomeDists != NULL)
	{
	delete m_pGenomeDists;
	m_pGenomeDists = NULL;
	}

if(m_pFeatCntDists != NULL)
	{
	delete m_pFeatCntDists;
	m_pFeatCntDists = NULL;
	}
m_PMode = ePMdefault; 
m_Strand = eStrandDflt;
m_RegRegionLen = cDfltRegLen;
m_NumEls = 0;
m_NumChromDists = 0;
}

void
Init(void)
{
m_pBiobed = NULL;
m_pCSVFile = NULL;
m_pBEDFile = NULL;
m_pChromDists = NULL;
m_pGenomeDists = NULL;
m_pFeatCntDists = NULL;
m_hRsltsFile = -1;
Reset();
}

#ifdef USETHISCODE
// MapLoci2Features
// Maps element loci to features and writes mapping for each element loci to pszRsltsFile
// Assumes that the bed file containing features (m_pBiobed) has been opened and that all element loci have been parsed into m_pHypers
int
MapLoci2Features(char *pszRsltsFile)
{
int Rslt;
char szLineBuff[0x07fff];
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
	if((m_hRsltFile = open(pszRsltsFile, O_RDWR | O_CREAT |O_TRUNC, S_IREAD | S_IWRITE))==-1)
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

printf("\nAssociating element %8.8d",1);
for(ElID = 1; ElID <= m_NumEls; ElID++)
	{
	if(!(ElID % 100000))
		printf("\b\b\b\b\b\b\b\b%8.8d",ElID);
	pEl = m_pHypers->GetElement(ElID);
	if(pszChrom == NULL || PrevChromID != pEl->ChromID)
		{
		pszChrom = m_pHypers->GetChrom(pEl->ChromID);
		// some old datasets may be referencing ChrM as mitochondria, or ChrC as chloroplast
		// so need to check for these
		if(!stricmp(pszChrom,"chloroplast"))
			pszChrom = (char *)"ChrC";
		else
			if(!stricmp(pszChrom,"mitochondria"))
				pszChrom = (char *)"ChrM";
		PrevChromID = pEl->ChromID;
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
	SrcID = pEl->SrcID;
	Features = pEl->Features;
	switch(m_Strand) {
		case eStrandDflt:
			break;
		case eStrandWatson:
			m_pBiobed->SetStrand(Strand);
			break;
		case eStrandCrick:
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

		if(Rslt > 0)
			{
			if(MaxChromID < ChromID)
				MaxChromID = ChromID;

			CoreStartLoci = StartLoci;
			CoreEndLoci = EndLoci;

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
					Features = m_pBiobed->GetFeatureBits(ChromID,			// feature is on which chromsome
									 CoreStartLoci,							// feature must end on or after Start
									 CoreEndLoci,							// and start on or before End
									cRegionFeatBits,
									m_RegRegionLen);
					Features |= m_pBiobed->GetSpliceSiteBits(ChromID,CoreStartLoci,CoreEndLoci,cMinSpliceOverlap);
					AccumFeatures |= Features;
					pCurFeatCntDist = &m_pFeatCntDists[FeatID-1];
					for(FeatMsk = 0x01, FeatIdx = 0; FeatIdx < 8; FeatIdx++, FeatMsk <<= 1)
						if(Features & FeatMsk)
							pCurFeatCntDist->RegionCnts[FeatIdx] += 1;
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
				
				AccumFeatures = m_pBiobed->GetFeatureBits(ChromID,	  // feature is on which chromsome
									 CoreStartLoci,       // feature must end on or after Start
									 CoreEndLoci,		  // and start on or before End
									cRegionFeatBits,
									m_RegRegionLen);
				
			
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
					SrcID,pszElType,pszRefSpecies,pszChrom,StartLoci,EndLoci,Len,Strand,AccumFeatures);
	
	if(!AccumFeatures)
		m_FeatureDist[ChromID-1][8] += 1;
	else
		for(FeatMsk = 0x01, FeatIdx = 0; FeatIdx < 8; FeatIdx++, FeatMsk <<= 1)
			if(AccumFeatures & FeatMsk)
				m_FeatureDist[ChromID-1][FeatIdx] += 1;

	if(m_hRsltFile != -1 && ((BuffIdx + 200) > sizeof(szLineBuff)))
		{
		CUtility::SafeWrite(m_hRsltFile,szLineBuff,BuffIdx);
		BuffIdx = 0;
		}
	}
printf("\b\b\b\b\b\b\b\b%8.8d",ElID-1);

if(m_hRsltFile != -1)
	{
	if(BuffIdx > 0)
		CUtility::SafeWrite(m_hRsltFile,szLineBuff,BuffIdx);
	close(m_hRsltFile);
	m_hRsltFile = -1;
	}

return(eBSFSuccess);
}
#else

int
MapLoci2Features(char *pszRsltsFile)
{
return(0);
}
#endif

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
		pSfxA = &pszFeatA[NameLen-1];
		}
	}
else
	ChrA = '\0';

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
		pSfxB = &pszFeatB[NameLen-1];
		}
	}
else
	ChrB = '\0';
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

#ifdef USETHISCODE
int
MapFeatures2Loci(char *pszRsltsFile)
{
char szLineBuff[0x07fff];
int BuffIdx;
int TotNumFeatures;
int FeatID;
int FeatIdx;
double RPKM;

tsFeatCntDist *pCurFeatCntDist;		// to hold currently being processed feature count distribution
tsFeatCntDist *pMaxRPKM;		// best RPKM isoform instance
tsFeatCntDist *pMaxExonReads;		// best reads isoform instance
bool bSameMaxRPKMFeature;		// true if processing same feature, different isoforms, for maximal RPKMs
bool bSameMaxExonReadsFeature;  // true if processing same feature,  different isoforms, for maximal exonic reads

if(m_hRsltFile != -1)		// shouldn't be open but let's make sure...
	{
	close(m_hRsltFile);
	m_hRsltFile = -1;
	}
if(m_pBiobed != NULL && pszRsltsFile != NULL && pszRsltsFile[0] != '\0')
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
	if((m_hRsltFile = open(pszRsltsFile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
	if((m_hRsltFile = open(pszRsltsFile,O_RDWR | O_CREAT |  O_TRUNC, S_IREAD | S_IWRITE))==-1)
#endif
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create %s - %s",pszRsltsFile,strerror(errno));
		Reset();
		return(eBSFerrCreateFile);
		}
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Features output file created/truncated: '%s'",pszRsltsFile);
	}
else
	TotNumFeatures = 0;



if(m_hRsltFile != -1 && m_pFeatCntDists != NULL)
	{
	// sort by feature name 
	qsort(m_pFeatCntDists,TotNumFeatures,sizeof(tsFeatCntDist),CompareFeatName);
	
	// generate RPKMs and ExonReads over all features and flag if feature is an isoform
	pCurFeatCntDist = m_pFeatCntDists;
	pMaxRPKM = NULL;
	pMaxExonReads = NULL;
	for(FeatID=1;FeatID<=TotNumFeatures;FeatID++,pCurFeatCntDist++)
		{
		pCurFeatCntDist = &m_pFeatCntDists[FeatID-1];
		pCurFeatCntDist->NumExonReads = pCurFeatCntDist->RegionCnts[0] + pCurFeatCntDist->RegionCnts[1] + pCurFeatCntDist->RegionCnts[2];
		if(pCurFeatCntDist->NumExonReads)	// if at least 1 read mapping to feature...
			{
			RPKM = (double)(pCurFeatCntDist->NumExonReads * 1000);
			RPKM /= (double)pCurFeatCntDist->TranscribedLen;
			RPKM *= (double)1000000/m_NumEls;
			}
		else
			RPKM = 0.0f;
		pCurFeatCntDist->RPKM = RPKM;

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


	BuffIdx = sprintf(szLineBuff,"\"FeatID\",\"Feature\",\"GeneLen\",\"TransLen\",\"CDS\",\"5'UTR\",\"3'UTR\",\"Introns\",\"5'upstream\",\"3'downstream\",\"intron3'/5'exon\",\"exon3'/5'intron\",\"RPKM\",\"ExonReads\"");
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
		BuffIdx += sprintf(&szLineBuff[BuffIdx],",%1.2f,%d",pCurFeatCntDist->RPKM,pCurFeatCntDist->NumExonReads);
		if((BuffIdx + 200) > sizeof(szLineBuff))
			{
			CUtility::SafeWrite(m_hRsltFile,szLineBuff,BuffIdx);
			BuffIdx = 0;
			}
		}
	if(BuffIdx)
		CUtility::SafeWrite(m_hRsltFile,szLineBuff,BuffIdx);
	close(m_hRsltFile);
	m_hRsltFile = -1;
	}
if(m_pFeatCntDists != NULL)
	{
	delete m_pFeatCntDists;
	m_pFeatCntDists = NULL;
	}
return(eBSFSuccess);
}
#else
int
MapFeatures2Loci(char *pszRsltsFile)
{
return(0);
}
#endif




int
BuildReadCoverage(char *pszChrom,int StartLoci,int EndLoci,int ElLen)
{
tsChromDist *pChromDist;
tsChromDist *pPrevChromDist;

if(m_pGenomeDists == NULL)
	{
	m_pGenomeDists = new tsChromDist;
	memset(m_pGenomeDists,0,sizeof(tsChromDist));
	strcpy(m_pGenomeDists->szChrom,(const char *)"GenomeWide");
	m_pGenomeDists->MinLen = ElLen;
	m_pGenomeDists->MaxLen = ElLen;
	}

pChromDist = m_pChromDists;
pPrevChromDist = NULL;
while(pChromDist != NULL)
	{
	if(!stricmp(pszChrom,pChromDist->szChrom))
		break;
	pPrevChromDist = pChromDist;
	pChromDist = pChromDist->pNext;
	}
if(pChromDist == NULL)
	{
	pChromDist = new tsChromDist;
	memset(pChromDist,0,sizeof(tsChromDist));
	pChromDist->ChromID = ++m_NumChromDists;
	strcpy(pChromDist->szChrom,pszChrom);
	if(pPrevChromDist != NULL)
		pPrevChromDist->pNext = pChromDist;
	else
		m_pChromDists = pChromDist;
	pChromDist->MinLen = ElLen;
	pChromDist->MaxLen = ElLen;
	}
pChromDist->TotCnts += 1;
if(ElLen < pChromDist->MinLen)
	pChromDist->MinLen = ElLen;
if(ElLen > pChromDist->MaxLen)
	pChromDist->MaxLen = ElLen;
pChromDist->Dist[ElLen] += 1;

m_pGenomeDists->TotCnts += 1;
if(ElLen < m_pGenomeDists->MinLen)
	m_pGenomeDists->MinLen = ElLen;
if(ElLen > m_pGenomeDists->MaxLen)
	m_pGenomeDists->MaxLen = ElLen;
m_pGenomeDists->Dist[ElLen] += 1;

return(eBSFSuccess);
}

int
ProcessElsCSV(char Strand,					// process for this strand only
		   int MinLen,					// elements must be of at least this length
		   int MaxLen,					// elements must be no longer than this length
		   char *pszInputFile)			// load elements from this file
{
int Rslt;
int NumEls;
int NumFields;
char *pszChrom;
int StartLoci;
int EndLoci;
char *pszStrand;
int NumFiltStrand;
int	ElLen;
int NumAccepted;
int UnderLength;
int OverLength;

if((m_pCSVFile = new CCSVFile)==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CCSVFile");
	Reset();
	return(eBSFerrObj);
	}

if((Rslt=m_pCSVFile->Open(pszInputFile)) !=eBSFSuccess)
	{
	while(m_pCSVFile->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,m_pCSVFile->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open '%s' for processing",pszInputFile);
	Reset();
	return(eBSFerrOpnFile);
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing elements from '%s'",pszInputFile);
NumEls = 0;
NumFiltStrand = 0;
UnderLength = 0;
OverLength = 0;
NumAccepted = 0;
while((Rslt=m_pCSVFile->NextLine()) > 0)	// onto next line containing fields
	{
	NumFields = m_pCSVFile->GetCurFields();
	if(NumFields < 7)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expected at least 7 fields in '%s', GetCurFields() returned '%d'",pszInputFile,NumFields);
		return(eBSFerrFieldCnt);
		}
	NumEls += 1;
	if(Strand != '*')
		{
		if(NumFields < 8)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Strand filtering, expected at least 8 fields in '%s', GetCurFields() returned '%d'",pszInputFile,NumFields);
			return(eBSFerrFieldCnt);
			}
		m_pCSVFile->GetText(8,&pszStrand);
		if(*pszStrand != Strand)
			{
			NumFiltStrand += 1;
			continue;
			}
		}

	m_pCSVFile->GetText(4,&pszChrom);
	m_pCSVFile->GetInt(5,&StartLoci);
	m_pCSVFile->GetInt(6,&EndLoci);
	m_pCSVFile->GetInt(7,&ElLen);
	if(ElLen < MinLen)
		{
		UnderLength += 1;
		continue;
		}
	if(ElLen > MaxLen)
		{
		OverLength += 1;
		continue;
		}

	if((Rslt=BuildReadCoverage(pszChrom,StartLoci,EndLoci,ElLen))<eBSFSuccess)
		break;
	NumAccepted += 1;
	}
if(Rslt >= eBSFSuccess)
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"From %d elements there were %d accepted, %d filtered out because of strand, %d under length and %d over length",
					NumEls,NumAccepted,NumFiltStrand,UnderLength,OverLength);
return(Rslt >= eBSFSuccess ? NumAccepted : Rslt);
}

int
ProcessElsBED(char Strand,				// process for this strand only
		   int MinLen,					// elements must be of at least this length
		   int MaxLen,					// elements must be no longer than this length
		   char *pszInputFile)			// load elements from this file
{
int Rslt;
int CurFeatureID;
int StartLoci;
int EndLoci;
int Score;
char szChrom[128];
char szFeatName[128];
char FeatStrand;
int RefID;
int IntergenicStart;
int NumEls;
int NumFiltStrand;
int	ElLen;
int NumAccepted;
int UnderLength;
int OverLength;

if((m_pBEDFile = new CBEDfile)==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CBEDfile");
	Reset();
	return(eBSFerrObj);
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Starting to load reads %s",pszInputFile);
if((Rslt=m_pBEDFile->Open(pszInputFile,eBTAnyBed)) !=eBSFSuccess)
	{
	while(m_pBEDFile->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,m_pBEDFile->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open '%s' for processing",pszInputFile);
	Reset();
	return(eBSFerrOpnFile);
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Load completed, iterating elements...");

// now iterate over the elements 
CurFeatureID = 0;
IntergenicStart = 0;
RefID = 0;
Rslt = eBSFSuccess;
NumEls = 0;
NumFiltStrand = 0;
UnderLength = 0;
OverLength = 0;
NumAccepted = 0;

while(Rslt == eBSFSuccess && (CurFeatureID = m_pBEDFile->GetNextFeatureID(CurFeatureID)) > 0)
	{
	NumEls += 1;
	m_pBEDFile->GetFeature(CurFeatureID,	// feature instance identifier
				szFeatName,				// where to return feature name
				szChrom,				// where to return chromosome name
				&StartLoci,				// where to return feature start on chromosome (0..n) 
				&EndLoci,				// where to return feature end on chromosome
 				&Score,					// where to return score
 				&FeatStrand);				// where to return strand

	if(Strand != '*')
		{
		if(FeatStrand != Strand)
			{
			NumFiltStrand += 1;
			continue;
			}
		}

	ElLen = 1 + EndLoci - StartLoci;
	if(ElLen < MinLen)
		{
		UnderLength += 1;
		continue;
		}
	if(ElLen > MaxLen)
		{
		OverLength += 1;
		continue;
		}

	if((Rslt=BuildReadCoverage(szChrom,StartLoci,EndLoci,ElLen)) < eBSFSuccess)
		break;
	NumAccepted += 1;
	}
if(Rslt >= eBSFSuccess)
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"From %d elements there were %d accepted, %d filtered out because of strand, %d under length and %d over length",
					NumEls,NumAccepted,NumFiltStrand,UnderLength,OverLength);
return(Rslt >= eBSFSuccess ? NumAccepted : Rslt);
}

int
ReportChromDists(int MinLength,				// minimum element length
		   int MaxLength)				// maximum element length
{
tsChromDist *pDist;
int CntIdx;
char szBuffer[16000];
int BuffOfs;

// firstly the heading..
BuffOfs = sprintf(szBuffer,",\"%s\"",m_pGenomeDists->szChrom);
pDist = m_pChromDists;
while(pDist != NULL)
	{
	BuffOfs += sprintf(&szBuffer[BuffOfs],",\"%s\"",pDist->szChrom);
	pDist = pDist->pNext;
	if(BuffOfs + 250 > sizeof(szBuffer))
		{
		CUtility::SafeWrite(m_hRsltsFile,szBuffer,BuffOfs);
		BuffOfs = 0;
		}
	}

// now for the counts
for(CntIdx = 0; CntIdx <= MaxLength; CntIdx++)
	{
	BuffOfs += sprintf(&szBuffer[BuffOfs],"\n%d,%d",CntIdx,m_pGenomeDists->Dist[CntIdx]);
	pDist = m_pChromDists;
	while(pDist != NULL)
		{
		BuffOfs += sprintf(&szBuffer[BuffOfs],",%d",pDist->Dist[CntIdx]);
		pDist = pDist->pNext;
		if(BuffOfs + 250 > sizeof(szBuffer))
			{
			CUtility::SafeWrite(m_hRsltsFile,szBuffer,BuffOfs);
			BuffOfs = 0;
			}
		}
	}
if(BuffOfs > 0)
	CUtility::SafeWrite(m_hRsltsFile,szBuffer,BuffOfs);
return(0);
}



int
ReportDist(int MinLength,				// minimum element length
		   int MaxLength,				// maximum element length
		   char *pszRsltsFile)
{
int Rslt;
if(m_pChromDists == NULL)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"No count distributions to report on");
	return(-1);
	}
#ifdef _WIN32
if((m_hRsltsFile = open(pszRsltsFile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
if((m_hRsltsFile = open(pszRsltsFile, O_RDWR | O_CREAT | O_TRUNC, S_IREAD | S_IWRITE ))==-1)
#endif
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create results file: %s - %s",pszRsltsFile,strerror(errno));
	Reset();
	return(eBSFerrCreateFile);
	}
Rslt = ReportChromDists(MinLength,MaxLength);
if(m_hRsltsFile != -1)
	{
	close(m_hRsltsFile);
	m_hRsltsFile = -1;
	}
return(Rslt);
}



int Process(etPMode PMode,				// processing mode
			char Strand,				// which element strand to filter (retain) on
			char *pszInLociFile,		// input CSV or BED loci file
			char *pszInBEDFile,			// input BED feature/gene file
			char *pszRsltsFile,			// output loci file
			int RegRegionLen,			// regulatory region length
			int MinLength,				// minimum element length
			int MaxLength)				// maximum element length
{
int Rslt;
Init();
m_PMode = PMode;
m_Strand = Strand;
m_RegRegionLen = RegRegionLen;

if(pszInBEDFile != NULL && pszInBEDFile[0] != '\0')
	{
	if((m_pBiobed = new CBEDfile) == NULL)
		{
		printf("\nUnable to instantiate CBEDfile gene/exon file '%s'",pszInBEDFile);
		Reset();
		return(eBSFerrObj);
		}

	if((Rslt=m_pBiobed->Open(pszInBEDFile,eBTGeneExons))!=eBSFSuccess)
		{
		while(m_pBiobed->NumErrMsgs())
			gDiagnostics.DiagOut(eDLFatal,gszProcName,m_pBiobed->GetErrMsg());
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open feature/gene file '%s'",pszInBEDFile);
		Reset();
		return(eBSFerrObj);
		}
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Opened features file '%s'",pszInBEDFile);
	}
else
	m_pBiobed = NULL;

// now load csv or bed input element file
switch(PMode) {
	case ePMdefault:		// default is for CSV loci processing
		Rslt = ProcessElsCSV(Strand,MinLength,MaxLength,pszInLociFile);
		break;

	case ePMBed:			// UCSC BED processing
		Rslt = ProcessElsBED(Strand,MinLength,MaxLength,pszInLociFile);
		break;

	default:
		break;
	}

// now report the distributions

Rslt = ReportDist(MinLength,				// minimum element length
		   MaxLength,				// maximum element length
		   pszRsltsFile);
Reset();
return(Rslt);
}


