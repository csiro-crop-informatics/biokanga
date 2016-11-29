// gensampler.cpp : Defines the entry point for the console application.
// Given a csv file containing loci + len + region will sample the specified genome assembly
// and generate a csv file containing randomly selected loci with same lengths and region as the original
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

const unsigned int cProgVer = 203;		// increment with each release


const int cMaxRandAttempts = 10000000; // max attempts at randomly selecting chromosomal locus in same functional region
const int cAllocSeqLen = 100000;	 // allocate for sequences in this sized chunks

typedef enum eProcMode {
	eProcModeStandard = 0				// default processing
} etProcMode;

int ParseRegions(char *pszRegionList);
char *Regions2Txt(int Regions);
int TrimQuotes(char *pszTxt);

int Process(int Mode,				// processing mode
		int FilterRegionsIn,		// retain any of these (exclusive) regions
		int RegLen,					// length of regulatory US/DS regions
		char *szInFile,				// input csv file containing loci + lengths
		char *pszInputBiobedFile,	// characterise regions from this file
		char *pszInSeqFile,			// use sequences in this file
		char *pszOutputFile);		// output to this file

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

int iMode;
char szOutputFile[_MAX_PATH];
char szInputFile[_MAX_PATH];
char szInSeqFile[_MAX_PATH];
char szInputBiobedFile[_MAX_PATH];
int iRegionsIn;
int iRegionsOut;
char szRegionsIn[128];		// regions to retain
char szRegionsOut[128];		// regions to exclude
int iRegLen;

// command line args
struct arg_lit  *help    = arg_lit0("hH","help",                "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_int *ScreenLogLevel=arg_int0("S", "ScreenLogLevel",	"<int>","Level of diagnostics written to logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",	"diagnostics log file");

struct arg_int  *Mode = arg_int0("m","procmode","<int>",	"processing mode 0:default");
struct arg_int  *RegLen = arg_int0("L","updnstream","<int>",	"length of 5'up or 3'down  stream regulatory region length (default = 2000) 0..1000000");
struct arg_file *InFile = arg_file1("i",NULL,"<file>",			"input from loci csv file");
struct arg_file *OutFile = arg_file1("o",NULL,"<file>",			"output to multifasta file");
struct arg_file *InSeqFile = arg_file1("a","assemb","<file>","get sequences fron this file");
struct arg_file *InBedFile = arg_file1("b","bed","<file>",		"characterise regions from biobed file");
struct arg_str  *RegionsOut = arg_str0("R","regionsout","<string>","Filter out random samples in any of these regions (space or comma delimit), 1: Intergenic, 2: US, 3: 5'UTR, 4: CDS, 5: Intron, 6: 3'UTR, 7: DS");
struct arg_str  *RegionsIn = arg_str0("r","regionsin","<string>","Accept random samples in any of these exclusive regions (space or comma delimit), 1: Intergenic, 2: US, 3: 5'UTR, 4: CDS, 5: Intron, 6: 3'UTR, 7: DS");

struct arg_end *end = arg_end(20);

void *argtable[] = {help,version,FileLogLevel,ScreenLogLevel,LogFile,
					Mode,RegLen,InFile,OutFile,InSeqFile,InBedFile,RegionsIn,RegionsOut,
					end};

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

	iMode = Mode->count ? Mode->ival[0] : eProcModeStandard;
	if(iMode < eProcModeStandard || iMode > eProcModeStandard)
		{
		printf("\nError: Requested processing mode '-x%d' not supported",iMode);
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
		}
	else
		{
		szRegionsIn[0] = '\0';
		iRegionsIn = (cFeatBitIG | cFeatBitUpstream | cFeatBit5UTR | cFeatBitCDS | cFeatBitIntrons | cFeatBit3UTR | cFeatBitDnstream);
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
	iRegionsIn &= ~iRegionsOut;
	if(!iRegionsIn)
		{
		printf("Error: no regions to retain");
		exit(1);
		}

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

	strcpy(szInputFile,InFile->filename[0]);
	strcpy(szOutputFile,OutFile->filename[0]);
	strcpy(szInputBiobedFile,InBedFile->filename[0]);
	strcpy(szInSeqFile,InSeqFile->filename[0]);

		// now that command parameters have been parsed then initialise diagnostics log system
	if(!gDiagnostics.Open(szLogFile,(etDiagLevel)iScreenLogLevel,(etDiagLevel)iFileLogLevel,true))
		{
		printf("\nError: Unable to start diagnostics subsystem.");
		if(szLogFile[0] != '\0')
			printf(" Most likely cause is that logfile '%s' can't be opened/created",szLogFile);
		exit(1);
		}

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Version: %d.%2.2d Processing parameters:",cProgVer/100,cProgVer%100);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Loci file (.csv) file to process: '%s'",szInputFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"bio assembly (.seq) file to process: '%s'",szInSeqFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"where to write out fasta with random sequences: '%s'",szOutputFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"biobed file containing regional features: '%s'",szInputBiobedFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Processing mode: %s",iMode == eProcModeStandard ? "standard" : iMode == eProcModeStandard ? "extended" : "summary");
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Accept random samples in any of these regions: '%s'",Regions2Txt(iRegionsIn));

#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif

	Rslt = Process(iMode,iRegionsIn,iRegLen,szInputFile,szInputBiobedFile,szInSeqFile,szOutputFile);
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


// ParseRegions
// Parses space or comma delimited list of regions in which
// 1 == Intergenic, 2 == US, 3 == 5'UTR, 4 == CDS, 5 == Intron, 6 == 3'UTR, 7 == DS
//
// Returns bitmap of regions or -1 if parse errors
// If no region specified then assumes all regions are selected
int
ParseRegions(char *pszRegionList)
{
// parse out region list
char Chr;
int Region = 0;
if(pszRegionList == NULL || *pszRegionList == '\0')
	return(cFeatBitIG | cFeatBitUpstream | cFeatBit5UTR | cFeatBitCDS | cFeatBitIntrons | cFeatBit3UTR | cFeatBitDnstream);

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
return(szRegions);
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



// OpenBedfile
// Attempts to open specified bedfile
// Returns ptr to opened bedfile or NULL
CBEDfile *
OpenBedfile(char *pToOpen)
{
int Rslt;
CBEDfile *pBed;
if(pToOpen != NULL && pToOpen[0] != '\0')
	{
	if((pBed = (CBEDfile *)new CBEDfile())==NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CBEDfile '%s'",pToOpen);
		return(NULL);
		}

	if((Rslt = pBed->Open(pToOpen,eBTAnyBed))!=eBSFSuccess)
		{
		while(pBed->NumErrMsgs())
			gDiagnostics.DiagOut(eDLFatal,gszProcName,pBed->GetErrMsg());
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open biobed file '%s'",pToOpen);
		delete pBed;
		return(NULL);
		}
	return(pBed);
	}
return(NULL);
}

int
MapBEDRegion(CBEDfile *pBed,int BEDChromID,int Start,int End,int UpDnStreamLen)
{
int FeatureBits;
int RegionIdx;
int BitMsk;
int FeatIdx;
FeatureBits = pBed->GetFeatureBits(BEDChromID,Start,End,cRegionFeatBits,UpDnStreamLen);
RegionIdx = 0;		// default to intergenic if no feature bits set
if(FeatureBits)		
	{
	BitMsk = cFeatBitCDS;
	for(FeatIdx = 1; FeatIdx < 7; FeatIdx++,BitMsk <<= 1)
		{
		if(BitMsk & FeatureBits)
			{
			RegionIdx = FeatIdx;
			break;
			}
		}
	switch(RegionIdx) {
		case 0: case 2: case 4: case 6:		// IG,5'UTR, Introns and 3'DS
			FeatureBits = RegionIdx;
			break;
		case 1:								// CDS
			FeatureBits = 3;
			break;
		case 3:								// 3'UTR
			FeatureBits = 5;
			break;
		case 5:								// 5'US
			FeatureBits = 1;
			break;
			
		}
	}
return(FeatureBits);
}

int
Process(int Mode,	// processing mode
		int FilterRegionsIn,		// retain any of these (exclusive) regions
		int RegLen,					// length of regulatory US/DS regions
		char *szInFile,				// input csv file containing loci + lengths
		char *pszInputBiobedFile,	// characterise regions from this file
		char *pszInSeqFile,			// use sequences in this file
		char *pszOutputFile)		// output to this file
{
int Rslt;
int NumEls;
CBEDfile *pBed;
CBioSeqFile *pSeq;
CCSVFile *pCSV;
char *pszSpecies;
char *pszChrom;
char szSpecies[64];
char szCurChrom[64];
INT64 Now;
int RandSeed;
int NumFields;
int CoreID;
int SeqStart;
int SeqLen;
int SeqChromLen;
int SeqChromID;
int BEDChromID;
int BEDStart;
int BEDEnd;
int BEDRegion;

int FiltRegion;

int Attempts;

int AllocSeqLen;
etSeqBase *pSeqBuffer;
etSeqBase *pBase;
int Cnt;

char szLineBuff[512];
int LineLen;
int hRsltsFile;

if((pBed = (CBEDfile *)new CBEDfile())==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"\nUnable to instantiate CBEDfile '%s'",pszInputBiobedFile);
	return(eBSFerrObj);
	}

if((Rslt = pBed->Open(pszInputBiobedFile,eBTAnyBed))!=eBSFSuccess)
	{
	while(pBed->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,pBed->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open biobed file %s",pszInputBiobedFile);
	delete pBed;
	return(Rslt);
	}

if((pSeq = (CBioSeqFile *)new CBioSeqFile())==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"\nUnable to instantiate CBioSeqFile '%s'",pszInSeqFile);
	delete pBed;
	return(eBSFerrObj);
	}

if((Rslt = pSeq->Open(pszInSeqFile,cBSFTypeSeq))!=eBSFSuccess)
	{
	while(pSeq->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,pSeq->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open bioseq file '%s'",pszInSeqFile);
	delete pBed;
	delete pSeq;
	return(Rslt);
	}

if((pCSV = new CCSVFile())==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"\nUnable to instantiate CCSVFile '%s'",szInFile);
	delete pBed;
	delete pSeq;
	return(eBSFerrObj);
	}
if((Rslt=pCSV->Open(szInFile))!=eBSFSuccess)
	{
	while(pCSV->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,pCSV->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open CSV loci file: '%s'",szInFile);
	delete pBed;
	delete pSeq;
	delete pCSV;
	return(Rslt);
	}

#ifdef _WIN32
if((hRsltsFile = open(pszOutputFile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
if((hRsltsFile = open(pszOutputFile, O_RDWR | O_CREAT | O_TRUNC, S_IREAD | S_IWRITE ))==-1)
#endif
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create %s - %s",pszOutputFile,strerror(errno));
		delete pBed;
		delete pSeq;
		delete pCSV;
		return(eBSFerrCreateFile);
		}

// all input files now opened


#ifdef _WIN32
QueryPerformanceCounter((LARGE_INTEGER *)&Now);
#else
struct tms Times;
Now = (INT64)times(&Times);
#endif
RandSeed = (int)(Now & 0x07fffffff);
CRandomMersenne RandBase(RandSeed);

pSeqBuffer = NULL;
AllocSeqLen = 0;

szCurChrom[0] = '\0';
// no headers, straight into the rows which are expected to contain at least 7 fields
NumEls = 0;
while((Rslt=pCSV->NextLine()) > 0)	// onto next line containing fields
	{
	NumFields = pCSV->GetCurFields();
	if(NumFields < 7)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expected at least 7 fields in '%s', GetCurFields() returned '%d'",szInFile,NumFields);
		close(hRsltsFile);
		delete pBed;
		delete pSeq;
		delete pCSV;
		return(eBSFerrFieldCnt);
		}
	if(!NumEls && pCSV->IsLikelyHeaderLine())
		continue;

	NumEls += 1;
	pCSV->GetInt(1,&CoreID);
	pCSV->GetText(3,&pszSpecies);
	strcpy(szSpecies,pszSpecies);
	pCSV->GetText(4,&pszChrom);
	pCSV->GetInt(5,&SeqStart);
	pCSV->GetInt(7,&SeqLen);

	// get chromome details when chromosome changes
	if(szCurChrom[0] == '\0' || stricmp(szCurChrom,pszChrom))
		{
		if((SeqChromID = pSeq->LocateEntryIDbyName(pszChrom)) <= 0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Can't locate expected chromosome '%s' in sequence file '%s'",pszChrom,pszInSeqFile);
			close(hRsltsFile);
			delete pBed;
			delete pSeq;
			delete pCSV;
			return(eBSFerrFieldCnt);
			}

		if((BEDChromID = pBed->LocateChromIDbyName(pszChrom)) <= 0)
			{
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"Can't locate expected chromosome '%s' in bed file '%s'",pszChrom,pszInputBiobedFile);
			BEDChromID = 0;
			}
		SeqChromLen = pSeq->GetDataLen(SeqChromID);
		strcpy(szCurChrom,pszChrom);
		}

	if(pSeqBuffer == NULL || SeqLen > AllocSeqLen)
		{
		if(pSeqBuffer != NULL)
			delete pSeqBuffer;
		if((pSeqBuffer = new unsigned char [SeqLen + cAllocSeqLen])==NULL)
			break;
		AllocSeqLen = SeqLen + cAllocSeqLen;
		}

	// randomly select start ofs on chromosome until requested genomic region located
	for(Attempts = 0; Attempts < cMaxRandAttempts; Attempts++)
		{
		BEDStart = RandBase.IRandom(0, SeqChromLen - SeqLen - 1);		// randomly select initial start on genome chromosome
		BEDEnd = BEDStart + SeqLen - 1;
		if(BEDChromID > 0)
			{
			BEDRegion = MapBEDRegion(pBed,BEDChromID,BEDStart,BEDEnd,RegLen);

			switch(BEDRegion) {
				case 0:	// intergenic
					FiltRegion = cFeatBitIG;
					break;
				case 1:	// US
					FiltRegion = cFeatBitUpstream;
					break;
				case 2:	// 5'UTR
					FiltRegion = cFeatBit5UTR;
					break;
				case 3:	// CDS
					FiltRegion = cFeatBitCDS;
					break;
				case 4:	// Intron
					FiltRegion = cFeatBitIntrons;
					break;
				case 5: // 3'UTR
					FiltRegion = cFeatBit3UTR;
					break;
				case 6:	// DS
					FiltRegion = cFeatBitDnstream;
					break;
				}
			}
		else
			FiltRegion = 0;

		if(BEDChromID == 0 || (FiltRegion & FilterRegionsIn))
			{
			// check if sequence contains any undefined bases, if so then need to find another sequence
			if((pSeq->GetData(SeqChromID,eSeqBaseType,BEDStart,pSeqBuffer,SeqLen)) < 0)
				continue;

			pBase = pSeqBuffer;
			Cnt = SeqLen;
			while(Cnt--)
				if(*pBase++ > eBaseT)
					break;
			if(Cnt < 0)
				break;
			}
		}

	if(Attempts == cMaxRandAttempts)
		{
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Unable to get random genomic region on chrom %s",szCurChrom);
		continue;
		}

	// have seq in requested genomic region
	LineLen = sprintf(szLineBuff,"%d,\"%s\",\"%s\",\"%s\",%d,%d,%d,\"%s\",%d\n",
		CoreID,"randcore",szSpecies,szCurChrom,BEDStart,BEDEnd,SeqLen,"RandomSelection",BEDRegion);
	CUtility::SafeWrite(hRsltsFile,szLineBuff,LineLen);
	}

if(hRsltsFile != -1)
	close(hRsltsFile);

if(pSeqBuffer != NULL)
	delete pSeqBuffer;
if(pBed != NULL)
	delete pBed;
if(pSeq != NULL)
	delete pSeq;
if(pCSV != NULL)
	delete pCSV;
return(0);
}


