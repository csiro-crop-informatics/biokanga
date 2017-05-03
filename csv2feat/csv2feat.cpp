// csv2feat.cpp : Defines the entry point for the console application.
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

const unsigned int cProgVer = 110;		// increment with each release
const int cDfltMinLengthRange = 4;		// default is to accept sequence lengths from this min length up
const int cMaxLengthRange = 1000000000;	// max feature length accepted

const int cMaxOutBuff = 32000;			// max number of chars to buffer in output 

// processing modes
typedef enum eProcMode {
	ePMStandard = 0			  // standard processing
} etProcMode;

int 
Process(etProcMode Mode,			// processing mode
		int MinLength,				// probe elements must be of at least this length
		int MaxLength,				// and no longer than this length
		int MinOverlap,				// overlap onto features must be at least this length
		char *pszInLociFile,		// CSV file containing elements to be mapped onto feature
		char *pszInFeatFile,		// file containing features to map elements onto
		char *pszRsltsFile);		// file to write results into

int								// eBSFSuccess, eBSFerrFeature or eBSFerrChrom
Loci2Gene(etProcMode Mode,		// processing mode
		char Strand,			// loci are on this strand
		CBEDfile *pBED,			// BED file containing features or genes to map loci onto
		char *pszChrom,			// chromosome
		int LociStart,			// loci start
		int LociEnd,			// loci end
		int MinOverlap,			// minimum number of bases required to overlap
		char *pszGene);			// where to return gene loci maps onto

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

int iMode;					// processing mode 0:standard
int iMinLength;				// probe elements must be of at least this length
int iMaxLength;				// and no longer than this length
int iMinOverlap;			// overlaps onto elements must be of at least this length

char szInLociFile[_MAX_PATH];	// input element loci from this file
char szInFeatFile[_MAX_PATH];	// input bioseq file features
char szRsltsFile[_MAX_PATH];	// output stats to this file


// command line args
struct arg_lit  *help    = arg_lit0("hH","help",                "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_int *ScreenLogLevel=arg_int0("S", "ScreenLogLevel",	"<int>","Level of diagnostics written to logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"diagnostics log file");

struct arg_int  *Mode = arg_int0("m","procmode","<int>",		"processing mode 0:standard");
struct arg_file *InLociFile = arg_file1("i","inloci","<file>",	"element loci CSV file");
struct arg_file *InFeatFile = arg_file1("I","feat","<file>",	"bioseq feature file");
struct arg_file *RsltsFile = arg_file1("o","output","<file>",	"output file");
struct arg_int  *MinLength = arg_int0("l","minlength","<int>",	"minimum element length (default 4)");
struct arg_int  *MaxLength = arg_int0("L","maxlength","<int>",	"maximum element length (default 1000000000)");
struct arg_int  *MinOverlap = arg_int0("M","minoverlap","<int>","minimum feature overlap (default 1)");
struct arg_end *end = arg_end(20);

void *argtable[] = {help,version,FileLogLevel,ScreenLogLevel,LogFile,
					Mode,InLociFile,InFeatFile,RsltsFile,MinLength,MaxLength,MinOverlap,
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
			printf("\n%s Version: %d.%2.2d\n",gszProcName,cProgVer/100,cProgVer%100);
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

	iMode = Mode->count ? Mode->ival[0] : ePMStandard;
	if(iMode < ePMStandard || iMode > ePMStandard)
		{
		printf("\nError: Requested processing mode '-m%d' not supported",iMode);
		exit(1);
		}

	iMinLength = MinLength->count ? MinLength->ival[0] : cDfltMinLengthRange;
	if(iMinLength < 1 || iMinLength > cMaxLengthRange)
		{
		printf("Error: Mininum element length '-l%d' is not in range 1..%d",iMinLength,cMaxLengthRange);
		exit(1);
		}

	iMaxLength = MaxLength->count ? MaxLength->ival[0] : cMaxLengthRange;
	if(iMaxLength < iMinLength || iMaxLength > cMaxLengthRange)
		{
		printf("Error: Maximum element length '-L%d' is not in range %d..%d",iMaxLength,iMinLength,cMaxLengthRange);
		exit(1);
		}


	iMinOverlap = MinOverlap->count ? MinOverlap->ival[0] : 1;
	if(iMinOverlap < 1 || iMinOverlap > iMinLength)
		{
		printf("Error: Minimum feature overlap length '-M%d' is not in range 1..%d",iMinOverlap,iMinLength);
		exit(1);
		}

	strncpy(szInLociFile,InLociFile->filename[0],_MAX_PATH);
	szInLociFile[_MAX_PATH-1] = '\0';
	strncpy(szInFeatFile,InFeatFile->filename[0],_MAX_PATH);
	szInFeatFile[_MAX_PATH-1] = '\0';
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
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Version: %d.%2.2d Processing parameters:",cProgVer/100,cProgVer%100);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Processing mode: Standard");
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Input CSV element loci file: '%s'",szInLociFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Input bioseq feature file: '%s'",szInFeatFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Output to file: '%s'",szRsltsFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Minimum element length: %d",iMinLength);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Maximum element length: %d",iMaxLength);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Minimum feature overlap length: %d",iMinOverlap);


#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	// processing here...
	Rslt = Process((etProcMode)iMode,iMinLength,iMaxLength,iMinOverlap,szInLociFile,szInFeatFile,szRsltsFile);

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
}

int 
Process(etProcMode Mode,			// processing mode
		int MinLength,				// core elements must be of at least this length
		int MaxLength,				// and no longer than this length
		int MinOverlap,				// must overlap by at least this number of bases onto feature
		char *pszInLociFile,		// CSV file containing elements
		char *pszInFeatFile,		// file containing features
		char *pszRsltsFile)			// file to write fasta into
{
int NumFields;

int Rslt;
int SrcID;
char *pszChrom;
char *pszElType;
char *pszRefSpecies;
char *pszFieldVal;
char szFeatHit[80];
char szLineBuff[cMaxOutBuff];
int BuffLen;

int StartLoci;
int EndLoci;
int Len;

int NumAccepted;
int NumProcessed;
int NumUnmapped;
int NumUnderLen;
int NumOverLen;

int hRsltFile = -1;
CBEDfile *pFeatFile = NULL;

CCSVFile *pCSV = new CCSVFile;
if(pCSV == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CCSVfile");
	return(eBSFerrObj);
	}

if((Rslt=pCSV->Open(pszInLociFile))!=eBSFSuccess)
	{
	while(pCSV->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,pCSV->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open file: %s",pszInLociFile);
	delete pCSV;
	return(Rslt);
	}


if((pFeatFile = new CBEDfile) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CBEDfile object");
	delete pCSV;
	return(eBSFerrObj);
	}

if((Rslt = pFeatFile->Open(pszInFeatFile))!=eBSFSuccess)
	{
	while(pFeatFile->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,pFeatFile->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open bioseq feature file '%s'",pszInFeatFile);
	delete pCSV;
	delete pFeatFile;
	return(Rslt);
	}

#ifdef _WIN32
if((hRsltFile = open(pszRsltsFile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
if((hRsltFile = open(pszRsltsFile, O_RDWR | O_CREAT |O_TRUNC, S_IREAD | S_IWRITE))==-1)
#endif
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create %s - %s",pszRsltsFile,strerror(errno));
	delete pCSV;
	delete pFeatFile;
	return(eBSFerrCreateFile);
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Output file created/truncated: '%s'",pszRsltsFile);

NumUnmapped = 0;
NumUnderLen = 0;
NumOverLen = 0;
NumAccepted = 0;
NumProcessed = 0;
BuffLen = 0;
while((Rslt=pCSV->NextLine()) > 0)	// onto next line containing fields
	{
	NumFields = pCSV->GetCurFields();
	if(NumFields < 7)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expected 7+ fields in '%s', GetCurFields() returned '%d'",pszInLociFile,NumFields);
		Rslt = eBSFerrFieldCnt;
		break;
		}
	if(!NumProcessed && pCSV->IsLikelyHeaderLine())
		continue;

	NumProcessed += 1;

	pCSV->GetInt(7,&Len);
	if(Len < MinLength)
		{
		NumUnderLen += 1;
		continue;
		}
		
	if(Len > MaxLength)
		{
		NumOverLen += 1;
		continue;
		}

	pCSV->GetInt(1,&SrcID);
	pCSV->GetText(2,&pszElType);
	pCSV->GetText(3,&pszRefSpecies);
	pCSV->GetText(4,&pszChrom);
	pCSV->GetInt(5,&StartLoci);
	pCSV->GetInt(6,&EndLoci);

	Rslt = Loci2Gene(Mode,'*',pFeatFile,pszChrom,StartLoci,EndLoci,MinOverlap,szFeatHit);
	if(Rslt == eBSFerrFeature)		// doesn't map onto any feature
		{
		NumUnmapped += 1;
		continue;
		}
	if(Rslt < eBSFSuccess)
		break;

	BuffLen += sprintf(&szLineBuff[BuffLen],"%d,\"%s\",\"%s\",\"%s\",%d,%d,%d",
		SrcID,pszElType,pszRefSpecies,pszChrom,StartLoci,EndLoci,Len);
	if(BuffLen > sizeof(szLineBuff)/2)
		{
		CUtility::SafeWrite(hRsltFile,szLineBuff,BuffLen);
		BuffLen = 0;
		}
	for(int FieldIdx = 8; FieldIdx <= NumFields; FieldIdx++)
		{
		pCSV->GetText(FieldIdx,&pszFieldVal);
		if(pCSV->GetQuoted(FieldIdx))
			BuffLen += sprintf(&szLineBuff[BuffLen],",\"%s\"",pszFieldVal);
		else
			BuffLen += sprintf(&szLineBuff[BuffLen],",%s",pszFieldVal);
		if(BuffLen > sizeof(szLineBuff)/2)
			{
			CUtility::SafeWrite(hRsltFile,szLineBuff,BuffLen);
			BuffLen = 0;
			}
		}
	BuffLen += sprintf(&szLineBuff[BuffLen],",\"%s\"\n",szFeatHit);
	
	if(BuffLen > sizeof(szLineBuff)/2)
		{
		CUtility::SafeWrite(hRsltFile,szLineBuff,BuffLen);
		BuffLen = 0;
		}
	NumAccepted += 1;
	}

if(Rslt >= eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Elements accepted: %d, Processed: %d, Unmapped: %d,UnderLen: %d, OverLen: %d",
			NumAccepted,NumProcessed,NumUnmapped,NumUnderLen,NumOverLen);
	}

if(BuffLen)
	CUtility::SafeWrite(hRsltFile,szLineBuff,BuffLen);

close(hRsltFile);
delete pCSV;
delete pFeatFile;
return(Rslt < 0 ? NumAccepted : Rslt);
}


int								// eBSFSuccess, eBSFerrFeature or eBSFerrChrom
Loci2Gene(etProcMode Mode,		// processing mode
		char Strand,			// loci are on this strand
		CBEDfile *pBED,			// BED file containing features or genes to map loci onto
		char *pszChrom,			// chromosome
		int LociStart,			// loci start
		int LociEnd,			// loci end
		int MinOverlap,			// minimum number of bases required to overlap
		char *pszGene)			// where to return gene loci maps onto
{
int Ith;
int ChromID;
int FeatID;
int InGeneFeatID;

int GeneStart;
int GeneEnd;
int InGeneStart;
int InGeneEnd;

int OverlapStart;
int OverlapEnd;

// length must be at least the min overlap!
if((1 + LociEnd - LociStart) < MinOverlap)
	return(eBSFerrFeature);

*pszGene = '\0';

// is the chromosome known?
if((ChromID = pBED->LocateChromIDbyName(pszChrom)) < 1)
	return(eBSFerrFeature);

// processing may be strand specific...
pBED->SetStrand(Strand);

FeatID = 0;
InGeneFeatID = 0;

// if loci is intragenic then...
if(pBED->InAnyFeature(ChromID,LociStart,LociEnd))
	{
	// iterate over all genes loci maps onto and use shortest gene completely containing loci
	Ith = 1;
	do {
		FeatID = pBED->LocateFeatureIDinRangeOnChrom(ChromID,LociStart,LociEnd,Ith++);
		if(FeatID < 1)
			continue;
		pBED->GetFeature(FeatID,NULL,NULL,&GeneStart,&GeneEnd);
				
		if((1 + GeneEnd - GeneStart) < MinOverlap)	// feature must be at least min overlap length
			continue;

		// check if at least minimum overlap onto feature
		if(GeneStart < LociStart)
			OverlapStart = LociStart;
		else
			OverlapStart = GeneStart;

		if(GeneEnd > LociEnd)
			OverlapEnd = LociEnd;
		else
			OverlapEnd = GeneEnd;

		if((1 + OverlapEnd - OverlapStart) < MinOverlap)
			continue;

		// go for smallest overlapped feature
		if(InGeneFeatID == 0)	// if 1st feature to overlap
			{
			InGeneFeatID = FeatID;
			InGeneStart = GeneStart;
			InGeneEnd = GeneEnd;
			}
		else
			{
			if(GeneStart <= LociStart && GeneEnd >= LociEnd)
				{
				if(GeneStart >= InGeneStart && GeneEnd <= InGeneEnd)
					{
					InGeneFeatID = FeatID;
					InGeneStart = GeneStart;
					InGeneEnd = GeneEnd;
					}
				}
			}
		}
	while(FeatID > 0);
	}

if(InGeneFeatID)
	FeatID = InGeneFeatID;
else
	return(eBSFerrFeature);	// unable to associate	

pBED->GetFeature(FeatID,pszGene,NULL,&GeneStart,&GeneEnd);
return(eBSFSuccess);
}
