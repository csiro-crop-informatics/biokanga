// genelementseq.cpp : Defines the entry point for the console application.
// Processes an input CSV file containing element loci, and outputs an extended CSV file
// with original loci plus the nucleotide sequence at that loci
// Version 2.xxxx
// Added modes whereby output consists of a fasta file containing either all
// element sequences simply concatenated into a single sequence, or as individual
// fasta sequences
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

const unsigned int cProgVer = 206;		// increment with each release

const int cMaxLengthRange = 1000000;	// maximal element length
const int cSeqAllocLen = 1000000;		// can handle sequences of upto this maximal length


typedef enum TAG_eRsltsFormat {
	eRsltsFCSV = 0,		// CSV output
	eRsltsFasta,		// all sequences concatenated as single fasta record
	eRsltsMultifasta	// sequences as multifasta
} teRsltsFormat;



int 
Process(teCSVFormat CSVFormat,	// expected input CSV format
		teRsltsFormat RsltsFormat,	//processing mode (default 0) 0: extended CSV file, 1: concatenated fasta, 2: multifasta 
		bool bFeatures,		// if true then generate feature detail from pszInBEDFile (file containing gene id's + loci)
		int RegRegionLen,	// up/down stream regulatory region length (only applies if bFeatures == true)
		int MinLength,		// only process elements of at least this length
		int MaxLength,		// only process elements which are no longer than this length
		int MinIdentity,	// minimum identity (0..100)
		int MaxIdentity,	// maximum identity (MinIdentity..100)
		char *pszInLociFile,// CSV file containing elements
		char *pszInBEDFile,	// file containing gene id's + loci
		char *pszInSeqFile,  // file containing genome assembly sequences
		char *pszRsltsFile);// file to write results into

char *RsltsFormat2Text(teRsltsFormat Format);
char *CSVFormat2Text(teCSVFormat Format);

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
	return _T("genelementseq");
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
int iCSVFormat;
int iRsltsFormat;
int iMinLength;				// core elements must be of at least this length
int iMaxLength;				// and no longer than this length
int iMinIdentity;			// out species must align with at least this identity
int iMaxIdentity;			// and no more than this identity
bool bFeatures;					// true if overlapping features to be located
int iRegLen;					// up/down stream regulatory region length
char szInLociFile[_MAX_PATH];	// input element loci from this file
char szInBEDFile[_MAX_PATH];	// input bed file containing gene identifiers
char szInSeqFile[_MAX_PATH];	// input bioseq file containing assembly
char szRsltsFile[_MAX_PATH];	// output loci + sequences or gene identifiers to this file


// command line args
struct arg_lit  *help    = arg_lit0("hH","help",                "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_int *ScreenLogLevel=arg_int0("S", "ScreenLogLevel",	"<int>","Level of diagnostics written to logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"diagnostics log file");

struct arg_int  *CSVFormat = arg_int0("c","informat","<int>",	"input CSV file type 0:Loci 1:locsfx probe 2:locsfx target 3: Genhyperconserved x2 4: GenMAlignScore m3");
struct arg_file *InLociFile = arg_file1("i","inloci","<file>",	"element CSV file");
struct arg_file *InBEDFile = arg_file0("I","inbed","<file>",	"gene biobed BED file");
struct arg_file *InSeqFile = arg_file1("a","assembly","<file>",	"genome assembly bioseq file");
struct arg_file *RsltsFile = arg_file1("o","output","<file>",	"output file");
struct arg_int  *RsltsFormat = arg_int0("p","outformat","<int>","results output as 0: extended CSV file, 1: concatenated fasta, 2: multifasta ");
struct arg_int  *MinLength = arg_int0("m","minlength","<int>",	"minimum element length (default 0)");
struct arg_int  *MaxLength = arg_int0("M","maxlength","<int>",	"maximum element length (default 1000000)");
struct arg_int *MinIdentity = arg_int0("d", "minident","<int>",	"input type 3 minimum out species alignment identity (default 0)");
struct arg_int *MaxIdentity = arg_int0("D", "maxident","<int>",	"input type 3 maximum out species alignment identity (default 100)");

struct arg_int  *RegLen = arg_int0("L","updnstream","<int>",	"length of 5'up or 3'down  stream regulatory region length (default = 2000) 0..1000000");
struct arg_lit  *Features = arg_lit0("x","features",			"regenerate any existing feature mappings");

struct arg_end *end = arg_end(20);

void *argtable[] = {help,version,FileLogLevel,ScreenLogLevel,LogFile,
					CSVFormat,InLociFile,InBEDFile,InSeqFile,RsltsFile,RsltsFormat,MinLength,MaxLength,MinIdentity,MaxIdentity,
					RegLen,Features,
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

	iCSVFormat = CSVFormat->count ? CSVFormat->ival[0] : eCSVFdefault;
	if(iCSVFormat < eCSVFdefault || iCSVFormat > eCSVFAlignM3)
		{
		printf("\nError: expected input CSV format specified '-c%d' must be in range 0..4",iCSVFormat);
		exit(1);
		}

	iRsltsFormat = RsltsFormat->count ? RsltsFormat->ival[0] : eRsltsFCSV;
	if(iRsltsFormat < eRsltsFCSV || iRsltsFormat > eRsltsMultifasta)
		{
		printf("\nError: RsltsFormat specified '-p%d' must be in range 0..2",iRsltsFormat);
		exit(1);
		}

	iMinLength = MinLength->count ? MinLength->ival[0] : 0;
	if(iMinLength < 0 || iMinLength > cMaxLengthRange)
		{
		printf("Error: Mininum element length '-m%d' is not in range 0..%d",iMinLength,cMaxLengthRange);
		exit(1);
		}

	iMaxLength = MaxLength->count ? MaxLength->ival[0] : cMaxLengthRange;
	if(iMaxLength < iMinLength || iMaxLength > cMaxLengthRange)
		{
		printf("Error: Maximum element length '-M%d' is not in range %d..%d",iMaxLength,iMinLength,cMaxLengthRange);
		exit(1);
		}

	if(iCSVFormat == eCSVFhyperX2)
		{
		iMinIdentity = MinIdentity->count ? MinIdentity->ival[0] : 0;
		if(iMinIdentity < 0 || iMinIdentity > 100)
			{
			printf("\nError: Minimum identity specified as '-d%d' must be in range 0..100\n",iMinIdentity);
			exit(1);
			}

		iMaxIdentity = MaxIdentity->count ? MaxIdentity->ival[0] : 100;
		if(iMaxIdentity < iMinIdentity || iMaxIdentity > 100)
			{
			printf("\nError: Maximum identity specified as '-D%d' must be in range %d..100\n",iMaxIdentity,iMinIdentity);
			exit(1);
			}
		szInBEDFile[0] = '\0';
		bFeatures = false;
		iRegLen = 0;
		}
	else
		{
		iMinIdentity = 0;
		iMaxIdentity = 100;
		}

	if(InBEDFile->count)
		{
		strncpy(szInBEDFile,InBEDFile->filename[0],_MAX_PATH);
		szInBEDFile[_MAX_PATH-1] = '\0';
		bFeatures = Features->count ? true : false;
		if(bFeatures)
			{
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
			iRegLen = 0;
		}
	else
		{
		szInBEDFile[0] = '\0';
		bFeatures = false;
		iRegLen = 0;
		}

	strncpy(szInLociFile,InLociFile->filename[0],_MAX_PATH);
	szInLociFile[_MAX_PATH-1] = '\0';
	strncpy(szInSeqFile,InSeqFile->filename[0],_MAX_PATH);
	szInSeqFile[_MAX_PATH-1] = '\0';
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
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Input CSV element loci file: '%s'",szInLociFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Expecting CSV format as: %s",CSVFormat2Text((teCSVFormat)iCSVFormat));
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Input bioseq genome assembly file: '%s'",szInSeqFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Output to file: '%s'",szRsltsFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Generating Results as: %s",RsltsFormat2Text((teRsltsFormat)iRsltsFormat));
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Minimum ref element length: %d",iMinLength);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Maximum ref element length: %d",iMaxLength);

	if(iCSVFormat == eCSVFhyperX2)
		{
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Minimum out species identity: %d",iMinIdentity);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Maximum out species identity: %d",iMaxIdentity);
		}

	if(szInBEDFile[0] != '\0')
		{
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Input biobed BED gene loci file: '%s'",szInBEDFile);
		if(bFeatures)
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Regenerate feature mappings with regulatory region length: %d",iRegLen);
		}

#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	// processing here...
	Rslt = Process((teCSVFormat)iCSVFormat,(teRsltsFormat)iRsltsFormat,bFeatures,iRegLen,iMinLength,iMaxLength,iMinIdentity,iMaxIdentity,szInLociFile,szInBEDFile,szInSeqFile,szRsltsFile);

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

char *
RsltsFormat2Text(teRsltsFormat Format)
{
switch(Format) {
	case eRsltsFCSV:	// CSV output
		return((char *)"CSV with sequence");
	case eRsltsFasta:	// all sequences concatenated as single fasta record
		return((char *)"Sequences concatenated into single Fasta record");
	case eRsltsMultifasta:	// sequences as multifasta
		return((char *)"Sequences as multiple fasta records");
	default:
		break;
	}
return((char *)"Unsupported");
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
	case eCSVFhyperX2:
		return((char *)"Genhyperconserved with -x2 option");
	case eCSVFAlignM3:
		return((char *)"GenMAlignScore with -m3 option");
	default:
		break;
	}
return((char *)"Unsupported");
}

int 
Process(teCSVFormat CSVFormat,		// expected input CSV format
		teRsltsFormat RsltsFormat,	// processing mode (default 0) 0: extended CSV file, 1: concatenated fasta, 2: multifasta 
		bool bFeatures,		// if true then generate feature detail from pszInBEDFile (file containing gene id's + loci)
		int RegRegionLen,	// up/down stream regulatory region length (only applies if bFeatures == true)
		int MinLength,				// core elements must be of at least this length
		int MaxLength,				// and no longer than this length
		int MinIdentity,			// out species must align with at least this identity
		int MaxIdentity,			// and no more than this identity
		char *pszInLociFile,// CSV file containing elements
		char *pszInBEDFile,	// file containing gene id's + loci
		char *pszInSeqFile,  // file containing genome assembly sequences
		char *pszRsltsFile)// file to write results into
{
int NumFields;
int NumElsRead;
int NumEls = 0;
int Rslt;
int SrcID;
char *pszChrom;
int ChromID;
char *pszElType;
char *pszRefSpecies;
char *pszRelSpecies;
int StartLoci;
int EndLoci;
int Len;
int Features;

int OGUnalignedBases;
int OGMatchCnt;
int OGMismatchCnt;
int OGInDelCnt;

int Identity;

char szGeneName[100];
int FeatID;
int NxtFeatID;
int NxtFeatStart;
int NxtFeatEnd;
int PrvFeatID;
int PrvFeatStart;
int PrvFeatEnd;

int NumCols;
int SeqOfs;
int NxtFastaCol;
bool bFastaStart;
char *pszLineBuff;
int BuffLen;
etSeqBase *pElSeqBuff;
int hRsltFile = -1;
CBEDfile *pBiobed = NULL;
CBioSeqFile *pBioseq = NULL;

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

if(pszInBEDFile != NULL && pszInBEDFile[0] != '\0')
	{
	if((pBiobed = new CBEDfile) == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CBEDfile object");
		delete pCSV;
		return(eBSFerrObj);
		}

	if((Rslt=pBiobed->Open(pszInBEDFile,eBTGeneExons))!=eBSFSuccess)
		{
		while(pBiobed->NumErrMsgs())
			gDiagnostics.DiagOut(eDLFatal,gszProcName,pBiobed->GetErrMsg());
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open gene file '%s'",pszInBEDFile);
		delete pCSV;
		return(eBSFerrObj);
		}
	}
else
	pBiobed = NULL;

if((pBioseq = new CBioSeqFile) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CBioSeqFile object");
	delete pCSV;
	if(pBiobed != NULL)
		delete pBiobed;
	return(eBSFerrObj);
	}
if((Rslt = pBioseq->Open(pszInSeqFile))!=eBSFSuccess)
	{
	while(pBioseq->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,pBioseq->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open assembly sequence file '%s'",pszInSeqFile);
	delete pCSV;
	if(pBiobed != NULL)
		delete pBiobed;
	delete pBioseq;
	return(Rslt);
	}

if((pElSeqBuff = (etSeqBase *)new unsigned char[cSeqAllocLen])==NULL)
	{
	printf("\nUnable to allocate %d bytes as a sequence buffer",cSeqAllocLen);
	delete pCSV;
	if(pBiobed != NULL)
		delete pBiobed;
	delete pBioseq;
	return(eBSFerrMem);
	}

if((pszLineBuff = (char *)new char[(cSeqAllocLen * 3)/2])==NULL) // need to account for EOLs plus fasta descriptors etc
	{
	printf("\nUnable to allocate %d bytes as a line buffer",(cSeqAllocLen * 3)/2);
	delete pCSV;
	if(pBiobed != NULL)
		delete pBiobed;
	delete pBioseq;
	delete pElSeqBuff;
	return(eBSFerrMem);
	}

#ifdef _WIN32
if((hRsltFile = open(pszRsltsFile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
if((hRsltFile = open(pszRsltsFile, O_RDWR | O_CREAT |O_TRUNC,S_IREAD | S_IWRITE))==-1)
#endif
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create %s - %s",pszRsltsFile,strerror(errno));
	delete pCSV;
	if(pBiobed != NULL)
		delete pBiobed;
	delete pBioseq;
	delete pElSeqBuff;
	return(eBSFerrCreateFile);
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Output file created/truncated: '%s'",pszRsltsFile);
bFastaStart = true;
NumElsRead = 0;
while((Rslt=pCSV->NextLine()) > 0)	// onto next line containing fields
	{
	NumFields = pCSV->GetCurFields();
	if(NumFields > 6 && (!NumElsRead && pCSV->IsLikelyHeaderLine()))
		continue;
	NumElsRead+=1;
	switch(CSVFormat) {
		case eCSVFdefault:
			if(NumFields < 9)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expected 9+ fields in '%s', GetCurFields() returned '%d'",pszInLociFile,NumFields);
				Rslt = eBSFerrFieldCnt;
				break;
				}

			pCSV->GetInt(7,&Len);
			if(Len < MinLength || Len > MaxLength)
				continue;
			pCSV->GetInt(1,&SrcID);
			pCSV->GetText(2,&pszElType);
			pCSV->GetText(3,&pszRefSpecies);
			pCSV->GetText(4,&pszChrom);
			pCSV->GetInt(5,&StartLoci);
			pCSV->GetInt(6,&EndLoci);

			pCSV->GetText(8,&pszRelSpecies);
			if(!bFeatures)
				pCSV->GetInt(9,&Features);			// use existing features
			else
				Features = 0;
			break;

		case eCSVFprobe:							// use probe loci
			if(NumFields < 7)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expected 7+ fields in '%s', GetCurFields() returned '%d'",pszInLociFile,NumFields);
				Rslt = eBSFerrFieldCnt;
				break;
				}

			pCSV->GetInt(7,&Len);
			if(Len < MinLength || Len > MaxLength)
				continue;
			pCSV->GetInt(1,&SrcID);
			pCSV->GetText(2,&pszElType);
			pCSV->GetText(3,&pszRefSpecies);
			pCSV->GetText(4,&pszChrom);
			pCSV->GetInt(5,&StartLoci);
			pCSV->GetInt(6,&EndLoci);
			pCSV->GetText(8,&pszRelSpecies);
			Features = 0;
			break;

		case eCSVFtarget:				// use target loci
			if(NumFields < 12)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expected 12+ fields in '%s', GetCurFields() returned '%d'",pszInLociFile,NumFields);
				Rslt = eBSFerrFieldCnt;
				break;
				}

			pCSV->GetInt(7,&Len);
			if(Len < MinLength || Len > MaxLength)
				continue;
			pCSV->GetInt(1,&SrcID);
			pCSV->GetText(2,&pszElType);
			pCSV->GetText(8,&pszRefSpecies);
			pCSV->GetText(9,&pszChrom);
			pCSV->GetInt(10,&StartLoci);
			pCSV->GetInt(11,&EndLoci);
			pCSV->GetText(3,&pszRelSpecies);
			Features = 0;
			break;
	
		case eCSVFAlignM3:			// GenMAlignScore with -m3 option
			if(NumFields < 7)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expected 7+ fields in '%s', GetCurFields() returned '%d'",pszInLociFile,NumFields);
				Rslt = eBSFerrFieldCnt;
				break;
				}

			pCSV->GetInt(2,&Len);
			if(Len < MinLength || Len > MaxLength)
				continue;
			pCSV->GetInt(1,&SrcID);
			pszElType = (char *)"BlockID";
			pCSV->GetText(3,&pszRefSpecies);
			pCSV->GetText(4,&pszChrom);
			pCSV->GetInt(6,&StartLoci);
			pCSV->GetInt(7,&EndLoci);
			pszRelSpecies = pszRefSpecies;
			Features = 0;
			break;

		case eCSVFhyperX2:
			if(NumFields < 13)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expected at least 13 fields in '%s', GetCurFields() returned '%d'",pszInLociFile,NumFields);
				return(eBSFerrFieldCnt);
				}

			pCSV->GetInt(7,&Len);
			if(Len < MinLength || Len > MaxLength)
				continue;

			pCSV->GetInt(11,&OGMatchCnt);
			Identity = (int)(((100.0 * OGMatchCnt) / (double)Len) + 0.001);
			if(Identity < MinIdentity || Identity > MaxIdentity)
				continue;

			pCSV->GetInt(1,&SrcID);
			pCSV->GetText(2,&pszElType);
			pCSV->GetText(3,&pszRefSpecies);
			pCSV->GetText(4,&pszChrom);
			pCSV->GetInt(5,&StartLoci);
			pCSV->GetInt(6,&EndLoci);

			pCSV->GetText(8,&pszRelSpecies);
			if(!bFeatures)
				pCSV->GetInt(9,&Features);			// use existing features
			else
				Features = 0;

			pCSV->GetInt(10,&OGUnalignedBases);
			pCSV->GetInt(12,&OGMismatchCnt);
			pCSV->GetInt(13,&OGInDelCnt);
			break;

		}

	if(pBiobed != NULL)
		{
		if((Rslt= ChromID = pBiobed->LocateChromIDbyName(pszChrom)) > 0)
			{
			// see if overlapping feature
			FeatID=pBiobed->LocateFeatureIDinRangeOnChrom(ChromID,	  // feature is on which chromsome
									 StartLoci,       // feature must end on or after Start
									 EndLoci,		  // and start on or before End 
									 1);		  // Ith instance to return (1..n)
			if(FeatID <= 0)
				{
				NxtFeatID = pBiobed->LocateFeatureAfter(ChromID,// feature is on this chromosome
							 EndLoci);					         // feature starts on or immediately after this offset
				if(NxtFeatID > 0)
					pBiobed->GetFeature(NxtFeatID,		// feature instance identifier
							 NULL,			// where to return feature name
							 NULL,					// where to return chromosome name
							 &NxtFeatStart,		// where to return feature start on chromosome (0..n) 
							 &NxtFeatEnd);		// where to return feature end on chromosome

				PrvFeatID = pBiobed->LocateFeatureBefore(ChromID,	// feature is on this chromosome
							 StartLoci);			// feature ends on or immediately before this offset
				if(PrvFeatID > 0)
					pBiobed->GetFeature(PrvFeatID,		// feature instance identifier
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
						if((NxtFeatStart - EndLoci) < (StartLoci - PrvFeatEnd))
							FeatID = NxtFeatID;
						else
							FeatID = PrvFeatID;
						}
					}
				}

			if(bFeatures)					// override any existing feature detail?
				{
				Features = pBiobed->GetFeatureBits(ChromID,	  // feature is on which chromsome
									 StartLoci,       // feature must end on or after Start
									 EndLoci,		  // and start on or before End
									cRegionFeatBits,
									RegRegionLen);
				if(Features != 0)			// if not intergenic then check if overlaying splice sites
					Features |= pBiobed->GetSpliceSiteBits(ChromID,StartLoci,EndLoci,cMinSpliceOverlap);
				}	

			if(FeatID >= 1)
				{
				pBiobed->GetFeature(FeatID,		// feature instance identifier
							 szGeneName);	// where to return feature name
				}
			else
				strcpy(szGeneName,"NONE"); 

			}
		else
			strcpy(szGeneName,"NONE"); 
		}
	else
		strcpy(szGeneName,"NONE"); 

// now get the sequence
	if((Rslt= ChromID = pBioseq->LocateEntryIDbyName(pszChrom))<=0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to locate chrom '%s' in assembly file '%s'",pszChrom,pszInSeqFile);
		break;
		}

	if((Rslt=pBioseq->GetData(ChromID,eSeqBaseType,StartLoci,pElSeqBuff,Len)) != Len)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Loading sequence failed from chrom: %s loci %d len: %d file: '%s'",pszChrom,StartLoci,Len,pszInSeqFile);
		break;
		}

	// processing mode specific
	switch(RsltsFormat) {
		case eRsltsFCSV:	// write element to disk with sequence as CSV
			BuffLen = sprintf(pszLineBuff,"%d,\"%s\",\"%s\",\"%s\",%d,%d,%d,\"%s\",%d,\"%s\",\"",
				SrcID,pszElType,pszRefSpecies,pszChrom,StartLoci,EndLoci,Len,pszRelSpecies,Features,szGeneName);
			CSeqTrans::MapSeq2Ascii(pElSeqBuff,Len,&pszLineBuff[BuffLen]);
			BuffLen+=Len;
			BuffLen+=sprintf(&pszLineBuff[BuffLen],"\"\n");
			break;
		case eRsltsFasta:
			if(bFastaStart)
				{
				BuffLen = sprintf(pszLineBuff,">concatenated\n");
				bFastaStart = false;
				NxtFastaCol = 0;
				}
			SeqOfs = 0;
			while(Len)
				{
				NumCols = Len > 70 ? 70 : Len;
				if((NumCols + NxtFastaCol) > 70)
					NumCols = 70 - NxtFastaCol;
				CSeqTrans::MapSeq2Ascii(&pElSeqBuff[SeqOfs],NumCols,&pszLineBuff[BuffLen]);
				BuffLen += NumCols;
				NxtFastaCol += NumCols;
				if(NxtFastaCol >= 70)
					{
					BuffLen += sprintf(&pszLineBuff[BuffLen],"\n");
					NxtFastaCol = 0;
					}
				Len -= NumCols;
				SeqOfs += NumCols;
				}

			break;
		case eRsltsMultifasta:
			BuffLen = sprintf(pszLineBuff,">%s|%s|%d|%d|%d|%d %s\n",
					pszRefSpecies,pszChrom,StartLoci,EndLoci,Len,Features,szGeneName);
			SeqOfs = 0;
			while(Len)
				{
				NumCols = Len > 70 ? 70 : Len;
				CSeqTrans::MapSeq2Ascii(&pElSeqBuff[SeqOfs],NumCols,&pszLineBuff[BuffLen]);
				BuffLen += NumCols;
				BuffLen += sprintf(&pszLineBuff[BuffLen],"\n");
				Len -= NumCols;
				SeqOfs += NumCols;
				}
			break;
		}
	CUtility::SafeWrite(hRsltFile,pszLineBuff,BuffLen);
	BuffLen = 0;
	NumEls++;
	}
close(hRsltFile);
delete pCSV;
if(pBiobed != NULL)
	delete pBiobed;
delete pBioseq;
delete pElSeqBuff;
delete pszLineBuff;
return(Rslt < 0 ? NumEls : Rslt);
}

