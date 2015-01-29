// genbiomultialign.cpp : Defines the entry point for the console application.
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

const int cMARawBuffSize = 0x03ffffff;
const int cMALineSize    = 0x07fffff;

CStopWatch gStopWatch;
CDiagnostics gDiagnostics;				// for writing diagnostics messages to log file
char gszProcName[_MAX_FNAME];			// process name

CAlignValidate AlignValidate;

const unsigned int cProgVer = 305;		// increment with each release



typedef struct TAG_sProcParams 
	{
	bool bIsAXT;			// true if source files are assumed to be AXT, otherwise (false) they are MAFs
	bool bSwapRefRel;		// exchange reference and relative sequences
	char *pszRefSpecies;	// reference species name
	char *pszRelSpecies;	// relative species name
	char *pRawBuffer;		// pts to raw char buffer from which MAF/AXT lines are processed
	char *pszLineBuffer;	// pts to line buffer
	char *pszAXTBuffer;		// pts to line buffer for use by AXT processing 
	CMAlignFile *pAlignFile; // bioseq multialignment file
	} tsProcParams; 

int CreateMAlignFile(char *pszChromLens,char *pszSrcDirPath,char *pszDestAlignFile,char *pszDescr,char *pszTitle,char *pszRefSpecies,char *pszRelSpecies,bool bIsAXT,bool bSwapRefRel);
int ProcessMAFline(int LineLen, tsProcParams *pProcParams);
int ProcessAXTline(int LineLen, tsProcParams *pProcParams);
char *StripWS(char *pChr);

int GenbioDataPointsFile(char *pszMAF,char *pszDataPointsFile,char *pszDescr,char *pszTitle);

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

bool bIsAXT;
bool bSwapRefRel;
int Rslt;
int iMode;

char szOutputFileSpec[_MAX_PATH];
char szInputFileSpec[_MAX_PATH];
char szChromLens[_MAX_PATH];
char szRefSpecies[80];
char szRelSpecies[80];
char szDescription[cMBSFFileDescrLen];
char szTitle[cMBSFShortFileDescrLen];


// command line args
struct arg_lit  *help    = arg_lit0("hH","help",                "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_int *ScreenLogLevel=arg_int0("S", "ScreenLogLevel",	"<int>","Level of diagnostics written to logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","FileLogLevel","<file>",	"diagnostics log file");

struct arg_lit  *IsAXT   = arg_lit0("x","axt",					"MAF (default) or AXT two species source files");
struct arg_file *ChromLens = arg_file0("c","chromlens","<file>","input file containing species.chrom lengths, required for .axt processing");
struct arg_int  *Mode = arg_int0("m","mode","<int>",			"processing mode: 0 (default) MAF or AXT -> multialign, 1 multialign -> datapoints");
struct arg_file *InFile = arg_file1("i",NULL,"<file>",			"input from .maf or axt files");
struct arg_file *OutFile = arg_file1("o",NULL,"<file>",			"output as multialign biosequence file");
struct arg_str *Descr = arg_str1("d","descr","<string>",		"full description");
struct arg_str *Title = arg_str1("t","title","<string>",		"short title");
struct arg_str *RefSpecies = arg_str1("r","ref","<string>",		"reference species ");
struct arg_str *RelSpecies = arg_str0("R","rel","<string>",		"relative species (axt only) ");
struct arg_lit  *SwapRefRel = arg_lit0("X","exchange",			"exchange ref and rel species - only applies to AXT alignments");


struct arg_end *end = arg_end(20);

void *argtable[] = {help,version,FileLogLevel,ScreenLogLevel,LogFile,SwapRefRel,IsAXT,ChromLens,Mode,InFile,OutFile,Descr,Title,RefSpecies,RelSpecies,end};

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

	bIsAXT = IsAXT->count ? true : false;
	bSwapRefRel = SwapRefRel->count ? true : false;

	strcpy(szInputFileSpec,InFile->filename[0]);
	strcpy(szOutputFileSpec,OutFile->filename[0]);
	if(bIsAXT)
		{
		if(ChromLens->count)
			{
			strncpy(szChromLens,ChromLens->filename[0],_MAX_PATH);
			szChromLens[_MAX_PATH-1] = '\0';
			}
		else
			{
			printf("\nError: Processing .axt files but no chromosome length file '-c<chromlen_file>' specified");
			exit(1);
			}
		}
	else
		szChromLens[0] = '\0';

	iMode = Mode->count ? Mode->ival[0] : 0;

	if(Title->count < 1)
		{
		printf("No short title specified with '-t' parameter\n");
		exit(1);
		}
	if(Descr->count < 1)
		{
		printf("No descriptive text specified with '-d' parameter\n");
		exit(1);
		}

	strncpy(szTitle,Title->sval[0],sizeof(szTitle));
	szTitle[cMBSFShortFileDescrLen-1] = '\0';
	strncpy(szDescription,Descr->sval[0],sizeof(szDescription));
	szDescription[cMBSFFileDescrLen-1] = '\0';


	if(RefSpecies->count)
		strcpy(szRefSpecies,RefSpecies->sval[0]);
	else
		szRefSpecies[0] = '\0';

	if(RelSpecies->count)
		strcpy(szRelSpecies,RelSpecies->sval[0]);
	else
		szRelSpecies[0] = '\0';

	if(iMode == 0)
		{
		if(szRefSpecies[0] == '\0')
			{
			printf("\nError: In Mode 0 (creating biosequence multiple alignment file) reference species must be specified as '-r RefSpecies'");
			exit(1);
			}

		if(bIsAXT && szRelSpecies[0] == '\0')
			{
			printf("\nError: In Mode 0 (creating biosequence multiple alignment file) and processing .AXT files, relative species must be specified as '-R RelSpecies'");
			exit(1);
			}
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


	gStopWatch.Start();
#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	switch(iMode) {
		case 0:	// creating bioseq multialignment file
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Version: %d.%2.2d Processing parameters:",cProgVer/100,cProgVer%100);
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Mode: 0 (creating bioseq .ALGN multiple alignment file from .AXT or MAF files)");
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Output file: '%s'",szOutputFileSpec);
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Source files: '%s'",szInputFileSpec);
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Source files are expected to be %s alignments",bIsAXT ? "AXT" : "MAF");
			if(bIsAXT)
				gDiagnostics.DiagOutMsgOnly(eDLInfo,"Chromosome lengths file: '%s'",szChromLens);
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Reference Species: %s",szRefSpecies);
			if(bIsAXT)
				gDiagnostics.DiagOutMsgOnly(eDLInfo,"Relative Species: %s",szRelSpecies);
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Exchange Ref/Rel alignments: %s",bSwapRefRel?"yes":"no");
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Title text: %s",szTitle);
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Descriptive text: %s",szDescription);
			Rslt = CreateMAlignFile(szChromLens,szInputFileSpec,szOutputFileSpec,szDescription,szTitle,szRefSpecies,szRelSpecies,bIsAXT,bSwapRefRel);
			break;

		case 1:	// creating bioseq data points file from bioseq multialignment file
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Version: %d.%2.2d Processing parameters:",cProgVer/100,cProgVer%100);
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Mode: 1 (creating bioseq data points .DPS file from bioseq .ALGN file)");
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Output file: '%s'",szOutputFileSpec);
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Source file: '%s'",szInputFileSpec);
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Title text: %s",szTitle);
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Descriptive text: %s",szDescription);
			Rslt = GenbioDataPointsFile(szInputFileSpec,szOutputFileSpec,szDescription,szTitle);
			break;

		default:
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"requested processing mode %d not supported\n",iMode);
			exit(1);
		}

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
exit(Rslt);
}


bool
IsBase(char Base)
{
if(Base == 'a' || Base == 'A' ||
	Base == 'c' || Base == 'C' ||
	Base == 'g' || Base == 'G' ||
	Base == 't' || Base == 'T')
	return(true);
return(false);
} 

char *
StripWS(char *pChr)
{
while(*pChr != '\0' && isspace(*pChr))
	pChr++;
return(pChr);
}

int
GenbioDataPointsFile(char *pszMAF,char *pszDataPointsFile,char *pszDescr,char *pszTitle)
{
int Rslt;
CMAlignFile *pAlignFile;
if((pAlignFile = new CMAlignFile)==NULL)
	return(eBSFerrObj);
Rslt = pAlignFile->MultiAlignAll2DataPts(pszMAF,pszDataPointsFile,pszDescr,pszTitle);
if(Rslt < 0)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"GenbioDataPointsFile: Errors creating '%s' from '%s'\n",pszDataPointsFile,pszMAF);
			while(pAlignFile->NumErrMsgs())
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error-> %s",pAlignFile->GetErrMsg());
	}

delete pAlignFile;
return(Rslt);
}

// ProcessThisFile
// Parse input multialignment file. Format is expected to be either MAF or AXT
// as specified by (tsProcParams *)pParams.
int
ProcessThisFile(char *pszFile,
					void *pParams)			// will be tsProcParams
{
int Rslt;
int hFile;
char *pNxtChr;
char Chr;
int iBuffLen;							// number of buffered chrs in szBuffer
int LineLen;							// current number of chars in pszLineBuffer;
tsProcParams *pProcParams = (tsProcParams *)pParams;
bool bIsAXT = pProcParams->bIsAXT;
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing File: %s",pszFile);

#ifdef _WIN32
hFile = open(pszFile, _O_RDONLY | _O_SEQUENTIAL);
#else
hFile = open(pszFile,O_RDWR);
#endif
if(hFile == -1)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessThisFile: Unable to open input file for processing - '%s' - %s", 
				pszFile,strerror(errno));
	return(eBSFerrOpnFile);
	}
#ifdef USENEWCODE
Rslt = StartAlignBlocks(pszFile,pProcParams->pAlignFile);
#else
Rslt = pProcParams->pAlignFile->StartFile(pszFile);
#endif
if(Rslt != eBSFSuccess)
	return(Rslt);

LineLen = 0;
while((iBuffLen = read(hFile,pProcParams->pRawBuffer,cMARawBuffSize)) > 0)
	{
	pNxtChr = pProcParams->pRawBuffer;
	int BuffIdx = 0;
	while(BuffIdx++ < iBuffLen)
		{
		Chr = *pNxtChr++;
		if(Chr == '\r')		// slough CR, use NL as EOL 
			continue;
		if(Chr == '\n')		// EOL
			{
			if(!LineLen)
				continue;
			pProcParams->pszLineBuffer[LineLen]=0;
			if(bIsAXT)
				Rslt=ProcessAXTline(LineLen,pProcParams);
			else
				Rslt=ProcessMAFline(LineLen,pProcParams);
			LineLen = 0;
			continue;
			}
		if(!LineLen)		// strip any leading space
			{
			if(isspace(Chr))
				continue;
			}
		pProcParams->pszLineBuffer[LineLen++] = Chr;
		}
	}
if(LineLen)
	{
	pProcParams->pszLineBuffer[LineLen]=0;
	if(bIsAXT)
		Rslt=ProcessAXTline(LineLen,pProcParams);
	else
		Rslt=ProcessMAFline(LineLen,pProcParams);
	}
pProcParams->pAlignFile->EndAlignBlock();
pProcParams->pAlignFile->EndFile();
gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessThisFile: Total blocks: %d Chroms: %d",pProcParams->pAlignFile->GetNumBlocks(),pProcParams->pAlignFile->GetNumChroms());
return(Rslt);
}

// ProcessMAFline
// Process line from assumed MAF formated file
//
//Each sequence line begins with "s" and contains the following required fields: 
//src -- Name of one of the source sequences included in the alignment. For sequences that are resident in a browser assembly, the form database.chromosome allows automatic creation of links to other assemblies. Non-browser sequences are typically referenced by the species name alone. Species names must not contain spaces: concatenate multi-word names or replace the space with an underscore. 
//start -- Start of the aligning region in the source sequence, using zero-based position coordinates. If the strand value is "-", this field defines the start relative to the reverse-complemented source sequence. 
//size -- Size of the aligning region in the source sequence. This is equal to the number of non-dash characters in the text field (see below). 
//strand -- "+" or "-". If the value is "-", the sequence aligns to the reverse-complemented source. 
//srcSize -- Size of the entire source sequence, not just the portions involved in the alignment. 
//text -- Nucleotides (or amino acids) in the alignment are represented by upper-case letters; repeats are shown in lower case. Insertions are indicated by "-". 
//
// Notes:
// MAF start coordinates are zero based (0..n) but when displayed in the UCSC browser they are one based (1..n+1)
int
ProcessMAFline(int LineLen, tsProcParams *pProcParams)
{
static long Score;						// alignment line score or pass value
static int RefChromID;					// reference chromosome
static bool bRefChromNxt;				// true if next chromosome is the first in block - assume this is the reference
char szSpecies[80];
char szChrom[cMaxDatasetSpeciesChrom];
unsigned long iStart;
unsigned long iLen;
unsigned long iAlignLen;
long ChromLen;
char cStrand;
int Psn;
int Cnt;
int Rslt = eBSFSuccess;
char *pszLine = pProcParams->pszLineBuffer;
switch(*pszLine) {
	case '#':					// parameter line or header
		if(pszLine[1] == '#')	// '##' is a header line, should have "##maf varname=value ..." format
			{
			pszLine = StripWS(pszLine+6);
			if(*pszLine != '\0')
				Rslt = pProcParams->pAlignFile->AddHeader(pszLine);
			}
		else						// '#' is a parameter line containing parameters used to run the alignment program...
			{
			pszLine = StripWS(pszLine+1);
			if(*pszLine != '\0')
				Rslt = pProcParams->pAlignFile->AddParameters(pszLine);
			}
		return(Rslt);

	case 'a':		// new alignment block starting
		pszLine = StripWS(pszLine+1);
		if(*pszLine != '\0')
			{
			if(!sscanf(pszLine,"score=%ld",&Score))	// note: only interested in integral part of score (scores are floating points)
				Score = 0;
			}
		else
			Score = 0;

		Rslt = pProcParams->pAlignFile->StartAlignBlock(Score);		// start new block, any existing will be autoclosed
		if(Rslt < eBSFSuccess)
			{
			gDiagnostics.DiagOut(eDLWarn,gszProcName,"ProcessMAFline: Unable to start alignment block - %d - %s",Rslt,CErrorCodes::ErrText((teBSFrsltCodes)Rslt));
			break;
			}
		bRefChromNxt = true;
		break;

	case 's':		// species sequence
		pszLine = StripWS(pszLine+1);

		Cnt=sscanf(pszLine,"%s %ld %ld %c %ld %n",
						szSpecies,&iStart,&iLen,&cStrand,&ChromLen,&Psn);
		if(Cnt != 5)
			{
			pszLine[25]='\0';
			gDiagnostics.DiagOut(eDLWarn,gszProcName,"ProcessMAFline: Unable to parse species sequence params: %s",pszLine);
			return(eBSFerrParse);
			}

		char *pChr = szSpecies;
		char Chr;
		while(Chr = *pChr++)
			{
			switch(Chr) {
				case '.':
				case '_':
					pChr[-1] = '\0';
					strncpy(szChrom,pChr,sizeof(szChrom));
					szChrom[cMaxDatasetSpeciesChrom-1] = '\0';
					break;
				default:
					continue;
				}
			break;
			}
		szSpecies[cMaxDatasetSpeciesChrom-1] = '\0';
		pszLine = StripWS(pszLine+Psn);
		iAlignLen = (unsigned long)strlen(pszLine);
		if((Rslt = pProcParams->pAlignFile->AddAlignSeq(szSpecies,szChrom,ChromLen,iStart,iAlignLen,cStrand,pszLine))!= eBSFSuccess)
			{
			gDiagnostics.DiagOut(eDLWarn,gszProcName,"ProcessMAFline: Unable to add sequence for %s.%s ofs: %d len: %d\n",szSpecies,szChrom,iStart,iAlignLen);
			while(pProcParams->pAlignFile->NumErrMsgs())
				gDiagnostics.DiagOut(eDLWarn,gszProcName,"ProcessMAFline: Error-> %s",pProcParams->pAlignFile->GetErrMsg());
			}
		break;
		}
return(Rslt);
}


void
RevSeq(int SeqLenInclInDels,char *pszSeq)
{
int Cnt;
char tmp;
char *pRight = pszSeq + SeqLenInclInDels - 1;
for(Cnt = 0; Cnt < SeqLenInclInDels/2; Cnt++)
	{
	tmp = *pszSeq;
	*pszSeq++ = *pRight;
	*pRight-- = tmp;
	}
}

// ReverseComplement
// Inplace reverse complement pszSeq
bool
ReverseComplement(int SeqLenInclInDels,char *pszSeq)
{
int Cnt;
RevSeq(SeqLenInclInDels,pszSeq);
for(Cnt = 0; Cnt < SeqLenInclInDels; Cnt++, pszSeq++)
	{
	switch(*pszSeq) {
		case 'a': 
			*pszSeq = 't';
			break;
		case 'A': 
			*pszSeq = 'T';
			break;
		case 'c': 
			*pszSeq = 'g';
			break;
		case 'C': 
			*pszSeq = 'G';
			break;
		case 'g': 
			*pszSeq = 'c';
			break;
		case 'G': 
			*pszSeq = 'C';
			break;
		case 't': 
			*pszSeq = 'a';
			break;
		case 'T': 
			*pszSeq = 'A';
			break;
		case 'u': 
			*pszSeq = 'a';
			break;
		case 'U': 
			*pszSeq = 'A';
			break;
		case '\0':
			return(false);
		default:
			break;
			}
	}
return(true);
}

//AXT file Structure
//Each alignment block in an axt file contains three lines: a summary line and 2 sequence lines. Blocks are separated from one another by blank lines. 

//1. Summary line 

// 0 chr19 3001012 3001075 chr11 70568380 70568443 - 3500
// The summary line contains chromosomal position and size information about the alignment.
// It consists of 9 required fields: 
// Alignment number -- The alignment numbering starts with 0 [0..numalignments-1]
// Chromosome (primary organism) 
// Alignment start (primary organism) -- The first base is numbered 1 (in MAF co-ordinates are 0 based!). 
// Alignment end (primary organism) -- The end base is included. 
// Chromosome (aligning organism) 
// Alignment start (aligning organism) 
// Alignment end (aligning organism) 
// Strand (aligning organism)
// Blastz score -- Different blastz scoring matrices are used for different organisms
//
// Note:  -- If the strand value is "-", the values of the aligning organism's start and end fields
//           are relative to the reverse-complemented coordinates of its chromosome. 
int
ProcessAXTline(int LineLen, tsProcParams *pProcParams)
{
int Rslt;
int SeqLenInclInDels;
int RefChromID;
int RelChromID;
int RefChromLen;
int RelChromLen;

static bool b1stAlign = true;
static long iAlignNum;
static long iStart1;
static long iEnd1;
static long iStart2;
static long iEnd2;
static long iLen1;
static long iLen2;
static char cStrand;
static char szChrom1[cMaxDatasetSpeciesChrom];
static char szChrom2[cMaxDatasetSpeciesChrom];
static long Score;
char *pszLine = pProcParams->pszLineBuffer;
static int BlockID = 0;

	// alignment number as a numeric starts a new block
if(isdigit(*pszLine))
	{
			// close any currently opened block
	b1stAlign = true;
	sscanf(pszLine,"%*d %s %ld %ld %s %ld %ld %c %ld", // note: only interested in integral part of score (scores are floating points)
					szChrom1,&iStart1,&iEnd1,szChrom2,&iStart2,&iEnd2,&cStrand,&Score);
	iLen1 = iEnd1 - iStart1 + 1;
	iLen2 = iEnd2 - iStart2 + 1;

//----- hack to get around the factor that there could be zillions of scafolds in some assemblies which will exceed the limit on number of chromosomes (currently 20000)
#ifdef _WIN32
	if((szChrom2[0]=='s' || szChrom2[0]=='S') && !strnicmp("SCAFFOLD",szChrom2,8)) // if chromosome name starts with "SCAFFOLD" then make all these names the same
#else
	if((szChrom2[0]=='s' || szChrom2[0]=='S') && !strncasecmp("SCAFFOLD",szChrom2,8)) // if chromosome name starts with "SCAFFOLD" then make all these names the same

#endif
		szChrom2[8] = '\0';				// by simply truncating
// ----- end of hack
	if((Rslt=pProcParams->pAlignFile->StartAlignBlock(Score)) < eBSFSuccess)	// start a new block
		return(Rslt);

	BlockID = Rslt;
	return(eBSFSuccess);
	}

	// sequence base starts alignment
Rslt = eBSFSuccess;
if(IsBase(*pszLine))
	{
	SeqLenInclInDels = (int)strlen(pszLine);
	if(b1stAlign)
		{
		strcpy(pProcParams->pszAXTBuffer,pszLine);
		b1stAlign = false;
		return(eBSFSuccess);
		}
	// not 1st alignment so pszLine pts to 2nd alignment sequence
	RefChromID = AlignValidate.GetChromID(pProcParams->pszRefSpecies,szChrom1);
	RefChromLen = AlignValidate.GetChromLen(RefChromID);
	if(RefChromLen < 0)
		RefChromLen = 0;
	RelChromID = AlignValidate.GetChromID(pProcParams->pszRelSpecies,szChrom2);
	RelChromLen = AlignValidate.GetChromLen(RelChromID);
	if(RelChromLen < 0)
		RelChromLen = 0;

	if(pProcParams->bSwapRefRel)	// need to exchange ref and rel?
		{
		if(cStrand == '-')
			{
			ReverseComplement(SeqLenInclInDels,pszLine);
			ReverseComplement(SeqLenInclInDels,pProcParams->pszAXTBuffer);
			}

		Rslt = pProcParams->pAlignFile->AddAlignSeq(pProcParams->pszRelSpecies,szChrom2,RelChromLen,iStart2-1,SeqLenInclInDels,'+',pszLine);
		Rslt = pProcParams->pAlignFile->AddAlignSeq(pProcParams->pszRefSpecies,szChrom1,RefChromLen,iStart1-1,SeqLenInclInDels,cStrand,pProcParams->pszAXTBuffer); // note that AXT have 1..n, we use 0..n locus co-ordinates

	}
	else	// not exchanging....
		{

		Rslt = pProcParams->pAlignFile->AddAlignSeq(pProcParams->pszRefSpecies,szChrom1,RefChromLen,iStart1-1,SeqLenInclInDels,'+',pProcParams->pszAXTBuffer); // note that AXT have 1..n, we use 0..n locus co-ordinates
		Rslt = pProcParams->pAlignFile->AddAlignSeq(pProcParams->pszRelSpecies,szChrom2,RelChromLen,iStart2-1,SeqLenInclInDels,cStrand,pszLine);
		
	}

	if(Rslt < 0)
		{
		gDiagnostics.DiagOut(eDLWarn,gszProcName,"Unable to add AXT Alignment!");
		while(pProcParams->pAlignFile->NumErrMsgs())
			gDiagnostics.DiagOut(eDLWarn,gszProcName,"ProcessAXTline: Error-> %s",pProcParams->pAlignFile->GetErrMsg());
		}
	}
return(Rslt);
}

// CreateMAlignFile
// Create a multialignment file from files (either axt or mfa format) in specified source directory
int
CreateMAlignFile(char *pszChromLens,char *pszSrcDirPath,char *pszDestAlignFile,
				 char *pszDescr,char *pszTitle,
				 char *pszRefSpecies,char *pszRelSpecies,bool bIsAXT,bool bSwapRefRel)
{
int Rslt;
tsProcParams ProcParams;

memset(&ProcParams,0,sizeof(tsProcParams));

ProcParams.pRawBuffer = new char[cMARawBuffSize];
if(ProcParams.pRawBuffer == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate memory (%d bytes) for raw buffer",cMARawBuffSize);
	return(eBSFerrMem);
	}

if((ProcParams.pszLineBuffer = new char[cMALineSize])==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate memory (%d bytes) for line buffer",cMALineSize);
	delete ProcParams.pRawBuffer;
	return(eBSFerrMem);
	}

if(bIsAXT)
	{
	if((ProcParams.pszAXTBuffer = new char[cMALineSize])==NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate memory (%d bytes) for AXT buffers",cMALineSize);
		delete ProcParams.pRawBuffer;
		delete ProcParams.pszLineBuffer;
		return(eBSFerrMem);
		}

	if((Rslt = AlignValidate.ProcessChroms(pszChromLens))!=eBSFSuccess)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to process file containing chromosome lengths from %s",pszChromLens);
		delete ProcParams.pRawBuffer;
		delete ProcParams.pszLineBuffer;
		if(ProcParams.pszAXTBuffer == NULL)
			delete ProcParams.pszAXTBuffer;
		return(Rslt);
		}
	}

ProcParams.bSwapRefRel = bSwapRefRel;
ProcParams.bIsAXT = bIsAXT;
ProcParams.pszRefSpecies = pszRefSpecies;
ProcParams.pszRelSpecies = pszRelSpecies;
if((ProcParams.pAlignFile = new CMAlignFile)== NULL)
	{
	delete ProcParams.pRawBuffer;
	delete ProcParams.pszLineBuffer;
	if(ProcParams.pszAXTBuffer == NULL)
		delete ProcParams.pszAXTBuffer;
	return(eBSFerrObj);
	}

if((Rslt=ProcParams.pAlignFile->Open(pszDestAlignFile,eMAPOCreate,pszRefSpecies))!=eBSFSuccess)
	{
	while(ProcParams.pAlignFile->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,ProcParams.pAlignFile->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Failed to create output file '%s'",pszDestAlignFile);
	}
if(Rslt==eBSFSuccess)
	{
	ProcParams.pAlignFile->SetDescription(pszDescr);
	ProcParams.pAlignFile->SetTitle(pszTitle);
	CSimpleGlob glob(SG_GLOB_FULLSORT);
	if (glob.Add(pszSrcDirPath) >= SG_SUCCESS)
		{
		for (int n = 0; Rslt >= eBSFSuccess &&  n < glob.FileCount(); ++n)
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"Will process in this order: %d '%s", n+1,glob.File(n));

	    for (int n = 0; Rslt >= eBSFSuccess &&  n < glob.FileCount(); ++n)
			Rslt = ProcessThisFile(glob.File(n),&ProcParams);
		}
	else
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to glob '%s",pszSrcDirPath);
		Rslt = eBSFerrOpnFile;	// treat as though unable to open file
	    }
	ProcParams.pAlignFile->Close();
	}
delete ProcParams.pAlignFile;
AlignValidate.Reset();
delete ProcParams.pRawBuffer;
delete ProcParams.pszLineBuffer;
if(ProcParams.pszAXTBuffer == NULL)
	delete ProcParams.pszAXTBuffer;
return(Rslt);
}







