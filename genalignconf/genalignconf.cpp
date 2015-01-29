// genalignconf.cpp : Defines the entry point for the console application.


#include "stdafx.h"
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#if _WIN32

#include "../libbiokanga/commhdrs.h"
#else
#include "../libbiokanga/commhdrs.h"
#endif
#include "./AlignSubSeqs.h"

const unsigned int cProgVer = 202;		// increment with each release


const int cMALineSize    = 0x03fffff; // max seq length that can be buffered in concatenated seqs

typedef enum eProcMode {
	eProcModeStandard = 0,				// default processing
	eProcModeSummary					// summary processing
} etProcMode;

int
Process(int iMode,	// processing mode
				  char *pszRefSpecies,		// reference species
				  char *pszRelSpecies,		// relative species
				  char *pszPrimeInputFile,  // alignment file to treat as being the primary
				  char *pszSecInputFile,	// alignment file to treat as being the secondary
				  char *pszRsltsFile,		// output results to this file
				  char *pszBiobedFile,		// bed file containing features
				  bool bMultipleFeatBits,	// true if overlapping featurebits allowed
				  int RegLen,				// regulatory u/d region length
				  // following only used in mode==1
				  char *pszChrom,			// only process this chromosome - "*" for all
				  int Region,				// only process blocks which are completely contained in specified region
				  int MinBlockLen,			// blocks must of at least this length
				  int MaxBlockLen,			// and no longer than this
				  int MinNumSpecies,		// blocks must have at least this number of species
				  int MaxNumSpecies,		// and no more than this
				  char *pszFilePrefix,		// file prefix 
				  char *pszTmpDir,			// temp dir to use
				  char *pszExternExeDir);	// path to extern aligners

int
ProcessAlignments(int AlignID,				// identifies which alignment is being processed
				  char *pszRefSpecies,		// reference species
				  char *pszRelSpecies,		// relative species
				  char *pszAlignmentFile);  // alignment file

int AddSubSeq(int AlignID,			// alignment identifier
			  int BlockID,			// alignment block containing subsequence
			  int RefChromID,		// reference chromosome identifier
			  char RefStrand,		// reference strand
			  int RefChromOfs,		// reference offset
			  int RelChromID,		// relative chromosome identifier
			  char RelStrand,		// relative strand
			  int RelChromOfs,		// relative offset
			  int SubSeqLen,		// subsequence length
			  int NumIdentities);	// number of identical matches in subsequence

char *TrimWhitespace(char *pszText);

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
	return _T("genalignconf");
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
int iMode;
bool bMultipleFeatBits;

char szOutputFile[_MAX_PATH];
char szPrimeInputFile[_MAX_PATH];
char szSecondaryInputFile[_MAX_PATH];

char szInputBiobedFile[_MAX_PATH];
char szRefSpecies[50];
char szRelSpecies[50];
int iRegLen;
int iRegion;
char szChrom[50];

int iMinBlockLen;								// blocks must of at least this length
int iMaxBlockLen;								// and no longer than this
int iMinNumSpecies;								// blocks must have at least this number of species
int iMaxNumSpecies;								// and no more than this
char szFilePrefix[50];							// file prefix 
char szTmpDir[_MAX_PATH];						// temp dir to use
char szExternExeDir[_MAX_PATH];					// path to extern aligners

// command line args
struct arg_lit  *help    = arg_lit0("hH","help",                "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_int *ScreenLogLevel=arg_int0("S", "ScreenLogLevel",	"<int>","Level of diagnostics written to logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",	"diagnostics log file");

struct arg_str  *RefSpecies = arg_str1("r","refspecies","<string>","reference species");
struct arg_str  *RelSpecies = arg_str1("R","relspecies","<string>","relative species");

struct arg_int  *RegLen = arg_int0("L","updnstream","<int>",	"length of 5'up or 3'down  stream regulatory region length (default = 2000) 0..1000000");
struct arg_file *PrimeInputFile = arg_file1("i",NULL,"<file>",		"input from primary (reference) .algn file");
struct arg_file *SecondaryInputFile = arg_file0("I",NULL,"<file>","input from secondary .algn file");
struct arg_file *OutFile = arg_file1("o",NULL,"<file>",			"output to statistics file as CSV");
struct arg_file *InBedFile = arg_file0("b","bed","<file>",		"characterise regions from biobed file");
struct arg_lit  *MultipleFeatBits  = arg_lit0("q","multiplefeatbits",	"single featbit (default) or multiple featbits allowed");
struct arg_lit  *ChromPer = arg_lit0("c","chromper",		"generate stats for each chromosome -= default is for complete genome");
struct arg_int  *Mode = arg_int0("m","procmode","<int>",	"processing mode 0:default, 1:extended, 2:summary");

struct arg_str  *Chrom = arg_str0("C","chrom","<string>",	"only process blocks with ref chromsome");
struct arg_int  *Region = arg_int0("k","region","<int>",	"only process blocks completely contained in this region (0:IG..6:3'DS)");
struct arg_int  *MinNumSpecies = arg_int0("z","minspecies", "<int>",	"minimum number of species in aligned block (default 2)");
struct arg_int  *MaxNumSpecies = arg_int0("Z","maxspecies","<int>",	"maximum number of species in aligned block (default 50)");
struct arg_int  *MinBlockLen = arg_int0("x","minblocklen","<int>",	"minimum block length (default = 100)");
struct arg_int  *MaxBlockLen = arg_int0("X","maxblocklen","<int>",	"maximum block length (default = 100000)");

struct arg_str *FilePrefix = arg_str0("p","prefix","<string>","file prefix to use");
struct arg_file *TmpDir = arg_file0("d","intermed","<file>","write intermediate files in this dir");
struct arg_file *ExternExeDir = arg_file0("D","externalign","<file>","path to extern aligners - clustalw.exe,muscle.exe");

struct arg_end *end = arg_end(20);

void *argtable[] = {help,version,FileLogLevel,ScreenLogLevel,LogFile,
					Mode,
					RegLen,PrimeInputFile,SecondaryInputFile,OutFile,InBedFile,
					MultipleFeatBits,ChromPer,
					RefSpecies,RelSpecies,
					Chrom,Region,MinNumSpecies,MaxNumSpecies,MinBlockLen,MaxBlockLen,
					FilePrefix,TmpDir,ExternExeDir,
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
	if(iMode < eProcModeStandard || iMode > eProcModeSummary)
		{
		printf("\nError: Requested processing mode '-m%d' not supported",iMode);
		exit(1);
		}

	bMultipleFeatBits = MultipleFeatBits->count ? true : false;

	strncpy(szRefSpecies,RefSpecies->sval[0],sizeof(szRefSpecies));
	szRefSpecies[sizeof(szRefSpecies)-1] = '\0';

	if(iMode == 0)
		{
		if(!RelSpecies->count)
			{
			printf("\nError: Relative species must be specified with '-R<species>' if '-m%d' mode requested",iMode);
			exit(1);
			}
		strncpy(szRelSpecies,RelSpecies->sval[0],sizeof(szRelSpecies));
		szRelSpecies[sizeof(szRelSpecies)-1] = '\0';
		if(!SecondaryInputFile->count)
			{
			printf("\nError: Secondary input file must be specified with '-I<file>' if '-m%d' mode requested",iMode);
			exit(1);
			}

		strcpy(szSecondaryInputFile,SecondaryInputFile->filename[0]);
		szFilePrefix[0] = '\0';
		szTmpDir[0] = '\0';
		szExternExeDir[0] = '\0';
		iRegion = 0;
		}
	else
		{
		szRelSpecies[0] = '\0';
		szSecondaryInputFile[0] = '\0';

		if(Chrom->count)
			{
			strncpy(szChrom,Chrom->sval[0],sizeof(szChrom));
			szChrom[sizeof(szChrom)-1] = '\0';
			}
		else
			strcpy(szChrom,"*");

		if(!Region->count)
			{
			if(iMode > 1)
				{
				printf("\nError: Containing region must be specified with '-k<region>' if '-m%d' mode requested",iMode);
				exit(1);
				}
			else
				iRegion = -1;
			}
		else
			{
			iRegion = Region->ival[0];
			if(iRegion < 0 || iRegion > 6)
				{
				printf("\nError: Region specified as '-k%d' must be in range 0..6",iRegion);
				exit(1);
				}
			}

		if(!TmpDir->count)
			{
			printf("\nError: Intermediate directory must be specified with '-d<file>' if '-m%d' mode requested",iMode);
			exit(1);
			}
		strcpy(szTmpDir,TmpDir->filename[0]);

		if(!ExternExeDir->count)
			{
			printf("\nError: Path to directory containing external aligners must be specified with '-D<file>' if '-m%d' mode requested",iMode);
			exit(1);
			}
		strcpy(szExternExeDir,ExternExeDir->filename[0]);

		if(!FilePrefix->count)
			{
			printf("\nError: Filename prefix must be specified with '-p<prefix>' if '-m%d' mode requested",iMode);
			exit(1);
			}
		strcpy(szFilePrefix,FilePrefix->sval[0]);

		if(MinNumSpecies->count)
			{
			iMinNumSpecies = MinNumSpecies->ival[0];
			if(iMinNumSpecies < 2 || iMinNumSpecies > 50)
				{
				printf("\nError: Minimum number species '-z%d' must be in range 2-50",iMinNumSpecies);
				exit(1);
				}
			}
		else
			iMinNumSpecies = 2;

		if(MaxNumSpecies->count)
			{
			iMaxNumSpecies = MaxNumSpecies->ival[0];
			if(iMaxNumSpecies < iMinNumSpecies || iMaxNumSpecies > 50)
				{
				printf("\nError: Minimum number species '-z%d' must be in range %d-50",iMaxNumSpecies,iMinNumSpecies);
				exit(1);
				}
			}
		else
			iMaxNumSpecies = 50;

		if(MinBlockLen->count)
			{
			iMinBlockLen = MinBlockLen->ival[0];
			if(iMinBlockLen < 10 || iMinBlockLen > 100000)
				{
				printf("\nError: Minimum block length '-x%d' must be in range 10-100000",iMinBlockLen);
				exit(1);
				}
			}
		else
			iMinBlockLen = 100;


		if(MaxBlockLen->count)
			{
			iMaxBlockLen = MaxBlockLen->ival[0];
			if(iMaxBlockLen < iMinBlockLen || iMaxBlockLen > 100000)
				{
				printf("\nError: Maximum block length '-X%d' must be in range %d-100000",iMinBlockLen,iMaxBlockLen);
				exit(1);
				}
			}
		else
			iMaxBlockLen = 100000;
		}

	strcpy(szPrimeInputFile,PrimeInputFile->filename[0]);
	strcpy(szOutputFile,OutFile->filename[0]);

	if(InBedFile->count)
		strcpy(szInputBiobedFile,InBedFile->filename[0]);
	else
		szInputBiobedFile[0] = '\0';

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
	

			// now that command parameters have been parsed then initialise diagnostics log system
	if(!gDiagnostics.Open(szLogFile,(etDiagLevel)iScreenLogLevel,(etDiagLevel)iFileLogLevel,true))
		{
		printf("\nError: Unable to start diagnostics subsystem.");
		if(szLogFile[0] != '\0')
			printf(" Most likely cause is that logfile '%s' can't be opened/created",szLogFile);
		exit(1);
		}

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Version: %d.%2.2d Processing parameters:",cProgVer/100,cProgVer%100);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Processing mode: %s",iMode == eProcModeStandard ? "standard" : iMode == eProcModeStandard ? "extended" : "summary");

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Reference species: '%s'",szRefSpecies);
	if(iMode == 0)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Relative species: '%s'",szRelSpecies);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Primary multialignment (.algn) file to process: '%s'",szPrimeInputFile);
	if(iMode == 0)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Secondary multialignment (.algn) file to process: '%s'",szSecondaryInputFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"where to write out stats: '%s'",szOutputFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"biobed file containing regional features: '%s'",szInputBiobedFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Multiple featurebits allowed: %s",bMultipleFeatBits ? "yes" : "no");
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Regulatory region size: %d",iRegLen);
	if(iMode == 1)
		{
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"reference species chromsome : %s",szChrom);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"minimum number of species in aligned block : %d",iMinNumSpecies);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"maximum number of species in aligned block : %d",iMaxNumSpecies);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"minimum block length : %d",iMinBlockLen);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"maximum block length: %d",iMaxBlockLen);

		gDiagnostics.DiagOutMsgOnly(eDLInfo,"File prefix to use: '%s'",szFilePrefix);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Write intermediate files in this dir: '%s'",szTmpDir);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Path to extern aligners: '%s'",szExternExeDir);
		}
	gStopWatch.Start();

#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif

	Rslt = Process(iMode,szRefSpecies,szRelSpecies,szPrimeInputFile,szSecondaryInputFile,szOutputFile,
		szInputBiobedFile,bMultipleFeatBits,iRegLen,
		szChrom,iRegion,
		iMinBlockLen,iMaxBlockLen,iMinNumSpecies,iMaxNumSpecies,
		szFilePrefix,szTmpDir,szExternExeDir);
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


int
Process(int Mode,	// processing mode
				  char *pszRefSpecies,		// reference species
				  char *pszRelSpecies,		// relative species
				  char *pszPrimeInputFile,  // alignment file to treat as being the primary
				  char *pszSecInputFile,	// alignment file to treat as being the secondary
				  char *pszRsltsFile,		// output results to this file
				  char *pszBiobedFile,		// bed file containing features
				  bool bMultipleFeatBits,	// true if overlapping featurebits allowed
				  int RegLen,				// regulatory u/d region length
				  				  // following only used in mode==1
  				  char *pszChrom,			// only process this chromosome - "*" for all
				  int Region,				// only process blocks which are completely contained in specified region
				  int MinBlockLen,			// blocks must of at least this length
				  int MaxBlockLen,			// and no longer than this
				  int MinNumSpecies,		// blocks must have at least this number of species
				  int MaxNumSpecies,		// and no more than this
				  char *pszFilePrefix,		// file prefix 
				  char *pszTmpDir,			// temp dir to use
				  char *pszExternExeDir)	// path to extern aligners
{
int Rslt;
CAlignSubSeqs *pAlignSubSeqs;

pAlignSubSeqs = new CAlignSubSeqs();

if((Rslt=pAlignSubSeqs->InitResults(Mode,pszRsltsFile,pszBiobedFile,bMultipleFeatBits,RegLen))!= eBSFSuccess)
	return(Rslt);

switch(Mode) {
	case 0:
		// process all subsequences in the secondary alignment into sorted array of subseqs 
		Rslt = pAlignSubSeqs->ProcessAlignment(true,pszRefSpecies,pszRelSpecies,pszSecInputFile);
		if(Rslt >= eBSFSuccess)		// any errors?
			{
			// now iterate over all subsequences in the primary alignment, looking for reciprocal
			// alignments in the previously processed secondary alignment 
			Rslt = pAlignSubSeqs->ProcessAlignment(true,pszRefSpecies,pszRelSpecies,pszPrimeInputFile);
			}
		pAlignSubSeqs->EndResults();
		break;
	case 1:
		Rslt = pAlignSubSeqs->ProcessCompareAlignment(
							  pszChrom,			// only process this chromosome - "*" for all
								Region, // blocks must be completely contained within this region
								   MinBlockLen,	// blocks must of at least this length
									   MaxBlockLen,			// and no longer than this
									   MinNumSpecies,		// blocks must have at least this number of species
									   MaxNumSpecies,		// and no more than this
										pszFilePrefix,	// file prefix 
										pszTmpDir, // temp dir to use
									    pszExternExeDir, // path to extern aligners
									   pszRefSpecies,// reference species
					   				   pszPrimeInputFile);	// alignment file
		break;
	}

if(pAlignSubSeqs != NULL)
	delete pAlignSubSeqs;

return(Rslt);
}

