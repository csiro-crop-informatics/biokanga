// CSVMerge.cpp : Defines the entry point for the console application.
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

#include "./ChromMaps.h"

const char *cpszProgVer = "1.1.1";		// increment with each release

const int cDfltMinLength = 4;			// default minimum element length
const int cDfltMaxLength = 1000000;		// default maximum element length
const int cMaxLengthRange = 10000000;	// maximal element length
const int cDfltJoinOverlap = 0;			// join output elements which are separated by at most this many bases
const int cMaxJoinOverlap  = 1000000;	// max allowed join overlap
const int cMaxExtendLength = 1000000;	// max allowed flank extension


const int cRsltsLineLen = 4096;			// buffer size for results loci lines

// processing mode
typedef enum eProcMode {
	ePMElIntersect = 0,				// Ref and Rel intersect (ref,rel must both be present)
	ePMElRefExclusive,				// Ref exclusive	(ref only)
	ePMElRelExclusive,				// Rel exclusive	(rel only)
	ePMElRefRelUnion,				// Ref and Rel union  (ref or rel)
	ePMElRefNotRefRel				// Ref or Rel not present
} etProcMode;



typedef struct TAG_sProcParams 
	{
	etProcMode ProcMode;	// processing mode
	int MinLength;			// only process input ref/rel elements after any flanking of at least this length
	int MaxLength;			// only process input ref/rel elements after any flanking which are no longer than this length
	char *pszRefFile;		// CSV file containing ref elements
	char *pszRelFile;		// CSV file containing rel elements
	char *pszOutLociFile;	// write element loci into this CSV file
	char *pszRefSpecies;	// replace element's ref species with this in OutLociFile
	char *pszRelSpecies;	// replace element's rel species with this in OutLociFile
	char *pszElType;		// replace element's type with this in OutLociFile

	int hRsltsLociFile;		// file handle for results loci file
	int NumRsltLines;		// number of result lines generated 
	int RsltLineBuffOfs;	// offset in szRsltLineBuff at which to write next results line
	char szRsltLineBuff[cRsltsLineLen]; // used to buffer output results loci lines

	int RefExtend;			// extend ref left+right flanks by this many bases
	int RelExtend;			// extend rel left+right flanks by this many bases

	// when outputing the merged elements then the following applies
	int JoinDistance;		// if > 0 then join elements which only differ by at most this distance beween end of element i and start of element i+1
	int MinMergeLength;		// only output merged elements of at least this length
	int MaxMergeLength;		// only output merged elements which are no longer than this length

	char szCurChrom[cMaxDatasetSpeciesChrom];	// used to hold chromosome when joining adjacent elements
	int CurStartLoci;		// used to hold element start when joining adjacent elements
	int CurEndLoci;			// used to hold element end when joining adjacent elements

	int NumUnmergedEls;		// number of elements after processing but not yet merged
	int NumMergedEls;		// number of merged elements before filtering for length
	int NumMergedFiltUnderLen;	// number of merged elements filtered out because of underlength
	int NumMergedFiltOverLen;	// number of merged elements filtered out because of overlength
	int NumMergedAccepted;  // number of merged elements accepted and output 
} tsProcParams; 

bool CleanupResources(tsProcParams *pProcParams);

int 
Process(etProcMode ProcMode,	// processing mode
		int MinLength,			// only process ref/rel elements of at least this length
		int MaxLength,			// only process ref/rel elements which are no longer than this length
		int RefExtend,			// extend ref left+right flanks by this many bases
		int RelExtend,			// extend rel left+right flanks by this many bases
		int JoinDistance,		// if > 0 then join elements which only differ by at most this distance beween end of element i and start of element i+1
		int MinMergeLength,		// only output merged elements of at least this length
		int MaxMergeLength,		// only output merged elements which are no longer than this length
		char *pszRefFile,		// CSV file containing ref elements
		char *pszRelFile,		// CSV file containing rel elements
		char *pszOutLociFile,	// write element loci into this CSV file
		char *pszRefSpecies,	// replacement ref species in output loci file
		char *pszRelSpecies,	// replacement rel species in output loci file
		char *pszElType);		// replacement element type in output loci file

char *ProcMode2Txt(etProcMode ProcMode);
int OutputLoci(char *pszChrom,int StartLoci,int EndLoci,void *pParm);
int WriteLoci(char *pszChrom,int StartLoci,int EndLoci,tsProcParams *pParams);

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
	return _T("CSVMerge");
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
int iProcMode;
int iMinLength;
int iMaxLength;

int iMinMergeLength;
int iMaxMergeLength;

char szRefFile[_MAX_PATH];	// process ref hypers from this file
char szRelFile[_MAX_PATH];	// process rel hypers from this file
char szOutLociFile[_MAX_PATH];	// write loci to this file

char szRefSpecies[cMaxDatasetSpeciesChrom];	// use this species as the ref species in generated szOutLociFile
char szRelSpecies[cMaxDatasetSpeciesChrom];	// use this species/list as the rel species in generated szOutLociFile
char szElType[cMaxDatasetSpeciesChrom];		// use this as the element type in generated szOutLociFile

int iRefExtend;			// extend ref element lengths left+right by this many bases
int iRelExtend;			// extend rel element lengths left+right by this many bases
int iJoinDistance;		// if > 0 then join elements which only differ by at most this distance beween end of element i and start of element i+1


// command line args
struct arg_lit  *help    = arg_lit0("h","help",                 "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to screen and logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"diagnostics log file");

struct arg_file *RefFile = arg_file1("i","reffile","<file>",	"reference hyper element CSV file");
struct arg_file *RelFile = arg_file0("I","relfile","<file>",	"relative hyper element CSV file");
struct arg_file *OutLociFile = arg_file1("o",NULL,"<file>",		"output loci to file as CSV");

struct arg_str  *RefSpecies = arg_str1("r","refspecies","<string>","output loci file ref species");
struct arg_str  *RelSpecies = arg_str1("R","relspecies","<string>","output loci file rel species");
struct arg_str  *ElType = arg_str0("t","eltype","<string>","output loci file element type");

struct arg_int  *ProcMode = arg_int0("p","mode","<int>",		 "processing mode: 0:Intersect (Ref & Rel)\n\t\t1:Ref exclusive (Ref & !Rel)\n\t\t2:Rel exclusive (!Ref & Rel)\n\t\t3:Union (Ref | Rel)\n\t\t4:Neither (!(Ref | Rel))");
struct arg_int  *MinLength = arg_int0("l","minlength","<int>",   "minimum input ref/rel element length (default 4)");
struct arg_int  *MaxLength = arg_int0("L","maxlength","<int>",   "maximum input ref/rel element length (default 1000000)");

struct arg_int  *MinMergeLength = arg_int0("m","minmergelength","<int>","minimum merged output element length (default 4)");
struct arg_int  *MaxMergeLength = arg_int0("M","maxmergelength","<int>","maximum merged output element length (default 1000000)");

struct arg_int  *RefExtend = arg_int0("e","refextend","<int>",	 "extend ref element flanks left+right by this many bases (default 0)");
struct arg_int  *RelExtend = arg_int0("E","relextend","<int>",	 "extend rel element flanks left+right by this many bases (default 0)");

struct arg_int  *JoinDistance = arg_int0("j","join","<int>",     "merge output elements which are only separated by this number of bases (default 0)");

struct arg_end *end = arg_end(20);

void *argtable[] = {help,version,FileLogLevel,LogFile,
					ProcMode,
					RefFile,RelFile,OutLociFile,
					MinLength,MaxLength,RefExtend,RelExtend,JoinDistance,MinMergeLength,MaxMergeLength,
					RefSpecies,RelSpecies,ElType,
					end};

char **pAllArgs;
int argerrors;
argerrors = CUtility::arg_parsefromfile(argc,(char **)argv,&pAllArgs);
if(argerrors >= 0)
	argerrors = arg_parse(argerrors,pAllArgs,argtable);

    /* special case: '--help' takes precedence over error reporting */
if (help->count > 0)
        {
		printf("\n%s CSV Merge Elements, Version %s\nOptions ---\n", gszProcName,cpszProgVer);
        arg_print_syntax(stdout,argtable,"\n");
        arg_print_glossary(stdout,argtable,"  %-25s %s\n");
		printf("\nNote: Parameters can be entered into a parameter file, one parameter per line.");
		printf("\n      To invoke this parameter file then precede it's name with '@'");
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


	iProcMode = ProcMode->count ? ProcMode->ival[0] : ePMElIntersect;
	if(iProcMode < ePMElIntersect || iProcMode > ePMElRefNotRefRel)
		{
		printf("Error: Processing mode '-p%d' is not in range 0..4",iProcMode);
		exit(1);
		}

	strncpy(szOutLociFile,OutLociFile->filename[0],_MAX_PATH);
	szOutLociFile[_MAX_PATH-1] = '\0';

	strncpy(szRefSpecies,RefSpecies->sval[0],sizeof(szRefSpecies));
	szRefSpecies[sizeof(szRefSpecies)-1] = '\0';
	strncpy(szRelSpecies,RelSpecies->sval[0],sizeof(szRelSpecies));
	szRelSpecies[sizeof(szRelSpecies)-1] = '\0';
	if(ElType->count)
		{
		strncpy(szElType,ElType->sval[0],sizeof(szElType));
		szElType[sizeof(szElType)-1] = '\0';
		}
	else
		strcpy(szElType,"merged");

	iMinLength = MinLength->count ? MinLength->ival[0] : cDfltMinLength;
	if(iMinLength < 0 || iMinLength > cMaxLengthRange)
		{
		printf("Error: Minimum element length '-l%d' is not in range 0..%d",iMinLength,cMaxLengthRange);
		exit(1);
		}

	iMaxLength = MaxLength->count ? MaxLength->ival[0] : cDfltMaxLength;
	if(iMaxLength < iMinLength || iMaxLength > cMaxLengthRange)
		{
		printf("Error: Maximum element length '-L%d' is not in range %d..%d",iMaxLength,iMinLength,cMaxLengthRange);
		exit(1);
		}

	iMinMergeLength = MinMergeLength->count ? MinMergeLength->ival[0] : cDfltMinLength;
	if(iMinMergeLength < 0 || iMinMergeLength > cMaxLengthRange)
		{
		printf("Error: Minimum output merged element length '-m%d' is not in range 0..%d",iMinMergeLength,cMaxLengthRange);
		exit(1);
		}

	iMaxMergeLength = MaxMergeLength->count ? MaxMergeLength->ival[0] : cDfltMaxLength;
	if(iMaxMergeLength < iMinMergeLength || iMaxMergeLength > cMaxLengthRange)
		{
		printf("Error: Maximum element length '-M%d' is not in range %d..%d",iMaxMergeLength,iMinMergeLength,cMaxLengthRange);
		exit(1);
		}

	iJoinDistance = JoinDistance->count ? JoinDistance->ival[0] : cDfltJoinOverlap;
	if(iJoinDistance < 0 || iJoinDistance > cMaxJoinOverlap)
		{
		printf("Error: Join separation length '-j%d' is not in range %d..%d",iJoinDistance,0,cMaxJoinOverlap);
		exit(1);
		}

	iRefExtend = RefExtend->count ? RefExtend->ival[0] : 0;
	if(iRefExtend < (-1 * cMaxExtendLength) || iRefExtend > cMaxExtendLength)
		{
		printf("Error: Ref Extension length '-e%d' is not in range %d..%d",iRefExtend,(-1 * cMaxExtendLength),cMaxExtendLength);
		exit(1);
		}

	iRelExtend = RelExtend->count ? RelExtend->ival[0] : 0;
	if(iRelExtend < (-1 * cMaxExtendLength) || iRelExtend > cMaxExtendLength)
		{
		printf("Error: Rel Extension length '-E%d' is not in range %d..%d",iRelExtend,(-1 * cMaxExtendLength),cMaxExtendLength);
		exit(1);
		}


	strncpy(szRefFile,RefFile->filename[0],_MAX_PATH);
	szRefFile[_MAX_PATH-1] = '\0';

	if(RelFile->count)
		{
		strncpy(szRelFile,RelFile->filename[0],_MAX_PATH);
		szRelFile[_MAX_PATH-1] = '\0';
		}
	else
		{
		if(iProcMode == ePMElRefExclusive || iProcMode == ePMElRefRelUnion)
			szRelFile[0] = '\0';
		else
			{
			printf("Error: Rel loci file must be specified in processing mode '-p%d' (%s)",iProcMode,ProcMode2Txt((etProcMode)iProcMode));
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
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Version: %s Processing parameters:",cpszProgVer);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Processing Mode: %d (%s)",iProcMode,ProcMode2Txt((etProcMode)iProcMode));
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Reference CSV file: '%s'",szRefFile);
	if(szRelFile[0] != '\0')
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Relative CSV file: '%s'",szRelFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Output processed loci into CSV file: '%s'",szOutLociFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Output loci file ref species: '%s'",szRefSpecies);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Output loci file rel species: '%s'",szRelSpecies);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Output loci file element type: '%s'",szElType);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Minimum input element length: %d",iMinLength);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Maximum input element length: %d",iMaxLength);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Ref element flank extension length: %d",iRefExtend);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Rel element flank extension length: %d",iRelExtend);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Merge output elements separated by at most this many bases: %d",iJoinDistance);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Minimum output merged element length: %d",iMinMergeLength);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Maximum output merged element length: %d",iMaxMergeLength);


	// processing here...
	gStopWatch.Start();
#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif	
	Rslt = Process((etProcMode)iProcMode,iMinLength,iMaxLength,iRefExtend,iRelExtend,iJoinDistance,iMinMergeLength,iMaxMergeLength,szRefFile,szRelFile,szOutLociFile,
		szRefSpecies,szRelSpecies,szElType);
	gStopWatch.Stop();
	Rslt = Rslt >=0 ? 0 : 1;
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Exit code: %d Total processing time: %s",Rslt,gStopWatch.Read());
	exit(Rslt);
	}
else
	{
	printf("\n%s CSV Merge Elements, Version %s\n",gszProcName,cpszProgVer);
	arg_print_errors(stdout,end,gszProcName);
	arg_print_syntax(stdout,argtable,"\nUse '-h' to view option and parameter usage\n");
	exit(1);
	}
}


char *
ProcMode2Txt(etProcMode ProcMode)
{
switch(ProcMode) {
	case ePMElIntersect:
		return((char *)"Ref and Rel intersect (ref,rel must both be present)");
	case ePMElRefExclusive:				
		return((char *)"Ref exclusive (ref only)");
	case ePMElRelExclusive:
		return((char *)"Rel exclusive (rel only)");
	case ePMElRefRelUnion:
		return((char *)"Ref and Rel union (ref or rel)");
	case ePMElRefNotRefRel:
		return((char *)"Neither Ref or Rel present");
	default:
		break;
	}
return((char *)"Unrecognised");
}

int 
Process(etProcMode ProcMode,	// processing mode
		int MinLength,			// only process ref/rel elements of at least this length
		int MaxLength,			// only process ref/rel elements which are no longer than this length
		int RefExtend,			// extend ref left+right flanks by this many bases
		int RelExtend,			// extend rel left+right flanks by this many bases
		int JoinDistance,		// if > 0 then join elements which only differ by at most this distance beween end of element i and start of element i+1
		int MinMergeLength,		// only output merged elements of at least this length
		int MaxMergeLength,		// only output merged elements which are no longer than this length
		char *pszRefFile,		// CSV file containing ref elements
		char *pszRelFile,		// CSV file containing rel elements
		char *pszOutLociFile,	// write element loci into this CSV file
		char *pszRefSpecies,	// replacement ref species in output loci file
		char *pszRelSpecies,	// replacement rel species in output loci file
		char *pszElType)		// replacement element type in output loci file
{
int Rslt;
CChromMaps *pChromMaps;

tsProcParams ProcParams;

memset(&ProcParams,0,sizeof(ProcParams));

ProcParams.MinLength = MinLength;
ProcParams.MaxLength = MaxLength;
ProcParams.MinMergeLength = MinMergeLength;
ProcParams.MaxMergeLength = MaxMergeLength;
ProcParams.pszOutLociFile = pszOutLociFile;
ProcParams.pszRefFile=pszRefFile;
ProcParams.pszRelFile=pszRelFile;
ProcParams.pszRefSpecies=pszRefSpecies;
ProcParams.pszRelSpecies=pszRelSpecies;
ProcParams.pszElType=pszElType;
ProcParams.RefExtend = RefExtend;
ProcParams.RelExtend = RelExtend;
ProcParams.JoinDistance = JoinDistance;

pChromMaps = new CChromMaps;
if((Rslt = pChromMaps->Load(false,MinLength,MaxLength,RefExtend,pszRefFile,RelExtend,pszRelFile))<0)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to load files");
	return(Rslt);
	}

#ifdef _WIN32
if((ProcParams.hRsltsLociFile = open(pszOutLociFile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
if((ProcParams.hRsltsLociFile = open(pszOutLociFile, O_RDWR | O_CREAT | O_TRUNC, S_IREAD | S_IWRITE ))==-1)
#endif
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create loci file: %s - %s",pszOutLociFile,strerror(errno));
	CleanupResources(&ProcParams);
	return(eBSFerrCreateFile);
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Loci CSV file created/truncated: '%s'",pszOutLociFile);


switch(ProcMode) {
	case ePMElIntersect:					
		Rslt = pChromMaps->Process(eElIntersect,OutputLoci,&ProcParams);
		break;
	case ePMElRefExclusive:				
		Rslt = pChromMaps->Process(eElRefExclusive,OutputLoci,&ProcParams);
		break;

	case ePMElRelExclusive:				
		Rslt = pChromMaps->Process(eElRelExclusive,OutputLoci,&ProcParams);
		break;

	case ePMElRefRelUnion:				
		Rslt = pChromMaps->Process(eElRefRelUnion,OutputLoci,&ProcParams);
		break;

	case ePMElRefNotRefRel:				
		Rslt = pChromMaps->Process(eElRefNotRefRel,OutputLoci,&ProcParams);
		break;
	}

if(Rslt >= 0 && ProcParams.hRsltsLociFile != -1 && ProcParams.JoinDistance > 0 && ProcParams.szCurChrom[0] != '\0')
	Rslt = WriteLoci(ProcParams.szCurChrom,ProcParams.CurStartLoci,ProcParams.CurEndLoci,&ProcParams);

if(Rslt >= 0 && ProcParams.hRsltsLociFile != -1 && ProcParams.RsltLineBuffOfs)
	{
	CUtility::SafeWrite(ProcParams.hRsltsLociFile,ProcParams.szRsltLineBuff,ProcParams.RsltLineBuffOfs);
	ProcParams.RsltLineBuffOfs = 0;
	}

if(Rslt >= 0 && ProcParams.hRsltsLociFile != -1)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Els Output: %d PreMerge: %d PostMerge: %d FilterUnderLen: %d FilterOverLen: %d",
		ProcParams.NumMergedAccepted,ProcParams.NumUnmergedEls,ProcParams.NumMergedEls,ProcParams.NumMergedFiltUnderLen,ProcParams.NumMergedFiltOverLen);
	}
CleanupResources(&ProcParams);
return(Rslt > 0 ? 0 : Rslt);
}


bool 
CleanupResources(tsProcParams *pProcParams)
{
if(pProcParams->hRsltsLociFile != -1)
	{
	close(pProcParams->hRsltsLociFile);
	pProcParams->hRsltsLociFile = -1;
	}
return(true);
}

int
OutputLoci(char *pszChrom,int StartLoci,int EndLoci,void *pParm)
{
int Rslt;

tsProcParams *pParams = (tsProcParams *)pParm;

if(pParams==NULL || pParams->hRsltsLociFile == -1)
	return(eBSFerrParams);

pParams->NumUnmergedEls += 1;

if(!pParams->JoinDistance)
	return(WriteLoci(pszChrom,StartLoci,EndLoci,pParams));

if(pParams->szCurChrom[0] == '\0')
	{
	strncpy(pParams->szCurChrom,pszChrom,sizeof(pParams->szCurChrom));
	pParams->szCurChrom[sizeof(pParams->szCurChrom)-1] = '\0';
	pParams->CurStartLoci = StartLoci;
	pParams->CurEndLoci = EndLoci;
	return(eBSFSuccess);
	}

if(stricmp(pszChrom,pParams->szCurChrom))
	{
	Rslt = WriteLoci(pParams->szCurChrom,pParams->CurStartLoci,pParams->CurEndLoci,pParams);
	strcpy(pParams->szCurChrom,pszChrom);
	pParams->CurStartLoci = StartLoci;
	pParams->CurEndLoci = EndLoci;
	return(eBSFSuccess);
	}

if((StartLoci - pParams->CurEndLoci - 1) > pParams->JoinDistance)
	{
	Rslt = WriteLoci(pParams->szCurChrom,pParams->CurStartLoci,pParams->CurEndLoci,pParams);
	pParams->CurStartLoci = StartLoci;
	pParams->CurEndLoci = EndLoci;
	return(eBSFSuccess);
	}

pParams->CurEndLoci = EndLoci;
return(eBSFSuccess);
}

int
WriteLoci(char *pszChrom,int StartLoci,int EndLoci,tsProcParams *pParams)
{
int Len = EndLoci - StartLoci + 1;

pParams->NumMergedEls += 1;
if(Len < pParams->MinMergeLength)
	{
	pParams->NumMergedFiltUnderLen += 1;
	return(eBSFSuccess);
	}
if(Len > pParams->MaxMergeLength)
	{
	pParams->NumMergedFiltOverLen += 1;
	return(eBSFSuccess);
	}
pParams->NumMergedAccepted += 1;
pParams->RsltLineBuffOfs += sprintf(&pParams->szRsltLineBuff[pParams->RsltLineBuffOfs],"%d,\"%s\",\"%s\",\"%s\",%d,%d,%d,\"%s\",0\n",
				(pParams->NumRsltLines+=1),pParams->pszElType,pParams->pszRefSpecies,pszChrom,StartLoci,EndLoci,EndLoci - StartLoci + 1,pParams->pszRelSpecies);
if(pParams->RsltLineBuffOfs > ((sizeof(pParams->szRsltLineBuff) * 4) / 5))
	{
	CUtility::SafeWrite(pParams->hRsltsLociFile,pParams->szRsltLineBuff,pParams->RsltLineBuffOfs);
	_commit(pParams->hRsltsLociFile);
	pParams->RsltLineBuffOfs = 0;
	}
return(eBSFSuccess);
}



