
#include "stdafx.h"
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#if _WIN32

#include "../libbiokanga/commhdrs.h"
#else
#include "../libbiokanga/commhdrs.h"
#endif

const unsigned int cProgVer = 305;				// increment with each release
const unsigned int cMaxNumRelFiles = 50;		// max number of relative files handled
const unsigned int cAllocNumRsltItems = 50000;	// allocate results items in this many elements

typedef enum eProcMode {
	eProcModeStandard = 0,				// default processing is Identity = Match/(Match+Mismatch)
	eProcModeIdentity,					// Identity = Match/(CoreLength)
	eProcModeAligned,					// Aligned = (Match+Mismatch)/CoreLength
	eProcModeScore						// use scores
} etProcMode;


typedef struct TAG_sRelCounts {
	int NumRels;						// number of times this set of counts was updated (0 if never accessed)
	int Unaligned;
	int Matches;
	int Mismatches;
	int InDels;
	int Score;
} tsRelCounts;

typedef struct TAG_sRsltsItem {
	int RefID;						// identifier as parsed from field 1 int reference CVS file
	char szSpecies[cMaxDatasetSpeciesChrom];				// species from field3
	char szChrom[cMaxDatasetSpeciesChrom];	// chromosome from field 4
	int StartLoci;					// start loci from field 5
	int EndLoci;					// end loci from field 6
	int Length;						// length from field 7
	int Region;						// region from field 9
	tsRelCounts Counts[cMaxNumRelFiles];	// counts of interest from each relative CVS files (unaligned,matches,mismatches,indels)
} tsRsltsItem;

// processing parameters
typedef struct TAG_sProcParams 
	{
	etProcMode ProcMode;			// processing mode
	char *pszRefFile;				// reference CSV file
	char *pszRelFile;				// relative filespec
	char *pszRsltsFile;				// results file
	int hRsltsFile;					// opened handle for results file otherwise -1

	int CurFileIdx;					// current relative file index-1 (0 if none)

	char szCurRelFiles[cMaxNumRelFiles][_MAX_FNAME];	// relative file names
	int MinLen;						// minimum accepted length
	int MaxLen;						// maximum accepted length

	int NumRsltItems;				// number of items 
	int AllocdRsltItems;			// number of items allocated
	tsRsltsItem *pRsltItems;		// pts to allocated result items

	CFilterRefIDs *pFilterRefIDs;	// used when filtering by RefID
	} tsProcParams; 


int Process(etProcMode iProcMode,	// processing mode 0: default
			char *pszRefFile,		// reference csv file to process
			char *pszFilterRefIDFile, // exclude any RefIDs in this filter file
			char *pszRelFile,		// relative csv file(s) to processs
			char *pszRsltfile,		// results file to generate
			int MinLen,				// minimum accepted length
			int MaxLen);			// maximum accepted length

tsRsltsItem *LocateRsltsItem(int RefID,tsProcParams *pParams);
int static SortRsltItems(const void *arg1, const void *arg2);

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
	return _T("ProcessCSVfiles");
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

int iProcMode;			// processing mode
char szRefFile[_MAX_PATH];	// nput from reference csv file
char szFilterRefIDFile[_MAX_PATH];  // exclude any RefIDs in this filter file
char szRelFile[_MAX_PATH];	// input from relative csv file(s)
char szRsltfile[_MAX_PATH];	// output to results CSV file
int iMinLen;			// minimum accepted length
int iMaxLen;			// maximum accepted length


// command line args
struct arg_lit  *help    = arg_lit0("hH","help",                "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_int *ScreenLogLevel=arg_int0("S", "ScreenLogLevel",	"<int>","Level of diagnostics written to logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"diagnostics log file");
struct arg_int *ProcMode = arg_int0("m", "mode","<int>",		"processing mode (default 0) 0: Identities (Match/Match+Mismatches), 1: Identities (Match/CoreLength), 2: Aligned Length 3: Score");

struct arg_file *RefFile = arg_file1("i","in","<file>",			"input from reference csv file");
struct arg_file *FilterRefIDFile = arg_file0("X",NULL,"<file>",	"filter out any ref or rel loci with RefIDs in this filter file");
struct arg_file *RelFile = arg_file1("I","in","<file>",			"input from relative csv file(s) - can use wildcards");
struct arg_file *Rsltfile= arg_file1("o","out","<file>",		"output to results CSV file");
struct arg_int *MinLen=arg_int0("l", "MinLen",	"<int>","minimum accepted length");
struct arg_int *MaxLen=arg_int0("L", "MaxLen",	"<int>","maximum accepted length");

struct arg_end *end = arg_end(20);
void *argtable[] = {help,version,FileLogLevel,ScreenLogLevel,LogFile,
			ProcMode,RefFile,FilterRefIDFile,RelFile,Rsltfile,
			MinLen,MaxLen,
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

	iProcMode = ProcMode->count ? ProcMode->ival[0] : eProcModeStandard;
	if(iProcMode < eProcModeStandard || iProcMode > eProcModeScore)
		{
		printf("\nError: Requested processing mode '-x%d' not supported",iProcMode);
		exit(1);
		}

	iMinLen = MinLen->count ? MinLen->ival[0] : 0;
	if(iMinLen < 0)
		{
		printf("\nError: Requested minimum length '-l%d' is negative",iMinLen);
		exit(1);
		}

	iMaxLen = MaxLen->count ? MaxLen->ival[0] : 1000000000;
	if(iMaxLen < iMinLen)
		{
		printf("\nError: Requested maximum ength '-l%d' must be >= minimum %d",iMaxLen,iMinLen);
		exit(1);
		}



	strcpy(szRefFile,RefFile->filename[0]);
	strcpy(szRelFile,RelFile->filename[0]);
	strcpy(szRsltfile,Rsltfile->filename[0]);

	if(FilterRefIDFile->count)
		{
		strncpy(szFilterRefIDFile,FilterRefIDFile->filename[0],sizeof(szFilterRefIDFile));
		szFilterRefIDFile[sizeof(szFilterRefIDFile)-1] = '\0';
		}
	else
		szFilterRefIDFile[0] = '\0';


			// now that command parameters have been parsed then initialise diagnostics log system
	if(!gDiagnostics.Open(szLogFile,(etDiagLevel)iScreenLogLevel,(etDiagLevel)iFileLogLevel,true))
		{
		printf("\nError: Unable to start diagnostics subsystem.");
		if(szLogFile[0] != '\0')
			printf(" Most likely cause is that logfile '%s' can't be opened/created",szLogFile);
		exit(1);
		}

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Version: %d.%2.2d Processing parameters:",cProgVer/100,cProgVer%100);
	switch(iProcMode) {
		case eProcModeStandard:
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Processing mode:  Identity = Match/(Match+Mismatch)");
			break;
		case eProcModeIdentity:
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Processing mode:  Identity = Match/(CoreLength)");
			break;
		case eProcModeAligned:
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Processing mode:  Aligned = (Match+Mismatch)/CoreLength");
			break;
		case eProcModeScore:
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Processing mode: score");
			break;
		}

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"reference csv file to process: '%s'",szRefFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"relative csv file(s) to processs: '%s'",szRelFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"results file to generate: '%s'",szRsltfile);
	if(szFilterRefIDFile[0])
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Exclude any RefIDs in this filter file: '%s'",szFilterRefIDFile);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"minimum length: %d",iMinLen);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"maximum length: %d",iMaxLen);

	gStopWatch.Start();
#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	Rslt = Process((etProcMode)iProcMode,		// processing mode 0: default
					szRefFile,		// reference csv file to process
					szFilterRefIDFile, // exclude any RefIDs in this filter file
					szRelFile,		// relative csv file(s) to processs
					szRsltfile,  	// results file to generate
					iMinLen,		// minimum accepted ref length
					iMaxLen);		// maximum accepted ref length

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
return 0;
}


// ProcessThisFile
bool ProcessThisFile(char *pszFile,
				 void *pParams)			
{
int Rslt;
int LineLen;
int NumFields;
int MinFields;
int SrcID;
char *pszChrom;
char *pszElType;
char *pszRefSpecies;
char *pszRelSpecies;
int StartLoci;
int EndLoci;
int Len;
int Features;
int Unaligned;
int Matches;
int Mismatches;
int InDels;
int Score;

int CountIdx;
int NumElsRead;		// number of elements before filtering
int NumElsAccepted;	// number of elements accepted after filtering
int NumFiltRefIDs;	// number of elements filtered out because of matching RefID
int NumFiltLen;		// number of elements filtered out because of length
int NumUnlocated;	// number of elements which couldn't be located

tsRsltsItem *pRsltItem;


tsProcParams *pProcParams = (tsProcParams *)pParams;	// saves a little type casting later!

if(pProcParams->CurFileIdx == cMaxNumRelFiles)			// ensure we don't process more relative files than configured for
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to process file '%d' as compiled to handle a maximum of %d files",pszFile,cMaxNumRelFiles);
	return(false);
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing file: %s",pszFile);
CountIdx = pProcParams->CurFileIdx++;
strcpy(pProcParams->szCurRelFiles[CountIdx],pszFile);

CCSVFile *pCSV = new CCSVFile;
if(pCSV == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CCSVfile");
	return(false);
	}

if((Rslt=pCSV->Open(pszFile))!=eBSFSuccess)
	{
	while(pCSV->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,pCSV->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open file: %s",pszFile);
	delete pCSV;
	return(false);
	}


if(pProcParams->ProcMode == eProcModeScore)
	MinFields = 14;
else
	MinFields = 13;

NumElsRead =0;		// number of elements before filtering
NumElsAccepted =0;	// number of elements accepted after filtering
NumFiltRefIDs =0;	// number of elements filtered out because of matching RefID
NumFiltLen=0;		// number of elements filtered out because of length
NumUnlocated=0;	// number of elements which couldn't be located
LineLen = 0;
while((Rslt=pCSV->NextLine()) > 0)	// onto next line containing fields
	{
	NumFields = pCSV->GetCurFields();
	if(NumFields < MinFields)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expected at least %d fields in '%s', GetCurFields() returned '%d'",MinFields,pszFile,NumFields);
		delete pCSV;
		return(false);
		}
	NumElsRead += 1;
	pCSV->GetInt(1,&SrcID);
	if(pProcParams->pFilterRefIDs != NULL && pProcParams->pFilterRefIDs->Locate(SrcID))
		{	
		NumFiltRefIDs += 1;
		continue;
		}

	pCSV->GetText(2,&pszElType);
	pCSV->GetText(3,&pszRefSpecies);
	pCSV->GetText(4,&pszChrom);
	pCSV->GetInt(5,&StartLoci);
	pCSV->GetInt(6,&EndLoci);
	pCSV->GetInt(7,&Len);
	if(Len < pProcParams->MinLen || Len > pProcParams->MaxLen)
		{
		NumFiltLen += 1;
		continue;
		}

	pCSV->GetText(8,&pszRelSpecies);
	pCSV->GetInt(9,&Features);

	pCSV->GetInt(10,&Unaligned);
	pCSV->GetInt(11,&Matches);
	pCSV->GetInt(12,&Mismatches);
	pCSV->GetInt(13,&InDels);

	if(NumFields >= 14)
		pCSV->GetInt(14,&Score);
	else
		Score = 0;
	if((pRsltItem = LocateRsltsItem(SrcID,pProcParams))==NULL)
		{
		if(NumUnlocated++ < 10)
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"Couldn't locate element %d",SrcID);
		continue;
		}
	
	pRsltItem->Counts[CountIdx].NumRels += 1;
	pRsltItem->Counts[CountIdx].Unaligned += Unaligned;
	pRsltItem->Counts[CountIdx].Matches += Matches;
	pRsltItem->Counts[CountIdx].Mismatches += Mismatches;
	pRsltItem->Counts[CountIdx].InDels += InDels;

	pRsltItem->Counts[CountIdx].Score += Score;
	NumElsAccepted++;
	}
delete pCSV;

if(Rslt >= 0 && !NumElsAccepted)
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"No elements parsed from '%s'",pszFile);
else
	if(Rslt >= 0)
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"%s: Read: %d, Accepted: %d, FiltRefIDs: %d, FiltLen: %d, Unlocated: %d",pszFile,NumElsRead,NumElsAccepted,NumFiltRefIDs,NumFiltLen,NumUnlocated);
return(Rslt >= 0 ? true : false);
}

tsRsltsItem *
AddRsltsItem(int RefID, char *pszSpecies, char *pszChrom, int StartLoci, int EndLoci, int Length, int Region,tsProcParams *pParams)
{
tsRsltsItem *pItem;
if(pParams->pRsltItems == NULL || pParams->NumRsltItems == pParams->AllocdRsltItems)
	{
	pItem = new tsRsltsItem [pParams->AllocdRsltItems + cAllocNumRsltItems];
	if(pParams->pRsltItems != NULL && pParams->NumRsltItems > 0)
		{
		memcpy(pItem,pParams->pRsltItems,pParams->NumRsltItems * sizeof(tsRsltsItem));
		delete pParams->pRsltItems;
		}
	else
		{
		pParams->AllocdRsltItems = 0;
		pParams->NumRsltItems = 0;
		}
	pParams->AllocdRsltItems += cAllocNumRsltItems;
	pParams->pRsltItems = pItem;
	}
pItem = &pParams->pRsltItems[pParams->NumRsltItems++];
memset(pItem,0,sizeof(tsRsltsItem));
pItem->RefID = RefID;						// identifier as parsed from field 1 int reference CVS file
strcpy(pItem->szSpecies,pszSpecies);		// species from field3
strcpy(pItem->szChrom,pszChrom);	// chromosome from field 4
pItem->StartLoci = StartLoci;					// start loci from field 5
pItem->EndLoci = EndLoci;					// end loci from field 6
pItem->Length = Length;						// length from field 7
pItem->Region = Region;						// region from field 9
return(pItem);
}

int
LoadRefFile(tsProcParams *pProcParams)
{
int Rslt;
int LineLen;
int NumFields;
int SrcID;
char *pszChrom;
char *pszElType;
char *pszRefSpecies;
char *pszRelSpecies;
int StartLoci;
int EndLoci;
int Len;
int Features;

int NumElsRead;		// number of elements before filtering
int NumElsAccepted;	// number of elements accepted after filtering
int NumFiltRefIDs;	// number of elements filtered out because of matching RefID
int NumFiltLen;		// number of elements filtered out because of length

CCSVFile *pCSV = new CCSVFile;
if(pCSV == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CCSVfile");
	return(eBSFerrObj);
	}

if((Rslt=pCSV->Open(pProcParams->pszRefFile))!=eBSFSuccess)
	{
	while(pCSV->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,pCSV->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open file: %s",pProcParams->pszRefFile);
	delete pCSV;
	return(Rslt);
	}

NumElsRead =0;		// number of elements before filtering
NumElsAccepted =0;	// number of elements accepted after filtering
NumFiltRefIDs =0;	// number of elements filtered out because of matching RefID
NumFiltLen=0;		// number of elements filtered out because of length
LineLen = 0;
while((Rslt=pCSV->NextLine()) > 0)	// onto next line containing fields
	{
	NumFields = pCSV->GetCurFields();
	if(NumFields < 9)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expected at least 9 fields in '%s', GetCurFields() returned '%d'",pProcParams->pszRefFile,NumFields);
		delete pCSV;
		return(eBSFerrFieldCnt);
		}

	NumElsRead += 1;
	pCSV->GetInt(1,&SrcID);
	if(pProcParams->pFilterRefIDs != NULL && pProcParams->pFilterRefIDs->Locate(SrcID))
		{	
		NumFiltRefIDs += 1;
		continue;
		}

	pCSV->GetText(2,&pszElType);
	pCSV->GetText(3,&pszRefSpecies);
	pCSV->GetText(4,&pszChrom);
	pCSV->GetInt(5,&StartLoci);
	pCSV->GetInt(6,&EndLoci);
	pCSV->GetInt(7,&Len);
	if(Len < pProcParams->MinLen || Len > pProcParams->MaxLen)
		{
		NumFiltLen += 1;
		continue;
		}

	pCSV->GetText(8,&pszRelSpecies);
	pCSV->GetInt(9,&Features);

	if(AddRsltsItem(SrcID,pszRefSpecies,pszChrom,StartLoci,EndLoci,Len,Features,pProcParams)==NULL)
		{
		Rslt = eBSFerrMem;
		break;
		}
	NumElsAccepted++;
	}
delete pCSV;

if(Rslt >= 0 && !pProcParams->NumRsltItems)
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"No elements parsed from '%s'",pProcParams->pszRefFile);

if(Rslt >= 0 && pProcParams->NumRsltItems > 1)
	qsort(pProcParams->pRsltItems,pProcParams->NumRsltItems,sizeof(tsRsltsItem),SortRsltItems);

// check for duplicate refs
tsRsltsItem *pItem = pProcParams->pRsltItems;
int CurRefID = -1;
int PrevDupRefID = -1;
int NumDuplicates = 0;
for(int Idx = 0; Idx < pProcParams->NumRsltItems; Idx++,pItem++)
	{
	if(pItem->RefID == CurRefID && CurRefID != PrevDupRefID)
		{
		gDiagnostics.DiagOut(eDLWarn,gszProcName,"Duplicated RefID %d",CurRefID);
		PrevDupRefID = CurRefID;
		NumDuplicates += 1;
		}
	CurRefID = pItem->RefID;
	}

if(Rslt >= 0)
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"%s: Read: %d, Accepted: %d, FiltRefIDs: %d, FiltLen: %d, Duplicates: %d",
				pProcParams->pszRefFile,NumElsRead,NumElsAccepted,NumFiltRefIDs,NumFiltLen,NumDuplicates);

return(Rslt >= 0 ? pProcParams->NumRsltItems : Rslt);
}

tsRsltsItem *
LocateRsltsItem(int RefID,tsProcParams *pParams)
{
tsRsltsItem *pProbe;
if(RefID <= 0)
	return(NULL);
int Lo,Mid,Hi;	// search limits
Lo = 0; Hi = pParams->NumRsltItems-1;
while(Hi >= Lo) {
	Mid = (Hi + Lo)/2;
	pProbe = &pParams->pRsltItems[Mid];
	if(pProbe->RefID == RefID)
		return(pProbe);
	if(RefID < pProbe->RefID)	
		{
		Hi = Mid - 1;
		continue;
		}
	Lo = Mid + 1;
	}
return(NULL);
}



int
GenResults(char *pszRsltsFile,tsProcParams *pParams)
{
char szLineBuff[2048];
int Len;
tsRsltsItem *pCurRslt;
tsRelCounts *pCounts;
int RsltIdx;
int SpeciesIdx;
double Identity;

#ifdef _WIN32
if((pParams->hRsltsFile = open(pszRsltsFile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
if((pParams->hRsltsFile = open(pszRsltsFile, O_RDWR | O_CREAT | O_TRUNC,S_IREAD | S_IWRITE))==-1)
#endif
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create output results file: %s - %s",pszRsltsFile,strerror(errno));
	return(eBSFerrCreateFile);
	}

Len = sprintf(szLineBuff,"\"RefID\",\"Chrom\",\"Start\",\"End\",\"Length\"");
for(SpeciesIdx = 0; SpeciesIdx < pParams->CurFileIdx; SpeciesIdx++)
	{
	char szFname[_MAX_FNAME];
#ifdef _WIN32
	_splitpath(pParams->szCurRelFiles[SpeciesIdx],NULL,NULL,szFname,NULL);
#else
	CUtility::splitpath(pParams->szCurRelFiles[SpeciesIdx],NULL,szFname);
#endif
	Len += sprintf(&szLineBuff[Len],",\"%s\"",szFname);
	}
Len += sprintf(&szLineBuff[Len],"\n");
CUtility::SafeWrite(pParams->hRsltsFile,szLineBuff,Len);

pCurRslt = pParams->pRsltItems;
for(RsltIdx = 0; RsltIdx < pParams->NumRsltItems; RsltIdx++,pCurRslt++)
	{
	Len = sprintf(szLineBuff,"%d,\"%s\",%d,%d,%d",pCurRslt->RefID,pCurRslt->szChrom,pCurRslt->StartLoci,pCurRslt->EndLoci,pCurRslt->Length);
	pCounts = &pCurRslt->Counts[0];
	for(SpeciesIdx = 0; SpeciesIdx < pParams->CurFileIdx; SpeciesIdx++,pCounts++)
		{
		if(pCounts->NumRels > 1)
			gDiagnostics.DiagOut(eDLWarn,gszProcName,"Rel file '%s' referenced RefID %d multiple (%d) times",
											pParams->szCurRelFiles[SpeciesIdx],
											pCurRslt->RefID,
											pCounts->NumRels);
		switch(pParams->ProcMode) {
			case eProcModeStandard:
				if(!(pCounts->Matches + pCounts->Mismatches))
					Identity = 0.00;
				else
					Identity = (pCounts->Matches * 100.0)/(double)(pCounts->Matches + pCounts->Mismatches);
				break;

			case eProcModeIdentity:
				Identity = (pCounts->Matches * 100.0)/(double)pCurRslt->Length;
				break;

			case eProcModeAligned:
				Identity = ((pCounts->Matches  + pCounts->Mismatches) * 100.0)/(double)pCurRslt->Length;
				if(Identity > 100.0)				// could have been more than 100% if core had InDels relative to outspecies
					Identity = 100.0;
				break;

			case eProcModeScore:
				Identity = (double)pCounts->Score/10.0;
				break;
			}

		Len += sprintf(&szLineBuff[Len],",%2.3f",Identity);
		}
	Len += sprintf(&szLineBuff[Len],"\n");
	CUtility::SafeWrite(pParams->hRsltsFile,szLineBuff,Len);
	}

if(pParams->hRsltsFile != -1)
	{
	close(pParams->hRsltsFile);
	pParams->hRsltsFile = -1;
	}
return(eBSFSuccess);
}


int Process(etProcMode ProcMode,	// processing mode 0: default
			char *pszRefFile,		// reference csv file to process
			char *pszFilterRefIDFile, // exclude any RefIDs in this filter file
			char *pszRelFile,		// relative csv file(s) to processs
			char *pszRsltfile,		// results file to generate
			int MinLen,			// minimum accepted length
			int MaxLen)			// maximum accepted length
{
int Rslt;
tsProcParams ProcParams;

memset(&ProcParams,0,sizeof(ProcParams));
ProcParams.ProcMode = ProcMode;
ProcParams.pszRefFile   = pszRefFile;
ProcParams.pszRsltsFile = pszRsltfile;
ProcParams.pszRelFile = pszRelFile;
ProcParams.MinLen = MinLen;
ProcParams.MaxLen = MaxLen;
ProcParams.hRsltsFile = -1;

if(pszFilterRefIDFile != NULL && 
   pszFilterRefIDFile[0] != '\0')
	{
	ProcParams.pFilterRefIDs = new CFilterRefIDs;
	if((Rslt=ProcParams.pFilterRefIDs->Open(pszFilterRefIDFile)) < 0)
		{
		while(ProcParams.pFilterRefIDs->NumErrMsgs())
			gDiagnostics.DiagOut(eDLFatal,gszProcName,ProcParams.pFilterRefIDs->GetErrMsg());
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open/parse input filter RefID file %s",pszFilterRefIDFile);
		return(Rslt);
		}
	}
 
// load reference file
if((Rslt = LoadRefFile(&ProcParams)) <= 0)
	{
	if(ProcParams.pFilterRefIDs != NULL)
		delete ProcParams.pFilterRefIDs;
	return(Rslt);
	}

CSimpleGlob glob(SG_GLOB_FULLSORT);
if (glob.Add(pszRelFile) >= SG_SUCCESS)
	{
	for (int n = 0; Rslt >= eBSFSuccess &&  n < glob.FileCount(); ++n)
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Will process in this order: %d '%s", n+1,glob.File(n));

    for (int n = 0; Rslt >= eBSFSuccess &&  n < glob.FileCount(); ++n)
		Rslt = ProcessThisFile(glob.File(n),&ProcParams);
	}
else
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to glob '%s",pszRelFile);
	Rslt = eBSFerrOpnFile;	// treat as though unable to open file
	}
if(Rslt >= eBSFSuccess)
	GenResults(pszRsltfile,&ProcParams);

if(ProcParams.hRsltsFile != -1)
	close(ProcParams.hRsltsFile);

if(ProcParams.pRsltItems != NULL)
	delete ProcParams.pRsltItems;
if(ProcParams.pFilterRefIDs != NULL)
	delete ProcParams.pFilterRefIDs;
return(Rslt);
}

int 
SortRsltItems(const void *arg1, const void *arg2)
{
tsRsltsItem *pEl1 = (tsRsltsItem *)arg1;
tsRsltsItem *pEl2 = (tsRsltsItem *)arg2;
if(pEl1->RefID < pEl2->RefID)
	return(-1);
if(pEl1->RefID > pEl2->RefID)
	return(1);
return(0);
}
