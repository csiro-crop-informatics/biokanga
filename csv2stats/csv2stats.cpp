// csv2stats.cpp : Defines the entry point for the console application.
//


#include "stdafx.h"
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "../libbiokanga/commhdrs.h"

const char *cpszProgVer = "1.0.0";		// increment with each release
const int cDfltMinLengthRange = 4;		// accept sequence lengths from this min length up
const int cMaxLengthRange = 500000000;	// max sequence length accepted

int 
Process(bool bSkipFirst,			// true if first line contains header and should be skipped
		int MinLength,				// core elements must be of at least this length
		int MaxLength,				// and no longer than this length
		char *pszInLociFile,		// CSV file containing elements
		char *pszInSeqFile,			// file containing genome assembly sequences
		char *pszRsltsFile);		// file to write results into

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
bool bSkipFirst;			// true if first line contains header and should be skipped
int iMinLength;				// core elements must be of at least this length
int iMaxLength;				// and no longer than this length
char szInLociFile[_MAX_PATH];	// input element loci from this file
char szInSeqFile[_MAX_PATH];	// input bioseq file containing assembly
char szRsltsFile[_MAX_PATH];	// output stats to this file


// command line args
struct arg_lit  *help    = arg_lit0("hH","help",                "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_int *ScreenLogLevel=arg_int0("S", "ScreenLogLevel",	"<int>","Level of diagnostics written to logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"diagnostics log file");

struct arg_file *InLociFile = arg_file1("i","inloci","<file>",	"element loci CSV file");
struct arg_file *InSeqFile = arg_file1("I","assembly","<file>",	"genome assembly bioseq file");
struct arg_file *RsltsFile = arg_file1("o","output","<file>",	"output file");
struct arg_lit  *SkipFirst    = arg_lit0("x","skipfirst",       "skip first line of CSV - header line");
struct arg_int  *MinLength = arg_int0("l","minlength","<int>",	"minimum element length (default 10)");
struct arg_int  *MaxLength = arg_int0("L","maxlength","<int>",	"maximum element length (default 1000000000)");
struct arg_end *end = arg_end(20);

void *argtable[] = {help,version,FileLogLevel,ScreenLogLevel,LogFile,
					InLociFile,InSeqFile,RsltsFile,SkipFirst,MinLength,MaxLength,
					end};

char **pAllArgs;
int argerrors;
argerrors = CUtility::arg_parsefromfile(argc,(char **)argv,&pAllArgs);
if(argerrors >= 0)
	argerrors = arg_parse(argerrors,pAllArgs,argtable);

    /* special case: '--help' takes precedence over error reporting */
if (help->count > 0)
        {
		printf("\n%s csv2stats, Version %s\nOptions ---\n", gszProcName,cpszProgVer);
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

	bSkipFirst = SkipFirst->count ? true : false;

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
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Version: %s Processing parameters:",cpszProgVer);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Input CSV element loci file: '%s'",szInLociFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Input bioseq genome assembly file: '%s'",szInSeqFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Output to file: '%s'",szRsltsFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"First line contains header: %s",bSkipFirst ? "yes" : "no");
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Minimum element length: %d",iMinLength);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Maximum element length: %d",iMaxLength);

#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	// processing here...
	Rslt = Process(bSkipFirst,iMinLength,iMaxLength,szInLociFile,szInSeqFile,szRsltsFile);

	gStopWatch.Stop();
	Rslt = Rslt >=0 ? 0 : 1;
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Exit code: %d Total processing time: %s",Rslt,gStopWatch.Read());
	exit(Rslt);
	}
else
	{
	printf("\n%s csv2stats, Version %s\n",gszProcName,cpszProgVer);
	arg_print_errors(stdout,end,gszProcName);
	arg_print_syntax(stdout,argtable,"\nUse '-h' to view option and parameter usage\n");
	exit(1);
	}
}

typedef struct TAG_sBaseCompStats {
	int TotBases;					// total number of bases (canonical and non-canonical)
	int CBases;						// total number of all canonical bases
	int NBases;						// total of all non-canonical bases
	int MonoCnts[4];				// raw counts for a,c,g,t
	int DiCnts[4*4];				// dinucleotide counts indexed by aa,ac,ag,at,ca,...,ta,tc,tg,tt
	int TriCnts[4*4*4];				// trinucleotide counts indexed by aaa,aac,aag,aat,tta,ttc,ttg,ttt
} tsBaseCompStats;


char *
Idx2Seq(int NMer,		  // number of bases encoded by EncodedBases (1..3)
		int EncodedBases) // bases encoded as index into BaseCnts[]...TrinucCnts[]
{
static char szBases[80];
char *pBase;
int Base;
int Idx;
pBase = &szBases[NMer-1];
for(Idx = 0; Idx < NMer; Idx++,pBase--)
	{
	Base = EncodedBases & 0x03;
	EncodedBases >>= 2;
	switch(Base) {
		case eBaseA:
			*pBase = 'A';
			break;
		case eBaseC:
			*pBase = 'C';
			break;
		case eBaseG:
			*pBase = 'G';
			break;
		case eBaseT:
			*pBase = 'T';
			break;
		}
	}
szBases[NMer] = '\0';
return(szBases);
}

int
ReportStats(int hFile,tsBaseCompStats *pStats)
{
char szBuff[10000];
int Len;
int Idx;
Len = sprintf(szBuff,"0,\"TotalBases\",%d\n",pStats->TotBases);
Len += sprintf(&szBuff[Len],"1,\"TotalCanonical\",%d\n",pStats->CBases);
Len += sprintf(&szBuff[Len],"2,\"TotalNoncanonical\",%d\n",pStats->NBases);

for(Idx = 0; Idx < 4; Idx++)
	Len += sprintf(&szBuff[Len],"%d,\"%s\",%d\n",3+Idx,Idx2Seq(1,Idx),pStats->MonoCnts[Idx]);

for(Idx = 0; Idx < (4*4); Idx++)
	Len += sprintf(&szBuff[Len],"%d,\"%s\",%d\n",7+Idx,Idx2Seq(2,Idx),pStats->DiCnts[Idx]);

for(Idx = 0; Idx < (4*4*4); Idx++)
	Len += sprintf(&szBuff[Len],"%d,\"%s\",%d\n",23+Idx,Idx2Seq(3,Idx),pStats->TriCnts[Idx]);
return(write(hFile,szBuff,Len));
}

int
GenStats(int SeqLen,etSeqBase *pSeq,tsBaseCompStats *pStats)
{
etSeqBase Base;
int DiIdx;
int TriIdx;
int DiCnts;
int TriCnts;
int Cnt;

DiIdx = TriIdx = 0;
DiCnts = TriCnts = 0;
for(Cnt = 0; Cnt < SeqLen; Cnt++)
	{
	pStats->TotBases += 1;
	Base = *pSeq++ & ~cRptMskFlg;
	switch(Base) {
		case eBaseA:
		case eBaseC:
		case eBaseG:
		case eBaseT:
			pStats->MonoCnts[(int)Base] += 1;
			pStats->CBases += 1;
			break;

		default:
			pStats->NBases += 1;
			DiCnts = TriCnts = 0;
			continue;
		}
	

	DiIdx <<= 2;
	DiIdx &= 0x0f;
	DiIdx |= (int)Base;
	if(DiCnts++ >= 1)
		pStats->DiCnts[DiIdx] += 1;

	TriIdx <<= 2;
	TriIdx &= 0x03f;
	TriIdx |= (int)Base;
	if(TriCnts++ >= 2)
		pStats->TriCnts[TriIdx] += 1;
	}
return(SeqLen);
}

int 
Process(bool bSkipFirst,			// true if first line contains header and should be skipped
		int MinLength,				// core elements must be of at least this length
		int MaxLength,				// and no longer than this length
		char *pszInLociFile,		// CSV file containing elements
		char *pszInSeqFile,			// file containing genome assembly sequences
		char *pszRsltsFile)			// file to write fasta into
{
int NumFields;

int Rslt;
int SrcID;
char *pszChrom;
int ChromID;
char *pszElType;
char *pszRefSpecies;

int StartLoci;
int EndLoci;
int Len;

int ReqAllocLen;
int AllocSeqBuffLen = 0;
int AllocLineBuffLen = 0;

int NumAccepted;
int NumProcessed;
int NumUnderLen;
int NumOverLen;

tsBaseCompStats CompStats;

etSeqBase *pElSeqBuff = NULL;
int hRsltFile = -1;
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


if((pBioseq = new CBioSeqFile) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CBioSeqFile object");
	delete pCSV;
	return(eBSFerrObj);
	}

if((Rslt = pBioseq->Open(pszInSeqFile))!=eBSFSuccess)
	{
	while(pBioseq->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,pBioseq->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open assembly sequence file '%s'",pszInSeqFile);
	delete pCSV;
	delete pBioseq;
	return(Rslt);
	}

#ifdef _WIN32
if((hRsltFile = open(pszRsltsFile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
if((hRsltFile = open(pszRsltsFile, O_RDWR | O_CREAT,S_IREAD | S_IWRITE))==-1)
#endif
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create %s - %s",pszRsltsFile,strerror(errno));
	delete pCSV;
	delete pBioseq;
	return(eBSFerrCreateFile);
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Output file created/truncated: '%s'",pszRsltsFile);

memset(&CompStats,0,sizeof(CompStats));

NumUnderLen = 0;
NumOverLen = 0;
NumAccepted = 0;
NumProcessed = 0;

while((Rslt=pCSV->NextLine()) > 0)	// onto next line containing fields
	{
	NumFields = pCSV->GetCurFields();
	if(NumFields < 7)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expected 7+ fields in '%s', GetCurFields() returned '%d'",pszInLociFile,NumFields);
		Rslt = eBSFerrFieldCnt;
		break;
		}

	if(bSkipFirst || (!NumProcessed && pCSV->IsLikelyHeaderLine()))
		{
		bSkipFirst = false;
		continue;
		}

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

	ReqAllocLen = Len;
	if(pElSeqBuff == NULL || ReqAllocLen > AllocSeqBuffLen)
		{
		if(pElSeqBuff != NULL)
			delete pElSeqBuff;
		ReqAllocLen += 10000;		  // a little extra to reduce number of potential subsequent reallocations
		if((pElSeqBuff = (etSeqBase *)new unsigned char[ReqAllocLen])==NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate %s bytes as a sequence buffer",ReqAllocLen);
			Rslt = eBSFerrMem;
			break;
			}
		AllocSeqBuffLen = ReqAllocLen;
		}

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

	GenStats(Len,pElSeqBuff,&CompStats);
	NumAccepted += 1;
	}
if(Rslt >= eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Elements accepted: %d, Processed: %d, UnderLen: %d, OverLen: %d",
					NumAccepted,NumProcessed,NumUnderLen,NumOverLen);

	ReportStats(hRsltFile,&CompStats);
	}

close(hRsltFile);
delete pCSV;
delete pBioseq;
if(pElSeqBuff != NULL)
	delete pElSeqBuff;


return(Rslt < 0 ? NumAccepted : Rslt);
}


