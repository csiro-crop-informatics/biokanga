// FindApproxMatches.cpp : Defines the entry point for the console application.
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "stdafx.h"

#if _WIN32

#include "../conservlib/commhdrs.h"

#else
#include "../libbiokanga/commhdrs.h"
#endif

CStopWatch gStopWatch;
CDiagnostics gDiagnostics;				// for writing diagnostics messages to log file
char gszProcName[_MAX_FNAME];			// process name

const unsigned int cProgVer = 201;		// increment with each release


int gMatchID = 0;						// used to keep a count of all matches generated

const int cFMBuffChunk = 0x0fffff;	// allocate FMIndex image size plus this chunk to reduce the number of subsequent allocations

const int cMaxFAMProbeLen = 10000;		   // maximum allowed probe length

typedef struct TAG_sFAMProbe {
	struct TAG_sFAMProbe *pNext;				// probes are linked
	int ProbID;								// uniquely identifies this probe
	char szDescr[80];						// as parsed from fasta descriptor
	int TotPlusHits;						// total hits on '+' strand
	int TotMinHits;							// total hits on '-' strand
	int ProbeLen;							// number of bases in Bases[]
	etSeqBase Bases[1];						// to hold probe sequence
} sFAMProbe;

// holds approximate matches as linked list prior to sorting when deduping
typedef struct TAG_sFAMMatch {
	struct TAG_sFAMMatch *pNext;			// matches are linked
	int MatchID;							// uniquely identifies this match
	int MatchLen;							// number of bases in Bases[]
	etSeqBase Bases[1];						// to hold match sequence
} sFAMMatch;


int Process(char *pszInputFile,char *pszProbeFile,char *pszOutputFile,char Strand,int Identity,int FlankLen,bool bCount,bool bGlobCounts);
int LoadProbes(char *pszProbesFile,struct TAG_sFAMProbe **ppProbes);
void DeleteProbes(struct TAG_sFAMProbe *pProbes);
int ApproxMatch(bool bCounts,bool bGlobCounts,sFAMProbe *pProbe,int hRsltsFile,CFMIndex *pFMIndex,char *pszChrom,char Strand,char *pszDescr,unsigned char *Pattern, int PatternLen, UINT32 Identity);
int Extract(int hRsltsFile,CFMIndex *pFMIndex,char *pszChrom,UINT32 position, UINT32 nchars);
int Display(int hRsltsFile,CFMIndex *pFMIndex,char *pszChrom,char Strand,unsigned char *Pattern, int PatternLen,UINT32 nFlankChars);
int Locate(bool bCounts,bool bGlobCounts,sFAMProbe *pProbe,int hRsltsFile,CFMIndex *pFMIndex,char *pszChrom,char Strand,char *pszDescr,unsigned char *Pattern,int PatternLen); 
int Count(bool bGlobCounts,sFAMProbe *pProbe,int hRsltsFile,CFMIndex *pFMIndex,char *pszChrom,char Strand,char *pszDescr,unsigned char *Pattern,int PatternLen); 
int SortMatches(const void *arg1, const void *arg2);


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
int iIdentity;				// minimum matching identity
int iFlankLen;				// flanking sequence length to return with probe matches
char cStrand;				// which strand '+' or '-' or '*'
bool bCounts;				// true if match counts per chromosome required
bool bGlobCounts;				// true if match counts over whole genome required


char szOutputFile[_MAX_PATH];
char szInputFile[_MAX_PATH];
char szProbesFile[_MAX_PATH];

// command line args
struct arg_lit  *help    = arg_lit0("hH","help",                "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_int *ScreenLogLevel=arg_int0("S", "ScreenLogLevel",	"<int>","Level of diagnostics written to screen 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"diagnostics log file");

struct arg_file *InFile = arg_file1("i","inindex","<file>",			"input from bioseq FMIndex file");
struct arg_file *ProbesFile = arg_file1("I","inprobes","<file>",	"input from multifasta probe file");
struct arg_file *OutFile = arg_file1("o","outmatches","<file>",		"output matches to this file as CSV");
struct arg_int *Identity=arg_int0("M", "identity",	"<int>",		"Minimum matching identity (70..100)");
struct arg_int *MaxMismatches=arg_int0("m", "maxmissmatches","<int>","Maximum number of mismatches (0..10) overides -I<identity>");
struct arg_str *Strand=arg_str0("s", "strand",	"<string>",			"Strand '+' or '-', default is '*' for both");
struct arg_lit  *counts = arg_lit0("c","counts",					"output hit counts only (per chromosome)");
struct arg_lit  *globcounts = arg_lit0("g","counts",				"output hit counts only (over whole genome)");
struct arg_int *FlankLen=arg_int0("f", "len",	"<int>",			"flanking sequence length");
struct arg_end *end = arg_end(20);

void *argtable[] = {help,version,FileLogLevel,ScreenLogLevel,LogFile,
					InFile,OutFile,ProbesFile,Identity,MaxMismatches,counts,globcounts,FlankLen,Strand,
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

	iFlankLen = FlankLen->count ? FlankLen->ival[0] : 0;
	if(iFlankLen < 0)
		iFlankLen = 0;

	// NOTE: identity will always be overridden by number of mismatches if mismatches are specified
	iIdentity = Identity->count ? Identity->ival[0] : 100;
	if(iIdentity < 70)
		iIdentity = 70;
	else
		if(iIdentity > 100)
			iIdentity = 100;

	if(MaxMismatches->count)
		{
		iIdentity = MaxMismatches->ival[0];
		if(iIdentity <= 0)
			iIdentity = 100;
		else
			if(iIdentity > 10)
				iIdentity = 10;
		}
	
	if(Strand->count)
		{
		cStrand = Strand->sval[0][0];
		if(cStrand != '*' && cStrand != '+' && cStrand != '-' && cStrand != '1' && cStrand != '0')
			{
			printf("\nError: Requested strand '-s%c' is invalid, must be '-s+' or '-s-' or '-s*'",cStrand);
			exit(1);
			}
		if(cStrand == '0')
			cStrand = '+';
		else
			if(cStrand == '0')
				cStrand = '-';
		}
	else
		cStrand = '*';

	if(InFile->count)
		strcpy(szInputFile,InFile->filename[0]);
	else
		strcpy(szInputFile,"in.fmi");


	if(OutFile->count)
		strcpy(szOutputFile,OutFile->filename[0]);
	else
		strcpy(szOutputFile,"out.csv");

	if(ProbesFile->count)
		strcpy(szProbesFile,ProbesFile->filename[0]);
	else
		strcpy(szProbesFile,"probes.fa");

	bGlobCounts = globcounts->count ? true : false;
	if(!bGlobCounts)
		bCounts = counts->count ? true : false;
	else
		bCounts = false;

		// now that command parameters have been parsed then initialise diagnostics log system
	if(!gDiagnostics.Open(szLogFile,(etDiagLevel)iScreenLogLevel,(etDiagLevel)iFileLogLevel,true))
		{
		printf("\nError: Unable to start diagnostics subsystem.");
		if(szLogFile[0] != '\0')
			printf(" Most likely cause is that logfile '%s' can't be opened/created",szLogFile);
		exit(1);
		}

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Version: %d.%2.2d Processing parameters:",cProgVer/100,cProgVer%100);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Input from bioseq FMIndex file: %s",szInputFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Probes from multifasta file: %s",szProbesFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Output as CSV file: '%s'",szOutputFile);
	if(bGlobCounts)
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Only hit counts per genome will be generated");
	else
		if(bCounts)
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Only hit counts per chromosome will be generated");
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Target strand: '%c'",cStrand);
	if(iIdentity <= 10)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Max mismatches: %d",iIdentity);
	else
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Minimum identity: %d",iIdentity);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Flanking sequence length: %d",iFlankLen);

#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	gStopWatch.Start();
	Rslt = Process(szInputFile,szProbesFile,szOutputFile,cStrand,iIdentity,iFlankLen,bCounts,bGlobCounts);
	gStopWatch.Stop();
	Rslt = Rslt >=0 ? 0 : 1;
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Exit code: %d Total processing time: %s",Rslt,gStopWatch.Read());
	exit(Rslt);
	}
else
	{
	arg_print_errors(stdout,end,gszProcName);
	arg_print_syntax(stdout,argtable,"\n [@<parameterfile>]\nend of help\n");
	exit(1);
	}
}


int 
Process(char *pszInputFile,char *pszProbesFile,char *pszOutputFile,char Strand,int Identity,int FlankLen,bool bCounts,bool bGlobCounts)
{
CBioSeqFile *pBioSeq;
CBEDfile *pBEDFile = NULL;
CFMIndex *pFMIndex;						// NOTE: instantiated and deleted for each new FMIndex image because I don't yet have confidence that all resources are being released by this class
int hOutFile;
sFAMProbe *pProbes;
sFAMProbe *pCurProbe;

int NumProbes;
char szChromName[cBSFSourceSize+1];		// to hold current chromosome name
tBSFEntryID CurEntry;					// current entry being processed

int FMIndexLen;							// length of FMIndex image being processed
UINT32 FMImageLen;						// image length as stored in bioseq entry
unsigned char *pFMBuff;					// buffer allocated to hold the FMIndex image
int FMBuffLen;							// current FMIndex image buffer size 
int Rslt;								// used to hold processing result
int FMRslt;								// result as returned by pFMINdex methods

unsigned char RevCplBuff[0x7fff];		// used to reverse complement the probes
unsigned char *pPattern;

int RsltsLen;



if((NumProbes=Rslt=LoadProbes(pszProbesFile,&pProbes))<=eBSFSuccess)
	return(Rslt);

// open bioseq file containing FMIndexes
if((pBioSeq = new CBioSeqFile()) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CBioSeqFile");
	DeleteProbes(pProbes);
	return(eBSFerrObj);
	}
if((Rslt=pBioSeq->Open(pszInputFile,cBSFTypeFMIndexes))!=eBSFSuccess)
	{
	while(pBioSeq->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,pBioSeq->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open input bioseq FMIndex file '%s'",pszInputFile);
	DeleteProbes(pProbes);
	delete pBioSeq;
	return(Rslt);
	}

// open and truncate CSV output file
#ifdef _WIN32
hOutFile = open(pszOutputFile, _O_WRONLY | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE );
#else
hOutFile = open(pszOutputFile, O_WRONLY | O_CREAT | O_TRUNC, S_IREAD | S_IWRITE );
#endif
if(hOutFile == -1)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open output CSV file '%s', error: %s",pszOutputFile,strerror(errno));
	DeleteProbes(pProbes);
	delete pBioSeq;
	return(eBSFerrOpnFile);
	}

	// iterate over each FMIndex
pFMIndex = NULL;
pFMBuff = NULL;
FMBuffLen = 0;
CurEntry = 0;
Rslt = eBSFSuccess;
while((CurEntry = pBioSeq->Next(CurEntry)) > 0)
	{
		// get entry name - assume it will be a chromosome/contig name
	pBioSeq->GetName(CurEntry,sizeof(szChromName),szChromName);
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing: %s",szChromName);

		// need to ensure buffer can hold the FMIndex
	FMIndexLen = pBioSeq->GetDataLen(CurEntry);
	if(pFMBuff == NULL || FMIndexLen > FMBuffLen)
		{
		if(pFMBuff != NULL)
			delete pFMBuff;
		FMBuffLen = FMIndexLen + cFMBuffChunk;
		pFMBuff = (unsigned char *)new unsigned char [FMBuffLen];
		if(pFMBuff == NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate memory (%d bytes) for holding FMIndex image",FMBuffLen);
			Rslt = eBSFerrMem;
			break;
			}
		}

	pBioSeq->GetData(CurEntry,eBinDataType,0,(unsigned char *)&FMImageLen,sizeof(UINT32));
	assert((FMImageLen + sizeof(UINT32)) == FMIndexLen);
	pBioSeq->GetData(CurEntry,eBinDataType,sizeof(UINT32),pFMBuff,FMImageLen);

	// instantiate FMIndex - when this class is known to release all resources correctly then can move outside of iterating loop
	if((pFMIndex = new CFMIndex()) == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CFMIndex");
		Rslt = eBSFerrObj;
		break;
		}
	if((FMRslt = pFMIndex->load_index_mem(pFMBuff,FMImageLen))!=0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to load image into FMIndex: %s",pFMIndex->GetErrText(FMRslt));
		Rslt = eBSFerrInternal;
		break;
		}

	pCurProbe = pProbes;
	while(pCurProbe != NULL)
		{
		switch(Strand) {
			case '-':
				memcpy(RevCplBuff,pCurProbe->Bases,pCurProbe->ProbeLen);
				CSeqTrans::ReverseComplement(pCurProbe->ProbeLen,RevCplBuff);
				pPattern = RevCplBuff;
				break;
			case '+':
			case '*':
				pPattern = pCurProbe->Bases;
				break;
			default:
				break;
			}

		if(Identity < 100)
			ApproxMatch(bCounts,bGlobCounts,pCurProbe,hOutFile,pFMIndex,szChromName,Strand == '*' ? '+' : Strand,pCurProbe->szDescr,pPattern,pCurProbe->ProbeLen,Identity);
		else
			Locate(bCounts,bGlobCounts,pCurProbe,hOutFile,pFMIndex,szChromName,Strand == '*' ? '+' : Strand,pCurProbe->szDescr,pPattern,pCurProbe->ProbeLen);

		if(Strand == '*')
			{
			memcpy(RevCplBuff,pCurProbe->Bases,pCurProbe->ProbeLen);
			CSeqTrans::ReverseComplement(pCurProbe->ProbeLen,RevCplBuff);
			pPattern = RevCplBuff;

			if(Identity < 100)
				ApproxMatch(bCounts,bGlobCounts,pCurProbe,hOutFile,pFMIndex,szChromName,'-',pCurProbe->szDescr,pPattern,pCurProbe->ProbeLen,Identity);
			else
				Locate(bCounts,bGlobCounts,pCurProbe,hOutFile,pFMIndex,szChromName,'-',pCurProbe->szDescr,pPattern,pCurProbe->ProbeLen);
			}
		pCurProbe = pCurProbe->pNext;
		}

	// get ready for next iteration
	delete pFMIndex;
	pFMIndex = NULL;
	}

if(bGlobCounts)
	{
	gMatchID = 0;
	pCurProbe = pProbes;
	while(pCurProbe != NULL)
		{
		if(pCurProbe->TotPlusHits || pCurProbe->TotMinHits)
			{
			RsltsLen = sprintf((char *)RevCplBuff,"%d,%d,\"%s\",%d,%d\n",++gMatchID,pCurProbe->ProbID,pCurProbe->szDescr,pCurProbe->TotPlusHits,pCurProbe->TotMinHits);
			write(hOutFile,RevCplBuff,RsltsLen);
			}
		pCurProbe = pCurProbe->pNext;
		}
	}
DeleteProbes(pProbes);
if(hOutFile!=-1)
	close(hOutFile);

if(pFMBuff!=NULL)
	delete pFMBuff;
if(pFMIndex != NULL)
	delete pFMIndex;
if(pBioSeq != NULL)
	delete pBioSeq;
if(pBEDFile != NULL)
	delete pBEDFile;
return(Rslt);
}



void
DeleteProbes(sFAMProbe *pProbes)
{
sFAMProbe *pCurProbe;
if(pProbes == NULL)
	return;
while((pCurProbe = pProbes) != NULL) 
	{
	pProbes = pCurProbe->pNext;
	delete pCurProbe;
	}
}

int
LoadProbes(char *pszProbesFile,struct TAG_sFAMProbe **ppProbes)
{
int Rslt;
sFAMProbe *pCurProbe;
sFAMProbe *pNxtProbe;
int NumProbes;
sFAMProbe *pProbes;
CFasta *pFasta;

char szDescription[cBSFDescriptionSize];
etSeqBase SeqBuff[cMaxFAMProbeLen];
int SeqLen;
bool bDescriptor;
bool bFirst;
int TruncCnt;

NumProbes = 0;
pProbes = NULL;
*ppProbes = NULL;

// open fasta file containing probes
if((pFasta = new CFasta()) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CFasta");
	return(eBSFerrObj);
	}
if((Rslt=pFasta->Open(pszProbesFile))!=eBSFSuccess)
	{
	while(pFasta->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,pFasta->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open input fasta file containing probes '%s'",pszProbesFile);
	delete pFasta;
	return(Rslt);
	}

// read each probe sequence into memory array
// note that only the first 10Knt of each probe is processed
bDescriptor = false;
bFirst = true;
TruncCnt = 0;
while((Rslt = SeqLen = pFasta->ReadSequence(SeqBuff,cMaxFAMProbeLen)) > eBSFSuccess)
	{
	if(SeqLen == eBSFFastaDescr)		// just read a descriptor line
		{
		NumProbes++;
		pFasta->ReadDescriptor(szDescription,cBSFDescriptionSize);
		TruncCnt = 0;
		bDescriptor = true;
		bFirst = true;					
		continue;
		}
	else
		if(!bDescriptor)					// if there was no descriptor then dummy up one...
			{
			if(!bFirst)						// only process 1st 10Knt
				{
				if(!TruncCnt++)
					gDiagnostics.DiagOut(eDLInfo,gszProcName,"LoadProbes [%s] truncated probe %s at %d length",pszProbesFile,szDescription, cMaxFAMProbeLen);
				continue;
				}
			NumProbes++;
			sprintf(szDescription,"Probe%d",NumProbes);
			TruncCnt = 0;
			bDescriptor = true;
			}
	bFirst = false;
	pNxtProbe = (sFAMProbe *)new unsigned char [ sizeof(sFAMProbe) + SeqLen];
	pNxtProbe->pNext = NULL;
	pNxtProbe->ProbID = NumProbes;
	pNxtProbe->ProbeLen = SeqLen;
	pNxtProbe->TotMinHits = 0;
	pNxtProbe->TotPlusHits = 0;
	strncpy(pNxtProbe->szDescr,szDescription,sizeof(pNxtProbe->szDescr)-1);
	pNxtProbe->szDescr[sizeof(pNxtProbe->szDescr)-1] = '\0';
	CSeqTrans::RemoveMasking(SeqBuff,SeqLen);
	memcpy(pNxtProbe->Bases,SeqBuff,SeqLen);
	if(pProbes == NULL)
		{
		pProbes = pCurProbe = pNxtProbe;
		}
	else
		{
		pCurProbe->pNext = pNxtProbe;
		pCurProbe = pNxtProbe;
		}
	bDescriptor = false;
	}
if(Rslt < eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessFastaFile [%s] %s",pFasta->ErrText((teBSFrsltCodes)Rslt),pFasta->GetErrMsg());
	return(Rslt);
	}
delete pFasta;
if(Rslt >= eBSFSuccess)
	{
	*ppProbes = pProbes;
	return(NumProbes);
	}
return(Rslt);
}

int 
Extract(int hRsltsFile,CFMIndex *pFMIndex,char *pszChrom,UINT32 position, UINT32 nchars) {

int error;
UINT32 NumBasesExtracted;
unsigned char *pBases;
char *pszAsciiBases;

error = pFMIndex->extract(position, position+nchars-1, &pBases, &NumBasesExtracted);
if(error >= 0 && NumBasesExtracted > 0)
	{
	pszAsciiBases = CSeqTrans::MapSeq2Ascii(pBases,NumBasesExtracted);
	fprintf(stdout, "Chrom: %s %lu to %lu\n   %s\n",pszChrom,position, position+nchars-1,pszAsciiBases);
	delete pBases;
	}
return(error >= 0 ? NumBasesExtracted : error);
}


int Display(int hRsltsFile,CFMIndex *pFMIndex,char *pszChrom,char Strand,unsigned char *Pattern, int PatternLen, UINT32 nFlankChars) {

	UINT32 numocc, i, *length_snippet, p, len;
	unsigned char *snippet_text;
	int error;
	char *pszAsciiBases;
	len = (UINT32)PatternLen + 2*nFlankChars;

	error=	pFMIndex->display (Pattern, (UINT32)PatternLen, nFlankChars, &numocc, 
				    	 &snippet_text, &length_snippet);
	if(error >= 0 && numocc > 0)
		{
		pszAsciiBases = CSeqTrans::MapSeq2Ascii(Pattern,PatternLen);
		fprintf (stdout, "Chrom: '%s' Strand: '%c' Pattern: '%s' # Snippets: %lu\n\n\n",pszChrom, Strand,pszAsciiBases, numocc);
		for (i = 0, p = 0; i < numocc; i++, p+=len) {
			pszAsciiBases = CSeqTrans::MapSeq2Ascii(&snippet_text[p],length_snippet[i]);
			fprintf (stdout, "length: %lu Bases: %s\n", length_snippet[i],pszAsciiBases);
			}
		}

	if(length_snippet != NULL)
		delete length_snippet;
	if(snippet_text != NULL)
		delete snippet_text;
	return(error >= 0 ? numocc : error);
}	


int ApproxMatch(bool bCounts,bool bGlobCounts,sFAMProbe *pProbe,int hRsltsFile,CFMIndex *pFMIndex,char *pszChrom,char Strand,char *pszDescr,unsigned char *pPattern, int PatternLen, UINT32 Identity) 
{
UINT32 numocc, i, *length_snippet, p, len;
UINT32 *occs;
unsigned char *snippet_text;
int error;
char *pszAsciiBases;
int MaxMismatches;
int Mismatches;
int WindowSize;
int FlankChars;
int Left;
etSeqBase *pPutativeSeq;
etSeqBase *pPutativeBase;
etSeqBase *pPatternBase;
int Idx;
int PutativeMismatches;
sFAMMatch *pMatches;
sFAMMatch *pNxtMatch;
sFAMMatch *pCurMatch;
int NumMatches;

pMatches = NULL;
pNxtMatch = NULL;
pCurMatch = NULL;
NumMatches = 0;

sFAMMatch **pSortMatches;
int NumUnique;
char szRsltsBuff[0x7fff];
int RsltsLen;

// determine max number of mismatches
if(Identity > 10)
	MaxMismatches = PatternLen - ((PatternLen * Identity)/100);
else
	MaxMismatches = Identity;
// assuming mismatches are uniformly distributed then determine the window size within at least one window 
// will never contain a mismatch

length_snippet = NULL;
snippet_text = NULL;
Left = 0;
Mismatches = MaxMismatches;
while(Left < PatternLen) {
	WindowSize = ((PatternLen - Left) + Mismatches) / (Mismatches+1);
	FlankChars = max(Left,PatternLen - (Left + WindowSize));
	error =	pFMIndex->display(&pPattern[Left],(UINT32)WindowSize, FlankChars, &numocc, 
			    	 &snippet_text, &length_snippet);
	if(error >= 0 && numocc > 0)
		{
		len = (UINT32)WindowSize + (2 * FlankChars);
		for (i = 0, p = 0; i < numocc; i++, p+=len) 
			{
			if(length_snippet[i] < len)
				continue;
			pPutativeSeq = &snippet_text[p + FlankChars - Left];
			pPutativeBase = pPutativeSeq;
			PutativeMismatches = MaxMismatches;
			pPatternBase = pPattern;
			for(Idx = 0; Idx < PatternLen && PutativeMismatches >= 0; Idx++)
				if(*pPutativeBase++ != *pPatternBase++)
					PutativeMismatches--;
			if(PutativeMismatches >= 0)
				{
				pNxtMatch = (sFAMMatch *)new unsigned char [sizeof(sFAMMatch) + PatternLen];
				pNxtMatch->MatchLen = PatternLen;
				pNxtMatch->MatchID = ++NumMatches;
				pNxtMatch->pNext = NULL;
				memcpy(pNxtMatch->Bases,pPutativeSeq,PatternLen);
				if(pMatches == NULL)
					pMatches = pCurMatch = pNxtMatch;
				else
					{
					pCurMatch->pNext = pNxtMatch;
					pCurMatch = pNxtMatch;
					}

				// we have struck lucky!!
	//			pszAsciiBases = CSeqTrans::MapSeq2Ascii(pPutativeSeq,PatternLen);
	//			fprintf (stdout, "Snippet length: %lu Bases: %s\n", PatternLen,pszAsciiBases);
				}
			}
		}
	Mismatches--;
	Left += WindowSize;
	if(length_snippet != NULL)
		{
		delete length_snippet;
		length_snippet = NULL;
		}
	if(snippet_text != NULL)
		{
		delete snippet_text;
		snippet_text = NULL;
		}
	}
if(length_snippet != NULL)
	delete length_snippet;
if(snippet_text != NULL)
	delete snippet_text;
if(!NumMatches || error < 0)
	return(error >= 0 ? 0 : error);

pCurMatch = pMatches;
pSortMatches = (sFAMMatch **)new sFAMMatch *[NumMatches];
for(Idx = 0; Idx < NumMatches; Idx++)
	{
	pSortMatches[Idx] = pCurMatch;
	pCurMatch = pCurMatch->pNext;
	}

NumUnique = 1;
if(NumMatches > 1)
	{
	qsort(pSortMatches,NumMatches,sizeof(sFAMMatch *),SortMatches);

	// now that matches are sorted then can remove any duplicates
	pNxtMatch = pSortMatches[0];
	for(Idx=1;Idx < NumMatches; Idx++)
		{
		pCurMatch = pSortMatches[Idx];
		if(!memcmp(pCurMatch->Bases,pNxtMatch->Bases,pCurMatch->MatchLen))
			continue;
		pSortMatches[NumUnique++] = pCurMatch;
		pNxtMatch = pCurMatch;
		}
	}

// deduped, can now get their loci or just output the counts
if(bCounts || bGlobCounts)
	{
	if(!bGlobCounts)
		{
		RsltsLen = sprintf(szRsltsBuff,"%d,\"%s\",\"%s\",\"%c\",%lu\n",++gMatchID,pszDescr,pszChrom,Strand,NumUnique);
		write(hRsltsFile,szRsltsBuff,RsltsLen);
		}
	if(Strand == '-')
		pProbe->TotMinHits += numocc;
	else
		pProbe->TotPlusHits += numocc;
	}
else
	{
	for(Idx = 0; Idx < NumUnique; Idx++)
		{
		pCurMatch = pSortMatches[Idx];
		error = pFMIndex->locate (pCurMatch->Bases, (UINT32)pCurMatch->MatchLen, &occs, &numocc);
		if(error >= 0 && numocc > 0)
			{
			pszAsciiBases = CSeqTrans::MapSeq2Ascii(pCurMatch->Bases,pCurMatch->MatchLen);
			for (i = 0; i < numocc; i++)
				{
				RsltsLen = sprintf(szRsltsBuff,"%d,\"%s\",\"%s\",\"%c\",%lu,%d,\"%s\"\n",++gMatchID,pszDescr,pszChrom,Strand,occs[i],pCurMatch->MatchLen,pszAsciiBases);
				write(hRsltsFile,szRsltsBuff,RsltsLen);
				}
			}
		if(occs != NULL) 
			delete occs;
		}
	}

// cleanup
delete pSortMatches;
while((pCurMatch = pMatches)!=NULL)
	{
	pMatches = pCurMatch->pNext;
	delete pCurMatch;
	}
return(error >= 0 ? numocc : error);
}	


int Locate(bool bCounts,bool bGlobCounts,sFAMProbe *pProbe,int hRsltsFile,CFMIndex *pFMIndex,char *pszChrom,char Strand,char *pszDescr,unsigned char *Pattern,int PatternLen) 
{ 
char szRsltsBuff[0x7fff];
int RsltsLen;
UINT32 *occs, numocc, Idx;
int error;
char *pszAsciiBases;

if(bCounts || bGlobCounts)
	return(Count(bGlobCounts,pProbe,hRsltsFile,pFMIndex,pszChrom,Strand,pszDescr,Pattern,PatternLen));
error = pFMIndex->locate (Pattern, (UINT32)PatternLen, &occs, &numocc);
if(error >= 0 && numocc > 0)
	{
	pszAsciiBases = CSeqTrans::MapSeq2Ascii(Pattern,PatternLen);
	for (Idx = 0; Idx < numocc; Idx++)
		{
		RsltsLen = sprintf(szRsltsBuff,"%d,\"%s\",\"%s\",\"%c\",%lu,%d,\"%s\"\n",++gMatchID,pszDescr,pszChrom,Strand,occs[Idx],PatternLen,pszAsciiBases);
		write(hRsltsFile,szRsltsBuff,RsltsLen);
		}
	}
if(occs != NULL) 
	delete occs;
return(error >= 0 ? numocc : error);
}	


int Count(bool bGlobCounts,sFAMProbe *pProbe,int hRsltsFile,CFMIndex *pFMIndex,char *pszChrom,char Strand,char *pszDescr,unsigned char *Pattern,int PatternLen) 
{
	char szRsltsBuff[0x7fff];
	UINT32 numocc;
	int error;
	int RsltsLen;

	error = pFMIndex->count (Pattern, (UINT32)PatternLen, &numocc);
	if(error >= 0 && numocc > 0)
		{
		if(!bGlobCounts)
			{
			RsltsLen = sprintf(szRsltsBuff,"%d,\"%s\",\"%s\",\"%c\",%lu\n",++gMatchID,pszDescr,pszChrom,Strand,numocc);
			write(hRsltsFile,szRsltsBuff,RsltsLen);
			}
		if(Strand == '-')
			pProbe->TotMinHits += numocc;
		else
			pProbe->TotPlusHits += numocc;
		}
	return(error >= 0 ? numocc : error);
}

int 
SortMatches(const void *arg1, const void *arg2)
{
sFAMMatch *pEl1 = *(sFAMMatch **)arg1;
sFAMMatch *pEl2 = *(sFAMMatch **)arg2;
return(memcmp(pEl1->Bases,pEl2->Bases,pEl1->MatchLen));
}




