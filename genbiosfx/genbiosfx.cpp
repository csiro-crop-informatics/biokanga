// genbiosfx.cpp : Defines the entry point for the console application.
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

const unsigned int cProgVer = 400;		// increment with each release
const int cMaxAllocBuffChunk = 0x00ffffff;	// buffer for fasta sequences is realloc'd in this sized chunks

CStopWatch gStopWatch;
CDiagnostics gDiagnostics;				// for writing diagnostics messages to log file
char gszProcName[_MAX_FNAME];			// process name


int
CreateBioseqSuffixFile(bool bTargDeps,			// true if process only if any independent src files newer than target
					   int iMode,char *pszSrcDirPath,char *pszDestSfxFile,char *pszRefSpecies,char *pszDescr,char *pszTitle);


CSfxArrayV2 *m_pSfxFile;				// suffix array file being created

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
bool bTargDeps;						// true if only process if independent files newer than target

char szOutputFileSpec[_MAX_PATH];
char szInputFileSpec[_MAX_PATH];
char szDescription[cMBSFFileDescrLen];
char szTitle[cMBSFShortFileDescrLen];
char szRefSpecies[cMaxDatasetSpeciesChrom];
int iMode;									// processing mode

// command line args
struct arg_lit  *help    = arg_lit0("hH","help",                "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_int *ScreenLogLevel=arg_int0("S", "ScreenLogLevel",	"<int>","Level of diagnostics written to screen 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"diagnostics log file");

struct arg_int *Mode=arg_int0("m", "mode",	"<int>",			"Input file type, 0=bioseq, 1=fasta");
struct arg_file *InFile = arg_file1("i",NULL,"<file>",			"input from wildcarded bioseq or fasta files");
struct arg_file *OutFile = arg_file1("o",NULL,"<file>",			"output as suffix array V2 file");
struct arg_str *Descr = arg_str0("d","descr","<string>",		"full description");
struct arg_str *Title = arg_str0("t","title","<string>",		"short title");
struct arg_str *RefSpecies = arg_str1("r","ref","<string>",		"reference species");
struct arg_lit  *TargDeps = arg_lit0("D","TargDep",				"Generate target file only if missing or older than any of the independent source files");

struct arg_end *end = arg_end(20);

void *argtable[] = {help,version,FileLogLevel,ScreenLogLevel,LogFile,Mode,InFile,OutFile,RefSpecies,Descr,Title,TargDeps,end};

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
		printf("%s Version: %d.%2.2d\n",gszProcName,cProgVer/100,cProgVer%100);
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

	iMode = Mode->count ? Mode->ival[0] : 0;
	if(iMode < 0)
		iMode = 0;
	else
		if(iMode > 1)
			iMode = 1;

	if(InFile->count)
		strcpy(szInputFileSpec,InFile->filename[0]);
	else
		strcpy(szInputFileSpec,"in.fa");


	if(OutFile->count)
		strcpy(szOutputFileSpec,OutFile->filename[0]);
	else
		strcpy(szOutputFileSpec,"out.seq");

	strncpy(szRefSpecies,RefSpecies->sval[0],cMaxDatasetSpeciesChrom);
	szRefSpecies[cMaxDatasetSpeciesChrom-1] = '\0';

	if(!Title->count)
		strcpy(szTitle,szRefSpecies);
	else
		{
		strncpy(szTitle,Title->sval[0],cMBSFShortFileDescrLen);
		szTitle[cMBSFShortFileDescrLen-1] = '\0';
		}

	if(!Descr->count)
		strcpy(szDescription,szRefSpecies);
	else
		{
		strncpy(szDescription,Descr->sval[0],cMBSFFileDescrLen);
		szDescription[cMBSFFileDescrLen-1] = '\0';
		}

	if (TargDeps->count > 0)
		bTargDeps = true;
	else
		bTargDeps = false;


			// now that command parameters have been parsed then initialise diagnostics log system
	if(!gDiagnostics.Open(szLogFile,(etDiagLevel)iScreenLogLevel,(etDiagLevel)iFileLogLevel,true))
		{
		printf("\nError: Unable to start diagnostics subsystem.");
		if(szLogFile[0] != '\0')
			printf(" Most likely cause is that logfile '%s' can't be opened/created",szLogFile);
		exit(1);
		}

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Version: %d.%2.2d Processing parameters:",cProgVer/100,cProgVer%100);
	if(bTargDeps)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Generate target file only if missing or older than any of the independent source files");
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Process input files as: '%s'",iMode == 0 ? "bioseq" : "fasta");
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Input source file spec: '%s'",szInputFileSpec);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Output to suffix array file: '%s'",szOutputFileSpec);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Reference species: '%s'",szRefSpecies);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Title text: '%s'",szTitle);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Descriptive text: '%s'",szDescription);
#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	gStopWatch.Start();
	Rslt = CreateBioseqSuffixFile(bTargDeps,iMode,szInputFileSpec,szOutputFileSpec,szRefSpecies,szDescription,szTitle);
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


void Init(void)
{
m_pSfxFile = NULL;
}

void Reset(void)
{
if(m_pSfxFile != NULL)
	{
	delete m_pSfxFile;
	m_pSfxFile = NULL;
	}
}

// ProcessBioseqFile
// Process input biosequence file into suffix file
bool ProcessBioseqFile(char *pszFile)
{
CBioSeqFile BioSeqFile;
etSeqBase *pSeqBuff = NULL;
int AllocLen = 0;
char szSource[cBSFSourceSize];
char szDescription[cBSFDescriptionSize];
int SeqLen;
int Rslt;
tBSFEntryID CurEntryID;

if((Rslt=BioSeqFile.Open(pszFile,cBSFTypeSeq,false))!=eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessBioseqFile - error %d %s",Rslt,BioSeqFile.GetErrMsg());
	return(false);
	}

CurEntryID = 0;
while((Rslt = CurEntryID = BioSeqFile.Next(CurEntryID)) > eBSFSuccess)
	{
	BioSeqFile.GetNameDescription(CurEntryID,cBSFSourceSize-1,(char *)&szSource,
											cBSFDescriptionSize-1,(char *)&szDescription);
	SeqLen = BioSeqFile.GetDataLen(CurEntryID);
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing %s|%s",szSource,szDescription);

	if(!SeqLen)
		continue;

	if(AllocLen < (SeqLen + 1))
		{
		if(pSeqBuff != NULL)
			delete pSeqBuff;
		AllocLen = SeqLen;
		pSeqBuff = new unsigned char [AllocLen];
		if(pSeqBuff == NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessBioseqFile - Unable to alloc %d bytes memory for pSeqBuff",AllocLen);
			Rslt = eBSFerrMem;
			break;
			}
		}

	if((Rslt = BioSeqFile.GetData(CurEntryID,eSeqBaseType,0,pSeqBuff,SeqLen)) != SeqLen)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessBioseqFile - error %d %s",Rslt,BioSeqFile.GetErrMsg());
		break;
		}



	// remove any repeat masking flags so that sorts can actually sort
	// if run of more than 25 Ns and at least 5 Ns to end of buffer then randomly mutate
	// every 13th N
	//	e.g <25Ns>r<12Ns>r<12Ns> where r is a pseudorandom base 
	etSeqBase *pMskBase = pSeqBuff;
	int SeqNs = 0;
	for(int MskIdx = 0; MskIdx < SeqLen; MskIdx++,pMskBase++)
		{
		*pMskBase &= ~cRptMskFlg;
		if(*pMskBase == eBaseN && (MskIdx+5) < SeqLen)
			{
			if(++SeqNs > 25 && 
				pMskBase[1] == eBaseN &&
				pMskBase[2] == eBaseN &&
				pMskBase[3] == eBaseN &&
				pMskBase[4] == eBaseN)
				{
				if(!(SeqNs % 13))	// mutate every 13th
					*pMskBase = rand() % 4;
				}
			}
		else
			SeqNs = 0;
		}



	if((Rslt=m_pSfxFile->AddEntry(szSource,pSeqBuff,SeqLen))<= 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessBioseqFile - error %d %s",Rslt,m_pSfxFile->GetErrMsg());
		break;
		}
	}
if(pSeqBuff != NULL)
	delete pSeqBuff;
BioSeqFile.Close();
return(Rslt ==  eBSFSuccess ? true : false);
}

// ProcessFastaFile
// Parse input fasta format file into a biosequence suffix array file
bool ProcessFastaFile(char *pszFile)
{
CFasta Fasta;
unsigned char *pSeqBuff;
unsigned char *pMskBase;
int MskIdx;
int BuffOfs;
int AllocdBuffSize;
int AvailBuffSize;
char szName[cBSFSourceSize];
char szDescription[cBSFDescriptionSize];
int SeqLen;
int Descrlen;
bool bFirstEntry;
bool bEntryCreated;
int Rslt;
int SeqID;

AllocdBuffSize = cMaxAllocBuffChunk * 16;
// note malloc is used as can then simply realloc to expand as may later be required
if((pSeqBuff = (unsigned char *)malloc(AllocdBuffSize)) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessFastaFile:- Unable to allocate memory (%d bytes) for sequence buffer",AllocdBuffSize);
	return(false);
	}
AvailBuffSize = AllocdBuffSize;

gDiagnostics.DiagOut(eDLInfo,gszProcName,"ProcessFastaFile:- Adding %s..",pszFile);
if((Rslt=Fasta.Open(pszFile,true))!=eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessFastaFile: Unable to open '%s' [%s] %s",pszFile,Fasta.ErrText((teBSFrsltCodes)Rslt),Fasta.GetErrMsg());
	free(pSeqBuff);
	return(false);
	}
bFirstEntry = true;
bEntryCreated = false;
SeqID = 0;
BuffOfs = 0;
while((Rslt = SeqLen = Fasta.ReadSequence(&pSeqBuff[BuffOfs],AvailBuffSize,true,false)) > eBSFSuccess)
	{
	if(SeqLen == eBSFFastaDescr)		// just read a descriptor line
		{
		SeqID++;
		if(bEntryCreated)				// add any previous entry
			{
			if((Rslt=m_pSfxFile->AddEntry(szName,pSeqBuff,BuffOfs)) <= 0)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessFastaFile - error %d %s",Rslt,m_pSfxFile->GetErrMsg());
				break;
				}
			}
		Descrlen = Fasta.ReadDescriptor(szDescription,cBSFDescriptionSize);
		// An assumption - will one day bite real hard - is that the
		// fasta descriptor line starts with some form of unique identifier.
		// Use this identifier as the entry name.
		if(sscanf(szDescription," %s[ ,]",szName)!=1)
			sprintf(szName,"%s.%d",pszFile,++SeqID);

		bFirstEntry = false;
		bEntryCreated = true;
		BuffOfs = 0;
		continue;
		}
	else
		if(bFirstEntry)	// if there was no descriptor then dummy up one...
			{
			SeqID++;
			sprintf(szName,"%s.%d",pszFile,SeqID);
			strcpy(szDescription,"No Description provided");
			bFirstEntry = false;
			bEntryCreated = true;
			}

	// remove any repeat masking flags so that sorts can actually sort
	// if run of more than 25 Ns and at least 5 Ns to end of buffer then randomly mutate
	// every 13th N
	//	e.g <25Ns>r<12Ns>r<12Ns> where r is a pseudorandom base 
	pMskBase = &pSeqBuff[BuffOfs];
	int SeqNs = 0;
	for(MskIdx = 0; MskIdx < SeqLen; MskIdx++,pMskBase++)
		{
		*pMskBase &= ~cRptMskFlg;
		if(*pMskBase == eBaseN && (MskIdx+5) < SeqLen)
			{
			if(++SeqNs > 25 && 
				pMskBase[1] == eBaseN &&
				pMskBase[2] == eBaseN &&
				pMskBase[3] == eBaseN &&
				pMskBase[4] == eBaseN)
				{
				if(!(SeqNs % 13))	// mutate every 13th
					*pMskBase = rand() % 4;
				}
			}
		else
			SeqNs = 0;
		}

	BuffOfs += SeqLen;
	AvailBuffSize -= SeqLen;
	if(AvailBuffSize < (cMaxAllocBuffChunk / 8))
		{
		int NewSize = cMaxAllocBuffChunk + AllocdBuffSize;
		unsigned char *pTmp;
		if((pTmp = (unsigned char *)realloc(pSeqBuff,NewSize))==NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessFastaFile:- Unable to reallocate memory (%d bytes) for sequence buffer",NewSize);
			return(false);
			}
		pSeqBuff = pTmp;
		AllocdBuffSize = NewSize;
		AvailBuffSize = AllocdBuffSize - BuffOfs;
		}
	}

if(Rslt >= eBSFSuccess && bEntryCreated && BuffOfs > 0)			// close entry
	if((Rslt=m_pSfxFile->AddEntry(szName,pSeqBuff,BuffOfs)) <= 0)
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessFastaFile - error %d %s",Rslt,m_pSfxFile->GetErrMsg());
	else
		Rslt = eBSFSuccess;

if(pSeqBuff != NULL)
	delete pSeqBuff;
return(Rslt >=  eBSFSuccess ? true : false);
}



int
CreateBioseqSuffixFile(bool bTargDeps,			// true if process only if any independent src files newer than target
					   int Mode,char *pszSrcDirPath,char *pszDestSfxFile,char *pszRefSpecies,char *pszDescr,char *pszTitle)
{
int Rslt;
bool bCreate;

Init();

CSimpleGlob glob(SG_GLOB_FULLSORT);
if (glob.Add(pszSrcDirPath) >= SG_SUCCESS)
	{
	bCreate = bTargDeps ? false : true;
	for (int n = 0; n < glob.FileCount(); ++n)
		{
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Will process in this order: %d '%s", n+1,glob.File(n));
		if(!bCreate && CUtility::ChkTargDepend(NULL,0,pszDestSfxFile,glob.File(n),NULL)!= 0)
			bCreate = true;
		}
	}
else
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to glob '%s",pszSrcDirPath);
	Reset();
	return(eBSFerrOpnFile);	
    }
if(!bCreate)
	{
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Target '%s' is already updated",pszDestSfxFile);
	Reset();
	return(eBSFSuccess);
	}

m_pSfxFile = new CSfxArrayV2;
if(m_pSfxFile == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"CreateBioseqSuffixFile - error Unable to create instance of CSfxArray");
	Reset();
	return(eBSFerrMem);
	}
if((Rslt=m_pSfxFile->Open(pszDestSfxFile,true))!=eBSFSuccess)
	{
	while(m_pSfxFile->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,m_pSfxFile->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"CreateBioseqSuffixFile, unable to create '%s' - error %d %s",pszDestSfxFile,Rslt,m_pSfxFile->GetErrMsg());
	Reset();
	return(Rslt);
	}

if((Rslt=m_pSfxFile->SetDescription(pszDescr)) != eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"CreateBioseqSuffixFile: Unable to set description '%s' into %s",pszDescr,pszDestSfxFile);
	Reset();
	return(Rslt);
	}
if((Rslt=m_pSfxFile->SetTitle(pszTitle)) != eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"CreateBioseqSuffixFile: Unable to set title '%s' into %s",pszTitle,pszDestSfxFile);
	Reset();
	return(Rslt);
	}

if((Rslt = m_pSfxFile->SetDatasetName(pszRefSpecies)) != eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"CreateBioseqSuffixFile: Unable to set dataset %s",pszRefSpecies);
	Reset();
	return(Rslt);
	}

Rslt = eBSFSuccess;
for (int n = 0; Rslt >= eBSFSuccess &&  n < glob.FileCount(); ++n)
	{
	if(Mode == 0)
		Rslt = ProcessBioseqFile(glob.File(n));
	else
		Rslt = ProcessFastaFile(glob.File(n));
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"CreateBioseqSuffixFile: finalising...");
m_pSfxFile->Finalise();
gDiagnostics.DiagOut(eDLInfo,gszProcName,"CreateBioseqSuffixFile: completed...");
m_pSfxFile->Reset(false);
Reset();
return(Rslt);
}


