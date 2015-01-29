// kangas.cpp : Defines the entry point for the console application.
// Generates biosequence files

#include "stdafx.h"

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#if _WIN32
#include "../libbiokanga/commhdrs.h"
#else
#include "../libbiokanga/commhdrs.h"
#endif

CStopWatch gStopWatch;
CDiagnostics gDiagnostics;				// for writing diagnostics messages to log file
char gszProcName[_MAX_FNAME];			// process name

const char *cpszProgVer = "1.12.0";		// increment with each release

const int cMaxAllocBuffChunk = 0x00ffffff;	// buffer for fasta sequences is realloc'd in this sized chunks


typedef struct TAG_sProcParams 
	{
	bool bCapsSoftMask;		// treat uppercase as softmasked bases
	CBioSeqFile *pSeqFile;	// preopened bioseq file
	} tsProcParams; 


int
CreateBioseqFastaFile(bool bTargDeps,			// true if process only if any independent src files newer than target
					  char *pszSrcDirPath,char *pszDestBioseqFile,char *pszRefSpecies,
						char *pszDescr,char *pszTitle,
						bool bCapsSoftMask);

int ValidateFastaFile(char *pszFastaFile,
					   char *pszBioSeqFile,
						unsigned int EntryID);			// initial entry identifier


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

bool bCapsSoftmask;

char szOutputFileSpec[_MAX_PATH];
char szInputFileSpec[_MAX_PATH];
char szDescription[cMBSFFileDescrLen];
char szTitle[cMBSFShortFileDescrLen];
char szRefSpecies[80];

// command line args
struct arg_lit  *help    = arg_lit0("h","help",                 "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to screen and logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"diagnostics log file");

struct arg_file *InFile = arg_file1("i",NULL,"<file>",			"input from multifasta files");
struct arg_file *OutFile = arg_file1("o",NULL,"<file>",			"output as bioseq file");
struct arg_str *RefSpecies = arg_str0("r","ref","<string>",		"reference species ");
struct arg_str *Descr = arg_str0("d","descr","<string>",		"full description");
struct arg_str *Title = arg_str0("t","title","<string>",		"short title");
struct arg_lit  *CapsSoftmask    = arg_lit0("c","capsmask",     "caps represent softmasked bases (default is lowercase)");
struct arg_lit  *TargDeps = arg_lit0("D","TargDep",				"Generate target file only if missing or older than any of the independent source files");

struct arg_end *end = arg_end(20);

void *argtable[] = {help,version,FileLogLevel,LogFile,InFile,OutFile,RefSpecies,Descr,Title,CapsSoftmask,
					TargDeps,end};

char **pAllArgs;
int argerrors;
argerrors = CUtility::arg_parsefromfile(argc,(char **)argv,&pAllArgs);
if(argerrors >= 0)
	argerrors = arg_parse(argerrors,pAllArgs,argtable);

/* special case: '--help' takes precedence over error reporting */
if (help->count > 0)
        {
		printf("\n%s Kanga Biosequence File Generator, Version %s\nOptions ---\n", gszProcName,cpszProgVer);
        arg_print_syntax(stdout,argtable,"\n");
        arg_print_glossary(stdout,argtable,"  %-25s %s\n");
		printf("\nNote: Parameters can be entered into a parameter file, one parameter per line.");
		printf("\n      To invoke this parameter file then precede it's name with '@'");
		printf("\n      e.g. %s @myparams.txt\n",gszProcName);
		printf("\nPlease report any issues regarding usage of %s to stuart.stephen@csiro.au\n\n",gszProcName);
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

	bCapsSoftmask = CapsSoftmask->count ? true : false;
	if(InFile->count)
		strcpy(szInputFileSpec,InFile->filename[0]);
	else
		strcpy(szInputFileSpec,"in.fa");


	if(OutFile->count)
		strcpy(szOutputFileSpec,OutFile->filename[0]);
	else
		strcpy(szOutputFileSpec,"out.seq");

	if(RefSpecies->count)
		{
		strncpy(szRefSpecies,RefSpecies->sval[0],cMBSFShortFileDescrLen);
		szRefSpecies[cMBSFShortFileDescrLen-1] = '\0';
		}
	else
		strcpy(szRefSpecies,"Dataset1");
	strcpy(szTitle,szRefSpecies);

	if(Title->count)
		{
		strncpy(szTitle,Title->sval[0],cMBSFShortFileDescrLen);
		szTitle[cMBSFShortFileDescrLen-1] = '\0';
		}
	else
		strcpy(szTitle,szRefSpecies);

	if(Descr->count)
		{
		strncpy(szDescription,Descr->sval[0],cMBSFFileDescrLen);
		szDescription[cMBSFFileDescrLen-1] = '\0';
		}
	else
		strcpy(szDescription,szTitle);

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

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Version: %s Processing parameters:",cpszProgVer);
	if(bTargDeps)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Generate target file only if missing or older than any of the independent source files");
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Processing Mode: %s","Normal");
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Input from multifasta files: %s",szInputFileSpec);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Output as bioseq file: '%s'",szOutputFileSpec);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Species: '%s'",szRefSpecies);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Title text: '%s'",szTitle);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Descriptive text: '%s'",szDescription);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Caps represent softmasked bases: '%s'",bCapsSoftmask ? "yes" : "no");

#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	gStopWatch.Start();
	Rslt = CreateBioseqFastaFile(bTargDeps,szInputFileSpec,szOutputFileSpec,szRefSpecies,szDescription,szTitle,
						bCapsSoftmask);
	gStopWatch.Stop();
	Rslt = Rslt >=0 ? 0 : 1;
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Exit code: %d Total processing time: %s",Rslt,gStopWatch.Read());
	exit(Rslt);
	}
else
	{
	printf("\n%s Kanga Biosequence File Generator, Version %s\n",gszProcName,cpszProgVer);
	arg_print_errors(stdout,end,gszProcName);
	arg_print_syntax(stdout,argtable,"\nUse '-h' to view option and parameter usage\n");
	exit(1);
	}
}

// ValidateFastaFile
// Validates that bioseq file starting with specified entryid contains sequences in specified fasta file
int 
ValidateFastaFile(char *pszFastaFile,
					   char *pszBioSeqFile,
						unsigned int EntryID)			// initial entry identifier
{
int Rslt;
CFasta *pFasta;
CBioSeqFile *pSeqFile;
unsigned char *pFastaSeqBuff;
unsigned char *pBioSeqBuff;
unsigned int BioEntryDataLen;
unsigned int CurSeqOfs;
unsigned int TestBuffSize;

int FastaSeqLen;
int BioSeqLen;

if((pFasta = new CFasta)==NULL)
	return(-1);

if((Rslt=pFasta->Open(pszFastaFile,true))!=eBSFSuccess)
	{
	while(pFasta->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,pFasta->GetErrMsg());
	delete pFasta;
	return(Rslt);
	}

if((pSeqFile = new CBioSeqFile)==NULL)
	{
	delete pFasta;
	return(eBSFerrObj);
	}

if((Rslt=pSeqFile->Open(pszBioSeqFile,eSeqBaseType))!=eBSFSuccess)
	{
	while(pSeqFile->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,pSeqFile->GetErrMsg());
	pFasta->Close();
	delete pFasta;
	delete pSeqFile;
	return(Rslt);
	}

if((Rslt=pSeqFile->Exists(EntryID))!= eBSFSuccess)
	{
	pFasta->Close();
	pSeqFile->Close();
	delete pFasta;
	delete pSeqFile;
	return(Rslt);
	}


if((pFastaSeqBuff = new unsigned char [cMaxAllocBuffChunk])==NULL)
	{
	pFasta->Close();
	pSeqFile->Close();
	delete pFasta;
	delete pSeqFile;
	return(eBSFerrMem);
	}

if((pBioSeqBuff = new unsigned char [cMaxAllocBuffChunk])==NULL)
	{
	delete pFastaSeqBuff;
	pFasta->Close();
	pSeqFile->Close();
	delete pFasta;
	delete pSeqFile;
	return(eBSFerrMem);
	}

BioEntryDataLen = pSeqFile->GetDataLen(EntryID);
CurSeqOfs = 0;
TestBuffSize = 12345;
while((FastaSeqLen = pFasta->ReadSequence(pFastaSeqBuff,TestBuffSize,true)) > eBSFSuccess)
	{
	if(FastaSeqLen == eBSFFastaDescr)		// just read a descriptor line
		continue;
	// see what the bioseq looks like - must be identical!
	BioSeqLen = pSeqFile->GetData(EntryID,eSeqBaseType,CurSeqOfs,pBioSeqBuff,FastaSeqLen);
	if(BioSeqLen != FastaSeqLen)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Sequence length error!");
		continue;
		}
	if(memcmp(pFastaSeqBuff,pBioSeqBuff,BioSeqLen))
		{
		for(int Idx = 0; Idx < BioSeqLen; Idx++)
			if(pFastaSeqBuff[Idx] != pBioSeqBuff[Idx])
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Mismatch error at %d!",Idx);
				break;
				}
		}
	CurSeqOfs += BioSeqLen;
	BioEntryDataLen -= BioSeqLen;
	// make buffsize a pseudorandom value between 1 and (cMaxAllocBuffChunk -1)
	TestBuffSize *= 125031;	// guess this must be a prime!
	TestBuffSize %= cMaxAllocBuffChunk;
	if(!TestBuffSize)
		TestBuffSize = 1;
	}
delete pFastaSeqBuff;
delete pBioSeqBuff;
pFasta->Close();
pSeqFile->Close();
delete pFasta;
delete pSeqFile;
gDiagnostics.DiagOut(eDLInfo,gszProcName,"No errors!");
return(eBSFSuccess);
}


// ProcessFastaFile
// Parse input fasta format file into a biosequence file
bool ProcessFastaFile(char *pszFile,
				 void *pParams)			
{
CFasta Fasta;
unsigned char *pSeqBuff;
char szName[cBSFSourceSize];
char szDescription[cBSFDescriptionSize];
int SeqLen;
int Descrlen;
bool bFirstEntry;
bool bEntryCreated;
tBSFEntryID EntryID;
int Rslt;
int SeqID;
bool bCapsSoftMask;

if((pSeqBuff = new unsigned char [cMaxAllocBuffChunk]) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessFastaFile:- Unable to allocate memory (%d bytes) for sequence buffer",cMaxAllocBuffChunk);
	return(false);
	}

CBioSeqFile *pSeqFile = ((tsProcParams *)pParams)->pSeqFile;
bCapsSoftMask = ((tsProcParams *)pParams)->bCapsSoftMask;

gDiagnostics.DiagOut(eDLInfo,gszProcName,"ProcessFastaFile:- Adding %s..",pszFile);
if((Rslt=Fasta.Open(pszFile,true))!=eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open '%s' [%s] %s",pszFile,Fasta.ErrText((teBSFrsltCodes)Rslt),Fasta.GetErrMsg());
	delete pSeqBuff;
	return(false);
	}
bFirstEntry = true;
bEntryCreated = false;
SeqID = 0;
while((Rslt = SeqLen = Fasta.ReadSequence(pSeqBuff,cMaxAllocBuffChunk,true,bCapsSoftMask)) > eBSFSuccess)
	{
	if(SeqLen == eBSFFastaDescr)		// just read a descriptor line
		{
		SeqID++;
		if(bEntryCreated)	// close any previous entry
			pSeqFile->SealEntry();
		Descrlen = Fasta.ReadDescriptor(szDescription,cBSFDescriptionSize);
		// An assumption - will one day bite real hard - is that the
		// fasta descriptor line starts with some form of unique identifier.
		// Use this identifier as the entry name.
		if(sscanf(szDescription," %s[ ,]",szName)!=1)
			sprintf(szName,"%s.%d",pszFile,++SeqID);
		EntryID = pSeqFile->CreateEntry(szName, szDescription,eSeqBaseType,bCapsSoftMask);
		if(EntryID < eBSFSuccess)
			{
			delete pSeqBuff;
			return(false);
			}
		bFirstEntry = false;
		bEntryCreated = true;
		continue;
		}
	else
		if(bFirstEntry)	// if there was no descriptor then dummy up one...
			{
			SeqID++;
			sprintf(szName,"%s.%d",pszFile,SeqID);
			EntryID =  pSeqFile->CreateEntry(szName, (char *)"No Description provided",eSeqBaseType,bCapsSoftMask);
			if(EntryID < eBSFSuccess)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unble to create entry '%s' [%s] %s",szName,Fasta.ErrText((teBSFrsltCodes)EntryID),Fasta.GetErrMsg());
				delete pSeqBuff;
				return(false);
				}
			bFirstEntry = false;
			bEntryCreated = true;
			}
	if((Rslt=pSeqFile->AddData(SeqLen,pSeqBuff)) < eBSFSuccess) // save sequence
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to add data [%s] %s",Fasta.ErrText((teBSFrsltCodes)Rslt),Fasta.GetErrMsg());
		delete pSeqBuff;
		return(false);
		}
	}
if(Rslt < eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessFastaFile [%s] %s",Fasta.ErrText((teBSFrsltCodes)Rslt),Fasta.GetErrMsg());
	delete pSeqBuff;
	return(false);
	}
if(bEntryCreated)
	pSeqFile->SealEntry();
delete pSeqBuff;
return(true);
}



int
CreateBioseqFastaFile(bool bTargDeps,			// true if process only if any independent src files newer than target
					  char *pszSrcDirPath,char *pszDestBioseqFile,char *pszRefSpecies,
						char *pszDescr,char *pszTitle,
						bool bCapsSoftMask)
{
int Rslt;
tsProcParams ProcParams;
bool bCreate;

CSimpleGlob glob(SG_GLOB_FULLSORT);
if (glob.Add(pszSrcDirPath) >= SG_SUCCESS)
	{
	bCreate = bTargDeps ? false : true;
	for (int n = 0; n < glob.FileCount(); ++n)
		{
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Will process in this order: %d '%s", n+1,glob.File(n));
		if(!bCreate && CUtility::ChkTargDepend(NULL,0,pszDestBioseqFile,glob.File(n),NULL)!= 0)
			bCreate = true;
		}
	}
else
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to glob '%s",pszSrcDirPath);
	return(eBSFerrOpnFile);	
    }
if(!bCreate)
	{
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Target '%s' is already updated",pszDestBioseqFile);
	return(eBSFSuccess);
	}

CBioSeqFile *pSeqFile = new CBioSeqFile;
if(pSeqFile == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"CreateBioseqFastaFile: Unable to instantiate CBioSeqFile");
	return(eBSFerrObj);
	}
if((Rslt=pSeqFile->Open(pszDestBioseqFile,cBSFTypeSeq,true)) != eBSFSuccess)
	{
	while(pSeqFile->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,pSeqFile->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"CreateBioseqFastaFile: Unable to create/truncate %s",pszDestBioseqFile);
	delete pSeqFile;
	return(Rslt);
	}
if((Rslt = pSeqFile->SetDescription(pszDescr)) != eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"CreateBioseqFastaFile: Unable to set description '%s' into %s",pszDescr,pszDestBioseqFile);
	delete pSeqFile;
	return(Rslt);
	}
if((Rslt=pSeqFile->SetTitle(pszTitle)) != eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"CreateBioseqFastaFile: Unable to set title '%s' into %s",pszTitle,pszDestBioseqFile);
	delete pSeqFile;
	return(Rslt);
	}

if((Rslt=pSeqFile->SetDatasetName(pszRefSpecies)) != eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"CreateBioseqFastaFile: Unable to set dataset %s",pszRefSpecies);
	delete pSeqFile;
	return(Rslt);
	}
ProcParams.bCapsSoftMask = bCapsSoftMask;
ProcParams.pSeqFile = pSeqFile;

srand(0);	// rand() used when processing for FMIndex generation to randomise long sequences of eBaseN otherwise ds_sort takes too long

for (int n = 0; Rslt >= eBSFSuccess &&  n < glob.FileCount(); ++n)
	Rslt=ProcessFastaFile(glob.File(n),&ProcParams);

pSeqFile->Close();
delete pSeqFile;
return(Rslt);
}






