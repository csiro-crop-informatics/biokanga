// kangarg.cpp : Defines the entry point for the console application
// Kanga genome assembly kmer randomiser
//
#include "stdafx.h"

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#if _WIN32
#include <process.h>
#include "../libbiokanga/commhdrs.h"
#else
#include <sys/mman.h>
#include <pthread.h>
#include "../libbiokanga/commhdrs.h"
#endif

const char *cpszProgVer = "1.0.3";			// increment with each release


const int cMaxSupportedModes = 1;		// currently just the one processing mode supported

const int cDfltKMerLen = 1;				// default length kmer frequency to maintain when randomising genome 
const int cMaxKMerLen = 15;				// max length kmer frequency to maintain when randomising genome 

const int cMaxInFileSpecs = 10;			// allow at most a total of this many wildcarded control or experiment input read alignment loci files
const int cSeqBuffLen = 10000;		// buffer sequences in blocks of this length

int GenerateRandFasta(int Mode,					    // processing mode
				int KMerLen,						// maintain occurance frequency composition of this K-mer length 
				int NumInputControlSpecs,			// number of input file specs 
				char **pszInControlFiles,			// input genome or chromosome fasta files
				char *pszOutputFile,				// write randomised genome to this file
				int RandSeed);						// use this random seed

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
	return _T("kangade");
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
int Idx;
int iMode = 0;			// processing mode
int iRandSeed = 0;		// random base seed to use (if < 0 then current time used to seed generator)

int KMerLen;			// what length kmer to generate for 1..cMaxKMerLen

int NumInFileSpecs;			// number of input file specs 
char *pszInFiles[cMaxInFileSpecs];			// input control aligned reads files

char szOutputFile[_MAX_PATH];

// command line args
struct arg_lit  *help    = arg_lit0("hH","help",                "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to screen and logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"diagnostics log file");

struct arg_int  *Mode = arg_int0("m","mode","<int>",			"processing mode - 0 randomise genome");
struct arg_int  *kmerlen = arg_int0("k","kmerlen","<int>",		"maintain frequency composition of K-mer length (default = 1, range 1..15)");

struct arg_file *InFiles = arg_filen("i",NULL,"<file>",1,cMaxInFileSpecs, "input genome assembly multifasta files (s) to randomise");
struct arg_file *OutFile= arg_file1("o",NULL,"<file>",			"output randomised assembly to this file as multifasta");
struct arg_int *RandSeed = arg_int0("s","randseed","<int>",		"random seed to use (default is use system time)");
struct arg_end *end = arg_end(20);

void *argtable[] = {help,version,FileLogLevel,LogFile,Mode,kmerlen,InFiles,OutFile,RandSeed,end};

char **pAllArgs;
int argerrors;
argerrors = CUtility::arg_parsefromfile(argc,(char **)argv,&pAllArgs);
if(argerrors >= 0)
	argerrors = arg_parse(argerrors,pAllArgs,argtable);

/* special case: '--help' takes precedence over error reporting */
if (help->count > 0)
        {
		printf("\n%s Kanga randomise genome K-mers, Version %s\nOptions ---\n", gszProcName,cpszProgVer);
        arg_print_syntax(stdout,argtable,"\n");
        arg_print_glossary(stdout,argtable,"  %-25s %s\n");
		printf("\nNote: Parameters can be entered into a parameter file, one parameter per line.");
		printf("\n      To invoke this parameter file then precede its name with '@'");
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
		printf("\nError: FileLogLevel '-f%d' specified but no logfile '-F<logfile>\n'",FileLogLevel->ival[0]);
		exit(1);
		}

	iScreenLogLevel = iFileLogLevel = FileLogLevel->count ? FileLogLevel->ival[0] : eDLInfo;
	if(iFileLogLevel < eDLNone || iFileLogLevel > eDLDebug)
		{
		printf("\nError: FileLogLevel '-l%d' specified outside of range %d..%d\n",iFileLogLevel,eDLNone,eDLDebug);
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

	// now that log parameters have been parsed then initialise diagnostics log system
	if(!gDiagnostics.Open(szLogFile,(etDiagLevel)iScreenLogLevel,(etDiagLevel)iFileLogLevel,true))
		{
		printf("\nError: Unable to start diagnostics subsystem\n");
		if(szLogFile[0] != '\0')
			printf(" Most likely cause is that logfile '%s' can't be opened/created\n",szLogFile);
		exit(1);
		}

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Version: %s",cpszProgVer);

	iMode = Mode->count ? Mode->ival[0] : 0;
	if(iMode < 0 || iMode >= cMaxSupportedModes)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Unsupported Mode '-m%d' requested",iMode);
		exit(1);
		}

	KMerLen = 0;
	switch(iMode) {
		case 0:					// generate random species fasta sequence
			for(NumInFileSpecs=Idx=0;NumInFileSpecs < cMaxInFileSpecs && Idx < InFiles->count; Idx++)
				{
				pszInFiles[Idx] = NULL;
				if(pszInFiles[NumInFileSpecs] == NULL)
					pszInFiles[NumInFileSpecs] = new char [_MAX_PATH];
				strncpy(pszInFiles[NumInFileSpecs],InFiles->filename[Idx],_MAX_PATH);
				pszInFiles[NumInFileSpecs][_MAX_PATH-1] = '\0';
				CUtility::TrimQuotedWhitespcExtd(pszInFiles[NumInFileSpecs]);
				if(pszInFiles[NumInFileSpecs][0] != '\0')
					NumInFileSpecs++;
				}

			if(!NumInFileSpecs)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: After removal of whitespace, no input file(s) specified with '-i<filespec>' option)\n");
				exit(1);
				}


			KMerLen = kmerlen->count ? kmerlen->ival[0] : cDfltKMerLen;
			if(KMerLen < 1 || KMerLen > cMaxKMerLen)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: K-mer length specified with '-k%d' is outside of range 1..10",KMerLen);
				exit(1);
				}
			

			if(!OutFile->count || OutFile->filename[0][0] == '\0')
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: No output fasta file specified with '-o<filename>'");
				exit(1);
				}
			strncpy(szOutputFile,OutFile->filename[0],_MAX_PATH);
			szOutputFile[_MAX_PATH] = '\0';

			if(!RandSeed->count)
				iRandSeed = -1;
			else
				{
				iRandSeed = RandSeed->ival[0];
				if(iRandSeed < 0 || iRandSeed > 32767)
					{
					gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Random seed specified as '-s%d' must be between 0 and 32767",iRandSeed);
					exit(1);
					}
				}
			break;

		}

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing parameters:");



	switch(iMode) {
		case 0:
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Mode: 0 (randomise genome K-mers)");
			for(Idx=0; Idx < NumInFileSpecs; Idx++)
				gDiagnostics.DiagOutMsgOnly(eDLInfo,"Use frequency compositions from these genome multifasta file(s) (%d): '%s'",Idx+1,pszInFiles[Idx]);
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Maintain K-mer composition of length: %d",KMerLen);
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Genome output file: '%s'",szOutputFile);
			if(iRandSeed >= 0)
				gDiagnostics.DiagOutMsgOnly(eDLInfo,"Random seed: %d",iRandSeed);
			else
				gDiagnostics.DiagOutMsgOnly(eDLInfo,"Random seed: will use current time as seed");
			break;
		}
	
	gStopWatch.Start();
#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	Rslt = GenerateRandFasta(iMode,KMerLen,NumInFileSpecs,pszInFiles,szOutputFile,iRandSeed);
	Rslt = Rslt >=0 ? 0 : 1;
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Exit code: %d Total processing time: %s",Rslt,gStopWatch.Read());
	exit(Rslt);
	}
else
	{
	printf("\n%s Kanga randomise genome K-mers, Version %s\n",gszProcName,cpszProgVer);
	arg_print_errors(stdout,end,gszProcName);
	arg_print_syntax(stdout,argtable,"\nUse '-h' to view option and parameter usage\n");
	exit(1);
	}
return 0;
}

#pragma pack(1)
typedef struct TAGKMerCnts {
	UINT32 OrigCnts;		// number of original counts for this kmer
	UINT32 CurCnts;			// temp counts
	UINT32 GenCnts;			// number of counts generated for this kmer
} tsKMerCnts;
#pragma pack()

int m_SeqBuffLen;	// m_pSeq can hold at most this many characters
UINT8 *m_pSeq;		// to hold generated random sequence of bases a,c,g and t

UINT32 m_KmerCntInsts;				// number of tsKMerCnts in m_pKmerCnts
size_t m_AllocdKmerCntsSize;
tsKMerCnts *m_pKmerCnts;		// to hold observed KMer counts for current chromosome 

void
Reset(void)
{
if(m_pKmerCnts != NULL)
	{
#ifdef _WIN32
	free(m_pKmerCnts);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pKmerCnts != MAP_FAILED)
		munmap(m_pKmerCnts,m_AllocdKmerCntsSize);
#endif	
	m_pKmerCnts = NULL;
	}
if(m_pSeq != NULL)
	delete m_pSeq;
}

char *
KmerToText(int KMerLen,int KMerBases)
{
static char SeqTxt[100];
char *pSeqBase;
int Idx;
int RandSelBase;
pSeqBase = (char *)&SeqTxt[KMerLen];
*pSeqBase-- = '\0';
for(Idx = 0; Idx < KMerLen; Idx++, pSeqBase--)
	{
	RandSelBase = KMerBases & 0x03;
	KMerBases >>= 2;
	switch(RandSelBase) {
		case 0:
			*pSeqBase = 'a';
			break;
		case 1:
			*pSeqBase = 'c';
			break;
		case 2:
			*pSeqBase = 'g';
			break;
		case 3:
			*pSeqBase = 't';
			break;
		case 4:
			*pSeqBase = 'n';
			break;
		}
	}
return(SeqTxt);
}

int
GenRandChrom(char *pszDescription,
			int KMerLen,			// what length kmer to generate for - 1..15
			UINT32 MaxAllocBaseCnts,	// cnts allocated
			tsKMerCnts *pBaseCnts,UINT32 ChromLen,CFasta *pOutFasta,UINT32 SeqBuffLen,UINT8 *pSeqBuff,CRandomMersenne *pRandomiser)
{
UINT32 TotBaseCnts;		// to hold total base counts or current chromosome length
UINT32 TotKMerCnts;	    // total number of cnts for current length KMer-1
tsKMerCnts *pBaseCnt;
etSeqBase RandSelBase;	// base randomly selected from CompDist[]
UINT32 BaseIdx;
UINT32 KMerIdx;
UINT32 KMerMsk;

char *pSeqBase;			// used to access pSeq[]

int CompSeqLen;
UINT32 Idx;

KMerMsk = 0;
for(Idx = 0; Idx < (UINT32)KMerLen; Idx++)
	{
	KMerMsk <<= 2;
	KMerMsk |= 0x03;
	}
KMerMsk &= 0xfffffffc;

pOutFasta->WriteDescriptor(pszDescription);
TotBaseCnts = 0;
pBaseCnt = pBaseCnts;
for(BaseIdx = 0; BaseIdx < MaxAllocBaseCnts; BaseIdx++,pBaseCnt++)
	{
	TotBaseCnts += pBaseCnt->OrigCnts;
	pBaseCnt->CurCnts = pBaseCnt->OrigCnts; 
	pBaseCnt->GenCnts = 0;
	}

CompSeqLen = 0;
// randomly select an initial k-mer, small chance that this could be a 0 instance kmer so better check!
do {
	KMerIdx = pRandomiser->IRandom(0, MaxAllocBaseCnts-1);
	pBaseCnt = &pBaseCnts[KMerIdx];
	}
while(pBaseCnt->CurCnts < 1);


// now output the initial kmer sequence
int KMerBases = KMerIdx;
pBaseCnt->GenCnts += 1;
if(KMerLen == 1)
	pBaseCnt->CurCnts -= 1;
pSeqBase = (char *)&pSeqBuff[KMerLen-1];
for(Idx = 0; Idx < (UINT32)KMerLen; Idx++, pSeqBase--)
	{
	RandSelBase = KMerBases & 0x03;
	KMerBases >>= 2;
	switch(RandSelBase) {
		case 0:
			*pSeqBase = 'a';
			break;
		case 1:
			*pSeqBase = 'c';
			break;
		case 2:
			*pSeqBase = 'g';
			break;
		case 3:
			*pSeqBase = 't';
			break;
		case 4:
			*pSeqBase = 'n';
			break;
		}
	}


// use the KMerLen - 1 prefix bases as the index to select the the group from which to next select the emitted base
pSeqBase = (char *)&pSeqBuff[KMerLen];
CompSeqLen = KMerLen;
UINT32 CntIdx;
UINT32 AccumCnts;

for(; Idx < ChromLen; Idx++,pSeqBase++)
	{
	KMerIdx <<= 2;
	KMerIdx &= KMerMsk;
	pBaseCnt = &pBaseCnts[KMerIdx];
	TotKMerCnts = 0;
	for(CntIdx = 0; CntIdx < 4; CntIdx++,pBaseCnt++)
		TotKMerCnts += pBaseCnt->CurCnts;
	if(TotKMerCnts == 0)		// shouldn't happen but if no KMers starting with the current KMerLen-1 prefix then generate a random KMer
		{
		RandSelBase = pRandomiser->IRandom(0,3);
		}
	else
		{
		do {
			CntIdx = pRandomiser->IRandom(0, TotKMerCnts-1);
			AccumCnts = 0;
			pBaseCnt = &pBaseCnts[KMerIdx];
			for(RandSelBase = 0; RandSelBase < 4; RandSelBase++)
				{
				if(pBaseCnt->CurCnts >= 1)
					{
					AccumCnts += pBaseCnt->CurCnts;
					if(CntIdx < (int)AccumCnts)
						break;
					}
				pBaseCnt += 1;
				}
			}
		while(RandSelBase >= 4);
		}
	if(KMerLen == 1)						// can easily get the exact same frequency distribution for monomers
		pBaseCnt->CurCnts -= 1;
	
	pBaseCnt->GenCnts += 1;
	KMerIdx |= RandSelBase;
	switch(RandSelBase) {
		case 0:
			*pSeqBase = 'a';
			break;
		case 1:
			*pSeqBase = 'c';
			break;
		case 2:
			*pSeqBase = 'g';
			break;
		case 3:
			*pSeqBase = 't';
			break;
		case 4:
			*pSeqBase = 'n';
			break;
		}
	CompSeqLen += 1;
	if(CompSeqLen == SeqBuffLen-1)
		{
		pOutFasta->Write((char *)pSeqBuff,CompSeqLen);
		pSeqBase = (char *)pSeqBuff;
		CompSeqLen = 0;
		}
	}
if(CompSeqLen > 0)
	{
	pOutFasta->Write((char *)pSeqBuff,CompSeqLen);
	CompSeqLen = 0;
	}

if(KMerLen <= 5)
	{
	pBaseCnt = pBaseCnts;
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"\tKMer,ObsCnt,GenCnt");
	for	(Idx = 0; Idx < (UINT32)MaxAllocBaseCnts; Idx++,pBaseCnt++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"\t%s,%d,%d",KmerToText(KMerLen,Idx),pBaseCnt->OrigCnts,pBaseCnt->GenCnts);
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Sequence randomisation completed");
return(eBSFSuccess);
}



int GenerateRandFasta(int Mode,					    // processing mode
				int KMerLen,						// maintain occurance frequency composition of this K-mer length 
				int NumInFileSpecs,					// number of input file specs 
				char **pszInFiles,					// input genome or chromosome fasta files
				char *pszOutputFile,				// write randomised genome to this file
				int RandSeed)						// use this random seed
{
int BaseIdx;
UINT8 *pSeqBase;	// used to access m_pSeq[]
CFasta *pInFasta = NULL;
CFasta *pOutFasta = NULL;
char szDescription[cBSFDescriptionSize+1];
int ChromLen;
int ChromOfs;
int CompSeqLen;
bool bDescriptor;
int NumDescriptors;
int Idx;
int Rslt;

int SeqInLen;
int NumNs;
int KMerIdx;
UINT32 KMerMsk;
INT64 Now;

m_pKmerCnts = NULL;
m_AllocdKmerCntsSize = 0;
m_pSeq = NULL;
m_SeqBuffLen = cSeqBuffLen;

KMerMsk = 0;
m_KmerCntInsts = 1;
for(Idx = 0; Idx < KMerLen; Idx++)
	{
	m_KmerCntInsts *= 4;
	KMerMsk <<= 2;
	KMerMsk |= 0x03;
	}
KMerMsk &= 0xfffffffc;

if((m_pSeq = new UINT8 [m_SeqBuffLen]) == NULL)			
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate memory for sequence buffering %d bytes",m_SeqBuffLen);
	Reset();
	return(eBSFerrMem);
	}

pOutFasta = new CFasta();
if(pOutFasta == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CFasta object");
		Reset();
		return(eBSFerrObj);
		}

if((Rslt = pOutFasta->Open(pszOutputFile,false)) != eBSFSuccess)
	{
	while(pOutFasta->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,pOutFasta->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create fasta file '%s'",pszOutputFile);
	Reset();
	return(Rslt);
	}

if(RandSeed < 0)
	{
#ifdef _WIN32
	QueryPerformanceCounter((LARGE_INTEGER *)&Now);
#else
	struct tms Times;
	Now = (INT64)times(&Times);
#endif
	RandSeed = (int)(Now & 0x07fffffff);
	}
CRandomMersenne RandomiserBase(RandSeed);


m_AllocdKmerCntsSize = m_KmerCntInsts * sizeof(tsKMerCnts);
#ifdef _WIN32
m_pKmerCnts = (tsKMerCnts *) malloc((size_t)m_AllocdKmerCntsSize);	
if(m_pKmerCnts == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"GenerateRandFasta: Memory allocation of %lld bytes for K-mer count instances failed",(INT64)m_AllocdKmerCntsSize);
	m_AllocdKmerCntsSize = 0;
	Reset();
	delete pOutFasta;
	return(eBSFerrMem);
	}
#else
// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
m_pKmerCnts = (tsKMerCnts *)mmap(NULL,(size_t)m_AllocdKmerCntsSize, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
if(m_pKmerCnts == MAP_FAILED)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"GenerateRandFasta: Memory allocation of %lld bytes through mmap()  for K-mer count instances failed",(INT64)m_AllocdKmerCntsSize,strerror(errno));
	m_pKmerCnts = NULL;
	m_AllocdKmerCntsSize = 0;
	Reset();
	delete pOutFasta;
	return(eBSFerrMem);
	}
#endif


int FileIdx;
int NumInFilesProcessed;
char *pszInputFile;
CSimpleGlob glob(SG_GLOB_FULLSORT);
NumInFilesProcessed = 0;
for(FileIdx = 0; FileIdx < NumInFileSpecs; FileIdx++)
	{
	glob.Init();
	if(glob.Add(pszInFiles[FileIdx]) < SG_SUCCESS)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to glob '%s",pszInFiles[FileIdx]);
		Reset();
		delete pOutFasta;
		return(eBSFerrOpnFile);	// treat as though unable to open file
		}
	if(glob.FileCount() <= 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to locate any input file matching '%s",pszInFiles[FileIdx]);
		continue;
		}

	Rslt = eBSFSuccess;
	for (int FileID = 0; Rslt >= eBSFSuccess &&  FileID < glob.FileCount(); ++FileID)
		{
		pszInputFile = glob.File(FileID);
		NumInFilesProcessed += 1;
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"GenerateRandFasta: Loading fasta sequences from file: %s",pszInputFile);

		if((pInFasta = new CFasta()) == NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CFasta object");
			Reset();
			delete pOutFasta;
			return(eBSFerrObj);
			}

		if((Rslt=pInFasta->Open(pszInputFile,true)) < eBSFSuccess)
			{
			while(pInFasta->NumErrMsgs())
				gDiagnostics.DiagOut(eDLFatal,gszProcName,pInFasta->GetErrMsg());
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open fasta file '%s'",pszInputFile);
			Reset();
			delete pInFasta;
			delete pOutFasta;
			return(Rslt);
			}
		strcpy(szDescription,"ChrSim");		// default in case no fasta descriptor
		bDescriptor = false;
		NumDescriptors = 0;
		ChromLen = 0;
		ChromOfs = 0;
		CompSeqLen = 0;
		KMerIdx = 0;
		NumNs = 0;
		memset(m_pKmerCnts,0,m_AllocdKmerCntsSize);
		while((Rslt = SeqInLen = pInFasta->ReadSequence(m_pSeq,m_SeqBuffLen)) > eBSFSuccess)	
			{
			if(Rslt == eBSFFastaDescr)		// just read a descriptor line
				{
				NumDescriptors += 1;
				if(ChromLen >= KMerLen && NumDescriptors >= 2)		// if already read at least one descriptor then this terminates the previous chrom 
					GenRandChrom(szDescription,KMerLen,m_KmerCntInsts,m_pKmerCnts,ChromLen,pOutFasta,m_SeqBuffLen,m_pSeq,&RandomiserBase);
				pInFasta->ReadDescriptor(szDescription,cBSFDescriptionSize);
				memset(m_pKmerCnts,0,m_AllocdKmerCntsSize);
				ChromLen = 0;
				KMerIdx = 0;
				ChromOfs = 0;
				NumNs = 0;
				bDescriptor = true;
				gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing fasta sequence: '%s'",szDescription);
				continue;
				}

			if(ChromLen == 0 && SeqInLen < KMerLen)
				{
				gDiagnostics.DiagOut(eDLInfo,gszProcName,"Sloughing undersized sequence (%d bases) for '%s'",SeqInLen,szDescription);
				continue;
				}

			ChromLen += SeqInLen;
			pSeqBase = m_pSeq;
			for(Idx = 0; Idx < SeqInLen; Idx++,pSeqBase++)
				{
				BaseIdx = *pSeqBase  & 0x07;
				if(BaseIdx >= eBaseN)		// better safe than sorry - treat any non-canonical base as if it were indeterminate
					BaseIdx = RandomiserBase.IRandom(0, 3);	// use a random base
				KMerIdx <<= 2;
				KMerIdx &= KMerMsk;
				KMerIdx |= BaseIdx;
				if(++ChromOfs >= KMerLen)
					m_pKmerCnts[KMerIdx].OrigCnts += 1;
				}
			}

		if(ChromLen > 0)
			GenRandChrom(szDescription,KMerLen,m_KmerCntInsts,m_pKmerCnts,ChromLen,pOutFasta,m_SeqBuffLen,m_pSeq,&RandomiserBase);

		if(pInFasta != NULL)
			{
			delete pInFasta;
			pInFasta = NULL;
			}
		}
	}


if(pOutFasta != NULL)
	{
	if(CompSeqLen > 0 && m_pSeq != NULL)
		pOutFasta->Write((char *)m_pSeq,CompSeqLen);
	pOutFasta->Close();
	delete pOutFasta;
	pOutFasta = NULL;
	}

if(pInFasta != NULL)
	delete pInFasta;

Reset();
return(eBSFSuccess);
}

