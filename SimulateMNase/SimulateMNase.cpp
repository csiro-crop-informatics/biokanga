#include "stdafx.h"
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#if _WIN32

#include "../libbiokanga/commhdrs.h"
#else
#include "../libbiokanga/commhdrs.h"
#endif

const unsigned int cProgVer = 100;		// increment with each release

const int cDfltNumReads = 4000000;		// default number of reads
const int cMinNumReads =  10000;		// minimum number of reads
const int cMaxNumReads =  500000000;	// maximum number of reads
const int cDfltReadLen =  36;			// default read length
const int cMinReadLen  =  20;			// minimum read length
const int cMinCutLen   =  20;			// minimum cut length
const int cDfltCutRange = 60;			// targeted at nucleosomes, add to user entered min cut length to derive default max cut length
const int cDfltMinCutLen = (147 - cDfltCutRange/2);	// default minimum cut length
const int cMaxCutLen   =  2000;			// maximum cut length

// processing modes
typedef enum TAG_ePMode {		
	ePMMNaseRand,					// default - MNase start with random end sites
	ePMMNaseMNase,					// MNase start and end sites
	ePMMRandRand,					// random start and end
	ePMplaceholder					// used to set the enumeration range
	} etPMode;


int
Process(int PMode,			// processing mode
		char Strand,		// generate for this strand '+' or '-' or for both '*'
		int NumReads,		// number of reads required
		int ReadLen,		// read lengths
		int CutMin,			// min cut length
		int CutMax,			// max cut length
		char *pszInFile,	// input from this bioseq assembly
		char *pszMNaseFile,	// input from this MNase site preferences file
		char *pszOutFile);	// output simulated reads to this file

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
	return _T("usimreads");
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
int Rslt = 0;   			// function result code >= 0 represents success, < 0 on failure

int PMode;					// processing mode
char Strand;				// generate for this strand '+' or '-' or for both '*'
int NumReads;				// number of reads required
int ReadLen;				// read lengths
int CutMin;					// min cut length
int CutMax;					// max cut length


char szInFile[_MAX_PATH];	// input from this bioseq assembly
char szMNaseFile[_MAX_PATH];// input from this MNase site preferences file
char szOutFile[_MAX_PATH];	// output simulated reads to this file


// command line args
struct arg_lit  *help    = arg_lit0("hH","help",                "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_int *ScreenLogLevel=arg_int0("S", "ScreenLogLevel",	"<int>","Level of diagnostics written to screen 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"diagnostics log file");

struct arg_int *pmode = arg_int0("m","mode","<int>",		    "processing mode: 0 - MNase start with random end sites, 1 - MNase start and end sites, 2 - random start and end (0 - default)");
struct arg_str *strand=arg_str0("s", "strand","<str>",          "generate for this strand '+' or '-' only (default is both)");
struct arg_int *numreads = arg_int0("n","numreads","<int>",	    "number of reads required (default = 4000000)");
struct arg_int *readlen = arg_int0("l","length","<int>",	    "read lengths (default = 36)");
struct arg_int *cutmin = arg_int0("c","cutmin","<int>",		    "min cut length (default = 122)");
struct arg_int *cutmax = arg_int0("C","cutmax","<int>",		    "max cut length (default = 172)");
struct arg_file *infile = arg_file1("i","in","<file>",			"input from this bioseq assembly");
struct arg_file *inmnase = arg_file0("I","in","<file>",			"input from this MNase site preferences file");
struct arg_file *outfile = arg_file1("o","out","<file>",		"output simulated reads to this file");
struct arg_end *end = arg_end(20);

void *argtable[] = {help,version,FileLogLevel,ScreenLogLevel,LogFile,
					pmode,strand,numreads,readlen,cutmin,cutmax,infile,inmnase,outfile,
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

	PMode = (etPMode)(pmode->count ? pmode->ival[0] : 0);
	if(PMode < 0 || PMode >= ePMplaceholder)
		{
		printf("\nError: Processing mode '-m%d' specified outside of range %d..%d",PMode,0,(int)ePMplaceholder-1);
		exit(1);
		}

	Strand = strand->count ? *(char *)strand->sval[0] : '*';
	if(!(Strand == '+' || Strand == '-' || Strand == '*'))
		{
		printf("\nError: Strand specified '-s%c' must be one of '+', '-' or '*'",Strand);
		exit(1);
		}

	NumReads = numreads->count ? numreads->ival[0] : cDfltNumReads;
	if(NumReads < cMinNumReads || NumReads >= cMaxNumReads)
		{
		printf("\nError: Number of reads '-a%d' specified outside of range %d..%d",NumReads,cMinNumReads,cMaxNumReads);
		exit(1);
		}

	ReadLen = readlen->count ? readlen->ival[0] : cDfltReadLen;
	if(ReadLen < cMinReadLen || ReadLen > cMaxReadLen)
		{
		printf("\nError: Read length '-a%d' specified outside of range %d..%d",ReadLen,cMinReadLen,cMaxReadLen);
		exit(1);
		}

	CutMin = cutmin->count ? cutmin->ival[0] : cDfltMinCutLen;
	if(CutMin < cMinCutLen || CutMin > cMaxCutLen)
		{
		printf("\nError: Minimum cut length '-a%d' must be in range %d..%d",CutMin,cMinCutLen,cMaxCutLen);
		exit(1);
		}

	CutMax = cutmax->count ? cutmax->ival[0] : CutMin + cDfltCutRange;
	if(CutMax < CutMin || CutMax >= cMaxCutLen)
		{
		printf("\nError: Maximum cut length '-a%d' must be in range %d..%d",CutMax,CutMin,cMaxCutLen);
		exit(1);
		}

	if(PMode != ePMMRandRand)
		{
		if(!inmnase->count)
			{
			printf("\nError: No MNase site preference file specified with '-I<file>'");
			exit(1);
			}
		strncpy(szMNaseFile,inmnase->filename[0],_MAX_PATH);
		szMNaseFile[_MAX_PATH-1] = '\0';
		}
	else
		szMNaseFile[0] = '\0';

	strcpy(szInFile,infile->filename[0]);
	
	strcpy(szOutFile,outfile->filename[0]);

				// now that command parameters have been parsed then initialise diagnostics log system
	if(!gDiagnostics.Open(szLogFile,(etDiagLevel)iScreenLogLevel,(etDiagLevel)iFileLogLevel,true))
		{
		printf("\nError: Unable to start diagnostics subsystem.");
		if(szLogFile[0] != '\0')
			printf(" Most likely cause is that logfile '%s' can't be opened/created",szLogFile);
		exit(1);
		}
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Version: %d.%2.2d Processing parameters:",cProgVer/100,cProgVer%100);

	const char *pszDescr;
	switch(PMode) {
		case ePMMNaseRand:				
			pszDescr = "MNase start and random end sites";
			break;
			case ePMMNaseMNase:				
			pszDescr = "MNase start and MNase end sites";
			break;
		case ePMMRandRand:				
			pszDescr = "random start and random end sites";
			break;
		}
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"processing mode is : '%s'",pszDescr);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Simulate for this strand : '%c'",Strand);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"number of reads required: %d",NumReads);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"read lengths: %d",ReadLen);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"min cut length: %d",CutMin);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"max cut length: %d",CutMax);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"input bioseq assembly file: '%s'",szInFile);
	if(PMode != ePMMRandRand)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"input MNase site preferences file: '%s'",szMNaseFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"output simulated reads file: '%s'",szOutFile);


#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	gStopWatch.Start();
	Rslt = Process(PMode,Strand,NumReads,ReadLen,CutMin,CutMax,szInFile,szMNaseFile,szOutFile);
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
return 0;

	return 0;
}

typedef struct TAG_sChromSeq {
	int ChromID;
	int ScaledLen;			// chrom length scaled such that the sum of all chrom scaled lengths is less than INT_MAX
	char szChromName[cMaxDatasetSpeciesChrom];
	int Len;
	etSeqBase *pSeq;
} tsChromSeq;

int m_hOutFile;					// output results file handle
CCSVFile *m_pMNaseCSV;				// used to hold DNase site preferences whilst loading into m_pMNaseSel 
double *m_pMNaseSel;				// allocated array of MNase site selection preferences (0.0..1.0) indexed by sequence octamers

CBioSeqFile *m_pBioSeqFile;		// genome assembly
int m_NumChromSeqs;				// number of chromosomes loaded
int m_AllocdChromSeqs;			// number allocated
tsChromSeq *m_pChromSeqs;		// pts to chromseqs array
int m_GenomeScaledLen;			// sum of all chrom scaled lengths, will always be less than INT_MAX

void
Init(void)
{
m_pBioSeqFile = NULL;
m_pChromSeqs = NULL;
m_NumChromSeqs= 0;
m_AllocdChromSeqs = 0;
m_pMNaseSel = NULL;			// allocated array of MNase site selection preferences (0.0..1.0) indexed by sequence octamers
m_pMNaseCSV = NULL;			// used to hold DNase site preferences whilst loading into m_pMNaseSel 
m_hOutFile = -1;
}

void
Reset(void)
{
int Idx;
tsChromSeq *pChromSeq;
if(m_hOutFile != -1)
	{
	close(m_hOutFile);
	m_hOutFile = -1;
	}
if(m_pBioSeqFile != NULL)
	{
	delete m_pBioSeqFile;
	m_pBioSeqFile = NULL;
	}
if(m_pMNaseSel != NULL)
	{
	delete m_pMNaseSel;
	m_pMNaseSel = NULL;
	}
if(m_pMNaseCSV != NULL)
	{
	delete m_pMNaseCSV;
	m_pMNaseCSV = NULL;
	}
if(m_pChromSeqs != NULL)
	{
	pChromSeq = m_pChromSeqs;
	for(Idx = 0; Idx < m_NumChromSeqs; Idx++, pChromSeq++)
		{
		if(pChromSeq->pSeq != NULL)
			delete pChromSeq->pSeq;
		}
	delete m_pChromSeqs;
	m_pChromSeqs = NULL;
	}
}

// generate '+' strand index
int
GenPSeqIdx(int SeqLen,etSeqBase *pSeq)
{
int Idx;
int SeqIdx;
int Base;
for(Idx=SeqIdx=0; Idx < SeqLen; Idx++,pSeq++)
	{
	Base = *pSeq & ~cRptMskFlg;
	if(Base > eBaseT)
		return(-1);
	SeqIdx <<= 2;
	SeqIdx |= Base;
	}
return(SeqIdx);
}

// generate '-' strand index
int
GenMSeqIdx(int SeqLen,etSeqBase *pSeq)
{
int Idx;
int SeqIdx;
int Base;

pSeq += SeqLen-1;
for(Idx=SeqIdx=0; Idx < SeqLen; Idx++,pSeq--)
	{
	Base = *pSeq & ~cRptMskFlg;
	if(Base > eBaseT)
		return(-1);
	switch(Base) {
		case 0:
			Base = 3;
			break;
		case 1:
			Base = 2;
			break;
		case 2:
			Base = 1;
			break;
		case 3:
			Base = 0;
			break;
		}
	SeqIdx <<= 2;
	SeqIdx |= Base;
	}
return(SeqIdx);
}

char *
StepIdx2Seq(int SeqLen,int SeqIdx)
{
static char szSeqBuff[256];
char *pChr;
int Base;
int Idx;

szSeqBuff[SeqLen] = '\0';
pChr = &szSeqBuff[SeqLen-1];
for(Idx=0; Idx < SeqLen; Idx++,pChr--)
	{
	Base = SeqIdx & 0x03;
	switch(Base) {
		case 0: *pChr = 'a'; break;
		case 1: *pChr = 'c'; break;
		case 2: *pChr = 'g'; break;
		case 3: *pChr = 't'; break;
		}
	SeqIdx >>= 2;
	}
return(szSeqBuff);
}

int
InitMNaseSitePrefs(char *pszInMNaseFile)	// read from this MNase site selectivity file (generated by MNaseSitePred process)
{
int Rslt;
int NumProcessed;
int NumFields;
int OctIdx;
char *pszOctamer;
etSeqBase Octamer[9];
double SitePref;

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Loading MNase site selection preferences from file: %s",pszInMNaseFile);

m_pMNaseSel = (double *) new double [0x010000];	// to hold 4^8 (octamer) site preferences
if(m_pMNaseSel == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadReads: Memory allocation for MNase site preferences of 64k doubles failed");
	Reset();
	return(eBSFerrMem);
	}
for(OctIdx = 0; OctIdx < 0x010000; OctIdx++)
	m_pMNaseSel[OctIdx] = 0.0f;

if((m_pMNaseCSV = new CCSVFile) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CCSVfile");
	Reset();
	return(eBSFerrObj);
	}

if((Rslt=m_pMNaseCSV->Open(pszInMNaseFile))!=eBSFSuccess)
	{
	while(m_pMNaseCSV->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,m_pMNaseCSV->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open file: %s",pszInMNaseFile);
	Reset();
	return(Rslt);
	}

NumProcessed = 0;
while((Rslt=m_pMNaseCSV->NextLine()) > 0)	// onto next line containing fields
	{
	NumFields = m_pMNaseCSV->GetCurFields();
	if(NumFields < 4)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"file: %s contains % fields, expected at least 4",pszInMNaseFile,NumFields);
		Reset();
		return(eBSFerrParams);
		}
	if(!NumProcessed && m_pMNaseCSV->IsLikelyHeaderLine())
		continue;
	NumProcessed += 1;
	m_pMNaseCSV->GetText(1,&pszOctamer);
	memcpy(Octamer,CSeqTrans::MapAscii2Sense(pszOctamer),8);
	OctIdx = GenPSeqIdx(8,Octamer);
	m_pMNaseCSV->GetDouble(4,&SitePref);
	m_pMNaseSel[OctIdx] = SitePref;
	}

delete m_pMNaseCSV;
m_pMNaseCSV = NULL;

return(eBSFSuccess);
}

int
LoadGenome(char *pszBioSeqFile)
{
int Rslt;
int Len;
size_t TotLen;
int ChromID;
tsChromSeq *pChromSeq;
etSeqBase *pSeq;
double LenSCF;			// length scaling factor

if((m_pBioSeqFile = new CBioSeqFile) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CBioSeqFile object");
	return(eBSFerrObj);
	}
if((Rslt = m_pBioSeqFile->Open(pszBioSeqFile))!=eBSFSuccess)
	{
	while(m_pBioSeqFile->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,m_pBioSeqFile->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open assembly sequence file '%s'",pszBioSeqFile);
	Reset();
	return(Rslt);
	}

m_AllocdChromSeqs = m_pBioSeqFile->NumEntries();
if((m_pChromSeqs = new tsChromSeq [m_AllocdChromSeqs])==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadAssembly: unable to allocate memory (%d bytes) for %d ChromSeqs",sizeof(tsChromSeq) * m_NumChromSeqs,m_NumChromSeqs);
	Reset();
	return(eBSFerrMem);
	}
memset(m_pChromSeqs,0,sizeof(tsChromSeq) * m_AllocdChromSeqs);

TotLen = 0;
m_NumChromSeqs = 0;
ChromID = 0;

while((ChromID = m_pBioSeqFile->Next(ChromID)) > 0) {
	pChromSeq = &m_pChromSeqs[ChromID-1];
	pChromSeq->ChromID = ChromID;
	m_pBioSeqFile->GetName(ChromID,sizeof(pChromSeq->szChromName),pChromSeq->szChromName);
	Len = m_pBioSeqFile->GetDataLen(ChromID);
	TotLen += Len;
	if((pSeq = new etSeqBase [Len+1]) == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadAssembly: unable to allocate memory (%d bytes) for '%s' sequence",Len,pChromSeq->szChromName);
		Reset();
		return(eBSFerrMem);
		}
	m_pBioSeqFile->GetData(ChromID,eSeqBaseType,0,pSeq,Len);
	pChromSeq->Len = Len;
	pChromSeq->pSeq = pSeq;
	m_NumChromSeqs += 1;

	}
delete m_pBioSeqFile;
m_pBioSeqFile = NULL;


LenSCF = (double)TotLen/(double)INT_MAX;
m_GenomeScaledLen = TotLen >= (INT64)INT_MAX ? (long)(TotLen * LenSCF) : (long)TotLen;
pChromSeq = &m_pChromSeqs[0];
for(ChromID = 0; ChromID < m_NumChromSeqs; ChromID++,pChromSeq++)
	pChromSeq->ScaledLen = TotLen >= (INT64)INT_MAX ? (int)(pChromSeq->Len * LenSCF) : pChromSeq->Len;
return(eBSFSuccess);
}

int
Process(int PMode,			// processing mode
		char Strand,		// generate for this strand '+' or '-' or for both '*'
		int NumReads,		// number of reads required
		int ReadLen,		// read lengths
		int CutMin,			// min cut length
		int CutMax,			// max cut length
		char *pszInFile,	// input from this bioseq assembly
		char *pszMNaseFile,	// input from this MNase site preferences file
		char *pszOutFile)	// output simulated reads to this file
{
int Rslt;
char szLineBuff[0x07fff];
int LineLen;
int NumCols;
int ChromIdx;
int RandChrom;
int RandCutSite1;
int RandCutSite2;
int RandStrand;				// 0 if '+', 1 if '-'
int RandCutLen;
int OctIdx;
double RandCutProb;
double MNaseScore;
int NumReadsGenerated;
tsChromSeq *pChromSeq;
etSeqBase *pSeq;
etSeqBase *pSeqSite;

if(PMode != ePMMRandRand)
	{
	if((Rslt=InitMNaseSitePrefs(pszMNaseFile))!=eBSFSuccess)
		{
		Reset();
		return(Rslt);
		}
	}

if((Rslt=LoadGenome(pszInFile))!=eBSFSuccess)	// need to load complete assembly into memory
	{
	Reset();
	return(Rslt);
	}

#ifdef _WIN32
m_hOutFile = open(pszOutFile,O_CREATETRUNC );
#else
if((m_hOutFile = open(pszOutFile,O_RDWR | O_CREAT,S_IREAD | S_IWRITE))!=-1)
    if(ftruncate(m_hOutFile,0)!=0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate %s - %s",pszOutFile,strerror(errno));
			return(eBSFerrCreateFile);
			}
#endif

if(m_hOutFile < 0)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate output file '%s'",pszOutFile);
	return(eBSFerrCreateFile);
	}

TRandomCombined<CRandomMother,CRandomMersenne> RG((int)time(0));

NumReadsGenerated = 0;
printf("\nGenerating read %1.9d",NumReadsGenerated);
LineLen = 0;

while(NumReadsGenerated < NumReads)
	{
	// first choose chromosome
	RandChrom = (int)RG.IRandom(1,m_GenomeScaledLen);
	pChromSeq = &m_pChromSeqs[0];
	for(ChromIdx = 0; ChromIdx < (m_NumChromSeqs-1); ChromIdx++, pChromSeq++)
		{
		if(pChromSeq->ScaledLen >= RandChrom)
			break;
		RandChrom -= pChromSeq->ScaledLen;
		}
	if(pChromSeq->Len < (CutMax + 10))	// skip any extremely short chromosomes - could be a contig?
		continue;

	// randomly choose strand?
	if(Strand == '*')
		RandStrand = (int)RG.IRandom(0,1);
	else
		RandStrand = Strand == '-' ? 1 : 0;

	// randomly choose cut length
	RandCutLen = (int)RG.IRandom(CutMin,CutMax);
	// randomly choose initial cut site
	// note that specified range is such that '+' and '-' start/end loci after allowing for cut lengths will always be on the chromosome
	RandCutSite1 = (int)RG.IRandom(5,pChromSeq->Len - (RandCutLen + 5));

	if(PMode != ePMMRandRand)
		{
		// randomly choose site1 cut prob
		RandCutProb = RG.Random();
		pSeq = pChromSeq->pSeq;
		pSeqSite = &pSeq[RandCutSite1-4];
		if(!RandStrand)
			OctIdx = GenPSeqIdx(8,pSeqSite);
		else
			OctIdx = GenMSeqIdx(8,pSeqSite);
		if(OctIdx < 0)		// -1 if pSeqSite contained 'n'
			continue;
		MNaseScore = m_pMNaseSel[OctIdx];
		if(MNaseScore < RandCutProb)
			continue;
		}
	
	pSeq = pChromSeq->pSeq;
	RandCutSite2 = RandCutSite1 + RandCutLen;

	if(PMode == ePMMNaseMNase)
		{
		// randomly choose site2 cut prob
		RandCutProb = RG.Random();
		pSeqSite = &pSeq[RandCutSite2-4];
		if(!RandStrand)
			OctIdx = GenPSeqIdx(8,pSeqSite);
		else
			OctIdx = GenMSeqIdx(8,pSeqSite);
		if(OctIdx < 0)		// -1 if seq contained 'n'
			continue;
		MNaseScore = m_pMNaseSel[OctIdx];
		if(MNaseScore < RandCutProb)
			continue;
		}

	// ensure that this sequence will be exactly alignable - slough if any contained 'n's
	pSeqSite = &pSeq[RandCutSite1];
	for(OctIdx = RandCutSite1; OctIdx <= RandCutSite2; OctIdx++,pSeqSite++)
		if((*pSeqSite & ~cRptMskFlg) > eBaseT)
			break;
	if(OctIdx <= RandCutSite2)
		continue;
	
	
	// we have a cut MNase cut sequence starting at RandCutSite1 and ending at RandCutSite2-1 which is of length RandCutLen !!!!!!
	NumReadsGenerated += 1;
	if(!(NumReadsGenerated % 10000))
		printf("\b\b\b\b\b\b\b\b\b%1.9d",NumReadsGenerated);

	LineLen+=sprintf(&szLineBuff[LineLen],"%s|simmnase %1.9d %s|%d|%d|%d|%c\n",
			NumReadsGenerated?">lcl":"\n>lcl",NumReadsGenerated,pChromSeq->szChromName,RandCutSite1,RandCutSite2,RandCutLen,RandStrand ? '-' : '+');

	if((LineLen + 200) > sizeof(szLineBuff))
		{
		if(write(m_hOutFile,szLineBuff,LineLen) != LineLen)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Errors whilst writing %d bytes - '%s' - %s",LineLen, pszOutFile, strerror(errno));
			Reset();
			return(eBSFerrFileAccess);
			}

		LineLen=0;
		}

	if(RandStrand)		// if on '-' strand
		RandCutSite1 = 1 + RandCutSite2 - ReadLen;
	if(RandCutSite1 < 0)
		printf("\nA problem!");

	RandCutLen = ReadLen;
	if(RandStrand)
		{
		RandCutSite2 = RandCutSite1;
		CSeqTrans::ReverseComplement(ReadLen,&pSeq[RandCutSite2]);
		}

	while(RandCutLen)
		{
		NumCols = RandCutLen > 70 ? 70 : RandCutLen;
		CSeqTrans::MapSeq2Ascii(&pSeq[RandCutSite1],NumCols,&szLineBuff[LineLen]);
		LineLen += NumCols;
		LineLen += sprintf(&szLineBuff[LineLen],"\n");
		RandCutLen -= NumCols;
		RandCutSite1 += NumCols;
		}
	if((LineLen + 200) > sizeof(szLineBuff))
		{
		if(write(m_hOutFile,szLineBuff,LineLen) != LineLen)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Errors whilst writing %d bytes - '%s' - %s",LineLen, pszOutFile, strerror(errno));
			Reset();
			return(eBSFerrFileAccess);
			}
		LineLen=0;
		}
	if(RandStrand)
		CSeqTrans::ReverseComplement(ReadLen,&pSeq[RandCutSite2]);
	}
if(LineLen && write(m_hOutFile,szLineBuff,LineLen) != LineLen)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Errors whilst writing %d bytes - '%s' - %s",LineLen, pszOutFile, strerror(errno));
	Reset();
	return(eBSFerrFileAccess);
	}

printf("\b\b\b\b\b\b\b\b%1.8d",NumReadsGenerated);

close(m_hOutFile);
m_hOutFile = -1;
Reset();
return(eBSFSuccess);
}

