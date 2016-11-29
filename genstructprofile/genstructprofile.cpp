// genstructprofile.cpp : Defines the entry point for the console application.
// Loads a set of fasta reads and attempts to find a structural profile
// common to a consensus of reads for a specified conformational characteristic

#include "stdafx.h"
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#if _WIN32
#include "../libbiokanga/commhdrs.h"
#else
#include "../libbiokanga/commhdrs.h"
#endif

const unsigned int cProgVer = 132;		// increment with each release

const int cMaxProcSeqLen = 1000;		// max length sequence processed - longer are truncated to this length
const int cTruncFastaDescrLen = 70;		// truncate fasta descriptors if this length or longer

const double cMinBkgndGroove = 11.10f; // min background minor groove width
const double cDfltBkgndGroove = 11.20f; // default background minor groove width
const double cMaxBkgndGroove = 11.40f; // max background minor groove width
const double cDeltaBkgndGroove = 0.025f; // incr background minor groove width by this when iterating from cMinBkgndGroove to cMaxBkgndGroove

const double cMinDyadratio = 1.000f;	// min dyad grooves must be at least this ratio to background
const double cDfltDyadratio = 1.030f;	// default dyad grooves must be at least this ratio to background
const double cMaxDyadratio = 1.050f;	// max dyad grooves must be at least this ratio to background
const double cDeltaDyadratio = 0.005f;	// incr dyad grooves by this when iterating from cMinDyadratio to cMaxDyadratio

const double cMinDyad2ratio = 1.000f;	// min immediately flanking grooves must be at least this ratio to background
const double cDfltDyad2ratio = 1.020f;	// default immediately flanking grooves must be at least this ratio to background
const double cMaxDyad2ratio = 1.040f;	// max immediately flanking grooves must be at least this ratio to background
const double cDeltaDyad2ratio = 0.005f;	// incr immediately flanking grooves by this when iterating from cMinDyad2ratio to cMaxDyad2ratio

const double cMinDyad3ratio = 1.000f;	// min remainder of flanking grooves must be at least this ration to background
const double cDfltDyad3ratio = 1.015f;	// default remainder of flanking grooves must be at least this ration to background
const double cMaxDyad3ratio = 1.030f;	// max remainder of flanking grooves must be at least this ration to background
const double cDeltaDyad3ratio = 0.005f;	// incr flanking grooves by this when iterating from cMinDyad3ratio to cMaxDyad3ratio

const int cMinProfileLen = 125;			// lower limit - just enough for minor groove to be determined at dyad and at 6pts left and right (total 13) every 360 degrees accumulated twist 
const int cDfltProfileLen = 147;		// assume looking for nucleosome profiles

const int cDfltNSamples = 100;			// default number of seleted or random samples - purely arbitary
const int cMaxNSamples = 4000000;		// upper limit on selected or random samples - purely arbitary

// processing modes
typedef enum TAG_ePMode {		
	ePMdefault,							// default processing mode is to process all reads
	ePMNsel,							// select 1st N reads to be processed
	ePMrandsel,							// randomly select N reads to be processed
	ePMNallOpt,							// process all reads for parameter optimisation
	ePMNselOpt,							// select 1st N reads reads for parameter optimisation
	ePMrandselOpt,						// randomly select N reads for parameter optimisation
	ePMplaceholder						// used to set the enumeration range
	} etPMode;

#pragma pack(1)
typedef struct TAG_sSeqValues {
	int SeqID;							// uniquely identifies this sequence (1..m_NumSeqValues)
	int Sizeof;							// total size, in bytes, of this instance
	UINT8 RandSel:1;					// currently reserved for use to indicate that this sequence was randomly selected
	int SeqLen;							// sequence length
	int LenDescr;						// length of fasta descriptor (appended after conformational values)
	int DyadIdx;						// step index (0..NumSteps-1) for putative dyad at which maximal profile score was obtained for this sequence (-1 if none)
	int Score;							// maximal score for this sequence
	UINT8 SeqDescr[1];				    // sequence followed by appended '\0' terminated fasta descriptor
} tsSeqValues;

typedef struct TAG_sConfRange {
	int Min;
	int Max;
	int Range;
} tsConfRange;
#pragma pack()


const int cAllocSeqSize = 100000 * (sizeof(tsSeqValues) + (2 * cMaxProcSeqLen) + cTruncFastaDescrLen + 1);	// allocation block size for holding tsSeqValues

etPMode m_PMode;						// processing mode

CFasta *m_pFasta;						// fasta file containing nucleosome wrap sequences 
CTwister *m_pTwister;					// contains hunter group dsDNA conformational characteristic values
int m_hProfileFile;						// file handle for output profile file

tsConfRange m_ConfRange;				// min/max/range value limits for selected minor groove conformational parameter
tsConfRange m_ConfTwistRange;			// min/max/range value limits specifically for twist

typedef struct TAG_sSeqSet {
	char szFastaFile[_MAX_PATH];		// file from which these sequences were read
	int NumSeqsProc;					// number of sequences processed
	int NumShortSeqs;					// number of sequences not accepted because they are shorter than ProfLen
	int NumTruncSeqs;					// number of sequences that were truncated as longer than cMaxProcSeqLen
	tsSeqValues *pSeqValues;			// memory allocated to hold loaded sequences
	size_t AllocdMemSeqValues;			// how many bytes of memory have been allocated to hold the loaded sequences
	size_t UsedMemSeqValues;			// how many allocated bytes of memory have been used to hold the loaded sequences
	int NumSeqValues;					// number of sequences loaded into memory
	int NSamples;						// sample this number of sequences
	int NumDyads;						// number of sequences with at least one dyad
	tsSeqValues *pIterSeqValues;		// marks last iterated sequence
	}	tsSeqSet;

tsSeqSet m_ExprSet;			// holds the loaded sequences deemed as being the experimental sequences
tsSeqSet m_CtrlSet;			// holds the loaded sequences deemed as being the control sequences

int 
LoadSequences(etPMode PMode,			// processing mode
			  tsSeqSet *pSeqSet,		// load into this sequence set
	  		  int NSamples,				// sample this number of sequences if non-zero
			  char *pszFastaFile);		// fasta reads

int
Process(etPMode PMode,				// processing mode
		int TruncLength,			// truncate sequences to be no  longer than this length after any offseting
		int OfsStart,				// offset start of sequences by this many bases
		int NSamples,				// sample this number of sequences
		double BkgndGroove,			// background minor groove (if 0 then determine background from the mean of all sequences) 
		double DyadratioThres,		// dyad grooves must be at least this ratio to background
		double Dyad2ratioThres,		// immediately flanking grooves must be at least this ratio to background
		double Dyad3ratioThres,		// remainder of flanking grooves must be at least this ration to background
		char *pszFastaFile,			// input reads fastas
		char *szFastaCtrlFile,		// input reads fastas for control parameterisation objective function mode input
		char *pszStructParamsFile,	// file containing conformation structural parameters
		char *pszProfileFile);		// write profiles to this file

int 
LocateDyads(etPMode PMode,			// processing mode
			tsSeqSet *pSeqSet,		// interested in dyads for this sequence set
			int OfsStart,			// offset start of sequences by this many bases
			int TruncLength,		// truncate sequences to be no  longer than this length after any offseting
			double BkgndGroove,			// dyads are relative to this background groove
			double DyadratioThres,		// dyad grooves must be at least this ratio to background
			double Dyad2ratioThres,		// immediately flanking grooves must be at least this ratio to background
			double Dyad3ratioThres);	// remainder of flanking grooves must be at least this ration to background

tsSeqValues *
IterNext(tsSeqSet *pCtrlSet,	// iterate this set of sequences
		 bool bFirst,		// true if first sequence to be iterated, subsequent calls must have bFirst set false
		 bool bRandSels);	// true if only sequences marked as RandSel to be returned

bool NxtOfsCombination(int ProfLen);	// returns true if more combinations, false if all combinations tried
int MarkRandoms(tsSeqSet *pSeqSet,int NSamples,bool bReset); // Randomly selects NSamples sequences and marks these as having been randomly selected

void Init(void);
void Reset(void);

CStopWatch gStopWatch;
CDiagnostics gDiagnostics;				// for writing diagnostics messages to log file
char gszProcName[_MAX_FNAME];			// process name
bool gbActivity;						// used to determine if activity messages vi printf's can be used - output activity if eDLInfo or less
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

etPMode PMode;				// processing mode
int NSamples;				// process this number of sequences
int TruncLength;			// truncate sequences to be no  longer than this length after any offseting
int OfsStart;				// offset start of sequences by this many bases
double BkgndGroove;			// background minor groove (if not specified then defaults to 11.12)
double Dyadratio;
double Dyad2ratio;
double Dyad3ratio;

char szFastaFile[_MAX_PATH];	// sequences from this file
char szFastaCtrlFile[_MAX_PATH];	// sequences from this file when in parameterisation objective function mode

char szStructParamsFile[_MAX_PATH]; // file containing conformation structural parameters
char szProfileFile[_MAX_PATH];	// output profile to this file

// command line args
struct arg_lit  *help    = arg_lit0("hH","help",                "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_int *ScreenLogLevel=arg_int0("S", "ScreenLogLevel",	"<int>","Level of diagnostics written to screen 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"diagnostics log file");

struct arg_int *pmode = arg_int0("m","mode","<int>",		    "processing mode: 0 - process all sequences, 1 - select 1st N samples, 2 - randomly select N sample sequences, 3-5 for iterative parameter optimisation");
struct arg_int *nsamples = arg_int0("n", "nsamples","<int>",	"number of samples in processing mode 1 and 2");
struct arg_int  *trunclength = arg_int0("T","truncatelength",   "<int>","truncate sequences to be no  longer than this length after any offseting, 125..1000 (default 300)");
struct arg_int  *ofsstart   = arg_int0("u","ofsstart","<int>",	"offset start of sequences by this many bases, 0 - 1000 (default 0)");
struct arg_dbl *bkgndgroove=arg_dbl0("b", "bkgndgroove",	"<double>","background minor groove (default 11.1200)");
struct arg_dbl *dyadratio=arg_dbl0("d", "dyadratio",	"<double>","dyad minor grooves must be at least this ratio to background (default 1.030");
struct arg_dbl *dyad2ratio=arg_dbl0("D", "dyad2ratio",	"<double>","immediately flanking minor grooves must be at least this ratio to background (default 1.020");
struct arg_dbl *dyad3ratio=arg_dbl0("e", "dyad3ratio",	"<double>","remainder flanking minor grooves must be at least this ratio to background (default 1.015");
struct arg_file *infile = arg_file1("i","in","<file>",				"sequences from this multifasta file");
struct arg_file *inctrlfile = arg_file0("I","inctrlfile","<file>",	"when iterating for parameter optimisation use this fasta file as the control");
struct arg_file *StructParams = arg_file1("p","params","<file>","file containing conformation structural parameters");
struct arg_file *outfile = arg_file1("o","out","<file>",		"output profile to this file");

struct arg_end *end = arg_end(20);
void *argtable[] = {help,version,FileLogLevel,ScreenLogLevel,LogFile,
					pmode,ofsstart,trunclength,bkgndgroove,dyadratio,dyad2ratio,dyad3ratio,nsamples,StructParams,infile,inctrlfile,outfile,
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

	PMode = (etPMode)(pmode->count ? pmode->ival[0] : ePMdefault);
	if(PMode < 0 || PMode >= ePMplaceholder)
		{
		printf("\nError: Processing mode '-m%d' specified outside of range %d..%d",PMode,ePMdefault,(int)ePMplaceholder-1);
		exit(1);
		}

	OfsStart = ofsstart->count ? ofsstart->ival[0] : 0;
	if(OfsStart < 0 || OfsStart > 1000)
		{
		printf("\nError: Sequence offset '-u%d' specified outside of range 0..1000",OfsStart);
		exit(1);
		}

	TruncLength = trunclength->count ? trunclength->ival[0] : 0;
	if(TruncLength != 0 && (TruncLength < 125 || TruncLength > 1000))
		{
		printf("\nError: Truncate length mode '-T%d' specified outside of range 125..1000",TruncLength);
		exit(1);
		}


	if(PMode != ePMdefault && PMode != ePMNallOpt)
		{
		NSamples = nsamples->count ? nsamples->ival[0] : cDfltNSamples;
		if(NSamples < 1 || NSamples > cMaxNSamples)
			{
			printf("\nError: number of samples '-n%d' specified outside of range %d..%d",NSamples,1,cMaxNSamples);
			exit(1);
			}
		}
	else
		NSamples = 0;

	if(PMode <= ePMrandsel)
		{
		BkgndGroove = bkgndgroove->count ? bkgndgroove->dval[0] : cDfltBkgndGroove;
		if(BkgndGroove != 0.0f && (BkgndGroove < cMinBkgndGroove || BkgndGroove > cMaxBkgndGroove))
			{
			printf("\nError: Bacgound minor groove '-b%1.4f' must be either 0.0 or in range %1.4f to %1.4f",BkgndGroove,cMinBkgndGroove,cMaxBkgndGroove);
			exit(1);
			}

		Dyadratio = dyadratio->count ? dyadratio->dval[0] : cDfltDyadratio;
		if(Dyadratio < cMinDyadratio || Dyadratio > cMaxDyadratio)
			{
			printf("\nError: Centrall dyad threshold ratio '-d%1.4f' must be in range %1.4f to %1.4f",Dyadratio,cMinDyadratio,cMaxDyadratio);
			exit(1);
			}


		Dyad2ratio = dyad2ratio->count ? dyad2ratio->dval[0] : cDfltDyad2ratio;
		if(Dyad2ratio < cMinDyad2ratio || Dyadratio > cMaxDyad2ratio)
			{
			printf("\nError: Centrall dyad threshold ratio '-d%1.4f' must be in range %1.4f to %1.4f",Dyad2ratio,cMinDyad2ratio,cMaxDyad2ratio);
			exit(1);
			}


		Dyad3ratio = dyad3ratio->count ? dyad3ratio->dval[0] : cDfltDyad3ratio;
		if(Dyad3ratio < cMinDyad2ratio || Dyadratio > cMaxDyad2ratio)
			{
			printf("\nError: Centrall dyad threshold ratio '-d%1.4f' must be in range %1.4f to %1.4f",Dyad3ratio,cMinDyad3ratio,cMaxDyad3ratio);
			exit(1);
			}
		}
	else
		{
		Dyadratio = cMinDyadratio;
		Dyad2ratio = cMinDyad2ratio;
		Dyad3ratio = cMinDyad3ratio;
		BkgndGroove = cMinBkgndGroove;
		}

	strcpy(szFastaFile,infile->filename[0]);

	if(PMode >= ePMNallOpt && inctrlfile->count)
		strcpy(szFastaCtrlFile,inctrlfile->filename[0]);
	else
		szFastaCtrlFile[0] = '\0';

	strcpy(szStructParamsFile,StructParams->filename[0]);
	strcpy(szProfileFile,outfile->filename[0]);

	// now that command parameters have been parsed then initialise diagnostics log system
	if(!gDiagnostics.Open(szLogFile,(etDiagLevel)iScreenLogLevel,(etDiagLevel)iFileLogLevel,true))
		{
		printf("\nError: Unable to start diagnostics subsystem.");
		if(szLogFile[0] != '\0')
			printf(" Most likely cause is that logfile '%s' can't be opened/created",szLogFile);
		exit(1);
		}

// show user current resource limits
#ifndef _WIN32
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Resources: %s",CUtility::ReportResourceLimits());
#endif

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Version: %d.%2.2d Processing parameters:",cProgVer/100,cProgVer%100);
	const char *pszDescr;
	switch(PMode) {
		case ePMdefault:
			pszDescr = "generate conformational profile for all reads";
			break;

		case ePMNsel:					
			pszDescr = "conformational profile for 1st N reads";
			break;

		case ePMrandsel:				
			pszDescr = "conformational profile for randomly selected N reads";
			break;

		case ePMNallOpt:				
			pszDescr = "process all reads iterating for parameter optimisation";
			break;

		case ePMNselOpt:		
			pszDescr = "process 1st N reads reads for parameter optimisation";
			break;
		case ePMrandselOpt:				
			pszDescr = "process randomly selected N reads for parameter optimisation";
			break;
		}
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"processing mode is : '%s'",pszDescr);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"offset sequence starts by : %dnt",OfsStart);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"truncate sequence after offset to be : %dnt",TruncLength);

	if(PMode != ePMdefault && PMode != ePMNallOpt)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"sample this number of sequences : '%d'",NSamples);

	if(PMode <= ePMrandsel)
		{
		if(BkgndGroove > 0.0f)
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"background minor groove: %1.4f",BkgndGroove);
		else
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"background minor groove: use mean of all sequences");

		gDiagnostics.DiagOutMsgOnly(eDLInfo,"dyad minor groove threshold: %1.4f",Dyadratio);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"dyad immediately flanking minor groove threshold: %1.4f",Dyad2ratio);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"dyad remainder flanking minor groove threshold: %1.4f",Dyad3ratio);
		}

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"input reads file: '%s'",szFastaFile);
	if(szFastaCtrlFile[0] != '\0')
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"parameterisation objective function mode input reads file : '%s'",szFastaCtrlFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"conformation struct param file: '%s'",szStructParamsFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"output profile file: '%s'",szProfileFile);

#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	gStopWatch.Start();
	Rslt = Process(PMode,TruncLength,OfsStart,NSamples,BkgndGroove,Dyadratio,Dyad2ratio,Dyad3ratio,szFastaFile,szFastaCtrlFile,szStructParamsFile,szProfileFile);
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
}

void
Init(void)
{
m_pFasta = NULL;
m_pTwister = NULL;
m_ExprSet.pSeqValues = NULL;
m_CtrlSet.pSeqValues = NULL;
m_ExprSet.pIterSeqValues = NULL;
m_CtrlSet.pIterSeqValues = NULL;
m_hProfileFile = -1;
Reset();
}

void
Reset(void)
{
if(m_hProfileFile != -1)
	{
	close(m_hProfileFile);
	m_hProfileFile = -1;
	}
if(m_pFasta != NULL)					// fasta file containing nucleosome wrap sequences 
	{
	delete m_pFasta;
	m_pFasta = NULL;
	}

if(m_pTwister != NULL)					// contains hunter group dsDNA conformational characteristic values
	{
	delete m_pTwister;
	m_pTwister = NULL;
	}

if(m_ExprSet.pSeqValues != NULL)			// loaded sequences
	{
	free(m_ExprSet.pSeqValues);				// malloc/realloc was used to alloc memory
	m_ExprSet.pSeqValues = NULL;
	}
if(m_CtrlSet.pSeqValues != NULL)			// loaded sequences
	{
	free(m_CtrlSet.pSeqValues);				// malloc/realloc was used to alloc memory
	m_CtrlSet.pSeqValues = NULL;
	}

memset(&m_ExprSet,0,sizeof(m_ExprSet));
memset(&m_CtrlSet,0,sizeof(m_CtrlSet));
m_PMode = ePMdefault;					// processing mode
memset(&m_ConfRange,0,sizeof(m_ConfRange));
memset(&m_ConfTwistRange,0,sizeof(m_ConfRange));
}

int
Process(etPMode PMode,				// processing mode
		int TruncLength,			// truncate sequences to be no  longer than this length after any offseting
		int OfsStart,				// offset start of sequences by this many bases
		int NSamples,				// sample this number of sequences
		double BkgndGroove,			// background minor groove (if 0 then determine background from the mean of all sequences) 
		double DyadratioThres,		// dyad grooves must be at least this ratio to background
		double Dyad2ratioThres,		// immediately flanking grooves must be at least this ratio to background
		double Dyad3ratioThres,		// remainder of flanking grooves must be at least this ration to background
		char *pszFastaFile,			// structural profile for these fastas
		char *szFastaCtrlFile,		// input reads fastas for control parameterisation objective function mode input
		char *pszStructParamsFile,	// file containing conformation structural parameters
		char *pszProfileFile)		// write profiles to this file
{
int Rslt;
bool bFirst;
tsSeqValues *pSeqValues;
int Ident;
char szSpecies[70];
char szChrom[70];
int StartLoci;
int EndLoci;
int ElLen;
char Strand;
int NumReads;
char *pszDescr;
int NumFields;
int NumDyads;
char szBuff[0x7fff];
int BuffOfs;
Init();
m_PMode = PMode;

// load the reads into memory, note that the loaded reads will be sorted by ReadID
if((Rslt = LoadSequences(PMode,&m_ExprSet,NSamples,pszFastaFile)) < eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to load sequences");
	return(Rslt);
	}

if((PMode == ePMrandsel || PMode == ePMrandselOpt) && NSamples < m_ExprSet.NumSeqValues)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Randomly marking %d sequences for processing...",NSamples);
	MarkRandoms(&m_ExprSet,NSamples,true);
	m_ExprSet.NSamples = NSamples;
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Random marking completed");
	}
else
	m_ExprSet.NSamples = m_ExprSet.NumSeqValues;

if(PMode >= ePMNallOpt && szFastaCtrlFile[0] != '\0')
	{
	// load the control reads into memory, note that the loaded reads will be sorted by ReadID
	if((Rslt = LoadSequences(PMode,&m_CtrlSet,NSamples,szFastaCtrlFile)) < eBSFSuccess)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to load sequences");
		return(Rslt);
		}

	if(PMode == ePMrandselOpt && NSamples < m_CtrlSet.NumSeqValues)
		{
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Randomly marking %d sequences for processing...",NSamples);
		MarkRandoms(&m_CtrlSet,NSamples,true);
		m_CtrlSet.NSamples = NSamples;
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Random marking completed");
		}
	else
		m_CtrlSet.NSamples = m_CtrlSet.NumSeqValues;
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Loading conformational characteristic values from file: '%s'",pszStructParamsFile);

if((m_pTwister = new CTwister)==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,"LoadSequences","Unable to create CTwister object");
	Reset();
	return(eBSFerrObj);
	}

if((Rslt = m_pTwister->LoadStructParams(pszStructParamsFile))  < eBSFSuccess)
	{
	while(m_pTwister->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,m_pTwister->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,"LoadSequences","LoadStructParams(%s) failed",pszStructParamsFile);
	Reset();
	return(Rslt);
	}

m_ConfRange.Min = min(m_pTwister->m_StructParamStats[eSSminorgroove].Max,m_pTwister->m_StructParamStats[eSSminorgroove].Min);
m_ConfRange.Max = max(m_pTwister->m_StructParamStats[eSSminorgroove].Max,m_pTwister->m_StructParamStats[eSSminorgroove].Min);
m_ConfRange.Range = abs(m_pTwister->m_StructParamStats[eSSminorgroove].Max - m_pTwister->m_StructParamStats[eSSminorgroove].Min);
m_ConfTwistRange.Min = min(m_pTwister->m_StructParamStats[eSStwist].Max,m_pTwister->m_StructParamStats[eSStwist].Min);
m_ConfTwistRange.Max = max(m_pTwister->m_StructParamStats[eSStwist].Max,m_pTwister->m_StructParamStats[eSStwist].Min);
m_ConfTwistRange.Range = abs(m_pTwister->m_StructParamStats[eSStwist].Max - m_pTwister->m_StructParamStats[eSStwist].Min);

if(PMode>= ePMNallOpt)
	{
#ifdef _WIN32
	if((m_hProfileFile = open(pszProfileFile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
	if((m_hProfileFile = open(pszProfileFile, O_RDWR | O_CREAT | O_TRUNC, S_IREAD | S_IWRITE ))==-1)
#endif
		{
		gDiagnostics.DiagOut(eDLFatal,"Process","Unable to create or truncate results profile file %s error: %s",pszProfileFile,strerror(errno));
		Reset();
		return(eBSFerrCreateFile);
		}
	int Iteration = 0;
	printf("Starting iteration %1.8d",Iteration);

	BuffOfs = sprintf(szBuff,"\"BkgndGroove\",\"DyadratioThres\",\"Dyad2ratioThres\",\"Dyad3ratioThres\",\"ExprSet.NSamples\",\"ExprSet.NumDyads\"");
	if(szFastaCtrlFile[0] != '\0')	
		BuffOfs += sprintf(&szBuff[BuffOfs],",\"CtrlSet.NSamples\",\"CtrlSet.NumDyads\"");
	CUtility::SafeWrite(m_hProfileFile,szBuff,BuffOfs);
	BuffOfs = 0;
	for(BkgndGroove = cMinBkgndGroove; BkgndGroove <= cMaxBkgndGroove; BkgndGroove += cDeltaBkgndGroove)
		{
		for(DyadratioThres = cMinDyadratio; DyadratioThres <= cMaxDyadratio; DyadratioThres += cDeltaDyadratio)
			{
			for(Dyad2ratioThres = cMinDyad2ratio; Dyad2ratioThres <= cMaxDyad2ratio; Dyad2ratioThres += cDeltaDyad2ratio)
				{
				for(Dyad3ratioThres = cMinDyad3ratio;Dyad3ratioThres <= cMaxDyad3ratio; Dyad3ratioThres += cDeltaDyad3ratio)
					{
					printf("\b\b\b\b\b\b\b\b%1.8d",++Iteration);

					LocateDyads(PMode,&m_ExprSet,OfsStart,TruncLength,BkgndGroove,DyadratioThres,Dyad2ratioThres,Dyad3ratioThres);
					if(szFastaCtrlFile[0] != '\0')
						LocateDyads(PMode,&m_CtrlSet,OfsStart,TruncLength,BkgndGroove,DyadratioThres,Dyad2ratioThres,Dyad3ratioThres);
					BuffOfs = sprintf(szBuff,"\n%1.4f,%1.4f,%1.4f,%1.4f,%d,%d",
						BkgndGroove,DyadratioThres,Dyad2ratioThres,Dyad3ratioThres,m_ExprSet.NSamples,m_ExprSet.NumDyads);
					if(szFastaCtrlFile[0] != '\0')
						BuffOfs += sprintf(&szBuff[BuffOfs],",%d,%d,%4.4f",
										m_CtrlSet.NSamples,m_CtrlSet.NumDyads,(float)m_ExprSet.NumDyads/m_CtrlSet.NumDyads);
					CUtility::SafeWrite(m_hProfileFile,szBuff,BuffOfs);
					BuffOfs = 0;
					}
				}
			}
		}
	close(m_hProfileFile);
	m_hProfileFile = -1;
 
	Reset();
	return(Rslt);
	}

LocateDyads(PMode,&m_ExprSet,OfsStart,TruncLength,BkgndGroove,DyadratioThres,Dyad2ratioThres,Dyad3ratioThres);
gDiagnostics.DiagOut(eDLInfo,"ProcessFastaStruct","From %d sequences there are %d predicted dyads, %d with no dyads",m_ExprSet.NSamples,m_ExprSet.NumDyads, m_ExprSet.NSamples - m_ExprSet.NumDyads);

#ifdef _WIN32
if((m_hProfileFile = open(pszProfileFile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
if((m_hProfileFile = open(pszProfileFile, O_RDWR | O_CREAT | O_TRUNC, S_IREAD | S_IWRITE ))==-1)
#endif
	{
	gDiagnostics.DiagOut(eDLFatal,"Process","Unable to create or truncate results profile file %s error: %s",pszProfileFile,strerror(errno));
	Reset();
	return(eBSFerrCreateFile);
	}

NumDyads = 0;
BuffOfs = 0;
bFirst = true;
while((pSeqValues = IterNext(&m_ExprSet,bFirst,PMode == ePMrandsel ? true : false))!=NULL)
	{
	bFirst = false;
	if(pSeqValues->DyadIdx == -1)
		continue;
	NumDyads += 1;

	pszDescr = (char *)&pSeqValues->SeqDescr[pSeqValues->SeqLen];
    NumFields = sscanf(pszDescr,"lcl|%d %70[^|]|%70[^|]|%d|%d|%d|%c|%d",&Ident,szSpecies,szChrom,&StartLoci,&EndLoci,&ElLen,&Strand,&NumReads);
	if(NumFields < 7)
		continue;

	if(NumFields == 7)
		NumReads = 1;

	StartLoci += pSeqValues->DyadIdx;
	StartLoci -= 74;
	EndLoci = StartLoci + 147;
	ElLen = 147;
	BuffOfs += sprintf(&szBuff[BuffOfs],"%d,\"%s\",\"%s\",%d,%d,%d,%d,%d\n",Ident,szSpecies,szChrom,StartLoci,EndLoci,ElLen,pSeqValues->Score,NumReads);
	if((BuffOfs + 200) > sizeof(szBuff))
		{	
		CUtility::SafeWrite(m_hProfileFile,szBuff,BuffOfs);
		BuffOfs = 0;
		}
	}
if(BuffOfs)
	CUtility::SafeWrite(m_hProfileFile,szBuff,BuffOfs);
gDiagnostics.DiagOut(eDLInfo,"ProcessFastaStruct","From %d sequences there are %d predicted dyads, %d with no dyads",m_ExprSet.NSamples,NumDyads, m_ExprSet.NSamples - NumDyads);

close(m_hProfileFile);
m_hProfileFile = -1;
 
Reset();
return(Rslt);
}

int 
LocateDyads(etPMode PMode,			// processing mode
			tsSeqSet *pSeqSet,		// interested in dyads for this sequence set
			int OfsStart,			// offset start of sequences by this many bases
			int TruncLength,		// truncate sequences to be no  longer than this length after any offseting
			double BkgndGroove,		// dyads are relative to this background groove
			double DyadratioThres,	// dyad grooves must be at least this ratio to background
			double Dyad2ratioThres,	// immediately flanking grooves must be at least this ratio to background
			double Dyad3ratioThres)	// remainder of flanking grooves must be at least this ration to background
{
int Rslt;
bool bFirst;
int Step;
tsSeqValues *pSeqValues;
int *pTwist;

double BaseLineAv;

double DyadRatio;
double Dyad2Ratio;
double Dyad3Ratio;

int TotNumSteps;
int NumSteps;
int *pValues;
int *pConfGroove;
int *pConfTwist;
int ChkGroove[13];		// to hold dyad (ChkGroove[6]) and +/- 6 at approx decimer (depends on twist) offsets 
int DecIdx;				// index into ChkGroove, incr/decr every 360 degree twist
int GrooveCnt;			// to hold number of groove values contributing to current ChkGroove[DecIdx] so average can be calculated
int AccumTwist;			// to hold accumulated twist relative to dyad
int ChkTwist;			// AccumTwist % 360 used to determine if minor groove back on same plane as at dyad

int ConfValues[cMaxProcSeqLen];
int TwistValues[cMaxProcSeqLen];	

int SeqIdx;
int DyadFirstOfs;
int	DyadLastOfs;
int NumTooTruncated;

if(BkgndGroove == 0.0f)
	{
	// bit slow - iterate each sequence to determine the baseline average
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Determining background baseline mean minor groove...");
	BaseLineAv = 0.0f;
	bFirst = true;
	TotNumSteps = 0;
	while((pSeqValues = IterNext(pSeqSet,bFirst,PMode == ePMrandsel ? true : false))!=NULL)
		{
		bFirst = false;
		if((Rslt = m_pTwister->GetSequenceConformation(eSSminorgroove,	// process for this conformational parameter
						  0,								// initial starting offset (0..n) in pSeq
						  0,				                // number of steps (0 for all) to process starting at pSeq[iStartPsn]|pSeq[iStartPsn+1]
  						  pSeqValues->SeqLen,				// number of nucleotides
						  pSeqValues->SeqDescr,				// sequence to be processed
						  ConfValues))!=eBSFSuccess)		// where to return conformational values
			{
			gDiagnostics.DiagOut(eDLFatal,"LoadSequences","ProcessSequence failed");
			Reset();
			return(Rslt);
			}

		pValues = ConfValues;
		for(Step = 0; Step < (pSeqValues->SeqLen-1); Step++,pValues++)
			BaseLineAv += (double)*pValues/10000.0f;
		TotNumSteps += pSeqValues->SeqLen-1;
		}
	BaseLineAv /= TotNumSteps;
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Background baseline mean minor groove: %1.6f",BaseLineAv);
	}
else
	BaseLineAv = BkgndGroove;


pSeqSet->NumDyads = 0;
NumTooTruncated = 0;
bFirst = true;
while((pSeqValues = IterNext(pSeqSet,bFirst,(PMode == ePMrandsel || PMode == ePMrandselOpt) ? true : false))!=NULL)
	{
	bFirst = false;
	if((Rslt = m_pTwister->GetSequenceConformation(eSSminorgroove,	// process for this conformational parameter
					  0,								// initial starting offset (0..n) in pSeq
					  0,				                // number of steps (0 for all) to process starting at pSeq[iStartPsn]|pSeq[iStartPsn+1]
  					  pSeqValues->SeqLen,				// number of nucleotides
					  pSeqValues->SeqDescr,				// sequence to be processed
					  ConfValues))!=eBSFSuccess)		// where to return conformational values
					{
					gDiagnostics.DiagOut(eDLFatal,"LoadSequences","ProcessSequence failed");
					Reset();
					return(Rslt);
					}

	if((Rslt = m_pTwister->GetSequenceConformation(eSStwist,	// process for this conformational parameter
						  0,									// initial starting offset (0..n) in pSeq
						  0,					                // number of steps (0 for all) to process starting at pSeq[iStartPsn]|pSeq[iStartPsn+1]
  						  pSeqValues->SeqLen,					// number of nucleotides
						  pSeqValues->SeqDescr,					// sequence to be processed
						  TwistValues))!=eBSFSuccess)			// where to return conformational values
						{
						gDiagnostics.DiagOut(eDLFatal,"LoadSequences","ProcessSequence failed");
						Reset();
						return(Rslt);
						}
    NumSteps = pSeqValues->SeqLen - 1;
	pValues = ConfValues;
	pTwist = TwistValues;
	pSeqValues->DyadIdx = -1;
	pSeqValues->Score = 0;
	DyadFirstOfs = OfsStart + 65;				// allows for twist as 6 360 degrees of twist can require 65 nt
	if(TruncLength > 0 && (NumSteps - OfsStart) > TruncLength)
		DyadLastOfs = OfsStart + TruncLength - 65;
	else
		DyadLastOfs = NumSteps - 65;
	if(DyadLastOfs < DyadFirstOfs)
		{
		NumTooTruncated += 1;
		continue;
		}
	memset(ChkGroove,0,sizeof(ChkGroove));
	for(SeqIdx = DyadFirstOfs; SeqIdx <= DyadLastOfs; SeqIdx++)
		{
		DecIdx = 6;
		pConfGroove = &pValues[SeqIdx];
		ChkGroove[DecIdx++] = *pConfGroove++;
		DyadRatio = (double)ChkGroove[6]/(BaseLineAv * 10000.0f);
		if(DyadRatio < DyadratioThres)
			continue;
		pConfTwist =  &pTwist[SeqIdx+1];
		AccumTwist = *pConfTwist;
		ChkGroove[DecIdx] = 0;
		GrooveCnt = 0;
		
		// iterate over bases to right of putative dyad and every rotation of the dsDNA get the minor groove
		int Bases = 1;
		while(DecIdx <= 12)
			{
			Bases += 1;
			pConfTwist += 1;
			pConfGroove += 1;
			AccumTwist += *pConfTwist;
			ChkTwist = AccumTwist % 3600000;
			if(ChkTwist >= 3300000 || ChkTwist <= 300000)
				{
				ChkGroove[DecIdx] += *pConfGroove;
				GrooveCnt += 1;
				}
			else
				{
				if(GrooveCnt > 0)
					{
					ChkGroove[DecIdx] /= GrooveCnt;
					GrooveCnt = 0;
					if(DecIdx++ < 12)
						ChkGroove[DecIdx]= 0;
					}
				}
			}
		// now iterate over bases to left of putative dyad and every rotation of the dsDNA get the minor groove
		DecIdx = 5;
		pConfGroove =  &pValues[SeqIdx-1];
		pConfTwist =  &pTwist[SeqIdx-1];
		AccumTwist = *pConfTwist;
		ChkGroove[DecIdx] = 0;
		GrooveCnt = 0;
		Bases = 1;
		while(DecIdx >= 0)
			{
			Bases += 1;
			pConfTwist -= 1;
			pConfGroove -= 1;
			AccumTwist += *pConfTwist;
			ChkTwist = AccumTwist % 3600000;
			if(ChkTwist >= 3300000 || ChkTwist <= 300000)
				{
				ChkGroove[DecIdx] += *pConfGroove;
				GrooveCnt += 1;
				}
			else
				{
				if(GrooveCnt > 0)
					{
					ChkGroove[DecIdx] /= GrooveCnt;
					GrooveCnt = 0;
					if(DecIdx-- > 0)
						ChkGroove[DecIdx] = 0;
					}
				}
			}

		Dyad2Ratio = (double)(ChkGroove[5] + ChkGroove[7])/(2*BaseLineAv*10000.0f);
		Dyad3Ratio = (double)(ChkGroove[0] + ChkGroove[1] + ChkGroove[2] + ChkGroove[3] + ChkGroove[4] +
					ChkGroove[8] + ChkGroove[9] + ChkGroove[10] + ChkGroove[11] + ChkGroove[12])/(10*BaseLineAv*10000.0f);

		
		if(Dyad2Ratio < Dyad2ratioThres || Dyad3Ratio < Dyad3ratioThres)
			continue;

			// how to really score these dyads? One day real scores will be associated...
		int Score = (int)(1000 * ((DyadRatio - 1.0f) + ((Dyad2Ratio - 1.0f) * 0.85) + ((Dyad3Ratio - 1.0f) * 0.75)));
		if(Score > pSeqValues->Score)
			{
			pSeqValues->Score = Score;
			if(pSeqValues->DyadIdx == -1)
				pSeqSet->NumDyads += 1;
			pSeqValues->DyadIdx = SeqIdx;
			}
		}
	}

if(NumTooTruncated)
	gDiagnostics.DiagOut(eDLInfo,"Process","Sloughed %d as after offseting and truncating these sequences were too short to process",NumTooTruncated);
return(pSeqSet->NumDyads);
}

int
LoadSequences(etPMode PMode,			// processing mode
  			  tsSeqSet *pSeqSet,		// load into this sequence set
	  		  int NSamples,				// sample this number of sequences if non-zero
			  char *pszFastaFile)		// fasta reads
{
int Rslt;
tsSeqValues *pSeqValues;
tsSeqValues *pTmpSeqValues;

int ElSize;
etSeqBase Sequence[cMaxProcSeqLen + 10];
  
int SeqLen;
int NumSteps;
char szDescr[cMaxFastaDescrLen];
int DescrLen;
bool bDescriptor;
bool bTrunc;

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Loading sequences from multifasta file: '%s'",pszFastaFile);
strncpy(pSeqSet->szFastaFile,pszFastaFile,sizeof(pSeqSet->szFastaFile));
pSeqSet->szFastaFile[sizeof(pSeqSet->szFastaFile)-1] = '\0';
if((m_pFasta = new CFasta())==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,"LoadSequences","Unable to create CFasta object");
	Reset();
	return(eBSFerrObj);
	}

if((Rslt = m_pFasta->Open(pszFastaFile,true)) != eBSFSuccess)
	{
	while(m_pFasta->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,m_pFasta->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,"LoadSequences","Unable to open fasta file '%s'",pszFastaFile);
	Reset();
	return(Rslt);
	}


gDiagnostics.DiagOut(eDLInfo,gszProcName,"Allocating memory for sequences");

if((pSeqSet->pSeqValues = (tsSeqValues *)malloc(cAllocSeqSize))==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate memory (%d bytes) for holding sequences",cAllocSeqSize);
	Reset();
	return(eBSFerrMem);
	}
pSeqSet->AllocdMemSeqValues = cAllocSeqSize;
pSeqValues = pSeqSet->pSeqValues;

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Loading sequences");

SeqLen= 0;
NumSteps = 0;
bTrunc = false;
while((Rslt = SeqLen = m_pFasta->ReadSequence(Sequence,sizeof(Sequence))) > eBSFSuccess)
	{
	if(SeqLen == eBSFFastaDescr)		// just read a descriptor line
		{
		m_pFasta->ReadDescriptor(szDescr,sizeof(szDescr));
		szDescr[cTruncFastaDescrLen-1] = '\0';
		pSeqSet->NumSeqsProc += 1;
		bDescriptor = true;
		bTrunc = false;
		continue;
		}
	if(!bDescriptor)					// if never seen a descriptor then dummy up one...
		{
		sprintf(szDescr,"Probe1");
		bDescriptor = true;
		pSeqSet->NumSeqsProc += 1;
		bTrunc = false;
		}
	if(bTrunc)
		continue;

	// sequence long enough to process?
	if(SeqLen < cMinProfileLen)
		{
		if(!pSeqSet->NumShortSeqs++)
			gDiagnostics.DiagOut(eDLInfo,"LoadSequences","Skipping at least one sequence shorter than %d, SeqLen = %d",cMinProfileLen,SeqLen);
		continue;
		}

	if(SeqLen > cMaxProcSeqLen)
		{
		if(!pSeqSet->NumTruncSeqs++)
			gDiagnostics.DiagOut(eDLInfo,"LoadSequences","Truncating at least one long sequence longer than %d",cMaxProcSeqLen);
		NumSteps = cMaxProcSeqLen - 1;
		bTrunc = true;
		}
	else
		NumSteps = SeqLen - 1;

	// remove any repeat masking and randomly substitute bases for eBaseN's - not expecting too many of these say's he hopefully!
	etSeqBase *pSeq = Sequence;
	for(int SeqIdx = 0; SeqIdx < SeqLen; SeqIdx++,pSeq++)
		if((*pSeq &= ~cRptMskFlg) > eBaseT)
			*pSeq = rand() % 4;


	DescrLen = (int)strlen(szDescr);
	ElSize = sizeof(tsSeqValues) + SeqLen + DescrLen;	
	if((pSeqSet->UsedMemSeqValues + ElSize) > pSeqSet->AllocdMemSeqValues)
		{
		if((pTmpSeqValues = (tsSeqValues *)realloc(pSeqSet->pSeqValues,pSeqSet->AllocdMemSeqValues + cAllocSeqSize))==NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate memory (%d bytes) for holding sequence conformational values",pSeqSet->AllocdMemSeqValues + cAllocSeqSize);
			Reset();
			return(eBSFerrMem);
			}
		pSeqSet->pSeqValues = pTmpSeqValues;
		pSeqSet->AllocdMemSeqValues += cAllocSeqSize;
		}
	pSeqValues = (tsSeqValues *)((char *)pSeqSet->pSeqValues + pSeqSet->UsedMemSeqValues);
	pSeqValues->SeqID = ++pSeqSet->NumSeqValues;
	pSeqValues->Sizeof = ElSize;
	pSeqValues->LenDescr = DescrLen; 
	pSeqValues->SeqLen = SeqLen;
	pSeqValues->DyadIdx = -1;
	pSeqValues->Score = 0;
	pSeqValues->RandSel = 0;
	memcpy(pSeqValues->SeqDescr,Sequence,SeqLen);
	strcpy((char *)&pSeqValues->SeqDescr[SeqLen],szDescr);
	pSeqSet->UsedMemSeqValues += ElSize;
	
	if((PMode == ePMNsel || PMode == ePMNselOpt) && NSamples == pSeqSet->NumSeqValues)
		break;
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Sequences (%d) loaded ready for profile alignment processing...",pSeqSet->NumSeqValues);
delete m_pFasta;
m_pFasta = NULL;
delete m_pTwister;
m_pTwister = NULL;
return(pSeqSet->NumSeqValues);
}

// MarkRandoms
// Randomly selects NSamples sequences and marks these as having been randomly selected
// Note that any previously marked sequences will be unmarked if bReset set true
int
MarkRandoms(tsSeqSet *pSeqSet,int NSamples,bool bReset)
{
bool bFirst;
tsSeqValues *pSeqValues;
tsSeqValues **ppIndex;
int Idx;

if((ppIndex = new tsSeqValues *[pSeqSet->NumSeqValues])==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate memory (%d bytes) for holding sequence index",NSamples * sizeof(tsSeqValues *));
	Reset();
	return(eBSFerrMem);
	}

Idx = 0;
bFirst = true;
while((pSeqValues=IterNext(pSeqSet,bFirst,false))!=NULL)
	{
	bFirst = false;
	if(bReset)
		pSeqValues->RandSel = 0;
	ppIndex[Idx++] = pSeqValues;
	}

TRandomCombined<CRandomMother,CRandomMersenne> RG((int)time(0));

while(NSamples--)
	{
	Idx = (int)RG.IRandom(0,pSeqSet->NumSeqValues-1);
	ppIndex[Idx]->RandSel = 1;
	}

delete ppIndex;
return(eBSFSuccess);
}

// IterNext
// returns next sequence values after current, or NULL if last sequence previously iterated
tsSeqValues *
IterNext(tsSeqSet *pCtrlSet,	// iterate this set of sequences
		 bool bFirst,		// true if first sequence to be iterated, subsequent calls must have bFirst set false
		 bool bRandSels)	// true if only sequences marked as RandSel to be returned
{
if(bFirst || pCtrlSet->pIterSeqValues == NULL)
	{
	pCtrlSet->pIterSeqValues = pCtrlSet->pSeqValues;
	if(bRandSels)
		{
		while(pCtrlSet->pIterSeqValues->SeqID != pCtrlSet->NumSeqValues && !pCtrlSet->pIterSeqValues->RandSel)
			pCtrlSet->pIterSeqValues = (tsSeqValues *)((char *)pCtrlSet->pIterSeqValues + pCtrlSet->pIterSeqValues->Sizeof);
		if(!pCtrlSet->pIterSeqValues->RandSel)
			return(NULL);
		}
	return(pCtrlSet->pIterSeqValues);
	}
if(pCtrlSet->pIterSeqValues->SeqID == pCtrlSet->NumSeqValues)
	return(NULL);
pCtrlSet->pIterSeqValues = (tsSeqValues *)((char *)pCtrlSet->pIterSeqValues + pCtrlSet->pIterSeqValues->Sizeof);
if(bRandSels)
	{
	while(pCtrlSet->pIterSeqValues->SeqID != pCtrlSet->NumSeqValues && !pCtrlSet->pIterSeqValues->RandSel)
		pCtrlSet->pIterSeqValues = (tsSeqValues *)((char *)pCtrlSet->pIterSeqValues + pCtrlSet->pIterSeqValues->Sizeof);
	if(!pCtrlSet->pIterSeqValues->RandSel)
		return(NULL);
	}
return(pCtrlSet->pIterSeqValues);
}