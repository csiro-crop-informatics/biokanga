// DNAseqSitePotential.cpp : Defines the entry point for the console application.
// Used to generate genome wide read start site potentials
#include "stdafx.h"
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#if _WIN32

#include "../libbiokanga/commhdrs.h"
#else
#include "../libbiokanga/commhdrs.h"
#endif

const char *cpszProgVer = "1.0.2";		// increment with each release

// processing modes
typedef enum TAG_ePMode {		
	ePMdefault,					// default processing mode
	ePMplaceholder				// used to set the enumeration range
	} etPMode;

typedef enum TAG_eReadsSortMode {
		eRSMReadID,				// index by ascending ReadID
		eRSMHitMatch,			// index by ascending chrom, loci, strand, level
		eRSMplaceholder			// used to limit the enumeration range
} etReadsSortMode;

#pragma pack(1)
typedef struct TAG_sAlignHit {
	UINT32 AlignHitIdx;			// current read hit index + 1 for this read
	char szChromName[cMaxDatasetSpeciesChrom];	// identifies hit chromosome
	UINT8 Strand;				// hit strand - '+' or '-'
	UINT32 Loci;				// offset on chromosome of hit
	UINT16 MatchLen;			// match length
	etSeqBase Seq[8];			// will hold 4nt 5' to read start and 4nt 3' to read start
} tsAlignHit;
#pragma pack()

int
Process(etPMode PMode,					// processing mode
		char Strand,					// process for this strand '+' or '-' or for both '*'
		char *pszReadsFile,				// input from this alignreads generated CSV file
		char *pszGenomeFile,			// input from this genome bioseq file
		char *pszPotentialsFile);		// output read start site potentials to this file

tsAlignHit *IterSortedReads(tsAlignHit *pCurAlignHit, bool bForward); //iterate over sorted reads
int SortAlignHits(void);		// index read hits

static int SortHitMatch(const void *arg1, const void *arg2);

CStopWatch gStopWatch;
CDiagnostics gDiagnostics;				// for writing diagnostics messages to log file
char gszProcName[_MAX_FNAME];			// process name


const int cReadsInitalAlloc   = 10000000;						 // initial allocation to hold this many reads
const int cReadsReAlloc = 5000000;	     // realloc allocation in this sized increments
tsAlignHit *m_pAlignHits = NULL;	// memory allocated to hold reads, reads are written contiguously into this memory
UINT32 m_AllocdAlignHits = 0;		// how instances of tsAlignHit have been allocated
UINT32 m_NumAlignHits = 0;			// m_pAlignHits contains this many reads

tsAlignHit **m_ppAlignHitsIdx = NULL;	// memory allocated to hold array of ptrs to read hits in m_pAlignHits - usually sorted by some critera
UINT32 m_AllocdAlignHitsIdx = 0;		// how many elements for m_pAlignHitsIdx have been allocated
etReadsSortMode	m_CurReadsSortMode;	    // sort mode last used on m_ppAlignHitsIdx

CBioSeqFile *m_pBioseq;

etSeqBase *m_pSeq;						// to hold chromosome sequences
int m_AllocdSeqLen;						// allocated to m_pSeq

int *m_pGenomeOct;						// to hold genome octamer cnts
int *m_pSiteOct;						// to hold read start site octamer cnts

int m_hRsltFile;						// to hold results file handle

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
	return _T("DNAseqSitePotential");
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

etPMode PMode;				// processing mode
int Strand;					// process for this strand '+' or '-' or for both '*'
char szReadsFile[_MAX_PATH];	// alignreads generated CSV file
char szGenomeFile[_MAX_PATH];	// genome bioseq file
char szPotentialFile[_MAX_PATH];	// output read start site potentials to this file


// command line args
struct arg_lit  *help    = arg_lit0("hH","help",                "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"diagnostics log file");

struct arg_int *pmode = arg_int0("m","mode","<int>",		    "processing mode: 0 - output read start site potentials (default = 0)");
struct arg_str *strand=arg_str0("s", "strand","<str>",          "process for this strand '+' or '-' only (default is both)");
struct arg_file *infile = arg_file1("i","in","<file>",			"input from this alignreads generated CSV file");
struct arg_file *ingfile = arg_file1("I","in","<file>",			"input from this genome bioseq file");
struct arg_file *outfile = arg_file1("o","out","<file>",		"output read start potentials to this file");

struct arg_end *end = arg_end(20);

void *argtable[] = {help,version,FileLogLevel,LogFile,
					pmode,strand,infile,ingfile,outfile,
					end};

char **pAllArgs;
int argerrors;
argerrors = CUtility::arg_parsefromfile(argc,(char **)argv,&pAllArgs);
if(argerrors >= 0)
	argerrors = arg_parse(argerrors,pAllArgs,argtable);

if (help->count > 0)
        {
		printf("\n%s Whole genome Read Start Site Potentials, Version %s\nOptions ---\n", gszProcName,cpszProgVer);
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

	strcpy(szPotentialFile,outfile->filename[0]);
	strcpy(szReadsFile,infile->filename[0]);
	strcpy(szGenomeFile,ingfile->filename[0]);

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

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Version: %s Processing parameters:",cpszProgVer);
	const char *pszDescr;
	switch(PMode) {
		case ePMdefault:				
			pszDescr = "Default - generate read start site potentials";
			break;
		}

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"processing mode is : '%s'",pszDescr);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Process for this strand : '%c'",Strand);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"input from this alignreads generated CSV file: '%s'",szReadsFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"input from this genome bioseq file: '%s'",szGenomeFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"output read start site potentials to this file: '%s'",szPotentialFile);

#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	gStopWatch.Start();
	Rslt = Process(PMode,Strand,szReadsFile,szGenomeFile,szPotentialFile);
	gStopWatch.Stop();
	Rslt = Rslt >=0 ? 0 : 1;
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Exit code: %d Total processing time: %s",Rslt,gStopWatch.Read());
	exit(Rslt);
	}
else
	{
	printf("\n%s Predict whole genome read start site potentials from aligned reads, Version %s\n",gszProcName,cpszProgVer);
	arg_print_errors(stdout,end,gszProcName);
	arg_print_syntax(stdout,argtable,"\nUse '-h' to view option and parameter usage\n");
	exit(1);
	}
return(0);
}

void
Init(void)
{
m_pAlignHits = NULL;
m_ppAlignHitsIdx = NULL;
m_pBioseq = NULL;
m_pSeq = NULL;
m_AllocdSeqLen = 0;
m_AllocdAlignHitsIdx = 0;
m_AllocdAlignHits = 0;
m_NumAlignHits = 0;
m_hRsltFile = -1;
}

void
Reset(void)
{
if(m_pAlignHits != NULL)
	{
	delete m_pAlignHits;
	m_pAlignHits = NULL;
	}
if(m_ppAlignHitsIdx != NULL)
	{
	delete m_ppAlignHitsIdx;
	m_ppAlignHitsIdx = NULL;
	}
if(m_pBioseq != NULL)
	{
	delete m_pBioseq;
	m_pBioseq = NULL;
	}
if(m_pSeq != NULL)
	{
	delete m_pSeq;
	m_pSeq = NULL;
	}
if(m_pGenomeOct != NULL)
	{
	delete m_pGenomeOct;
	m_pGenomeOct = NULL;
	}
if(m_pSiteOct != NULL)
	{
	delete m_pSiteOct;
	m_pSiteOct = NULL;
	}
if(m_hRsltFile != -1)
	{
	close(m_hRsltFile);
	m_hRsltFile = -1;
	}
m_AllocdSeqLen = 0;
m_AllocdAlignHitsIdx = 0;
m_AllocdAlignHits = 0;
m_NumAlignHits = 0;
}

int
LoadReadsBED(char *pszInFile,
			char FiltStrand)
{
int Rslt;
CBEDfile *pBEDFile;
int CurFeatureID;
int StartLoci;
int EndLoci;
int Score;
char szChrom[128];
char szFeatName[128];
char Strand;
tsAlignHit *pAlignHit;

if((pBEDFile = new CBEDfile)==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CBEDfile");
	return(eBSFerrObj);
	}

if((Rslt=pBEDFile->Open(pszInFile,eBTAnyBed)) !=eBSFSuccess)
	{
	while(pBEDFile->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,pBEDFile->GetErrMsg());
	delete pBEDFile;
	return(eBSFerrOpnFile);
	}
CurFeatureID = 0;
while(Rslt == eBSFSuccess && (CurFeatureID = pBEDFile->GetNextFeatureID(CurFeatureID)) > 0)
	{
	pBEDFile->GetFeature(CurFeatureID,// feature instance identifier
				szFeatName,				// where to return feature name
				szChrom,				// where to return chromosome name
				&StartLoci,				// where to return feature start on chromosome (0..n) 
				&EndLoci,				// where to return feature end on chromosome
 				&Score,					// where to return score
 				&Strand);				// where to return strand

	if(FiltStrand != '*' && Strand != FiltStrand)
		continue;

			// need to allocate more memory?
	if(m_NumAlignHits >= m_AllocdAlignHits)
		{
		pAlignHit = (tsAlignHit *) new tsAlignHit [m_AllocdAlignHits + cReadsReAlloc];	
		if(pAlignHit == NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadReads: Memory re-allocation to %d bytes - %s",(m_AllocdAlignHits + cReadsReAlloc) * sizeof(tsAlignHit),strerror(errno));
			delete pBEDFile;
			return(eBSFerrMem);
			}

		memcpy(pAlignHit,m_pAlignHits,m_NumAlignHits*sizeof(tsAlignHit));
		delete m_pAlignHits;
		m_pAlignHits = pAlignHit;
		m_AllocdAlignHits += cReadsReAlloc;
		}

	pAlignHit = &m_pAlignHits[m_NumAlignHits++];
	memset(pAlignHit,0,sizeof(tsAlignHit));
	strcpy(pAlignHit->szChromName,szChrom);
	pAlignHit->Loci = StartLoci;
	pAlignHit->MatchLen = 1 + EndLoci - StartLoci;
	pAlignHit->Strand = Strand;
	}
delete pBEDFile;
return(Rslt);
}




int
LoadReadsCSV(char *pszInFile,			// load reads from this file
		  char Strand)			    // process for this strand '+' or '-' or for both '*'
{
int Rslt;
int NumProcessed;
int ReadID;
int MatchLen;
char *pszTargSpecies;
char *pszChromName;
int Loci;
char *pszStrand;
int NumExclReads;

int NumFields;

tsAlignHit *pAlignHit;					// current read hit
int BuffLen = 0;
int BuffOfs = 0;

m_NumAlignHits = 0;
CCSVFile *pCSVAligns = NULL;


// load into memory
if((pCSVAligns = new CCSVFile) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CCSVfile");
	return(eBSFerrObj);
	}
pCSVAligns->SetMaxFields(14);	// expecting at least 12, upto 14, fields in reads alignment CSV file
if((Rslt=pCSVAligns->Open(pszInFile))!=eBSFSuccess)
	{
	while(pCSVAligns->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,pCSVAligns->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open file: %s",pszInFile);
	delete pCSVAligns;
	return(Rslt);
	}

NumExclReads = 0;
NumProcessed = 0;
while((Rslt=pCSVAligns->NextLine()) > 0)	// onto next line containing fields
	{
	if(!(m_NumAlignHits % 100000))
		{
		if(!m_NumAlignHits)
			printf("\n     processing match %8.8d",m_NumAlignHits);
		else
			printf("\b\b\b\b\b\b\b\b%8.8d",m_NumAlignHits);
		}

	NumFields = pCSVAligns->GetCurFields();
	if(NumFields < 8)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"file: %s contains % fields, expected at least 8",pszInFile,NumFields);
		delete pCSVAligns;
		return(eBSFerrParams);
		}

	if(!NumProcessed && pCSVAligns->IsLikelyHeaderLine())
		continue;
	NumProcessed += 1;
	pCSVAligns->GetInt(1,&ReadID);
	pCSVAligns->GetText(3,&pszTargSpecies);
	pCSVAligns->GetText(4,&pszChromName);
	pCSVAligns->GetInt(5,&Loci);
	pCSVAligns->GetInt(7,&MatchLen);
	pCSVAligns->GetText(8,&pszStrand);
	if(Strand != '*' && Strand == *pszStrand)
		continue;

		// need to allocate more memory?
	if(m_NumAlignHits >= m_AllocdAlignHits)
		{
		pAlignHit = (tsAlignHit *) new tsAlignHit [m_AllocdAlignHits + cReadsReAlloc];	
		if(pAlignHit == NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadReads: Memory re-allocation to %d bytes - %s",(m_AllocdAlignHits + cReadsReAlloc) * sizeof(tsAlignHit),strerror(errno));
			delete pCSVAligns;
			return(eBSFerrMem);
			}

		memcpy(pAlignHit,m_pAlignHits,m_NumAlignHits*sizeof(tsAlignHit));
		delete m_pAlignHits;
		m_pAlignHits = pAlignHit;
		m_AllocdAlignHits += cReadsReAlloc;
		}

	pAlignHit = &m_pAlignHits[m_NumAlignHits++];
	memset(pAlignHit,0,sizeof(tsAlignHit));
	strcpy(pAlignHit->szChromName,pszChromName);
	pAlignHit->Loci = Loci;
	pAlignHit->MatchLen = MatchLen;
	pAlignHit->Strand = *pszStrand;
	}
delete pCSVAligns;
return(Rslt);
}

int
LoadReads(char *pszInFile,			// load reads from this file
		  char Strand)			    // process for this strand '+' or '-' or for both '*'
{
int Rslt;

if(m_pAlignHits == NULL)		// initial allocation required?
	{
	m_pAlignHits = (tsAlignHit *) new tsAlignHit [cReadsInitalAlloc];	// initial allocation
	if(m_pAlignHits == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadReads: Memory allocation of %d bytes - %s",cReadsInitalAlloc*sizeof(tsAlignHit),strerror(errno));
		return(eBSFerrMem);
		}
	m_AllocdAlignHits = cReadsInitalAlloc;
	m_NumAlignHits = 0;
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Loading read alignments from file: %s",pszInFile);

// first try loading aligned reads as BED, if that fails then try as CSV
if((Rslt = LoadReadsBED(pszInFile,Strand)) == eBSFerrOpnFile)
	Rslt = LoadReadsCSV(pszInFile,Strand);

if(Rslt >= 0)
	{
	printf("\b\b\b\b\b\b\b\b%8.8d",m_NumAlignHits);
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"LoadReads: accepted %d aligned reads on strand '%c'",m_NumAlignHits,Strand);

	// finally, create sorted index by chrom, loci, strand over loaded reads
	Rslt = SortAlignHits();
	}

return(Rslt);
}

int
GenSeqIdx(int SeqLen,etSeqBase *pSeq)
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
Process(etPMode PMode,					// processing mode
		 char Strand,				    // process for this strand '+' or '-' or for both '*'
		char *pszReadsFile,				// input from this alignreads generated CSV file
		char *pszGenomeFile,			// input from this genome bioseq file
		char *pszPotentialFile)				// output read start site potentials to this file
{
int Rslt;
int Ofs;
int ChromID;
int SeqLen;
int SeqIdx;
int OctIdx;
etSeqBase *pSeq;
tsAlignHit *pRead;
char szCurChrom[cMaxDatasetSpeciesChrom];

if((m_pGenomeOct = new int [0x10000])==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to alloc memory (%d requested) for holding genome octamer cnts",0x0ffff);
	Reset();
	return(eBSFerrMem);
	}

if((m_pSiteOct = new int [0x10000])==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to alloc memory (%d requested) for holding read start site octamer cnts",0x0ffff);
	Reset();
	return(eBSFerrMem);
	}

// load the read loci
Init();
if((Rslt = LoadReads(pszReadsFile,Strand)) < eBSFSuccess)
	{
	Reset();
	return(Rslt);
	}

if((m_pBioseq = new CBioSeqFile) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CBioSeqFile object");
	Reset();
	return(eBSFerrObj);
	}
if((Rslt = m_pBioseq->Open(pszGenomeFile))!=eBSFSuccess)
	{
	while(m_pBioseq->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,m_pBioseq->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open assembly sequence file '%s'",pszGenomeFile);
	Reset();
	return(Rslt);
	}

// next load the sequences bracketing the read start cut sites
memset(m_pSiteOct,0,sizeof(int) * 0x010000);
pRead = NULL;
while((pRead = IterSortedReads(pRead,true))!=NULL)
	{
	if((Rslt= ChromID = m_pBioseq->LocateEntryIDbyName(pRead->szChromName))<=0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to locate chrom '%s' in assembly file '%s'",pRead->szChromName,pszGenomeFile);
		Reset();
		return(Rslt);
		}
	SeqLen = m_pBioseq->GetDataLen(ChromID);

	if(pRead->Strand == '-')
		Ofs = (int)(pRead->Loci + pRead->MatchLen) - 4;
	else
		Ofs = (int)pRead->Loci - 4;
	if(Ofs < 0 || (Ofs + 8) >= SeqLen)
		{
		pRead->MatchLen = 0;
		continue;
		}
	
	if((Rslt=m_pBioseq->GetData(ChromID,eSeqBaseType,Ofs,pRead->Seq,8)) != 8)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Loading sequence failed from chrom: %s loci %d len: %d file: '%s'",pRead->szChromName,Ofs,8,pszGenomeFile);
		Reset();
		return(Rslt);
		}
//	if(pRead->Strand == '-')
//		CSeqTrans::ReverseComplement((unsigned int)8,pRead->Seq);
	OctIdx = GenSeqIdx(8,pRead->Seq);
	if(OctIdx > 0x0ffff)
		printf("Problem!");
	if(OctIdx >= 0)
		m_pSiteOct[OctIdx] += 1;
	}


// now get the genome wide distribution of octamers
memset(m_pGenomeOct,0,sizeof(int) * 0x010000);
ChromID = 0;
while((ChromID = m_pBioseq->Next(ChromID))>0)
	{
	m_pBioseq->GetName(ChromID,sizeof(szCurChrom),szCurChrom);
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing %s genome wide background octamer counts...",szCurChrom);
	SeqLen = m_pBioseq->GetDataLen(ChromID);
	if(m_pSeq == NULL || SeqLen > m_AllocdSeqLen)
		{
		if(m_pSeq != NULL)
			{
			delete m_pSeq;
			m_pSeq = NULL;
			}
		if((m_pSeq = new unsigned char [SeqLen+0x07fff])==NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to alloc memory (%d requested) for holding raw sequence data",SeqLen+0x07fff);
			Reset();
			return(eBSFerrMem);
			}
		m_AllocdSeqLen = SeqLen+0x07fff;
		}

	if((Rslt=m_pBioseq->GetData(ChromID,eSeqBaseType,0,m_pSeq,SeqLen)) != SeqLen)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Loading sequence failed from chrom: %s loci %d len: %d file: '%s'",pRead->szChromName,Ofs,8,pszGenomeFile);
		Reset();
		return(Rslt);
		}

	pSeq = m_pSeq;
	for(SeqIdx = 0; SeqIdx < SeqLen - 8; SeqIdx++)
		{
		OctIdx = GenSeqIdx(8,pSeq++);
		if(OctIdx > 0x0ffff)
			printf("Problem!");

		if(OctIdx >= 0)
			m_pGenomeOct[OctIdx] += 1;
		}
	}

// for now, we just output the raw counts
char szBuff[4096];
int BuffIdx;

#ifdef _WIN32
if((m_hRsltFile = open(pszPotentialFile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
if((m_hRsltFile = open(pszPotentialFile, O_RDWR | O_CREAT |O_TRUNC,S_IREAD | S_IWRITE))==-1)
#endif
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create %s - %s",pszPotentialFile,strerror(errno));
	Reset();
	return(eBSFerrCreateFile);
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Output file created/truncated: '%s'",pszPotentialFile);
BuffIdx = 0;
double Ratio;

for(OctIdx = 0; OctIdx <= 0x0ffff; OctIdx++)
	{
	if(m_pGenomeOct[OctIdx] == 0)
		{
		Ratio = 0.0;
		if(m_pSiteOct[OctIdx] != 0)
			printf("\nWhy does m_pSiteOct[%d] == %d",OctIdx,m_pSiteOct[OctIdx]);
		}
	else
		Ratio = (double)m_pSiteOct[OctIdx]/m_pGenomeOct[OctIdx];
	BuffIdx += sprintf(&szBuff[BuffIdx],"\"%s\",%d,%d,%8.8f\n",StepIdx2Seq(8,OctIdx),m_pGenomeOct[OctIdx],m_pSiteOct[OctIdx],Ratio);
	if(BuffIdx + 100 > sizeof(szBuff))
		{
		CUtility::SafeWrite(m_hRsltFile,szBuff,BuffIdx);
		BuffIdx = 0;
		}
	}
if(BuffIdx)
	CUtility::SafeWrite(m_hRsltFile,szBuff,BuffIdx);
Reset();
return(Rslt);
}


// IterSortedReads
// use to iterate over sorted reads returning ptr to -
// bForward == true : next read following the current read
// bForward == false : read previous to the current read
// To start from first (bForward == true) read then pass in NULL as pCurAlignHit
// To start from last (bForward == false) read then pass in NULL as pCurAlignHit
// Returns NULL if all read hits have been iterated
tsAlignHit *
IterSortedReads(tsAlignHit *pCurAlignHit, bool bForward)
{
tsAlignHit *pNxtAlignHit = NULL;
if(bForward)
	{
	if(pCurAlignHit == NULL)
		pNxtAlignHit = m_ppAlignHitsIdx[0];
	else
		if(pCurAlignHit->AlignHitIdx < m_NumAlignHits)
			pNxtAlignHit = m_ppAlignHitsIdx[pCurAlignHit->AlignHitIdx];
	}
else
	{
	if(pCurAlignHit == NULL)
		pNxtAlignHit = m_ppAlignHitsIdx[m_NumAlignHits-1];
	else
		if(pCurAlignHit->AlignHitIdx > 1)
			pNxtAlignHit = m_ppAlignHitsIdx[pCurAlignHit->AlignHitIdx-2];
	}
return(pNxtAlignHit);
}

int
SortAlignHits(void)
{
tsAlignHit *pAlignHit;
int Idx;

if(m_ppAlignHitsIdx != NULL && m_AllocdAlignHitsIdx >= m_NumAlignHits)
	return(eBSFSuccess);

if(m_ppAlignHitsIdx == NULL || m_AllocdAlignHitsIdx < m_NumAlignHits)
	{
	if(m_ppAlignHitsIdx != NULL)
		{
		delete m_ppAlignHitsIdx;
		m_ppAlignHitsIdx = NULL;
		m_AllocdAlignHitsIdx = 0;
		}
	if((m_ppAlignHitsIdx = (tsAlignHit **) new tsAlignHit * [m_NumAlignHits])==NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"SortAlignHits: Memory reads index allocation for %d ptrs - %s",m_NumAlignHits,strerror(errno));
		return(eBSFerrMem);
		}
	m_AllocdAlignHitsIdx = m_NumAlignHits;
	}

pAlignHit = NULL;
tsAlignHit **pIdx = m_ppAlignHitsIdx;

for(Idx = 0; Idx < (int)m_NumAlignHits; Idx++)
	*pIdx++ = &m_pAlignHits[Idx];

qsort(m_ppAlignHitsIdx,m_NumAlignHits,sizeof(tsAlignHit *),SortHitMatch);

for(Idx = 0; Idx < (int)m_NumAlignHits; Idx++)
	m_ppAlignHitsIdx[Idx]->AlignHitIdx = Idx + 1;
return(eBSFSuccess);
}


// SortHitmatch
// Sort by ascending chrom, loci, strand
static int
SortHitMatch(const void *arg1, const void *arg2)
{
tsAlignHit *pEl1 = *(tsAlignHit **)arg1;
tsAlignHit *pEl2 = *(tsAlignHit **)arg2;
int Rslt;

if(pEl1->szChromName[0] < pEl2->szChromName[0] )
		return(-1);
if(pEl1->szChromName[0] > pEl2->szChromName[0] )
	return(1);
if((Rslt = stricmp(pEl1->szChromName,pEl2->szChromName))!=0)
	return(Rslt);

if(pEl1->Loci < pEl2->Loci )
		return(-1);
if(pEl1->Loci > pEl2->Loci )
	return(1);
if(pEl1->Strand < pEl2->Strand )
		return(-1);
if(pEl1->Strand > pEl2->Strand )
	return(1);

return(0);
}

