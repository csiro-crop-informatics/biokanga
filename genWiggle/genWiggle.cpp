// genWiggle.cpp : Defines the entry point for the console application.
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

const char *cpszProgVer = "1.2.0";		// increment with each release
const int cChromSeqReAlloc = 5000000;	// realloc chrom sequence size

// processing modes
typedef enum TAG_ePMode {		
	ePMdefault,					// CSV loci
	ePMBed,						// UCSC BED format
	ePMNase,					// MNase digestion
	ePMconf,					// Conformational characteristics
	ePMhamming,					// Hamming distances
	ePMplaceholder				// used to set the enumeration range
	} etPMode;

typedef enum eBEDRegion {
	eMEGRAny = 0,				// process any region
	eMEGRIntergenic,			// only process intergenic
	eMEGRExons,					// only process exons
	eMEGRIntrons,				// only process introns
	eMEGRCDS,					// only process CDSs
	eMEGUTR,					// only process UTRs
	eMEG5UTR,					// only process 5'UTRs
	eMEG3UTR					// only process 3'UTRs
} etBEDRegion;

#pragma pack(1)
// Hamming specific structures
const int cMaxHammingChroms = 200;	// can handle at most this many chromosomes with hammings
typedef struct TAG_sHamChrom {
	UINT32 ChromID;					// uniquely identifies this chromosome
	UINT8  szChrom[cMaxDatasetSpeciesChrom];	// chrom name
	UINT32 NumEls;					// number of subsequences with hammings on this chrom
	UINT8 Dists[1];					// array, in ascending loci order, of hamming distances
} tsHamChrom;

typedef struct TAG_sHamHdr {
	UINT8 Magic[4];		        // magic chars 'bham' to identify this file as a biosequence file containing hamming edit distances
	UINT32 Version;				// structure version 
	INT32 Len;					// file length, also allocation size required when loading hammings into memory
	UINT16 NumChroms;		    // number of chromosomes with Hammings
	UINT32 ChromOfs[cMaxHammingChroms];	// offsets to each chromosomes respective tsHamChrom
} tsHamHdr;

#pragma pack()

int
Process(etPMode PMode,					// processing mode
		etBEDRegion Region,				// regions of interest if processing BED files
		char Strand,					// process for this strand only
		teOctStructStats ConfParam,		// selected conformation
		int Limit,						// limit (0 if no limit) processing to this many bases total
		char *pszInFile,				// input MNase or conformational file
		char *pszInBioseqFile,			// bioseq genome file
		char *pszRsltsFile);			// output UCSC wiggle file

char *Conf2Txt(teOctStructStats ConfParam); // returns short description of specified conformational characteristic

CStopWatch gStopWatch;
CDiagnostics gDiagnostics;				// for writing diagnostics messages to log file
char gszProcName[_MAX_FNAME];			// process name

int m_hRsltsFile;					// handle for opened results file

CTwister *m_pTwister;				// used if processing structural conformations
int *m_pConfValues;					// to hold conformational values for current chromosome being processed
int m_AllocdConfValues;				// how many conformational values can be stored in m_pConfValues
CCSVFile *m_pMNaseCSV;				// used to hold DNase site preferences whilst loading into m_pMNaseSel 
CBioSeqFile *m_pBioSeqFile;			// holds instantiated genome assembly sequences, used if MNase scoring
etSeqBase *m_pChromSeq;				// holds current assembly chromosome sequence, used if MNase scoring
int m_AllocdChromSeq;				// allocated assembly chromosome sequence length
int m_ChromSeqLen;					// current assembly chromosome sequence length
char m_szCurChrom[cMaxDatasetSpeciesChrom]; // for this chromosome
int m_CurChromID;					// for this chromosome
double *m_pMNaseSel;				// allocated array of MNase site selection preferences (0.0..1.0) indexed by sequence octamers

tsHamHdr *m_pHamHdr;			// header for binary format hamming edit distances 
tsHamChrom *m_pCurHamChrom;		// pts to current chromosome specific binary hammings

CCSVFile *m_pCSVFile;			// used if processing input CSV file for coverage
CBEDfile *m_pBEDFile;			// used if processing input BED files for coverage

const int cMaxChromCov = 200;			// can handle at most this many chromosomes
const int cAllocCovCnts = 0x07fffff;	// allocate for chrom coverage cnts in this sized increments
typedef struct TAG_sChromCnts {
	char szChrom[cMaxDatasetSpeciesChrom+1];	// coverage is on this chromosome
	int AllocCovCnts;							// allocated to hold cnts for this sized chromosome
	int StartOfs;								// pCovCnts[offset] of first coverage cnt
	int EndOfs;									// pCovCnts[offset] of last coverage cnt
	UINT16 *pCovCnts;							// coverage counts
} tsChromCnts;

tsChromCnts m_ChromCnts[cMaxChromCov];
int m_NumChromsCov = 0;		

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
	return _T("genWiggle");
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

etPMode PMode;				// processing mode
int Rslt = 0;   			// function result code >= 0 represents success, < 0 on failure
int Strand;					// filter for this strand
teOctStructStats ConfParam;	// selected conformation

char szRsltsFile[_MAX_PATH];	// output wiggle file
char szInFile[_MAX_PATH];		// input MNase, conformation, or Hammings 
char szGenomeFile[_MAX_PATH];	// bioseq genome file

int Limit;						// 0 if no limit, otherwise only process for at most this number of bases

// command line args
struct arg_lit  *help    = arg_lit0("hH","help",                "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"diagnostics log file");
struct arg_int  *strand = arg_int0("s","strand","<int>",		"filter for this strand: 0 - any, 1 - Watson '+', 2 - Crick '-' (default is any)");
struct arg_int *pmode = arg_int0("m","mode","<int>",		    "processing mode: 0 - CSV, 1 - BED, 2 - MNase, 3 - Conformation, 4 - Hamming edit distances (default = 0)");
struct arg_int *confparam = arg_int0("c","conf","<int>",		"conformational characteristic: 0-energy,1-minorgroove,2-majorgroove,3-twist,4-roll,5-tilt,6-rise,7-slide,8-shift,9-rmsd,10-ORChid (default 1)");
struct arg_file *infile = arg_file1("i","in","<file>",			"input from this CSV loci, BED, MNase, conformational, or Hamming file");
struct arg_file *gfile = arg_file0("I","sfx","<file>",			"input bioseq genome file");
struct arg_file *outfile = arg_file1("o","out","<file>",		"output UCSC wiggle format to this file");
struct arg_int *limit = arg_int0("l","limit","<int>",		    "limit number of bases or features processed whilst debugging (0 == no limit)");


struct arg_end *end = arg_end(20);

void *argtable[] = {help,version,FileLogLevel,LogFile,
					pmode,strand,confparam,infile,gfile,outfile,limit,
					end};

char **pAllArgs;
int argerrors;
argerrors = CUtility::arg_parsefromfile(argc,(char **)argv,&pAllArgs);
if(argerrors >= 0)
	argerrors = arg_parse(argerrors,pAllArgs,argtable);

    /* special case: '--help' takes precedence over error reporting */
if (help->count > 0)
        {
		printf("\n%s a component of the K-mer Adaptive Next Generation Aligner toolset, Version %s\nOptions ---\n", gszProcName,cpszProgVer);
        arg_print_syntax(stdout,argtable,"\n");
        arg_print_glossary(stdout,argtable,"  %-25s %s\n");
		printf("\nNote: Parameters can be entered into a parameter file, one parameter per line.");
		printf("\n      To invoke this parameter file then precede it's name with '@'");
		printf("\n      e.g. %s @myparams.txt\n",gszProcName);
		printf("\nPlease report any issues regarding usage of %s at https://github.com/csiro-crop-informatics/biokanga/issues\n\n",gszProcName);
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

	if(PMode == ePMconf)
		{
		ConfParam = (teOctStructStats)(confparam->count ? confparam->ival[0] : eSSminorgroove);
		if(ConfParam < eSSenergy ||ConfParam >= eSSNumStatParams)
			{
			printf("\nError: conformation parameter '-c%d' specified outside of range %d..%d",ConfParam,eSSenergy,(int)eSSNumStatParams-1);
			exit(1);
			}
		}
	else
		ConfParam = eSSenergy;

	Strand = strand->count ? strand->ival[0] : 0;
	if(Strand < 0 || Strand > 2)
		{
		printf("\nError: Strand '-s%d' specified outside of range 0..2",Strand);
		exit(1);
		}

	if(Strand != 0 && PMode >= ePMNase)
		{
		printf("\nError: Strand '-s%d' only processed in -m0 (CSV) and -m1 (BED) modes",Strand);
		exit(1);
		}
	switch(Strand) {
		case 1: Strand = (int)'+'; break;
		case 2: Strand = (int)'-'; break;
		case 0: Strand = (int)'*'; break;
		}


	Limit = limit->count ? limit->ival[0] : 0;
	if(Limit < 0)
		Limit = 0;

	strcpy(szInFile,infile->filename[0]);
	strcpy(szRsltsFile,outfile->filename[0]);
	if(PMode == ePMNase || PMode == ePMconf)
		{
		if(gfile->count)
			strcpy(szGenomeFile,gfile->filename[0]);
		else
			{
			printf("\nError: no bioseq genome assembly file specified");
			exit(1);
			}
		}
	else
		szGenomeFile[0] = '\0';

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
			pszDescr = "CSV loci";
			break;
		case ePMBed:			// UCSC BED format
			pszDescr = "UCSC BED format";
			break;
		case ePMNase:				
			pszDescr = "MNase prefered sites";
			break;
		case ePMconf:				
			pszDescr = "Conformational characteristics";
			break;
		case ePMhamming:				
			pszDescr = "Hamming edit distances";
			break;

		}
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"processing mode is : '%s'",pszDescr);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"process for this strand only: '%c'",(char)Strand);
	if(PMode == ePMconf)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"conformational characteristic: '%s'",Conf2Txt(ConfParam));
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"input %s file: '%s'",pszDescr, szInFile);
	if(PMode == ePMNase || PMode == ePMconf)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"input bioseq genome file: '%s'",szGenomeFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"output UCSC Wiggle format file: '%s'",szRsltsFile);
	if(Limit > 0)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"limit processing to first %d bases or elements",Limit);

#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	gStopWatch.Start();
	// processing here...
	Rslt = Process(PMode,eMEGRAny,Strand,ConfParam,Limit,szInFile,szGenomeFile,szRsltsFile);

	gStopWatch.Stop();
	Rslt = Rslt >=0 ? 0 : 1;
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Exit code: %d Total processing time: %s",Rslt,gStopWatch.Read());
	exit(Rslt);
	}
else
	{
	printf("\n%s component of the the K-mer Adaptive Next Generation Aligner toolkit, Version %s\n",gszProcName,cpszProgVer);
	arg_print_errors(stdout,end,gszProcName);
	arg_print_syntax(stdout,argtable,"\nUse '-h' to view option and parameter usage\n");
	exit(1);
	}
}

char *
Conf2Txt(teOctStructStats ConfParam) // returns short description of specified conformational characteristic
{
char *pszDescr;
switch(ConfParam) {
	case eSSenergy:
		pszDescr = (char *)"minimal energy";
		break;
	case eSSminorgroove:				
		pszDescr = (char *)" minor groove";
		break;
	case eSSmajorgroove:				
		pszDescr = (char *)"major groove inferenced from octamer twist + rise";
		break;
	case eSStwist:					
		pszDescr = (char *)"twist";
		break;
	case eSSroll:					
		pszDescr = (char *)"roll";
		break;
	case eSStilt:					
		pszDescr = (char *)"tilt";
		break;
	case eSSrise:					
		pszDescr = (char *)"rise";
		break;
	case eSSslide:						
		pszDescr = (char *)"slide";
		break;
	case eSSshift:					
		pszDescr = (char *)"shift";
		break;
	case eSSrmsd:					
		pszDescr = (char *)"rmsd";
		break;
	case eSSORChidVal:				
		pszDescr = (char *)"hydroxyl radical cleavage value from ORChid dataset";
		break;
	default:
		pszDescr = (char *)"Unsupported";
		break;
	}
return(pszDescr);
}

void
Init(void)
{
m_pBioSeqFile = NULL;		// holds instantiated genome assembly sequences, used if MNase scoring
m_pChromSeq = NULL;			// holds current assembly chromosome sequence, used if MNase scoring
m_pMNaseSel = NULL;			// allocated array of MNase site selection preferences (0.0..1.0) indexed by sequence octamers
m_pMNaseCSV = NULL;			// used to hold DNase site preferences whilst loading into m_pMNaseSel 
m_pTwister = NULL;			// used if processing structural conformation
m_pConfValues = NULL;		// to hold conformational values for current chromosome being processed
m_pHamHdr = NULL;			// Hamming edit distances
m_pCurHamChrom = NULL;		// current Hamming chromosome

m_pCSVFile = NULL;
m_pBEDFile = NULL;

m_AllocdConfValues = 0;		// how many conformational values can be stored in m_pConfValues
m_hRsltsFile = -1;			// handle for opened results file
m_AllocdChromSeq = 0;		// allocated assembly chromosome sequence length
m_ChromSeqLen = 0;			// current assembly chromosome sequence length
m_CurChromID = 0;			// for this chromosome
m_szCurChrom[0] = '\0';

for(int ChromIdx = 0; ChromIdx < cMaxChromCov; ChromIdx++)
	{
	m_ChromCnts[ChromIdx].AllocCovCnts = 0;
	m_ChromCnts[ChromIdx].StartOfs = 0;
	m_ChromCnts[ChromIdx].EndOfs = 0;
	m_ChromCnts[ChromIdx].pCovCnts = NULL;
	m_ChromCnts[ChromIdx].szChrom[0] = '\0';
	}
m_NumChromsCov = 0;	
}

void
Reset(void)
{
if(m_hRsltsFile != -1)
	{
	close(m_hRsltsFile);
	m_hRsltsFile = -1;
	}
if(m_pBioSeqFile != NULL)
	{
	delete m_pBioSeqFile;
	m_pBioSeqFile = NULL;
	}
if(m_pTwister != NULL)
	{
	delete m_pTwister;
	m_pTwister = NULL;
	}
if(m_pConfValues != NULL)
	{
	delete m_pConfValues;
	m_pConfValues = NULL;
	}
if(m_pChromSeq != NULL)
	{
	delete m_pChromSeq;
	m_pChromSeq = NULL;
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
if(m_pHamHdr != NULL)
	{
	delete(m_pHamHdr);			
	m_pHamHdr = NULL;
	}
if(m_pBEDFile != NULL)
	{
	delete m_pBEDFile;
	m_pBEDFile = NULL;
	}
if(m_pCSVFile != NULL)
	{
	delete m_pCSVFile;
	m_pCSVFile = NULL;
	}
m_pCurHamChrom = NULL;
m_AllocdChromSeq = 0;				// allocated assembly chromosome sequence length
m_ChromSeqLen = 0;					// current assembly chromosome sequence length
m_CurChromID = 0;						// for this chromosome
m_szCurChrom[0] = '\0';
if(m_NumChromsCov > 0)
	{
	for(int ChromIdx = 0; ChromIdx < m_NumChromsCov; ChromIdx++)
		{
		m_ChromCnts[ChromIdx].StartOfs = 0;
		m_ChromCnts[ChromIdx].EndOfs = 0;
		if(m_ChromCnts[ChromIdx].pCovCnts != NULL)
			{
#ifdef _WIN32
			free(m_ChromCnts[ChromIdx].pCovCnts);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
			if(m_ChromCnts[ChromIdx].pCovCnts != MAP_FAILED)
				munmap(m_ChromCnts[ChromIdx].pCovCnts,m_ChromCnts[ChromIdx].AllocCovCnts * sizeof(UINT16));
#endif
			m_ChromCnts[ChromIdx].pCovCnts = NULL;
			}
		m_ChromCnts[ChromIdx].AllocCovCnts = 0;
		m_ChromCnts[ChromIdx].szChrom[0] = '\0';
		}
	m_NumChromsCov = 0;	
	}
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


teBSFrsltCodes
LoadHammings(char *pszHammings)
{
int hHamFile;
tsHamHdr HamHdr;

#ifdef _WIN32
hHamFile = open(pszHammings, O_READSEQ ); // file access will be sequential..
#else
hHamFile = open64(pszHammings, O_READSEQ ); // file access will be sequential..
#endif
if(hHamFile == -1)					// check if file open succeeded
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open %s - %s",pszHammings,strerror(errno));
	Reset();
	return(eBSFerrOpnFile);
	}

if(read(hHamFile,&HamHdr,sizeof(tsHamHdr))!=sizeof(tsHamHdr))
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to read %s - %s",pszHammings,strerror(errno));
	close(hHamFile);
	Reset();
	return(eBSFerrParse);
	}
if(HamHdr.Magic[0] != 'b' && HamHdr.Magic[0] != 'h' && HamHdr.Magic[0] != 'a' && HamHdr.Magic[0] != 'm')
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Not a binary format Hamming distance file - %s",pszHammings);
	close(hHamFile);
	Reset();
	return(eBSFerrFileType);
	}
if(HamHdr.NumChroms < 1)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Hamming distance file contains no edit distances - %s",pszHammings);
	close(hHamFile);
	Reset();
	return(eBSFerrNoEntries);
	}

if((m_pHamHdr = (tsHamHdr *)new UINT8 [HamHdr.Len])==NULL)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Unable to allocate memory (%d bytes) for holding Hamming distances loaded from - %s",HamHdr.Len,pszHammings);
	close(hHamFile);
	Reset();
	return(eBSFerrMem);
	}
memcpy(m_pHamHdr,&HamHdr,sizeof(tsHamHdr));
if(read(hHamFile,(UINT8 *)m_pHamHdr+sizeof(tsHamHdr),HamHdr.Len-sizeof(tsHamHdr))!=HamHdr.Len-sizeof(tsHamHdr))
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Unable to read all Hamming edit distances from - %s",pszHammings);
	close(hHamFile);
	Reset();
	return(eBSFerrFileAccess);
	}
close(hHamFile);
int ChromIdx;
for(ChromIdx = 0; ChromIdx < m_pHamHdr->NumChroms; ChromIdx++)
	{
	m_pCurHamChrom = (tsHamChrom *)((UINT8 *)m_pHamHdr + m_pHamHdr->ChromOfs[ChromIdx]);
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"chrom '%s' loaded with %d Hammings",m_pCurHamChrom->szChrom,m_pCurHamChrom->NumEls);
	}
m_pCurHamChrom = NULL;
return(eBSFSuccess);
}

int						// returned Hamming distance, 0 if no hammings loaded, -1 if not located
LocateHamming(char *pszChrom,UINT32 Loci)
{
int ChromIdx;
if(m_pHamHdr == NULL)
	return(0);
if(m_pCurHamChrom == NULL || stricmp((char *)m_pCurHamChrom->szChrom,pszChrom))
	{
	for(ChromIdx = 0; ChromIdx < m_pHamHdr->NumChroms; ChromIdx++)
		{
		m_pCurHamChrom = (tsHamChrom *)((UINT8 *)m_pHamHdr + m_pHamHdr->ChromOfs[ChromIdx]);
		if(!stricmp((char *)m_pCurHamChrom->szChrom,pszChrom))
			break;
		}
	if(ChromIdx == m_pHamHdr->NumChroms)
		{
		m_pCurHamChrom = NULL;
		return(-1);
		}
	}
if(m_pCurHamChrom->NumEls < Loci)
	return(-1);
return(m_pCurHamChrom->Dists[Loci]);
}

int
GenHammingWiggle(int Limit,				// limit (0 if no limit) processing to this many bases total
		char *pszHammFile,				// input binary format Hamming edit distance
		char *pszRsltsFile)				// output UCSC wiggle file
{
int Rslt;
int SeqIdx;
char szLineBuff[4096];
int BuffIdx;

gDiagnostics.DiagOutMsgOnly(eDLInfo,"Loading Hamming edit distances from file '%s'",pszHammFile);
if((Rslt = LoadHammings(pszHammFile)) < eBSFSuccess)
	{
	Reset();
	return(Rslt);
	}
gDiagnostics.DiagOutMsgOnly(eDLInfo,"Loading Hamming edit distances completed");

BuffIdx = sprintf(szLineBuff,"track type=wiggle_0 name=\"Hamming\" description=\"Hamming Edit Distances\"\n");
CUtility::SafeWrite(m_hRsltsFile,szLineBuff,BuffIdx);
BuffIdx = 0;

int ChromIdx;
int HammingDist;
for(ChromIdx = 0; ChromIdx < m_pHamHdr->NumChroms; ChromIdx++)
	{
	m_pCurHamChrom = (tsHamChrom *)((UINT8 *)m_pHamHdr + m_pHamHdr->ChromOfs[ChromIdx]);
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing chromosome: '%s'",m_pCurHamChrom->szChrom);
	BuffIdx += sprintf(&szLineBuff[BuffIdx],"fixedStep chrom=%s start=%d step=1\n",m_pCurHamChrom->szChrom,1);
	CUtility::SafeWrite(m_hRsltsFile,szLineBuff,BuffIdx);
	BuffIdx = 0;
	for(SeqIdx = 0; SeqIdx < (int)m_pCurHamChrom->NumEls; SeqIdx++)
		{
		if(Limit > 0 && SeqIdx > Limit)
			break;

		HammingDist = LocateHamming((char *)m_pCurHamChrom->szChrom,SeqIdx);
		BuffIdx += sprintf(&szLineBuff[BuffIdx],"%d\n",HammingDist);
		if(BuffIdx + 100 > sizeof(szLineBuff))
			{
			CUtility::SafeWrite(m_hRsltsFile,szLineBuff,BuffIdx);
			BuffIdx = 0;
			}
		}
	}
if(BuffIdx)
	CUtility::SafeWrite(m_hRsltsFile,szLineBuff,BuffIdx);
return(eBSFSuccess);
}


int
LoadNxtChrom(void)
{
int Rslt;

if(m_pBioSeqFile == NULL)
	return(0);
if((m_CurChromID = m_pBioSeqFile->Next(m_CurChromID)) < 1)
	{
	m_CurChromID = 0;						
	m_szCurChrom[0] = '\0';
	return(eBSFerrEntry);
	}

m_pBioSeqFile->GetName(m_CurChromID,sizeof(m_szCurChrom),m_szCurChrom);
m_ChromSeqLen = m_pBioSeqFile->GetDataLen(m_CurChromID);
if(m_pChromSeq == NULL || m_ChromSeqLen >= m_AllocdChromSeq)
	{
	if(m_pChromSeq != NULL)
		{
		delete m_pChromSeq;
		m_pChromSeq = NULL;
		m_AllocdChromSeq = 0;
		}
	if((m_pChromSeq = new etSeqBase [m_ChromSeqLen + cChromSeqReAlloc]) == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate memory for chromosome sequence of %d bytes",m_ChromSeqLen + cChromSeqReAlloc);
		return(eBSFerrMem);
		}
	m_AllocdChromSeq = m_ChromSeqLen + cChromSeqReAlloc;
	}
if((Rslt=m_pBioSeqFile->GetData(m_CurChromID,eSeqBaseType,0,m_pChromSeq,m_ChromSeqLen))!=m_ChromSeqLen)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to load chromosome sequence for '%s' of expected length",m_szCurChrom,m_ChromSeqLen);
	return(Rslt);
	}
return(m_ChromSeqLen);
}

int
InitMNaseSitePrefs(char *pszInMNaseFile)			// score from this MNase site selectivity file (generated by MNaseSitePred process)
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
	OctIdx = GenSeqIdx(8,Octamer);
	m_pMNaseCSV->GetDouble(4,&SitePref);
	m_pMNaseSel[OctIdx] = SitePref;
	}
delete m_pMNaseCSV;
m_pMNaseCSV = NULL;

return(Rslt);
}

int
GenMNaseWiggle(int Limit,				// limit (0 if no limit) processing to this many bases total
		char *pszInFile,				// input MNase
		char *pszInBioseqFile,			// bioseq genome file
		char *pszRsltsFile)				// output UCSC wiggle file
{
int Rslt;
etSeqBase *pSeq;
int SeqIdx;
int OctIdx;
double MNaseScore;
char szLineBuff[4096];
int BuffIdx;
bool bStartSect;

if((Rslt = InitMNaseSitePrefs(pszInFile)) < 0)
	{
	Reset();
	return(Rslt);
	}

BuffIdx = sprintf(szLineBuff,"track type=wiggle_0 name=\"MNase Prefsites\" description=\"MNase Cut Site Preferences\"\n");
CUtility::SafeWrite(m_hRsltsFile,szLineBuff,BuffIdx);
BuffIdx = 0;

while((Rslt = LoadNxtChrom()) > 0)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing chromosome: '%s'",m_szCurChrom);
	// have DNase and chromosome sequence, can now generate wiggle!
	bStartSect = true;
	pSeq = m_pChromSeq;
	for(SeqIdx = 0; SeqIdx+7 < m_ChromSeqLen; SeqIdx++,pSeq++)
		{
		OctIdx = GenSeqIdx(8,pSeq);
		if(OctIdx < 0)
			{
			bStartSect = true;
			continue;
			}
		if(bStartSect)
			{
			BuffIdx += sprintf(&szLineBuff[BuffIdx],"fixedStep chrom=%s start=%d step=1\n",m_szCurChrom,SeqIdx+1);
			CUtility::SafeWrite(m_hRsltsFile,szLineBuff,BuffIdx);
			BuffIdx = 0;
			bStartSect = false;
			}

		MNaseScore = m_pMNaseSel[OctIdx];
		BuffIdx += sprintf(&szLineBuff[BuffIdx],"%1.4f\n",MNaseScore);
		if(BuffIdx + 100 > sizeof(szLineBuff))
			{
			CUtility::SafeWrite(m_hRsltsFile,szLineBuff,BuffIdx);
			BuffIdx = 0;
			}
		if(Limit > 0 && SeqIdx > Limit)
			break;
		}
	}
if(BuffIdx)
	CUtility::SafeWrite(m_hRsltsFile,szLineBuff,BuffIdx);
return(eBSFSuccess);
}

int
initConformation(char *pszConfFile)		// input structural conformation characteristics file
{
int Rslt;
if((m_pTwister = new CTwister)==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,"genWiggle","Unable to create CTwister object");
	Reset();
	return(eBSFerrObj);
	}

if((Rslt = m_pTwister->LoadStructParams(pszConfFile))  < eBSFSuccess)
	{
	while(m_pTwister->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,m_pTwister->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,"ProcessBioseqStruct","LoadStructParams(%s) failed",pszConfFile);
	Reset();
	return(Rslt);
	}
return(eBSFSuccess);
}


int
BuildReadCoverage(char *pszChrom,		// coverage is onto this chrom
			  int StartOfs,				// coverage start at this offset 
			  int EndOfs,				// and ends at this offset inclusive
			  int Cnt)					// increment coverage by this
{
tsChromCnts *pChrom;
int ChromIdx;
int AllocCovCnts;
UINT16 *pCovCnts;
size_t ReallocTo;

if(pszChrom == NULL || pszChrom[0] == '\0')
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"BuildReadCoverage: No chromosome specified");
	Reset();
	return(eBSFerrChrom);
	}

// arbitary count clamping to be in range 1..10000 inclusive
if(Cnt < 1)
	Cnt = 1;
else
	if(Cnt > 10000)
		Cnt = 10000;

// ensure StartOfs and EndOfs are both >= 0
if(StartOfs < 0)
	StartOfs = 0;
if(EndOfs < 0)
	EndOfs = 0;

// ensure StartOfs <= EndOfs
if(StartOfs > EndOfs)
	{
	int TmpOfs = EndOfs;
	EndOfs = StartOfs;
	StartOfs = TmpOfs;
	}

// check if this is a new chrom or if coverage is onto an existing chrom
pChrom = &m_ChromCnts[0];
ChromIdx = 0;
if(m_NumChromsCov > 0)
	{
	for(ChromIdx = 0; ChromIdx < m_NumChromsCov; ChromIdx++,pChrom++)
		{
		if(!strnicmp(pszChrom,pChrom->szChrom,cMaxDatasetSpeciesChrom))
			break;
		}
	}
if(ChromIdx == m_NumChromsCov)	// if a new or first chrom
	{
	strncpy(pChrom->szChrom,pszChrom,cMaxDatasetSpeciesChrom);
	pChrom->szChrom[cMaxDatasetSpeciesChrom-1] = 0;
	pChrom->StartOfs = StartOfs;
	pChrom->EndOfs = EndOfs;
	AllocCovCnts = EndOfs + cAllocCovCnts;
	ReallocTo =  AllocCovCnts * sizeof(UINT16);
#ifdef _WIN32
	pChrom->pCovCnts = (UINT16 *) malloc(ReallocTo);	// initial and perhaps the only allocation
	if(pChrom->pCovCnts == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"BuildReadCoverage: Memory allocation of %lld bytes - %s",(INT64)ReallocTo,strerror(errno));
		Reset();
		return(eBSFerrMem);
		}
#else
	// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
	pChrom->pCovCnts = (UINT16 *)mmap(NULL,ReallocTo, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
	if(pChrom->pCovCnts == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"BuildReadCoverage: Memory allocation of %lld bytes through mmap()  failed - %s",(INT64)ReallocTo,strerror(errno));
		pChrom->pCovCnts = NULL;
		Reset();
		return(eBSFerrMem);
		}
#endif
	pChrom->AllocCovCnts = AllocCovCnts;
	memset(pChrom->pCovCnts,0,ReallocTo);
	m_NumChromsCov += 1;
	}

// check if chrom coverage cnts needs to be extended
if(EndOfs >= pChrom->AllocCovCnts)
	{
	AllocCovCnts = EndOfs + cAllocCovCnts;
	ReallocTo = AllocCovCnts * sizeof(UINT16);
#ifdef _WIN32
	pCovCnts = (UINT16 *) realloc(pChrom->pCovCnts,ReallocTo);
#else
	pCovCnts = (UINT16 *)mremap(pChrom->pCovCnts,pChrom->AllocCovCnts * sizeof(UINT16),ReallocTo,MREMAP_MAYMOVE);
	if(pCovCnts == MAP_FAILED)
		pCovCnts = NULL;
#endif
	if(pCovCnts == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"BuildReadCoverage: Memory re-allocation to %d bytes - %s",ReallocTo,strerror(errno));
		Reset();
		return(eBSFerrMem);
		}
	pChrom->pCovCnts = pCovCnts;
	memset(&pChrom->pCovCnts[pChrom->AllocCovCnts],0,(AllocCovCnts - pChrom->AllocCovCnts) * sizeof(UINT16));
	pChrom->AllocCovCnts = AllocCovCnts;
	}

if(EndOfs > pChrom->EndOfs)
	pChrom->EndOfs = EndOfs;
if(StartOfs < pChrom->StartOfs)
	pChrom->StartOfs = StartOfs;

pCovCnts = &pChrom->pCovCnts[StartOfs];
while(StartOfs++ <= EndOfs)
	{
	// clamp accumulated cnts to be no greater than 0x0fffe
	if((0x0fffe - *pCovCnts) > Cnt)
		*pCovCnts += Cnt;
	else
		*pCovCnts = 0x0fffe;
	pCovCnts += 1;
	}

return(eBSFSuccess);
}

int
WriteReadsWig(char *pszSrcFile)
{
int BuffIdx;
char szLineBuff[8096];
tsChromCnts *pChrom;
UINT16 *pCnts;
int ChromIdx;
int SeqIdx;
bool bStartRegion;
int EmptyRegionLen;

BuffIdx = sprintf(szLineBuff,"track type=wiggle_0 color=50,150,255 autoScale=off maxHeightPixels=128:32:8 name=\"Reads - %s\" description=\"Reads distribution for %s\"\n",pszSrcFile,pszSrcFile);
CUtility::SafeWrite(m_hRsltsFile,szLineBuff,BuffIdx);
BuffIdx = 0;
pChrom = &m_ChromCnts[0];
for(ChromIdx = 0 ; ChromIdx < m_NumChromsCov; ChromIdx++,pChrom++)
	{
	EmptyRegionLen = 0;
	bStartRegion = true;
	pCnts = &pChrom->pCovCnts[pChrom->StartOfs];
	for(SeqIdx = pChrom->StartOfs; SeqIdx < pChrom->EndOfs; SeqIdx++,pCnts++)
		{
		if(*pCnts == 0)
			{
			EmptyRegionLen += 1;
			if(EmptyRegionLen > 10)
				bStartRegion = true;
			continue;
			}
		else
			{
			if(!bStartRegion && EmptyRegionLen > 0 && EmptyRegionLen <= 10)
				{
				while(EmptyRegionLen--)
					BuffIdx += sprintf(&szLineBuff[BuffIdx],"0\n");
				}
			EmptyRegionLen = 0;
			}
		if(bStartRegion)
			{
			BuffIdx += sprintf(&szLineBuff[BuffIdx],"fixedStep chrom=%s start=%d step=1\n",pChrom->szChrom,SeqIdx+1);
			CUtility::SafeWrite(m_hRsltsFile,szLineBuff,BuffIdx);
			BuffIdx = 0;
			bStartRegion = false;
			}
		BuffIdx += sprintf(&szLineBuff[BuffIdx],"%d\n",*pCnts);
		if(BuffIdx + 100 > sizeof(szLineBuff))
			{
			CUtility::SafeWrite(m_hRsltsFile,szLineBuff,BuffIdx);
			BuffIdx = 0;
			}

		}
	}
if(BuffIdx)
	CUtility::SafeWrite(m_hRsltsFile,szLineBuff,BuffIdx);
return(eBSFSuccess);
}

int GenCSVWiggle(char FiltStrand,					// process for this strand only
				 int Limit,						// limit (0 if no limit) processing to this many elements total
				 char *pszInFile,				// CSV loci input file
				 char *pszRsltsFile)			// output to this UCSC wiggle file
{
int Rslt;
int NumEls;
int NumFields;
char *pszStrand;
char *pszChrom;
int StartLoci;
int EndLoci;

if((m_pCSVFile = new CCSVFile)==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CCSVFile");
	Reset();
	return(eBSFerrObj);
	}

if((Rslt=m_pCSVFile->Open(pszInFile)) !=eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open '%s' for processing",pszInFile);
	Reset();
	return(eBSFerrOpnFile);
	}
NumEls = 0;
while((Rslt=m_pCSVFile->NextLine()) > 0)	// onto next line containing fields
	{
	NumFields = m_pCSVFile->GetCurFields();
	if(NumFields < 7)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expected at least 7 fields in '%s', GetCurFields() returned '%d'",pszInFile,NumFields);
		return(eBSFerrFieldCnt);
		}

	if(FiltStrand != '*')
		{
		if(NumFields >= 8)					// check if strand has been specified
			{
			m_pCSVFile->GetText(8,&pszStrand);
			if(pszStrand[0] != '-')			// assume anything other than '-' is on the plus strand
				pszStrand[0] = '+';
			if(FiltStrand != pszStrand[0])
				continue;
			}
		else
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Strand filtering requested, expected at least 8 fields in '%s', GetCurFields() returned '%d'",pszInFile,NumFields);
			return(eBSFerrFieldCnt);
			}
		}
	NumEls += 1;
	m_pCSVFile->GetText(4,&pszChrom);
	m_pCSVFile->GetInt(5,&StartLoci);
	m_pCSVFile->GetInt(6,&EndLoci);

	if((Rslt=BuildReadCoverage(pszChrom,StartLoci,EndLoci,1))!=eBSFSuccess)
		break;
	if(Limit && NumEls > Limit)
		break;
	}
if(Rslt == eBSFSuccess)
	return(WriteReadsWig(pszInFile));
return(Rslt);
}

int GenBEDWiggle(char FiltStrand,					// process for this strand only
				 int Limit,						// limit (0 if no limit) processing to this many bases total
				 etBEDRegion Region,			// which regions are of interest
				 char *pszInFile,				// UCSC BED input file
				 char *pszRsltsFile)			// output to this UCSC wiggle file
{
int Rslt;
int CurFeatureID;
int NumExons;
int NumIntrons;
int CDSstart;
int CDSend;
int StartLoci;
int EndLoci;
int Score;
char szChrom[128];
char szPrevChrom[128];
char szFeatName[128];
char Strand;
int RefID;
int IntergenicStart;
int Idx;
int NumEls;

if((m_pBEDFile = new CBEDfile)==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CBEDfile");
	Reset();
	return(eBSFerrObj);
	}

if((Rslt=m_pBEDFile->Open(pszInFile,eBTAnyBed)) !=eBSFSuccess)
	{
	while(m_pBEDFile->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,m_pBEDFile->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open '%s' for processing",pszInFile);
	Reset();
	return(eBSFerrOpnFile);
	}

if(!m_pBEDFile->ContainsGeneDetail() && Region > eMEGRIntergenic)			// returns true if file contains gene detail (utr/cds/intron etc)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"bed file %s does not contain gene regions",pszInFile);
	Reset();
	return(eBSFerrFeature);
	}

// now iterate over the features, filtering as may be appropriate
szPrevChrom[0] = '\0';
CurFeatureID = 0;
IntergenicStart = 0;
RefID = 0;
Rslt = eBSFSuccess;
NumEls = 0;
while(Rslt == eBSFSuccess && (CurFeatureID = m_pBEDFile->GetNextFeatureID(CurFeatureID)) > 0)
	{

	m_pBEDFile->GetFeature(CurFeatureID,	// feature instance identifier
				szFeatName,				// where to return feature name
				szChrom,				// where to return chromosome name
				&StartLoci,				// where to return feature start on chromosome (0..n) 
				&EndLoci,				// where to return feature end on chromosome
 				&Score,					// where to return score
 				&Strand);				// where to return strand

	if(CurFeatureID == 1 || stricmp(szChrom,szPrevChrom))	// if new chromosome then reset IntergenicStart
		{
		strcpy(szPrevChrom,szChrom);
		IntergenicStart = 0;
		}

	if(FiltStrand != '*' && FiltStrand != Strand)
		continue;

	NumEls += 1;
	if(Limit && NumEls > Limit)
		break;

	if(Region != eMEGRAny)
		{
		NumExons = m_pBEDFile->GetNumExons(CurFeatureID);		// returns number of exons - includes UTRs + CDS
		NumIntrons = m_pBEDFile->GetNumIntrons(CurFeatureID);
		CDSstart = StartLoci + m_pBEDFile->GetCDSStart(CurFeatureID);		// returns relative start offset of CDS - NOTE add to '+' strand gene start, subtract on '-' strand gene start
		CDSend = StartLoci + m_pBEDFile->GetCDSEnd(CurFeatureID);			// returns relative end offset of CDS - NOTE add to '+' strand gene start, subtract on '-' strand gene start
		}
	Rslt = eBSFSuccess;

	switch(Region) {
		case eMEGRAny:					// process any region
			Rslt=BuildReadCoverage(szChrom,StartLoci,EndLoci,1);
			continue;

		case eMEGRIntergenic:	// only process intergenic
			if(IntergenicStart < StartLoci)
				Rslt=BuildReadCoverage(szChrom,IntergenicStart,StartLoci-1,1);
			if(IntergenicStart <= EndLoci)
				IntergenicStart = EndLoci+1;
			continue;

		case eMEGRExons:		// only process exons
			for(Idx = 1; Idx <= NumExons; Idx++)
				{
				StartLoci = m_pBEDFile->GetExonStart(CurFeatureID,Idx);
				EndLoci   = m_pBEDFile->GetExonEnd(CurFeatureID,Idx);
				if(StartLoci <= EndLoci)
					Rslt=BuildReadCoverage(szChrom,StartLoci,EndLoci,1);
				}
			continue;

		case eMEGRIntrons:		// only process introns
			for(Idx = 1; Idx <= NumIntrons; Idx++)
				{
				StartLoci = m_pBEDFile->GetIntronStart(CurFeatureID,Idx);
				EndLoci = m_pBEDFile->GetIntronEnd(CurFeatureID,Idx);
				if(StartLoci <= EndLoci)
					Rslt=BuildReadCoverage(szChrom,StartLoci,EndLoci,1);
				}
			continue;

		case eMEGRCDS:			// only process CDSs
			for(Idx = 1; Idx <= NumExons; Idx++)
				{
				StartLoci = m_pBEDFile->GetExonStart(CurFeatureID,Idx);
				EndLoci   = m_pBEDFile->GetExonEnd(CurFeatureID,Idx);
				if(EndLoci < CDSstart || StartLoci > CDSend)
					continue;
				if(StartLoci < CDSstart)
					StartLoci = CDSstart;
				if(EndLoci > CDSend)
					EndLoci = CDSend;
				if(StartLoci <= EndLoci)
					BuildReadCoverage(szChrom,StartLoci,EndLoci,1);
				}
			continue;

		case eMEGUTR:			// only process UTRs - single exon may have both 5' and 3' UTRs
			for(Idx = 1; Idx <= NumExons; Idx++)
				{
				StartLoci = m_pBEDFile->GetExonStart(CurFeatureID,Idx);
				EndLoci   = m_pBEDFile->GetExonEnd(CurFeatureID,Idx);
				if(EndLoci <= CDSend && StartLoci >= CDSstart) // is exon CDS only?
					continue;

				// check if 5' UTR 
				if(StartLoci < CDSstart)
					{
					if(EndLoci >= CDSstart)
						Rslt=BuildReadCoverage(szChrom,StartLoci,CDSstart-1,1);
					else
						Rslt=BuildReadCoverage(szChrom,StartLoci,EndLoci,1);
					}
					
				// check if 3'UTR
				if(EndLoci > CDSend)
					{
					if(StartLoci <= CDSend)
						StartLoci = CDSend+1;
					Rslt=BuildReadCoverage(szChrom,StartLoci,EndLoci,1);
					}
				}
			continue;

		case eMEG5UTR:			// only process 5'UTRs - strand sensitive
			for(Idx = 1; Idx <= NumExons; Idx++)
				{
				StartLoci = m_pBEDFile->GetExonStart(CurFeatureID,Idx);
				EndLoci   = m_pBEDFile->GetExonEnd(CurFeatureID,Idx);
				if(EndLoci <= CDSend && StartLoci >= CDSstart) // is exon CDS only?
					continue;

				if(Strand != '-')
					{
					// check if 5' UTR on '+' strand 
					if(StartLoci < CDSstart)
						{
						if(EndLoci >= CDSstart)
							EndLoci = CDSstart - 1;
						}
					}
				else 
					{	
					// check if 5'UTR on '-' strand
					if(EndLoci > CDSend)
						{
						if(StartLoci <= CDSend)
							StartLoci = CDSend+1;
						}
					}
				Rslt=BuildReadCoverage(szChrom,StartLoci,EndLoci,1);
				}
			continue;

		case eMEG3UTR:			// only process 3'UTRs  - strand sensitive
			for(Idx = 1; Idx <= NumExons; Idx++)
				{
				StartLoci = m_pBEDFile->GetExonStart(CurFeatureID,Idx);
				EndLoci   = m_pBEDFile->GetExonEnd(CurFeatureID,Idx);
				if(EndLoci <= CDSend && StartLoci >= CDSstart) // is exon CDS only?
					continue;

				if(Strand == '-')
					{
					// check if 3' UTR on '-' strand 
					if(StartLoci < CDSstart)
						{
						if(EndLoci >= CDSstart)
							EndLoci = CDSstart - 1;
						}
					}
				else 
					{	
					// check if 3'UTR on '+' strand
					if(EndLoci > CDSend)
						{
						if(StartLoci <= CDSend)
							StartLoci = CDSend+1;
						}
					}
				Rslt=BuildReadCoverage(szChrom,StartLoci,EndLoci,1);
				}
			continue;
		}
	}
if(Rslt == eBSFSuccess)
	return(WriteReadsWig(pszInFile));
return(Rslt);
}

int
GenConfWiggle(int Limit,						// limit (0 if no limit) processing to this many bases total
		teOctStructStats ConfParam,	// selected conformation
		char *pszInFile,				// input MNase
		char *pszInBioseqFile,			// bioseq genome file
		char *pszRsltsFile)				// output UCSC wiggle file
{
int Rslt;
int BuffIdx;
char szLineBuff[4096];
bool bStartSect;

char *pszConf;
etSeqBase *pSeq;
int *pConfValue;
double ConfValue;
int ConfIdx;

if((Rslt = initConformation(pszInFile)) < 0)
	{
	Reset();
	return(Rslt);
	}
pszConf = Conf2Txt(ConfParam);
BuffIdx = sprintf(szLineBuff,"track type=wiggle_0 name=\"Conformation - %s\" description=\"DNA Conformation for %s\"\n",pszConf,pszConf);
CUtility::SafeWrite(m_hRsltsFile,szLineBuff,BuffIdx);
BuffIdx = 0;

while((Rslt = LoadNxtChrom()) > 0)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing chromosome: '%s'",m_szCurChrom);
	bStartSect = true;
	if(m_pConfValues == NULL || m_ChromSeqLen > m_AllocdConfValues)
		{
		if(m_pConfValues != NULL)
			{
			delete m_pConfValues;
			m_pConfValues = NULL;
			m_AllocdConfValues = 0;
			}
		if((m_pConfValues = new int [m_ChromSeqLen + cChromSeqReAlloc]) == NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate memory (%d bytes) for holding conformation values",m_ChromSeqLen + cChromSeqReAlloc);
			Rslt = eBSFerrMem;
			break;
			}
		m_AllocdConfValues = m_ChromSeqLen + cChromSeqReAlloc;
		}

	pSeq = m_pChromSeq;
	m_pChromSeq[m_ChromSeqLen] = eBaseEOS;
	if((Rslt = m_pTwister->GetSequenceConformation(ConfParam,	// process for this conformational parameter
							  0,						// initial starting offset (0..n) in pSeq
							  0,						// number of steps (0 for all) to process starting at pSeq[iStartPsn]|pSeq[iStartPsn+1]
  							  m_ChromSeqLen,			// total length of sequence
							  m_pChromSeq,				// sequence to be processed
							  m_pConfValues))!=eBSFSuccess) // where to return step conformational values
					{
					gDiagnostics.DiagOut(eDLFatal,"genWiggle","GetSequenceConformation failed");
					break;
					}
	pConfValue = m_pConfValues;
	for(ConfIdx = 0; ConfIdx < m_ChromSeqLen; ConfIdx++,pConfValue++)
		{
		if(*pConfValue != INT_MIN)
			ConfValue = (double)*pConfValue/10000.0f;
		else
			{
			bStartSect = true;
			continue;
			}
		if(bStartSect)
			{
			BuffIdx += sprintf(&szLineBuff[BuffIdx],"fixedStep chrom=%s start=%d step=1\n",m_szCurChrom,ConfIdx+1);
			CUtility::SafeWrite(m_hRsltsFile,szLineBuff,BuffIdx);
			BuffIdx = 0;
			bStartSect = false;
			}
		BuffIdx += sprintf(&szLineBuff[BuffIdx],"%1.3f\n",ConfValue);
		if(BuffIdx + 100 > sizeof(szLineBuff))
			{
			CUtility::SafeWrite(m_hRsltsFile,szLineBuff,BuffIdx);
			BuffIdx = 0;
			}
		if(Limit > 0 && ConfIdx > Limit)
			break;
		}
	}
if(BuffIdx)
	CUtility::SafeWrite(m_hRsltsFile,szLineBuff,BuffIdx);
return(eBSFSuccess);
}

int
Process(etPMode PMode,					// processing mode
		etBEDRegion Region,				// regions of interest if processing BED files
		char Strand,					// process for this strand only
		teOctStructStats ConfParam,		// selected conformation
		int Limit,						// limit (0 if no limit) processing to this many bases total
		char *pszInFile,				// input CSV, BED, MNase, conformational file, or Hamming distances
		char *pszInBioseqFile,			// bioseq genome file
		char *pszRsltsFile)				// output UCSC wiggle file
{
int Rslt;

Init();

if(PMode == ePMNase || PMode == ePMconf)
	{
	if((m_pBioSeqFile = new CBioSeqFile) == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CBioSeqFile object");
		Reset();
		return(eBSFerrObj);
		}

	if((Rslt = m_pBioSeqFile->Open(pszInBioseqFile))!=eBSFSuccess)
		{
		while(m_pBioSeqFile->NumErrMsgs())
			gDiagnostics.DiagOut(eDLFatal,gszProcName,m_pBioSeqFile->GetErrMsg());
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open genome assembly sequence file '%s'",pszInBioseqFile);
		Reset();
		return(Rslt);
		}
	}

#ifdef _WIN32
if((m_hRsltsFile = open(pszRsltsFile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
if((m_hRsltsFile = open(pszRsltsFile, O_RDWR | O_CREAT | O_TRUNC, S_IREAD | S_IWRITE ))==-1)
#endif
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create UCSC Wiggle file: %s - %s",pszRsltsFile,strerror(errno));
	Reset();
	return(eBSFerrCreateFile);
	}

switch(PMode) {
	case ePMdefault:		// default is for CSV loci processing
		Rslt = GenCSVWiggle(Strand,Limit,pszInFile,pszRsltsFile);
		break;

	case ePMBed:			// UCSC BED processing
		Rslt = GenBEDWiggle(Strand,Limit,Region,pszInFile,pszRsltsFile);
		break;

	case ePMNase:			// MNase (default) processing mode
		Rslt = GenMNaseWiggle(Limit,pszInFile,pszInBioseqFile,pszRsltsFile);
		break;
	case ePMconf:			// Conformation processing mode
		Rslt = GenConfWiggle(Limit,ConfParam,pszInFile,pszInBioseqFile,pszRsltsFile);
		break;
	case ePMhamming:			// Hamming distance processing mode
		Rslt = GenHammingWiggle(Limit,pszInFile,pszRsltsFile);
		break;
	default:
		break;
	}
Reset();
return(Rslt);
}



