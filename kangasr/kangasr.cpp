// kangasr.cpp : Defines the entry point for the console application.
//

// usimreads.cpp : Defines the entry point for the console application
// simulates reads by randomly sampling from user specified target genome
// will then mutate these reads according to a user specified distribution profile

// 1.7.7   ensured that if no errors then output files _commit/fsync'd before close
// Release 1.7.8   Allows comment lines in parameter files
// Release 1.10.0	Restructured build directories to reflect application name
// Release 1.10.1   Added capability to add adaptor/linker sequence artifacts to 5' or 3' of generated reads
//                  Default read length increased to 100bp, default output format now to be fasta
// Release 1.11.0   public binary release
// Release 1.11.1   Bug fix for user entered parameters not being correctly defaulted when generating hammings
// Release 1.11.2   Bug fix for memory allocation segfault when processing input multifasta files with millions of entries
// Release 1.11.3   Check and warn if chrom lengths in input assembly are longer than the simulated read length requested
// Release 1.11.4   When generating paired ends, default min and max fragment lengths set to be 2x and 3x respectively the read length
// Release 1.11.6   Allow paired end reads to overlap enabling simulated datasets which can be processed by AllpathsLG
// Release 1.12.0   public binary release
// Release 1.12.1   changed _MAX_PATH to 260 to maintain compatibility with windows
// Release 1.12.2   increased upper limit of 100M sinulated reads to be 500M
// Release 1.12.4   increased simulated SNPs maximum rate to be 20000 per million bases
// Release 1.12.4   Updated illumina profile
// Release 1.12.5   Updated illumina profile to include small 5' end bias

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

const char *cpszProgVer = "1.12.5";		// increment with each release

const int cMaxArtefSeqs = 20;			// allow at most this number of artefact sequences
const int cMaxArtefSeqLen = 40;			// artefact sequences can be at most this length
const char *pszArtef5Seq = "ACACTCTTTCCCTACACGACGCTGTTCCATCT";	// default artifact seq for 5' simulated read ends (Illumina Single End Adapter 1)
const char *pszArtef3Seq = "ACACTCTTTCCCTACACGACGCTCTTCCGATCT"; // default artefact seq for 3' simulated read ends (Illumina Single End Sequencing Primer)


const int cDfltNumReads = 10000000;		// default number of reads
const int cMinNumReads =  100;			// minimum number of reads
const int cMaxNumReads =  500000000;	// maximum number of reads
const int cDfltReadLen =  100;			// default read length
const int cMinReadLen  =  20;			// minimum read length
const int cMaxSimReadLen = 1000;		// max simulated read length

const int cUpdnstream = 2000;			// when regional filtering this is the length of the up/down stream region

const int cMinCutLen   =  cMinReadLen;	// minimum cut length
const int cMaxCutLen   =  1000;			// maximum cut length

const int cMinPEFragLen     = 75;		// paired end minimum fragment length
const int cMaxPEFragLen     = 100000;	// paired end maximum fragment length

const int cMaxDistClusterLen = 300;		// max clustered read bin size

const int cMaxHammingBatchSize = 5000;	// max number of dynamic Hamming reads to simulate per thread before checkpointing to disk
const int cMaxProfileBatchSize = 500000;  // max number of end profiled reads to simulate per thread before checkpointing to disk
const int cMaxBatchSize = 5000000;		// max number of defaulted Hamming and non-profiled reads to simulate per thread before checkpointing to disk

#define NUCONLYMSK (~cRptMskFlg & 0x0f)	// hiorder bits used as attributes - bit 5 used to flag subsequence already selected
#define SSSELECTED 0x010				// used as an attribute to flag subsequence starting this loci already selected


const int cMaxWorkerThreads = 64;		// allow for at most 64 worker threads

// processing modes
typedef enum TAG_ePMode {
	ePMStandard,					// default - standard random start and end
	ePMProfRand,					// profiled start with random end sites
	ePMRandProf,					// random start with profiled end sites
	ePMProfProf,					// profiled start with profiled end sites
	ePMSampHamm,					// same as ePMStandard except hammings also generated in same format as 'uhamming' generated
	ePMplaceholder					// used to set the enumeration range
	} etPMode;

typedef enum eBEDRegion {
	eMEGRAny = 0,		// process any region
	eMEGRCDS,			// part of feature overlaps CDS
	eMEGR5UTR,			// part of feature overlaps 5'UTR
	eMEGR3UTR,			// part of feature overlaps 3'UTR
	eMEGRIntrons,		// part of feature overlaps Intron
	eMEGRUpstream,		// part of feature overlaps 5'upstream regulatory
	eMEGRDnstream,		// part of feature overlaps 3'downstream regulatory
	eMEGRIntergenic		// part of feature overlaps intergenic
} etBEDRegion;


char *Region2Txt(etBEDRegion Region);

// output format
typedef enum TAG_eFMode {
	eFMcsv,					// default - CSV loci only
	eFMcsvSeq,				// CSV loci + sequence
	eFMFasta,				// multifasta, sequences wrapped if longer than 79bp
	eFMNWFasta,				// multifasta, sequences non-wrapped even if longer than 79bp
	eFMSOLiD,				// SOLiD colorspace reads as csfasta
	eFMSOLiDbwa,			// SOLiD colorspace reads in double encoded basespace to suit BWA as fastq
	eFMplaceholder			// used to set the enumeration range
	} etFMode;

// simulated error rate mode selection
typedef enum TAG_eSEMode {
	eSEPnone,				// no simulated errors
	eSEPfixerrs,			// simulate fixed number of errors in each read
	eSEPstatic,				// use internal static profile
	eSEPdyn,				// dynamic according to '-z<rate>'
	eSEPplaceholder			// used to set the enumeration range
	} etSEMode;

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

int	TrimSeqChrs(char *pszTxt);	// trims quote marks, space/tabs and validates sequence as only containing acgtu

int
Process(etPMode PMode,		// processing mode
		etSEMode SEMode,	// induced sequencer error rate mode
		bool bPEgen,		// true if paired ends are to be generated
		int PEmin,			// PE minimum insert
		int PEmax,			// PE maximum insert
		double PropRandReads, // generate completely random reads at this rate
		int DistCluster,	// cluster generated reads into windows of this median length, 0 if no clustering
		double SeqErrRate,	// dynamic sequencing errors are to be induced into generated sequences at this rate
		bool bSeqErrProfile,// true if to generate composite sequencer errors with uniform profile (default is Illumina 3' skewed)
		int SNPrate,		// simulate SNPs at this rate per million bases
		int InDelSize,		// simulated InDel size range
		double InDelRate,	// simulated InDel rate per read
		bool bReadHamDist,	// true if hamming distributions from each sampled read to all other genome subsequences to be generated
		etFMode FMode,		// output format
		int NumThreads,		// number of worker threads to use
		char Strand,		// generate for this strand '+' or '-' or for both '*'
		int NumReads,		// number of reads required (will be doubled if paired end reads)
		int ReadLen,		// read lengths
		double Artef5Rate,			// rate (0..1) at which to insert 5' artefact sequences
		int NumArtef5Seqs,			// number of user specified 5' artefact sequences
		char *pszArtef5Seqs[], // 5' artefact sequences
		double Artef3Rate,			// rate (0..1) at which to insert 3' artefact sequences
		int NumArtef3Seqs,			// number of user specified 3' artefact sequences
		char *pszArtef3Seqs[], // 5' artefact sequences
		int CutMin,			// min cut length
		int CutMax,			// max cut length
		bool bDedupe,		// true if unique read sequences only to be generated
		int DfltHamming,	// if < 0 then dynamically determine Hamming distance otherwise use this value as the Hamming
		int Region,			// Process regions 0:ALL,1:CDS,2:5'UTR,3:3'UTR,4:Introns,5:5'US,6:3'DS,7:Intergenic (default = ALL)
		int UpDnRegLen,		// if processing regions then up/down stream regulatory length
		char *pszFeatFile,	// optionally generate transcriptome reads from features or genes in this BED file
		char *pszInFile,	// input from this bioseq assembly
		char *pszProfFile,	// input from this profile site preferences file
		char *pszHammFile,	// use Hamming edit distances from this file
		char *pszOutPEFile, // output partner paired end simulated reads to this file
		char *pszOutFile,	// output simulated reads to this file
		char *pszOutSNPs);   // output simulated SNP loci to this file

int
LocateRevCplSeq(int Len,etSeqBase *pProbe,int NumSortedReads);

int			// index of exactly matching probe or -1 if no match
LocateFirstExact(etSeqBase *pProbe,				// pts to probe sequence
				 int ProbeLen,					// probe length to exactly match over
				  int IdxLo,					// low index in m_pSimReads
				  int IdxHi);					// high index in m_pSimReads

static int SortSimReads(const void *arg1, const void *arg2);
static int SortSimLoci(const void *arg1, const void *arg2);

int		// returned minimum Hamming distance, checks for self-loci
MinHammingDistW(int MinHamming,	// initial minimum Hamming distance
			   int ReadLen,		// read length
			   int ChromID,     // from which chromosome was this read was derived
			   int ReadLoci,	// read start loci
			   etSeqBase *pRead); // read sequence

int		// returned minimum Hamming distance, no check for self-loci
MinHammingDistC(int MinHamming,	// initial minimum Hamming distance
			   int ReadLen,		// read length
			   int ChromID,     // from which chromosome was this read was derived
			   int ReadLoci,	// read start loci
			   etSeqBase *pRead); // read sequence

int		// returned minimum Watson Hamming
HammingDistCntsW(int ReadLen,		// read length
			 int ChromID,     // from which chromosome was this read was derived
			 int ReadLoci,	// read start loci
			 etSeqBase *pRead, // read sequence
			 UINT32 *pHammDist);	// where to return hamming dist counts (assumes at least ReadLen elements)

int		// returned minimum Crick Hamming
HammingDistCntsC(int ReadLen,		// read length
			   int ChromID,     // from which chromosome was this read was derived
			   int ReadLoci,	// read start loci
			   etSeqBase *pRead, // read sequence
   			 UINT32 *pHammDist);	// where to return hamming dist counts (assumes at least ReadLen elements)


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
	return _T("kangasr");
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

int Idx;
int LenReq;

int PMode;					// processing mode
int FMode;					// output format - csv loci, csv loci +seq or multifasta
char Strand;				// generate for this strand '+' or '-' or for both '*'
int NumReads;				// number of reads required
int ReadLen;				// read lengths
int CutMin;					// min length
int CutMax;					// max cut length
bool bDedupe;				// true if unique reads only to be generated
int DfltHamming;			// if >= 0 then the default Hamming edit distance to use
int NumberOfProcessors;		// number of installed CPUs
int NumThreads;				// number of threads (0 defaults to number of CPUs)
bool bReadHamDist;			// true if hamming distributions from each sampled read to all other genome subsequences to be generated
int SNPrate;				// generate SNPs at this rate per million bases

int Region;					// Process regions 0:ALL,1:CDS,2:5'UTR,3:3'UTR,4:Introns,5:5'US,6:3'DS,7:Intergenic (default = ALL)
int UpDnRegLen;				// if processing regions then up/down stream regulatory length

double PropRandReads;		// proportion of reads to be generated with completely random sequences

etSEMode SEMode;			// induced sequencer error rate mode
int DistCluster;			// cluster generated reads into windows of this median length, 0 if no clustering
double SeqErrRate;			// dynamic sequencer error rate
int InDelSize;				// simulated InDel size range
double InDelRate;			// simulated InDel rate per read
bool bSeqErrProfile;		// true if to generate composite sequencer errors with uniform profile (default is Illumina 3' skewed)
bool bPEgen;				// true if paired ends are to be generated
int PEmin;					// PE minimum insert
int PEmax;					// PE maximum insert
char szInFile[_MAX_PATH];	// input from this bioseq assembly
char szProfFile[_MAX_PATH];// input from this profile site preferences file
char szOutFile[_MAX_PATH];	// output simulated reads to this file
char szOutPEFile[_MAX_PATH];// output partner paired end simulated reads to this file
char szSNPFile[_MAX_PATH];	// output simulated SNPs to this file
char szHammFile[_MAX_PATH];	// hamming edit distance file
char szFeatFile[_MAX_PATH]; // optional BED feature or gene file if generating transcriptome reads

double Artef5Rate;			// rate (0..1) at which to insert 5' artefact sequences
int NumArtef5Seqs;			// number of user specified 5' artefact sequences
char *pszArtef5Seqs[cMaxArtefSeqs]; // 5' artefact sequences
double Artef3Rate;			// rate (0..1) at which to insert 3' artefact sequences
int NumArtef3Seqs;			// number of user specified 3' artefact sequences
char *pszArtef3Seqs[cMaxArtefSeqs]; // 5' artefact sequences

// command line args
struct arg_lit  *help    = arg_lit0("h","help",                "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to screen and logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"diagnostics log file");

struct arg_int *pmode = arg_int0("m","mode","<int>",		    "processing mode: 0 - random start and end, 1 - Profiled start with random end sites, 2 - random start with profiled end sites,  3 - profiled start with profiled end sites, 4 - same as standard with Hammings (0 - default)");
struct arg_int *fmode = arg_int0("M","format","<int>",		    "output format: 0 - CSV loci only, 1 - CSV loci with sequence, 2 - multifasta wrapped, 3 - multifasta non-wrapped, 4 - SOLiD colorspace, 5 - SOLiD for BWA (default: 3)");
struct arg_int *numreads = arg_int0("n","nreads","<int>",	    "number of reads required, minimum 100, maximum 500000000 (default = 10000000)");


struct arg_int *region    = arg_int0("G","genomicregion","<int>","Process regions 0:ALL,1:CDS,2:5'UTR,3:3'UTR,4:Introns,5:5'US,6:3'DS,7:Intergenic (default = ALL)");
struct arg_int *updnreglen    = arg_int0("u","updnreglen","<int>","Up/Dn stream regulatory region length, used if processing regions (default is 2000)");


struct arg_int *generrmode = arg_int0("g","generrmode","<int>", "simulate sequencer error modes: 0 - no errors, 1 - induce fixed num errors per read, 2 - static profile, 3 - dynamic according to '-z<rate>' (0 - default)");
struct arg_dbl *seqerrrate = arg_dbl0("z","seqerrs","<dbl>",	"simulate composite sequencer errors with induced error mean rate in range 0-0.20%");
struct arg_lit *seqerrprof = arg_lit0("Z","seqerrprofile",      "generate composite sequencer errors with uniform profile (default is Illumina 3' skewed)");

struct arg_int *indelsize = arg_int0("x","indelsize","<int>",   "simulate micro-InDel size range: 1..9 of this length maximum (default is 3)");
struct arg_dbl *indelrate = arg_dbl0("X","indelrate","<dbl>",	"simulate micro-InDels with mean rate per read in range 0 - 100% of all reads (default is 0)");

struct arg_str *strand=arg_str0("s", "strand","<str>",          "generate for this strand '+' or '-' only (default is both)");
struct arg_int *readlen = arg_int0("l","length","<int>",	    "read lengths (default = 100, max is 1000)");
struct arg_int *cutmin = arg_int0("c","cutmin","<int>",		    "min cut length (minimum = read length)");
struct arg_int *cutmax = arg_int0("C","cutmax","<int>",		    "max cut length (maximum = 1000)");
struct arg_file *infile = arg_file1("i","in","<file>",			"input from this raw multifasta or bioseq assembly");
struct arg_file *inmnase = arg_file0("I","in","<file>",			"input from this profile site preferences file");

struct arg_int *distcluster = arg_int0("D","distcluster","<int>","distribute generated reads as clusters into windows of this median length (0..300), 0 if no clustering");

struct arg_file *featfile = arg_file0("t","featfile","<file>",	"use features or genes in this BED file to generate transcriptome reads from target genome");

struct arg_file *hammfile = arg_file0("H","hammingfile","<file>","use distances in this Hamming edit distance file (.hmg) for targeted genome instead of dynamic generation");
struct arg_file *outfile = arg_file1("o","out","<file>",		"output simulated (or N/1 if paired) reads to this file");
struct arg_file *outpefile = arg_file0("O","outpe","<file>",	"output simulated (N/2) paired end reads to this file");
struct arg_file *outsnpfile = arg_file0("u","outsnp","<file>",	"output simulated SNP loci to this BED file");
struct arg_int *threads = arg_int0("T","threads","<int>",		"number of processing threads 0..n (defaults to 0 which sets threads to number of CPUs)");
struct arg_lit  *dedupe = arg_lit0("d","dedupe",                "generate unique read sequences only");
struct arg_int *hamming = arg_int0("e","hamming","<int>",		"if specified and < 0, then dynamically generate Hamming edit distances, otherwise use this static distance (default = static generation with Hamming 0)");
struct arg_lit  *readhamdist = arg_lit0("r","readhamdist",      "generate hamming distribution from each simulated read to all other subsequences of same length in genome");
struct arg_int *snprate = arg_int0("N","snprate","<int>",       "generate SNPs at the specified rate per Mb (default = 0, max = 20000)");

struct arg_lit  *pegen = arg_lit0("p","pegen",				    "generate paired end reads");
struct arg_int *pemin = arg_int0("j","pemin","<int>",           "generate paired end reads with minimum fragment lengths (default is 2x read length)");
struct arg_int *pemax = arg_int0("J","pemax","<int>",           "generate paired end reads with maximum  fragment lengths  (default is min fragment length plus read length)");

struct arg_dbl *proprandreads = arg_dbl0("R","randreads","<dbl>","proportion (0 to 0.9000) of reads to be generated with random sequences not likely to be in target genome (default = 0.0)");

struct arg_dbl *artif5rate = arg_dbl0("a","artif5rate","<dbl>",	"randomly induce sequencer adaptor/linker artefacts at 5' end of sequences rate (default 0.0, max 0.9)");
struct arg_str *artif5str = arg_strn("A","artif5str","<string>",0,cMaxArtefSeqs,"5' artefacts sequence(s) (default is 'ACACTCTTTCCCTACACGACGCTGTTCCATCT'");

struct arg_dbl *artif3rate = arg_dbl0("b","artif3rate","<dbl>",	"randomly induce sequencer adaptor/linker artefacts at 3' end of sequences rate (default 0.0, max 0.9)");
struct arg_str *artif3str = arg_strn("A","artif3str","<string>",0,cMaxArtefSeqs,"3' artefacts sequences (default is 'ACACTCTTTCCCTACACGACGCTCTTCCGATCT'");

struct arg_end *end = arg_end(20);

void *argtable[] = {help,version,FileLogLevel,LogFile,
					pmode,region,updnreglen,pegen,pemin,pemax,readhamdist,fmode,numreads,proprandreads,snprate,generrmode,distcluster,seqerrrate,seqerrprof,
					artif5rate,artif5str,artif3rate,artif3str,
					indelsize,indelrate,strand,readlen,cutmin,cutmax,dedupe,hamming,featfile,
					infile,inmnase,hammfile,outpefile,outfile,outsnpfile,threads,
					end};

char **pAllArgs;
int argerrors;
argerrors = CUtility::arg_parsefromfile(argc,(char **)argv,&pAllArgs);
if(argerrors >= 0)
	argerrors = arg_parse(argerrors,pAllArgs,argtable);

/* special case: '--help' takes precedence over error reporting */
if (help->count > 0)
        {
		printf("\n%s the Kanga Simulated reads generator, Version %s\nOptions ---\n", gszProcName,cpszProgVer);
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

		// now that log parameters have been parsed then initialise diagnostics log system
	if(!gDiagnostics.Open(szLogFile,(etDiagLevel)iScreenLogLevel,(etDiagLevel)iFileLogLevel,true))
		{
		printf("\nError: Unable to start diagnostics subsystem\n");
		if(szLogFile[0] != '\0')
			printf(" Most likely cause is that logfile '%s' can't be opened/created\n",szLogFile);
		exit(1);
		}
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Version: %s",cpszProgVer);

	szFeatFile[0] = '\0';
	szSNPFile[0] = '\0';

	PMode = (etPMode)(pmode->count ? pmode->ival[0] : ePMStandard);
	if(PMode < ePMStandard || PMode >= ePMplaceholder)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Processing mode '-m%d' specified outside of range %d..%d",PMode,(int)ePMStandard,(int)ePMplaceholder-1);
		exit(1);
		}

	bReadHamDist = readhamdist->count ? true : false;
	if(PMode == ePMSampHamm)
		bReadHamDist = true;

	if(bReadHamDist)
		{
		Region = region->count ? region->ival[0] : eMEGRAny;	// default as being any region
		if(Region < eMEGRAny)
			{
			printf("\nSpecified region '-G%d' < 0, assuming you meant ANY region",Region);
			Region = eMEGRAny;
			}
		else
			{
			if(Region > eMEGRIntergenic)
				{
				printf("\nSpecified region '-n%d' > %d, assuming you meant Intergenic",Region,eMEGRIntergenic);
				Region = eMEGRIntergenic;

				}
			}
		if(Region != eMEGRAny)
			{
			UpDnRegLen = updnreglen->count ? updnreglen->ival[0] : cUpdnstream;
			if(UpDnRegLen < 0 || UpDnRegLen > 100000)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Up/Dn stream regulatory region '-u%d' specified outside the range 0..100000",UpDnRegLen);
				exit(1);
				}
			}
		else
			UpDnRegLen = 0;
		}
	else
		{
		Region = eMEGRAny;
		UpDnRegLen = 0;
		}

	InDelRate = 0.0f;
	InDelSize = 0;
	if(!bReadHamDist)
		{
		bPEgen = pegen->count ? true : false;
		if(bPEgen && PMode != ePMStandard)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Processing mode '-m%d' not supported when generating paired ends",PMode);
			exit(1);
			}

		if(bPEgen && outpefile->count == 0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Paired end processing '-p' requested but no partner paired end output file specifed with '-O<file>'");
			exit(1);
			}

		FMode = (etFMode)(fmode->count ? fmode->ival[0] : eFMNWFasta);
		if(FMode < eFMcsv || FMode >= eFMplaceholder)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Output format '-m%d' specified outside of range %d..%d",PMode,(int)eFMcsv,(int)eFMplaceholder-1);
			exit(1);
			}

		if(bPEgen && !(FMode == eFMNWFasta || FMode == eFMFasta))
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Output format '-M%d' not supported (only fasta output supported) when generating paired ends",FMode);
			exit(1);
			}

		PropRandReads = proprandreads->count ? proprandreads->dval[0] : 0.0;
		if(PropRandReads < 0.0 || PropRandReads >= 0.9000)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Proportion of random reads '-R%f' must in range 0 to 0.9000",PropRandReads);
			exit(1);
			}

		SNPrate = snprate->count ? snprate->ival[0] : 0;
		if(SNPrate < 0 || SNPrate > 20000)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: simulated SNP rate '-N%d' specified outside of range 0..20000",SNPrate);
			exit(1);
			}

		if(outsnpfile->count)
			{
			strncpy(szSNPFile,outsnpfile->filename[0],_MAX_PATH);
			szSNPFile[_MAX_PATH-1] = '\0';
			if(SNPrate == 0)
				SNPrate = 1000;
			}
		else
			szSNPFile[0] = '\0';

		DistCluster = distcluster->count ? distcluster->ival[0] : 0;
		if(DistCluster < 0 || DistCluster > cMaxDistClusterLen)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: distributed cluster size '-D%d' must be in range 0 to %d",DistCluster,cMaxDistClusterLen);
			exit(1);
			}

		InDelRate = indelrate->count ? indelrate->dval[0] : 0.0f;
		if(InDelRate < 0.0 || InDelRate > 1.0)
			{
			printf("\nError: simulated InDel size range '-X%f' must be in range 0.0 to 1.0",InDelRate);
			exit(1);
			}
		if(InDelRate > 0.0)
			{
			InDelSize = indelsize->count ? indelsize->ival[0] :3;
			if(InDelSize < 1 || InDelSize > 9)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: simulated InDel size range '-x%d' specified outside of range 1..9",InDelSize);
				exit(1);
				}
			}
		else
			InDelSize = 0;

		SEMode = (etSEMode)(generrmode->count ? generrmode->ival[0] : eSEPnone);
		if(SEMode < eSEPnone || SEMode >= eSEPplaceholder)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: simulated error rate mode '-m%d' specified outside of range %d..%d",SEMode,(int)eSEPnone,(int)eSEPplaceholder-1);
			exit(1);
			}
		if(SEMode == eSEPdyn || SEMode == eSEPfixerrs)
			{
			SeqErrRate = SEMode == eSEPdyn ? 0.01 : 5.0;
			double MaxSeqErrRate = SEMode == eSEPdyn ? 0.20 : 30;
			SeqErrRate = seqerrrate->count ? seqerrrate->dval[0] : SeqErrRate;
			if(SeqErrRate < 0.0 || SeqErrRate > MaxSeqErrRate)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: sequencer error rate '-z%f' must be in range 0.0 to %f",SeqErrRate,MaxSeqErrRate);
				exit(1);
				}
			}
		else
			SeqErrRate = -1;

		if((FMode < eFMNWFasta) && SEMode != eSEPnone)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: sequencer error rate mode '-g%d' only allowed if generating multifasta output with '-M2/3/4'",SEMode);
			exit(1);
			}

		bSeqErrProfile = seqerrprof->count ? true : false;
		if(bSeqErrProfile && SEMode == eSEPnone)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Uniform profile '-Z' requested but no error rate mode specified with '-g<mode>'");
			exit(1);
			}

		}
	else
		{
		SEMode = eSEPnone;
		FMode = eFMcsv;
		bPEgen = false;
		DistCluster = 0;
		SNPrate = 0;
		InDelSize = 0;
		PropRandReads = 0.0;
		SeqErrRate = 0;
		bSeqErrProfile = false;
		}


	if(bReadHamDist)
		Strand = '+';
	else
		{
		Strand = strand->count ? *(char *)strand->sval[0] : '*';
		if(!(Strand == '+' || Strand == '-' || Strand == '*'))
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Strand specified '-s%c' must be one of '+', '-' or '*'",Strand);
			exit(1);
			}
		}

	NumReads = numreads->count ? numreads->ival[0] : cDfltNumReads;
	if(NumReads < cMinNumReads || NumReads > cMaxNumReads)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Number of reads '-n%d' specified outside of range %d..%d",NumReads,cMinNumReads,cMaxNumReads);
		exit(1);
		}

	ReadLen = readlen->count ? readlen->ival[0] : cDfltReadLen;
	if(ReadLen < cMinReadLen || ReadLen > cMaxSimReadLen)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Read length '-a%d' specified outside of range %d..%d",ReadLen,cMinReadLen,cMaxSimReadLen);
		exit(1);
		}

	if(!bReadHamDist) // if not generating Hammings to all other subsequences
		{
		DfltHamming = hamming->count ? hamming->ival[0] : 0;
		if(DfltHamming > ReadLen)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Hamming edit distance '-e%d' must be <= read length %d",DfltHamming,ReadLen);
			exit(1);
			}
		else
			if(DfltHamming < 0)
				DfltHamming = -1;

		if((DfltHamming == -1 || hamming->count == 0) && hammfile->count)
			{
			strncpy(szHammFile,hammfile->filename[0],_MAX_PATH);
			szHammFile[_MAX_PATH-1] = '\0';
			DfltHamming = -1;
			}
		else
			szHammFile[0] = '\0';

		CutMin = cutmin->count ? cutmin->ival[0] : ReadLen;
		if(CutMin < ReadLen || CutMin > cMaxCutLen)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Minimum cut length '-c%d' must be in range %d..%d",CutMin,ReadLen,cMaxCutLen);
			exit(1);
			}

		CutMax = cutmax->count ? cutmax->ival[0] : CutMin;
		if(CutMax < CutMin || CutMax > cMaxCutLen)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Maximum cut length '-C%d' must be in range %d..%d",CutMax,CutMin,cMaxCutLen);
			exit(1);
			}


		if(bPEgen)
			{
			PEmin = pemin->count ? pemin->ival[0] : ReadLen * 2;
			if(PEmin < cMinPEFragLen || PEmin > cMaxPEFragLen)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: paired end minimum fragment length '-j%d' specified outside of range %d..%d",PEmin,cMinPEFragLen,cMaxPEFragLen);
				exit(1);
				}
			PEmax = min(PEmin + ReadLen,cMaxPEFragLen);
			PEmax = pemax->count ? pemax->ival[0] : PEmax;
			if(PEmax < PEmin || PEmax > cMaxPEFragLen)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: paired end maximum fragment length '-J%d' specified outside of range %d..%d",PEmax,PEmin,cMaxPEFragLen);
				exit(1);
				}

			if(PEmin <= CutMin)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: paired end min fragment length (%d) must be more than min read length (%d)",PEmin,CutMin);
				exit(1);
				}

			}
		else
			{
			PEmin = 0;
			PEmax = 0;
			}
		}
	else
		{
		CutMin = ReadLen;
		CutMax = ReadLen;
		DfltHamming = 0;
		szHammFile[0] = '\0';
		PEmin = 0;
		PEmax = 0;
		}

	bDedupe = dedupe->count ? true : false;

	if(bPEgen && bDedupe)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Deduping not supported when generating paired ends");
		exit(1);
		}

	if(PMode != ePMStandard && inmnase->count)
		{
		strncpy(szProfFile,inmnase->filename[0],_MAX_PATH);
		szProfFile[_MAX_PATH-1] = '\0';
		}
	else
		szProfFile[0] = '\0';

	strcpy(szInFile,infile->filename[0]);

	strcpy(szOutFile,outfile->filename[0]);

	if(featfile->count)
		{
		strncpy(szFeatFile,featfile->filename[0],_MAX_PATH);
		szFeatFile[_MAX_PATH-1] = '\0';
		FMode = eFMcsv;
		}
	else
		szFeatFile[0] = '\0';

	if(Region != eMEGRAny && szFeatFile[0] == '\0')
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Filtering on region specified with '-G%d' but no BED feature file specified with '-t<file>'",Region);
		exit(1);
		}

	if(bPEgen)
		strcpy(szOutPEFile,outpefile->filename[0]);
	else
		szOutPEFile[0] = '\0';


	int MaxArtefLen = min(ReadLen - 5,cMaxArtefSeqLen);
	if(artif5rate->count ||  artif3rate->count)
		{
		if(!(FMode == eFMNWFasta || FMode == eFMFasta ) || bPEgen || Region != eMEGRAny || bDedupe || bReadHamDist)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Sorry artefact sequence processing not supported with currently selected options - only single end fasta is supported!");
			exit(1);
			}
		}


	if(FMode == eFMNWFasta || FMode == eFMFasta)
		{
		Artef5Rate = artif5rate->count ? artif5rate->dval[0] : 0.0;
		if(Artef5Rate < 0.0 || Artef5Rate > 1.0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: 5' artefact rate '-a%1.2f' must be in the range 0.0 to 1.0",Artef5Rate);
			exit(1);
			}
		}
	else
		Artef5Rate = 0.0;

	if(Artef5Rate > 0.0)
		{
		NumArtef5Seqs = artif5str->count;
		if(NumArtef5Seqs > 0)
			{
			for(Idx=0;Idx < artif5str->count; Idx++)
				{
				LenReq = (int)strlen(artif5str->sval[Idx]);
				pszArtef5Seqs[Idx] = new char [LenReq+1];
				strcpy(pszArtef5Seqs[Idx],artif5str->sval[Idx]);
				if((LenReq = TrimSeqChrs(pszArtef5Seqs[Idx])) < 1 || LenReq > MaxArtefLen)
					{
					if(LenReq == -1)
						gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: 5' artefact sequence '-A\"%s\"' contains non-sequence char",artif5str->sval[Idx]);
					else
						{
						if(LenReq == 0)
							gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: 5' artefact sequence '-A\"%s\"' after filtering is empty",artif5str->sval[Idx]);
						else
							gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: 5' artefact sequence '-A\"%s\"' is too long (must be <= %d bp)",pszArtef5Seqs[Idx],MaxArtefLen);
						}

					exit(1);
					}
				}
			}
		else
			{
			pszArtef5Seqs[0] = new char [cMaxArtefSeqLen+1];
			strncpy(pszArtef5Seqs[0],pszArtef5Seq,MaxArtefLen);
			pszArtef5Seqs[0][MaxArtefLen] = '\0';
			NumArtef5Seqs = 1;
			}
		}
	else
		{
		NumArtef5Seqs = 0;
		pszArtef5Seqs[0] = NULL;
		}

	if(FMode == eFMNWFasta || FMode == eFMFasta)
		{
		Artef3Rate = artif3rate->count ? artif3rate->dval[0] : 0.0;
		if(Artef3Rate < 0.0 || Artef3Rate > 1.0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: 5' artefact rate '-b%1.2f' must be in the range 0.0 to 1.0",Artef3Rate);
			exit(1);
			}
		}
	else
		Artef3Rate = 0.0;



	if(Artef3Rate > 0.0)
		{
		NumArtef3Seqs = artif3str->count;
		if(NumArtef3Seqs > 0)
			{
			for(Idx=0;Idx < artif3str->count; Idx++)
				{
				LenReq = (int)strlen(artif3str->sval[Idx]);
				pszArtef3Seqs[Idx] = new char [LenReq+1];
				strcpy(pszArtef3Seqs[Idx],artif3str->sval[Idx]);
				if((LenReq = TrimSeqChrs(pszArtef3Seqs[Idx])) < 1 || LenReq > MaxArtefLen)
					{
					if(LenReq == -1)
						gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: 3' artefact sequence '-A\"%s\"' contains non-sequence char",artif3str->sval[Idx]);
					else
						{
						if(LenReq == 0)
							gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: 3' artefact sequence '-A\"%s\"' after filtering is empty",artif3str->sval[Idx]);
						else
							gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: 3' artefact sequence '-A\"%s\"' is too long (must be <= %d bp)",pszArtef3Seqs[Idx],MaxArtefLen);
						}

					exit(1);
					}
				}
			}
		else
			{
			pszArtef3Seqs[0] = new char [cMaxArtefSeqLen+1];
			strncpy(pszArtef3Seqs[0],pszArtef3Seq,MaxArtefLen);
			pszArtef3Seqs[0][MaxArtefLen] = '\0';
			NumArtef3Seqs = 1;
			}
		}
	else
		{
		NumArtef3Seqs = 0;
		pszArtef3Seqs[0] = NULL;
		}

#ifdef _WIN32
	SYSTEM_INFO SystemInfo;
	GetSystemInfo(&SystemInfo);
	NumberOfProcessors = SystemInfo.dwNumberOfProcessors;
#else
	NumberOfProcessors = sysconf(_SC_NPROCESSORS_CONF);
#endif
	int MaxAllowedThreads = min(cMaxWorkerThreads,NumberOfProcessors);	// limit to be at most cMaxWorkerThreads
	if((NumThreads = threads->count ? threads->ival[0] : MaxAllowedThreads)==0)
		NumThreads = MaxAllowedThreads;
	if(NumThreads < 0 || NumThreads > MaxAllowedThreads)
		{
		gDiagnostics.DiagOut(eDLWarn,gszProcName,"Warning: Number of threads '-T%d' specified was outside of range %d..%d",NumThreads,1,MaxAllowedThreads);
		gDiagnostics.DiagOut(eDLWarn,gszProcName,"Warning: Defaulting number of threads to %d",MaxAllowedThreads);
		NumThreads = MaxAllowedThreads;
		}

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing parameters:");

	const char *pszDescr;
	switch(PMode) {
		case ePMStandard:
			pszDescr = "random start and random end sites";
			break;
		case ePMProfRand:
			pszDescr = "profiled start with random end sites";
			break;
		case ePMRandProf:
			pszDescr = "profiled start with random end sites";
			break;
		case ePMProfProf:
			pszDescr = "profiled start with profiled end sites";
			break;
		case ePMSampHamm:
			pszDescr = "random start and random end sites with Hamming generation in 'uhamming' format";
			break;
		}
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"processing mode is : '%s'",pszDescr);

	switch(FMode) {
		case eFMcsv:
			pszDescr = "CSV loci only";
			break;
		case eFMcsvSeq:
			pszDescr = "CSV loci + sequence";
			break;
		case eFMFasta:
			pszDescr = "Multifasta, wrapped read sequences";
			break;
		case eFMNWFasta:
			pszDescr = "Multifasta, non-wrapped read sequences";
			break;
		case eFMSOLiD:
			pszDescr = "SOLiD colorspace csfasta";
			break;
		case eFMSOLiDbwa:
			pszDescr = "SOLiD colorspace reads in double encoded basespace to suit BWA in fastq";
			break;
		}


	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Generate reads as : '%s'", bPEgen ? "paired ends" : "single ended");
	if(bPEgen)
		gDiagnostics.DiagOutMsgOnly(eDLInfo," Paired end fragment length: min %d, max %d", PEmin, PEmax);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Output format is : '%s'",pszDescr);

	switch(SEMode) {
		case eSEPnone:
			pszDescr = "No induced simulated errors";
			break;
		case eSEPfixerrs:
			pszDescr = "simulate fixed number of errors in each read at rate specified by '-z<rate>'";
			break;
		case eSEPstatic:
			pszDescr = "use internal static profile";
			break;
		case eSEPdyn:
			pszDescr = "dynamic according to '-z<rate>'";
			break;
		}
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Simulated sequencer error mode : '%s'",pszDescr);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Simulate SNPs at this rate per million bases : %d",SNPrate);
	if(szSNPFile[0] != '\0')
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Write simulated SNP loci to this BED file : '%s'",szSNPFile);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Generate random reads in this Proportion: %0.3f",PropRandReads);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Cluster reads into these sized windows : %d", DistCluster);

	if(SEMode == eSEPdyn)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Induce composite sequence errors at mean rate of : %f",SeqErrRate);
	else
		if(SEMode == eSEPfixerrs)
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Induce fixed number of sequence errors : %d",(int)SeqErrRate);
	if(SEMode != eSEPnone)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Induce sequence errors profile : '%s'",bSeqErrProfile ? "uniform" : "3' skewed");

	if(InDelRate == 0.0f)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Do not induce InDels");
	else
		{
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Induce InDels at mean rate per read of : %1.3f",InDelRate);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"InDels will have a size range of : 1 to %d",InDelSize);
		}
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Simulate for this strand : '%c'",Strand);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Number of %s reads required: %d",bPEgen ?  "paired end" : "single ended", NumReads);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Read lengths: %d",ReadLen);

	if(Artef5Rate > 0.0 && NumArtef5Seqs > 0)
		{
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Induce 5' artefacts at this rate: %1.2f",Artef5Rate);
		for(Idx = 0; Idx < NumArtef5Seqs; Idx++)
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"   5' artefact (%d) : '%s'",Idx+1,pszArtef5Seqs[Idx]);
		}
	else
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Induce 5' artefacts at this rate: %1.2f",Artef5Rate);

	if(Artef3Rate > 0.0 && NumArtef3Seqs > 0)
		{
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Induce 3' artefacts at this rate: %1.2f",Artef3Rate);
		for(Idx = 0; Idx < NumArtef3Seqs; Idx++)
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"   3' artefact (%d) : '%s'",Idx+1,pszArtef3Seqs[Idx]);
		}
	else
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Induce 3' artefacts at this rate: %1.2f",Artef3Rate);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Generate Hamming distributions for each simulated read: %s",bReadHamDist?"Yes":"No");

	if(!bReadHamDist)
		{
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"min cut length: %d",CutMin);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"max cut length: %d",CutMax);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"unique read sequences only: %s",bDedupe ? "yes" : "no");
		if(DfltHamming == -1)
			{
			if(szHammFile[0] == '\0')
				gDiagnostics.DiagOutMsgOnly(eDLInfo,"Dynamically generate Hamming edit distances (VERY SLOW!)");
			else
				gDiagnostics.DiagOutMsgOnly(eDLInfo,"Use Hamming edit distances from this file: '%s'",szHammFile);
			}
		else
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Default Hamming edit distances: %d",DfltHamming);
		}

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"input bioseq assembly file: '%s'",szInFile);
	if(PMode != ePMStandard)
		{
		if(szProfFile[0] == '\0')
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"static profiled site preferences derived from site octamer");
		else
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"input profiled site preferences file: '%s'",szProfFile);
		}
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"output simulated reads file: '%s'",szOutFile);
	if(bPEgen)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"output simulated paired read partners to file: %s", szOutPEFile);

	if(szFeatFile[0] != '\0')
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"input BED feature or gene file: '%s'",szFeatFile);

	if(Region != eMEGRAny)
		{
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Process for region: %s",Region2Txt((etBEDRegion)Region));
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Up/Dn stream regulatory length: %d",UpDnRegLen);
		}

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"number of threads : %d",NumThreads);
#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	gStopWatch.Start();
	Rslt = Process((etPMode)PMode,SEMode,bPEgen,PEmin,PEmax,PropRandReads,
			DistCluster,SeqErrRate,bSeqErrProfile,SNPrate,
						InDelSize,InDelRate,bReadHamDist,(etFMode)FMode,NumThreads,Strand,NumReads,
						ReadLen,Artef5Rate,NumArtef5Seqs,pszArtef5Seqs,Artef3Rate,NumArtef3Seqs,pszArtef3Seqs,
						CutMin,CutMax,bDedupe,DfltHamming,Region,UpDnRegLen,szFeatFile,szInFile,szProfFile,szHammFile,szOutPEFile,szOutFile,szSNPFile);
	gStopWatch.Stop();
	Rslt = Rslt >=0 ? 0 : 1;
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Exit code: %d Total processing time: %s",Rslt,gStopWatch.Read());
	exit(Rslt);
	}
else
	{
	printf("\n%s Kanga Simulated reads generator, Version %s\n",gszProcName,cpszProgVer);
	arg_print_errors(stdout,end,gszProcName);
	arg_print_syntax(stdout,argtable,"\nUse '-h' to view option and parameter usage\n");
	exit(1);
	}
return 0;
}

// TrimSeqChrs
// Removes any quotes and whitespace (space and tabs only) from pszTxt
// Also checks for legal base chars 'acgt'
// Returns -1 if illegal char, 0 if empty sequence, otherwise the sequence length
int
TrimSeqChrs(char *pszTxt)
{
char *pDst = pszTxt;
char Chr;
int Len = 0;
while((Chr = *pszTxt++))
	{
	if(Chr == '\0' || (Chr == '\t'  || Chr == ' '  || Chr == '"' || Chr == '\''))
		continue;
	switch(Chr) {
		case 'a': case 'A':
			Chr = 'a';
			break;
		case 'c': case 'C':
			Chr = 'c';
			break;
		case 'g': case 'G':
			Chr = 'g';
			break;
		case 't': case 'T': case 'u': case 'U':
			Chr = 't';
			break;
		default:
			return(-1);
		}
	*pDst++ = Chr;
	Len++;
	}
*pDst = Chr;
return(Len);
}

char *
Region2Txt(etBEDRegion Region)
{
switch(Region) {
	case eMEGRAny:		// process any region
		return((char *)"All");

	case eMEGRCDS:	// only process CDS
		return((char *)"CDS");

	case eMEGR5UTR:	// only process 5' UTR
		return((char *)"5' UTR");

	case eMEGR3UTR:	// only process 3' UTR
		return((char *)"3' UTR");

	case eMEGRIntrons:		// only process Introns
		return((char *)"Introns");

	case eMEGRUpstream:		// only process 5'upstream regulatory
		return((char *)"5'US");

	case eMEGRDnstream:		// only process 3'upstream regulatory
		return((char *)"3'DS");

	case eMEGRIntergenic:		// only process intergenic
		return((char *)"Intergenic");

	default:
		break;
	}
return((char *)"Unsupported");
}

typedef struct TAG_sChromSeq {
	int ChromSeqID;			// m_pChromSeqs[ChromSeqID-1] of this tsChromSeq instance
	int ChromID;
	char Strand;			// '*' if chrom seq is for genome assembly, '+' or '-' if seq is for a feature or gene
	int RelDensity;			// relative (0 to 1000), to other reads, compared with rand(0,1000) and if less then no read generated in current iteration
	int ScaledStart;		// scaled start of this chrom
	int ScaledLen;			// chrom length scaled such that the sum of all chrom scaled lengths is less than INT_MAX
	char szChromName[cMaxDatasetSpeciesChrom];
	int Len;				// actual length
	UINT32 SeqOfs;			// offset at which this chroms sequence starts
} tsChromSeq;


#pragma pack(4)
typedef struct TAG_sSimRead {
	int ChromSeqID;			// simulated read is on this m_pChromSeqs[]
	int Status;				// 0 if yet to be reported on, 1 if reported, 2 if not reported because dup read seq, 3 if not reported because from same loci
	int FlgPE2:1;			// if paired reads generation then set if this is the 3' end
	struct TAG_sSimRead *pPartner; // if paired reads generation then points to partner of this read
	int ChromID;			// chromosome identifier
	int Strand;				// on this strand (0 if '+', 1 if '-')
	int StartLoci;			// starts at this loci
	int EndLoci;			// ends at this loci
	int Len;				// and was originally of this length
	int InDelLen;			// and has had a deletion (<0) or an insertion (>0) of this length
	int Lenx;				// resulting after the InDel in the simulated read being of this length
	int HammingDist;		// read is at least this hamming distance from any other sequence of same length in targeted genome
	etSeqBase *pSeq;		// pts to start of this simulated reads sequence
	UINT32 *pHamDistFreq;	// hamming distributions from this read to all other genome subsequences of same length
	} tsSimRead;

typedef struct TAG_sWorkerPars {
#ifdef _WIN32
	HANDLE threadHandle;	// handle as returned by _beginthreadex()
	unsigned int threadID;	// identifier as set by _beginthreadex()
#else
	int threadRslt;			// result as returned by pthread_create ()
	pthread_t threadID;		// identifier as set by pthread_create ()
#endif
	int ThreadIdx;			// index of this thread 0..n
	int RandSeed;			// random generator seed for this worker
	etPMode PMode;			// processing mode
	int Region;				//  Process regions 0:ALL,1:CDS,2:5'UTR,3:3'UTR,4:Introns,5:5'US,6:3'DS,7:Intergenic (default = ALL)
	int UpDnRegLen;			// if processing regions then up/down stream regulatory length
	int NumReqReads;		// number of reads to be generated by this worker thread
	int NumGenReads;		// number of actually generated reads
	int DfltHamming;		// if < 0 then dynamically determine Hamming distance otherwise use this value as the Hamming
	bool bUseLocateHamming; // if true then call LocateHamming() using distances from file instead of dynamically generating Hammings
	bool bDedupe;			// if true then slough any reads with a Hamming of 0
	bool bReadHamDist;		// true if hamming distributions from each sampled read to all other genome subsequences to be generated
	int ReadLen;			// read length
	int CutMin;				// min cut length
	int CutMax;				// max cut length
	char Strand;			// generate for this strand '+' or '-' or for both '*'
	bool bPEgen;			// true if paired ends are to be generated
	bool bMaxIters;			// set true by thread if exhusted attempts to find chrom from which reads can be simulated
	int PEmin;				// PE minimum fragment
	int PEmax;				// PE maximum fragment
	tsSimRead *pReads;		// to hold all reads generated by this worker thread - will have been preallocated to hold NumReqReads
} tsWorkerPars;
#pragma pack()

etPMode m_PMode;				// processing mode
etFMode m_FMode;				// output format
int m_hOutFile;					// output results file handle
char *m_pszOutFile;				// output file name

int m_hOutPEFile;					// output results paired end file handle
char *m_pszOutPEFile;				// output paired end file name

CCSVFile *m_pProfCSV;			// used to hold profile site preferences whilst loading into m_pProfSel
double *m_pProfSel;				// allocated array of profile site selection preferences (0.0..1.0) indexed by sequence octamers

int m_PropRandReads;			// what proportion ( * 1 million) are to be generated as completely random reads
int m_DistCluster;				// cluster generated reads into windows of this median length, 0 if no clustering
int m_InDelSize;				// simulated InDel size range
double m_InDelRate;				// simulated InDel rate per read


double m_Artef5Rate;			// rate (0..1) at which to insert 5' artefact sequences
int m_NumArtef5Seqs;			// number of user specified 5' artefact sequences
char **m_ppszArtef5Seqs;		// 5' artefact sequences
double m_Artef3Rate;			// rate (0..1) at which to insert 3' artefact sequences
int m_NumArtef3Seqs;			// number of user specified 3' artefact sequences
char **m_ppszArtef3Seqs;		// 5' artefact sequences

int m_Artef5SeqLens[cMaxArtefSeqs];				// sequence lengths for each 5' artefact
etSeqBase m_Artef5Seqs[cMaxArtefSeqs][cMaxArtefSeqLen+1];	// each 5' artefact sequence
int m_Artef3SeqLens[cMaxArtefSeqs];				// sequence lengths for each 3' artefact
etSeqBase m_Artef3Seqs[cMaxArtefSeqs][cMaxArtefSeqLen+1];	// each 3' artefact sequence

CBioSeqFile *m_pBioSeqFile;		// genome assembly
CBEDfile *m_pBEDFile;			// optional if simulating transcript reads

char m_szSpecies[cMaxDatasetSpeciesChrom+1];		// species (title) as retrieved from bioseqfile
int m_NumChromSeqs;				// number of chromosomes loaded
int m_AllocdChromSeqs;			// number allocated
size_t m_AllocdChromSize;		// allocd mem size for m_pChromSeqs
tsChromSeq *m_pChromSeqs;		// pts to chromseqs array

UINT8 *m_pGenomeSeq;			// allocated to hold concatenated (separated by eBaseEOSs) chromosome sequences
size_t m_AllocdGenomeMem;		// allocd mem size for m_pGenomeSeq
INT64 m_GenomeLen;				// total genome length including concatenators
INT32 m_GenomeScaledLen;		// sum of all chrom scaled lengths, will always be less than INT_MAX

tsHamHdr *m_pHamHdr;			// header for binary format hamming edit distances
tsHamChrom *m_pCurHamChrom;		// pts to current chromosome specific binary hammings

tsSimRead *m_pSimReads;			// allocated to hold simulated reads - NOTE: mmap/malloc used because of GNU malloc limitations
INT64 m_AllocdMemReads;			// actual memory allocation size used when allocating m_pSimReads
int m_NumReadsAllocd;			// this many reads have been allocated for in m_pSimReads
int m_TotReqReads;				// number of reads required to be simulated - will be 2x user requested number if simulating paired end reads
int m_CurNumGenReads;			// current number of generated reads - updated every N reads generated by worker threads
UINT32 *m_pHamDistFreq;			// allocated to hold hamming distance counts from one read to all other genome subsequences

int m_MaxFastaLineLen;			// wrap sequences in multifasta output files if line is longer than this many bases

#ifdef _WIN32
HANDLE m_hMtxDedupe;
unsigned __stdcall WorkerThread(void *pThreadPars);
#else
pthread_mutex_t m_hMtxDedupe;
void *WorkerThread(void * pThreadPars);
#endif

TRandomCombined<TRanrotWGenerator,TRandomMersenne> RGseeds((int)time(0));


int
SimInDels(tsSimRead *pSimRead,int *pReadLen,etSeqBase *pRead)
{
int Idx;
etSeqBase InsertBases[20];
bool bInsert;
double Thres;
int InDelSize;
int InDelPsn;
if(m_InDelRate == 0.0 || m_InDelSize == 0)
	return(0);
Thres = RGseeds.Random();	// pick a number, any number, between 0 and 1 inclusive
if(Thres > m_InDelRate)		// InDel this read?
	return(0);

bInsert = RGseeds.IRandom(0,1) == 0 ? false : true;			// pick insert or delete
if(bInsert)
	{
	for(Idx = 0; Idx < m_InDelSize; Idx++)
		InsertBases[Idx] = (etSeqBase) RGseeds.IRandom(0,3);
	InDelSize = (int)RGseeds.IRandom(1,m_InDelSize);            // pick the insertion size
	InDelPsn = (int)RGseeds.IRandom(0,*pReadLen-InDelSize);     // pick the insert position
	memcpy(&pRead[InDelPsn+InDelSize],&pRead[InDelPsn],*pReadLen-(InDelPsn+InDelSize));
	memcpy(&pRead[InDelPsn],InsertBases,InDelSize);
	pSimRead->InDelLen = InDelSize;
	return(-1 * InDelSize);
	}
else
	{
	InDelSize = (int)RGseeds.IRandom(1,m_InDelSize);            // pick the deletion size
	InDelPsn = (int)RGseeds.IRandom(0,*pReadLen-InDelSize);     // pick the deletion position
	memcpy(&pRead[InDelPsn],&pRead[InDelPsn+InDelSize],*pReadLen-(InDelPsn+InDelSize));
	pSimRead->InDelLen = -1 * InDelSize;
	return(InDelSize);
	}
}

int
SimArtefacts(bool b3ArtefSeq,		// if false then 5' artefact else if true then 3' artefact
			int ReadLen,etSeqBase *pRead)
{
int NumArtSeqs;
double ArtefRate;

int ArtefSeqIdx;
int ArtefSeqLen;
int ArtefactLen;
etSeqBase *pArtefSeq;

double Thres;

if(b3ArtefSeq)
	{
	ArtefRate = m_Artef3Rate;
	NumArtSeqs = m_NumArtef3Seqs;
	}
else
	{
	ArtefRate = m_Artef5Rate;
	NumArtSeqs = m_NumArtef5Seqs;
	}

if(NumArtSeqs == 0 || ArtefRate == 0.0)
	return(ReadLen);

Thres = RGseeds.Random();	// pick a number, any number, between 0 and 1 inclusive
if(ArtefRate < 1.0 && Thres > ArtefRate)	// 5' artefact this read?
	return(0);

if(NumArtSeqs > 1)
	ArtefSeqIdx = RGseeds.IRandom(0,NumArtSeqs);			// pick artefact sequence
else
	ArtefSeqIdx = 0;
if(b3ArtefSeq)
	{
	pArtefSeq = m_Artef3Seqs[ArtefSeqIdx];
	ArtefSeqLen = m_Artef3SeqLens[ArtefSeqIdx];
	}
else
	{
	pArtefSeq = m_Artef5Seqs[ArtefSeqIdx];
	ArtefSeqLen = m_Artef5SeqLens[ArtefSeqIdx];
	}

ArtefactLen = min(ArtefSeqLen,ReadLen - 10);
ArtefactLen = RGseeds.IRandom(1,ArtefSeqLen);

if(b3ArtefSeq)
	memcpy(&pRead[ReadLen - ArtefactLen],pArtefSeq,ArtefactLen);
else
	{
	memcpy(&pRead[ArtefactLen],pRead,ReadLen - ArtefactLen);
	memcpy(pRead,&pArtefSeq[ArtefSeqLen - ArtefactLen],ArtefactLen);
	}

return(ReadLen);
}

// Induced errors will be heavily 3' biased simulating Illumina read substitution profiles
typedef struct TAG_sInduceErrProf {
	double Proportion;			// proportion of reads (0 - 100.0) with
	int NumSubs;				// this number of sequencer errors
} tsInduceErrProf;

// static profile - sequencing error rates as initially hardcoded
// updated (15th August 2013) to better represent the observed error rates with current Illumina 100bp read sets
// assumes mean substitution rate per 100bp of 0.01%
// =POISSON.DIST(NumSubs,MeanExpectedSubs=1,FALSE)
tsInduceErrProf StaticErrProfile[] =
	{
	{ 0.367879, 0},	// proportion of reads with no substitutions
	{ 0.367879, 1},	// proportion of reads  with 1 substitution
	{ 0.183944, 2},	// proportion of reads   with 2 substitutions
	{ 0.061313, 3},	// proportion of reads  with 3 substitutions
	{ 0.015328, 4},	// proportion of reads  with 4 substitutions
	{ 0.003066, 5},	// proportion of reads   with 5 substitutions
	{ 0.000511, 6},	// proportion of reads  with 6 substitutions
	{ 0.000073, 7},	// proportion of reads  with 7 substitutions
	{ -1.0f, 8 }	// remainder of reads  with 8 substitutions
	};
#define cNumProfEls (int)(sizeof(StaticErrProfile)/sizeof(tsInduceErrProf))

// dynamic profile - initialised according to user specified induced error rates
tsInduceErrProf DynErrProfile[cNumProfEls];


// Illumina cumulative error profile along length of read with moderate increase in subs at 5' start of reads but most subs are down at the 3' end of reads
int IlluminaSpatialDist[20] = { 40,55,64,72,80,88,96,104,112,121,131,142,156,174,197,228,270,325,400,500 };
#define cNumIlluminaSpatialDist (int)(sizeof(IlluminaSpatialDist)/sizeof(int))  

etSEMode m_SEMode;				// simulated sequencer error mode
double m_DynProbErr = 0.01;		// default as being a 1% error rate
bool m_bUniformDist = true;		// default as being a uniform distribution of induced errors
int m_PrvProfSeqLen = -1;		// set to read length for which DynErrProfile has been generated
int m_InducedErrDist[cNumProfEls];	// to hold induced error count distribution
int m_InducedErrPsnDist[2000];  // to hold read sequence psn induced error count


// CAUTION: uses the static profile for
// determining composite read distribution
int // number of substitutions inplace induced into this read
SimSeqErrors(int SeqLen,etSeqBase *pRead)
{
int Idx;
int RandSubs;
int SubBase;
int NumSubs2Induce;
etSeqBase *pSeq;
double Thres;

UINT8 Subd[2000];
if(m_SEMode == eSEPnone)
	return(0);
pSeq = pRead;
for(Idx = 0; Idx < SeqLen; Idx++,pSeq++)
	*pSeq = *pSeq & ~cRptMskFlg;

// determine how many errors are to be induced into the current read
tsInduceErrProf *pProfile;
switch(m_SEMode) {
	case eSEPdyn:
		if(m_PrvProfSeqLen != SeqLen)
			{
			pProfile = &DynErrProfile[0];
			double CurThres = pow(1.0-m_DynProbErr,(double)SeqLen);
			double AccThres = 0.0;
			for(Idx = 0; Idx < (cNumProfEls-1); Idx++,pProfile++)
				{
				pProfile->NumSubs = Idx;
				pProfile->Proportion = CurThres;
				AccThres += CurThres;
				CurThres = (1-AccThres)/2;
				}
			pProfile->NumSubs = cNumProfEls-1;
			pProfile->Proportion = -1.0;		// flags last
			m_PrvProfSeqLen = SeqLen;
			}
		pProfile = &DynErrProfile[0];
		break;
	case eSEPstatic:
		pProfile = StaticErrProfile;
		break;
	case eSEPfixerrs:
		if((int)m_DynProbErr <= 0)
			{
			m_InducedErrDist[0] += 1;
			return(0);
			}
		pProfile = NULL;
		break;
	}

if(pProfile != NULL)
	{
	Thres = RGseeds.Random();	// pick a number, any number, between 0 and 1 inclusive
	for(Idx = 0; Idx < cNumProfEls; Idx++, pProfile++)
		{
		if(pProfile->Proportion <= 0.0 || pProfile->Proportion >= Thres)
			break;
		Thres -= pProfile->Proportion;
		}
	if(!pProfile->NumSubs)
		{
		m_InducedErrDist[0] += 1;
		return(0);
		}
	}
if(pProfile != NULL)
	NumSubs2Induce = pProfile->NumSubs;
else
	NumSubs2Induce = (int)m_DynProbErr;
memset(Subd,0,SeqLen);
RandSubs = 0;

int Psn;
while(RandSubs < NumSubs2Induce)
	{
	if(m_bUniformDist)
		Idx = RGseeds.IRandom(0,SeqLen);
	else
		{
		int MinPsn;
		int MaxPsn;
		int DistIdx;
		Psn = RGseeds.IRandom(0,IlluminaSpatialDist[cNumIlluminaSpatialDist-1]);
		for(DistIdx = 0; DistIdx < cNumIlluminaSpatialDist-1; DistIdx++)
			{
			if(Psn <= IlluminaSpatialDist[DistIdx])
				break;
			}
		MinPsn = (DistIdx * SeqLen) / cNumIlluminaSpatialDist;
		if(DistIdx == cNumIlluminaSpatialDist - 1)
			MaxPsn = SeqLen-1;
		else
			MaxPsn = (MinPsn + (SeqLen / cNumIlluminaSpatialDist)) - 1;

		Psn = RGseeds.IRandom(MinPsn,MaxPsn);
		Idx=min(SeqLen-1,Psn);
		}
	if(Subd[Idx])
		continue;
	pSeq = &pRead[Idx];
	do {
		SubBase = RGseeds.IRandom(0,3);
		}
	while(SubBase == *pSeq);
	*pSeq = SubBase | cRptMskFlg;
	Subd[Idx] = 1;
	m_InducedErrPsnDist[Idx] += 1;
	RandSubs += 1;
	}
m_InducedErrDist[NumSubs2Induce] += 1;
return(NumSubs2Induce);
}

//
// Randomises all bases in this sequence making it highly unlikely that it can be aligned
int // number of substitutions inplace induced into this read
SimSeqRand(int SeqLen,etSeqBase *pRead)
{
int Idx;
int SubBase;
etSeqBase *pSeq;
if(m_SEMode == eSEPnone)
	return(0);
pSeq = pRead;
for(Idx = 0; Idx < SeqLen; Idx++,pSeq++)
	{
	*pSeq = *pSeq & ~cRptMskFlg;
	while(1)
		{
		SubBase = RGseeds.IRandom(0,3);
		if(SubBase != *pSeq)
			{
			*pSeq = SubBase;
			break;
			}
		}
	};
return(SeqLen);
}


void
Reset(bool bSync)
{
if(m_hOutFile != -1)
	{
	if(bSync)
#ifdef _WIN32
		_commit(m_hOutFile);
#else
		fsync(m_hOutFile);
#endif
	close(m_hOutFile);
	m_hOutFile = -1;
	}
if(m_hOutPEFile != -1)
	{
	if(bSync)
#ifdef _WIN32
		_commit(m_hOutPEFile);
#else
		fsync(m_hOutPEFile);
#endif
	close(m_hOutPEFile);
	m_hOutPEFile = -1;
	}
if(m_pBioSeqFile != NULL)
	{
	delete m_pBioSeqFile;
	m_pBioSeqFile = NULL;
	}
if(m_pBEDFile != NULL)
	{
	delete m_pBEDFile;
	m_pBEDFile = NULL;
	}
if(m_pProfSel != NULL)
	{
	delete m_pProfSel;
	m_pProfSel = NULL;
	}
if(m_pProfCSV != NULL)
	{
	delete m_pProfCSV;
	m_pProfCSV = NULL;
	}

if(m_pChromSeqs != NULL)
	{
#ifdef _WIN32
	free(m_pChromSeqs);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pChromSeqs != MAP_FAILED)
		munmap(m_pChromSeqs,m_AllocdChromSize);
#endif
	m_pChromSeqs = NULL;
	}

if(m_pSimReads != NULL)
	{
#ifdef _WIN32
	free(m_pSimReads);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pSimReads != MAP_FAILED)
		munmap(m_pSimReads,m_AllocdMemReads);
#endif
	m_pSimReads = NULL;
	}

if(m_pGenomeSeq != NULL)
	{
#ifdef _WIN32
	free(m_pGenomeSeq);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pGenomeSeq != MAP_FAILED)
		munmap(m_pGenomeSeq,m_AllocdGenomeMem);
#endif
	m_pGenomeSeq = NULL;
	}
if(m_pHamHdr != NULL)
	{
	delete(m_pHamHdr);
	m_pHamHdr = NULL;
	}

if(m_pHamDistFreq != NULL)
	{
	delete m_pHamDistFreq;
	m_pHamDistFreq = NULL;
	}

m_pCurHamChrom = NULL;
m_AllocdMemReads = 0;
m_NumReadsAllocd = 0;

m_GenomeLen = 0;
m_AllocdGenomeMem = 0;
m_GenomeScaledLen = 0;
m_szSpecies[0] = '\0';
m_NumChromSeqs = 0;
m_AllocdChromSeqs = 0;
m_AllocdChromSize = 0;

m_DynProbErr = 0.01;	// default as being a 1% error rate
m_bUniformDist = false;		// default as being a uniform distribution of induced errors
m_InDelSize = 3;				// simulated InDel size range
m_InDelRate = 0.0;				// simulated InDel rate per read
m_PrvProfSeqLen = -1;  // set to read length for which DynErrProfile has been generated
m_TotReqReads = 0;
m_CurNumGenReads = 0;
m_MaxFastaLineLen = 79;
memset(m_InducedErrDist,0,sizeof(m_InducedErrDist));	// to hold induced error count distribution
memset(m_InducedErrPsnDist,0,sizeof(m_InducedErrPsnDist));	// to hold read sequence psn induced error count
}

void
Init(void)
{
m_pBioSeqFile = NULL;
m_pBEDFile = NULL;
m_pChromSeqs = NULL;
m_pProfSel = NULL;			// allocated array of profile site selection preferences (0.0..1.0) indexed by sequence octamers
m_pProfCSV = NULL;			// used to hold DNase site preferences whilst loading into m_pProfSel
m_hOutFile = -1;
m_hOutPEFile = -1;
m_pSimReads = NULL;			// allocated to hold simulated reads
m_pGenomeSeq = NULL;		// allocated to hold concatenated (separated by eBaseEOSs) chromosome sequences
m_pHamHdr = NULL;
m_pCurHamChrom = NULL;
m_pHamDistFreq = NULL;
Reset(false);
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
	Reset(false);
	return(eBSFerrOpnFile);
	}

if(read(hHamFile,&HamHdr,sizeof(tsHamHdr))!=sizeof(tsHamHdr))
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to read %s - %s",pszHammings,strerror(errno));
	close(hHamFile);
	Reset(false);
	return(eBSFerrParse);
	}
if(HamHdr.Magic[0] != 'b' && HamHdr.Magic[0] != 'h' && HamHdr.Magic[0] != 'a' && HamHdr.Magic[0] != 'm')
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Not a binary format Hamming distance file - %s",pszHammings);
	close(hHamFile);
	Reset(false);
	return(eBSFerrFileType);
	}
if(HamHdr.NumChroms < 1)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Hamming distance file contains no edit distances - %s",pszHammings);
	close(hHamFile);
	Reset(false);
	return(eBSFerrNoEntries);
	}

if((m_pHamHdr = (tsHamHdr *)new UINT8 [HamHdr.Len])==NULL)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Unable to allocate memory (%d bytes) for holding Hamming distances loaded from - %s",HamHdr.Len,pszHammings);
	close(hHamFile);
	Reset(false);
	return(eBSFerrMem);
	}
memcpy(m_pHamHdr,&HamHdr,sizeof(tsHamHdr));
if(read(hHamFile,(UINT8 *)m_pHamHdr+sizeof(tsHamHdr),HamHdr.Len-sizeof(tsHamHdr))!=HamHdr.Len-sizeof(tsHamHdr))
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Unable to read all Hamming edit distances from - %s",pszHammings);
	close(hHamFile);
	Reset(false);
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

int				// returned Hamming distance, -1 if no hammings loaded, -2 if chrom not matched, -3 if loci outside of range
LocateHamming(char *pszChrom,UINT32 Loci)
{
int ChromIdx;
tsHamChrom *pCurHamChrom;
if(m_pHamHdr == NULL)
	return(-1);
for(ChromIdx = 0; ChromIdx < m_pHamHdr->NumChroms; ChromIdx++)
	{
	pCurHamChrom = (tsHamChrom *)((UINT8 *)m_pHamHdr + m_pHamHdr->ChromOfs[ChromIdx]);
	if(!stricmp((char *)pCurHamChrom->szChrom,pszChrom))
		break;
	}
if(ChromIdx == m_pHamHdr->NumChroms)
	return(-2);
if(Loci >= pCurHamChrom->NumEls)
	return(-3);
return(pCurHamChrom->Dists[Loci]);
}


// generate '+' strand index from K-mer of length SeqLen starting at pSeq
int
GenPSeqIdx(int SeqLen,etSeqBase *pSeq)
{
int Idx;
int SeqIdx;
int Base;
for(Idx=SeqIdx=0; Idx < SeqLen; Idx++,pSeq++)
	{
	Base = *pSeq & NUCONLYMSK;
	if(Base > eBaseT)
		return(-1);
	SeqIdx <<= 2;
	SeqIdx |= Base;
	}
return(SeqIdx);
}

// generate '-' strand index from K-mer of length SeqLen starting at pSeq
int
GenMSeqIdx(int SeqLen,etSeqBase *pSeq)
{
int Idx;
int SeqIdx;
int Base;

pSeq += SeqLen-1;
for(Idx=SeqIdx=0; Idx < SeqLen; Idx++,pSeq--)
	{
	Base = *pSeq & NUCONLYMSK;
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



int
InitProfSitePrefs(char *pszInProfFile)	// read from this profile site selectivity file (for MNase, generated by MNaseSitePred process), if NULL then static profiling
{
int Rslt;
int NumProcessed;
int NumFields;
int OctIdx;
char *pszOctamer;
etSeqBase Octamer[9];
double SitePref;

if(pszInProfFile == NULL || pszInProfFile[0] == '\0')
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Static profile site selection preferences based on site octamer");
else
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Loading profile site selection preferences from file: %s",pszInProfFile);

m_pProfSel = (double *) new double [0x010000];	// to hold 4^8 (octamer) site preferences
if(m_pProfSel == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadReads: Memory allocation for profile site preferences of 64k doubles failed");
	Reset(false);
	return(eBSFerrMem);
	}

if(pszInProfFile == NULL || pszInProfFile[0] == '\0')
	{
	for(OctIdx = 0; OctIdx < 0x010000; OctIdx++)
		m_pProfSel[OctIdx] = max(0.0000001f,(double)OctIdx/(double)0x0ffff);
	return(eBSFSuccess);
	}

for(OctIdx = 0; OctIdx < 0x010000; OctIdx++)
	m_pProfSel[OctIdx] = 0.0f;

if((m_pProfCSV = new CCSVFile) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CCSVfile");
	Reset(false);
	return(eBSFerrObj);
	}

if((Rslt=m_pProfCSV->Open(pszInProfFile))!=eBSFSuccess)
	{
	while(m_pProfCSV->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,m_pProfCSV->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open file: %s",pszInProfFile);
	Reset(false);
	return(Rslt);
	}

NumProcessed = 0;
while((Rslt=m_pProfCSV->NextLine()) > 0)	// onto next line containing fields
	{
	NumFields = m_pProfCSV->GetCurFields();
	if(NumFields < 4)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"file: %s contains % fields, expected at least 4",pszInProfFile,NumFields);
		Reset(false);
		return(eBSFerrParams);
		}
	if(!NumProcessed && m_pProfCSV->IsLikelyHeaderLine())
		continue;
	NumProcessed += 1;
	m_pProfCSV->GetText(1,&pszOctamer);
	memcpy(Octamer,CSeqTrans::MapAscii2Sense(pszOctamer),8);
	OctIdx = GenPSeqIdx(8,Octamer);
	m_pProfCSV->GetDouble(4,&SitePref);
	m_pProfSel[OctIdx] = SitePref;
	}

delete m_pProfCSV;
m_pProfCSV = NULL;

return(eBSFSuccess);
}

const int cMaxAllocBuffChunk = 0x07fffff;			// buffer fasta
const int cAllocNumChroms = 4096;					// allocate for chrom seqs in this many increments

int
LoadFasta(size_t *pTotLen,int MinChromLen,char *pszFastaFile)
{
CFasta Fasta;
unsigned char *pSeqBuff;
char szName[cBSFSourceSize];
char szDescription[cBSFDescriptionSize];
int SeqLen;
int Descrlen;
bool bFirstEntry;
bool bEntryCreated;
int Rslt;
int SeqID;
tsChromSeq *pChromSeq;
UINT8 *pSeq;
int ChromID;
int SeqOfs;
int NumLenWarnings;
size_t TotLen;
tsChromSeq *pTmp;
size_t memreq;

m_GenomeLen = 0;
m_AllocdGenomeMem = 0;
TotLen = 0;
m_NumChromSeqs = 0;
m_GenomeLen = 0;
ChromID = 0;

m_AllocdChromSize = sizeof(tsChromSeq) * cAllocNumChroms;
#ifdef _WIN32
m_pChromSeqs = (tsChromSeq *) malloc((size_t)m_AllocdChromSize);
if(m_pChromSeqs == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadFasta: Memory allocation of %lld bytes failed",(INT64)m_AllocdChromSize);
	m_AllocdChromSize = 0;
	Reset(false);
	return(eBSFerrMem);
	}
#else
// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
m_pChromSeqs = (tsChromSeq *)mmap(NULL,(size_t)m_AllocdChromSize, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
if(m_pChromSeqs == MAP_FAILED)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadFasta: Memory allocation of %lld bytes through mmap()  failed",(INT64)m_AllocdChromSize,strerror(errno));
	m_pChromSeqs = NULL;
	m_AllocdChromSize = 0;
	Reset(false);
	return(eBSFerrMem);
	}
#endif
m_AllocdChromSeqs = cAllocNumChroms;
memset(m_pChromSeqs,0,m_AllocdChromSize);

m_AllocdGenomeMem = cMaxAllocBuffChunk;		// an initial allocation to hold the assembly sequence, will be extended as may be required
#ifdef _WIN32
m_pGenomeSeq = (UINT8 *) malloc((size_t)m_AllocdGenomeMem);
if(m_pGenomeSeq == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadFasta: Memory allocation of %lld bytes failed",(INT64)m_AllocdGenomeMem);
	Reset(false);
	return(eBSFerrMem);
	}
#else
// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
m_pGenomeSeq = (UINT8 *)mmap(NULL,(size_t)m_AllocdGenomeMem, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
if(m_pGenomeSeq == MAP_FAILED)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadFasta: Memory allocation of %lld bytes through mmap()  failed",(INT64)m_AllocdGenomeMem,strerror(errno));
	m_pGenomeSeq = NULL;
	Reset(false);
	return(eBSFerrMem);
	}
#endif


if((pSeqBuff = new unsigned char [cMaxAllocBuffChunk]) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessFastaFile:- Unable to allocate memory (%d bytes) for sequence buffer",cMaxAllocBuffChunk);
	Reset(false);
	return(eBSFerrMem);
	}

if((Rslt=Fasta.Open(pszFastaFile,true))!=eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open '%s' [%s] %s",pszFastaFile,Fasta.ErrText((teBSFrsltCodes)Rslt),Fasta.GetErrMsg());
	Reset(false);
	delete pSeqBuff;
	return(Rslt);
	}
bFirstEntry = true;
bEntryCreated = false;
SeqID = 0;
ChromID = 0;
TotLen = 0;

#ifdef _WIN32
_ASSERTE( _CrtCheckMemory());
#endif
NumLenWarnings = 0;
while((Rslt = SeqLen = Fasta.ReadSequence(pSeqBuff,cMaxAllocBuffChunk-1,true,false)) > eBSFSuccess)
	{
	if(bFirstEntry || SeqLen == eBSFFastaDescr)		// just read a descriptor line
		{
		if(!bFirstEntry)
			{
			pChromSeq->Len = SeqOfs;
			if(pChromSeq->Len < MinChromLen && NumLenWarnings++ < 10)
				gDiagnostics.DiagOut(eDLWarn,gszProcName,"LoadFasta: Chrom '%s' of length %d may be too short from which to reliably sample simulated reads...",pChromSeq->szChromName,pChromSeq->Len);
			m_pGenomeSeq[pChromSeq->SeqOfs+SeqOfs] = eBaseEOS;
			m_GenomeLen += SeqOfs + 1;
			}

		if(ChromID == m_AllocdChromSeqs)
			{
			memreq = sizeof(tsChromSeq) * (m_AllocdChromSeqs + cAllocNumChroms);
#ifdef _WIN32
			pTmp = (tsChromSeq *) realloc(m_pChromSeqs,memreq);
#else
			pTmp = (tsChromSeq *)mremap(m_pChromSeqs,m_AllocdChromSize,memreq,MREMAP_MAYMOVE);
			if(pTmp == MAP_FAILED)
				pTmp = NULL;
#endif
			if(pTmp == NULL)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadFasta: Memory re-allocation to %lld bytes - %s",(INT64)(memreq),strerror(errno));
				delete pSeqBuff;
				return(eBSFerrMem);
				}
			m_pChromSeqs = pTmp;
			m_AllocdChromSize = memreq;
			m_AllocdChromSeqs += cAllocNumChroms;
			}

     	pChromSeq = &m_pChromSeqs[ChromID++];
		pChromSeq->ChromSeqID = ChromID;
		pChromSeq->ChromID = ChromID;
		m_NumChromSeqs = ChromID;
		Descrlen = Fasta.ReadDescriptor(szDescription,cBSFDescriptionSize);
		if(SeqLen != eBSFFastaDescr || sscanf(szDescription," %s[ ,]",szName)!=1)
			sprintf(szName,"%s.%d",pszFastaFile,ChromID);
		strncpy(pChromSeq->szChromName,szName,sizeof(pChromSeq->szChromName));
		pChromSeq->Strand = '*';
		pChromSeq->RelDensity = 1000;
		pChromSeq->SeqOfs = (UINT32)m_GenomeLen;
		bEntryCreated = true;
		bFirstEntry = false;
		SeqOfs = 0;
		if(SeqLen == eBSFFastaDescr)
			continue;
		}

	// realloc m_pGenomeSeq as and when required
	if((m_GenomeLen + SeqOfs + SeqLen + 100) >= (INT64)m_AllocdGenomeMem)
		{
		memreq = m_AllocdGenomeMem + SeqLen + cMaxAllocBuffChunk;
#ifdef _WIN32
		pSeq = (UINT8 *) realloc(m_pGenomeSeq,memreq);
#else
		pSeq = (UINT8 *)mremap(m_pGenomeSeq,m_AllocdGenomeMem,memreq,MREMAP_MAYMOVE);
		if(pSeq == MAP_FAILED)
			pSeq = NULL;
#endif
		if(pSeq == NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadFasta: Memory re-allocation to %lld bytes - %s",(INT64)(memreq),strerror(errno));
			delete pSeqBuff;
			return(eBSFerrMem);
			}
		m_pGenomeSeq = pSeq;
		m_AllocdGenomeMem = memreq;
		}

#ifdef _WIN32
_ASSERTE( _CrtCheckMemory());
#endif

	pSeq = &m_pGenomeSeq[pChromSeq->SeqOfs+SeqOfs];
	memcpy(pSeq,pSeqBuff,SeqLen);
	SeqOfs += SeqLen;
	TotLen += SeqLen;
	}

if(Rslt < eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessFastaFile [%s] %s",Fasta.ErrText((teBSFrsltCodes)Rslt),Fasta.GetErrMsg());
	delete pSeqBuff;
	return(false);
	}
if(bEntryCreated)
	{
	pChromSeq->Len = SeqOfs;
	if(pChromSeq->Len < MinChromLen && NumLenWarnings++ < 10)
		gDiagnostics.DiagOut(eDLWarn,gszProcName,"LoadFasta: Chrom '%s' of length %d may be too short from which to reliably sample simulated reads...",pChromSeq->szChromName,pChromSeq->Len);

	m_pGenomeSeq[pChromSeq->SeqOfs+SeqOfs] = eBaseEOS;
	m_GenomeLen += SeqOfs + 1;
	}
delete pSeqBuff;
*pTotLen = TotLen;
return(eBSFSuccess);
}

int
LoadBioseq(size_t *pTotLen,int MinChromLen,char *pszBioSeqFile)
{
int Rslt;
size_t TotLen;
int ChromID;
etSeqBase *pSeq;
tsChromSeq *pChromSeq;
int Len;
int NumLenWarnings;

*pTotLen = 0;
if((m_pBioSeqFile = new CBioSeqFile) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CBioSeqFile object");
	return(eBSFerrObj);
	}
if((Rslt = m_pBioSeqFile->Open(pszBioSeqFile))!=eBSFSuccess)
	{
	if(Rslt == eBSFerrFileAccess)
		return(Rslt);

	while(m_pBioSeqFile->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,m_pBioSeqFile->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open assembly sequence file '%s'",pszBioSeqFile);
	Reset(false);
	return(Rslt);
	}

m_pBioSeqFile->GetTitle(sizeof(m_szSpecies),m_szSpecies);

m_GenomeLen = 0;
m_AllocdGenomeMem = 0;
TotLen = 0;
m_NumChromSeqs = 0;
ChromID = 0;
m_AllocdChromSeqs = m_pBioSeqFile->NumEntries();

m_AllocdChromSize = sizeof(tsChromSeq) * m_AllocdChromSeqs;
#ifdef _WIN32
m_pChromSeqs = (tsChromSeq *) malloc((size_t)m_AllocdChromSize);
if(m_pChromSeqs == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadAssembly: Memory allocation of %lld bytes failed",(INT64)m_AllocdChromSize);
	Reset(false);
	return(eBSFerrMem);
	}
#else
// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
m_pChromSeqs = (tsChromSeq *)mmap(NULL,(size_t)m_AllocdChromSize, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
if(m_pChromSeqs == MAP_FAILED)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadAssembly: Memory allocation of %lld bytes through mmap()  failed",(INT64)m_AllocdChromSize,strerror(errno));
	m_pChromSeqs = NULL;
	Reset(false);
	return(eBSFerrMem);
	}
#endif

memset(m_pChromSeqs,0,m_AllocdChromSize);

// determine total sequence length required for allocating to hold genome as one concatenated sequence
NumLenWarnings = 0;
while((ChromID = m_pBioSeqFile->Next(ChromID)) > 0) {
	pChromSeq = &m_pChromSeqs[ChromID-1];
	pChromSeq->ChromID = ChromID;
	m_pBioSeqFile->GetName(ChromID,sizeof(pChromSeq->szChromName),pChromSeq->szChromName);
	pChromSeq->Strand = '*';
	pChromSeq->RelDensity = 1000;
	Len = m_pBioSeqFile->GetDataLen(ChromID);
	if(Len < MinChromLen && NumLenWarnings++ < 10)
		gDiagnostics.DiagOut(eDLWarn,gszProcName,"LoadAssembly: Chrom '%s' of length %d may be too short from which to reliably sample simulated reads...",pChromSeq->szChromName,Len);
	m_GenomeLen += (INT64)Len + 1;		// allow for a concatenation separator
	pChromSeq->Len = Len;
	pChromSeq->SeqOfs = 0;
	pChromSeq->ChromSeqID = ++m_NumChromSeqs;
	}
m_AllocdGenomeMem = (size_t)m_GenomeLen;
#ifdef _WIN32
m_pGenomeSeq = (UINT8 *) malloc((size_t)m_AllocdGenomeMem);
if(m_pGenomeSeq == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadAssembly: Memory allocation of %lld bytes failed",(INT64)m_AllocdGenomeMem);
	Reset(false);
	return(eBSFerrMem);
	}
#else
// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
m_pGenomeSeq = (UINT8 *)mmap(NULL,(size_t)m_AllocdGenomeMem, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
if(m_pGenomeSeq == MAP_FAILED)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadAssembly: Memory allocation of %lld bytes through mmap()  failed",(INT64)m_AllocdGenomeMem,strerror(errno));
	m_pGenomeSeq = NULL;
	Reset(false);
	return(eBSFerrMem);
	}
#endif

ChromID = 0;
TotLen = 0;
pSeq = m_pGenomeSeq;
while((ChromID = m_pBioSeqFile->Next(ChromID)) > 0) {
	pChromSeq = &m_pChromSeqs[ChromID-1];
	pChromSeq->SeqOfs = (UINT32)(UINT64)(pSeq - m_pGenomeSeq);
	TotLen += pChromSeq->Len;
	m_pBioSeqFile->GetData(ChromID,eSeqBaseType,0,pSeq,pChromSeq->Len);
	pSeq += pChromSeq->Len;
	*pSeq++ = eBaseEOS;		// concatenator separator
	}
delete m_pBioSeqFile;
m_pBioSeqFile = NULL;
*pTotLen = TotLen;
return(eBSFSuccess);
}

int
LoadGenome(int MinChromLen,			// warn if loaded chromosomes are less than this length; may be too short to sample reads from
			char *pszBioSeqFile)
{
int Rslt;
size_t TotLen;
int ChromID;
tsChromSeq *pChromSeq;
double LenSCF;			// length scaling factor
int CurScaledStart;

// assume a bioseq file, if not a bioseq file then try to load as a raw multifasta
if((Rslt = LoadBioseq(&TotLen,MinChromLen,pszBioSeqFile)) != eBSFSuccess)
	{
	if(Rslt != eBSFerrFileAccess)
		return(Rslt);

	if((Rslt = LoadFasta(&TotLen,MinChromLen,pszBioSeqFile)) != eBSFSuccess)
		return(Rslt);
	}

// assembly or multifasta sequences now loaded
LenSCF = (double)INT_MAX/(double)TotLen;
m_GenomeScaledLen = TotLen >= (INT64)INT_MAX ? (long)(TotLen * LenSCF) : (long)TotLen;
CurScaledStart = 0;
pChromSeq = &m_pChromSeqs[0];
for(ChromID = 0; ChromID < m_NumChromSeqs; ChromID++,pChromSeq++)
	{
	pChromSeq->ScaledLen = TotLen >= (INT64)INT_MAX ? (int)(pChromSeq->Len * LenSCF) : pChromSeq->Len;
	pChromSeq->ScaledStart = CurScaledStart;
	CurScaledStart += pChromSeq->ScaledLen;
	}
return(eBSFSuccess);
}


int
LoadTranscriptome(char *pszBioSeqFile,			// genome assembly
				char *pszFeatFile)				// BED file containing features or genes
{
int Rslt;
bool bTranscribed;
int Len;
size_t TotLen;
int ChromID;
int FeatID;
tsChromSeq *pChromSeq;
etSeqBase *pSeq;
double LenSCF;			// length scaling factor
int CurScaledStart;

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
	Reset(false);
	return(Rslt);
	}

m_pBioSeqFile->GetTitle(sizeof(m_szSpecies),m_szSpecies);

// open feature file
if((m_pBEDFile = new CBEDfile) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CBEDfile object");
	return(eBSFerrObj);
	}
if((Rslt = m_pBEDFile->Open(pszFeatFile))!=eBSFSuccess)
	{
	while(m_pBEDFile->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,m_pBEDFile->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open feature/gene BED file '%s'",pszFeatFile);
	Reset(false);
	return(Rslt);
	}

// check if BED contains features as genes (transcribed lengths) or simply features as being regions in the genome with no transcriptional information
bTranscribed = m_pBEDFile->ContainsGeneDetail();

// alloc for total number of features, sum of transcribed length if genes or sum of feature lengths if regions
m_AllocdChromSeqs = m_pBEDFile->GetNumFeatures();
if((m_pChromSeqs = new tsChromSeq [m_AllocdChromSeqs])==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadTranscriptome: unable to allocate memory (%d bytes) for %d Feature Seqs",sizeof(tsChromSeq) * m_AllocdChromSeqs,m_AllocdChromSeqs);
	Reset(false);
	return(eBSFerrMem);
	}
memset(m_pChromSeqs,0,sizeof(tsChromSeq) * m_AllocdChromSeqs);

// determine total sequence length required for allocating to hold all features or genes as one concatenated sequence
m_GenomeLen = 0;
TotLen = 0;
m_NumChromSeqs = 0;

int CurFeatScore;
int MaxFeatScore;
char CurStrand;

MaxFeatScore = 0;
FeatID = 0;
while((FeatID = m_pBEDFile->GetNextFeatureID(FeatID)) > 0) {
	pChromSeq = &m_pChromSeqs[FeatID-1];
	pChromSeq->ChromID = FeatID;
	m_pBEDFile->GetFeature(FeatID,pChromSeq->szChromName,NULL,NULL,NULL,&CurFeatScore,&CurStrand);
	pChromSeq->Strand = CurStrand;
	CurFeatScore = min(CurFeatScore,1000);
	MaxFeatScore = max(MaxFeatScore,CurFeatScore);
	if(bTranscribed)
		Len = m_pBEDFile->GetTranscribedLen(FeatID);			// get transcribed length for gene
	else
		Len = m_pBEDFile->GetFeatLen(FeatID);
	m_GenomeLen += (INT64)Len + 1;		// allow for a concatenation separator
	pChromSeq->Len = Len;
	pChromSeq->SeqOfs = 0;
	pChromSeq->ChromSeqID = ++m_NumChromSeqs;
	}

FeatID = 0;
while((FeatID = m_pBEDFile->GetNextFeatureID(FeatID)) > 0) {
	pChromSeq = &m_pChromSeqs[FeatID-1];
	m_pBEDFile->GetFeature(FeatID,NULL,NULL,NULL,NULL,&CurFeatScore);
	CurFeatScore = min(1000,CurFeatScore);
	if(MaxFeatScore == 0)			// if no scores for any feature then assume all are to be equally sampled
		pChromSeq->RelDensity = 1000;
	else
		pChromSeq->RelDensity = CurFeatScore == 0 ? 0 : 1 + ((1000 * CurFeatScore) / MaxFeatScore);
	}

// now know the feature sequence lengths, alloc to hold all sequences
#ifdef _WIN32
m_pGenomeSeq = (UINT8 *) malloc((size_t)m_GenomeLen);
if(m_pGenomeSeq == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadTranscriptome: Memory allocation of %lld bytes failed",(INT64)m_GenomeLen);
	Reset(false);
	return(eBSFerrMem);
	}
#else
// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
m_pGenomeSeq = (UINT8 *)mmap(NULL,(size_t)m_GenomeLen, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
if(m_pGenomeSeq == MAP_FAILED)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadAssembly: Memory allocation of %lld bytes through mmap()  failed",(INT64)m_GenomeLen,strerror(errno));
	m_pGenomeSeq = NULL;
	Reset(false);
	return(eBSFerrMem);
	}
#endif

int CurRegionLen;
int NumExons;
int StartLoci;
int EndLoci;
int Idx;
int ExonLen;
int NumWarns = 0;
int CurFeatureID = 0;
int PrevFeatureChromID = 0;
TotLen = 0;
pSeq = m_pGenomeSeq;

char szGenomeChromName[128];
int FeatureChromID;
int GenomeChromID;

while((CurFeatureID = m_pBEDFile->GetNextFeatureID(CurFeatureID)) > 0) {
	FeatureChromID = m_pBEDFile->GetFeatureChromID(CurFeatureID);
	if(FeatureChromID != PrevFeatureChromID)
		{
		m_pBEDFile->GetChromosome(FeatureChromID,szGenomeChromName);
		if((GenomeChromID = m_pBioSeqFile->LocateEntryIDbyName(szGenomeChromName)) <= 0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadTranscriptome: Unable to load locate feature '%s' chrom '%s'",pChromSeq->szChromName,szGenomeChromName);
			Reset(false);
			return(eBSFerrChrom);
			}
		PrevFeatureChromID = FeatureChromID;
		}

	pChromSeq = &m_pChromSeqs[CurFeatureID-1];
	pChromSeq->SeqOfs = (UINT32)(UINT64)(pSeq - m_pGenomeSeq);
	TotLen += pChromSeq->Len;
	CurRegionLen = 0;
	if(pChromSeq->Len > 0)
		{
		if(bTranscribed)
			NumExons = m_pBEDFile->GetNumExons(CurFeatureID);
		else
			NumExons = 1;
		for(Idx = 1; Idx <= NumExons; Idx++)
			{
			StartLoci = m_pBEDFile->GetExonStart(CurFeatureID,Idx);
			EndLoci   = m_pBEDFile->GetExonEnd(CurFeatureID,Idx);
			ExonLen = 1 + EndLoci - StartLoci;
			CurRegionLen += ExonLen;
			if((Rslt=m_pBioSeqFile->GetData(GenomeChromID,eSeqBaseType,StartLoci,pSeq,ExonLen))<=0)
				{
				if(NumWarns++ < 10)
					gDiagnostics.DiagOut(eDLWarn,gszProcName,"LoadTranscriptome: Unable to load sequence for feature '%s' from chrom '%s'",pChromSeq->szChromName,szGenomeChromName);
				else
					{
					Reset(false);
					return(eBSFerrMem);
					}
				}
			pSeq += ExonLen;
			}
		}

	if(pChromSeq->Len != CurRegionLen)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadTranscriptome: Expected sum of exon lengths (%d) to equal transcript length (%d) for element/gene '%s'",CurRegionLen,pChromSeq->Len,pChromSeq->szChromName);
		Reset(false);
		return(eBSFerrInternal);
		}
	*pSeq++ = eBaseEOS;		// concatenator separator
	}

delete m_pBioSeqFile;
m_pBioSeqFile = NULL;

LenSCF = (double)INT_MAX/(double)TotLen;
m_GenomeScaledLen = TotLen >= (INT64)INT_MAX ? (long)(TotLen * LenSCF) : (long)TotLen;
CurScaledStart = 0;
pChromSeq = &m_pChromSeqs[0];
for(ChromID = 0; ChromID < m_NumChromSeqs; ChromID++,pChromSeq++)
	{
	pChromSeq->ScaledLen = TotLen >= (INT64)INT_MAX ? (int)(pChromSeq->Len * LenSCF) : pChromSeq->Len;
	pChromSeq->ScaledStart = CurScaledStart;
	CurScaledStart += pChromSeq->ScaledLen;
	}
return(eBSFSuccess);
}


// SimulateSNPs
// Simulate SNPs in the genome at the specified rate per million bases
// Reports their loci
int
SimulateSNPs(char *pszOutSNPs,   // output simulated SNP loci to this file
			int SNPrate)		 // required SNPs per million bases
{
int hFile;
int BuffOfs;
char szLineBuff[8196];
char szChromName[128];
int Ofs;
int SNPiD;
int NumChromSNPs;
etSeqBase *pSeq;
etSeqBase PrevBase;
etSeqBase SNPbase;

int ChromID;
tChromID AssembChromID;
int ChromLoci;

tsChromSeq *pChromSeq;

if(SNPrate == 0)
	return(eBSFSuccess);

if(pszOutSNPs != NULL && pszOutSNPs[0] != '\0')
	{
#ifdef _WIN32
	hFile = open(pszOutSNPs,O_CREATETRUNC );
#else
	if((hFile = open(pszOutSNPs,O_RDWR | O_CREAT,S_IREAD | S_IWRITE))!=-1)
		if(ftruncate(hFile,0)!=0)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate %s - %s",pszOutSNPs,strerror(errno));
				return(eBSFerrCreateFile);
				}
#endif

	if(hFile < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate output file '%s'",pszOutSNPs);
		Reset(false);
		return(eBSFerrCreateFile);
		}
	BuffOfs = sprintf(szLineBuff,"track type=bed name=\"SimSNPs\" description=\"Simulated SNPS\"\n");
	CUtility::SafeWrite(hFile,szLineBuff,BuffOfs);
	}
else
	hFile = -1;

TRandomCombined<TRanrotWGenerator,TRandomMersenne> RG((int)time(0));
BuffOfs = 0;
SNPiD = 0;
for(ChromID = 0; ChromID < m_NumChromSeqs; ChromID++,pChromSeq++)
	{
	pChromSeq = &m_pChromSeqs[ChromID];
	pSeq = &m_pGenomeSeq[pChromSeq->SeqOfs];
	NumChromSNPs = 1 + (int)(pChromSeq->Len * ((double)SNPrate/1000000));
	while(NumChromSNPs)
		{
		Ofs = (int)RG.IRandom(0,pChromSeq->Len-1);
		pSeq = &m_pGenomeSeq[pChromSeq->SeqOfs+Ofs];
		if((PrevBase = *pSeq) > eBaseT)
			continue;
		if((SNPbase = (int)RG.IRandom(0,3))==PrevBase)
			continue;
		*pSeq = SNPbase;
		NumChromSNPs -= 1;

		if(hFile != -1)
			{
			// if transcriptome read generation then remap the transcript relative SNP loci back to the assembly chrom loci
			if(m_pBEDFile != NULL)
				{
				m_pBEDFile->MapTransOfs2Loci(pChromSeq->ChromID,Ofs,NULL,&AssembChromID,&ChromLoci);
				m_pBEDFile->GetChromosome(AssembChromID,szChromName);
				}
			else
				{
				ChromLoci = Ofs;
				strcpy(szChromName,pChromSeq->szChromName);
				}

			SNPiD += 1;
			BuffOfs+=sprintf(&szLineBuff[BuffOfs],"%s\t%d\t%d\tSNP_%d\t%d\t+\n",
					szChromName,ChromLoci,ChromLoci+1,SNPiD,0);

			if(BuffOfs + 200 > sizeof(szLineBuff))
				{
				CUtility::SafeWrite(hFile,szLineBuff,BuffOfs);
				BuffOfs = 0;
				}
			}
		}
	}


if(hFile != -1 && BuffOfs > 0)
	{
	CUtility::SafeWrite(hFile,szLineBuff,BuffOfs);
	}
if(hFile != -1)
	close(hFile);
return(eBSFSuccess);
}

// CmpSeqs
// Returns the number of differences between two sequences of length Len
int
CmpSeqs(int Len,etSeqBase *pSeq1,etSeqBase *pSeq2)
{
int Diffs = 0;
etSeqBase Base1;
etSeqBase Base2;
int Idx;
for(Idx=0; Idx < Len; Idx++)
	if((Base1 = (*pSeq1++ & NUCONLYMSK)) != (Base2 = (*pSeq2++ & NUCONLYMSK)))
		Diffs += 1;
return(Diffs);
}

bool
IsDupSeq(int Len,etSeqBase *pSeq1,etSeqBase *pSeq2)
{
return(CmpSeqs(Len,pSeq1,pSeq2) == 0 ? true : false);
}



// both the forward and reverse complements are compared for equivilence
bool
IsDupSeqStrands(int Len,etSeqBase *pReadSeq1,etSeqBase *pReadSeq2)
{
etSeqBase *pSeq1 = pReadSeq1;
etSeqBase *pSeq2 = pReadSeq2;
etSeqBase Base1;
etSeqBase Base2;
int Idx;

// first assume both reads will ultimately align onto the forward (watson) strand
pSeq1 = pReadSeq1;
pSeq2 = pReadSeq2;
for(Idx=0; Idx < Len; Idx++)
	if((Base1 = (*pSeq1++ & NUCONLYMSK)) != (Base2 = (*pSeq2++ & NUCONLYMSK)))
		break;
if(Idx == Len)
	return(true);

// now assume that pReadSeq1 stays watson but that pReadSeq2 will be reverse complemented to align onto the crick strand
pSeq1 = pReadSeq1;
pSeq2 = &pReadSeq2[Len-1];
for(Idx=0; Idx < Len; Idx++)
	{
	Base1 = (*pSeq1++ & NUCONLYMSK);
	Base2 = (*pSeq2-- & NUCONLYMSK);
	switch(Base2) {
		case eBaseA: Base2=eBaseT; break;
		case eBaseC: Base2=eBaseG; break;
		case eBaseG: Base2=eBaseC; break;
		case eBaseT: Base2=eBaseA; break;
		}
	if(Base1 != Base2)
		break;
	}
if(Idx == Len)
	return(true);

// finally assume that pReadSeq2 stays watson but that pReadSeq1 will be reverse complemented to align onto the crick strand
pSeq1 = &pReadSeq1[Len-1];
pSeq2 = pReadSeq2;
for(Idx=0; Idx < Len; Idx++)
	{
	Base1 = (*pSeq1-- & NUCONLYMSK);
	Base2 = (*pSeq2++ & NUCONLYMSK);
	switch(Base1) {
		case eBaseA: Base1=eBaseT; break;
		case eBaseC: Base1=eBaseG; break;
		case eBaseG: Base1=eBaseC; break;
		case eBaseT: Base1=eBaseA; break;
		}
	if(Base1 != Base2)
		return(false);
	}
return(true);
}

//
void
ShowDupSeqs(int ReadLen,tsSimRead *pRead1,tsSimRead *pRead2)
{
int Idx;
etSeqBase *pSeq;
UINT8 Seq[2000];
char szSeq[2000];
pSeq = pRead1->pSeq;
for(Idx = 0; Idx < ReadLen; Idx++,pSeq++)
	Seq[Idx] = *pSeq & NUCONLYMSK;
CSeqTrans::MapSeq2Ascii(Seq,ReadLen,szSeq);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"chrom: %1.2d, loci: %1.9d HD: %1.2d seq: '%s'",
			 pRead1->ChromID,pRead1->StartLoci,pRead1->HammingDist,szSeq);
pSeq = pRead2->pSeq;
for(Idx = 0; Idx < ReadLen; Idx++,pSeq++)
	Seq[Idx] = *pSeq & NUCONLYMSK;
CSeqTrans::MapSeq2Ascii(Seq,ReadLen,szSeq);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"chrom: %1.2d, loci: %1.9d HD: %1.2d seq: '%s'",
			 pRead2->ChromID,pRead2->StartLoci,pRead2->HammingDist,szSeq);
}


// TransformToColorspace
// transformation of basespace and colorspace (SOLiD)
//
UINT8 CvrtSOLiDmap[5][7] = {
	{0,1,2,3,4,(UINT8)eBaseEOS,(UINT8)eBaseEOG},	// PrevBase==etBaseA
	{1,0,3,2,4,(UINT8)eBaseEOS,(UINT8)eBaseEOG},	// PrevBase==etBaseC
	{2,3,0,1,4,(UINT8)eBaseEOS,(UINT8)eBaseEOG},	// PrevBase==etBaseG
	{3,2,1,0,4,(UINT8)eBaseEOS,(UINT8)eBaseEOG},	// PrevBase==etBaseT
	{0,1,2,3,4,(UINT8)eBaseEOS,(UINT8)eBaseEOG},	// PrevBase==etBaseN
	};

int
TransformToColorspace(UINT8 *pSrcBases,		// basespace sequence to transform
					UINT32 SeqLen,			// sequence length
					char *pDstColorspace,	// where to write transformed
					etSeqBase Seed)			// use this base as the seed base immediately preceding the first base
{
etSeqBase PrvBase;
etSeqBase CurBase;
UINT8 ColorSpace;
if(pSrcBases == NULL || SeqLen == 0 || pDstColorspace == NULL || Seed > eBaseN)
	return(eBSFerrParams);
PrvBase = Seed;
while(SeqLen--)
	{
	CurBase = *pSrcBases++ & 0x0f;						// conversion is on the low nibble
	if(!(CurBase == eBaseEOS || CurBase == eBaseEOG) && PrvBase <= 4)
		{
		CurBase &= ~cRptMskFlg;
		if(CurBase >= eBaseN)
			CurBase = eBaseN;
		ColorSpace = CvrtSOLiDmap[PrvBase][CurBase];
		*pDstColorspace++ = 0x30 + ColorSpace;
		}
	else
		*pDstColorspace++ = 0x30 + CurBase;
	PrvBase = CurBase;
	}
*pDstColorspace = '\0';
return(eBSFSuccess);
}

// TransformToColorspaceBase
// BWA requires that the colorspace be represented as double encoded basespace - UGH!!!
int
TransformToColorspaceBase(UINT8 *pSrcBases,		// basespace sequence to transform
					UINT32 SeqLen,			// sequence length
					char *pDstColorspace,	// where to write transformed
					etSeqBase Seed)			// use this base as the seed base immediately preceding the first base
{
etSeqBase PrvBase;
etSeqBase CurBase;
UINT8 ColorSpace;
UINT8 ColorSpaceBase;
if(pSrcBases == NULL || SeqLen == 0 || pDstColorspace == NULL || Seed > eBaseN)
	return(eBSFerrParams);
PrvBase = Seed;
while(SeqLen--)
	{
	CurBase = *pSrcBases++ & 0x0f;						// conversion is on the low nibble
	if(!(CurBase == eBaseEOS || CurBase == eBaseEOG) && PrvBase <= 4)
		{
		CurBase &= ~cRptMskFlg;				// remove any repeat masking
		if(CurBase > eBaseN)
			CurBase = eBaseN;

		ColorSpace = CvrtSOLiDmap[PrvBase][CurBase];
		switch(ColorSpace) {
			case 0: ColorSpaceBase = 'A'; break;
			case 1: ColorSpaceBase = 'C'; break;
			case 2: ColorSpaceBase = 'G'; break;
			case 3: ColorSpaceBase = 'T'; break;
			default: ColorSpaceBase = 'N'; break;
			}

		*pDstColorspace++ = ColorSpaceBase;
		}
	else
		*pDstColorspace++ = 'N';
	PrvBase = CurBase;
	}
*pDstColorspace = '\0';
return(eBSFSuccess);
}

//
typedef struct TAG_sReadRprt {
	etSeqBase FwdSeq[cMaxReadLen+1];
	etSeqBase RevSeq[cMaxReadLen+1];
	char szColorspace[cMaxReadLen+2];
	char szLineBuff[0x03fff];
	int LineLen;
	int NumCols;
} sReadRprt;

int
ReportReads(bool bPEgen,			// true if paired end simulated reads being simulated
		    int Region,			// Process regions 0:ALL,1:CDS,2:5'UTR,3:3'UTR,4:Introns,5:5'US,6:3'DS,7:Intergenic (default = ALL)
			int NomReadLen,			// reads nominal length
			int NumPrevReported,	// number of reads previously reported
			int MaxReads,			// maximum number of reads remaining to be reported
			bool bDedupe,			// true if reads are to be deduped
			bool bReadHamDist,	    // true if hamming distributions from each sampled read to all other genome subsequences to be generated
			int AvailReads)			// number of reads to be reported on from m_pSimReads[]
{
tsSimRead *pRead;
tsSimRead *pPairRead;
tsSimRead *pPrevRead;
tsChromSeq *pChromSeq;
etSeqBase *pReadSeq;
char *pszIsRand;
char *pQualChr;
int ReadOfs;
int ReadsIdx;
int ReadLenRem;

int NumOfDups;

etSeqBase CSseedBase;

sReadRprt PEreads[2];			// to hold each paired end read


int NumReported;

// check if any reads to report
if(m_pSimReads == NULL || MaxReads < 1 || AvailReads < 1)
	return(0);

if(bDedupe && AvailReads > 1)
	{
	NumOfDups = 0;
	// for detection of the duplicate read sequences, the sequences need to be sorted
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Deduping %d read sequences, sorting...",AvailReads);
	qsort(m_pSimReads,AvailReads,sizeof(tsSimRead),SortSimReads);
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Sorting completed, checking for duplicates");
	// firstly find and mark all watson '+' strand duplicates
	pPrevRead = m_pSimReads;
	pRead = &m_pSimReads[1];
	for(ReadsIdx = 1; ReadsIdx < AvailReads; ReadsIdx++, pRead++)
		{
		if(pPrevRead->Len == pRead->Len && IsDupSeq(pRead->Len,pPrevRead->pSeq,pRead->pSeq))
			{
			if(pRead->Status > 1)
				pPrevRead->Status = 2;
			pRead->Status = 2;
			if(pRead->pPartner != NULL)
				{
				pRead->pPartner->Status = 2;
				}
			NumOfDups+=1;
			}
		else
			{
			if(pPrevRead->Len == pRead->Len)
				{
				int Hamming = CmpSeqs(pRead->Len,pPrevRead->pSeq,pRead->pSeq);
				if(Hamming < pPrevRead->HammingDist || Hamming < pRead->HammingDist)
					{
					if(pRead->Status > 1)
						pPrevRead->Status = 2;
					pRead->Status = 2;		// currently treat these as if duplicates
					if(pRead->pPartner != NULL)
						pRead->pPartner->Status = 2;
					NumOfDups+=1;
					}
				}
			pPrevRead = pRead;
			}
		}
	// next, for those sequences not already marked as duplicates, see if after revcpl into the crick they match one still on watson
	// iterate and revcpl each non-dup sequence
	// search for this crick sequence
	// if none found this is a truely unique sequence
	// if found then mark all instances as being duplicates
	int RevIdx;
	pRead = m_pSimReads;
	for(ReadsIdx = 0; ReadsIdx < AvailReads; ReadsIdx++, pRead++)
		{
		if((RevIdx=LocateRevCplSeq(pRead->Len,pRead->pSeq,AvailReads)) >= 0)
			{
			ShowDupSeqs(pRead->Len,&m_pSimReads[RevIdx],pRead);
			pRead->Status = 2;
			if(pRead->pPartner != NULL)
				pRead->pPartner->Status = 2;
			NumOfDups+=1;
			}
		}
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Marked %d reads out of %d as duplicate sequences",NumOfDups,AvailReads);
	}

// if transcriptome read generation then remap the transcript relative loci back to the assembly chrom loci
if(Region == eMEGRAny && m_pBEDFile != NULL)
	{
	tChromID ChromID;
	int ChromLoci;

	pRead = m_pSimReads;
	NumReported = 0;
	for(ReadsIdx = 0; ReadsIdx < AvailReads; ReadsIdx++, pRead++)
		{
		if(pRead->Status > 0)	// slough if read was determined as being a duplicate
			continue;
		pChromSeq = &m_pChromSeqs[pRead->ChromSeqID-1];
		m_pBEDFile->MapTransOfs2Loci(pRead->ChromID,pRead->StartLoci,NULL,&ChromID,&ChromLoci);
		m_pBEDFile->GetChromosome(ChromID,pChromSeq->szChromName);
		pRead->StartLoci = ChromLoci;
		m_pBEDFile->MapTransOfs2Loci(pRead->ChromID,pRead->EndLoci,NULL,NULL,&ChromLoci);
		pRead->ChromID = (int)ChromID;
		pRead->EndLoci = ChromLoci;
		}
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Writing simulated reads to file: '%s'",m_pszOutFile);


memset(PEreads,0,sizeof(PEreads));
pRead = m_pSimReads;
NumReported = 0;
if(m_PMode == ePMSampHamm)
	PEreads[0].LineLen = sprintf(PEreads[0].szLineBuff,"\"Chrom\",\"Loci\",\"Hamming\"\n");

for(ReadsIdx = 0; ReadsIdx < AvailReads; ReadsIdx++, pRead++)
	{
	if(pRead->Status > 0)	// slough if read already reported on, or if a duplicate
		continue;
	int Idx;
	int Pair1ReadLen;
	int Pair2ReadLen;

	pReadSeq = pRead->pSeq;
	Pair1ReadLen = pRead->Len;
	for(Idx = 0; Idx < (Pair1ReadLen + m_InDelSize); Idx++,pReadSeq++)	// extra read length is in case InDels being simulated and upto 10 bases deleted
		{
		if((*pReadSeq & NUCONLYMSK) > eBaseN)
			break;
		PEreads[0].FwdSeq[Idx] = *pReadSeq & NUCONLYMSK;
		}
	if(Idx != (Pair1ReadLen + m_InDelSize))
		{
		pRead->Status = 2;
		if(bPEgen == true)
			pRead->pPartner->Status = 2;
		continue;
		}

	if(bPEgen == true)
		{
		pPairRead =  pRead->pPartner;
		Pair2ReadLen = pPairRead->Len;
		pReadSeq = pPairRead->pSeq;
		for(Idx = 0; Idx < (Pair2ReadLen + m_InDelSize); Idx++,pReadSeq++)	// extra read length is in case InDels being simulated and upto 10 bases deleted
			{
			if((*pReadSeq & NUCONLYMSK) > eBaseN)
				break;
			PEreads[1].FwdSeq[Idx] = *pReadSeq & NUCONLYMSK;
			}
		if(Idx != (Pair2ReadLen + m_InDelSize))
			{
			pPairRead->Status = 2;
			pRead->Status = 2;
			continue;
			}
		}

	pReadSeq = PEreads[0].FwdSeq;

	NumReported += 1;
	pChromSeq = &m_pChromSeqs[pRead->ChromSeqID-1];
	switch(m_FMode) {
		case eFMcsv:
		case eFMcsvSeq:
			if(m_PMode == ePMSampHamm)
				{
				PEreads[0].LineLen += sprintf(&PEreads[0].szLineBuff[PEreads[0].LineLen],"\"%s\",%d,%d\n",
					pChromSeq->szChromName,pRead->StartLoci,pRead->HammingDist);
				break;
				}

			PEreads[0].LineLen += sprintf(&PEreads[0].szLineBuff[PEreads[0].LineLen],"%d,\"%s\",\"%s\",\"%s\",%d,%d,%d,\"%c\",%d",
					NumPrevReported+NumReported,"usimreads",m_szSpecies,pChromSeq->szChromName,pRead->StartLoci,pRead->EndLoci,pRead->Len,pRead->Strand ? '-' : '+',pRead->HammingDist);
			if(m_FMode == eFMcsvSeq)
				{
				PEreads[0].LineLen += sprintf(&PEreads[0].szLineBuff[PEreads[0].LineLen],",\"");
				pReadSeq = PEreads[0].FwdSeq;

				if(pRead->Strand)
					{
					memcpy(PEreads[0].RevSeq,pReadSeq,Pair1ReadLen);
					CSeqTrans::ReverseComplement(Pair1ReadLen,PEreads[0].RevSeq);
					pReadSeq = PEreads[0].RevSeq;
					}
				CSeqTrans::MapSeq2Ascii(pReadSeq,Pair1ReadLen,&PEreads[0].szLineBuff[PEreads[0].LineLen]);
				PEreads[0].LineLen += Pair1ReadLen;
				PEreads[0].LineLen += sprintf(&PEreads[0].szLineBuff[PEreads[0].LineLen],"\"%s\"",pReadSeq);
				}
			if(bReadHamDist)
				{
				int HamDistIdx;
				UINT32 *pCnt = pRead->pHamDistFreq;
				for(HamDistIdx=0; HamDistIdx <= pRead->Len; HamDistIdx++,pCnt++)
					PEreads[0].LineLen += sprintf(&PEreads[0].szLineBuff[PEreads[0].LineLen],",%d",*pCnt);
				}
			PEreads[0].LineLen += sprintf(&PEreads[0].szLineBuff[PEreads[0].LineLen],"\n");
			break;

		case eFMFasta:
		case eFMNWFasta:
			{
			int NumSubs;
			int InDelSize;

			if(m_FMode == eFMFasta)
				m_MaxFastaLineLen = 79;
			else
				m_MaxFastaLineLen = cMaxSimReadLen;

			pReadSeq = PEreads[0].FwdSeq;
			if(pRead->Strand)
				{
				memcpy(PEreads[0].RevSeq,pReadSeq,Pair1ReadLen+m_InDelSize);
				CSeqTrans::ReverseComplement(Pair1ReadLen+m_InDelSize,PEreads[0].RevSeq);
				pReadSeq = PEreads[0].RevSeq;
				}

			SimArtefacts(false,Pair1ReadLen,pReadSeq);
			SimArtefacts(true,Pair1ReadLen,pReadSeq);

			if(!m_PropRandReads || (m_PropRandReads <= RGseeds.IRandom(0,1000000)))
				{
				InDelSize = SimInDels(pRead,&Pair1ReadLen,pReadSeq);
				NumSubs = SimSeqErrors(Pair1ReadLen,pReadSeq);
				pszIsRand = (char *)"lcl";
				}
			else
				{
				InDelSize = 0;
				NumSubs = SimSeqRand(Pair1ReadLen,pReadSeq);
				pszIsRand = (char *)"lcr";
				}

			PEreads[0].LineLen+=sprintf(&PEreads[0].szLineBuff[PEreads[0].LineLen],">%s|usimreads|%1.8d|%s|%d|%d|%d|%c|%d|%d|%d\n",pszIsRand,
					NumPrevReported+NumReported,pChromSeq->szChromName,pRead->StartLoci,pRead->EndLoci+InDelSize,Pair1ReadLen,pRead->Strand ? '-' : '+',pRead->HammingDist,NumSubs,InDelSize);

			ReadOfs = 0;
			ReadLenRem = Pair1ReadLen;
			while(ReadLenRem)
				{
				PEreads[0].NumCols = ReadLenRem > m_MaxFastaLineLen ? m_MaxFastaLineLen : ReadLenRem;
				CSeqTrans::MapSeq2UCAscii(&pReadSeq[ReadOfs],PEreads[0].NumCols,&PEreads[0].szLineBuff[PEreads[0].LineLen]);
				PEreads[0].LineLen += PEreads[0].NumCols;
				PEreads[0].LineLen += sprintf(&PEreads[0].szLineBuff[PEreads[0].LineLen],"\n");
				ReadLenRem -= PEreads[0].NumCols;
				ReadOfs += PEreads[0].NumCols;
				}

			if(bPEgen == true)
				{
				NumReported += 1;
				pReadSeq = PEreads[1].FwdSeq;
				if(pPairRead->Strand)
					{
					memcpy(PEreads[1].RevSeq,pReadSeq,Pair2ReadLen+m_InDelSize);
					CSeqTrans::ReverseComplement(Pair2ReadLen+m_InDelSize,PEreads[1].RevSeq);
					pReadSeq = PEreads[1].RevSeq;
					}

				SimArtefacts(false,Pair2ReadLen,pReadSeq);
				SimArtefacts(true,Pair2ReadLen,pReadSeq);

				if(!m_PropRandReads || (m_PropRandReads <= RGseeds.IRandom(0,1000000)))
					{
					InDelSize = SimInDels(pPairRead,&Pair2ReadLen,pReadSeq);
					NumSubs = SimSeqErrors(Pair2ReadLen,pReadSeq);
					pszIsRand = (char *)"lcl";
					}
				else
					{
					InDelSize = 0;
					NumSubs = SimSeqRand(Pair2ReadLen,pReadSeq);
					pszIsRand = (char *)"lcr";
					}

				PEreads[1].LineLen+=sprintf(&PEreads[1].szLineBuff[PEreads[1].LineLen],">%s|usimreads|%1.8d|%s|%d|%d|%d|%c|%d|%d|%d\n",pszIsRand,
						NumPrevReported+NumReported,pChromSeq->szChromName,pPairRead->StartLoci,pPairRead->EndLoci+InDelSize,Pair2ReadLen,pPairRead->Strand ? '-' : '+',pPairRead->HammingDist,NumSubs,InDelSize);

				ReadOfs = 0;
				ReadLenRem = Pair2ReadLen;
				while(ReadLenRem)
					{
					PEreads[1].NumCols = ReadLenRem > m_MaxFastaLineLen ? m_MaxFastaLineLen : ReadLenRem;
					CSeqTrans::MapSeq2UCAscii(&pReadSeq[ReadOfs],PEreads[1].NumCols,&PEreads[1].szLineBuff[PEreads[1].LineLen]);
					PEreads[1].LineLen += PEreads[1].NumCols;
					PEreads[1].LineLen += sprintf(&PEreads[1].szLineBuff[PEreads[1].LineLen],"\n");
					ReadLenRem -= PEreads[1].NumCols;
					ReadOfs += PEreads[1].NumCols;
					}
				}
			break;
			}

		case eFMSOLiD:
		case eFMSOLiDbwa:
			{
			int NumSubs;
			pReadSeq = PEreads[0].FwdSeq;
			if(pRead->Strand)
				{
				memcpy(PEreads[0].RevSeq,pReadSeq,Pair1ReadLen);
				if(m_FMode == eFMSOLiD)
					CSeqTrans::ReverseComplement(Pair1ReadLen,PEreads[0].RevSeq);
				else
					CSeqTrans::ReverseSeq(Pair1ReadLen,PEreads[0].RevSeq);
				pReadSeq = PEreads[0].RevSeq;
				}

			int InDelSize;

			if(!m_PropRandReads || (m_PropRandReads < RGseeds.IRandom(0,1000000)))
				{
				InDelSize = SimInDels(pRead,&Pair1ReadLen,pReadSeq);
				NumSubs = SimSeqErrors(Pair1ReadLen,pReadSeq);
				pszIsRand = (char *)"lcl";
				}
			else
				{
				InDelSize = 0;
				NumSubs = SimSeqRand(Pair1ReadLen,pReadSeq);
				pszIsRand = (char *)"lcr";
				}

			CSseedBase = (UINT8)RGseeds.IRandom(0,3);
			PEreads[0].szColorspace[0] = CSeqTrans::MapBase2Ascii(CSseedBase);
			if(m_FMode == eFMSOLiD)
				TransformToColorspace(pReadSeq,Pair1ReadLen,&PEreads[0].szColorspace[1],CSseedBase);
			else
				TransformToColorspaceBase(pReadSeq,Pair1ReadLen,&PEreads[0].szColorspace[1],CSseedBase);

			if(m_FMode == eFMSOLiD)
				{
				PEreads[0].LineLen+=sprintf(&PEreads[0].szLineBuff[PEreads[0].LineLen],">%s|usimreads|%1.8d|%s|%d|%d|%d|%c|%d|%d|%d\n",pszIsRand,
					NumPrevReported+NumReported,pChromSeq->szChromName,pRead->StartLoci,pRead->EndLoci+InDelSize,Pair1ReadLen,pRead->Strand ? '-' : '+',pRead->HammingDist,NumSubs,InDelSize);
				PEreads[0].LineLen+=sprintf(&PEreads[0].szLineBuff[PEreads[0].LineLen],"%s\n",PEreads[0].szColorspace);
				}
			else
				{
				PEreads[0].LineLen+=sprintf(&PEreads[0].szLineBuff[PEreads[0].LineLen],"@%s|usimreads|%1.8d|%s|%d|%d|%d|%c|%d|%d|%d\n",pszIsRand,
					NumPrevReported+NumReported,pChromSeq->szChromName,pRead->StartLoci,pRead->EndLoci+InDelSize,Pair1ReadLen,pRead->Strand ? '-' : '+',pRead->HammingDist,NumSubs,InDelSize);
				PEreads[0].LineLen+=sprintf(&PEreads[0].szLineBuff[PEreads[0].LineLen],"%s\n",PEreads[0].szColorspace);
				PEreads[0].LineLen+=sprintf(&PEreads[0].szLineBuff[PEreads[0].LineLen],"+ descriptor line\n");
				pQualChr = &PEreads[0].szLineBuff[PEreads[0].LineLen];
				*pQualChr++ = '!';			// dummy quality scores
				for(Idx = 0; Idx < Pair1ReadLen; Idx++)
					*pQualChr++ = (char)(((pRead->StartLoci + Idx) % 20) + 65);
				*pQualChr++ = '\n';
				*pQualChr = '\0';
				PEreads[0].LineLen += Pair1ReadLen + 2;
				}

			if(bPEgen == true)
				{
				NumReported += 1;
				pReadSeq = PEreads[1].FwdSeq;
				if(pPairRead->Strand)
					{
					memcpy(PEreads[1].RevSeq,pReadSeq,Pair2ReadLen);
					if(m_FMode == eFMSOLiD)
						CSeqTrans::ReverseComplement(Pair2ReadLen,PEreads[1].RevSeq);
					else
						CSeqTrans::ReverseSeq(Pair2ReadLen,PEreads[1].RevSeq);
					pReadSeq = PEreads[1].RevSeq;
					}

				if(!m_PropRandReads || (m_PropRandReads < RGseeds.IRandom(0,1000000)))
					{
					InDelSize = SimInDels(pPairRead,&Pair2ReadLen,pReadSeq);
					NumSubs = SimSeqErrors(Pair2ReadLen,pReadSeq);
					pszIsRand = (char *)"lcl";
					}
				else
					{
					InDelSize = 0;
					NumSubs = SimSeqRand(Pair2ReadLen,pReadSeq);
					pszIsRand = (char *)"lcr";
					}

				CSseedBase = (UINT8)RGseeds.IRandom(0,3);
				PEreads[1].szColorspace[0] = CSeqTrans::MapBase2Ascii(CSseedBase);

				if(m_FMode == eFMSOLiD)
					TransformToColorspace(pReadSeq,Pair2ReadLen,&PEreads[1].szColorspace[1],CSseedBase);
				else
					TransformToColorspaceBase(pReadSeq,Pair2ReadLen,&PEreads[1].szColorspace[1],CSseedBase);

				if(m_FMode == eFMSOLiD)
					{
					PEreads[1].LineLen+=sprintf(&PEreads[1].szLineBuff[PEreads[1].LineLen],">%s|usimreads|%1.8d|%s|%d|%d|%d|%c|%d|%d|%d\n",pszIsRand,
						NumPrevReported+NumReported,pChromSeq->szChromName,pPairRead->StartLoci,pPairRead->EndLoci+InDelSize,Pair2ReadLen,pPairRead->Strand ? '-' : '+',pPairRead->HammingDist,NumSubs,InDelSize);
					PEreads[1].LineLen+=sprintf(&PEreads[1].szLineBuff[PEreads[1].LineLen],"%s\n",PEreads[1].szColorspace);
					}
				else
					{
					PEreads[1].LineLen+=sprintf(&PEreads[1].szLineBuff[PEreads[1].LineLen],">%s|usimreads|%1.8d|%s|%d|%d|%d|%c|%d|%d|%d\n",pszIsRand,
					NumPrevReported+NumReported,pChromSeq->szChromName,pPairRead->StartLoci,pPairRead->EndLoci+InDelSize,Pair2ReadLen,pPairRead->Strand ? '-' : '+',pPairRead->HammingDist,NumSubs,InDelSize);
					PEreads[1].LineLen+=sprintf(&PEreads[1].szLineBuff[PEreads[1].LineLen],"%s\n",PEreads[1].szColorspace);
					PEreads[1].LineLen+=sprintf(&PEreads[1].szLineBuff[PEreads[1].LineLen],"+ descriptor line\n");
					pQualChr = &PEreads[1].szLineBuff[PEreads[1].LineLen];
					*pQualChr++ = '!';
					for(Idx = 0; Idx < Pair2ReadLen; Idx++)
						*pQualChr++ = (char)(((pPairRead->StartLoci + Idx) % 20) + 65);
					*pQualChr++ = '\n';
					pQualChr = '\0';
					PEreads[1].LineLen += Pair2ReadLen + 2;
					}
				}
			break;
			}
		}

	if((PEreads[0].LineLen + 4000) > sizeof(PEreads[0].szLineBuff))
		{
		if(write(m_hOutFile,PEreads[0].szLineBuff,PEreads[0].LineLen) != PEreads[0].LineLen)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Errors whilst writing %d bytes - '%s' - %s",PEreads[0].LineLen, m_pszOutFile, strerror(errno));
			Reset(false);
			return(eBSFerrFileAccess);
			}
		PEreads[0].LineLen=0;
		}
	if(bPEgen == true)
		{
		if((PEreads[1].LineLen + 4000) > sizeof(PEreads[1].szLineBuff))
			{
			if(write(m_hOutPEFile,PEreads[1].szLineBuff,PEreads[1].LineLen) != PEreads[1].LineLen)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Errors whilst writing %d bytes - '%s' - %s",PEreads[1].LineLen, m_pszOutPEFile, strerror(errno));
				Reset(false);
				return(eBSFerrFileAccess);
				}
			PEreads[1].LineLen=0;
			}
		}

	pRead->Status = 1;						// mark this read as having being reported
	if(bPEgen)
		pPairRead->Status = 1;

	if(NumPrevReported + NumReported == MaxReads)
		break;
	}

if(PEreads[0].LineLen && write(m_hOutFile,PEreads[0].szLineBuff,PEreads[0].LineLen) != PEreads[0].LineLen)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Errors whilst writing %d bytes - '%s' - %s",PEreads[0].LineLen, m_pszOutFile, strerror(errno));
	Reset(false);
	return(eBSFerrFileAccess);
	}
#ifdef _WIN32
_commit(m_hOutFile);
#else
fsync(m_hOutFile);
#endif

if(bPEgen == true)
	{
	if(PEreads[1].LineLen && write(m_hOutPEFile,PEreads[1].szLineBuff,PEreads[1].LineLen) != PEreads[1].LineLen)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Errors whilst writing %d bytes - '%s' - %s",PEreads[1].LineLen, m_pszOutPEFile, strerror(errno));
		Reset(false);
		return(eBSFerrFileAccess);
		}
	#ifdef _WIN32
	_commit(m_hOutPEFile);
	#else
	fsync(m_hOutPEFile);
	#endif
	}

return(NumReported);
}

int
Process(etPMode PMode,		// processing mode
		etSEMode SEMode,	// induced sequencer error rate mode
		bool bPEgen,		// true if paired ends are to be generated
		int PEmin,			// PE minimum fragment length
		int PEmax,			// PE maximum fragment length
		double PropRandReads, // generate completely random reads at this rate
		int DistCluster,	// cluster generated reads into windows of this median length, 0 if no clustering
		double SeqErrRate,	// dynamic sequencing errors are to be induced into generated sequences at this rate
		bool bSeqErrProfile,// true if to generate composite sequencer errors with uniform profile (default is Illumina 3' skewed)
		int SNPrate,		// simulate SNPs at this rate per million bases
		int InDelSize,		// simulated InDel size range
		double InDelRate,	// simulated InDel rate per read
		bool bReadHamDist,	// true if hamming distributions from each sampled read to all other genome subsequences to be generated
		etFMode FMode,		// output format
		int NumThreads,		// number of worker threads to use
		char Strand,		// generate for this strand '+' or '-' or for both '*'
		int NumReads,		// number of reads required (will be 2x this number if generating paired ends)
		int ReadLen,		// read lengths
		double Artef5Rate,			// rate (0..1) at which to insert 5' artefact sequences
		int NumArtef5Seqs,			// number of user specified 5' artefact sequences
		char *pszArtef5Seqs[], // 5' artefact sequences
		double Artef3Rate,			// rate (0..1) at which to insert 3' artefact sequences
		int NumArtef3Seqs,			// number of user specified 3' artefact sequences
		char *pszArtef3Seqs[], // 5' artefact sequences
		int CutMin,			// min cut length
		int CutMax,			// max cut length
		bool bDedupe,		// true if unique read sequences only to be generated
		int DfltHamming,	// if < 0 then Hamming distance to be dynamically calculated otherwise default Hamming to this value
		int Region,			// Process regions 0:ALL,1:CDS,2:5'UTR,3:3'UTR,4:Introns,5:5'US,6:3'DS,7:Intergenic (default = ALL)
		int UpDnRegLen,		// if processing regions then up/down stream regulatory length
		char *pszFeatFile,	// optionally generate transcriptome reads from features or genes in this BED file
		char *pszInFile,	// input from this raw multifasta or bioseq assembly
		char *pszProfFile,	// input from this profile site preferences file
		char *pszHammFile,	// use Hamming edit distances from this file
		char *pszOutPEFile, // output partner paired end simulated reads to this file
		char *pszOutFile,	// output simulated reads to this file
		char *pszOutSNPs)   // output simulated SNP loci to this file
{
int Rslt;
bool bUseLocateHamming;
int ReadsOfs;
int NumReadsReq;
int ReadsCnt;
tsWorkerPars WorkerThreads[cMaxWorkerThreads];
tsWorkerPars *pCurThread;
int ReadsPerThread;
int ThreadIdx;
bool bFirst;
int MaxReadsPerBatch;		// process at most this number of simulated reads per batch
int ReportedReads;			// number of reads reported on by last call to ReportReads()
int TotReportedReads;		// total number of reads reported on by all batches processed by ReportReads()
int ExhustedChroms;			// number of threads which exhusted attempts to find chroms from which reads can be simulated
int CurNumGenReads;
int PrevNumGenReads;
int MinChromLen;

Init();
m_PMode = PMode;
m_FMode = FMode;
m_SEMode = SEMode;
m_PropRandReads = (int)(PropRandReads * 1000000.0);
m_DynProbErr = SeqErrRate;
m_bUniformDist = bSeqErrProfile;
m_InDelSize = InDelSize;
m_InDelRate = InDelRate;
m_DistCluster = DistCluster;
m_TotReqReads = bPEgen ? NumReads * 2 : NumReads;

m_Artef5Rate = Artef5Rate;				// rate (0..1) at which to insert 5' artefact sequences
m_NumArtef5Seqs = NumArtef5Seqs;		// number of user specified 5' artefact sequences
m_ppszArtef5Seqs = pszArtef5Seqs;		// 5' artefact sequences
m_Artef3Rate = Artef3Rate;				// rate (0..1) at which to insert 3' artefact sequences
m_NumArtef3Seqs = NumArtef3Seqs;		// number of user specified 3' artefact sequences
m_ppszArtef3Seqs = pszArtef3Seqs;		// 5' artefact sequences
m_Artef5SeqLens[0] = 0;
m_Artef3SeqLens[0] = 0;
for(int Idx = 0; Idx < m_NumArtef5Seqs; Idx++)
	{
	m_Artef5SeqLens[Idx] = (int)strlen(pszArtef5Seqs[Idx]);
	CSeqTrans::MapAscii2Sense(pszArtef5Seqs[Idx],m_Artef5SeqLens[Idx],m_Artef5Seqs[Idx]);
	}
for(int Idx = 0; Idx < m_NumArtef3Seqs; Idx++)
	{
	m_Artef3SeqLens[Idx] = (int)strlen(pszArtef3Seqs[Idx]);
	CSeqTrans::MapAscii2Sense(pszArtef3Seqs[Idx],m_Artef3SeqLens[Idx],m_Artef3Seqs[Idx]);
	}

bUseLocateHamming = (pszHammFile != NULL && pszHammFile[0] != '\0') ? true : false;

if(PMode != ePMStandard)
	{
	if((Rslt=InitProfSitePrefs(pszProfFile))!=eBSFSuccess)
		{
		Reset(false);
		return(Rslt);
		}
	}

// estimate minimum length chroms from which reads can be reliably sampled
if(bPEgen)
	MinChromLen = PEmax;
else
	MinChromLen = CutMax;
MinChromLen += 20;

if(Region != eMEGRAny || pszFeatFile == NULL || pszFeatFile[0] == '\0')
	{
	if((Rslt=LoadGenome(MinChromLen,pszInFile))!=eBSFSuccess)	// need to load complete assembly into memory
		{
		Reset(false);
		return(Rslt);
		}
	}
else
	{
	if((Rslt=LoadTranscriptome(pszInFile,pszFeatFile))!=eBSFSuccess)	// need to load transcriptome into memory
		{
		Reset(false);
		return(Rslt);
		}
	}

// if filtering by region...
if(Region != eMEGRAny)
	{
	// open feature file
	if((m_pBEDFile = new CBEDfile) == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CBEDfile object");
		return(eBSFerrObj);
		}

	if((Rslt = m_pBEDFile->Open(pszFeatFile))!=eBSFSuccess)
		{
		while(m_pBEDFile->NumErrMsgs())
			gDiagnostics.DiagOut(eDLFatal,gszProcName,m_pBEDFile->GetErrMsg());
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open feature/gene BED file '%s'",pszFeatFile);
		Reset(false);
		return(Rslt);
		}

	 if(!m_pBEDFile->ContainsGeneDetail())
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"BED file '%s' does not contain gene feature regions",pszFeatFile);
		Reset(false);
		return(eBSFerrEntry);
		}
	}

if(SNPrate > 0)				// simulate SNPs?
	{
	if((Rslt = SimulateSNPs(pszOutSNPs,SNPrate))!=eBSFSuccess)
		{
		Reset(false);
		return(Rslt);
		}
	}

#ifdef _WIN32
m_hOutFile = open(pszOutFile,O_CREATETRUNC );
#else
if((m_hOutFile = open(pszOutFile,O_RDWR | O_CREAT,S_IREAD | S_IWRITE))!=-1)
    if(ftruncate(m_hOutFile,0)!=0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate %s - %s",pszOutFile,strerror(errno));
			Reset(false);
			return(eBSFerrCreateFile);
			}
#endif

if(m_hOutFile < 0)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate output file '%s'",pszOutFile);
	Reset(false);
	return(eBSFerrCreateFile);
	}
m_pszOutFile = pszOutFile;

if(bPEgen)
	{
#ifdef _WIN32
	m_hOutPEFile = open(pszOutPEFile,O_CREATETRUNC );
#else
	if((m_hOutPEFile = open(pszOutPEFile,O_RDWR | O_CREAT,S_IREAD | S_IWRITE))!=-1)
		if(ftruncate(m_hOutPEFile,0)!=0)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate %s - %s",pszOutPEFile,strerror(errno));
				Reset(false);
				return(eBSFerrCreateFile);
				}
#endif

	if(m_hOutPEFile < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate output file '%s'",pszOutPEFile);
		Reset(false);
		return(eBSFerrCreateFile);
		}
	m_pszOutPEFile = pszOutPEFile;
	}
else
	{
	m_pszOutPEFile = NULL;
	m_hOutPEFile = -1;
	}

if(pszHammFile != NULL && pszHammFile[0] != '\0')
	{
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Loading Hamming edit distances from file '%s'",pszHammFile);
	if((Rslt = LoadHammings(pszHammFile)) < eBSFSuccess)
		{
		Reset(false);
		return(Rslt);
		}
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Loading Hamming edit distances completed");
	}

// Don't bother checkpointing (write to file in batches every N reads generated) unless dynamically generating Hamming or profiling
// Checkpointing is only used if the processing is expected to take some time to complete
if(DfltHamming < 0 && (pszHammFile == NULL || pszHammFile[0] == '\0'))
	MaxReadsPerBatch = cMaxHammingBatchSize;
else
	{
	if(PMode != ePMStandard)
		MaxReadsPerBatch = min(m_TotReqReads,cMaxProfileBatchSize);
	else
		MaxReadsPerBatch = min(m_TotReqReads,cMaxBatchSize);
	}
MaxReadsPerBatch *= NumThreads;


// Allocate to hold all reads if bDedupe is TRUE even though they will be checkpointed, and written to disk, every cChkNumReads. This is because
// when deduping the reads the deduping needs to be over all reads and not just the reads in the current checkpointed batch
if(bDedupe)			// if will be deduping then generate some extra reads assuming that a few will be dups and thus be removed
	{
	NumReadsReq = (int)(((INT64)m_TotReqReads * 120)/100);
	m_NumReadsAllocd = NumReadsReq + 100;
	}
else
	{
	NumReadsReq = m_TotReqReads;
	m_NumReadsAllocd = MaxReadsPerBatch + 100;
	}

m_AllocdMemReads = (INT64)m_NumReadsAllocd * sizeof(tsSimRead);
#ifdef _WIN32
m_pSimReads = (tsSimRead *) malloc((size_t)m_AllocdMemReads);
if(m_pSimReads == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: Memory allocation of %lld bytes for simulated reads failed",(INT64)m_AllocdMemReads);
	Reset(false);
	return(eBSFerrMem);
	}
#else
// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
m_pSimReads = (tsSimRead *)mmap(NULL,(size_t)m_AllocdMemReads, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
if(m_pSimReads == MAP_FAILED)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: Memory allocation of %lld bytes for simulated reads through mmap()  failed",(INT64)m_AllocdMemReads,strerror(errno));
	m_pSimReads = NULL;
	Reset(false);
	return(eBSFerrMem);
	}
#endif

memset(m_pSimReads,0,(size_t)m_AllocdMemReads);
if(bReadHamDist) {
	tsSimRead *pRead;
	int ReadIdx;
	m_pHamDistFreq = new UINT32 [m_NumReadsAllocd * (ReadLen+1)]; // need to allow for distances ranging 0..ReadLen inclusive
	if(m_pHamDistFreq == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: Memory allocation of %lld bytes for simulated reads failed",(INT64)m_AllocdMemReads);
		Reset(false);
		return(eBSFerrMem);
		}

	pRead = m_pSimReads;
	for(ReadIdx = 0; ReadIdx < m_NumReadsAllocd; ReadIdx++,pRead++)
		pRead->pHamDistFreq = &m_pHamDistFreq[ReadIdx * (ReadLen+1)];
	}

#ifdef _WIN32
if((m_hMtxDedupe = CreateMutex(NULL,false,NULL))==NULL)
#else
if(pthread_mutex_init (&m_hMtxDedupe,NULL)!=0)
#endif
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: unable to create mutex");
	Reset(false);
	return(eBSFerrInternal);
	}

CurNumGenReads = 0;
PrevNumGenReads = 0;
ReadsOfs = 0;
TotReportedReads = 0;
bFirst =true;
TRandomCombined<TRanrotWGenerator,TRandomMersenne> RGseeds((int)time(0));
do {
	// initialise all worker thread parameters and start the threads
	if(!bDedupe)
		{
		ReadsOfs = 0;
		ReadsCnt = min((m_TotReqReads - TotReportedReads),MaxReadsPerBatch);
		}
	else
		{
		if(bUseLocateHamming)
			ReadsCnt = min(1000+(m_TotReqReads - TotReportedReads),MaxReadsPerBatch);
		else
			ReadsCnt = min((int)(((INT64)(m_TotReqReads - TotReportedReads) * 105)/100),MaxReadsPerBatch);
		if(ReadsCnt < 10000 && (ReadsOfs + 10000) < m_NumReadsAllocd)
			ReadsCnt = 10000;
		}
	NumReadsReq = 0;
	ExhustedChroms = 0;
	for(ThreadIdx = 0; ThreadIdx < NumThreads; ThreadIdx++)
		{
		pCurThread = &WorkerThreads[ThreadIdx];
		ReadsPerThread = ReadsCnt/(NumThreads-ThreadIdx);
		ReadsCnt -= ReadsPerThread;
		memset(pCurThread,0,sizeof(tsWorkerPars));
		pCurThread->RandSeed = RGseeds.IRandom(1,INT_MAX);
		pCurThread->ThreadIdx = ThreadIdx;
		pCurThread->DfltHamming = DfltHamming;
		pCurThread->bUseLocateHamming = bUseLocateHamming;
		pCurThread->bDedupe = bDedupe;
		pCurThread->bReadHamDist = bReadHamDist;
		pCurThread->ReadLen = ReadLen;
		pCurThread->Region = Region;
		pCurThread->UpDnRegLen = UpDnRegLen;
		pCurThread->CutMax = CutMax;
		pCurThread->CutMin = CutMin;

		pCurThread->bPEgen = bPEgen;
		pCurThread->PEmin = PEmin;
		pCurThread->PEmax = PEmax;
		pCurThread->NumGenReads=0;
		pCurThread->NumReqReads=ReadsPerThread;
		NumReadsReq += ReadsPerThread;
		pCurThread->pReads = &m_pSimReads[ReadsOfs];
		ReadsOfs += ReadsPerThread;
		pCurThread->PMode=PMode;
		pCurThread->Strand=Strand;
		pCurThread->bMaxIters = false;

	#ifdef _WIN32
		pCurThread->threadHandle = (HANDLE)_beginthreadex(NULL,0x0fffff,WorkerThread,pCurThread,0,&pCurThread->threadID);
	#else
		pCurThread->threadRslt =	pthread_create(&pCurThread->threadID , NULL , WorkerThread , pCurThread );
	#endif
		}
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Randomly selecting up to %s%d reads...",bFirst?" ":" another ",NumReadsReq);
	bFirst = false;


	// wait for all threads to terminate - be patient, could be a long, long wait if dynamic Hamming dist determinations
    PrevNumGenReads = 0;
	for(ThreadIdx = 0; ThreadIdx < NumThreads; ThreadIdx++)
		{
		pCurThread = &WorkerThreads[ThreadIdx];

		// report on number of reads generated every 60 secs assuming that at least 5 reads were generated in that time period
#ifdef _WIN32
		while(WAIT_TIMEOUT == WaitForSingleObject( pCurThread->threadHandle, 60000))
			{
			WaitForSingleObject(m_hMtxDedupe,INFINITE);
			CurNumGenReads = m_CurNumGenReads;
			ReleaseMutex(m_hMtxDedupe);
			if(CurNumGenReads > (PrevNumGenReads+1))
				{
				gDiagnostics.DiagOut(eDLInfo,gszProcName,"Progress: %u reads generated",CurNumGenReads);
				PrevNumGenReads = CurNumGenReads;
				}
			}
		CloseHandle(pCurThread->threadHandle);
#else
		struct timespec ts;
		int JoinRlt;
		clock_gettime(CLOCK_REALTIME, &ts);
		ts.tv_sec += 60;
		while((JoinRlt = pthread_timedjoin_np(pCurThread->threadID, NULL, &ts)) != 0)
			{
			pthread_mutex_lock(&m_hMtxDedupe);
			CurNumGenReads = m_CurNumGenReads;
			pthread_mutex_unlock(&m_hMtxDedupe);
			if(CurNumGenReads > (PrevNumGenReads+1))
				{
				gDiagnostics.DiagOut(eDLInfo,gszProcName,"Progress: %u reads generated",CurNumGenReads);
				PrevNumGenReads = CurNumGenReads;
				}
			ts.tv_sec += 60;
			}

#endif
		if(pCurThread->bMaxIters)
			ExhustedChroms += 1;
		}

	if(ExhustedChroms < NumThreads)
		{
		ReportedReads = ReportReads(bPEgen,		// true if paired end simulated reads being simulated
		        Region,						// Hamming regional processing?
				ReadLen,					// read length
			    TotReportedReads,			// number of reads thus far reported on
				m_TotReqReads,				// maximum number of reads required to be reported on
				bDedupe,					// true if reads are to be deduped
				bReadHamDist,				// true if hamming distributions from each sampled read
				ReadsOfs);					// number of reads to report on from this batch
		TotReportedReads += ReportedReads;
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Current block processing: %1.7d total: %1.8d",ReportedReads,TotReportedReads);
		}
	else
		{
		gDiagnostics.DiagOut(eDLWarn,gszProcName,"Unable to generate simulated reads as assembly chromosomes not of sufficent length from which to randomly select reads");
		break;
		}
	}
while(TotReportedReads < m_TotReqReads && ReadsOfs < m_NumReadsAllocd);

#ifdef _WIN32
CloseHandle(m_hMtxDedupe);
#else
pthread_mutex_destroy(&m_hMtxDedupe);
#endif
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Total simulated reads generated: %1.8d",TotReportedReads);

if(m_hOutFile != -1)
	{
#ifdef _WIN32
	_commit(m_hOutFile);
#else
	fsync(m_hOutFile);
#endif
	close(m_hOutFile);
	m_hOutFile = -1;
	}
if(m_hOutPEFile != -1)
	{
#ifdef _WIN32
	_commit(m_hOutPEFile);
#else
	fsync(m_hOutPEFile);
#endif
	close(m_hOutPEFile);
	m_hOutPEFile = -1;
	}

if(Region != eMEGRAny)
	{
	if(m_pBEDFile != NULL)
		{
		delete m_pBEDFile;
		m_pBEDFile = NULL;
		}
	}

// show distributions
gDiagnostics.DiagOut(eDLInfo,gszProcName,"\"ReadOfs\",\"InducedSubs\"");
int Idx;
for(Idx=0;Idx<ReadLen;Idx++)
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"%d,%d",Idx,m_InducedErrPsnDist[Idx]);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"\"SubDist\",\"Cnt\"");

for(Idx=0;Idx<cNumProfEls;Idx++)
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"%d,%d",Idx,m_InducedErrDist[Idx]);
Reset(true);
return(eBSFSuccess);
}


// Workerthread
#ifdef _WIN32
unsigned __stdcall WorkerThread(void * pThreadPars)
#else
void *WorkerThread(void * pThreadPars)
#endif
{
tsWorkerPars *pWorkerPars = (tsWorkerPars *)pThreadPars;
etSeqBase ReadSeq[cMaxReadLen+1];
etSeqBase *pReadSeq;
etSeqBase Nuc;
tsSimRead *pRead;
int HammingDist;
int HammingDistW;
int RandChrom;

int BEDChromID;
int Features;

int RandCutSite1;
int RandCutLen;
int RandCutSite2;
int RandStrand;				// 0 if '+', 1 if '-'

int PEFragSize;
int PERandCutSite1;
int PERandCutLen;
int PERandCutSite2;
int PERandStrand;

int OctIdx;
int CurNumGenReads;

double RandCutProb;
double ProfScore;
tsChromSeq *pChromSeq;
etSeqBase *pSeq;
etSeqBase *pSeqSite;

int TargPsn;
int IdxLo;
int IdxHi;
UINT64 NumChromIters;

CurNumGenReads = 0;
NumChromIters = 0;
pRead = pWorkerPars->pReads;
pRead += pWorkerPars->NumGenReads;
pWorkerPars->bMaxIters = false;
TRandomCombined<TRanrotWGenerator,TRandomMersenne> RG(pWorkerPars->RandSeed);
while(pWorkerPars->NumGenReads < pWorkerPars->NumReqReads)
	{
	if(NumChromIters++ > ((UINT64)m_NumChromSeqs * 50))
		{
		pWorkerPars->bMaxIters = true;			// flag this thread is terminating because it has exhusted attempts to find chrom from which read can be simulated
		break;
		}
		// first randomly choose chromosome
	RandChrom = (int)RG.IRandom(1,m_GenomeScaledLen);

	IdxLo = 0;
	IdxHi = m_NumChromSeqs-1;
	do {
		TargPsn = (IdxLo + IdxHi) / 2;
		pChromSeq = &m_pChromSeqs[TargPsn];
		if(pChromSeq->ScaledStart <= RandChrom && (pChromSeq->ScaledStart + pChromSeq->ScaledLen) >= RandChrom)
			break;

		if(pChromSeq->ScaledStart > RandChrom)
			IdxHi = TargPsn - 1;
		else
			IdxLo = TargPsn+1;
		}
	while(IdxHi >= IdxLo);

	if(pWorkerPars->Strand != '*' && pChromSeq->Strand != '*' && pWorkerPars->Strand != pChromSeq->Strand)
		continue;

	// randomly choose cut length?
	if(pWorkerPars->CutMin != pWorkerPars->CutMax)
		{
		RandCutLen = (int)RG.IRandom(pWorkerPars->CutMin,pWorkerPars->CutMax);
		if(pWorkerPars->bPEgen)
			PERandCutLen = (int)RG.IRandom(pWorkerPars->CutMin,pWorkerPars->CutMax);
		}
	else
		{
		RandCutLen = pWorkerPars->CutMin;
		PERandCutLen = RandCutLen;
		}

		// if paired end then choose the fragment length
	if(pWorkerPars->bPEgen)
		{
		PEFragSize = (int)RG.IRandom(pWorkerPars->PEmin,pWorkerPars->PEmax);
		if(PEFragSize < (min(RandCutLen,PERandCutLen) + 1))		// allow for user simulating overlapped paired reads as required for AllpathsLG
			continue;
		}
	else
		PEFragSize = 0;

	// skip any extremely short chromosomes - could be a contig?
	if(pWorkerPars->bPEgen)
		{
		if(pChromSeq->Len < (PEFragSize + 20))
			continue;
		}
	else
		if(pChromSeq->Len < (pWorkerPars->CutMax + 20))
			continue;

	if(pChromSeq->RelDensity < 1)				// don't bother to generate reads if relative density too low
		continue;

	if(pChromSeq->RelDensity < 1000)				// if less than 1000 then don't accept this putative read if RelDensity < rand(0,999)
		{
		if(pChromSeq->RelDensity < (int)RG.IRandom(0,999))
			continue;
		}

	// try and replicate, very crude attempt!, both the variance in transcript levels and the clumped distribution of reads within transcripts
	// as observed with RNASeq sequenced read distributions
	// firstly the transcript level is simply a linear function of the first 8bp at the start of the transcript sequence
	if(m_DistCluster > 0 && pChromSeq->RelDensity == 1000)
		{
		int TransLev = 0;
		TransLev = GenPSeqIdx(8,&m_pGenomeSeq[pChromSeq->SeqOfs]);
		// overall transcript level now determined
		if(RG.IRandom(0,0x010000) >= TransLev)
			continue;
		}

	// randomly choose strand?
	if(pWorkerPars->Strand == '*')
		RandStrand = (int)RG.IRandom(0,1);
	else
		RandStrand = pWorkerPars->Strand == '-' ? 1 : 0;
	if(pWorkerPars->bPEgen)
		PERandStrand = RandStrand == 0 ? 1 : 0;

	// randomly choose initial cut sites
	// note that specified range is such that '+' and '-' start/end loci after allowing for cut lengths will always be on the chromosome
	if(pWorkerPars->bPEgen)
		{
		if(RandStrand == 0)
			{
			RandCutSite1 = (int)RG.IRandom(7,pChromSeq->Len - (PEFragSize + 7));
			PERandCutSite1 = RandCutSite1 + PEFragSize - PERandCutLen;
			}
		else
			{
			RandCutSite1 = (int)RG.IRandom(7 + PEFragSize - RandCutLen,pChromSeq->Len - RandCutLen - 7);
			PERandCutSite1 = RandCutSite1 - (PEFragSize - RandCutLen) - 1;
			}
		}
	else
		{
		RandCutSite1 = (int)RG.IRandom(7,pChromSeq->Len - (RandCutLen + 7));
		PERandCutSite1 = RandCutSite1;
		}

	RandCutSite2 = RandCutSite1 + RandCutLen;
	PERandCutSite2 = PERandCutSite1 + PERandCutLen;


	// filter by region here
	if(pWorkerPars->Region != eMEGRAny)
		{
		if(m_pBEDFile != NULL)
			{
			if((BEDChromID = m_pBEDFile->LocateChromIDbyName(pChromSeq->szChromName)) < 1)
				continue;
			Features = m_pBEDFile->GetFeatureBits(BEDChromID,RandCutSite1+(RandCutLen/2),RandCutSite1+(RandCutLen/2),cRegionFeatBits,pWorkerPars->UpDnRegLen);
			switch(pWorkerPars->Region) {
				case eMEGRCDS:			// part of feature overlaps CDS
					if(!(Features & cFeatBitCDS))
						continue;
					break;
				case eMEGR5UTR:			// part of feature overlaps 5'UTR
					if(!(Features & cFeatBit5UTR) || (Features & cFeatBitCDS))
						continue;
					break;
				case eMEGR3UTR:			// part of feature overlaps 3'UTR
					if(!(Features & cFeatBit3UTR) || (Features & (cFeatBitCDS | cFeatBit5UTR)))
						continue;
					break;
				case eMEGRIntrons:		// part of feature overlaps Intron
					if(!(Features & cFeatBitIntrons) || (Features & (cFeatBitCDS | cFeatBit5UTR | cFeatBit3UTR)))
						continue;
					break;
				case eMEGRUpstream:		// part of feature overlaps 5'upstream regulatory
					if(!(Features & cFeatBitUpstream) || (Features & (cFeatBitCDS | cFeatBit5UTR | cFeatBit3UTR | cFeatBitIntrons)))
						continue;
					break;
				case eMEGRDnstream:		// part of feature overlaps 3'downstream regulatory
					if(!(Features & cFeatBitDnstream) || (Features & (cFeatBitCDS | cFeatBit5UTR | cFeatBit3UTR | cFeatBitIntrons | cFeatBitUpstream)))
						continue;
					break;
				case eMEGRIntergenic:	// part of feature overlaps intergenic
					if(Features & (cFeatBitCDS | cFeatBit5UTR | cFeatBit3UTR | cFeatBitIntrons | cFeatBitUpstream | cFeatBitDnstream))
						continue;
					break;
				}
			}
		}

		// now for the clumping factor
		// clumps are distributed along the length of the transcript into non-overlapping clustering bins of width m_DistCluster
	if(m_DistCluster && pChromSeq->Len > m_DistCluster)
		{
		int ClumpProb;
		int ClumpOfs;
		ClumpOfs = (RandCutSite1 / m_DistCluster) * m_DistCluster; // ClumpOfs is the start of bin, of m_DistCluster len, which would contain RandCutSite
		if(((ClumpOfs/m_DistCluster) % 5) == 2)					   // every 5th bin starting at the second is the bin containing the clustered reads
			{
			pSeq = &m_pGenomeSeq[pChromSeq->SeqOfs  + ClumpOfs];	// bases at the start of clustering bin used to determine the prob of reads in that bin
			ClumpProb = GenPSeqIdx(5,pSeq);
			}
		else
			ClumpProb = 0;											// reads outside of clustering bins are given a no chance...
		if(RG.IRandom(1,1024) > ClumpProb)
				continue;
		}

	if(pWorkerPars->PMode == ePMProfRand || pWorkerPars->PMode == ePMProfProf)
		{
		// randomly choose site1 cut prob
		RandCutProb = RG.Random();
		pSeq =  &m_pGenomeSeq[pChromSeq->SeqOfs];
		pSeqSite = &pSeq[RandCutSite1-4];
		if(!RandStrand)
			OctIdx = GenPSeqIdx(8,pSeqSite);
		else
			OctIdx = GenMSeqIdx(8,pSeqSite);
		if(OctIdx < 0)		// -1 if pSeqSite contained 'n'
			continue;
		ProfScore = m_pProfSel[OctIdx];
		if(ProfScore < RandCutProb)
			continue;
		}

	pSeq = &m_pGenomeSeq[pChromSeq->SeqOfs];

	if(pWorkerPars->PMode == ePMProfProf || pWorkerPars->PMode == ePMRandProf)
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
		ProfScore = m_pProfSel[OctIdx];
		if(ProfScore < RandCutProb)
			continue;
		}

		// ensure that this sequence will be exactly alignable - slough if any contained 'n's
	pSeqSite = &pSeq[RandCutSite1];
	pReadSeq = ReadSeq;
	for(OctIdx = RandCutSite1; OctIdx < RandCutSite2; OctIdx++,pSeqSite++)
		{
		if((Nuc = (*pSeqSite & NUCONLYMSK)) > eBaseT)
			break;
		*pReadSeq++ = Nuc;
		}
	if(OctIdx < RandCutSite2)
		continue;
	*pReadSeq = eBaseEOS;

	// if generating paired ends then check if partner read contains any 'n's - if so then slough
	if(pWorkerPars->bPEgen)
		{
		pSeqSite = &pSeq[PERandCutSite1];
		for(OctIdx = PERandCutSite1; OctIdx < PERandCutSite2; OctIdx++,pSeqSite++)
			{
			if((Nuc = (*pSeqSite & NUCONLYMSK)) > eBaseT)
				break;
			}
		if(OctIdx < PERandCutSite2)
			continue;
		}

	if(pWorkerPars->bReadHamDist)
		{
		if(m_PMode == ePMSampHamm)
			{
			HammingDistW = MinHammingDistW(pWorkerPars->ReadLen,pWorkerPars->ReadLen,pChromSeq->ChromID,RandCutSite1,ReadSeq);
			CSeqTrans::ReverseComplement(RandCutLen,ReadSeq);
			HammingDist = MinHammingDistC(HammingDistW,pWorkerPars->ReadLen,pChromSeq->ChromID,RandCutSite1,ReadSeq);
			}
		else
			{
			memset(pRead->pHamDistFreq,0,sizeof(UINT32) * (pWorkerPars->ReadLen + 1));
			HammingDistW = HammingDistCntsW(pWorkerPars->ReadLen,pChromSeq->ChromID,RandCutSite1,ReadSeq,pRead->pHamDistFreq);
			CSeqTrans::ReverseComplement(RandCutLen,ReadSeq);
			HammingDist = HammingDistCntsC(pWorkerPars->ReadLen,pChromSeq->ChromID,RandCutSite1,ReadSeq,pRead->pHamDistFreq);
			if(HammingDistW < HammingDist)
				HammingDist = HammingDistW;
			}
		}
	else
		{
		if(pWorkerPars->bUseLocateHamming)
			{
			HammingDist = LocateHamming(pChromSeq->szChromName,RandCutSite1);
			if(HammingDist < 0)	// need to tighten up on this: -1 returned if no such chromosome or no hamming for subsequence at specified loci
				continue;		// in this experimental release these chroms are simply sloughed
			}
		else
			{
			if(pWorkerPars->DfltHamming < 0)
				{
				// now determine the minimal Hamming distance of this sequence from any other sequence in the targeted genome
				// on both the watson and crick strands
				// Dynamic Hamming distance determination is an extremely slow process, how to speed it up?????
				// small speedup by setting initial hamming to be readlen/2. With length 36, currently the largest actual hamming observed
				// was around 14
				HammingDistW = MinHammingDistW(pWorkerPars->ReadLen/2,pWorkerPars->ReadLen,pChromSeq->ChromID,RandCutSite1,ReadSeq);
				CSeqTrans::ReverseComplement(RandCutLen,ReadSeq);
				HammingDist = MinHammingDistC(HammingDistW,pWorkerPars->ReadLen,pChromSeq->ChromID,RandCutSite1,ReadSeq);
				}
			else
				HammingDist = pWorkerPars->DfltHamming;
			}
		if(pWorkerPars->bDedupe && HammingDist == 0)
			continue;
		}


	// we have a sequence starting at RandCutSite1 and ending at RandCutSite2-1 which is of length RandCutLen
	// if non-duplicates required then mark subsequence as selected
	if(pWorkerPars->bDedupe)
		{
		pSeqSite = &pSeq[RandCutSite1];
#ifdef _WIN32
		WaitForSingleObject(m_hMtxDedupe,INFINITE);
#else
		pthread_mutex_lock(&m_hMtxDedupe);
#endif
		if(*pSeqSite & SSSELECTED)	// check if this subsequence already selected
			{
#ifdef _WIN32
			ReleaseMutex(m_hMtxDedupe);
#else
			pthread_mutex_unlock(&m_hMtxDedupe);
#endif
			continue;
			}
		*pSeqSite |= SSSELECTED;
#ifdef _WIN32
		ReleaseMutex(m_hMtxDedupe);
#else
		pthread_mutex_unlock(&m_hMtxDedupe);
#endif
		}

	pRead->Status = 0;
	pRead->ChromSeqID = pChromSeq->ChromSeqID;
	pRead->ChromID = pChromSeq->ChromID;
	pRead->Strand = RandStrand;
	pRead->StartLoci = RandCutSite1;
	pRead->EndLoci = RandCutSite2 - 1;
	pRead->Len = RandCutLen;
	pRead->HammingDist = HammingDist;
	pRead->pSeq = &pSeq[RandCutSite1];
	pWorkerPars->NumGenReads += 1;
	pRead->FlgPE2 = 0;
	if(pWorkerPars->bPEgen)
		pRead->pPartner = &pRead[1];
	else
		pRead->pPartner = NULL;
	pRead += 1;

	if(pWorkerPars->bPEgen)
		{
		pRead->pPartner = &pRead[-1];
		pRead->FlgPE2 = 1;
		pRead->Status = 0;
		pRead->ChromSeqID = pChromSeq->ChromSeqID;
		pRead->ChromID = pChromSeq->ChromID;
		pRead->Strand = PERandStrand;
		pRead->StartLoci = PERandCutSite1;
		pRead->EndLoci = PERandCutSite2 - 1;
		pRead->Len = PERandCutLen;
		pRead->HammingDist = HammingDist;
		pRead->pSeq = &pSeq[pRead->StartLoci];
		pWorkerPars->NumGenReads += 1;
		pRead += 1;
		}

	CurNumGenReads += 1;
	NumChromIters = 0;

	// time to let main thread know that some progress is being made?
	if(pWorkerPars->bReadHamDist == true || CurNumGenReads >= 500)
		{
#ifdef _WIN32
		WaitForSingleObject(m_hMtxDedupe,INFINITE);
#else
		pthread_mutex_lock(&m_hMtxDedupe);
#endif

		m_CurNumGenReads += CurNumGenReads;
		CurNumGenReads = 0;
#ifdef _WIN32
		ReleaseMutex(m_hMtxDedupe);
#else
		pthread_mutex_unlock(&m_hMtxDedupe);
#endif
		}
	}
#ifdef _WIN32
_endthreadex(0);
return(eBSFSuccess);
#else
pthread_exit(NULL);
#endif
}



int		// returned minimum Hamming distance of Watson sequence - note the check for self-loci otherwise hamming would always be 0!
MinHammingDistW(int MinHamming,	// initial minimum Hamming distance
			   int ReadLen,		// read length
			   int ChromID,     // from which chromosome was this read was derived
			   int ReadLoci,	// read start loci
			   etSeqBase *pRead) // read sequence
{
int ChromIdx;
int SeqIdx;
int ChromSeqIdx;
int CurHamming;
int EndIdx;
tsChromSeq *pChrom;
etSeqBase *pChromBase;
etSeqBase *pChromSeq;
etSeqBase *pReadBase;

if(MinHamming <= 0)
	return(0);

pChrom = &m_pChromSeqs[0];
for(ChromIdx = 0; ChromIdx < m_NumChromSeqs; ChromIdx++, pChrom++)
	{
	pChromSeq =  &m_pGenomeSeq[pChrom->SeqOfs];
	EndIdx = pChrom->Len - ReadLen;
	for(ChromSeqIdx=0;ChromSeqIdx <= EndIdx; ChromSeqIdx++,pChromSeq++)
		{
		if(ChromSeqIdx == ReadLoci && ChromID == pChrom->ChromID)
			continue;

		pChromBase = pChromSeq;
		pReadBase = pRead;
		CurHamming = 0;
		for(SeqIdx=0;SeqIdx < ReadLen; SeqIdx++,pChromBase++,pReadBase++)
			if((*pChromBase & NUCONLYMSK) != (*pReadBase & NUCONLYMSK) &&
				++CurHamming >= MinHamming)
					break;
		if(CurHamming < MinHamming)
			{
			MinHamming = CurHamming;
			if(MinHamming == 0)
				return(MinHamming);
			}
		}
	}
return(MinHamming);
}

int		// returned minimum Hamming distance of Crick (reverse complement of Watson) sequence - note no check for self-loci
MinHammingDistC(int MinHamming,	// initial minimum Hamming distance
			   int ReadLen,		// read length
			   int ChromID,     // from which chromosome was this read was derived
			   int ReadLoci,	// read start loci
			   etSeqBase *pRead) // read sequence
{
int ChromIdx;
int SeqIdx;
int ChromSeqIdx;
int CurHamming;
int EndIdx;
tsChromSeq *pChrom;
etSeqBase *pChromBase;
etSeqBase *pChromSeq;
etSeqBase *pReadBase;

if(MinHamming <= 0)
	return(0);

pChrom = &m_pChromSeqs[0];
for(ChromIdx = 0; ChromIdx < m_NumChromSeqs; ChromIdx++, pChrom++)
	{
	pChromSeq =  &m_pGenomeSeq[pChrom->SeqOfs];
	EndIdx = pChrom->Len - ReadLen;
	for(ChromSeqIdx=0;ChromSeqIdx <= EndIdx; ChromSeqIdx++,pChromSeq++)
		{
		pChromBase = pChromSeq;
		pReadBase = pRead;
		CurHamming = 0;
		for(SeqIdx=0;SeqIdx < ReadLen; SeqIdx++,pChromBase++,pReadBase++)
			if((*pChromBase & NUCONLYMSK) != (*pReadBase & NUCONLYMSK) &&
				++CurHamming >= MinHamming)
					break;
		if(CurHamming < MinHamming)
			{
			MinHamming = CurHamming;
			if(MinHamming == 0)
				return(MinHamming);
			}
		}
	}
return(MinHamming);
}


int		// minimum Hamming for Watson strand
HammingDistCntsW(int ReadLen,		// read length
			 int ChromID,     // from which chromosome was this read was derived
			 int ReadLoci,	// read start loci
			 etSeqBase *pRead, // read sequence
			 UINT32 *pHammDist)	// where to return hamming dist counts (assumes at least ReadLen elements)
{
int MinHamming;
int ChromIdx;
int SeqIdx;
int ChromSeqIdx;
int CurHamming;
int EndIdx;
tsChromSeq *pChrom;
etSeqBase *pChromBase;
etSeqBase *pChromSeq;
etSeqBase *pReadBase;

MinHamming = ReadLen;
pChrom = &m_pChromSeqs[0];
for(ChromIdx = 0; ChromIdx < m_NumChromSeqs; ChromIdx++, pChrom++)
	{
	pChromSeq =  &m_pGenomeSeq[pChrom->SeqOfs];
	EndIdx = pChrom->Len - ReadLen;
	for(ChromSeqIdx=0;ChromSeqIdx <= EndIdx; ChromSeqIdx++,pChromSeq++)
		{
		if(ChromSeqIdx == ReadLoci && ChromID == pChrom->ChromID)
			continue;

		pChromBase = pChromSeq;
		pReadBase = pRead;
		CurHamming = 0;
		for(SeqIdx=0;SeqIdx < ReadLen; SeqIdx++,pChromBase++,pReadBase++)
			if((*pChromBase & NUCONLYMSK) != (*pReadBase & NUCONLYMSK))
				CurHamming++;
		pHammDist[CurHamming] += 1;
		if(MinHamming > CurHamming)
			MinHamming = CurHamming;
		}
	}
return(MinHamming);
}

int		// minimum Hamming for Crick strand
HammingDistCntsC(int ReadLen,		// read length
			   int ChromID,     // from which chromosome was this read was derived
			   int ReadLoci,	// read start loci
			   etSeqBase *pRead, // read sequence
   			 UINT32 *pHammDist)	// where to return hamming dist counts (assumes at least ReadLen elements)
{
int MinHamming;
int ChromIdx;
int SeqIdx;
int ChromSeqIdx;
int CurHamming;
int EndIdx;
tsChromSeq *pChrom;
etSeqBase *pChromBase;
etSeqBase *pChromSeq;
etSeqBase *pReadBase;

MinHamming = ReadLen;
pChrom = &m_pChromSeqs[0];
for(ChromIdx = 0; ChromIdx < m_NumChromSeqs; ChromIdx++, pChrom++)
	{
	pChromSeq = &m_pGenomeSeq[pChrom->SeqOfs];
	EndIdx = pChrom->Len - ReadLen;
	for(ChromSeqIdx=0;ChromSeqIdx <= EndIdx; ChromSeqIdx++,pChromSeq++)
		{
		pChromBase = pChromSeq;
		pReadBase = pRead;
		CurHamming = 0;
		for(SeqIdx=0;SeqIdx < ReadLen; SeqIdx++,pChromBase++,pReadBase++)
			if((*pChromBase & NUCONLYMSK) != (*pReadBase & NUCONLYMSK))
				CurHamming+=1;
		pHammDist[CurHamming] += 1;
		if(MinHamming > CurHamming)
			MinHamming = CurHamming;
		}
	}
return(MinHamming);
}

int
LocateRevCplSeq(int Len,etSeqBase *pProbe,int NumSortedReads)
{
int Rslt;
etSeqBase RevSeq[cMaxReadLen+1];
for(int Idx = 0; Idx < Len; Idx++,pProbe++)
	RevSeq[Idx] = *pProbe & NUCONLYMSK;
CSeqTrans::ReverseComplement(Len,RevSeq);
Rslt = LocateFirstExact(RevSeq,Len,0,NumSortedReads-1);
return(Rslt);
}



int			// index of exactly matching probe or -1 if no match
LocateFirstExact(etSeqBase *pProbe,				// pts to probe sequence
				 int ProbeLen,					// probe length to exactly match over
				  int IdxLo,					// low index in m_pSimReads
				  int IdxHi)					// high index in m_pSimReads
{
etSeqBase *pEl1;
etSeqBase *pEl2;

int CmpRslt;

int TargPsn;
int LowPsn;

pEl1 = pProbe;
do {
	TargPsn = (IdxLo + IdxHi) / 2;
	pEl2 = m_pSimReads[TargPsn].pSeq;
	CmpRslt = CmpSeqs(ProbeLen,pEl1,pEl2);

	if(!CmpRslt)	// if have a match but might not be the lowest indexed match
		{
		if(TargPsn == 0 || IdxLo == TargPsn) // check if lowest
			return(TargPsn);
		LowPsn = LocateFirstExact(pProbe,ProbeLen,IdxLo,TargPsn - 1);
		return(LowPsn < 0 ? TargPsn : LowPsn);
		}

	if(CmpRslt < 0)
		IdxHi = TargPsn - 1;
	else
		IdxLo = TargPsn+1;
	}
while(IdxHi >= IdxLo);

return(-1);	// unable to locate any instance of pProbe
}



// SortSimLoci
// Sort simulated reads by chrom,loci,len
int
SortSimLoci(const void *arg1, const void *arg2)
{
tsSimRead *pEl1 = (tsSimRead *)arg1;
tsSimRead *pEl2 = (tsSimRead *)arg2;
if(pEl1->ChromID < pEl2->ChromID)
	return(-1);
if(pEl1->ChromID > pEl2->ChromID)
	return(1);
if(pEl1->StartLoci < pEl2->StartLoci)
	return(-1);
if(pEl1->StartLoci > pEl2->StartLoci)
	return(1);
if(pEl1->Len < pEl2->Len)
	return(-1);
if(pEl1->Len > pEl2->Len)
	return(1);
return(0);
}


// SortSimReads
// Sort simulated reads by sequence
int
SortSimReads(const void *arg1, const void *arg2)
{
tsSimRead *pEl1 = (tsSimRead *)arg1;
tsSimRead *pEl2 = (tsSimRead *)arg2;
etSeqBase *pSeq1 = pEl1->pSeq;
etSeqBase *pSeq2 = pEl2->pSeq;
etSeqBase Base1;
etSeqBase Base2;
int Idx;
int Len = min(pEl1->Len,pEl2->Len);
for(Idx=0; Idx < Len; Idx++)
	if((Base1 = (*pSeq1++ & NUCONLYMSK)) != (Base2 = (*pSeq2++ & NUCONLYMSK)))
		return(Base1 - Base2);
return(0);
}

