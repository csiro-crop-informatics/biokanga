// gencentroidmetrics.cpp : Defines the entry point for the console application.
// Generates centroid metrics from multialignment file containing the reference, relative and outgroup species.
// An example of this would be to use hg18multiz17way alignments and process hg18 vs panTro1 with rheMac1 (outgroup).

#include "stdafx.h"
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#if _WIN32

#include "../libbiokanga/commhdrs.h"
#else
#include "../libbiokanga/commhdrs.h"
#endif

#include "./SNPsFile.h"

const unsigned int cProgVer = 300;		// increment with each release


const int cMARawBuffSize = 0x0ffffff;
const int cMALineSize    = 0x0fffff;
const int cMAFmaxSpecies =  3;    // maximum number of species that can be handled (ref + rel + optional outgroup)

const int cMinNMer = 3;			  // minimal length centroid oligo that can be handled
const int cDfltNMer = 5;		  // default ength centroid oligo that can be handled
const int cMaxNMer = 13;		  // maximal length centroid oligo that can be handled

const int cMaxProcMode = 2;		  // maximal processing supported mode
const int cMaxIncludeFiles = 10;  // maximun number of include region filter files
const int cMaxExcludeFiles = 10;  // maximun number of exclude region filter files
const int cRegionLen = 4;		  // number of ints required per region when capturing regional counts

const int cChromSeqLen = 0x0ffffff; // initial data allocation for holding complete chromosomes


typedef enum eProcMode {
	ePMCentroids = 0,			  // nMer centroid alignment stats  
	ePMnMerFreq,				  // nMer frequency distribution stats for genome
	ePMSNPs						  // Create SNP stats
} teProcMode;

CStopWatch gStopWatch;
CDiagnostics gDiagnostics;				// for writing diagnostics messages to log file
char gszProcName[_MAX_FNAME];			// process name



typedef struct TAG_sProcParams 	{
	char szOnlyChrom[cMaxDatasetSpeciesChrom]; // a single chromosome only processing - saves waiting for all!
	int OnlyChromID;				// chrom identifier for szOnlyChrom
	int	OnlyStart;					// if OnlyChromID and >= 0 then start offset along chromosome
	int	OnlyEnd;					// if OnlyChromID and >= 0 then end offset along chromosome
	char OnStrand;					// stats on which strand - if '*' then any strand otherwise '+' or '-' features only will be processed
	int RefSpeciesIdx;				// current index into pSeqs for the reference species
	char szSpecies[cMAFmaxSpecies][cMaxDatasetSpeciesChrom];	// species names of interest - other species are sloughed
									// only alignments with species included will be processed
									// first species is the reference species
	char szRefChrom[cMaxDatasetSpeciesChrom];
	int NumSpeciesAligned;			// actual number of species in current alignment
	int NumSpecies2Align;			// number of species to be aligned in current phase
	bool bFiltLoConfidence;			// true if low confidence subseqences to be filtered out
	bool bFilt1stLast;				// treat 1st and last subsequences as being low confidence
	int MinIdent;					// treat subsequences of less than this identity as being low confidence
	int MinSubSeqLen;				// subsequences of less than this length are treated as being low confidence
	bool bInferLoose;				// if true then if refbase == relbase accept refbase as ancestor even if outgroup differs
	int	 iWindowSize;				// sampling window size
	int NxtOutputOffset;			// when to next output results - i.e end of current window
	bool bOverlapNMers;				// if true then overlap NMers when generating genome statistics
	etSeqBase *pSeq;				// ptr to sequence
	int AllocdDataLen;				// how much memory (bytes) have been allocated to hold pSeq data
	int SeqOfs;						// offset into sequences
	int SeqLen;						// sequence length
	etSeqBase *pSeqs[cMAFmaxSpecies];  // ptrs to each sequence
	int NMer;						// what length oligos to process
	int FlankLen;					// number of nucleotides flanking each side of centroid in NMer
	int MerIdxMsk;
	int NumCnts;					// number of steps in pCntStepCnts
	int Regions;					// number of regions per step
	bool bCountsAvail;				// true if counts have been generated, reset when results have been output
	int *pCntStepCnts;
	CSNPsFile *pSNPs;				// pts to SNPs class
	CBEDfile *pBiobed;				// pts to biobed file containing feature classifications
	int BEDChromID;					// regional classification BED file chrom identifier
	int NumIncludes;				// number of biobed files containing regions to include
	int NumExcludes;				// number of biobed files containing regions to exclude
	CBEDfile *pIncludes[cMaxIncludeFiles];	// if opened biobed files for regions to include - all other regions are to be excluded
	CBEDfile *pExcludes[cMaxExcludeFiles];	// if opened biobed files for regions to exclude 
	int UpDnStreamLen;				// up/dn stream regional length when characterising
	bool bMultipleFeatBits;			// if false then stats only generated if a single feature bit is set - e.g if both CDS and introns overlapped then no stat generated
	bool bMultiPhase;				// false if single phase processing (ref+rel+outgroup in same alignment file)
	int CurPhase;					// if multiphase then the current phase (0..n)
	int hRsltsFile;					// write results into this file
	bool bPerChrom;					// true if stats to be generated per chromosome
	bool bGenEmptyRows;				// true if rows with 0 counts still to be generated in output
	bool bGenEmptyWindows;			// true if windows with no counts still to be generated in output
	} tsProcParams; 

int NormaliseInDelColumns(tsProcParams *pProcParams,int AlignLen,int NumSpecies);
bool ProcAlignBlock(int RefChromID,int RefChromOfs,int AlignLen,tsProcParams *pProcParams);
bool ProcessAlignment(int RefChromID,int SubRefOfs,int SubSeqLen,int SubSeqStartIdx,tsProcParams *pProcParams);
bool ProcessSubSeq(int RefChromID,int SubRefOfs,int SubSeqStartIdx,int MidBaseIdx,tsProcParams *pProcParams);
char *Idx2Bases(int StepIdx,int NumBases,char *pszBuffer);
int GetNumSubseqs(int AlignLen, tsProcParams *pProcParams);
int ProcessAlignments(char *pszMAF,tsProcParams *pProcParams);
int ParseNumSpecies(char *pszSpeciesList,tsProcParams *pProcParams);
int TrimQuotes(char *pszToTrim);
bool ChkOutputGenomeResults(char *pszChrom, int ChromOffset, tsProcParams *pProcParams);
bool ChkOutputResults(char *pszChrom, int ChromOffset, tsProcParams *pProcParams);
bool OutputGenomeResults(char *pszChrom,int ChromOfs,tsProcParams *pProcParams);
bool OutputResults(char *pszChrom,int ChromOfs,tsProcParams *pProcParams);
double CalcChiSqr(int Rows,int Cols,int *pCells);
bool CalcAllChiSqrs(tsProcParams *pProcParams);
int MapStep2Region(int Step);
bool ProcessSNPSequence(int RefChromID,int SubRefOfs,tsProcParams *pProcParams);


int
CreateAlignStats(char *pszRefRelFile,		// file containing ref + rel species, optionally the outgroup species
				 char *pszStatsFile,		// where to write out stats
				 char *pszBiobedFile,		// biobed file containing regional features
			 	 int NumIncludes,			// number of biobed files containing regions to include
  				 char **pszIncludeFiles,	// biobed files containing regions to include - default is to include everything
 				 int NumExcludes,			// number of biobed files containing regions to exclude
				 char **pszExcludeFiles,	// biobed files containing regions to exclude - default is to exclude nothing
				 int iWindowSize,			// sampling window size
				 char *pszSpeciesList,		// space or comma separated list of species, priority ordered 
				 int nMer,					// lengths of oligos to process
				 int RegLen,				// regulatory region length - up/dn stream of 5/3' 
				 bool MultipleFeatBits,		// if true then accept alignments in which multiple feature bits are set
				 bool bFiltLoConfidence,	// true if to filter out low confidence subsequences
			 	bool bFilt1stLast,			// treat 1st and last subsequences as being low confidence
				int MinSubSeqLen,			// subsequences of less than this length are treated as being low confidence
				int MinIdent,				// treat subsequences of less than this identity as being low confidence
				bool bInferLoose,			// if true then if refbase == relbase accept refbase as ancestor even if outgroup differs
				int Mode,					// processing mode
				bool bGenEmptyRows,			// true if rows with 0 counts still to be generated in output
				bool bGenEmptyWindows,		// true if windows with no counts still to be generated in output
				bool bPerChrom,				// true if counts to be generated per chromosome
				char OnlyStrand,			// any or '+','-' features only
				char *pszOnlyChrom,			// if not NULL then specifies single chromosome to process
				int	OnlyStart,				// if OnlyChromID and >= 0 then start offset along chromosome
				int	OnlyEnd);				// if OnlyChromID and >= 0 then end offset along chromosome

bool ProcessSequence(int RefChromID,int SubRefOfs,tsProcParams *pProcParams);

// CreateGenomeStats
// Create nMer distribution stats for complete genome
int
CreateGenomeStats(char *pszGenomeFile,		// file containing genome assembly
				 char *pszStatsFile,		// where to write out stats
				 char *pszBiobedFile,		// biobed file containing regional features
			 	 int NumIncludes,			// number of biobed files containing regions to include
  				 char **pszIncludeFiles,	// biobed files containing regions to include - default is to include everything
				 int NumExcludes,			// number of biobed files containing regions to exclude
				 char **pszExcludeFiles,	// biobed files containing regions to exclude - default is to exclude nothing
				 int iWindowSize,			// sampling window size
				 int nMer,					// subseq length
				 int RegLen,				// regulatory region length - up/dn stream of 5/3' 
				 bool bMultipleFeatBits,	// if true then accept alignments in which multiple feature bits are set
				 bool bOverlapNMers,		// if true then overlap NMers when generating genome statistics
				bool bGenEmptyRows,			// true if rows with 0 counts still to be generated in output
				bool bGenEmptyWindows,		// true if windows with no counts still to be generated in output
				bool bPerChrom,				// true if stats to be generated per chromosome
				char OnlyStrand,			// any or '+','-' features only
				char *pszOnlyChrom,			// if not NULL then specifies single chromosome to process
				int	OnlyStart,				// if OnlyChromID and >= 0 then start offset along chromosome
				int	OnlyEnd);				// if OnlyChromID and >= 0 then end offset along chromosome


// CreateSNPStats
// Create nMer distribution stats for complete genome
int
CreateSNPStats(char *pszGenomeFile,			// file containing genome assembly
				 char *pszSNPFile,			// file containing SNPs
				 char *pszStatsFile,		// where to write out stats
				 char *pszBiobedFile,		// biobed file containing regional features
			 	 int NumIncludes,			// number of biobed files containing regions to include
  				 char **pszIncludeFiles,	// biobed files containing regions to include - default is to include everything
				 int NumExcludes,			// number of biobed files containing regions to exclude
				 char **pszExcludeFiles,	// biobed files containing regions to exclude - default is to exclude nothing
				 int iWindowSize,			// sampling window size
				 int nMer,					// subseq length
				 int RegLen,				// regulatory region length - up/dn stream of 5/3' 
				 bool bMultipleFeatBits,	// if true then accept alignments in which multiple feature bits are set
				 bool bOverlapNMers,		// if true then overlap NMers when generating genome statistics
				bool bGenEmptyRows,			// true if rows with 0 counts still to be generated in output
				bool bGenEmptyWindows,		// true if windows with no counts still to be generated in output
				bool bPerChrom,				// true if stats to be generated per chromosome
				char OnlyStrand,			// any or '+','-' features only
				char *pszOnlyChrom,			// if not NULL then specifies single chromosome to process
				int	OnlyStart,				// if OnlyChromID and >= 0 then start offset along chromosome
				int	OnlyEnd);				// if OnlyChromID and >= 0 then end offset along chromosome


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

bool bMultipleFeatBits;
bool bFiltLoConfidence;
bool bInferLoose;
bool bPerChrom;
bool bGenEmptyRows;
bool bGenEmptyWindows;
bool bOverlapNMers;
int Rslt;
int Idx;
char szOutputFile[_MAX_PATH];
char szInputFile[_MAX_PATH];
char szInputBiobedFile[_MAX_PATH];
char szSNPFile[_MAX_PATH];
char szOnlyChrom[cMaxDatasetSpeciesChrom];
int NumIncludeFiles;
char *pszIncludeFiles[cMaxExcludeFiles];
int NumExcludeFiles;
char *pszExcludeFiles[cMaxIncludeFiles];
int iWindowSize;
char szSpeciesList[2048];
int iNumSpecies;
int iNMer;
int iMode;
int iRegLen;
int iStart;
int iEnd;
char cStrand;

// command line args
struct arg_lit  *help    = arg_lit0("hH","help",                "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_int *ScreenLogLevel=arg_int0("S", "ScreenLogLevel",	"<int>","Level of diagnostics written to logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",	"diagnostics log file");

struct arg_int  *Mode = arg_int0("m","mode","<int>",			"processing mode - 0 (default) = alignments, 1 = genomes, 2 = SNPs");
struct arg_int  *NMer = arg_int0("n","nmer","<int>",			"centroid oligo length to process, must be 3,5,7,9,11 or 13 (default=5)");
struct arg_int  *RegLen = arg_int0("L","updnstream","<int>",	"5'up or 3'down  stream regulatory region length (default = 2000) 0..1000000");
struct arg_file *InFile = arg_file1("i",NULL,"<file>",			"ref, rel and outgroup species alignment (.algn) file");
struct arg_file *SNPFile = arg_file0("p","snp","<file>",		"SNP input file when analysing SNPs in mode 2 processing");
struct arg_file *OutFile = arg_file1("o",NULL,"<file>",			"output to global results file as CSV");
struct arg_file *InBedFile = arg_file0("b","bed","<file>",		"characterise regions from biobed file");
struct arg_str  *SpeciesList = arg_str1("s","species","<string>","species list, 1st=ref 2nd = rel 3rd = outgroup");
struct arg_lit  *MultipleFeatBits  = arg_lit0("q","multiplefeatbits",	"single non-overlapping featbit (default) or multiple overlapping featbits allowed");
struct arg_lit  *FiltLoConfidence = arg_lit0("l","filt",		"filter out low confidence subsequences,( slough subsequence < 15mer and first/last subsequence, subsequence must start/end on identical, identity > 70)");
struct arg_str  *OnlyChrom = arg_str0("C","onlychrom","<string>","only gen stats for this chromosome");
struct arg_int  *Start = arg_int0("L","start","<int>",        "start from chromosome offset (0..n), -C <chrom> must also be specified");
struct arg_int  *End = arg_int0("E","end","<int>",            "end at chromosome offset (0..n), -C <chrom>, -O <start> must also be specified");
struct arg_file *ExcludeFile = arg_filen("X","exclude","<file>",0,cMaxExcludeFiles,	"exclude all regions in biobed file from processing ");
struct arg_file *IncludeFile = arg_filen("I","include","<file>",0,cMaxExcludeFiles,	"include all regions (unless specific regions excluded) in biobed file");
struct arg_int  *WindowSize = arg_int0("w","windowsize","<int>","if non-zero then sets fixed size window (in Knt) used to sample along genome");
struct arg_lit  *InferLoose = arg_lit0("g","inferloose",		"relax requirements on outgroup flanking bases - need not be identical");
struct arg_lit  *OverlapNMers = arg_lit0("z","overlapnmers",	"overlap NMers when generating stats in mode==1");
struct arg_lit  *PerChrom = arg_lit0("c","chromper",			"generate stats for each chromosome (default = complete genome)");
struct arg_int  *EmptyCounts = arg_int0("y","emptycount","<int>","1 == output empty rows, 2 == output empty windows, 3 == output both");
struct arg_str  *Strand = arg_str0("P","strand","<strand>",    "only process features on specified strand '+' or '-'");


struct arg_end *end = arg_end(20);

void *argtable[] = {help,version,FileLogLevel,ScreenLogLevel,LogFile,WindowSize,NMer,RegLen,InFile,SNPFile,OutFile,IncludeFile,ExcludeFile,InBedFile,SpeciesList,MultipleFeatBits,Mode,OverlapNMers,FiltLoConfidence,OnlyChrom,Strand,PerChrom,Start,End,EmptyCounts,InferLoose,end};

char **pAllArgs;
int argerrors;
argerrors = CUtility::arg_parsefromfile(argc,(char **)argv,&pAllArgs);
if(argerrors >= 0)
	argerrors = arg_parse(argerrors,pAllArgs,argtable);

    /* special case: '--help' takes precedence over error reporting */
if (help->count > 0)
        {
		printf("Usage: %s ", gszProcName);
        arg_print_syntax(stdout,argtable,"\n");
        arg_print_glossary(stdout,argtable,"  %-25s %s\n");
		exit(1);
        }

    /* special case: '--version' takes precedence error reporting */
if (version->count > 0)
        {
			printf("\n%s: Version: %d.%2.2d\n",gszProcName,cProgVer/100,cProgVer%100);
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


	bInferLoose = InferLoose->count ? true : false;
	bOverlapNMers = OverlapNMers->count ? true : false;
	bFiltLoConfidence = FiltLoConfidence->count ? true : false;
	bMultipleFeatBits = MultipleFeatBits->count ? true : false;
	bPerChrom = PerChrom->count ? true : false;
	bGenEmptyRows = false;
	bGenEmptyWindows = false;
	if(EmptyCounts->count)
		{
		bGenEmptyRows = EmptyCounts->ival[0] & 0x01 ? true : false ; 
		bGenEmptyWindows = EmptyCounts->ival[0] & 0x02 ? true : false ; 
		}
	iMode = Mode->count ? Mode->ival[0] : ePMCentroids;
	if(iMode < ePMCentroids || iMode > ePMSNPs)
		{
		printf("\nError: processing mode specified as '-m%d' must be in range of 0..2",iMode);
		exit(1);
		}
	
	if(Strand->count)
		{
		cStrand = *(Strand->sval[0]);
		if(cStrand != '*' && cStrand != '-' && cStrand != '+')
			{
			printf("\nStrand must be specified as being '*' (any) or '+','-' and not '%c'\n",cStrand);
			exit(1);
			}	
		}
	else
		cStrand = '*';
	strcpy(szInputFile,InFile->filename[0]);

	if(iMode == ePMSNPs)
		{
		if(SNPFile->count)
			strcpy(szSNPFile,SNPFile->filename[0]);
		else
			{
			printf("\nProcessing SNPs (mode === 2) but no SNP file specified!\n");
			exit(1);
			}
		}

	strcpy(szOutputFile,OutFile->filename[0]);
	if(InBedFile->count)
		strcpy(szInputBiobedFile,InBedFile->filename[0]);
	else
		szInputBiobedFile[0] = '\0';

	NumIncludeFiles = IncludeFile->count;
	for(Idx=0;Idx < IncludeFile->count; Idx++)
		pszIncludeFiles[Idx]=(char *)IncludeFile->filename[Idx];
	NumExcludeFiles = ExcludeFile->count;
	for(Idx=0;Idx < ExcludeFile->count; Idx++)
		pszExcludeFiles[Idx] = (char *)ExcludeFile->filename[Idx];

	iStart = Start->count ? Start->ival[0] : 0;
	iEnd = End->count ? End->ival[0] : 0;

	if(OnlyChrom->count)
		{
		strncpy(szOnlyChrom,OnlyChrom->sval[0],cMaxDatasetSpeciesChrom);
		szOnlyChrom[cMaxDatasetSpeciesChrom-1]='\0';
		bPerChrom = true;
		if(iStart < 0)
			iStart = 0;
		if(iEnd < iStart)
			{
			printf("\nUnable to continue, End is before Start\n");
			exit(1);
			}
		}
	else
		{
		iStart = -1;
		iEnd = -1;
		szOnlyChrom[0] = '\0';
		}

	iWindowSize = WindowSize->count ? WindowSize->ival[0] : 0;
	if(iWindowSize < 0)
		iWindowSize = 0;

	iRegLen = RegLen->count ? RegLen->ival[0] : cDfltRegLen;
	if(iRegLen < cMinRegLen)
		iRegLen = cMinRegLen;
	else
		if(iRegLen > cMaxRegLen)
			iRegLen = cMaxRegLen;

	iNMer = NMer->count ? NMer->ival[0] : cDfltNMer;
	if(iNMer != 3 && iNMer != 5 && iNMer != 7 && iNMer != 9 && iNMer != 11 && iNMer != 13)
		{
		printf("\nError: centriod oligo length specified as '-n%d' must be one of 3,5,7,9,11 or 13",iNMer);
		exit(1);
		}

	// if processing alignments then need 2 or 3 species to be specified, first is ref, 2nd is rel and 3rd is optional outgroup
	strcpy(szSpeciesList,SpeciesList->sval[0]);
	TrimQuotes(szSpeciesList);
	iNumSpecies = ParseNumSpecies(szSpeciesList,NULL);

	if(iMode == 0 && iNumSpecies < 2)
		{
		printf("Only %d species specified, need to specify at least Ref + Rel, Outgroup is optional, species in SpeciesList when processing alignments",iNumSpecies);
		exit(1);
		}
	if(iMode == 0 && iNumSpecies > 3)
		printf("%d species specified, additional species after Ref + Rel + Outgroup species in SpeciesList ignored",iNumSpecies);
	else
		if(iMode > 0 && iNumSpecies > 1)
			printf("%d species specified, additional species after initial species in SpeciesList ignored",iNumSpecies);
	
	if(iRegLen < cMinRegLen)
		iRegLen = cMinRegLen;
	else
		if(iRegLen > cMaxRegLen)
			iRegLen = cMaxRegLen;


			// now that command parameters have been parsed then initialise diagnostics log system
	if(!gDiagnostics.Open(szLogFile,(etDiagLevel)iScreenLogLevel,(etDiagLevel)iFileLogLevel,true))
		{
		printf("\nError: Unable to start diagnostics subsystem.");
		if(szLogFile[0] != '\0')
			printf(" Most likely cause is that logfile '%s' can't be opened/created",szLogFile);
		exit(1);
		}

	
	gStopWatch.Start();
	switch(iMode) {
		case ePMnMerFreq:
			// Create nMer distribution stats for complete genome
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"Version: %d.%2.2d Processing parameters:",cProgVer/100,cProgVer%100);
			gDiagnostics.DiagOutMsgOnly(eDLInfo,gszProcName,"Mode (nMer frequency distribution stats for complete genome): %d",iMode);

			gDiagnostics.DiagOutMsgOnly(eDLInfo,"file containing genome assembly: '%s'",szInputFile);
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"where to write out stats: '%s'",szOutputFile);
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"biobed file containing regional features: '%s'",szInputBiobedFile);
			for(Idx = 0; Idx < NumIncludeFiles; Idx++)
				gDiagnostics.DiagOutMsgOnly(eDLInfo,"biobed file containing regions to include: '%s'",pszIncludeFiles[Idx]); 
			for(Idx = 0; Idx < NumExcludeFiles; Idx++)
				gDiagnostics.DiagOutMsgOnly(eDLInfo,"biobed file containing regions to exclude: '%s'",pszExcludeFiles[Idx]); 
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"sampling window size: %d",iWindowSize);
			
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"subseq length: %d",iNMer);
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"regulatory region length: %d",iRegLen);
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"accept alignments in which multiple feature bits are set: %s",bMultipleFeatBits ? "yes" : "no");
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"overlap NMers when generating genome statistics: %s",bOverlapNMers ? "yes" : "no");
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"empty rows in statistics: %s",bGenEmptyRows ? "yes" : "no");
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"empty windows in statistics: %s",bGenEmptyWindows ? "yes" : "no");
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"statistics for each chromosome: %s",bPerChrom ? "yes" : "no");
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"process this strand: '%c'",cStrand);
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"process this chromosome: '%s'",szOnlyChrom);
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"start loci: %d",iStart);
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"end loci: %d",iEnd);
			Rslt = CreateGenomeStats(szInputFile,	// file containing genome assembly
					 szOutputFile,		// where to write out stats
					 szInputBiobedFile,	// biobed file containing regional features
					 NumIncludeFiles,	// number of include region files
					 pszIncludeFiles,	// biobed files containing regions to include - default is to exclude none
					 NumExcludeFiles,	// number of exclude region files
					 pszExcludeFiles,	// biobed file containing regions to exclude - default is to include all
					 iWindowSize,		// sampling window size
					 iNMer,				// subseq length
					 iRegLen,			// regulatory region length - up/dn stream of 5/3' 
					 bMultipleFeatBits,	// if true then accept alignments in which multiple feature bits are set
					 bOverlapNMers,	    // if true then overlap NMers when generating genome statistics, otherwise NMers are non-overlaping
					 bGenEmptyRows,
					 bGenEmptyWindows,
					 bPerChrom,			// generate stats for each chromosome
					 cStrand,			// any or '+','-' features only
					 szOnlyChrom,		// if not '\0' then specifies single chromosome to process
					 iStart,			// start offset along chromosome
					 iEnd);				// end offset along chromosome
			break;
		case ePMSNPs:
					// Create SNP stats
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"Version: %d.%2.2d Processing parameters:",cProgVer/100,cProgVer%100);
			gDiagnostics.DiagOutMsgOnly(eDLInfo,gszProcName,"Mode (SNP centroids): %d",iMode);
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"file containing genome assembly: '%s'",szInputFile);
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"biobed file containing SNPs: '%s'",szSNPFile);
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"where to write out stats: '%s'",szOutputFile);
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"biobed file containing regional features: '%s'",szInputBiobedFile);
			for(Idx = 0; Idx < NumIncludeFiles; Idx++)
				gDiagnostics.DiagOutMsgOnly(eDLInfo,"biobed file containing regions to include: '%s'",pszIncludeFiles[Idx]); 
			for(Idx = 0; Idx < NumExcludeFiles; Idx++)
				gDiagnostics.DiagOutMsgOnly(eDLInfo,"biobed file containing regions to exclude: '%s'",pszExcludeFiles[Idx]); 
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"sampling window size: %d",iWindowSize);
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"subseq length: %d",iNMer);
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"regulatory region length: %d",iRegLen);
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"accept alignments in which multiple feature bits are set: %s",bMultipleFeatBits ? "yes" : "no");
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"overlap NMers when generating genome statistics: %s",bOverlapNMers ? "yes" : "no");
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"empty rows in statistics: %s",bGenEmptyRows ? "yes" : "no");
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"empty windows in statistics: %s",bGenEmptyWindows ? "yes" : "no");
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"statistics for each chromosome: %s",bPerChrom ? "yes" : "no");
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"process this strand: '%c'",cStrand);
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"process this chromosome: '%s'",szOnlyChrom);
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"start loci: %d",iStart);
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"end loci: %d",iEnd);

			Rslt = CreateSNPStats(szInputFile,	// file containing genome assembly
					szSNPFile,			// biobed file containing SNPs
					szOutputFile,		// where to write out stats
					szInputBiobedFile,	// biobed file containing regional features
					NumIncludeFiles,	// number of include region files
					pszIncludeFiles,	// biobed files containing regions to include - default is to exclude none
					NumExcludeFiles,	// number of exclude region files
					pszExcludeFiles,	// biobed file containing regions to exclude - default is to include all
					iWindowSize,		// sampling window size
					iNMer,				// subseq length
					iRegLen,			// regulatory region length - up/dn stream of 5/3' 
					bMultipleFeatBits,	// if true then accept alignments in which multiple feature bits are set
					bOverlapNMers,	    // if true then overlap NMers when generating genome statistics, otherwise NMers are non-overlaping
					bGenEmptyRows,
					bGenEmptyWindows,
					bPerChrom,			// generate stats for each chromosome
					cStrand,			// any or '+','-' features only
					szOnlyChrom,		// if not NULL then specifies single chromosome to process
					iStart,				// start offset along chromosome
					iEnd);		
			break;

		case ePMCentroids:
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"Version: %d.%2.2d Processing parameters:",cProgVer/100,cProgVer%100);
				gDiagnostics.DiagOutMsgOnly(eDLInfo,gszProcName,"Mode (Alignment centroids): %d",iMode);
				gDiagnostics.DiagOutMsgOnly(eDLInfo,"file containing alignment: '%s'",szInputFile);
				gDiagnostics.DiagOutMsgOnly(eDLInfo,"where to write out stats: '%s'",szOutputFile);
				gDiagnostics.DiagOutMsgOnly(eDLInfo,"biobed file containing regional features: '%s'",szInputBiobedFile);
			for(Idx = 0; Idx < NumIncludeFiles; Idx++)
				gDiagnostics.DiagOutMsgOnly(eDLInfo,"biobed file containing regions to include: '%s'",pszIncludeFiles[Idx]); 
			for(Idx = 0; Idx < NumExcludeFiles; Idx++)
				gDiagnostics.DiagOutMsgOnly(eDLInfo,"biobed file containing regions to exclude: '%s'",pszExcludeFiles[Idx]); 
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"sampling window size: %d",iWindowSize);
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"species list: '%s'",szSpeciesList);
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"subseq length: %d",iNMer);
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"regulatory region length: %d",iRegLen);
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"accept alignments in which multiple feature bits are set: %s",bMultipleFeatBits ? "yes" : "no");

			gDiagnostics.DiagOutMsgOnly(eDLInfo,"filter out low confidence subsequences: %s",bFiltLoConfidence ? "yes" : "no");
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"treat 1st and last subsequences as being low confidence: %s",bFiltLoConfidence ? "yes" : "no");
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"subsequences of less than this length are treated as being low confidence: %d",bFiltLoConfidence ? 15 : 0);
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"treat subsequences of less than this identity as being low confidence: %d",bFiltLoConfidence? 70 : 0);
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"if refbase == relbase accept refbase as ancestor even if outgroup differs: %s",bInferLoose ? "yes" : "no");
		
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"empty rows in statistics: %s",bGenEmptyRows ? "yes" : "no");
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"empty windows in statistics: %s",bGenEmptyWindows ? "yes" : "no");
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"statistics for each chromosome: %s",bPerChrom ? "yes" : "no");
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"process this strand: '%c'",cStrand);
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"process this chromosome: '%s'",szOnlyChrom);
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"start loci: %d",iStart);
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"end loci: %d",iEnd);

			Rslt = CreateAlignStats(szInputFile,szOutputFile,szInputBiobedFile,
					NumIncludeFiles,	// number of include region files
					pszIncludeFiles,	// biobed files containing regions to include - default is to exclude none
					NumExcludeFiles,	// number of exclude region files
					pszExcludeFiles,	// biobed file containing regions to exclude - default is to include all
 					iWindowSize,		// sampling window size
					szSpeciesList,iNMer,iRegLen,bMultipleFeatBits,bFiltLoConfidence,
					bFiltLoConfidence, bFiltLoConfidence ? 15 : 0,bFiltLoConfidence? 70 : 0,
					bInferLoose,
					iMode,
  				 bGenEmptyRows,
				 bGenEmptyWindows,
 				    bPerChrom,		// generate stats for each chromosome
				 cStrand,			// any or '+','-' features only
				 szOnlyChrom,		// if not NULL then specifies single chromosome to process
				 iStart,			// start offset along chromosome
				 iEnd);				// end offset along chromosome
		}
	gStopWatch.Stop();
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Exit Code: % d Total processing time: %s",Rslt,gStopWatch.Read());
	exit(Rslt);
	}
else
	{
	arg_print_errors(stdout,end,gszProcName);
	arg_print_syntax(stdout,argtable,"\nend of help\n");
	exit(1);
	}
}


// ProcAlignBlock
// Process an alignment block which may contain aligned subsequences meeting processing requirements
// If bFiltLoConfidence then -
// A) First and last subsequences in any alignment are discarded
// B) Subsequences are trimed left and right to the first base which is not mutatated over all species 
// C) Remaining subsequences of less than 15mer are discarded
bool 
ProcAlignBlock(int RefChromID,	    // reference chromosome
			   int RefChromOfs,		// offset along reference chromosome 
			   int AlignLen,		// alignment length incl InDels
			   tsProcParams *pProcParams) // global processing parameters
{
int CurSeqLen;
int CurFiltSeqLen;	
int SubRefOfs;
int SubSeqStartIdx;
int NumIdents;
int SeqIdx;
int NumSubSeqs;
int CurNumSeqs;
int VIdx;
bool bInDel = false;
bool bAllIdentical; 
etSeqBase *pSeq;
etSeqBase RefBase;
etSeqBase SeqBase;
SubRefOfs = RefChromOfs;
bInDel = false;
if(pProcParams->bFiltLoConfidence && pProcParams->bFilt1stLast)
	NumSubSeqs = CUtility::GetNumSubseqs(AlignLen,pProcParams->NumSpecies2Align,pProcParams->pSeqs);
else
	NumSubSeqs = 0;
CurNumSeqs=0;
CurSeqLen = 0;
for(SeqIdx = 0; SeqIdx < AlignLen; SeqIdx++)
	{
	bAllIdentical = true;			// assume all identical until at least one mismatch
	for(VIdx = 0; VIdx < pProcParams->NumSpecies2Align; VIdx++)
		{
		pSeq = pProcParams->pSeqs[VIdx];
		SeqBase = pSeq[SeqIdx] & ~cRptMskFlg;
		if(SeqBase > eBaseT && (VIdx < 2 || !pProcParams->bInferLoose))		 // treat anything not a,c,g or t as an InDel
			bInDel = true;
		if(VIdx == pProcParams->RefSpeciesIdx)
			RefBase = SeqBase;
		else
			if(SeqBase != RefBase && (VIdx < 2 || !pProcParams->bInferLoose)) // is it a mismatch???
				bAllIdentical = false;
		}

	if(bInDel)						// if a InDel then check if subsequence to process
		{
		if(CurSeqLen)
			{
			CurNumSeqs+=1;
			if(CurFiltSeqLen > pProcParams->MinSubSeqLen)			
				{
				if(pProcParams->bFiltLoConfidence)					// when filtering simply slough any low confidence subsequences
					{
					if(pProcParams->bFilt1stLast && (CurNumSeqs==1 || CurNumSeqs == NumSubSeqs))	// don't process first and last subsequence
						CurFiltSeqLen = 0;
					if(CurFiltSeqLen &&  pProcParams->MinIdent > 0 && ((NumIdents * 100)/CurFiltSeqLen) < pProcParams->MinIdent)
						CurFiltSeqLen = 0;
					}
				if(CurFiltSeqLen > pProcParams->NMer)
					ProcessAlignment(RefChromID,SubRefOfs,CurFiltSeqLen,SubSeqStartIdx,pProcParams);
				}
			}
		CurSeqLen = 0;
		CurFiltSeqLen = 0;
		bInDel = false;
		if(RefBase != eBaseInDel)
			RefChromOfs++;
		continue;
		}

	// no InDels or 'N's in any of the required alignments, must be part of a subsequence
	if(CurSeqLen == 0)		// if first base of putative subsequence...
		{
		if(pProcParams->bFiltLoConfidence && !bAllIdentical) // when filtering, must start on an identical base
			continue;
		NumIdents = 0;
		SubSeqStartIdx = SeqIdx;			// mark where subsequence starts 
		SubRefOfs = RefChromOfs; 
		}
	CurSeqLen++;
	if(bAllIdentical)						// mark last identical...
		{
		NumIdents++;
		CurFiltSeqLen = CurSeqLen;
		}
	else
		CurFiltSeqLen = CurSeqLen;
	}

if(CurFiltSeqLen > pProcParams->MinSubSeqLen)			
	{
	if(pProcParams->bFiltLoConfidence) // when filtering simply slough any low confidence subsequences
		{
		if(pProcParams->bFilt1stLast && (!CurNumSeqs || CurNumSeqs == (NumSubSeqs-1)))	// don't process first and last subsequence
			CurFiltSeqLen = 0;
		if(CurFiltSeqLen &&  pProcParams->MinIdent > 0 && ((NumIdents * 100)/CurFiltSeqLen) < pProcParams->MinIdent)
			CurFiltSeqLen = 0;
		}
	if(CurFiltSeqLen > pProcParams->NMer)
			ProcessAlignment(RefChromID,SubRefOfs,CurFiltSeqLen,SubSeqStartIdx,pProcParams);
	}
return(true);
}

int
ProcessGenome(char *pszGenomeFile,tsProcParams *pProcParams)
{
int ReqStart;
int ReqLen;
int Rslt;
int CurEntryID;
int CurDataLen;
CBioSeqFile *pBioSeq;
pBioSeq = new CBioSeqFile;
if((Rslt=pBioSeq->Open(pszGenomeFile))<0)
	{
	while(pBioSeq->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,pBioSeq->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open bioseq file '%s' for processing",pszGenomeFile);
	return(Rslt);
	}

if((pProcParams->pSeq = new unsigned char [cChromSeqLen])==NULL)
	{
	pBioSeq->Close();
	delete pBioSeq;
	return(eBSFerrMem);
	}
pProcParams->AllocdDataLen = cChromSeqLen;

CurEntryID = 0;
pProcParams->NxtOutputOffset = pProcParams->iWindowSize;

// if only a single chromosome to be processed then ensure chromosme is represented in genome sequence file
if(pProcParams->szOnlyChrom[0]!='\0')
	{
	pProcParams->OnlyChromID=pBioSeq->LocateEntryIDbyName(pProcParams->szOnlyChrom); //returns chromosome identifier
	if(pProcParams->OnlyChromID < 1)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Chromosome '%s' specified as only chromosome to process is not present for species '%s' in bioseq file '%s'",
				pProcParams->szOnlyChrom,pProcParams->szSpecies[0],pszGenomeFile);
		delete pBioSeq;
		return(eBSFerrChrom);
		}
	}

while((CurEntryID = pBioSeq->Next(CurEntryID))>0)
	{
	if(pProcParams->OnlyChromID > 0 && CurEntryID != pProcParams->OnlyChromID)
		continue;
	CurDataLen = pBioSeq->GetDataLen(CurEntryID);
	if(pProcParams->OnlyChromID > 0 && CurEntryID == pProcParams->OnlyChromID)
		{
		if(pProcParams->OnlyEnd <= 0)
			pProcParams->OnlyEnd = CurDataLen-1;
		if(pProcParams->OnlyStart >= CurDataLen || pProcParams->OnlyStart > pProcParams->OnlyEnd)
			{
			gDiagnostics.DiagOut(eDLWarn,gszProcName,"Requested start/end (%d/%d) past end (%d) of chromosome %s",pProcParams->OnlyStart,pProcParams->OnlyEnd,CurDataLen,pProcParams->szOnlyChrom);
			break;
			}
		ReqStart = pProcParams->OnlyStart;
		ReqLen = pProcParams->OnlyEnd - pProcParams->OnlyStart + 1;
		}
	else
		{
		ReqStart = 0;
		ReqLen = CurDataLen;
		}

	if(ReqLen > pProcParams->AllocdDataLen)
		{
		delete pProcParams->pSeq;
		if((pProcParams->pSeq = new unsigned char [ReqLen+0x07fff])==NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to alloc memory (%d requested) for holding raw sequence data",ReqLen+0x07fff);
			Rslt = eBSFerrMem;
			break;
			}
		pProcParams->AllocdDataLen = ReqLen+0x07fff;
		}

	if((Rslt=pBioSeq->GetData(CurEntryID,eSeqBaseType,ReqStart,pProcParams->pSeq,ReqLen))<eBSFSuccess)
		continue;
	pBioSeq->GetName(CurEntryID,sizeof(pProcParams->szRefChrom),pProcParams->szRefChrom);
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing %s...",pProcParams->szRefChrom);
	if(pProcParams->pBiobed != NULL)
		{
		pProcParams->BEDChromID = pProcParams->pBiobed->LocateChromIDbyName(pProcParams->szRefChrom);
		if(pProcParams->BEDChromID < 1)
			{
			gDiagnostics.DiagOut(eDLWarn,gszProcName,"\nUnable to locate chromosome '%s' in biobed file, unable to generate regional stats...\n",pProcParams->szRefChrom);
			pProcParams->BEDChromID = 0;
			}
		}
	else
		pProcParams->BEDChromID = -1;
	pProcParams->SeqLen = ReqLen;
	ProcessSequence(CurEntryID,ReqStart,pProcParams);
	if(pProcParams->iWindowSize > 0)
		ChkOutputGenomeResults(pProcParams->szRefChrom,-1,pProcParams);	// output results as appropriate
	else
		if(pProcParams->bPerChrom)
			OutputGenomeResults(pProcParams->szRefChrom,ReqStart,pProcParams);
	}

if(pProcParams->iWindowSize == 0 && !pProcParams->bPerChrom)
	OutputGenomeResults((char *)"Genome",0,pProcParams);
if(pBioSeq != NULL)
	{
	pBioSeq->Close();
	delete pBioSeq;
	}
if(pProcParams->pSeq != NULL)
	{
	delete pProcParams->pSeq;
	pProcParams->pSeq = NULL;
	}
pProcParams->AllocdDataLen = 0;
return(Rslt);
}



int
ProcessSNPs(char *pszGenomeFile,char *pszSNPFile,tsProcParams *pProcParams)
{
int ReqStart;
int ReqLen;
int Rslt;
int CurEntryID;
int CurDataLen;
CBioSeqFile *pBioSeq;

if((pProcParams->pSNPs = new CSNPsFile)==NULL)
	return(eBSFerrObj);

if((Rslt=pProcParams->pSNPs->ProcessSNPsFile(pszSNPFile))<0)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open SNPs raw table file '%s' for processing",pszGenomeFile);
	delete pProcParams->pSNPs;
	pProcParams->pSNPs = NULL;
	return(Rslt);
	}

if((pBioSeq = new CBioSeqFile)==NULL)
	return(eBSFerrObj);

if((Rslt=pBioSeq->Open(pszGenomeFile))<0)
	{
	while(pBioSeq->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,pBioSeq->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open bioseq file '%s' for processing",pszGenomeFile);
	delete pProcParams->pSNPs;
	pProcParams->pSNPs = NULL;
	return(Rslt);
	}

if((pProcParams->pSeq = new unsigned char [cChromSeqLen])==NULL)
	{
	pBioSeq->Close();
	delete pBioSeq;
	delete pProcParams->pSNPs;
	pProcParams->pSNPs = NULL;
	return(eBSFerrMem);
	}
pProcParams->AllocdDataLen = cChromSeqLen;

CurEntryID = 0;
pProcParams->NxtOutputOffset = pProcParams->iWindowSize;

// if only a single chromosome to be processed then ensure chromosme is represented in
// both the alignment file and also that there is at least one SNP on that chromosome
if(pProcParams->szOnlyChrom[0]!='\0')
	{
	pProcParams->OnlyChromID=pBioSeq->LocateEntryIDbyName(pProcParams->szOnlyChrom); //returns chromosome identifier
	if(pProcParams->OnlyChromID < 1)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Chromosome '%s' specified as only chromosome to process is not present for species '%s' in bioseq file '%s'",
				pProcParams->szOnlyChrom,pProcParams->szSpecies[0],pszGenomeFile);
		delete pBioSeq;
		delete pProcParams->pSNPs;
		pProcParams->pSNPs = NULL;
		return(eBSFerrChrom);
		}

	if(pProcParams->pSNPs->LocateChrom(pProcParams->szOnlyChrom) < 1)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Chromosome '%s' specified as only chromosome to process but no SNPs on this chromosome\n",
				pProcParams->szOnlyChrom);

		delete pBioSeq;
		delete pProcParams->pSNPs;
		pProcParams->pSNPs = NULL;
		return(eBSFerrChrom);
		}
	}

while((CurEntryID = pBioSeq->Next(CurEntryID))>0)
	{
	if(pProcParams->OnlyChromID > 0 && CurEntryID != pProcParams->OnlyChromID)
		continue;
	CurDataLen = pBioSeq->GetDataLen(CurEntryID);
	if(pProcParams->OnlyChromID > 0 && CurEntryID == pProcParams->OnlyChromID)
		{
		if(pProcParams->OnlyEnd <= 0)
			pProcParams->OnlyEnd = CurDataLen-1;
		if(pProcParams->OnlyStart >= CurDataLen || pProcParams->OnlyStart > pProcParams->OnlyEnd)
			{
			gDiagnostics.DiagOut(eDLWarn,gszProcName,"\nRequested start/end (%d/%d) past end (%d) of chromosome %s",pProcParams->OnlyStart,pProcParams->OnlyEnd,CurDataLen,pProcParams->szOnlyChrom);
			break;
			}
		ReqStart = pProcParams->OnlyStart;
		ReqLen = pProcParams->OnlyEnd - pProcParams->OnlyStart + 1;
		}
	else
		{
		ReqStart = 0;
		ReqLen = CurDataLen;
		}

	if(ReqLen > pProcParams->AllocdDataLen)
		{
		delete pProcParams->pSeq;
		if((pProcParams->pSeq = new unsigned char [ReqLen+0x07fff])==NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to alloc memory (%d requested) for holding raw sequence data\n",ReqLen+0x07fff);
			Rslt = eBSFerrMem;
			break;
			}
		pProcParams->AllocdDataLen = ReqLen+0x07fff;
		}

	pBioSeq->GetName(CurEntryID,sizeof(pProcParams->szRefChrom),pProcParams->szRefChrom);
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing %s...",pProcParams->szRefChrom);
	if(pProcParams->pSNPs->LocateChrom(pProcParams->szRefChrom) < 1)
		{
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"no SNPs on chromosome '%s'\n",
				pProcParams->szRefChrom);
		continue;
		}

	if((Rslt=pBioSeq->GetData(CurEntryID,eSeqBaseType,ReqStart,pProcParams->pSeq,ReqLen))<eBSFSuccess)
		continue;
	if(pProcParams->pBiobed != NULL)
		{
		pProcParams->BEDChromID = pProcParams->pBiobed->LocateChromIDbyName(pProcParams->szRefChrom);
		if(pProcParams->BEDChromID < 1)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to locate chromosome '%s' in biobed file, unable to generate regional stats...\n",pProcParams->szRefChrom);
			pProcParams->BEDChromID = 0;
			}
		}
	else
		pProcParams->BEDChromID = -1;
	pProcParams->SeqLen = ReqLen;
	if(pProcParams->iWindowSize)
		ChkOutputResults(pProcParams->szRefChrom, ReqStart,pProcParams);
	ProcessSNPSequence(CurEntryID,ReqStart,pProcParams);
	if(pProcParams->iWindowSize > 0)
		ChkOutputGenomeResults(pProcParams->szRefChrom,-1,pProcParams);	// output results as appropriate
	else
		if(pProcParams->bPerChrom)
			OutputResults(pProcParams->szRefChrom,ReqStart,pProcParams);
	}
if(pProcParams->iWindowSize)
	ChkOutputResults(pProcParams->szRefChrom, -1,pProcParams);
else
	{
	if(pProcParams->bPerChrom)
		OutputResults(pProcParams->szRefChrom,0,pProcParams);
	else
		OutputResults((char *)"Genome",0,pProcParams);
	}
if(pBioSeq != NULL)
	{
	pBioSeq->Close();
	delete pBioSeq;
	}

if(pProcParams->pSNPs!=NULL)
	{
	delete pProcParams->pSNPs;
	pProcParams->pSNPs = NULL;
	}

if(pProcParams->pSeq != NULL)
	{
	delete pProcParams->pSeq;
	pProcParams->pSeq = NULL;
	}
pProcParams->AllocdDataLen = 0;
return(Rslt);
}



// OpenBedfile
// Attempts to open specified bedfile
// Returns ptr to opened bedfile or NULL
CBEDfile *
OpenBedfile(char *pToOpen,char OnStrand)
{
int Rslt;
CBEDfile *pBed;
if(pToOpen != NULL && pToOpen[0] != '\0')
	{
	if((pBed = (CBEDfile *)new CBEDfile())==NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CBEDfile!");
		return(NULL);
		}

	if((Rslt = pBed->Open(pToOpen,eBTAnyBed))!=eBSFSuccess)
		{
		while(pBed->NumErrMsgs())
			gDiagnostics.DiagOut(eDLFatal,gszProcName,pBed->GetErrMsg());
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open biobed file %s",pToOpen);
		delete pBed;
		return(NULL);
		}
	pBed->SetStrand(OnStrand);
	return(pBed);
	}
return(NULL);
}

//CloseBedfiles
//Closes and deletes all created and opened Biobed files
bool
CloseBedfiles(tsProcParams *pProcParams)
{
int Idx;
for(Idx=0; Idx < pProcParams->NumIncludes; Idx++)
	{
	if(pProcParams->pIncludes[Idx] != NULL)
		delete pProcParams->pIncludes[Idx];
	pProcParams->pIncludes[Idx] = NULL;
	}
for(Idx=0; Idx < pProcParams->NumExcludes; Idx++)
	{
	if(pProcParams->pExcludes[Idx] != NULL)
		delete pProcParams->pExcludes[Idx];
	pProcParams->pExcludes[Idx] = NULL;
	}
if(pProcParams->pBiobed != NULL)
	delete pProcParams->pBiobed;
pProcParams->pBiobed = NULL;
return(true);
}

// CreateAlignStats
// Create alignment stats
int
CreateAlignStats(char *pszRefRelFile,		// file containing ref + rel + outgroup species
				 char *pszStatsFile,		// where to write out stats
				 char *pszBiobedFile,		// biobed file containing regional features
				 int NumIncludeFiles,		// number of include files
  				 char **ppszIncludeFiles,	// biobed files containing regions to exclude - default is to exclude none
				 int NumExcludeFiles,		// number of exclude files
				 char **ppszExcludeFiles,	// biobed files containing regions to include - default is to include all
				 int  iWindowSize,			// sampling window size
				 char *pszSpeciesList,		// space or comma separated list of species, ref followed by rel followed by outgroup
				 int iNMer,					// oligo length to process
				 int RegLen,				// regulatory region length - up/dn stream of 5/3' 
				 bool MultipleFeatBits,		// if true then accept alignments in which multiple feature bits are set
				 bool bFiltLoConfidence,	// true if to filter out low confidence subsequences
				 bool bFilt1stLast,			// true if to treat 1st and last subsequences as being low confidence
			 	int MinSubSeqLen,			// subsequences of less than this length are treated as being low confidence
			 	int MinIdent,				// treat subsequences of less than this identity as being low confidence
				bool bInferLoose,			// if true then if refbase == relbase accept refbase as ancestor even if outgroup differs
 			    int Mode,					// processing mode
				bool bGenEmptyRows,			// true if rows with 0 counts still to be generated in output
				bool bGenEmptyWindows,		// true if windows with no counts still to be generated in output
				bool bPerChrom,				// generate stats for each chromosome
				char OnStrand,				// stats on which strand? '*' any or '+','-' features only
				char *pszOnlyChrom,			// if not NULL then specifies single chromosome to process
				int	OnlyStart,				// if OnlyChromID and >= 0 then start offset along chromosome
				int	OnlyEnd)				// if OnlyChromID and >= 0 then end offset along chromosome
{
int Rslt;
int hRsltsFile;
int Pwr;

int Idx;
int *pCntStepCnts;
int NumCnts;
int Regions;
int UpDnStreamLen;
tsProcParams ProcParams;
int NumSpecies;

memset(&ProcParams,0,sizeof(tsProcParams));
// parse out species list - interested in 1st 3 species: ref,rel and optional outgroup
NumSpecies = ParseNumSpecies(pszSpeciesList,&ProcParams);
if(NumSpecies < 2)
	return(eBSFerrParams);

for(Idx=0;Idx<NumIncludeFiles; Idx++)
	{
	if((ProcParams.pIncludes[Idx] = OpenBedfile(ppszIncludeFiles[Idx],OnStrand))==NULL)
		{
		CloseBedfiles(&ProcParams);
		return(eBSFerrObj);
		}
	ProcParams.NumIncludes++;
	}

for(Idx=0;Idx<NumExcludeFiles; Idx++)
	{
	if((ProcParams.pExcludes[Idx] = OpenBedfile(ppszExcludeFiles[Idx],OnStrand))==NULL)
		{
		CloseBedfiles(&ProcParams);
		return(eBSFerrObj);
		}
	ProcParams.NumExcludes++;
	}

if(pszBiobedFile != NULL && pszBiobedFile[0] != '\0')
	{
	if((ProcParams.pBiobed = OpenBedfile(pszBiobedFile,OnStrand))==NULL)
		{
		CloseBedfiles(&ProcParams);
		return(eBSFerrObj);
		}
	Regions = 7;	// intergenic,5'upstream,5'utr,CDS,introns,3'utr,3'dnstream
	UpDnStreamLen = RegLen;
	}
else
	{
	Regions = 1;
	UpDnStreamLen = 0;
	}

NumCnts= 4;
for(Pwr=1;Pwr < iNMer; Pwr++)	// NumCnts = 4^NMer as central step is what we are deriving stats on
	NumCnts *= 4;

if((pCntStepCnts = new int[NumCnts * cRegionLen * Regions])==NULL)
	{
	CloseBedfiles(&ProcParams);
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate memory for holding alignment statistics\n");
	return(eBSFerrMem);
	}
memset(pCntStepCnts,0,NumCnts * cRegionLen  * Regions * sizeof(int));

ProcParams.bFiltLoConfidence = bFiltLoConfidence;
ProcParams.bFilt1stLast = bFilt1stLast;
ProcParams.MinIdent = MinIdent;
ProcParams.MinSubSeqLen = MinSubSeqLen;
ProcParams.NumSpecies2Align = NumSpecies;
ProcParams.NumSpeciesAligned = 0;
ProcParams.RefSpeciesIdx = 0;
ProcParams.bInferLoose = bInferLoose;
ProcParams.bGenEmptyRows=bGenEmptyRows;				// true if rows with 0 counts still to be generated in output
ProcParams.bGenEmptyWindows=bGenEmptyRows;			// true if windows with no counts still to be generated in output
ProcParams.bPerChrom = bPerChrom;
ProcParams.OnStrand = OnStrand;

for(Idx = 0; Idx < NumSpecies; Idx++)
	{
	if((ProcParams.pSeqs[Idx] = new unsigned char [cMALineSize])==NULL)
		{
		delete pCntStepCnts;
		CloseBedfiles(&ProcParams);
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate memory for holding species sequences");
		return(eBSFerrMem);
		}
	}

ProcParams.NumCnts = NumCnts;
ProcParams.Regions = Regions;
ProcParams.NMer = iNMer;
if(iNMer < 3)
	ProcParams.FlankLen = 0;
else
	ProcParams.FlankLen = (iNMer-1)/2;
ProcParams.pCntStepCnts = pCntStepCnts;
ProcParams.UpDnStreamLen = UpDnStreamLen;
ProcParams.iWindowSize = iWindowSize;
ProcParams.bMultipleFeatBits = MultipleFeatBits;
if(pszOnlyChrom == NULL || pszOnlyChrom[0] == '\0')
	ProcParams.szOnlyChrom[0] = '\0';
else
	strcpy(ProcParams.szOnlyChrom,pszOnlyChrom);
ProcParams.OnlyEnd = OnlyEnd;
ProcParams.OnlyStart = OnlyStart;

#ifdef _WIN32
if((hRsltsFile = open(pszStatsFile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
if((hRsltsFile = open(pszStatsFile,O_RDWR | O_CREAT |O_TRUNC, S_IREAD | S_IWRITE ))==-1)
#endif
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create %s - %s",pszStatsFile,strerror(errno));
	CloseBedfiles(&ProcParams);
	return(eBSFerrCreateFile);
	}

ProcParams.hRsltsFile = hRsltsFile;
Rslt = ProcessAlignments(pszRefRelFile,&ProcParams);
if(pCntStepCnts != NULL)
	delete(pCntStepCnts);
if(hRsltsFile != -1)
	close(hRsltsFile);
CloseBedfiles(&ProcParams);

for(Idx = 0; Idx < NumSpecies; Idx++)
	{
	if(ProcParams.pSeqs[Idx] != NULL)
		delete ProcParams.pSeqs[Idx];
	}
return(Rslt);
}



// CreateGenomeStats
// Create nMer distribution stats for complete genome
int
CreateGenomeStats(char *pszGenomeFile,		// file containing genome assembly
				 char *pszStatsFile,		// where to write out stats
				 char *pszBiobedFile,		// biobed file containing regional features
				 int NumIncludeFiles,		// number of include files
  				 char **ppszIncludeFiles,	// biobed files containing regions to exclude - default is to exclude none
				 int NumExcludeFiles,		// number of exclude files
				 char **ppszExcludeFiles,	// biobed files containing regions to include - default is to include all
 				 int  iWindowSize,			// sampling window size
				 int iNMer,					// subseq length
				 int RegLen,				// regulatory region length - up/dn stream of 5/3' 
				 bool bMultipleFeatBits,	// if true then accept alignments in which multiple feature bits are set
				 bool bOverlapNMers,	    // if true then overlap NMers when generating genome statistics
				bool bGenEmptyRows,			// true if rows with 0 counts still to be generated in output
				bool bGenEmptyWindows,		// true if windows with no counts still to be generated in output
				bool bPerChrom,				// generate stats for each chromosome
				char OnStrand,				// stats on which strand? '*' any or '+','-' features only
				char *pszOnlyChrom,			// if not NULL then specifies single chromosome to process
				int	OnlyStart,				// if OnlyChromID and >= 0 then start offset along chromosome
				int	OnlyEnd)				// if OnlyChromID and >= 0 then end offset along chromosome
{
int Rslt;
int hRsltsFile;
int Pwr;
int Idx;
int *pCntStepCnts;
int NumCnts;
int Regions;
int UpDnStreamLen;
CBEDfile *pBiobed = NULL;
CBEDfile *pInclude = NULL;
CBEDfile *pExclude = NULL;

tsProcParams ProcParams;

memset(&ProcParams,0,sizeof(tsProcParams));

for(Idx=0;Idx<NumIncludeFiles; Idx++)
	{
	if((ProcParams.pIncludes[Idx] = OpenBedfile(ppszIncludeFiles[Idx],OnStrand))==NULL)
		{
		CloseBedfiles(&ProcParams);
		return(eBSFerrObj);
		}
	ProcParams.NumIncludes++;
	}

for(Idx=0;Idx<NumExcludeFiles; Idx++)
	{
	if((ProcParams.pExcludes[Idx] = OpenBedfile(ppszExcludeFiles[Idx],OnStrand))==NULL)
		{
		CloseBedfiles(&ProcParams);
		return(eBSFerrObj);
		}
	ProcParams.NumExcludes++;
	}

if(pszBiobedFile != NULL && pszBiobedFile[0] != '\0')
	{
	if((ProcParams.pBiobed = OpenBedfile(pszBiobedFile,OnStrand))==NULL)
		{
		CloseBedfiles(&ProcParams);
		return(eBSFerrObj);
		}

	Regions = 7;	// intergenic,5'upstream,5'utr,CDS,introns,3'utr,3'dnstream
	UpDnStreamLen = RegLen;
	}
else
	{
	pBiobed = NULL;
	Regions = 1;
	UpDnStreamLen = 0;
	}

NumCnts= 4;
for(Pwr=1;Pwr < iNMer; Pwr++)	// NumCnts = 4^NMer as central step is what we are deriving stats on
	NumCnts *= 4;

if((pCntStepCnts = new int[NumCnts * Regions])==NULL)
	{
	CloseBedfiles(&ProcParams);
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate memory for holding alignment statistics\n");
	return(eBSFerrMem);
	}
memset(pCntStepCnts,0,NumCnts * Regions * sizeof(int));
ProcParams.NumCnts = NumCnts;
ProcParams.Regions = Regions;
ProcParams.NMer = iNMer;
if(iNMer < 3)
	ProcParams.FlankLen = 0;
else
	ProcParams.FlankLen = (iNMer-1)/2;
ProcParams.MerIdxMsk = 0x03;
for(Idx = 1; Idx < iNMer; Idx++)
	{
	ProcParams.MerIdxMsk <<= 2;
	ProcParams.MerIdxMsk |= 0x03;
	}
ProcParams.pCntStepCnts = pCntStepCnts;
ProcParams.UpDnStreamLen = UpDnStreamLen;
ProcParams.bMultipleFeatBits = bMultipleFeatBits;
ProcParams.bOverlapNMers = bOverlapNMers;
ProcParams.bGenEmptyRows=bGenEmptyRows;				// true if rows with 0 counts still to be generated in output
ProcParams.bGenEmptyWindows=bGenEmptyRows;			// true if windows with no counts still to be generated in output
ProcParams.OnStrand = OnStrand;
ProcParams.bPerChrom = bPerChrom;
ProcParams.iWindowSize = iWindowSize;
if(pszOnlyChrom == NULL || pszOnlyChrom[0] == '\0')
	ProcParams.szOnlyChrom[0] = '\0';
else
	strcpy(ProcParams.szOnlyChrom,pszOnlyChrom);
ProcParams.OnlyEnd = OnlyEnd;
ProcParams.OnlyStart = OnlyStart;

#ifdef _WIN32
if((hRsltsFile = open(pszStatsFile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
if((hRsltsFile = open(pszStatsFile,O_RDWR | O_CREAT |O_TRUNC, S_IREAD | S_IWRITE ))==-1)
#endif
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create %s - %s",pszStatsFile,strerror(errno));
	CloseBedfiles(&ProcParams);
	return(eBSFerrCreateFile);
	}


ProcParams.hRsltsFile = hRsltsFile;
Rslt = ProcessGenome(pszGenomeFile,&ProcParams);
if(pCntStepCnts != NULL)
	delete(pCntStepCnts);
if(hRsltsFile != -1)
	close(hRsltsFile);
CloseBedfiles(&ProcParams);
if(ProcParams.pSeq != NULL)
	delete ProcParams.pSeq;
return(Rslt);
}



// CreateSNPStats
// Create nMer distribution stats for complete genome
int
CreateSNPStats(char *pszGenomeFile,		// file containing genome assembly
			     char *pszSNPFile,			// file containing SNPs
				 char *pszStatsFile,		// where to write out stats
				 char *pszBiobedFile,		// biobed file containing regional features
				 int NumIncludeFiles,		// number of include files
  				 char **ppszIncludeFiles,	// biobed files containing regions to exclude - default is to exclude none
				 int NumExcludeFiles,		// number of exclude files
				 char **ppszExcludeFiles,	// biobed files containing regions to include - default is to include all
 				 int  iWindowSize,			// sampling window size
				 int iNMer,					// subseq length
				 int RegLen,				// regulatory region length - up/dn stream of 5/3' 
				 bool bMultipleFeatBits,	// if true then accept alignments in which multiple feature bits are set
				 bool bOverlapNMers,	    // if true then overlap NMers when generating genome statistics
				bool bGenEmptyRows,			// true if rows with 0 counts still to be generated in output
				bool bGenEmptyWindows,		// true if windows with no counts still to be generated in output
				bool bPerChrom,				// generate stats for each chromosome
				char OnStrand,				// stats on which strand? '*' any or '+','-' features only
				char *pszOnlyChrom,			// if not NULL then specifies single chromosome to process
				int	OnlyStart,				// if OnlyChromID and >= 0 then start offset along chromosome
				int	OnlyEnd)				// if OnlyChromID and >= 0 then end offset along chromosome
{
int Rslt;
int hRsltsFile;
int Pwr;
int Idx;
int *pCntStepCnts;
int NumCnts;
int Regions;
int UpDnStreamLen;
CBEDfile *pBiobed = NULL;
CBEDfile *pInclude = NULL;
CBEDfile *pExclude = NULL;

tsProcParams ProcParams;

memset(&ProcParams,0,sizeof(tsProcParams));

for(Idx=0;Idx<NumIncludeFiles; Idx++)
	{
	if((ProcParams.pIncludes[Idx] = OpenBedfile(ppszIncludeFiles[Idx],OnStrand))==NULL)
		{
		CloseBedfiles(&ProcParams);
		return(eBSFerrObj);
		}
	ProcParams.NumIncludes++;
	}

for(Idx=0;Idx<NumExcludeFiles; Idx++)
	{
	if((ProcParams.pExcludes[Idx] = OpenBedfile(ppszExcludeFiles[Idx],OnStrand))==NULL)
		{
		CloseBedfiles(&ProcParams);
		return(eBSFerrObj);
		}
	ProcParams.NumExcludes++;
	}

if(pszBiobedFile != NULL && pszBiobedFile[0] != '\0')
	{
	if((ProcParams.pBiobed = OpenBedfile(pszBiobedFile,OnStrand))==NULL)
		{
		CloseBedfiles(&ProcParams);
		return(eBSFerrObj);
		}

	Regions = 7;	// intergenic,5'upstream,5'utr,CDS,introns,3'utr,3'dnstream
	UpDnStreamLen = RegLen;
	}
else
	{
	pBiobed = NULL;
	Regions = 1;
	UpDnStreamLen = 0;
	}

NumCnts= 4;
for(Pwr=1;Pwr < iNMer; Pwr++)	// NumCnts = 4^NMer as central step is what we are deriving stats on
	NumCnts *= 4;

if((pCntStepCnts = new int[NumCnts * cRegionLen * Regions])==NULL)
	{
	CloseBedfiles(&ProcParams);
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate memory for holding alignment statistics\n");
	return(eBSFerrMem);
	}
memset(pCntStepCnts,0,NumCnts * cRegionLen  * Regions * sizeof(int));
ProcParams.NumCnts = NumCnts;
ProcParams.Regions = Regions;
ProcParams.NMer = iNMer;
if(iNMer < 3)
	ProcParams.FlankLen = 0;
else
	ProcParams.FlankLen = (iNMer-1)/2;
ProcParams.MerIdxMsk = 0x03;
for(Idx = 1; Idx < iNMer; Idx++)
	{
	ProcParams.MerIdxMsk <<= 2;
	ProcParams.MerIdxMsk |= 0x03;
	}
ProcParams.pCntStepCnts = pCntStepCnts;
ProcParams.UpDnStreamLen = UpDnStreamLen;
ProcParams.bMultipleFeatBits = bMultipleFeatBits;
ProcParams.bOverlapNMers = bOverlapNMers;
ProcParams.bGenEmptyRows=bGenEmptyRows;				// true if rows with 0 counts still to be generated in output
ProcParams.bGenEmptyWindows=bGenEmptyRows;			// true if windows with no counts still to be generated in output
ProcParams.OnStrand = OnStrand;
ProcParams.bPerChrom = bPerChrom;
ProcParams.iWindowSize = iWindowSize;
if(pszOnlyChrom == NULL || pszOnlyChrom[0] == '\0')
	ProcParams.szOnlyChrom[0] = '\0';
else
	strcpy(ProcParams.szOnlyChrom,pszOnlyChrom);
ProcParams.OnlyEnd = OnlyEnd;
ProcParams.OnlyStart = OnlyStart;

#ifdef _WIN32
if((hRsltsFile = open(pszStatsFile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
if((hRsltsFile = open(pszStatsFile,O_RDWR | O_CREAT |O_TRUNC, S_IREAD | S_IWRITE ))==-1)
#endif
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create %s - %s",pszStatsFile,strerror(errno));
	CloseBedfiles(&ProcParams);
	return(eBSFerrCreateFile);
	}


ProcParams.hRsltsFile = hRsltsFile;
Rslt = ProcessSNPs(pszGenomeFile,pszSNPFile,&ProcParams);
if(pCntStepCnts != NULL)
	delete(pCntStepCnts);
if(hRsltsFile != -1)
	close(hRsltsFile);
CloseBedfiles(&ProcParams);
if(ProcParams.pSeq != NULL)
	delete ProcParams.pSeq;
return(Rslt);
}


char *
Idx2Bases(int StepIdx,int NumBases,char *pszBuffer)
{
char *pDst = pszBuffer;
int CurBase;
int BasMsk = 0x03 << ((NumBases - 1) * 2);
while(NumBases--) {
	CurBase = (StepIdx & BasMsk) >> (NumBases * 2);
	BasMsk >>= 2;
	switch(CurBase) {
		case 0:
			*pDst++ = 'a';
			break;
		case 1:
			*pDst++ = 'c';
			break;
		case 2:
			*pDst++ = 'g';
			break;
		case 3:
			*pDst++ = 't';
			break;
		}
	}
*pDst = '\0';
return(pszBuffer);
}

// IncludeFilter
// Returns true if subsequence can be processed or false if subsequence is not to be processed
// To be processed a subsequence -
// A) If no Include BED files specified then subsequence is assumed included unless excluded by following rule B).
//    If at least one Include BED file was specified then if subsequence not in any of the Include files then false is returned
// B) If no Exclude files specified then subsequence is returned as Ok (true) to process. If subsequence not in any of the Exclude files
//    then subsequence is returned as Ok (true) to process, otherwise false is returned
bool
IncludeFilter(int RefChromID,int SubRefOfs,int SubRefEndOfs,tsProcParams *pProcParams)
{
int Idx;
int BEDChromID;
if(pProcParams->NumIncludes)
	{
	for(Idx = 0; Idx < pProcParams->NumIncludes; Idx++)
		{
		if(pProcParams->pIncludes[Idx] == NULL) // should'nt ever be NULL but...
			continue;
		if((BEDChromID = pProcParams->pIncludes[Idx]->LocateChromIDbyName(pProcParams->szRefChrom))<1)
			continue;
		if(pProcParams->pIncludes[Idx]->InAnyFeature(BEDChromID,SubRefOfs,SubRefEndOfs))
			break;
		}
	if(Idx == pProcParams->NumIncludes)
		return(false);
	}

if(pProcParams->NumExcludes)
	{
	for(Idx = 0; Idx < pProcParams->NumExcludes; Idx++)
		{
		if(pProcParams->pExcludes[Idx] == NULL) // should'nt ever be NULL but...
			continue;
		if((BEDChromID = pProcParams->pExcludes[Idx]->LocateChromIDbyName(pProcParams->szRefChrom))<1)
			continue;
		if(pProcParams->pExcludes[Idx]->InAnyFeature(BEDChromID,SubRefOfs,SubRefEndOfs))
			return(false);
		}
	}
return(true);
}



// determine central step distributions of specified aligned subsequence
// 
bool
ProcessAlignment(int RefChromID,int SubRefOfs,int SubSeqLen,int SubSeqStartIdx,tsProcParams *pProcParams)
{
static int m_NumExcluded = 0;			// count of all putative subsequences of length NMer which were excluded
static int m_NumIncluded = 0;			// count of all putative subsequences of length NMer which were included
int MidBaseIdx = pProcParams->NMer/2;
while(SubSeqLen >= pProcParams->NMer)
	{
	// ensure the subsequence is in an included region and not part of an excluded region
	if(IncludeFilter(RefChromID,SubRefOfs,SubRefOfs+pProcParams->NMer-1,pProcParams))
		{
		m_NumIncluded++;
		if(ProcessSubSeq(RefChromID,SubRefOfs,SubSeqStartIdx,MidBaseIdx,pProcParams))
			{
			SubSeqStartIdx += MidBaseIdx;
			SubSeqLen-=MidBaseIdx;
			SubRefOfs+=MidBaseIdx;
			}
		}
	else
		m_NumExcluded++;
	SubSeqStartIdx += 1;
	SubSeqLen -= 1;
	SubRefOfs += 1;
	}
return(true);
}


bool
ProcessSubSeq(int RefChromID,int SubRefOfs,int SubSeqStartIdx,int MidBaseIdx,tsProcParams *pProcParams)
{
int Idx;
int VIdx;
int MerIdx;
etSeqBase *pVSeq;
etSeqBase VBase;
etSeqBase RefBase;		
etSeqBase CentBase[3];			// centroid bases - CentBase[0] = ref, CentBase[1] = rel, CentBase[2] = outgroup
int ConBases[4];
int FeatureBits;
int BitMsk;
int FeaturePsn;
int *pCnt;
ConBases[0] = ConBases[1] = ConBases[2] = ConBases[3] = 0;

MerIdx = 0;
for(Idx = 0; Idx < pProcParams->NMer ; Idx++)
	{
	for(VIdx = 0; VIdx < pProcParams->NumSpecies2Align; VIdx++)
		{
		pVSeq = pProcParams->pSeqs[VIdx];
		VBase = pVSeq[Idx+SubSeqStartIdx] & ~cRptMskFlg;
		
		if(VBase > eBaseT && (VIdx < 2 || !pProcParams->bInferLoose))		// not interested in any bases not a,c,g or t
				return(false);			
		
		if(!VIdx)					// if reference sequence
			RefBase = VBase;		// then record its base
		else						// else must be either rel or outgroup
			{
			if(Idx != MidBaseIdx)
				{
				if(VBase != RefBase && (VIdx < 2 || !pProcParams->bInferLoose))	// if flanking base then must match exactly 
					return(false);
				}
			}

		// if centroid base
		// at least one of ref or rel must be identical to outgroup before ancestral base can be infered
		if(Idx == MidBaseIdx)
			{
			CentBase[VIdx] = VBase;
			if(VIdx == (pProcParams->NumSpecies2Align-1))	// only when ref+rel+outgroup known then can determine if reasonable to infer outgroup as ancesteral base
				{
				if(VIdx == 2)								// if outgroup is being processed
					{
					if(!(CentBase[0] == CentBase[2] || CentBase[1] == CentBase[2])) // nothing can be infered if neither ref or rel match outgroup
						return(false);
					RefBase = CentBase[2];			// centroid will always be outgroup
					if(CentBase[0] == CentBase[1])		// if ref == rel then no mutation
						ConBases[CentBase[0]] = 1;
					else
						{
						if(CentBase[0] != CentBase[2])  // relative to the outgroup has ref?
							ConBases[CentBase[0]] = 1;	// ref has mutated
						else
							ConBases[CentBase[1]]= 1; // rel has mutated
						}
					}
				else									// no outgroup is being processed, only interested in mutational rate
					{
					if(CentBase[0] == CentBase[1])		// if ref == rel then no mutation
						ConBases[CentBase[0]] = 1;
					else
						ConBases[CentBase[1]]= 1;		// assume that rel has mutated
					}
				}
			}
		}
	MerIdx <<= 2;
	MerIdx += (int)RefBase;
	}

// have a sequence which we are interested in
// clasify it as to which region it is from
if(pProcParams->BEDChromID >= 0)
	{
	if(pProcParams->BEDChromID > 0)
		FeatureBits = pProcParams->pBiobed->GetFeatureBits(pProcParams->BEDChromID,SubRefOfs,SubRefOfs+pProcParams->NMer-1,cRegionFeatBits,pProcParams->UpDnStreamLen);
	else
		FeatureBits = 0;
	FeaturePsn = 0;		// default to intergenic if no feature bits set
	if(FeatureBits)		
		{
		BitMsk = cFeatBitCDS;
		for(Idx = 1; Idx < 7; Idx++,BitMsk <<= 1)
			{
			if(BitMsk & FeatureBits)
				{
				if(FeaturePsn)		// if already have feature
					return(true);	// although was sequence of interest, more than one feature bit so can't contribute to stats
				FeaturePsn = Idx * cRegionLen;
				if(pProcParams->bMultipleFeatBits)	// if multiple features allowed then don't check for any additional
					break;
				}
			}
		}
	pCnt = &pProcParams->pCntStepCnts[(MerIdx * cRegionLen * 7 ) + FeaturePsn];
	pCnt[0] +=ConBases[0];
	pCnt[1] +=ConBases[1];
	pCnt[2] +=ConBases[2];
	pCnt[3] +=ConBases[3];
		
	pProcParams->bCountsAvail = true;
	if(pProcParams->iWindowSize)
		ChkOutputResults(pProcParams->szRefChrom,SubRefOfs,pProcParams);	// output results as may be appropriate
	}
else
	{
	pCnt = &pProcParams->pCntStepCnts[(MerIdx * cRegionLen)];
	pCnt[0] +=ConBases[0];
	pCnt[1]+=ConBases[1];
	pCnt[2]+=ConBases[2];
	pCnt[3]+=ConBases[3];
	pProcParams->bCountsAvail = true;
	if(pProcParams->iWindowSize)
		ChkOutputResults(pProcParams->szRefChrom,SubRefOfs,pProcParams);	// output results as appropriate
	}
return(true);
}



bool
ProcessSequence(int RefChromID,int SubRefOfs,tsProcParams *pProcParams)
{
int Idx;
int Psn;
int MerIdx;
int SubSeqLen;
 
etSeqBase Base;		
int FeatureBits;
int BitMsk;
int FeaturePsn;
int nMer;
etSeqBase *pSeq;
if(pProcParams->iWindowSize)
	ChkOutputGenomeResults(pProcParams->szRefChrom,SubRefOfs,pProcParams);	// output results as appropriate

// iterate over all bases in sequence, slough any subsequences with non-nucleotides
// create mask for use when determining MerIdx
nMer = pProcParams->NMer;
SubSeqLen = 0;
MerIdx = 0;
pSeq = pProcParams->pSeq;
for(Psn = 0; Psn < pProcParams->SeqLen; Psn++,pSeq++,SubRefOfs++)
	{
	if((Base = (*pSeq & ~cRptMskFlg)) > eBaseT)
		{
		SubSeqLen = 0;
		MerIdx = 0;
		continue;
		}
	MerIdx <<= 2;
	MerIdx += (int)Base;
	MerIdx &= pProcParams->MerIdxMsk;
	if(++SubSeqLen < nMer)
		continue;

	// ensure the subsequence is in an included region and not part of an excluded region
	if(!IncludeFilter(RefChromID,SubRefOfs,SubRefOfs+nMer-1,pProcParams))
		{
		if(!pProcParams->bOverlapNMers)
			{
			SubSeqLen = 0;
			MerIdx = 0;
			}
		continue;
		}

	// have a sequence which we are interested in
	if(pProcParams->iWindowSize)
		ChkOutputGenomeResults(pProcParams->szRefChrom,SubRefOfs,pProcParams);	// output results as appropriate
	// classify it as to which region it is from
	if(pProcParams->BEDChromID >= 0)
		{
		if(pProcParams->BEDChromID > 0)
			FeatureBits = pProcParams->pBiobed->GetFeatureBits(pProcParams->BEDChromID,SubRefOfs,SubRefOfs+nMer-1,cRegionFeatBits,pProcParams->UpDnStreamLen);
		else
			FeatureBits = 0;
		FeaturePsn = 0;
		if(FeatureBits)
			{
			BitMsk = cFeatBitCDS;
			for(Idx = 1; Idx < 7; Idx++,BitMsk <<= 1)
				{
				if(BitMsk & FeatureBits)
					{
					if(FeaturePsn)
						break;	// although was sequence of interest, more than one feature bit so can't contribute to stats
					FeaturePsn = Idx;
					if(pProcParams->bMultipleFeatBits)
						break;
					}
				}

			if(Idx != 7 && !pProcParams->bMultipleFeatBits)
				{
				if(!pProcParams->bOverlapNMers)
					{
					SubSeqLen = 0;
					MerIdx = 0;
					}
				continue;
				}
			}
		pProcParams->pCntStepCnts[(MerIdx * 7) + FeaturePsn]++;
		}
	else
		pProcParams->pCntStepCnts[MerIdx]++;

	pProcParams->bCountsAvail = true;
	if(pProcParams->iWindowSize)
		ChkOutputGenomeResults(pProcParams->szRefChrom,SubRefOfs,pProcParams);	// output results as appropriate
	if(!pProcParams->bOverlapNMers)
		{
		SubSeqLen = 0;
		MerIdx = 0;
		}
	}
if(pProcParams->iWindowSize)
	ChkOutputGenomeResults(pProcParams->szRefChrom,SubRefOfs,pProcParams);	// output results as appropriate
return(true);
}



bool
ProcessSNPSequence(int RefChromID,int SubRefOfs,tsProcParams *pProcParams)
{
int Base;
int ChromOfs;
char Strand;
int Idx;
int MerIdx;
int SNPChromID;
int CurSNPID;
int RelOfs;		
int FeatureBits;
int BitMsk;
int FeaturePsn;
etSeqBase *pBase;
int *pCnt;

if(pProcParams->iWindowSize)
	ChkOutputResults(pProcParams->szRefChrom, SubRefOfs,pProcParams);

SNPChromID = pProcParams->pSNPs->LocateChrom(pProcParams->szRefChrom);
CurSNPID = 0;
// iterate over all SNPs on this chromosome
while((CurSNPID = pProcParams->pSNPs->GetNextSNP(CurSNPID,SNPChromID))>0)
	{
	Strand = pProcParams->pSNPs->GetStrand(CurSNPID);
	if(pProcParams->OnStrand != '*' && Strand != pProcParams->OnStrand)
		continue;
	ChromOfs = pProcParams->pSNPs->GetOffset(CurSNPID);
	if(ChromOfs > pProcParams->SeqLen - pProcParams->NMer)
		break;
	if(ChromOfs < pProcParams->FlankLen)
		continue;
	RelOfs = SubRefOfs + ChromOfs;
		// ensure the SNP is in an included region and not part of an excluded region
	if(!IncludeFilter(RefChromID,RelOfs,RelOfs,pProcParams))
		continue;
	
	MerIdx = 0;
	pBase = &pProcParams->pSeq[ChromOfs - pProcParams->FlankLen];
	for(Idx = 0; Idx < pProcParams->NMer; Idx++,pBase++)
		{
		Base = *pBase & ~cRptMskFlg;
		if(Base > eBaseT)
			break;
		MerIdx <<= 2;
		MerIdx += (int)Base;
		}
	if(Idx < pProcParams->NMer)		// can only handle subsequences which contain a,c,g and t
		continue;

	Strand = pProcParams->pSNPs->GetStrand(CurSNPID);
	for(Base =eBaseA; Base <= eBaseT; Base++)
		{
		if(pProcParams->pSNPs->IsPM(CurSNPID,(teSeqBases)Base))
			{
			// have a sequence which we are interested in
			// clasify it as to which region it is from
			if(pProcParams->BEDChromID >= 0)
				{
				if(pProcParams->BEDChromID > 0)
					FeatureBits = pProcParams->pBiobed->GetFeatureBits(pProcParams->BEDChromID,RelOfs,RelOfs,cRegionFeatBits,pProcParams->UpDnStreamLen);
				else
					FeatureBits = 0;
				FeaturePsn = 0;		// default to intergenic if no feature bits set
				if(FeatureBits)		
					{
					BitMsk = cFeatBitCDS;
					for(Idx = 1; Idx < 7; Idx++,BitMsk <<= 1)
						{
						if(BitMsk & FeatureBits)
							{
							if(FeaturePsn)		// if already have feature
								continue;		// although was sequence of interest, more than one feature bit so can't contribute to stats
							FeaturePsn = Idx * cRegionLen;
							if(pProcParams->bMultipleFeatBits)	// if multiple features allowed then don't check for any additional
								break;
							}
						}
					}
				pCnt = &pProcParams->pCntStepCnts[(MerIdx * cRegionLen * 7 ) + FeaturePsn];
				pCnt[Base]++;
				pProcParams->bCountsAvail = true;
				if(pProcParams->iWindowSize)
					ChkOutputResults(pProcParams->szRefChrom,RelOfs,pProcParams);	// output results as may be appropriate
				}
			else
				{
				pCnt = &pProcParams->pCntStepCnts[(MerIdx * cRegionLen)];
				pCnt[Base]++;
				pProcParams->bCountsAvail = true;
				if(pProcParams->iWindowSize)
					ChkOutputResults(pProcParams->szRefChrom,RelOfs,pProcParams);	// output results as appropriate
				}
			}
		}
	}
return(true);
}

// TrimQuotes
// Removes any leading and trailing quotes from pszTxt
int
TrimQuotes(char *pszTxt)
{
char *pDst = pszTxt;
char Chr;
int Len = 0;
while((Chr = *pszTxt++))	
	{
	if((!Len || *pszTxt == '\0') && (Chr == '"' || Chr == '\''))
		continue;
	*pDst++ = Chr;
	Len++;
	}
*pDst = Chr;
return(Len);
}


// ParseNumSpecies
// Initialises pProcParams with parsed species names in space or comma delimited list ptd at by pszSpeciesList
// Returns number of species
int
ParseNumSpecies(char *pszSpeciesList,tsProcParams *pProcParams)
{
// parse out species list
char Chr;
char *pSpecies;
int NumSpecies = 0;
bool InToken = false;
if(pszSpeciesList == NULL || *pszSpeciesList == '\0')
	return(0);

while(Chr = *pszSpeciesList++)
	{
	if(Chr == '"' || Chr == '\'') // change any single or double quotes into spaces
		Chr = ' ';
	if(isspace(Chr) || Chr==',')
		{
		if(!InToken)			// slough whitespace or ',' if not inside a token parse
			continue;
		InToken = false;
		pszSpeciesList[-1] = '\0';
		if(pProcParams != NULL)
			{
			strncpy(pProcParams->szSpecies[NumSpecies],pSpecies,cMaxDatasetSpeciesChrom);
			pProcParams->szSpecies[NumSpecies][cMaxDatasetSpeciesChrom-1] = '\0';
			}
		pszSpeciesList[-1] = Chr;
		NumSpecies++;
		if(NumSpecies >= cMAFmaxSpecies)
			break;
		continue;
		}
	if(!InToken)			// if not already inside token then start token 
		{
		pSpecies = pszSpeciesList-1;
		InToken = true;
		}
	}
if(InToken)
	{
	pszSpeciesList[-1] = '\0';
	if(pProcParams != NULL)
		{
		strncpy(pProcParams->szSpecies[NumSpecies],pSpecies,cMaxDatasetSpeciesChrom);
		pProcParams->szSpecies[NumSpecies][cMaxDatasetSpeciesChrom-1] = '\0';
		}
	pszSpeciesList[-1] = Chr;
	NumSpecies++;
	}
return(NumSpecies);
}


// ProcessAlignments
// Process ref+rel+outgroup in single pass
int 
ProcessAlignments(char *pszMAF,			 // source multialignment file
				  tsProcParams *pProcParams) // processing parameters
{
int RefSpeciesID;
int RefChromID;
int BEDChromID;
int PrevRefChromID;
int PrevDispRefChromID;
char *pszRefChrom;
int RefChromOfs;
char RefStrand;
int SpeciesIDs[cMAFmaxSpecies];
int Rslt;
etSeqBase *pSeq;
int RefAlignLen;
CMAlignFile *pAlignments;
int Idx;
int CurBlockID;

if((pAlignments = new CMAlignFile())==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create new instance of CMAlignFile\n");
	return(eBSFerrObj);
	}

if((Rslt=pAlignments->Open(pszMAF))!=eBSFSuccess)
	{
	while(pAlignments->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,pAlignments->GetErrMsg());

	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open MAF file %s",pszMAF);
	return(Rslt);
	}

// ensure all species are represented in multispecies alignment file plus get their species identifiers
for(Idx = 0; Idx < pProcParams->NumSpecies2Align; Idx++)
	{
	if((Rslt = SpeciesIDs[Idx] = pAlignments->LocateSpeciesID(pProcParams->szSpecies[Idx]))<1)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Species '%s' not represented in %s\nRepresented species are:",pProcParams->szSpecies[Idx],pszMAF);
		Rslt = pAlignments->GetNumSpecies();
		for(Idx=1;Idx<=Rslt;Idx++)
			gDiagnostics.DiagOut(eDLFatal,gszProcName," %s",pAlignments->GetSpeciesName(Idx));
		delete pAlignments;
		return(eBSFerrEntry);
		}
	if(!Idx)	// reference species is always the first species in the species name list
		RefSpeciesID = SpeciesIDs[0];
	}

CurBlockID = 0;
PrevRefChromID = 0;
PrevDispRefChromID = 0;
BEDChromID = 0;
pProcParams->RefSpeciesIdx = 0;		// reference sequence will always be 1st
pProcParams->NxtOutputOffset = pProcParams->iWindowSize;

// if only a single chromosome to be processed then ensure chromosme is represented in alignment file
if(pProcParams->szOnlyChrom[0]!='\0')
	{
	pProcParams->OnlyChromID=pAlignments->LocateChromID(RefSpeciesID,pProcParams->szOnlyChrom); //returns chromosome identifier
	if(pProcParams->OnlyChromID < 1)
		{

		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Chromosome '%s' specified as only chromosome to process is not present for species '%s' in alignment file '%s'\n",
				pProcParams->szOnlyChrom,pProcParams->szSpecies[0],pszMAF);
		delete pAlignments;
		return(eBSFerrChrom);
		}
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing chromosome '%s'...",pProcParams->szOnlyChrom);
	RefChromID = pProcParams->OnlyChromID;
	strcpy(pProcParams->szRefChrom,pProcParams->szOnlyChrom);
	if(pProcParams->pBiobed != NULL)
		{
		pProcParams->BEDChromID = pProcParams->pBiobed->LocateChromIDbyName(pProcParams->szRefChrom);
		if(pProcParams->BEDChromID < 1)
			{
			gDiagnostics.DiagOut(eDLWarn,gszProcName,"Unable to locate chrom '%s' in biobed file, no regional stats for this chromosome...\n",pProcParams->szRefChrom);
			pProcParams->BEDChromID = 0;
			}
		}
	else
		pProcParams->BEDChromID = -1;
		
	while((CurBlockID = pAlignments->NxtBlockChrom(CurBlockID,pProcParams->OnlyChromID))>0)
		{
				// not interested in alignments in which ref+rel+outgroup are not all represented
		if(pAlignments->GetNumSpecies(CurBlockID) < pProcParams->NumSpecies2Align)
			continue;
		if((RefAlignLen =  pAlignments->GetAlignLen(CurBlockID,RefSpeciesID)) < pProcParams->NMer)
			continue;
		if(pProcParams->bFiltLoConfidence && RefAlignLen < pProcParams->MinSubSeqLen)
			continue;

		// iterate over species aligned in current block, species are in priority order
		for(pProcParams->NumSpeciesAligned = Idx = 0; Idx < pProcParams->NumSpecies2Align; Idx++) 
			{
			// get each sequence
			if((pSeq = pAlignments->GetSeq(CurBlockID,SpeciesIDs[Idx]))==NULL) // some species are not present in all blocks
				{
				if(Idx < pProcParams->NumSpecies2Align)			// ref+rel plus outgroup must always be present
					break;
				continue;		
				}
			memcpy(pProcParams->pSeqs[pProcParams->NumSpeciesAligned++],pSeq,RefAlignLen);
			}
		// check that required number of species required were in alignment
		if(pProcParams->NumSpeciesAligned < pProcParams->NumSpecies2Align)	
			continue;

		// check that after InDel normalisation the seq length is still of interest
		if((RefAlignLen = NormaliseInDelColumns(pProcParams,RefAlignLen,pProcParams->NumSpeciesAligned))< pProcParams->NMer)
			continue;
		if(pProcParams->bFiltLoConfidence && RefAlignLen < pProcParams->MinSubSeqLen)
			continue;

		RefChromOfs= pAlignments->GetRelChromOfs(CurBlockID,RefSpeciesID);
		RefStrand = pAlignments->GetStrand(CurBlockID,RefSpeciesID);

		// this is where we can actually process the aligned block
		if(pProcParams->iWindowSize)
			ChkOutputResults(pProcParams->szRefChrom, RefChromOfs,pProcParams);
		ProcAlignBlock(RefChromID,RefChromOfs,RefAlignLen,pProcParams);
		}
	// requested chromosome processed
	if(pProcParams->iWindowSize)
		ChkOutputResults(pProcParams->szRefChrom, -1,pProcParams);
	else
		OutputResults(pProcParams->szRefChrom,0,pProcParams);
	}
else // iterate over reference blocks which are sorted by chrom then offset
	{
	while((CurBlockID = pAlignments->NxtBlock(CurBlockID))>0)
		{
		RefChromID = pAlignments->GetRelChromID(CurBlockID,RefSpeciesID);
		pszRefChrom = pAlignments->GetChromName(RefChromID);
		if(RefChromID != PrevDispRefChromID)
			{
			if(PrevDispRefChromID != 0)
				{
				if(pProcParams->iWindowSize)
					ChkOutputResults(pProcParams->szRefChrom, -1,pProcParams);
				else
					if(pProcParams->bPerChrom)
						OutputResults(pProcParams->szRefChrom,0,pProcParams);
				}
			PrevDispRefChromID = RefChromID;
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing chromosome %s",pszRefChrom);
			pProcParams->NxtOutputOffset = pProcParams->iWindowSize;
			strcpy(pProcParams->szRefChrom,pszRefChrom);
			if(pProcParams->pBiobed != NULL)
				{
				pProcParams->BEDChromID = pProcParams->pBiobed->LocateChromIDbyName(pProcParams->szRefChrom);
				if(pProcParams->BEDChromID < 1)
					{
					gDiagnostics.DiagOut(eDLInfo,gszProcName,"Unable to locate chrom '%s' in biobed file, no regional stats for this chromosome...\n",pProcParams->szRefChrom);
					pProcParams->BEDChromID = 0;
					}
				}
			else
				pProcParams->BEDChromID = -1;
			}

			// not interested in alignments in which ref+rel+outgroup are not all represented
		if(pAlignments->GetNumSpecies(CurBlockID) < pProcParams->NumSpecies2Align)
			continue;
		if((RefAlignLen =  pAlignments->GetAlignLen(CurBlockID,RefSpeciesID)) < pProcParams->NMer)
			continue;
		if(pProcParams->bFiltLoConfidence && RefAlignLen < pProcParams->MinSubSeqLen)
			continue;

			// iterate over species aligned in current block, species are in priority order
		for(pProcParams->NumSpeciesAligned = Idx = 0; Idx < pProcParams->NumSpecies2Align; Idx++) 
			{
			// get each sequence
			if((pSeq = pAlignments->GetSeq(CurBlockID,SpeciesIDs[Idx]))==NULL) // some species are not present in all blocks
				{
				if(Idx < pProcParams->NumSpeciesAligned)		// ref+rel+outgroup must always be present
					break;
				continue;		
				}
			memcpy(pProcParams->pSeqs[pProcParams->NumSpeciesAligned++],pSeq,RefAlignLen);
			}
		// check that required number of species required were in alignment
		if(pProcParams->NumSpeciesAligned < pProcParams->NumSpecies2Align)	
			continue;

				// check that after InDel normalisation the seq length is still of interest
		if((RefAlignLen = NormaliseInDelColumns(pProcParams,RefAlignLen,pProcParams->NumSpeciesAligned))< pProcParams->NMer)
			continue;
		if(pProcParams->bFiltLoConfidence && RefAlignLen < pProcParams->MinSubSeqLen)
			continue;


		if(RefChromID != PrevRefChromID)
			PrevRefChromID = RefChromID;

		RefChromOfs= pAlignments->GetRelChromOfs(CurBlockID,RefSpeciesID);
		RefStrand = pAlignments->GetStrand(CurBlockID,RefSpeciesID);

		// this is where we can actually process the aligned block
		if(pProcParams->iWindowSize)
			ChkOutputResults(pProcParams->szRefChrom, RefChromOfs,pProcParams);
		if(!ProcAlignBlock(RefChromID,RefChromOfs,RefAlignLen,pProcParams))
			break;
		}
	if(pProcParams->iWindowSize)
		ChkOutputResults(pProcParams->szRefChrom, -1,pProcParams);
	else
		{
		if(pProcParams->bPerChrom)
			OutputResults(pProcParams->szRefChrom,0,pProcParams);
		else
			OutputResults((char *)"Genome",0,pProcParams);
		}
	}
delete pAlignments;
return(eBSFSuccess);
}


// NormaliseInDelColumns
// Because multialignments may be the result of merged alignments resulting in InDels being generated back into the 
// existing reference sequence then if only a subset of species are being processed there could be all InDels in any column
// of the subset sequences.
// This function will delete all columns which only contain InDels
// Returns the subset sequence length after InDel columns have been deleted
int
NormaliseInDelColumns(tsProcParams *pProcParams,int AlignLen,int NumSpecies)
{
etSeqBase *pSeq;
etSeqBase *pSrcSeq;
etSeqBase SeqBase;
int NumVInDels;
int SeqIdx;
int VIdx;
int FirstInDel=-1;
int NormAlignLen=0;
for(SeqIdx = 0; SeqIdx < AlignLen; SeqIdx++)
	{
	NumVInDels = 0;					// assume no InDels in column
	for(VIdx = 0; VIdx < NumSpecies; VIdx++)
		{
		pSeq = pProcParams->pSeqs[VIdx];
		SeqBase = pSeq[SeqIdx] & ~cRptMskFlg;
		if(SeqBase == eBaseInDel)
			NumVInDels++;
		else
			break;				
		}

	if(NumVInDels == NumSpecies)	// all were InDels?
		{
		// mark col for deletion
		for(VIdx = 0; VIdx < NumSpecies; VIdx++)
			{
			pSeq = pProcParams->pSeqs[VIdx];
			pSeq[SeqIdx] = (etSeqBase)0x0ff;
			}
		if(FirstInDel == -1)		// note idx of first InDel
			FirstInDel = SeqIdx;
		}
	else
		NormAlignLen++;				// accept this column
	}
if(NormAlignLen == AlignLen)	// if no columns to delete
	return(AlignLen);

// have at least one column which is all InDels and is to be deleted
for(VIdx = 0; VIdx < NumSpecies; VIdx++)
	{
	pSeq = pProcParams->pSeqs[VIdx];
	pSeq = &pSeq[FirstInDel];
	pSrcSeq = pSeq;
	for(SeqIdx = FirstInDel; SeqIdx < AlignLen; SeqIdx++,pSrcSeq++)
		{
		if(*pSrcSeq != (etSeqBase)0x0ff)
			*pSeq++ = *pSrcSeq;
		}
	}
return(NormAlignLen);
}



// ChkOutputGenomeResults
// 
bool
ChkOutputGenomeResults(char *pszChrom, int ChromOffset, tsProcParams *pProcParams)
{
if(!pProcParams->iWindowSize)
	return(false);
while(ChromOffset < 0 || ChromOffset > pProcParams->NxtOutputOffset)
	{
	if(pProcParams->bGenEmptyWindows || pProcParams->bCountsAvail)
		OutputGenomeResults(pszChrom, pProcParams->NxtOutputOffset, pProcParams);
	pProcParams->bCountsAvail = false;
	pProcParams->NxtOutputOffset += pProcParams->iWindowSize;
	if(ChromOffset < 0)
		break;
	}
return(true);
}

// ChkOutputResults
bool
ChkOutputResults(char *pszChrom, int ChromOffset, tsProcParams *pProcParams)
{
if(!pProcParams->iWindowSize)
	return(false);
while(ChromOffset < 0 || ChromOffset > pProcParams->NxtOutputOffset)
	{
	if(pProcParams->bGenEmptyWindows || pProcParams->bCountsAvail)
		OutputResults(pszChrom, pProcParams->NxtOutputOffset, pProcParams);
	pProcParams->bCountsAvail = false;
	pProcParams->NxtOutputOffset += pProcParams->iWindowSize;
	if(ChromOffset < 0)
		break;
	}
return(true);
}


bool
OutputGenomeResults(char *pszChrom, int ChromOffset, tsProcParams *pProcParams)
{
static bool bOutputHdrFirst = true;
char szSubSeq[512];
char szLineBuff[2058];
int Idx;
int Steps;
int Instances;
int Len;
int *pStep;

if(bOutputHdrFirst)
	{
	bOutputHdrFirst = false;
	Len = sprintf(szLineBuff,"\"Area\",\"Offset\",\"SeqID\",\"centroid\",\"centroid3\",\"Seq\",\"TotInstances\"");
	if(pProcParams->Regions > 1)
		Len += sprintf(&szLineBuff[Len],",\"IGInstances\",\"US5Instances\",\"UTR5Instances\",\"CDSInstances\",\"IntronInstances\",\"UTR3Instances\",\"DS3Instances\"");
	CUtility::SafeWrite(pProcParams->hRsltsFile,szLineBuff,Len);
	}

pStep = pProcParams->pCntStepCnts;
for(Idx = 0; Idx < pProcParams->NumCnts; Idx++, pStep += pProcParams->Regions)
	{
	for(Instances = Steps = 0; Steps < pProcParams->Regions; Steps++)
		Instances += pStep[Steps];
	if(!pProcParams->bGenEmptyRows && Instances == 0)
		continue;
	Idx2Bases(Idx,pProcParams->NMer,szSubSeq);
	Len = sprintf(szLineBuff,"\n\"%s\",%d,%d,\"%c\",\"%c%c%c\",\"%s\",%d",
		pszChrom,ChromOffset,Idx,
		szSubSeq[pProcParams->FlankLen],			// gives centroids, makes it easier for manual sorting
		pProcParams->NMer >= 3 ? szSubSeq[pProcParams->FlankLen-1] : ' ',
		szSubSeq[pProcParams->FlankLen],
		pProcParams->NMer >= 2 ? szSubSeq[pProcParams->FlankLen+1] : ' ',
		szSubSeq,Instances);
	if(pProcParams->bGenEmptyRows && Instances == 0)
		{
		for(Steps = 0; Steps < pProcParams->Regions; Steps++)
			Len += sprintf(&szLineBuff[Len],",0");
		CUtility::SafeWrite(pProcParams->hRsltsFile,szLineBuff,Len);
		continue;
		}
	// for each region generate raw stats
	for(Steps = 0; Steps < pProcParams->Regions; Steps++)
		Len += sprintf(&szLineBuff[Len],",%d",pStep[MapStep2Region(Steps)]);
	CUtility::SafeWrite(pProcParams->hRsltsFile,szLineBuff,Len);
	}
memset(pProcParams->pCntStepCnts,0,pProcParams->NumCnts * pProcParams->Regions * sizeof(int));
pProcParams->bCountsAvail = false;
return(true);
}

bool	
OutputResults(char *pszChrom, int ChromOffset, tsProcParams *pProcParams)
{
static bool bOutputHdrFirst = true;
char szSubSeq[512];
char szLineBuff[2058];
int Idx;
int Steps;
int RegionIdx;
int Instances;
int Len;
int *pStep;
int *pCnt;

if(bOutputHdrFirst)
	{
	bOutputHdrFirst = false;
	if(pProcParams->Regions==1)
		Len = sprintf(szLineBuff,"\"Area\",\"Offset\",\"SeqID\",\"Seq\",\"Instances\",\"BaseA\",\"BaseC\",\"BaseG\",\"BaseT\"\n");
	else
		{
		Len = sprintf(szLineBuff,"\"Chrom\",\"Offset\",\"SeqID\",\"centroid\",\"centroid3\",\"Seq\",\"Instances\",");
		Len += sprintf(&szLineBuff[Len],"\"IGBaseA\",\"IGBaseC\",\"IGBaseG\",\"IGBaseT\",");
		Len += sprintf(&szLineBuff[Len],"\"USBaseA\",\"USBaseC\",\"USBaseG\",\"USBaseT\",");
		Len += sprintf(&szLineBuff[Len],"\"UTR5BaseA\",\"UTR5BaseC\",\"UTR5BaseG\",\"UTR5BaseT\",");
		Len += sprintf(&szLineBuff[Len],"\"CDSBaseA\",\"CDSBaseC\",\"CDSBaseG\",\"CDSBaseT\",");
		Len += sprintf(&szLineBuff[Len],"\"IntronBaseA\",\"IntronBaseC\",\"IntronBaseG\",\"IntronBaseT\",");
		Len += sprintf(&szLineBuff[Len],"\"UTR3BaseA\",\"UTR3BaseC\",\"UTR3BaseG\",\"UTR3BaseT\",");
		Len += sprintf(&szLineBuff[Len],"\"DS3BaseA\",\"DS3BaseC\",\"DS3BaseG\",\"DS3BaseT\"");
		}
	CUtility::SafeWrite(pProcParams->hRsltsFile,szLineBuff,Len);
	}

pStep = pProcParams->pCntStepCnts;
for(Idx = 0; Idx < pProcParams->NumCnts; Idx++, pStep += cRegionLen * pProcParams->Regions)
	{
	for(Instances = Steps = 0; Steps < (cRegionLen * pProcParams->Regions); Steps++)
		Instances += pStep[Steps];
	if(!pProcParams->bGenEmptyRows && Instances == 0)
		continue;
	Idx2Bases(Idx,pProcParams->NMer,szSubSeq);
	Len = sprintf(szLineBuff,"\n\"%s\",%d,%d,\"%c\",\"%c%c%c\",\"%s\",%d",
		pszChrom,ChromOffset,Idx,
		szSubSeq[pProcParams->FlankLen],			// gives centroids, makes it easier for manual sorting
		pProcParams->NMer >= 3 ? szSubSeq[pProcParams->FlankLen-1] : ' ',
		szSubSeq[pProcParams->FlankLen],
		pProcParams->NMer >= 3 ? szSubSeq[pProcParams->FlankLen+1] : ' ',
		szSubSeq,Instances);

	if(pProcParams->bGenEmptyRows && Instances == 0)
		{
		for(Steps = 0; Steps < pProcParams->Regions; Steps++)
			Len += sprintf(&szLineBuff[Len],",0,0,0,0");
		CUtility::SafeWrite(pProcParams->hRsltsFile,szLineBuff,Len);
		continue;
		}

	// for each region output the counts
	for(Steps = 0; Steps < pProcParams->Regions; Steps++)
		{
		RegionIdx = MapStep2Region(Steps);
		pCnt = &pStep[RegionIdx*cRegionLen];
		Len += sprintf(&szLineBuff[Len],",%d,%d,%d,%d",	pCnt[0],pCnt[1],pCnt[2],pCnt[3]);
		}
	CUtility::SafeWrite(pProcParams->hRsltsFile,szLineBuff,Len);
	}
memset(pProcParams->pCntStepCnts,0,pProcParams->NumCnts * cRegionLen * pProcParams->Regions * sizeof(int));
pProcParams->bCountsAvail = false;
return(true);
}


int
MapStep2Region(int Step)
{
int RegionIdx;
switch(Step) {
	case 0:
		RegionIdx = 0;	// IG
		break;
	case 1:
		RegionIdx = 5;	// 5'US
		break;
	case 2:
		RegionIdx = 2;	// 5'UTR
		break;
	case 3:
		RegionIdx = 1;	// CDS
		break;
	case 4:
		RegionIdx = 4;	// Intron
		break;
	case 5:
		RegionIdx = 3;	// 3'UTR
		break;
	case 6:				// 3'DS
		RegionIdx = 6;
		break;
	}
return(RegionIdx);
}


