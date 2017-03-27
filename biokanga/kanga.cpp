/*
 * CSIRO Open Source Software License Agreement (variation of the BSD / MIT License)
 * Copyright (c) 2017, Commonwealth Scientific and Industrial Research Organisation (CSIRO) ABN 41 687 119 230.
 * See LICENCE for the complete licence information (https://github.com/csiro-crop-informatics/biokanga/LICENSE)
 * Contact: Alex Whan <alex.whan@csiro.au>
 */

// kanga.cpp : Defines the entry point for the console application.

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

#include "../libbiokanga/SAMfile.h"
#include "biokanga.h"
#include "Aligner.h"


int
Process(etPMode PMode,					// processing mode
		UINT32 SampleNthRawRead,		// sample every Nth raw read, or read pair, for processing (1..10000)
		etFQMethod Quality,				// quality scoring for fastq sequence files
		bool bSOLiD,					// if true then processing in colorspace
		bool bBisulfite,				// if true then process for bisulfite methylation patterning
		etPEproc PEproc,				// paired reads alignment processing mode
		int PairMinLen,					// accept paired end alignments with observed insert size of at least this (default = 100)
		int PairMaxLen,					// accept paired end alignments with observed insert size of at most this (default = 1000)
		bool bPairStrand,				// accept paired ends if on same strand
		bool bPEcircularised,			// experimental - true if processing for PE spaning circularised fragments
		bool bPEInsertLenDist,			// true if stats file to include PE insert length distributions for each transcript
		eALStrand AlignStrand,			// align on to watson, crick or both strands of target
		int MinChimericLen,				// minimum chimeric length as a percentage (0 to disable, otherwise 50..99) of probe sequence length: negative if chimeric diagnostics to be reported
		bool bChimericRpt,				// report chimeric trimming detail for individual reads (default is not to report)
		int microInDelLen,				// microInDel length maximum
		int SpliceJunctLen,				// maximum splice junction length when aligning RNAseq reads
		int MinSNPreads,				// must be at least this number of reads covering any loci before processing for SNPs at this loci
		double QValue,					// QValue controlling FDR (Benjamini–Hochberg) SNP prediction
		double SNPNonRefPcnt,			// only process for SNP if more/equal than this percentage number of reads are non-ref at putative SNP loci (defaults to 25.0) 
		int MarkerLen,					// marker sequences of this length
		double MarkerPolyThres,			// maximum allowed marker sequence base polymorphism independent of centroid SNP (default 0.333, range 0.0 to 0.5)
		int PCRartefactWinLen,			// if >= 0 then window size to use when attempting to reduce the number of  PCR differential amplification artefacts (reads stacking to same loci)
		etMLMode MLMode,				// multimatch loci reads processing mode
		int MaxMLmatches,				// accept at most this number of multimatched alignments for any read
		bool bClampMaxMLmatches,		// accept as if MaxMLmatches even if more than this number of MaxMLmatches multimached alignments
		bool bLocateBestMatches,		// align for best matches, not just those 1H better than any other, upto at most MaxMLmatches
		int MaxNs,					    // allow at most this number of indeterminate eBaseNs in read before deeming as nonalignable
		int MinEditDist,				// any matches must have at least this edit distance to the next best match
		int MaxSubs,					// maximum number of substitutions allowed per 100bp of actual read length
		int Trim5,						// trim this number of bases from 5' end of reads when loading the reads
		int Trim3,						// trim this number of bases from 3' end of reads when loading the reads
		int MinAcceptReadLen,			// only accepting reads for alignment if at least this length after any end trimming
		int MaxAcceptReadLen,			// only accepting reads for alignment if no longer than this length after any end trimming
		int MinFlankExacts,				// trim matched reads on 5' and 3' flanks until at least this number of exactly matching bases in flanks
		int PCRPrimerCorrect,			// initially align with MaxSubs+PCRPrimerCorrect subs allowed but then attempt to correct 5' PCR primer artefacts until read within MaxSubs
		int MaxRptSAMChroms,			// report all SAM chroms or sequences if number of reference chroms <= this limit (defaults to 10000)
		etFMode FMode,					// output format mode
		teSAMFormat SAMFormat,			// if SAM output format then could be SAM, BAM or BAM compressed dependent on the file extension used
		int SitePrefsOfs,				// offset read start sites when processing  octamer preferencing, range -100..100
		int NumThreads,					// number of worker threads to use
		char *pszTrackTitle,			// track title if output format is UCSC BED
		int NumPE1InputFiles,			// number of input PE1 or single ended file specs
		char *pszPE1InputFiles[],		// names of input files (wildcards allowed unless processing paired ends) containing raw reads
		int NumPE2InputFiles,			// number of input PE2 file specs
		char *pszPE2InputFiles[],		// optional raw paired reads are in these files
		char *pszPriorityRegionFile,	// optional high priority BED file contains prioritised region loci
		bool bFiltPriorityRegions,		// true if non-priority alignments to be filtered out 
		char *pszOutFile,				// where to write alignments
		char *pszSNPFile,				// Output SNPs (CSV format) to this file (default is to output file name with '.snp' appended)
		char *pszMarkerFile,			// Output markers to this file
		char *pszSNPCentroidFile,		// Output SNP centorids (CSV format) to this file (default is for no centroid processing)
		char *pszSfxFile,				// target as suffix array
		char *pszStatsFile,				// aligner induced substitutions stats file
		char *pszMultiAlignFile,		// file to contain reads which are aligned to multiple locations
		char *pszNoneAlignFile,			// file to contain reads which were non-alignable
		char *pszSitePrefsFile,			// file to contain aligned reads octamer preferencing
		char *pszLociConstraintsFile,	// loci base constraints file
		char *pszContamFile,			// file containing contaminant sequences
		int	NumIncludeChroms,			// number of chromosome regular expressions to include
		char **ppszIncludeChroms,		// array of include chromosome regular expressions
		int	NumExcludeChroms,			// number of chromosome expressions to exclude
		char **ppszExcludeChroms);		// array of exclude chromosome regular expressions

#ifdef _WIN32
int kanga(int argc, char* argv[])
{
// determine my process name
_splitpath(argv[0],NULL,NULL,gszProcName,NULL);
#else
int
kanga(int argc, char** argv)
{
// determine my process name
CUtility::splitpath((char *)argv[0],NULL,gszProcName);
#endif

int iFileLogLevel;			// level of file diagnostics
int iScreenLogLevel;		// level of file diagnostics
char szLogFile[_MAX_PATH];	// write diagnostics to this file
int Rslt = 0;   			// function result code >= 0 represents success, < 0 on failure
int Idx;
int LenChromList;

int PMode;					// processing mode
UINT32 SampleNthRawRead;	// sample every Nth raw read (or read pair) for processing (1..10000)
bool bSOLiD;				// if true then process for colorspace (SOLiD)
bool bBisulfite;			// if true then process for bisulfite methylation patterning
int PCRartefactWinLen;		// if >= 0 then window size to use when attempting to reduce the number of  PCR differential amplification artefacts (reads stacking to same loci)
int FMode;					// format output mode
teSAMFormat SAMFormat;		// if SAM output format then could be either SAM or BAM compressed dependent on the file extension used

int NumberOfProcessors;		// number of installed CPUs
int NumThreads;				// number of threads (0 defaults to number of CPUs)
int Quality;				// quality scoring for fastq sequence files
int MinEditDist;			// any matches must have at least this edit distance to the next best match
int MaxSubs;				// maximum number of substitutions allowed per 100bp of read length
int MLMode;					// processing mode for multiple loci aligned reads
int MaxMLmatches;			// accept at most this number of multimached alignments for any read
bool bClampMaxMLmatches;	// accept as if MaxMLmatches even if more than this number of MaxMLmatches multimached alignments
int MaxNs;				    // allow at most this number of indeterminate eBaseNs in read before deeming as nonalignable
int MinFlankExacts;			// trim matched reads on 5' and 3' flanks until at least this number of exactly matching bases in flanks
int PCRPrimerCorrect;		// initially align with MaxSubs+PCRPrimerCorrect subs allowed but then attempt to correct 5' PCR primer artefacts until read within MaxSubs
int Trim5;					// trim this number of bases from 5' end of reads when loading the reads
int Trim3;					// trim this number of bases from 3' end of reads when loading the reads
int MinAcceptReadLen;				// only accepting reads for alignment if at least this length after any end trimming
int MaxAcceptReadLen;				// only accepting reads for alignment if no longer than this length after any end trimming
int MaxRptSAMSeqsThres;		// report all SAM chroms or sequences to SAM header if number of reference chroms <= this limit (defaults to 10000)
int AlignStrand;			// align on to watson, crick or both strands of target
int MinChimericLen;			// minimum chimeric length as a percentage (0 to disable, otherwise 50..99) of probe sequence length: negative if chimeric diagnostics to be reported
bool bChimericRpt;			// report chimeric trimming detail for individual reads (default is not to report)
int microInDelLen;			// microInDel length maximum
int SpliceJunctLen;			// maximum splice junction length when aligning RNAseq reads
int MinSNPreads;			// must be at least this number of reads covering any loci before processing for SNPs at this loci
double QValue;			    // used in Benjamini–Hochberg for SNP prediction
double SNPNonRefPcnt;		// only process for SNP if more/equal than this percentage number of reads are non-ref at putative SNP loci (defaults to 25.0) 

bool bLocateBestMatches;			// instead of unique matches then return multiple best matches

char szTrackTitle[cMaxDatasetSpeciesChrom];		// track title if output format is UCSC BED
char szRsltsFile[_MAX_PATH];			// results to this file
char szTargFile[_MAX_PATH];				// align against this target suffix array genome file

int NumPE1InputFiles;					// number of input PE1 or single ended file spe
char *pszPE1InputFiles[cMaxInFileSpecs];		// names of input files (wildcards allowed unless processing paired ends) containing raw reads
int NumPE2InputFiles;					// number of input PE2 file specs
char *pszPE2InputFiles[cMaxInFileSpecs];		// optional raw paired reads are in these files

int PEproc;								// paired end processing mode
int PairMinLen;							// accept paired end alignments with observed insert size of at least this (default = 100)
int PairMaxLen;							// accept paired end alignments with observed insert size of at most this (default = 1000)
bool bPEInsertLenDist;					// true if stats file to include PE insert length distributions for each transcript
bool bPairStrand;						// 5' and 3' are on same strand
bool bPEcircularised;					// experimental - true if processing for PE spaning circularised fragments

char szPriorityRegionFile[_MAX_PATH];	// optional exact match priority file contains prioritised region loci as a BED file
bool bFiltPriorityRegions;				// if priority regions requested then can also request that only alignments into these regions be reported

int MarkerLen;							// length of marker sequences to generate
double MarkerPolyThres;					// maximum allowed marker sequence base polymorphism independent of centroid SNP (default 0.333, range 0.0 to 0.5)
char szSNPFile[_MAX_PATH];				// Output SNPs (CSV format) to this file (default is to output file name with '.snp' appended)
char szMarkerFile[_MAX_PATH];			// Output markers (multifasta) to this file (name is SNP file with '.marker' appended)
char szSNPCentroidFile[_MAX_PATH];		// Output SNP centorids (CSV format) to this file (default is for no centroid processing)

char szLociConstraintsFile[_MAX_PATH];  // optional loci base constraints from this file
char szContamFile[_MAX_PATH];			// optional file containing contaminant sequences

char szStatsFile[_MAX_PATH];			// optional output basic distribution counts/stats to this file
char szMultiAlignFile[_MAX_PATH];		// optional output file to contain reads which are aligned to multiple locations
char szNoneAlignFile[_MAX_PATH];		// optional output file to contain reads which could not be aligned
char szSitePrefsFile[_MAX_PATH];		// optional output file to contain aligned reads start site octamer preferencing
int SitePrefsOfs;						// offset read start sites when processing  octamer preferencing, range -100..100


int NumIncludeChroms;
char *pszIncludeChroms[cMaxIncludeChroms];
int NumExcludeChroms;
char *pszExcludeChroms[cMaxExcludeChroms];

char szSQLiteDatabase[_MAX_PATH];	// results summaries to this SQLite file
char szExperimentName[cMaxDatasetSpeciesChrom+1];			// experiment name
char szExperimentDescr[1000];		// describes experiment

//
struct arg_lit  *help    = arg_lit0("h","help",                 "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to screen and logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"diagnostics log file");

struct arg_int *pmode = arg_int0("m","mode","<int>",		    "alignment processing mode: 0 - standard sensitivity, 1 - more sensitive (slower), 2 - ultra sensitive (slowest), 3 - less sensitive (quicker)");
struct arg_int *format = arg_int0("M","format","<int>",		    "output format: 0 - CSV loci only, 1 - CSV loci + match sequence, 2 - CSV loci + read sequence, 3 - CSV loci + read + match sequence, 4 - UCSC BED, 5 - SAM (BAM) format with accepted aligned, 6 - SAM (BAM) format with all reads (default: 5)");
struct arg_int *peproc = arg_int0("U","pemode","<int>",		    "paired end processing mode: 0 - none, 1 - paired end with recover orphan ends, 2 - paired end no orphan recovery, 3 - paired end with recover orphan ends treating orphan ends as SE, 4 - paired end no orphan recovery treating orphan ends as SE (default: 0)");

struct arg_lit  *bisulfite = arg_lit0("b","bisulfite",          "process for bisulfite methylation patterning");
struct arg_lit  *solid = arg_lit0("C","colorspace",             "process for colorspace (SOLiD)");
struct arg_int  *pcrartefactwinlen = arg_int0("k","pcrwin", "<int>",   "PCR differential amplification artefact reduction window length (default is for no artefact reduction, 0..250)");
struct arg_file *pe1inputfiles = arg_filen("i","in","<file>",0,cMaxInFileSpecs,"input from these raw sequencer read files (.gz allowed), wildcards allowed if single ended");
struct arg_file *pe2inputfiles = arg_filen("u","pair","<file>",0,cMaxInFileSpecs,"if raw paired end processing then input read pairs from these raw paired end files (.gz allowed)");

struct arg_file *priorityregionfile = arg_file0("B","priorityregionfile","<file>",	"prioritise exact match alignments to loci regions in this BED file");
struct arg_lit  *nofiltpriority = arg_lit0("V","nofiltpriority",  "do not filter priority region alignments");

struct arg_int  *pairminlen = arg_int0("d","pairminlen","<int>", "accept paired end alignments with observed insert sizes of at least this (default = 100)");
struct arg_int  *pairmaxlen = arg_int0("D","pairmaxlen","<int>", "accept paired end alignments with observed insert sizes of at most this (default = 1000)");
struct arg_lit  *pairstrand = arg_lit0("E","pairstrand",         "5' and 3' are on same strand");

struct arg_int  *alignstrand = arg_int0("Q","alignstrand","<int>", "align to this strand: 0 either, 1 Sense '+', 2 Antisense '-' (default is to align to either strand)");

struct arg_int  *samplenthrawread = arg_int0("#","samplenthrawread","<int>", "sample every Nth raw read or read pair for processing (default 1, range 1..10000)");

struct arg_int  *minchimericlen = arg_int0("c","minchimeric","<int>", "minimum chimeric length as a percentage of probe length (default is 0 to disable, otherwise 50..99)");
struct arg_lit  *chimericrpt = arg_lit0("0","chimericrpt",     "report chimeric trimming detail for individual reads (default is not to report)");


struct arg_int *qual = arg_int0("g","quality","<int>",		    "fastq quality scoring - 0 - Sanger or Illumina 1.8+, 1 = Illumina 1.3+, 2 = Solexa < 1.3, 3 = Ignore quality (default = 3)");
struct arg_file *sfxfile = arg_file1("I","sfx","<file>",		"align against this suffix array (kangax generated) file");
struct arg_file *outfile = arg_file1("o","out","<file>",		"output alignments to this file");

struct arg_int  *microindellen = arg_int0("a","microindellen","<int>", "accept microInDels inclusive of this length: 0 to 20 (default = 0 or no microIndels)");
struct arg_int  *splicejunctlen = arg_int0("A","splicejunctlen","<int>", "aligning RNA-seq, force flank trim, accept splice junctions separated by at most this distance: 25 to 100000 (default = 0 for DNA non-spliced aligning)");

struct arg_file *statsfile = arg_file0("O","stats","<file>",	"output aligner induced substitution distribution stats (not supported for '-M6' output mode) or paired end length distributions to this file");
struct arg_file *nonealignfile = arg_file0("j","nonealign","<file>",	"output unalignable reads to this file ");
struct arg_file *multialignfile = arg_file0("J","multialign","<file>",	"output multialigned reads to this file");

struct arg_file *siteprefsfile = arg_file0("8","siteprefs","<file>", "output aligned reads start site octamer preferencing to this file");
struct arg_int *siteprefsofs = arg_int0("9","siteprefsofs","<int>", "offset read start sites when processing site octamer preferencing, range -100..100 (default is -4)");

struct arg_file *lociconstraintsfile = arg_file0("5","lociconstraints","<file>", "Optional loci base constraints CSV file");
struct arg_file *contamsfile = arg_file0("H","contaminants","<file>", "Optional contaminant sequences multifasta file");


struct arg_file *snpfile = arg_file0("S","snpfile","<file>",	"Output SNPs (CSV or VCF if file name extension is '.vcf') to this file (default is to output as CSV to file name with '.snp' appended)");
struct arg_file *centroidfile = arg_file0("7","snpcentroid","<file>", "Output SNP centroid distributions (CSV format) to this file (default is for no centroid processing)");
struct arg_int *markerlen = arg_int0("K","markerlen","<int>", "output marker sequences of this length with centralised SNP, output to SNP file name with '.marker' apended (default no markers, range 25..500)");
struct arg_dbl *markerpolythres = arg_dbl0("G","markerpolythres","<dbl>", "maximum allowed marker sequence base polymorphism independent of centroid SNP (default 0.333, range 0.0 to 0.5)");

struct arg_int *mineditdist = arg_int0("e","editdelta","<int>",	"accepted matches must be at least this Hamming edit distance from the next best match (default is 1, max 2)");
struct arg_int *maxsubs = arg_int0("s","substitutions","<int>",	"accept up to this number of aligner induced substitutions per 100bp of individual read length (default is 10, range 0..15)");

struct arg_int *minflankexacts = arg_int0("x","minflankexacts","<int>",	"trim matching reads on 5' and 3' flanks until at least this number of exactly matching bases in flanks (default is 0 or no trimming)");
struct arg_int *pcrprimercorrect = arg_int0("6","pcrprimersubs","<int>",  "experimental - initially align with MaxSubs+pcrprimersubs subs allowed but then attempt to correct 5' PCR primer artefacts over first 12bp until read within MaxSubs");

struct arg_int *maxns = arg_int0("n","maxns","<int>",	        "maximum number, percentage if read length > 100, of indeterminate 'N's in reads before treating read as unalignable (default is 1, max 5)");

struct arg_str  *title = arg_str0("t","title","<string>",       "track title");

struct arg_str  *ExcludeChroms = arg_strn("Z","chromexclude",	"<string>",0,cMaxExcludeChroms,"high priority - regular expressions defining chromosomes to exclude");
struct arg_str  *IncludeChroms = arg_strn("z","chromeinclude",	"<string>",0,cMaxIncludeChroms,"low priority - regular expressions defining chromosomes to include");
struct arg_int *threads = arg_int0("T","threads","<int>",		"number of processing threads 0..128 (defaults to 0 which sets threads to number of CPU cores)");

struct arg_int *maxmlmatches = arg_int0("R","maxmulti","<int>",	"allow any read to match at most this many genome loci then process according to mlmode (default is 5)");
struct arg_lit *clampmaxmulti = arg_lit0("X","clampmaxmulti",	 "treat reads mapping to more than limit set with '-R<n>' as if exactly <n> matches (default is not to further process reads exceeding limit set with '-R<n>')");

struct arg_lit *bestmatches = arg_lit0("N","bestmatches",	    "accept best '-R<N>' multiple alignments with alignments having at most '-s<max subs>'");


struct arg_int *mlmode = arg_int0("r","mlmode","<int>",			"processing mode for reads mapped to multiple loci: 0 slough, 1 stats only, 2 rand, 3 cluster with uniques only, 4 cluster with uniques + other multi, 5 report all match loci up to '-R<limit>' (default is 0)");
struct arg_int *minsnpreads = arg_int0("p","snpreadsmin","<int>","minimum read coverage at loci before processing for SNP determination (default is 0 or no SNP processing)");
struct arg_dbl *snpnonrefpcnt = arg_dbl0("1","snpnonrefpcnt","<dbl>", "Min percentage non-ref bases at putative SNP loci (defaults to 25.0, range 0.1 to 35.0)");
struct arg_lit  *pecircularised    = arg_lit0("2","pecircularised",  "experimental - set true if processing for PE circularised fragments spanning");

struct arg_lit  *peinsertlendist    = arg_lit0("3","petranslendist",  "experimental - include PE length distributions for each transcript in output stats file");
struct arg_int *rptsamseqsthres = arg_int0("4","rptsamseqsthres","<int>", "if SAM format and refence sequences more than this threshold then only write sequences with reads aligned to SAM header (defaults to 10000)");

struct arg_dbl *qvalue = arg_dbl0("P","qvalue","<dbl>",		"QValue controlling FDR (Benjamini-Hochberg) SNP prediction (default is 0.05, range is <= 0.40)");
struct arg_int *trim5 = arg_int0("y","trim5","<int>",		"trim this number of bases from 5' end of reads when loading raw reads (default is 0)");
struct arg_int *trim3 = arg_int0("Y","trim3","<int>",		"trim this number of bases from 3' end of reads when loading raw reads (default is 0)");

struct arg_int *minacceptreadlen = arg_int0("l","minacceptreadlen","<int>",		"after any end trimming only accept read for further processing if read is at least this length (default is 50bp, range 15..2000)");
struct arg_int *maxacceptreadlen = arg_int0("L","maxacceptreadlen","<int>",		"after any end trimming only accept read for further processing if read is at no longer than this length (default is 500bp, range minacceptreadlen..2000)");

struct arg_file *summrslts = arg_file0("q","sumrslts","<file>",		"Output results summary to this SQLite3 database file");
struct arg_str *experimentname = arg_str0("w","experimentname","<str>",		"experiment name SQLite3 database file");
struct arg_str *experimentdescr = arg_str0("W","experimentdescr","<str>",	"experiment description SQLite3 database file");

struct arg_end *end = arg_end(200);

void *argtable[] = {help,version,FileLogLevel,LogFile,
					summrslts,experimentname,experimentdescr,
					pmode,samplenthrawread,alignstrand,minchimericlen,chimericrpt,pecircularised,peinsertlendist,microindellen,splicejunctlen,solid,pcrartefactwinlen,qual,mlmode,trim5,trim3,minacceptreadlen,maxacceptreadlen,maxmlmatches,rptsamseqsthres,clampmaxmulti,bisulfite,
					mineditdist,maxsubs,maxns,minflankexacts,pcrprimercorrect,minsnpreads,markerlen,markerpolythres,qvalue,snpnonrefpcnt,format,title,priorityregionfile,nofiltpriority,bestmatches,
					pe1inputfiles,peproc,pairminlen,pairmaxlen,pairstrand,pe2inputfiles,sfxfile,snpfile,centroidfile,
					outfile,nonealignfile,multialignfile,statsfile,siteprefsfile,siteprefsofs,lociconstraintsfile,contamsfile,ExcludeChroms,IncludeChroms,threads,
					end};

char **pAllArgs;
int argerrors;
argerrors = CUtility::arg_parsefromfile(argc,(char **)argv,&pAllArgs);
if(argerrors >= 0)
	argerrors = arg_parse(argerrors,pAllArgs,argtable);

/* special case: '--help' takes precedence over error reporting */
if (help->count > 0)
        {
		printf("\n%s %s %s, Version %s\nOptions ---\n", gszProcName,gpszSubProcess->pszName,gpszSubProcess->pszFullDescr,cpszProgVer);
        arg_print_syntax(stdout,argtable,"\n");
        arg_print_glossary(stdout,argtable,"  %-25s %s\n");
		printf("\nNote: Parameters can be entered into a parameter file, one parameter per line.");
		printf("\n      To invoke this parameter file then precede it's name with '@'");
		printf("\n      e.g. %s %s @myparams.txt\n",gszProcName,gpszSubProcess->pszName);
		printf("\nPlease report any issues regarding usage of %s at https://github.com/csiro-crop-informatics/biokanga/issues\n\n",gszProcName);
		return(1);
        }

    /* special case: '--version' takes precedence error reporting */
if (version->count > 0)
        {
		printf("\n%s %s Version %s\n",gszProcName,gpszSubProcess->pszName,cpszProgVer);
		return(1);
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

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Subprocess %s Version %s starting",gpszSubProcess->pszName,cpszProgVer);
	gExperimentID = 0;
	gProcessID = 0;
	gProcessingID = 0;
	szSQLiteDatabase[0] = '\0';
	szExperimentName[0] = '\0';
	szExperimentDescr[0] = '\0';


	if(experimentname->count)
		{
		strncpy(szExperimentName,experimentname->sval[0],sizeof(szExperimentName));
		szExperimentName[sizeof(szExperimentName)-1] = '\0';
		CUtility::TrimQuotedWhitespcExtd(szExperimentName);
		CUtility::ReduceWhitespace(szExperimentName);
		}
	else
		szExperimentName[0] = '\0';

	gExperimentID = 0;
	gProcessID = 0;
	gProcessingID = 0;
	szSQLiteDatabase[0] = '\0';
	szExperimentDescr[0] = '\0';

	if(summrslts->count)
		{
		strncpy(szSQLiteDatabase,summrslts->filename[0],sizeof(szSQLiteDatabase)-1);
		szSQLiteDatabase[sizeof(szSQLiteDatabase)-1] = '\0';
		CUtility::TrimQuotedWhitespcExtd(szSQLiteDatabase);
		if(strlen(szSQLiteDatabase) < 1)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: After removal of whitespace, no SQLite database specified with '-q<filespec>' option");
			return(1);
			}

		if(strlen(szExperimentName) < 1)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: After removal of whitespace, no SQLite experiment name specified with '-w<str>' option");
			return(1);
			}
		if(experimentdescr->count)
			{
			strncpy(szExperimentDescr,experimentdescr->sval[0],sizeof(szExperimentDescr)-1);
			szExperimentDescr[sizeof(szExperimentDescr)-1] = '\0';
			CUtility::TrimQuotedWhitespcExtd(szExperimentDescr);
			}
		if(strlen(szExperimentDescr) < 1)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: After removal of whitespace, no SQLite experiment description specified with '-W<str>' option");
			return(1);
			}

		gExperimentID = gSQLiteSummaries.StartExperiment(szSQLiteDatabase,false,true,szExperimentName,szExperimentName,szExperimentDescr);
		if(gExperimentID < 1)
			return(1);
		gProcessID = gSQLiteSummaries.AddProcess((char *)gpszSubProcess->pszName,(char *)gpszSubProcess->pszName,(char *)gpszSubProcess->pszFullDescr);
		if(gProcessID < 1)
			return(1);
		gProcessingID = gSQLiteSummaries.StartProcessing(gExperimentID,gProcessID,(char *)cpszProgVer);
		if(gProcessingID < 1)
			return(1);
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Initialised SQLite database '%s' for results summary collection",szSQLiteDatabase);
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"SQLite database experiment identifier for '%s' is %d",szExperimentName,gExperimentID);
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"SQLite database process identifier for '%s' is %d",(char *)gpszSubProcess->pszName,gProcessID);
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"SQLite database processing instance identifier is %d",gProcessingID);
		}
	else
		{
		szSQLiteDatabase[0] = '\0';
		szExperimentDescr[0] = '\0';
		}

	// ensure all filenames are initialised in case not user specified
	szRsltsFile[0] = '\0';
	szTargFile[0] = '\0';
	NumPE1InputFiles = 0;
	NumPE2InputFiles = 0;
	szSNPFile[0] = '\0';
	szSNPCentroidFile[0] = '\0';
	szMarkerFile[0] = '\0';

	szLociConstraintsFile[0] = '0';
	szContamFile[0] = '0';
	szStatsFile[0] = '\0';
	szMultiAlignFile[0] = '\0';
	szNoneAlignFile[0] = '\0';
	szSitePrefsFile[0] = '\0';
	MinChimericLen = 0;
	bChimericRpt = false;
	bLocateBestMatches = false;

	MinAcceptReadLen = cDfltMinAcceptReadLen;
	MaxAcceptReadLen = cDfltMaxAcceptReadLen;
	
	PMode = (etPMode)(pmode->count ? pmode->ival[0] : ePMdefault);
	if(PMode < ePMdefault || PMode >= ePMplaceholder)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Processing mode '-m%d' specified outside of range %d..%d\n",PMode,ePMdefault,(int)ePMplaceholder-1);
		exit(1);
		}

	SampleNthRawRead = samplenthrawread->count ? samplenthrawread->ival[0] : 1;
	if(SampleNthRawRead < 1)
		SampleNthRawRead = 1;
	else
		if(SampleNthRawRead > 10000)
			SampleNthRawRead = 10000;
	
	AlignStrand = (eALStrand)(alignstrand->count ? alignstrand->ival[0] : eALSboth);
	if(AlignStrand < eALSboth || AlignStrand >= eALSnone)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Aligned to strand '-Q%d' specified outside of range %d..%d\n",AlignStrand,eALSboth,(int)eALSnone-1);
		exit(1);
		}


	Quality = (etFQMethod)(qual->count ? qual->ival[0] : eFQIgnore);
	if(Quality < eFQSanger || Quality >= eFQplaceholder)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: fastq quality '-g%d' specified outside of range %d..%d\n",Quality,eFQSanger,(int)eFQplaceholder-1);
		exit(1);
		}


	MLMode = (etMLMode)(mlmode->count ? mlmode->ival[0] : 0);
	if(MLMode < eMLdefault || MLMode >= eMLplaceholder)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: multiple aligned reads processing mode '-r%d' specified outside of range 0..%d\n",MLMode,eMLplaceholder-1);
		exit(1);
		}

	FMode = (etFMode)(format->count ? format->ival[0] : eFMsam);
	if(FMode < eFMdefault || FMode >= eFMplaceholder)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Output format mode '-m%d' specified outside of range %d..%d\n",FMode,eFMdefault,(int)eFMplaceholder-1);
		exit(1);
		}


	MinFlankExacts = minflankexacts->count ? minflankexacts->ival[0] : 0;
	SpliceJunctLen = splicejunctlen->count ? splicejunctlen->ival[0] : 0;
	microInDelLen = microindellen->count ? microindellen->ival[0] : 0;
	bBisulfite = bisulfite->count ? true : false;
	bSOLiD = solid->count ? true : false;
	if(bBisulfite && bSOLiD)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: colorspace (SOLiD) '-C' and bisulfite methylation '-b' processing modes are mutually exclusive\n");
		exit(1);
		}

	PEproc = (etPEproc)(peproc->count ? peproc->ival[0] : 0);
	if(PEproc < ePEdefault || PEproc >= ePEplaceholder)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: PE processing mode '-U%d' specified outside of range %d..%d\n",FMode,ePEdefault,(int)ePEplaceholder-1);
		exit(1);
		}

	if(pe2inputfiles->count > 0 && PEproc == ePEdefault)	// if PE2 files specified then must also have PEproc mode specified
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: PE2 file(s) with '-u<files>' specified but PE processing mode '-U<mode>' not specified\n");
		exit(1);
		}

	if(PEproc > ePEdefault)			// only a limited subset of processing options/parameters avail if paired end processing
		{
		if(!(FMode == eFMdefault || FMode == eFMbed || FMode >= eFMsam))
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Sorry, currently output format '-M%d' not supported with paired end '-U%d' processing\n",MLMode,PEproc);
			exit(1);
			}

		if(bBisulfite)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Sorry, currently bisulfite processing '-b' not supported in paired end '-U%d' processing\n",PEproc);
			exit(1);
			}

		if(MLMode != eMLdefault)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Sorry, currently multiloci processing '-r%d' not supported in paired end '-U%d' processing\n",MLMode,PEproc);
			exit(1);
			}

		if(microInDelLen != 0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Sorry, currently microInDel processing '-a%d' not supported in paired end '-U%d' processing\n",microInDelLen,PEproc);
			exit(1);
			}

		if(SpliceJunctLen != 0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Sorry, currently RNA-seq splice junction processing '-A%d' not supported in paired end '-U%d' processing\n",SpliceJunctLen,PEproc);
			exit(1);
			}

		bPEcircularised = pecircularised->count ? true : false;
		bPEInsertLenDist = peinsertlendist->count ? true : false;
		}
	else
		{
		bPEcircularised = false;
		bPEInsertLenDist = false;
		}

	if(statsfile->count)
		{
		if(PEproc == ePEdefault && FMode == eFMsamAll)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Output statistics file '-O<file>' not available in '-M6' output mode\n");
			exit(1);
			}
		strncpy(szStatsFile,statsfile->filename[0],_MAX_PATH);
		szStatsFile[_MAX_PATH-1] = '\0';
		}
	else
		{
		if(bPEInsertLenDist)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: PE insert length distributions requested but no output statistics file '-O<file>' specified");
			exit(1);
			}
		}

	for(NumPE1InputFiles=Idx=0;NumPE1InputFiles < cMaxInFileSpecs && Idx < pe1inputfiles->count; Idx++)
		{
		pszPE1InputFiles[Idx] = NULL;
		if(pszPE1InputFiles[NumPE1InputFiles] == NULL)
			pszPE1InputFiles[NumPE1InputFiles] = new char [_MAX_PATH];
		strncpy(pszPE1InputFiles[NumPE1InputFiles],pe1inputfiles->filename[Idx],_MAX_PATH);
		pszPE1InputFiles[NumPE1InputFiles][_MAX_PATH-1] = '\0';
		CUtility::TrimQuotedWhitespcExtd(pszPE1InputFiles[NumPE1InputFiles]);
		if(pszPE1InputFiles[NumPE1InputFiles][0] != '\0')
			NumPE1InputFiles++;
		}

	if(!NumPE1InputFiles)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: After removal of whitespace, no input file(s) specified with '-i<filespec>' option)\n");
		exit(1);
		}

	if(PEproc != ePEdefault)
		{
		for(NumPE2InputFiles=Idx=0;NumPE2InputFiles < cMaxInFileSpecs && Idx < pe2inputfiles->count; Idx++)
			{
			pszPE2InputFiles[Idx] = NULL;
			if(pszPE2InputFiles[NumPE2InputFiles] == NULL)
				pszPE2InputFiles[NumPE2InputFiles] = new char [_MAX_PATH];
			strncpy(pszPE2InputFiles[NumPE2InputFiles],pe2inputfiles->filename[Idx],_MAX_PATH);
			pszPE2InputFiles[NumPE2InputFiles][_MAX_PATH-1] = '\0';
			CUtility::TrimQuotedWhitespcExtd(pszPE2InputFiles[NumPE2InputFiles]);
			if(pszPE2InputFiles[NumPE2InputFiles][0] != '\0')
				NumPE2InputFiles++;
			}

		if(!NumPE2InputFiles)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: After removal of whitespace, no input file(s) specified with '-u<filespec>' option)\n");
			exit(1);
			}

		if(NumPE1InputFiles != NumPE2InputFiles)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Number (%d) of PE1 files must match number of PE2 (%d) files\n",NumPE1InputFiles,NumPE2InputFiles);
			exit(1);
			}


		PairMinLen = pairminlen->count ? pairminlen->ival[0] : cDfltPairMinLen;
		if(PairMinLen < cPairMinLen || PairMinLen > cPairMaxLen)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: paired end apparent min length '-d%d' must be in range %d..%d\n",PairMinLen,cPairMinLen,cPairMaxLen);
			exit(1);
			}
		PairMaxLen = pairmaxlen->count ? pairmaxlen->ival[0] : max(cDfltPairMaxLen,PairMinLen);
		if(PairMaxLen < PairMinLen || PairMaxLen > cPairMaxLen)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: paired end apparent max length '-D%d' must be in range %d..%d\n",PairMaxLen,PairMinLen,cPairMaxLen);
			exit(1);
			}
		bPairStrand = pairstrand->count ? true : false;
		}
	else
		{
		pszPE2InputFiles[0] = NULL;
		PairMinLen = 0;
		PairMaxLen = 0;
		bPairStrand = 0;
		}

	MinChimericLen = minchimericlen->count ? minchimericlen->ival[0] : 0;
	if(MinChimericLen != 0 && !(abs(MinChimericLen) >= 50 && abs(MinChimericLen) <= 99))
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: minimum chimeric length percentage '-c%d' specified outside of range 50..99\n",abs(MinChimericLen));
		exit(1);
		}

	if(MinChimericLen > 0 && MLMode == eMLdefault)
		bChimericRpt = chimericrpt->count ? true : false;
	else
		bChimericRpt = false;

	if(MinChimericLen && (bBisulfite || bSOLiD))
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Chimeric trimming not supported when either colorspace (SOLiD) '-C' or bisulfite methylation processing\n");
		exit(1);
		}

	bLocateBestMatches = false;
	if(MLMode != eMLdefault)
		{
		MaxMLmatches = maxmlmatches->count ? maxmlmatches->ival[0] : cDfltMaxMultiHits;
		if(MLMode != eMLall)
			{
			if(MaxMLmatches < 2 || MaxMLmatches > cMaxMultiHits)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: multiple aligned reads '-R%d' specified outside of range 2..%d\n",MaxMLmatches,cMaxMultiHits);
				exit(1);
				}
			}
		else   
			{
			if(MaxMLmatches < 2 || MaxMLmatches > cMaxAllHits)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: multiple aligned reads '-R%d' specified outside of range 2..%d\n",MaxMLmatches,cMaxAllHits);
				exit(1);
				}
			}
		bLocateBestMatches = bestmatches->count ? true : false;
		}
	else
		MaxMLmatches = 1;

	if(MaxMLmatches > 1)
		bClampMaxMLmatches = clampmaxmulti->count ? true : false;
	else
		bClampMaxMLmatches = false;
	if(bLocateBestMatches)
		bClampMaxMLmatches = true;


	microInDelLen = microindellen->count ? microindellen->ival[0] : 0;
	if(microInDelLen < 0 || microInDelLen > cMaxMicroInDelLen)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: microInDel length maximum '-a%d' specified outside of range 0..%d\n",microInDelLen,cMaxMicroInDelLen);
		exit(1);
		}

	if(MLMode == eMLall && microInDelLen != 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: microInDels not supported when reporting multiloci alignments");
		exit(1);
		}

	if(abs(MinChimericLen) > 0 && (bSOLiD || bLocateBestMatches))		// chimeric read processing not currently supported for SOLiD colorspace processing or for reporting multiple best matches
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Sorry, chimeric read processing not supported in this release if either SOLiD or locating multiple best matches also requested");
		exit(1);
		}

	
	PCRartefactWinLen = pcrartefactwinlen->count ? pcrartefactwinlen->ival[0] : -1;
	if(PCRartefactWinLen != -1 && (PCRartefactWinLen < 0 || PCRartefactWinLen > 250))
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: PCR differential amplification artefacts window length '-k%d' specified outside of range 0..%d\n",PCRartefactWinLen,250);
		exit(1);
		}

	SpliceJunctLen = splicejunctlen->count ? splicejunctlen->ival[0] : 0;
	if(MLMode != eMLall)
		{
		if(SpliceJunctLen != 0 && (SpliceJunctLen < cMinJunctAlignSep || SpliceJunctLen > cMaxJunctAlignSep))
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: RNAseq maximum splice junction separation '-A%d' must be either 0 or in the range %d..%d\n",SpliceJunctLen,cMinJunctAlignSep,cMaxJunctAlignSep);
			exit(1);
			}
		}
	else
		{
		if(SpliceJunctLen != 0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: in report all multiloci mode '-r5', there is no splice junction processing..\n");
			exit(1);
			}
		}

	Trim5 = trim5->count ? trim5->ival[0] : 0;
	if(Trim5 < 0 || Trim5 > 50)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Trim 5' raw reads '-y%d' specified outside of range %d..%d\n",Trim5,0,50);
		exit(1);
		}
	Trim3 = trim3->count ? trim3->ival[0] : 0;
	if(Trim3 < 0 || Trim3 > 50)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Trim 3' raw reads '-y%d' specified outside of range %d..%d\n",Trim3,0,50);
		exit(1);
		}

	MinAcceptReadLen = minacceptreadlen->count ? minacceptreadlen->ival[0] : cDfltMinAcceptReadLen;
	if(MinAcceptReadLen < cMinSeqLen || MinAcceptReadLen > cMaxSeqLen)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: minimum accepted read length '-l%d' specified outside of range %d..%d\n",Trim5,cMinSeqLen,cMaxSeqLen);
		exit(1);
		}
	MaxAcceptReadLen = maxacceptreadlen->count ? maxacceptreadlen->ival[0] : cDfltMaxAcceptReadLen;
	if(MaxAcceptReadLen < MinAcceptReadLen || MaxAcceptReadLen > cMaxSeqLen)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: maximum accepted read length '-L%d' specified outside of range %d..%d\n",MaxAcceptReadLen,MinAcceptReadLen,cMaxSeqLen);
		exit(1);
		}

	MinEditDist = mineditdist->count ? mineditdist->ival[0] : 1;
	if(MinEditDist < 1 || MinEditDist > 2)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Minimum edit distance '-e%d' specified outside of range %d..%d\n",MinEditDist,1,2);
		exit(1);
		}

	MaxSubs = maxsubs->count ?  maxsubs->ival[0] : cDfltAllowedSubs;
	if(MaxSubs < 0 || MaxSubs > cMaxAllowedSubs)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Max allowed substitutions per 100bp read length '-s%d' specified outside of range %d..%d\n",MaxSubs,0,cMaxAllowedSubs);
		exit(1);
		}

	PCRPrimerCorrect = pcrprimercorrect->count ? pcrprimercorrect->ival[0] : 0;
	if(PCRPrimerCorrect < 0 || PCRPrimerCorrect > cPCRPrimerSubs)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: PCR primer correction subs '-6%d' specified outside of range 0..%d\n",PCRPrimerCorrect,cPCRPrimerSubs);
		exit(1);
		}

	if(PCRPrimerCorrect != 0 && MinChimericLen != 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: PCR primer correction subs not allowed when also specifying chimeric trimming\n");
		exit(1);
		}

	MaxNs = maxns->count ? maxns->ival[0] : cDfltMaxNs;
	if(MaxNs < 0 || MaxNs > cMaxNs)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Allowed number, percentage if read length > 100, of indeterminate bases in reads '-n%d' specified outside of range 0..%d\n",MaxNs,cMaxNs);
		exit(1);
		}

	if(MinFlankExacts < 0 || MinFlankExacts > cMaxAllowedSubs/2)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Max flank trimming '-x%d' specified outside of range 0..%d\n",MinFlankExacts,cMaxAllowedSubs/2);
		exit(1);
		}

	if(SpliceJunctLen > 0 && MinChimericLen == 0 && MinFlankExacts == 0)
		MinFlankExacts = MaxSubs;

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

	if(MLMode == eMLall && !(FMode == eFMdefault || FMode == eFMbed || FMode == eFMsam || FMode == eFMsamAll))
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Output format mode '-M%d' not supported when reporting all multihit read loci\n",FMode);
		exit(1);
		}

	if(bSOLiD && (FMode == eFMread || FMode == eFMreadmatch))
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Output format '-M%d' currently not supported in '-C' colorspace (SOLiD) Processing mode\n",FMode);
		exit(1);
		}

	strcpy(szTargFile,sfxfile->filename[0]);
	strcpy(szRsltsFile,outfile->filename[0]);

	SAMFormat = etSAMFformat;
	if(FMode >= eFMsam)
		{
		// if SAM output format then could be BGZF compressed BAM; use the extension to determine which...
		// if extension is '.bam' then BGZF compressed BAM, any other extension is for SAM
		int Len;
		Len = (int)strlen(szRsltsFile);

		if(Len > 5)
			{
			if(!stricmp(".bam",&szRsltsFile[Len-4]))
				SAMFormat = etSAMFBAM;
			}

		MaxRptSAMSeqsThres = rptsamseqsthres->count ? rptsamseqsthres->ival[0] : cMaxRptSAMSeqsThres;
		if(MaxRptSAMSeqsThres < 1)
			MaxRptSAMSeqsThres = 1;
		}
	else
		MaxRptSAMSeqsThres = 0;

	MarkerLen = 0;
	szMarkerFile[0] = '\0';
	MinSNPreads = 0;
	QValue = 0.0;
	szSNPFile[0] = '\0';
	szSNPCentroidFile[0] = '\0';
	MinSNPreads = minsnpreads->count ? minsnpreads->ival[0] : 0;
	if(MinSNPreads != 0 && (MinSNPreads < cMinSNPreads || MinSNPreads > cMaxSNPreads))
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Minimum read coverage at any loci '-p%d' must be in range %d..%d\n",MinSNPreads,cMinSNPreads,cMaxSNPreads);
		exit(1);
		}
	if(MinSNPreads > 0)
		{
		QValue = qvalue->count ? qvalue->dval[0] : 0.0;
		if(QValue < 0.0 || QValue > 0.40)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: QValue '-P%1.5f' for controlling SNP FDR (Benjamini-Hochberg) must be in range 0.0 to 0.4\n",QValue);
			exit(1);
			}
		if(QValue > 0.0 && MinSNPreads == 0)
			MinSNPreads = cDfltMinSNPreads;
		else
			if(QValue == 0.0 && MinSNPreads != 0)
				QValue = cDfltQValueSNP;

		SNPNonRefPcnt = snpnonrefpcnt->count ? snpnonrefpcnt->dval[0] : 25.0;
		if(SNPNonRefPcnt < 0.1 || SNPNonRefPcnt > 35.0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: SNP minimum non-ref '-1%f' for controlling SNP FDR must be in range 0.1 to 35.0\n",SNPNonRefPcnt);
			exit(1);
			}
		}
	else
		{
		QValue = 0.0;
		SNPNonRefPcnt = 0.0;
		}

	if(MinSNPreads > 0 && (FMode > eFMsam || bBisulfite == true))
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: SNP processing not currently supported if processing bisulfite reads\n");
		exit(1);
		}

	if(MinSNPreads > 0 && MLMode == eMLall)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: SNP processing not currently supported if reporting multiloci alignments\n");
		exit(1);
		}

	if(snpfile->count)
		{
		strncpy(szSNPFile,snpfile->filename[0],_MAX_PATH);
		szSNPFile[_MAX_PATH-1] = '\0';
		}
	else
		{
		strcpy(szSNPFile,szRsltsFile);
		strcat(szSNPFile,".snp");
		}

	MarkerLen = markerlen->count ? markerlen->ival[0] : 0;
	if(MarkerLen != 0 && MarkerLen < cMinMarkerLen || MarkerLen > cMaxMarkerLen)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Marker length specified with '-K%d' must be in range %d to %d",MarkerLen,cMinMarkerLen,cMaxMarkerLen);
		exit(1);
		}
	if(MarkerLen)
		{
		strcpy(szMarkerFile,szSNPFile);
		strcat(szMarkerFile,".markers");
		}
	else
		szMarkerFile[0] ='\0';

	if(MarkerLen)
		{
		MarkerPolyThres = markerpolythres->count ? markerpolythres->dval[0] : cDfltMinMarkerSNPProp;
		if(MarkerPolyThres < 0.0 || MarkerPolyThres > 0.50)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Max marker sequence base poymorphism specified with '-G%1.3f' must be in range 0.0 to 0.5",MarkerPolyThres);
			exit(1);
			}
		}
	else
		MarkerPolyThres = 0.0;

	if(centroidfile->count)
		{
		strncpy(szSNPCentroidFile,centroidfile->filename[0],_MAX_PATH);
		szSNPCentroidFile[_MAX_PATH-1] = '\0';
		}
	else
		szSNPCentroidFile[0] = '\0';

	szTrackTitle[0] = '\0';

	if(FMode == eFMbed)
		{
		if(title->count)
			{
			strncpy(szTrackTitle,title->sval[0],sizeof(szTrackTitle));
			szTrackTitle[sizeof(szTrackTitle)-1] = '\0';
			CUtility::TrimQuotedWhitespcExtd(szTrackTitle);
			CUtility::ReduceWhitespace(szTrackTitle);
			}
		if(szTrackTitle[0] == '\0')
			{
			gDiagnostics.DiagOut(eDLWarn,gszProcName,"Warning: output format requested to be UCSC BED but no track title with '-t<title' specified, defaulting to 'kanga'\n");
			strcpy(szTrackTitle,"kanga");
			}
		}

	NumIncludeChroms = IncludeChroms->count;
	for(Idx=0;Idx < IncludeChroms->count; Idx++)
		{
		LenChromList = (int)strlen(IncludeChroms->sval[Idx]);
		pszIncludeChroms[Idx] = new char [LenChromList+1];
		strcpy(pszIncludeChroms[Idx],IncludeChroms->sval[Idx]);
		CUtility::TrimQuotes(pszIncludeChroms[Idx]);
		}

	NumExcludeChroms = ExcludeChroms->count;
	for(Idx=0;Idx < ExcludeChroms->count; Idx++)
		{
		LenChromList = (int)strlen(ExcludeChroms->sval[Idx]);
		pszExcludeChroms[Idx] = new char [LenChromList+1];
		strcpy(pszExcludeChroms[Idx],ExcludeChroms->sval[Idx]);
		CUtility::TrimQuotes(pszExcludeChroms[Idx]);
		}

	if(lociconstraintsfile->count)
		{
		strncpy(szLociConstraintsFile,lociconstraintsfile->filename[0],_MAX_PATH);
		szLociConstraintsFile[_MAX_PATH-1] = '\0';
		}
	else
		szLociConstraintsFile[0] = '\0';

	if(contamsfile->count)
		{
		strncpy(szContamFile,contamsfile->filename[0],_MAX_PATH);
		szContamFile[_MAX_PATH-1] = '\0';
		CUtility::TrimQuotes(szContamFile);
		}
	else
		szContamFile[0] = '\0';

	if(statsfile->count)
		{
		if(FMode == eFMsamAll)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Output induced substitution mode '-O<file>' not available in '-M6' output mode\n");
			exit(1);
			}
		strncpy(szStatsFile,statsfile->filename[0],_MAX_PATH);
		szStatsFile[_MAX_PATH-1] = '\0';
		CUtility::TrimQuotes(szStatsFile);
		}
	else
		szStatsFile[0] = '\0';
	
	SitePrefsOfs = siteprefsofs->count ? siteprefsofs->ival[0] : cDfltRelSiteStartOfs;
	if(abs(SitePrefsOfs) > cMaxSitePrefOfs)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: offset read start sites '-9%d' when processing site octamer preferencing must be in range -100..100\n",SitePrefsOfs);
		exit(1);
		}

	if(siteprefsfile->count)
		{
		strncpy(szSitePrefsFile,siteprefsfile->filename[0],_MAX_PATH);
		szSitePrefsFile[_MAX_PATH-1] = '\0';
		CUtility::TrimQuotes(szSitePrefsFile);
		}
	else
		szSitePrefsFile[0] = '\0';

	if(MLMode != eMLall)
		{
		if(nonealignfile->count)
			{
			strncpy(szNoneAlignFile,nonealignfile->filename[0],_MAX_PATH);
			szNoneAlignFile[_MAX_PATH-1] = '\0';
			CUtility::TrimQuotes(szNoneAlignFile);
			}
		else
			szNoneAlignFile[0] = '\0';

		if(multialignfile->count)
			{
			strncpy(szMultiAlignFile,multialignfile->filename[0],_MAX_PATH);
			szMultiAlignFile[_MAX_PATH-1] = '\0';
			CUtility::TrimQuotes(szMultiAlignFile);
			}
		else
			szMultiAlignFile[0] = '\0';
		}
	else
		{
		szNoneAlignFile[0] = '\0';
		szMultiAlignFile[0] = '\0';
		}

	if(priorityregionfile->count)
		{
		strncpy(szPriorityRegionFile,priorityregionfile->filename[0],_MAX_PATH);
		szPriorityRegionFile[_MAX_PATH-1] = '\0';
		CUtility::TrimQuotes(szPriorityRegionFile);
		bFiltPriorityRegions = nofiltpriority->count ? false : true;
		}
	else
		{
		szPriorityRegionFile[0] = '\0';
		bFiltPriorityRegions = false;
		}

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing parameters:");

	const char *pszDescr;
	switch(PMode) {
		case ePMdefault:
			pszDescr = "Standard alignment sensitivity";
			break;
		case ePMMoreSens:
			pszDescr = "More sensitive alignment (slower)";
			break;
		case ePMUltraSens:
			pszDescr = "Ultra sensitive alignment (very slow)";
			break;
		case ePMLessSens:
		default:
			pszDescr = "Less sensitive alignment (quicker)";
			break;
		}

// show user current resource limits
#ifndef _WIN32
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Resources: %s",CUtility::ReportResourceLimits());
#endif

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Processing mode is : '%s'",pszDescr);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Raw read or paired reads sampling is : every %u",SampleNthRawRead);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Processing in %s mode",bSOLiD ? "colorspace (SOLiD)" : "standard basespace");

	if(bBisulfite)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Processing for: %s",bBisulfite ? "bisulfite methylation patterning" : "standard");

	switch(AlignStrand) {
		case eALSboth:
			pszDescr = "Either sense '+' or antisense '-' strands";
			break;
		case eALSWatson:
			pszDescr = "Sense '+' strand only";
			break;
		case eALSCrick:
			pszDescr = "Antisense '-' strand only";
			break;
		}

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"alignments are to : %s",pszDescr);


	if(PCRartefactWinLen >= 0)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Use this window length when reducing PCR differential amplification artefacts : %d",PCRartefactWinLen);
	else
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"No PCR differential amplification artefact reduction");
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"trim 5' ends raw reads by : %d",Trim5);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"trim 3' ends raw reads by : %d",Trim3);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"after any end trimming reads must be at least this length : %dbp",MinAcceptReadLen);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"after any end trimming reads must be no longer than this length : %dbp",MaxAcceptReadLen);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"maximum aligner induced substitutions : %d subs per 100bp of actual read length",MaxSubs);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"minimum Hamming edit distance : %d",MinEditDist);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"maximum number, precentage of length if read length > 100, of indeterminate 'N's : %d",MaxNs);
	if(MinFlankExacts >= 0)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"minimum 5' and 3' flank exacts : %d",MinFlankExacts);
	if(PCRPrimerCorrect > 0)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"correcting PCR 5' artefacts : %d",PCRPrimerCorrect);

	switch(FMode) {
		case eFMdefault:
			pszDescr = "CSV match loci only";
			break;
		case eFMmatch:
			pszDescr = "CSV loci + match sequence";
			break;
		case eFMread:
			pszDescr = "CSV loci + read sequence";
			break;
		case eFMreadmatch:
			pszDescr = "CSV loci + read + match sequence";
			break;
		case eFMbed:
			pszDescr = "UCSC BED";
			break;
		case eFMsam:
			pszDescr = "SAM Toolset Format, accepted aligned reads only";
			break;
		case eFMsamAll:
			pszDescr = "SAM toolset format, includes all reads";
			break;
		}

	const char *pszProcMode;
	switch(Quality) {
		case eFQSanger:
			pszProcMode = "Sanger or Illumina 1.8+ Phred";
			break;

		case eFQIllumia:
			pszProcMode = "Illumina 1.3+ Phred";
			break;

		case eFQSolexa:
			pszProcMode = "Solexa/Illumia pre-1.3";
			break;

		case eFQIgnore:
			pszProcMode = "Ignore";
			break;

		default:
			pszProcMode = "Unknown quality scoring, defaulting to ignore";
			Quality = eFQIgnore;
			break;
		};

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Raw read quality scores are : '%s'",pszProcMode);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"output format is : '%s'",pszDescr);
	switch(SAMFormat) {
		case etSAMFBAM:					// output as uncompressed BAM
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"SAM will be generated as BGZF compressed BAM");
			break;
		default:
			break;
		}

	if(MaxRptSAMSeqsThres > 0)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"If number of target sequences no more than this threshold then write all sequence names to SAM header: %d",MaxRptSAMSeqsThres);

	if(szTrackTitle[0] != '\0')
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"track title: '%s'",szTrackTitle);
	for(Idx=0; Idx < NumPE1InputFiles; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"input PE1 raw reads files (%d): '%s'",Idx+1,pszPE1InputFiles[Idx]);

	switch(PEproc) {
		case ePEdefault:
			pszProcMode = "Single ended reads";
			break;

		case ePEorphan:					// process for paired ends and if one end is orphaned because multiple aligned then try to locate unique alignment downstream
			pszProcMode = "Paired end reads with orphan partner processing";
			break;

		case ePEunique:					// process for paired ends but only accept putative if both ends uniquely aligned
			pszProcMode = "Paired end reads with both ends uniquely aligned within the targeted genome";
			break;

		case ePEorphanSE:					// process for paired ends and if one end is orphaned because multiple aligned then try to locate unique alignment downstream
			pszProcMode = "Paired end reads with orphan partner processing treating orphans as SE";
			break;

		case ePEuniqueSE:					// process for paired ends but only accept putative if both ends uniquely aligned
			pszProcMode = "Paired end reads with both ends uniquely aligned within the targeted genome treating orphans as SE";
			break;

		}

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"process for: '%s'",pszProcMode);

	for(Idx=0; Idx < NumPE2InputFiles; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"input PE2 raw reads files (%d): '%s'",Idx+1,pszPE2InputFiles[Idx]);

	if(szPriorityRegionFile[0] != '\0')
		{
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Prioritise exact matchs to regions in this BED file: '%s'",szPriorityRegionFile);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Filter out matchs in non-prioritorised regions: '%s'",bFiltPriorityRegions ? "Yes" : "No");
		}

	if(PEproc > ePEdefault)
		{
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Accept as paired if observed insert size is between %d and %d",PairMinLen,PairMaxLen);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Accept as paired if 5' and 3' are are same strand: '%s'",bPairStrand ? "Yes" : "No");
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Experimental: Output PE insert length distributions for each transcript or contig : '%s'",bPEInsertLenDist ? "Yes" : "No");
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Experimental: Process PEs for spanning of circularised fragments: '%s'",bPEcircularised ? "Yes" : "No");
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Output paired end sequence length distribution to file: '%s'",szStatsFile[0] == '\0' ? "none specified" : szStatsFile);
		}
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"input target sequence(s) suffix array file: '%s'",szTargFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"output results file: '%s'",szRsltsFile);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Output none-aligned reads to fasta file: '%s'",szNoneAlignFile[0] == '\0' ? "none specified" : szNoneAlignFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Output reads with multiple alignments to fasta file: '%s'",szMultiAlignFile[0] == '\0' ? "none specified" : szMultiAlignFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Output aligner induced substitution distributions to file: '%s'",szStatsFile[0] == '\0' ? "none specified" : szStatsFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Loci base constraints file: '%s'",szLociConstraintsFile[0] == '\0' ? "none specified" : szLociConstraintsFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Contaminant sequences file: '%s'",szContamFile[0] == '\0' ? "none specified" : szContamFile);


	switch(MLMode) {
		case eMLdefault:
			pszDescr = "slough all reads which match to multiple loci";
			break;
		case eMLdist:
			pszDescr = "accumulate distribution stats only, slough reads matching to multiple loci";
			break;
		case eMLrand:
			pszDescr = "randomly select one of the aligned loci";
			break;
		case eMLuniq:
			pszDescr = "cluster with reads which are uniquely aligned";
			break;
		case eMLcluster:
			pszDescr = "cluster with unique (high priority) and other multiloci alignments";
			break;
		case eMLall:
			pszDescr = "report all multiple loci to which reads align";
			break;
		}
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Process multiple alignment reads by: '%s'",pszDescr);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Offset read start sites when processing site octamer preferencing: %d",SitePrefsOfs);

	if(szSitePrefsFile[0] != '\0')
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Aligned read octamer site preferencing into this file: '%s'",szSitePrefsFile);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Allow microInDels of upto this inclusive length: %d",microInDelLen);

	if(MinChimericLen != 0)
		{
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Check for chimeric sequences in reads of at least this percentage length: %d",abs(MinChimericLen));
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Report on individual read sequence chirmeric trimming: %s",bChimericRpt ? "Yes" : "no");
		}

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Maximum RNA-seq splice junction separation distance: %d",SpliceJunctLen);
	if(MinSNPreads == 0)
		{
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Minimum read coverage at loci before processing for SNP: No SNP processing");
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"QValue controlling FDR (Benjamini-Hochberg) SNP prediction : No SNP processing");
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Min percentage non-ref bases at putative SNP loci : No SNP processing");
		
		}
	else
		{
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Minimum read coverage at loci before processing for SNP: %d",MinSNPreads);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"QValue controlling FDR (Benjamini-Hochberg) SNP prediction : %0.5f",QValue);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Min percentage non-ref bases at putative SNP loci : %1.3f",SNPNonRefPcnt);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"SNP predictions written to file : '%s'",szSNPFile);
		if(szSNPCentroidFile[0] != '\0')
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"SNP prediction centroid distributions written to file : '%s'",szSNPCentroidFile);
		if(MarkerLen)
			{
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Marker sequence length : %d",MarkerLen);
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Max marker sequence base polymorphism : %1.3f",MarkerPolyThres);
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Marker sequences written to file : '%s'",szMarkerFile);
			}
		}

	if(MaxMLmatches == 1)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Only accept reads which uniquely match a single loci");
	else
		{
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Allow at most any read to match this many loci and then process: %d",MaxMLmatches);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"If read aligns to more than %d loci then treat as if aligned to the first %d loci discovered : %s",MaxMLmatches, MaxMLmatches,bClampMaxMLmatches ? "Yes" : "No");
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Report best alignments: '%s'",bLocateBestMatches ? "Yes" : "No");
		}

	for(Idx = 0; Idx < NumIncludeChroms; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"reg expressions defining chroms to include: '%s'",pszIncludeChroms[Idx]);
	for(Idx = 0; Idx < NumExcludeChroms; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"reg expressions defining chroms to exclude: '%s'",pszExcludeChroms[Idx]);

	if(szExperimentName[0] != '\0')
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"This processing reference: %s",szExperimentName);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"number of threads : %d",NumThreads);

	if(gExperimentID > 0)
		{
		int ParamID;
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szLogFile),"log",szLogFile);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(PMode),"mode",&PMode);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(SampleNthRawRead),"samplenthrawread",&SampleNthRawRead);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(FMode),"format",&FMode);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(PEproc),"pemode",&PEproc);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTBool,(int)sizeof(bBisulfite),"bisulfite",&bBisulfite);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTBool,(int)sizeof(bSOLiD),"colorspace",&bSOLiD);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(PCRartefactWinLen),"pcrwin",&PCRartefactWinLen);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTBool,(int)sizeof(bFiltPriorityRegions),"nofiltpriority",&bFiltPriorityRegions);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(PairMinLen),"pairminlen",&PairMinLen);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(PairMaxLen),"pairmaxlen",&PairMaxLen);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTBool,(int)sizeof(bPairStrand),"pairstrand",&bPairStrand);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTBool,(int)sizeof(bPEcircularised),"pecircularised",&bPEcircularised);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTBool,(int)sizeof(bPEInsertLenDist),"peinsertlendist",&bPEInsertLenDist);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(AlignStrand),"alignstrand",&AlignStrand);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(Quality),"quality",&Quality);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(microInDelLen),"microindellen",&microInDelLen);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(MinChimericLen),"minchimericlen",&MinChimericLen);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(SpliceJunctLen),"splicejunctlen",&SpliceJunctLen);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(SitePrefsOfs),"siteprefsofs",&SitePrefsOfs);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(MinEditDist),"editdelta",&MinEditDist);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(MaxSubs),"substitutions",&MaxSubs);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(MinFlankExacts),"minflankexacts",&MinFlankExacts);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTBool,(int)sizeof(PCRPrimerCorrect),"pcrprimercorrect",&PCRPrimerCorrect);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(MaxNs),"maxns",&MaxNs);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(MaxMLmatches),"maxmulti",&MaxMLmatches);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTBool,(int)sizeof(bClampMaxMLmatches),"clampmaxmulti",&bClampMaxMLmatches);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(MLMode),"mlmode",&MLMode);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(MinSNPreads),"snpreadsmin",&MinSNPreads);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTDouble,(int)sizeof(QValue),"qvalue",&QValue);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTDouble,(int)sizeof(SNPNonRefPcnt),"snpnonrefpcnt",&SNPNonRefPcnt);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(MaxRptSAMSeqsThres),"rptsamseqsthres",&MaxRptSAMSeqsThres);

		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTBool,(int)sizeof(bLocateBestMatches),"bestmatches",&bLocateBestMatches);

		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(MarkerLen),"markerlen",&MarkerLen);
		
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTDouble,(int)sizeof(MarkerPolyThres),"markerpolythres",&MarkerPolyThres);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szMarkerFile),"markerfile",szMarkerFile);

		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(Trim5),"trim5",&Trim5);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(Trim3),"trim3",&Trim3);

		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(MinAcceptReadLen),"minacceptreadlen",&MinAcceptReadLen);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(MaxAcceptReadLen),"maxacceptreadlen",&MaxAcceptReadLen);
		
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szTrackTitle),"title",szTrackTitle);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(NumIncludeChroms),"NumIncludeChroms",&NumIncludeChroms);
		for(Idx = 0; Idx < NumIncludeChroms; Idx++)
			ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(pszIncludeChroms[Idx]),"chromeinclude",pszIncludeChroms[Idx]);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(NumExcludeChroms),"NumExcludeChroms",&NumExcludeChroms);
		for(Idx = 0; Idx < NumExcludeChroms; Idx++)
			ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(pszExcludeChroms[Idx]),"chromeexclude",pszExcludeChroms[Idx]);

		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(NumPE1InputFiles),"NumPE1InputFiles",&NumPE1InputFiles);
		for(Idx=0; Idx < NumPE1InputFiles; Idx++)
			ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(pszPE1InputFiles[Idx]),"in",pszPE1InputFiles[Idx]);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(NumPE2InputFiles),"NumPE2InputFiles",&NumPE2InputFiles);
		for(Idx=0; Idx < NumPE2InputFiles; Idx++)
			ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(pszPE2InputFiles[Idx]),"pair",pszPE2InputFiles[Idx]);

		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szPriorityRegionFile),"priorityregionfile",szPriorityRegionFile);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szTargFile),"sfx",szTargFile);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szRsltsFile),"out",szRsltsFile);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szStatsFile),"stats",szStatsFile);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szNoneAlignFile),"nonealign",szNoneAlignFile);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szMultiAlignFile),"multialign",szMultiAlignFile);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szSitePrefsFile),"siteprefs",szSitePrefsFile);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szSNPFile),"snpfile",szSNPFile);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szSNPCentroidFile),"snpcentroid",szSNPCentroidFile);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szLociConstraintsFile),"lociconstraints",szLociConstraintsFile);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szContamFile),"contamsfile",szContamFile);

		
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(NumThreads),"threads",&NumThreads);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(NumberOfProcessors),"cpus",&NumberOfProcessors);

		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szSQLiteDatabase),"sumrslts",szSQLiteDatabase);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szExperimentName),"experimentname",szExperimentName);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szExperimentDescr),"experimentdescr",szExperimentDescr);
		}


#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	gStopWatch.Start();
	Rslt = Process((etPMode)PMode,SampleNthRawRead,(etFQMethod)Quality,bSOLiD,bBisulfite,(etPEproc)PEproc,PairMinLen,PairMaxLen,bPairStrand,bPEcircularised,bPEInsertLenDist,
				    (eALStrand)AlignStrand,MinChimericLen,bChimericRpt,microInDelLen,SpliceJunctLen,
					MinSNPreads,QValue,SNPNonRefPcnt,MarkerLen,MarkerPolyThres,PCRartefactWinLen,(etMLMode)MLMode,
					MaxMLmatches,bClampMaxMLmatches,bLocateBestMatches,
					MaxNs,MinEditDist,MaxSubs,Trim5,Trim3,MinAcceptReadLen,MaxAcceptReadLen,MinFlankExacts,PCRPrimerCorrect, MaxRptSAMSeqsThres,
					(etFMode)FMode,SAMFormat,SitePrefsOfs,NumThreads,szTrackTitle,
					NumPE1InputFiles,pszPE1InputFiles,NumPE2InputFiles,pszPE2InputFiles,szPriorityRegionFile,bFiltPriorityRegions,szRsltsFile, szSNPFile, szMarkerFile, szSNPCentroidFile, szTargFile,
					szStatsFile,szMultiAlignFile,szNoneAlignFile,szSitePrefsFile,szLociConstraintsFile,szContamFile,NumIncludeChroms,pszIncludeChroms,NumExcludeChroms,pszExcludeChroms);
	Rslt = Rslt >=0 ? 0 : 1;
	if(gExperimentID > 0)
		{
		if(gProcessingID)
			gSQLiteSummaries.EndProcessing(gProcessingID,Rslt);
		gSQLiteSummaries.EndExperiment(gExperimentID);
		}
	gStopWatch.Stop();

#ifdef _DEBUG
#ifdef _WIN32
_ASSERTE( _CrtCheckMemory());
#endif
#endif

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Exit code: %d Total processing time: %s",Rslt,gStopWatch.Read());
	exit(Rslt);
	}
else
	{
    printf("\n%s %s %s, Version %s\n", gszProcName,gpszSubProcess->pszName,gpszSubProcess->pszFullDescr,cpszProgVer);
	arg_print_errors(stdout,end,gszProcName);
	arg_print_syntax(stdout,argtable,"\nUse '-h' to view option and parameter usage\n");
	exit(1);
	}
return 0;
}


int	(CSfxArrayV3::*m_pIterateExactsFn)(etSeqBase *,unsigned int,unsigned int,unsigned int *,unsigned int *);


int
Process(etPMode PMode,					// processing mode
		UINT32 SampleNthRawRead,		// sample every Nth raw read (or read pair) for processing (1..10000)
		etFQMethod Quality,				// quality scoring for fastq sequence files
		bool bSOLiD,					// if true then processing in colorspace
		bool bBisulfite,				// if true then process for bisulfite methylation patterning
		etPEproc PEproc,				// paired reads alignment processing mode
		int PairMinLen,					// accept paired end alignments with observed insert size of at least this (default = 100)
		int PairMaxLen,					// accept paired end alignments with observed insert size of at most this (default = 1000)
		bool bPairStrand,				// accept paired ends if on same strand
		bool bPEcircularised,			// experimental - true if processing for PE spaning circularised fragments
		bool bPEInsertLenDist,			// true if stats file to include PE insert length distributions for each transcript
		eALStrand AlignStrand,			// align on to watson, crick or both strands of target
		int MinChimericLen,				// minimum chimeric length as a percentage (0 to disable, otherwise 50..99) of probe sequence
		bool bChimericRpt,				// report chimeric trimming detail for individual reads (default is not to report)
		int microInDelLen,				// microInDel length maximum
		int SpliceJunctLen,				// maximum splice junction length when aligning RNAseq reads
		int MinSNPreads,				// must be at least this number of reads covering any loci before processing for SNPs at this loci
		double QValue,					// QValue controlling FDR (Benjamini–Hochberg) SNP prediction
		double SNPNonRefPcnt,			// only process for SNP if more/equal than this percentage number of reads are non-ref at putative SNP loci (defaults to 25.0) 
		int MarkerLen,					// marker sequences of this length
		double MarkerPolyThres,			// maximum allowed marker sequence base polymorphism independent of centroid SNP (default 0.333, range 0.0 to 0.5)
		int PCRartefactWinLen,			// if >= 0 then window size to use when attempting to reduce the number of  PCR differential amplification artefacts (reads stacking to same loci)
		etMLMode MLMode,				// multimatch loci reads processing mode
		int MaxMLmatches,				// accept at most this number of multimatched alignments for any read
		bool bClampMaxMLmatches,		// accept as if MaxMLmatches even if more than this number of MaxMLmatches multimached alignments
		bool bLocateBestMatches,		// align for best matches, not just those 1H better than any other, upto at most MaxMLmatches
		int MaxNs,					    // allow at most this number of indeterminate eBaseNs in read before deeming as nonalignable
		int MinEditDist,				// any matches must have at least this edit distance to the next best match
		int MaxSubs,					// maximum number of substitutions allowed per 100bp of actual read length
		int Trim5,						// trim this number of bases from 5' end of reads when loading the reads
		int Trim3,						// trim this number of bases from 3' end of reads when loading the reads
		int MinAcceptReadLen,			// only accepting reads for alignment if at least this length after any end trimming
		int MaxAcceptReadLen,			// only accepting reads for alignment if no longer than this length after any end trimming
		int MinFlankExacts,				// trim matched reads on 5' and 3' flanks until at least this number of exactly matching bases in flanks
		int PCRPrimerCorrect,			// initially align with MaxSubs+PCRPrimerCorrect subs allowed but then attempt to correct 5' PCR primer artefacts until read within MaxSubs
		int MaxRptSAMChroms,			// report all SAM chroms or sequences if number of reference chroms <= this limit (defaults to 10000)
		etFMode FMode,					// output format mode
		teSAMFormat SAMFormat,			// if SAM output format then could be SAM, BAM or BAM compressed dependent on the file extension used
		int SitePrefsOfs,				// offset read start sites when processing  octamer preferencing, range -100..100
		int NumThreads,					// number of worker threads to use
		char *pszTrackTitle,			// track title if output format is UCSC BED
		int NumPE1InputFiles,			// number of input PE1 or single ended file specs
		char *pszPE1InputFiles[],		// names of input files (wildcards allowed unless processing paired ends) containing raw reads
		int NumPE2InputFiles,			// number of input PE2 file specs
		char *pszPE2InputFiles[],		// optional raw paired reads are in these files
		char *pszPriorityRegionFile,	// optional high priority BED file contains prioritised region loci
		bool bFiltPriorityRegions,		// true if non-priority alignments to be filtered out 
		char *pszOutFile,				// where to write alignments
		char *pszSNPFile,				// Output SNPs (CSV format) to this file (default is to output file name with '.snp' appended)
		char *pszMarkerFile,			// Output markers to this file
		char *pszSNPCentroidFile,		// Output SNP centorids (CSV format) to this file (default is for no centroid processing)
		char *pszSfxFile,				// target as suffix array
		char *pszStatsFile,				// aligner induced substitutions stats file
		char *pszMultiAlignFile,		// file to contain reads which are aligned to multiple locations
		char *pszNoneAlignFile,			// file to contain reads which were non-alignable
		char *pszSitePrefsFile,			// file to contain aligned reads octamer preferencing
		char *pszLociConstraintsFile,	// loci base constraints file
		char *pszContamFile,			// file containing contaminant sequences
		int	NumIncludeChroms,			// number of chromosome regular expressions to include
		char **ppszIncludeChroms,		// array of include chromosome regular expressions
		int	NumExcludeChroms,			// number of chromosome expressions to exclude
		char **ppszExcludeChroms)		// array of exclude chromosome regular expressions
{
int Rslt;
CAligner *pAligner;

if((pAligner = new CAligner)==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: Unable to instantiate CAligner");
	return(eBSFerrObj);
	}

Rslt = pAligner->Align(PMode,			// processing mode
			SampleNthRawRead,			// sample every Nth raw read (or read pair) for processing (1..10000)
			Quality,					// quality scoring for fastq sequence files
			bSOLiD,						// if true then processing in colorspace
			bBisulfite,					// if true then process for bisulfite methylation patterning
			PEproc,						// paired reads alignment processing mode
			PairMinLen,					// accept paired end alignments with observed insert size of at least this (default = 100)
			PairMaxLen,					// accept paired end alignments with observed insert size of at most this (default = 1000)
			bPairStrand,				// accept paired ends if on same strand
			bPEcircularised,			// experimental - true if processing for PE spaning circularised fragments
			bPEInsertLenDist,			// experimental - true if stats file to include PE insert length distributions for each transcript
			AlignStrand,				// align on to watson, crick or both strands of target
			MinChimericLen,				// minimum chimeric length as a percentage (0 to disable, otherwise 50..99) of probe sequence
			bChimericRpt,				// report chimeric trimming detail for individual reads (default is not to report)
			microInDelLen,				// microInDel length maximum
			SpliceJunctLen,				// maximum splice junction length when aligning RNAseq reads
			MinSNPreads,				// must be at least this number of reads covering any loci before processing for SNPs at this loci
			QValue,						// QValue controlling FDR (Benjamini–Hochberg) SNP prediction
			SNPNonRefPcnt,				// only process for SNP if more/equal than this percentage number of reads are non-ref at putative SNP loci (defaults to 25) 
			MarkerLen,					// marker sequences of this length
			MarkerPolyThres,			// maximum allowed marker sequence base polymorphism independent of centroid SNP (default 0.333, range 0.0 to 0.5)
			PCRartefactWinLen,			// if >= 0 then window size to use when attempting to reduce the number of  PCR differential amplification artefacts (reads stacking to same loci)
			MLMode,						// multimatch loci reads processing mode
			MaxMLmatches,				// accept at most this number of multimatched alignments for any read
			bClampMaxMLmatches,			// accept as if MaxMLmatches even if more than this number of MaxMLmatches multimached alignments
			bLocateBestMatches,			// align for best matches, not just those 1H better than any other, upto at most MaxMLmatches
			MaxNs,					    // allow at most this number of indeterminate eBaseNs in read before deeming as nonalignable
			MinEditDist,				// any matches must have at least this edit distance to the next best match
			MaxSubs,					// maximum number of substitutions allowed per 100bp of actual read length
			Trim5,						// trim this number of bases from 5' end of reads when loading the reads
			Trim3,						// trim this number of bases from 3' end of reads when loading the reads
			MinAcceptReadLen,			// only accepting reads for alignment if at least this length after any end trimming
			MaxAcceptReadLen,			// only accepting reads for alignment if no longer than this length after any end trimming
			MinFlankExacts,				// trim matched reads on 5' and 3' flanks until at least this number of exactly matching bases in flanks
			PCRPrimerCorrect,			// initially align with MaxSubs+PCRPrimerCorrect subs allowed but then attempt to correct 5' PCR primer artefacts until read within MaxSubs
			MaxRptSAMChroms,			// report all SAM chroms or sequences if number of reference chroms <= this limit (defaults to 10000)
			FMode,						// output format mode
			SAMFormat,					// if SAM output format then could be SAM, BAM or BAM compressed dependent on the file extension used
			SitePrefsOfs,				// offset read start sites when processing  octamer preferencing, range -100..100
			NumThreads,					// number of worker threads to use
			pszTrackTitle,				// track title if output format is UCSC BED
			NumPE1InputFiles,			// number of input PE1 or single ended file specs
			pszPE1InputFiles,			// names of input files (wildcards allowed unless processing paired ends) containing raw reads
			NumPE2InputFiles,			// number of input PE2 file specs
			pszPE2InputFiles,			// optional raw paired reads are in these files
			pszPriorityRegionFile,		// optional high priority BED file contains prioritised region loci
			bFiltPriorityRegions,		// true if non-priority alignments to be filtered out 
			pszOutFile,					// where to write alignments
			pszSNPFile,					// Output SNPs (CSV format) to this file (default is to output file name with '.snp' appended)
			pszMarkerFile,				// Output markers to this file
			pszSNPCentroidFile,			// Output SNP centorids (CSV format) to this file (default is for no centroid processing)
			pszSfxFile,					// target as suffix array
			pszStatsFile,				// aligner induced substitutions stats file
			pszMultiAlignFile,			// file to contain reads which are aligned to multiple locations
			pszNoneAlignFile,			// file to contain reads which were non-alignable
			pszSitePrefsFile,			// file to contain aligned reads octamer preferencing
			pszLociConstraintsFile,		// loci base constraints file
			pszContamFile,				// file containing contaminant sequences
			NumIncludeChroms,			// number of chromosome regular expressions to include
			ppszIncludeChroms,			// array of include chromosome regular expressions
			NumExcludeChroms,			// number of chromosome expressions to exclude
			ppszExcludeChroms);			// array of exclude chromosome regular expressions);

delete pAligner;
return(Rslt);
}