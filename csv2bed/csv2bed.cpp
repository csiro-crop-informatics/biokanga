// csv2bed.cpp : Defines the entry point for the console application.
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

const unsigned int cProgVer = 306;		// increment with each release

const unsigned int cLineBuffSize = 10000; // lines are buffered up until this length before write to disk
const int cMinLenEl = 0;				// minimum length element goes down to 0
const int cDfltMinLenEl = 20;			// default minimum if none specified
const int cMaxLenEl = 100000000;		// max length element
const int cDfltMaxLenEl = 100000000;	// default max length if none specified

CStopWatch gStopWatch;
CDiagnostics gDiagnostics;				// for writing diagnostics messages to log file
char gszProcName[_MAX_FNAME];			// process name

int
Process(int MinLen,					// filter out elements of less than this length
		int MaxLen,					// or longer than this length
		char *pszInfile,			// csv file containing element loci
		char *pszFilterRefIDFile,	// RefIDs to filter
		char *pszTrackName,			// optional track name; if track name and description specified specified then BED header will be created
		char *pszDescr,				// optional track description; if track name and description specified specified then BED header will be created
		bool bAppend,				// false to create/truncate existing file, true to append
		char *pszOutfile);			// BED file to create


int
ProcessFindApproxMatches(int MinLen,// filter out elements of less than this length
		int MaxLen,					// or longer than this length
		char *pszInfile,			// csv file containing element loci
		char *pszFilterRefIDFile,	// RefIDs to filter
		char *pszOutfile);			// BED file to create

int
ProcessGenhyperconservedX2(int MinLen,// filter out core elements of less than this length
		int MaxLen,					// or longer than this length
		int MinIdentity,			// filter out if out species has less than this identity (0..100)
		int MaxIdentity,			// filter out if out species has more than this identity (MinIdentity..100)
		char *pszInfile,			// csv file containing element loci
		char *pszFilterRefIDFile,	// RefIDs to filter
		char *pszOutfile);			// BED file to create

int
ProcessGenMAlignScoreM3(int MinLen,	// filter out elements of less than this length
		int MaxLen,					// or longer than this length
		char *pszInfile,			// csv file containing element loci
		char *pszFilterRefIDFile,	// RefIDs to filter
		char *pszOutfile);			// BED file to create


int
ProcessPSL2CSV(int MinLen,			// filter out elements of less than this length
		int MaxLen,					// or longer than this length
		int MinInstances,			// must be at least this number of instances
		int MaxInstances,			// must be no more than this number of instances
		char *pszInfile,			// csv file containing element loci
		char *pszFilterRefIDFile,	// RefIDs to filter
		char *pszOutfile);			// BED file to create


int
ProcessBLAST2CSV(int MinLen,		// filter out elements of less than this length
		int MaxLen,					// or longer than this length
		int MinInstances,			// must be at least this number of instances
		int MaxInstances,			// must be no more than this number of instances
		char *pszInfile,			// csv file containing element loci
		char *pszFilterRefIDFile,	// RefIDs to filter
		char *pszOutfile);			// BED file to create


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
char szInfile[_MAX_PATH];
char szFilterRefIDFile[_MAX_PATH];  // exclude any RefIDs in this filter file
char szOutfile[_MAX_PATH];
char szDescr[80];
char szTrack[80];
char szChrom[80];
char szName[80];
int iEntry,iEntryfld,iStepfld,iNbins,iStartstep,iEndstep,iValuefield,iWindowsize,iRemap;
bool bAppend, bWiggle;
int iMode;

CCSV2BED *pCsv2Bed = NULL;
int iMinLenEl;			// core elements must be of at least this length
int iMaxLenEl;			// and no longer than this length
int iMinIdentity;		// out species must align with at least this identity
int iMaxIdentity;		// and no more than this identity

int iMinInstances;		// minimum number of instances with same RefID in mode 6
int iMaxInstances;		// maximum number of instances with same RefID in mode 6


// command line args
struct arg_lit  *help    = arg_lit0("hH","help",                "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_int *ScreenLogLevel=arg_int0("S", "ScreenLogLevel",	"<int>","Level of diagnostics written to screen 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",	"diagnostics log file");

struct arg_int *Mode = arg_int0("m", "mode","<int>","processing mode (default 0) 0: CSV loci input file, 1: Mutational probabilities file, 2: FindApproxMatches CSV, 3: Genhyperconserved x2 CSV  4: GenMAlignScore 5: GenMAlignScore m3 7: psl2csv CSV 8: blast2csv CSV");

struct arg_file *infile = arg_file1("i","in","<file>",			"input from csv file");
struct arg_file *FilterRefIDFile = arg_file0("X",NULL,"<file>",	"exclude any RefIDs in this filter file");
struct arg_file *outfile= arg_file1("o","out","<file>",			"output to BED or WIG file");
struct arg_int *MinLenEl = arg_int0("l", "minlenel","<int>",	"mode 0,2,3 minimum length element (default 20)");
struct arg_int *MaxLenEl = arg_int0("L", "maxlenel","<int>",	"mode 0,2,3 max length element (default 100000000)");

struct arg_int *MinInstances = arg_int0("y", "mininsts","<int>",	"mode 7 minimum number of instances with same RefID (default 1)");
struct arg_int *MaxInstances = arg_int0("Y", "maxinsts","<int>",	"mode 7 maximum number of instances with same RefID (default 100000000)");

struct arg_int *MinIdentity = arg_int0("d", "minident","<int>",	"mode 3 minimum out species alignment identity (default 0)");
struct arg_int *MaxIdentity = arg_int0("k", "maxident","<int>",	"mode 3 maximum out species alignment identity (default 100)");

struct arg_str *descr = arg_str0("D","descr","<string>",		"BED description");
struct arg_str *track = arg_str0("t","track","<string>",		"BED track name");
struct arg_str *chrom = arg_str0("c","chromosome","<string>",	"chromosome");
struct arg_str *name = arg_str0("n","name","<string>",			"entry name");
struct arg_int *nbins = arg_int0("N", "numbins","<int>",		"number of bins (default 10)");
struct arg_int *windowsize = arg_int0("w", "windowsize","<int>","windowsize to smooth over (default 5)");
struct arg_int *entry = arg_int0("e", "entry","<int>",			"only output these matching entries");
struct arg_int *entryfld = arg_int0("E", "entryfield","<int>",	"which field (1..n) contains entry value");
struct arg_int *stepfld = arg_int0("z", "stepfield","<int>",	"which field (1..n) contains step value");
struct arg_int *startstep = arg_int0("s", "startstep","<int>",	"start step (default is 1st)");
struct arg_int *endstep = arg_int0("x", "endstep","<int>",		"end step (default is last)");
struct arg_int *remap = arg_int0("r", "remap","<int>",			"remap steps in output by N (default is 0)");
struct arg_int *valuefield = arg_int0("V", "valuefield","<int>","which field (1..n) contains value");
struct arg_lit  *append  = arg_lit0("a","append",               "append to output file");
struct arg_lit  *wiggle  = arg_lit0("W","wiggle",               "create WIG wiggle file (default is BED)");


struct arg_end *end = arg_end(20);
void *argtable[] = {help,version,FileLogLevel,ScreenLogLevel,LogFile,
			Mode,MinLenEl,MaxLenEl,MinIdentity,MaxIdentity,MinInstances,MaxInstances,
			infile,FilterRefIDFile,outfile,descr,track,chrom,name,nbins,entry,entryfld,windowsize,stepfld,startstep,endstep,remap,valuefield,append,wiggle,end};

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

	iMode = Mode->count ? Mode->ival[0] : eCSVFdefault;
	if(iMode == eCSVFSeq || iMode < eCSVFdefault || iMode > eCSVFAlignBLASTcsv)
		{
		printf("Error: Mode specified '-m%d' not in range 0..5, or 7",iMode);
		exit(1);
		}
	
	if(iMode == eCSVFAlignPSLcsv || iMode == eCSVFAlignBLASTcsv)
		{
		iMinInstances = MinInstances->count ? MinInstances->ival[0] : 1;
		iMaxInstances = MaxInstances->count ? MaxInstances->ival[0] : 100000000;
		if(iMinInstances < 1 || iMinInstances > iMaxInstances || iMaxInstances > 100000000)
			{
			printf("Error: Min instances '-y%d' and max instances '-Y%d' not in range 1..100000000",iMinInstances,iMaxInstances);
			exit(1);
			}
		}
	else
		{
		iMinInstances = 0;
		iMaxInstances = 0;
		}

	if(iMode != eCSVFprobe)
		{
		iMinLenEl = MinLenEl->count ? MinLenEl->ival[0] : cDfltMinLenEl;
		if(iMinLenEl < 0 || iMinLenEl > cMaxLenEl)
			{
			printf("\nError: Minimum element length specified as '-l%d' must be in range 0..100000000\n",iMinLenEl);
			exit(1);
			}

		iMaxLenEl = MaxLenEl->count ? MaxLenEl->ival[0] : cDfltMaxLenEl;
		if(iMaxLenEl < iMinLenEl || iMaxLenEl > cMaxLenEl)
			{
			printf("\nError: Maximum element length specified as '-L%d' must be in range %d..100000000\n",iMaxLenEl,iMinLenEl);
			exit(1);
			}

		if(iMode == eCSVFhyperX2)
			{
			iMinIdentity = MinIdentity->count ? MinIdentity->ival[0] : 0;
			if(iMinIdentity < 0 || iMinIdentity > 100)
				{
				printf("\nError: Minimum identity specified as '-d%d' must be in range 0..100\n",iMinIdentity);
				exit(1);
				}

			iMaxIdentity = MaxIdentity->count ? MaxIdentity->ival[0] : 100;
			if(iMaxIdentity < iMinIdentity || iMaxIdentity > 100)
				{
				printf("\nError: Maximum identity specified as '-D%d' must be in range %d..100\n",iMaxIdentity,iMinIdentity);
				exit(1);
				}
			}

		if(iMode == eCSVFdefault)
			{
			if(append->count)
				bAppend = true;
			else
				bAppend = false;

			if(!track->count)
				{
				szTrack[0] = '\0';
				szDescr[0] = '\0';
				}
			else
				{
				strcpy(szTrack,track->sval[0]);
				if(descr->count)
					strcpy(szDescr,descr->sval[0]);
				else
					strcpy(szDescr,szTrack);
				}
			}
		}
	else
		{
		iMaxIdentity = 100;
		iMinIdentity = 0;
		iMaxLenEl = cDfltMaxLenEl;
		iMinLenEl = cDfltMinLenEl;
		if(wiggle->count)
			bWiggle = true;
		else
			bWiggle = false;

		if(append->count)
			bAppend = true;
		else
			bAppend = false;

		if(windowsize->count)
			iWindowsize = windowsize->ival[0];
		else
			iWindowsize = 5;

		if(entry->count)
			iEntry = entry->ival[0];
		else
			iEntry = 0;

		if(remap->count)
			iRemap = remap->ival[0];
		else
			iRemap = 0;

		if(entryfld->count)
			iEntryfld = entryfld->ival[0];
		else
			iEntryfld = 0;

		if(stepfld->count)
			iStepfld = stepfld->ival[0];
		else
			iStepfld = 0;

		if(startstep->count)
			iStartstep = startstep->ival[0];
		else
			iStartstep = 0;

		if(endstep->count)
			iEndstep = endstep->ival[0];
		else
			iEndstep = 0;

		if(valuefield->count)
			iValuefield = valuefield->ival[0];
		else
			iValuefield = 0;

		if(nbins->count)
			iNbins = nbins->ival[0];
		else
			iNbins = 10;


		if(chrom->count)
			strcpy(szChrom,chrom->sval[0]);
		else
			szChrom[0] = '\0';

		if(name->count)
			strcpy(szName,name->sval[0]);
		else
			szName[0] = '\0';

		if(!descr->count)
			{
			printf("\nError: No description specified");
			exit(1);
			}
		strcpy(szDescr,descr->sval[0]);
	
		if(!track->count)
			{
			printf("\nError: No track specified");
			exit(1);
			}
		strcpy(szTrack,track->sval[0]);
		pCsv2Bed = new CCSV2BED;
		}

	if(FilterRefIDFile->count)
		{
		strncpy(szFilterRefIDFile,FilterRefIDFile->filename[0],sizeof(szFilterRefIDFile));
		szFilterRefIDFile[sizeof(szFilterRefIDFile)-1] = '\0';
		}
	else
		szFilterRefIDFile[0] = '\0';

	strcpy(szInfile,infile->filename[0]);
	strcpy(szOutfile,outfile->filename[0]);

		// now that command parameters have been parsed then initialise diagnostics log system
	if(!gDiagnostics.Open(szLogFile,(etDiagLevel)iScreenLogLevel,(etDiagLevel)iFileLogLevel,true))
		{
		printf("\nError: Unable to start diagnostics subsystem.");
		if(szLogFile[0] != '\0')
			printf(" Most likely cause is that logfile '%s' can't be opened/created",szLogFile);
		exit(1);
		}

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Version: %d.%2.2d Processing parameters:",cProgVer/100,cProgVer%100);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"input CSV file: '%s'",szInfile);				// csv file containing structual values
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"BED file to create: '%s'", szOutfile);			// BED file to create
	if(szFilterRefIDFile[0])
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Exclude any RefIDs in this filter file: '%s'",szFilterRefIDFile);

	if(iMode!=eCSVFprobe)
		{
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Minimum element length: %d", iMinLenEl);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Maximum element length: %d", iMaxLenEl);
		if(iMode == eCSVFAlignPSLcsv || iMode == eCSVFAlignBLASTcsv)
			{
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Minimum RefID instances: %d", iMinInstances);
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Maximum RefID instances: %d", iMaxInstances);
			}
		else
			if(iMode == eCSVFhyperX2)
				{
				gDiagnostics.DiagOutMsgOnly(eDLInfo,"Minimum out group identity: %d", iMinIdentity);
				gDiagnostics.DiagOutMsgOnly(eDLInfo,"Maximum out group identity: %d", iMaxIdentity);
				}
		if(iMode == eCSVFdefault)
			{
			if(szTrack[0] != '\0' && szDescr[0] != '\0')
				{
				gDiagnostics.DiagOutMsgOnly(eDLInfo,"Track name: '%s'",szTrack);				// track name for BED file
				gDiagnostics.DiagOutMsgOnly(eDLInfo,"Description: '%s'",szDescr);				// description for BED file
				}
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Append to existing CSV file: %s",bAppend ? "yes" : "no");	// true to append on to existing file otherwise (false) truncate existing file
			}
		}
	else
		{
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Track name: '%s'",szTrack);				// track name for BED file
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Description: '%s'",szDescr);				// description for BED file
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Chromosome: '%s'",szChrom);				// chromosome for BED file
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Entry name: '%s'",szName);				// entry name for BED file
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Numb bins: %d", iNbins);				// number of classes into which scored values are to be bined
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Data smoothing widow size: %d",iWindowsize);			// window size to use for data mean smoothing (1 if no smoothing)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"CSV file entry field: %d",iEntryfld);			// which field in the csv file contains the entry identifier (1..n)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Match entry identifier: %d",iEntry);				// entry identifer to match 
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"CSV file step psn field: %d",iStepfld);			// which field in the csv file contains the step psn
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Start step: %d", iStartstep);			// starting step
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"End step: %d", iEndstep);				// last step
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Remap steps delta: %d", iRemap);				// remap steps by Remap delta
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"CSV file value field: %d",iValuefield);	// which field contains the value of interest
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Append to existing CSV file: %s",bAppend ? "yes" : "no");	// true to append on to existing file otherwise (false) truncate existing file
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Create as Wiggle or BED: %s",bWiggle ? "Wiggle" : "BED");	// true if to output as WIG file, default is as BED
		}
	gStopWatch.Start();
	switch(iMode) {
		case eCSVFdefault:
			Rslt = Process(iMinLenEl,iMaxLenEl,szInfile,szFilterRefIDFile,szTrack,szDescr,bAppend,szOutfile);
			break;

		case eCSVFprobe:
			Rslt = pCsv2Bed->process(szInfile,	// csv file containing structual values
						  szOutfile,			// BED file to create
						  szTrack,				// track name for BED file
						  szDescr,				// description for BED file
						  szChrom,				// chromosome for BED file
						  szName,				// entry name for BED file
						  iNbins,				// number of classes into which scored values are to be bined
						  iWindowsize,			// window size to use for data mean smoothing (1 if no smoothing)
						  iEntryfld,			// which field in the csv file contains the entry identifier (1..n)
						  iEntry,				// entry identifer to match 
						  iStepfld,				// which field in the csv file contains the step psn
						  iStartstep,			// starting step
						  iEndstep,				// last step
						  iRemap,				// remap steps by Remap delta
						  iValuefield,			// which field contains the value of interest
  						  bAppend,				// true to append on to existing file otherwise (false) truncate existing file
						  bWiggle);				// true if to output as WIG file, default is as BED
			break;

		case eCSVFtarget:
			Rslt = ProcessFindApproxMatches(iMinLenEl,iMaxLenEl,szInfile,szFilterRefIDFile,szOutfile);
			break;

		case eCSVFhyperX2:
			Rslt = ProcessGenhyperconservedX2(iMinLenEl,iMaxLenEl,iMinIdentity,iMaxIdentity,szInfile,szFilterRefIDFile,szOutfile);
			break;

		case eCSVFAlignM3:
			Rslt = ProcessGenMAlignScoreM3(iMinLenEl,iMaxLenEl,szInfile,szFilterRefIDFile,szOutfile);
			break;

		case eCSVFAlignPSLcsv:
			Rslt = ProcessPSL2CSV(iMinLenEl,iMaxLenEl,iMinInstances,iMaxInstances,szInfile,szFilterRefIDFile,szOutfile);
			break;

		case eCSVFAlignBLASTcsv:
			Rslt = ProcessBLAST2CSV(iMinLenEl,iMaxLenEl,iMinInstances,iMaxInstances,szInfile,szFilterRefIDFile,szOutfile);
			break;

		default:
			printf("\nError: Unsupported processing mode");
			Rslt = eBSFerrParams;
			break;
		}
	if(pCsv2Bed!=NULL)
		delete pCsv2Bed;
	gStopWatch.Stop();
	Rslt = Rslt < 0 ? 1 : 0;
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Exit Code: %d Total processing time: %s",Rslt,gStopWatch.Read());
	}
else
	{
	arg_print_errors(stdout,end,gszProcName);
	arg_print_syntax(stdout,argtable,"\nend of help\n");
	exit(1);
	}
return(0);
}

int
Process(int MinLen,					// filter out elements of less than this length
		int MaxLen,					// or longer than this length
		char *pszInfile,			// input csv file containing element loci
		char *pszFilterRefIDFile,	// RefIDs to filter
		char *pszTrackName,			// optional track name; if track name and description specified specified then BED header will be created
		char *pszDescr,				// optional track description; if track name and description specified specified then BED header will be created
		bool bAppend,				// false to create/truncate existing file, true to append
		char *pszOutfile)			// BED file to create
{
int NumFields;
int Rslt;
int SrcID;
char *pszChrom;
char *pszAliasChrom;
char *pszElType;
char *pszRefSpecies;
int StartLoci;
int EndLoci;
int Len;
char *pszStrand;
char Strand;
int openopts;
INT64 FileLen;

int LineLen;
int hFile;
char szLineBuff[cLineBuffSize];

int NumElsRead;		// number of elements before filtering
int NumElsAccepted;	// number of elements accepted after filtering
int NumFiltRefIDs;	// number of elements filtered out because of matching RefID
int NumFiltLen;		// number of elements filtered out because of length

CFilterRefIDs *pFilterRefIDs = NULL;

if(pszFilterRefIDFile != NULL && pszFilterRefIDFile[0] != '\0')
	{
	pFilterRefIDs = new CFilterRefIDs;
	if((Rslt=pFilterRefIDs->Open(pszFilterRefIDFile)) < 0)
		{
		while(pFilterRefIDs->NumErrMsgs())
			gDiagnostics.DiagOut(eDLFatal,gszProcName,pFilterRefIDs->GetErrMsg());
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open/parse input filter RefID file %s",pszFilterRefIDFile);
		delete pFilterRefIDs;
		return(Rslt);
		}
	}

CCSVFile *pCSV = new CCSVFile;
if(pCSV == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CCSVfile");
	if(pFilterRefIDs != NULL)
		delete pFilterRefIDs;
	return(eBSFerrObj);
	}

if((Rslt=pCSV->Open(pszInfile))!=eBSFSuccess)
	{
	while(pCSV->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,pCSV->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open file: %s",pszInfile);
	delete pCSV;
	if(pFilterRefIDs != NULL)
		delete pFilterRefIDs;
	return(Rslt);
	}

// create or open existing output file taking into account if the file is to be appended or truncated
#ifdef _WIN32
openopts = _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT;
openopts |= bAppend ? 0 : _O_TRUNC;
hFile = open(pszOutfile,openopts, _S_IREAD | _S_IWRITE );
#else
openopts = O_RDWR | O_CREAT;
openopts |= bAppend ? 0 : O_TRUNC;
hFile = open(pszOutfile,openopts, S_IREAD | S_IWRITE );
#endif
if(hFile == -1)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open %s - %s",pszOutfile,strerror(errno));
	delete pCSV;
	if(pFilterRefIDs != NULL)
		delete pFilterRefIDs;
	return(eBSFerrCreateFile);
	}
if(bAppend)
	FileLen = _lseeki64(hFile,0,SEEK_END);

if((!bAppend || FileLen == 0) && pszTrackName != NULL && pszTrackName[0] != '\0' && pszDescr != NULL && pszDescr[0] != '\0')
	{
	LineLen = sprintf(szLineBuff,"\ntrack name=%s description=\"%s\" useScore=0\n",pszTrackName,pszDescr);
	CUtility::SafeWrite(hFile,szLineBuff,LineLen);
	LineLen = 0;
	}
	

NumElsRead=0;		// number of elements before filtering
NumElsAccepted=0;	// number of elements accepted after filtering
NumFiltRefIDs=0;	// number of elements filtered out because of matching RefID
NumFiltLen=0;		// number of elements filtered out because of length
LineLen = 0;
while((Rslt=pCSV->NextLine()) > 0)	// onto next line containing fields
	{
	NumFields = pCSV->GetCurFields();
	if(NumFields < 7)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expected at least 7 fields in '%s', GetCurFields() returned '%d'",pszInfile,NumFields);
		close(hFile);
		delete pCSV;
		if(pFilterRefIDs != NULL)
			delete pFilterRefIDs;
		return(eBSFerrFieldCnt);
		}

	if(!NumElsRead && pCSV->IsLikelyHeaderLine())
		continue;

	NumElsRead += 1;
	pCSV->GetInt(1,&SrcID);
	if(pFilterRefIDs != NULL && pFilterRefIDs->Locate(SrcID))
		{	
		NumFiltRefIDs += 1;
		continue;
		}

	pCSV->GetText(2,&pszElType);
	pCSV->GetText(3,&pszRefSpecies);
	pCSV->GetText(4,&pszChrom);
// a hack because some datasets still reference chloroplast and mitochondria instead of ChrC and ChrM as expected by the UCSC Browser
	pszAliasChrom = pszChrom;
	if(!stricmp("chloroplast",pszChrom))
		pszAliasChrom = (char *)"ChrC";
	else
		if(!stricmp("mitochondria",pszChrom))
			pszAliasChrom = (char *)"ChrM";
	pCSV->GetInt(5,&StartLoci);
	pCSV->GetInt(6,&EndLoci);
	pCSV->GetInt(7,&Len);
	if(Len < MinLen || Len > MaxLen)
		{
		NumFiltLen += 1;
		continue;
		}

	Strand = '+';				// default to watson if not specified
	if(NumFields > 7)			// see if strand is specified
		{
		pCSV->GetText(8,&pszStrand);
		if(*pszStrand == '+' || *pszStrand == '-')
			Strand = *pszStrand;
		}

	LineLen += sprintf(&szLineBuff[LineLen],"%s\t%d\t%d\t%s_El%d\t%d\t%c\n",
		pszAliasChrom,StartLoci,StartLoci+Len,pszAliasChrom,SrcID,Len < 1000 ? Len : 1000,Strand);
	if(LineLen > (cLineBuffSize-200))
		{
		CUtility::SafeWrite(hFile,szLineBuff,LineLen);
		LineLen = 0;
		}
	NumElsAccepted += 1;
	}
if(LineLen)
	CUtility::SafeWrite(hFile,szLineBuff,LineLen);
close(hFile);
delete pCSV;
if(pFilterRefIDs != NULL)
	delete pFilterRefIDs;

gDiagnostics.DiagOut(eDLInfo,gszProcName,"%s: Read: %d, Accepted: %d, FiltRefIDs: %d, FiltLen: %d",pszInfile,NumElsRead,NumElsAccepted,NumFiltRefIDs,NumFiltLen);
return(Rslt < 0 ? NumElsAccepted : Rslt);
}

int
ProcessFindApproxMatches(int MinLen,// filter out elements of less than this length
		int MaxLen,					// or longer than this length
		char *pszInfile,			// csv file containing element loci
		char *pszFilterRefIDFile,	// RefIDs to filter
		char *pszOutfile)			// BED file to create
{
int NumFields;
int Rslt;
int SrcID;
char *pszChrom;
char *pszStrand;
int StartLoci;
int Len;

int LineLen;
int hFile;
char szLineBuff[cLineBuffSize];

int NumElsRead;		// number of elements before filtering
int NumElsAccepted;	// number of elements accepted after filtering
int NumFiltRefIDs;	// number of elements filtered out because of matching RefID
int NumFiltLen;		// number of elements filtered out because of length


CFilterRefIDs *pFilterRefIDs = NULL;

if(pszFilterRefIDFile != NULL && pszFilterRefIDFile[0] != '\0')
	{
	pFilterRefIDs = new CFilterRefIDs;
	if((Rslt=pFilterRefIDs->Open(pszFilterRefIDFile)) < 0)
		{
		while(pFilterRefIDs->NumErrMsgs())
			gDiagnostics.DiagOut(eDLFatal,gszProcName,pFilterRefIDs->GetErrMsg());
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open/parse input filter RefID file %s",pszFilterRefIDFile);
		delete pFilterRefIDs;
		return(Rslt);
		}
	}

CCSVFile *pCSV = new CCSVFile;
if(pCSV == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CCSVfile");
	if(pFilterRefIDs != NULL)
		delete pFilterRefIDs;
	return(eBSFerrObj);
	}

if((Rslt=pCSV->Open(pszInfile))!=eBSFSuccess)
	{
	while(pCSV->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,pCSV->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open file: %s",pszInfile);
	delete pCSV;
	if(pFilterRefIDs != NULL)
		delete pFilterRefIDs;
	return(Rslt);
	}

#ifdef _WIN32
if((hFile = open(pszOutfile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
if((hFile = open(pszOutfile, O_RDWR | O_CREAT | O_TRUNC, S_IREAD | S_IWRITE ))==-1)
#endif
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create output BED file: %s - %s",pszOutfile,strerror(errno));
	delete pCSV;
	if(pFilterRefIDs != NULL)
		delete pFilterRefIDs;
	return(eBSFerrCreateFile);
	}

NumElsRead=0;		// number of elements before filtering
NumElsAccepted=0;	// number of elements accepted after filtering
NumFiltRefIDs=0;	// number of elements filtered out because of matching RefID
NumFiltLen=0;		// number of elements filtered out because of length
LineLen = 0;

LineLen = 0;
while((Rslt=pCSV->NextLine()) > 0)	// onto next line containing fields
	{
	NumFields = pCSV->GetCurFields();
	if(NumFields < 6)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expected at least 6 fields in '%s', GetCurFields() returned '%d'",pszInfile,NumFields);
		delete pCSV;
		if(pFilterRefIDs != NULL)
			delete pFilterRefIDs;
		return(eBSFerrFieldCnt);
		}

	if(!NumElsRead && pCSV->IsLikelyHeaderLine())
		continue;

	NumElsRead += 1;
	pCSV->GetInt(6,&Len);
	if(Len < MinLen || Len > MaxLen)
		{
		NumFiltLen += 1;
		continue;
		}

	pCSV->GetInt(1,&SrcID);
	if(pFilterRefIDs != NULL && pFilterRefIDs->Locate(SrcID))
		{	
		NumFiltRefIDs += 1;
		continue;
		}

	pCSV->GetText(3,&pszChrom);
	pCSV->GetText(4,&pszStrand);
	pCSV->GetInt(5,&StartLoci);

	LineLen += sprintf(&szLineBuff[LineLen],"%s\t%d\t%d\tEl%d\t%d\t%s\n",
		pszChrom,StartLoci,StartLoci+Len,SrcID,Len < 1000 ? Len : 1000,pszStrand);
	if(LineLen > (cLineBuffSize-200))
		{
		CUtility::SafeWrite(hFile,szLineBuff,LineLen);
		LineLen = 0;
		}
	NumElsAccepted += 1;
	}
if(LineLen)
	CUtility::SafeWrite(hFile,szLineBuff,LineLen);
close(hFile);
delete pCSV;
if(pFilterRefIDs != NULL)
	delete pFilterRefIDs;
gDiagnostics.DiagOut(eDLInfo,gszProcName,"%s: Read: %d, Accepted: %d, FiltRefIDs: %d, FiltLen: %d",pszInfile,NumElsRead,NumElsAccepted,NumFiltRefIDs,NumFiltLen);
return(Rslt < 0 ? NumElsAccepted : Rslt);
}

int
ProcessGenhyperconservedX2(int MinLen,// filter out core elements of less than this length
		int MaxLen,					// or longer than this length
		int MinIdentity,			// filter out if out species has less than this identity (0..100)
		int MaxIdentity,			// filter out if out species has more than this identity (MinIdentity..100)
		char *pszInfile,			// csv file containing element loci
		char *pszFilterRefIDFile,	// RefIDs to filter
		char *pszOutfile)			// BED file to create
{
int NumFields;
int Rslt;
int SrcID;
char *pszChrom;
char *pszElType;
char *pszRefSpecies;
char *pszRelSpecies;
int StartLoci;
int EndLoci;
int Len;
int Features;

int OGUnalignedBases;
int OGMatchCnt;
int OGMismatchCnt;
int OGInDelCnt;

int Identity;

int LineLen;
int hFile;
char szLineBuff[cLineBuffSize];
int NumElsRead;		// number of elements before filtering
int NumElsAccepted;	// number of elements accepted after filtering
int NumFiltRefIDs;	// number of elements filtered out because of matching RefID
int NumFiltLen;		// number of elements filtered out because of length
int NumFiltIdentity;  // number of elements filtered out because of identity


CFilterRefIDs *pFilterRefIDs = NULL;

if(pszFilterRefIDFile != NULL && pszFilterRefIDFile[0] != '\0')
	{
	pFilterRefIDs = new CFilterRefIDs;
	if((Rslt=pFilterRefIDs->Open(pszFilterRefIDFile)) < 0)
		{
		while(pFilterRefIDs->NumErrMsgs())
			gDiagnostics.DiagOut(eDLFatal,gszProcName,pFilterRefIDs->GetErrMsg());
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open/parse input filter RefID file %s",pszFilterRefIDFile);
		delete pFilterRefIDs;
		return(Rslt);
		}
	}

CCSVFile *pCSV = new CCSVFile;
if(pCSV == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CCSVfile");
	if(pFilterRefIDs != NULL)
		delete pFilterRefIDs;
	return(eBSFerrObj);
	}

if((Rslt=pCSV->Open(pszInfile))!=eBSFSuccess)
	{
	while(pCSV->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,pCSV->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open file: %s",pszInfile);
	delete pCSV;
	if(pFilterRefIDs != NULL)
		delete pFilterRefIDs;
	return(Rslt);
	}

#ifdef _WIN32
if((hFile = open(pszOutfile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
if((hFile = open(pszOutfile, O_RDWR | O_CREAT | O_TRUNC, S_IREAD | S_IWRITE ))==-1)
#endif
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create output BED file: %s - %s",pszOutfile,strerror(errno));
	delete pCSV;
	if(pFilterRefIDs != NULL)
		delete pFilterRefIDs;
	return(eBSFerrCreateFile);
	}

NumElsRead=0;		// number of elements before filtering
NumElsAccepted=0;	// number of elements accepted after filtering
NumFiltRefIDs=0;	// number of elements filtered out because of matching RefID
NumFiltLen=0;		// number of elements filtered out because of length
NumFiltIdentity=0;  // number of elements filtered out because of identity

LineLen = 0;
while((Rslt=pCSV->NextLine()) > 0)	// onto next line containing fields
	{
	NumFields = pCSV->GetCurFields();
	if(NumFields < 13)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expected at least 13 fields in '%s', GetCurFields() returned '%d'",pszInfile,NumFields);
		delete pCSV;
		if(pFilterRefIDs != NULL)
			delete pFilterRefIDs;
		return(eBSFerrFieldCnt);
		}

	if(!NumElsRead && pCSV->IsLikelyHeaderLine())
		continue;

	NumElsRead += 1;
	pCSV->GetInt(7,&Len);
	if(Len < MinLen || Len > MaxLen)
		{
		NumFiltLen += 1;
		continue;
		}
	pCSV->GetInt(1,&SrcID);
	if(pFilterRefIDs != NULL && pFilterRefIDs->Locate(SrcID))
		{	
		NumFiltRefIDs += 1;
		continue;
		}


	pCSV->GetInt(11,&OGMatchCnt);
	Identity = (int)(((100.0 * OGMatchCnt) / (double)Len) + 0.001);
	if(Identity < MinIdentity || Identity > MaxIdentity)
		{
		NumFiltIdentity += 1;
		continue;
		}

	pCSV->GetText(2,&pszElType);
	pCSV->GetText(3,&pszRefSpecies);
	pCSV->GetText(4,&pszChrom);
	pCSV->GetInt(5,&StartLoci);
	pCSV->GetInt(6,&EndLoci);

	pCSV->GetText(8,&pszRelSpecies);
	pCSV->GetInt(9,&Features);

	pCSV->GetInt(10,&OGUnalignedBases);
	pCSV->GetInt(12,&OGMismatchCnt);
	pCSV->GetInt(13,&OGInDelCnt);


	LineLen += sprintf(&szLineBuff[LineLen],"%s\t%d\t%d\tEl%d\t%d\t+\n",
		pszChrom,StartLoci,StartLoci+Len,SrcID,Identity*10);
	if(LineLen > (cLineBuffSize-200))
		{
		CUtility::SafeWrite(hFile,szLineBuff,LineLen);
		LineLen = 0;
		}
	NumElsAccepted += 1;
	}
if(LineLen)
	CUtility::SafeWrite(hFile,szLineBuff,LineLen);
close(hFile);
delete pCSV;
if(pFilterRefIDs != NULL)
	delete pFilterRefIDs;
gDiagnostics.DiagOut(eDLInfo,gszProcName,"%s: Read: %d, Accepted: %d, FiltRefIDs: %d, FiltLen: %d, FiltIdentity",pszInfile,NumElsRead,NumElsAccepted,NumFiltRefIDs,NumFiltLen,NumFiltIdentity);
return(Rslt < 0 ? NumElsAccepted : Rslt);
}


int
ProcessGenMAlignScoreM3(int MinLen,	// filter out elements of less than this length
		int MaxLen,					// or longer than this length
		char *pszInfile,			// csv file containing element loci
		char *pszFilterRefIDFile,	// RefIDs to filter
		char *pszOutfile)			// BED file to create
{
int NumFields;
int Rslt;
int SrcID;
char *pszChrom;
char *pszElType;
char *pszRefSpecies;
int StartLoci;
int EndLoci;
int Len;

int LineLen;
int hFile;
char szLineBuff[cLineBuffSize];
int NumElsRead;		// number of elements before filtering
int NumElsAccepted;	// number of elements accepted after filtering
int NumFiltRefIDs;	// number of elements filtered out because of matching RefID
int NumFiltLen;		// number of elements filtered out because of length

CFilterRefIDs *pFilterRefIDs = NULL;

if(pszFilterRefIDFile != NULL && pszFilterRefIDFile[0] != '\0')
	{
	pFilterRefIDs = new CFilterRefIDs;
	if((Rslt=pFilterRefIDs->Open(pszFilterRefIDFile)) < 0)
		{
		while(pFilterRefIDs->NumErrMsgs())
			gDiagnostics.DiagOut(eDLFatal,gszProcName,pFilterRefIDs->GetErrMsg());
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open/parse input filter RefID file %s",pszFilterRefIDFile);
		delete pFilterRefIDs;
		return(Rslt);
		}
	}

CCSVFile *pCSV = new CCSVFile;
if(pCSV == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CCSVfile");
	if(pFilterRefIDs != NULL)
		delete pFilterRefIDs;
	return(eBSFerrObj);
	}

if((Rslt=pCSV->Open(pszInfile))!=eBSFSuccess)
	{
	while(pCSV->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,pCSV->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open file: %s",pszInfile);
	delete pCSV;
	if(pFilterRefIDs != NULL)
		delete pFilterRefIDs;
	return(Rslt);
	}

#ifdef _WIN32
if((hFile = open(pszOutfile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
if((hFile = open(pszOutfile, O_RDWR | O_CREAT | O_TRUNC, S_IREAD | S_IWRITE ))==-1)
#endif
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create output BED file: %s - %s",pszOutfile,strerror(errno));
	delete pCSV;
	if(pFilterRefIDs != NULL)
		delete pFilterRefIDs;
	return(eBSFerrCreateFile);
	}

NumElsRead=0;		// number of elements before filtering
NumElsAccepted=0;	// number of elements accepted after filtering
NumFiltRefIDs=0;	// number of elements filtered out because of matching RefID
NumFiltLen=0;		// number of elements filtered out because of length
LineLen = 0;
while((Rslt=pCSV->NextLine()) > 0)	// onto next line containing fields
	{
	NumFields = pCSV->GetCurFields();
	if(NumFields < 7)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expected at least 7 fields in '%s', GetCurFields() returned '%d'",pszInfile,NumFields);
		delete pCSV;
		if(pFilterRefIDs != NULL)
			delete pFilterRefIDs;
		return(eBSFerrFieldCnt);
		}

	if(!NumElsRead && pCSV->IsLikelyHeaderLine())
		continue;

	NumElsRead += 1;
	pCSV->GetInt(2,&Len);
	if(Len < MinLen || Len > MaxLen)
		{
		NumFiltLen += 1;
		continue;
		}
	pCSV->GetInt(1,&SrcID);
	if(pFilterRefIDs != NULL && pFilterRefIDs->Locate(SrcID))
		{	
		NumFiltRefIDs += 1;
		continue;
		}

	pszElType = (char *)"BlockID";
	pCSV->GetText(3,&pszRefSpecies);
	pCSV->GetText(4,&pszChrom);
	pCSV->GetInt(5,&StartLoci);
	pCSV->GetInt(6,&EndLoci);

	LineLen += sprintf(&szLineBuff[LineLen],"%s\t%d\t%d\tEl%d\t%d\t+\n",
		pszChrom,StartLoci,StartLoci+Len,SrcID,Len < 1000 ? Len : 1000);
	if(LineLen > (cLineBuffSize-200))
		{
		CUtility::SafeWrite(hFile,szLineBuff,LineLen);
		LineLen = 0;
		}
	NumElsAccepted += 1;
	}
if(LineLen)
	CUtility::SafeWrite(hFile,szLineBuff,LineLen);
close(hFile);
delete pCSV;
if(pFilterRefIDs != NULL)
	delete pFilterRefIDs;
if(Rslt >= 0)
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"%s: Read: %d, Accepted: %d, FiltRefIDs: %d, FiltLen: %d",pszInfile,NumElsRead,NumElsAccepted,NumFiltRefIDs,NumFiltLen);
return(Rslt < 0 ? NumElsAccepted : Rslt);
}

const int cRefIDInitalAllocNum = 10000;
const int cRefIDGrowAllocNum = 10000;

typedef struct TAG_sRefIDChrom {
	unsigned short  ChromID;			// uniquely identifies chromosome
	unsigned short Hash;				// hash on chromosome name - GenHash()
	char szChromName[cMaxDatasetSpeciesChrom];	// chromosome name
} tsRefIDChrom;

typedef struct TAG_sRefIDInstances {
	int SrcID;				// instance identifier
	int ChromID;			// on this chromosome
	char Strand;			// on this strand
	int StartLoci;			// between this start
	int EndLoci;			// and this end
	int InstCnt;			// number of instances with same instance identifier
} tsRefIDInstances;


tsRefIDChrom *gpChroms = NULL;
int gNumChromsAllocd = 0;
int gCurNumChroms = 0;
int gMRAChromID = 0;

tsRefIDInstances *gpRefIDInstances;
int gNumRefIDInstancesAllocd=0;
int gCurNumRefIDInstances=0;

// SortRefIDInstances
// Used to sort RefID instances by RefID-->ChromID --->Strand -->StartLoci ascending
static int 
SortRefIDInstances( const void *arg1, const void *arg2)
{
tsRefIDInstances *pEl1 = (tsRefIDInstances *)arg1;
tsRefIDInstances *pEl2 = (tsRefIDInstances *)arg2;
if(pEl1->SrcID < pEl2->SrcID)
	return(-1);
if(pEl1->SrcID > pEl2->SrcID)
	return(1);
if(pEl1->ChromID < pEl2->ChromID)
	return(-1);
if(pEl1->ChromID > pEl2->ChromID)
	return(1);
if(pEl1->Strand < pEl2->Strand)
	return(-1);
if(pEl1->Strand > pEl2->Strand)
	return(1);
if(pEl1->StartLoci < pEl2->StartLoci)
	return(-1);
if(pEl1->StartLoci > pEl2->StartLoci)
	return(1);
return(0);
}


void
ResetRefIDInstances(void)
{
if(gpChroms != NULL)
	{
	delete gpChroms;
	gpChroms = NULL;
	}
gNumChromsAllocd = 0;
gCurNumChroms = 0;
gMRAChromID = 0;
if(gpRefIDInstances != NULL)
	{
	delete gpRefIDInstances;
	gpRefIDInstances = NULL;
	}
gNumRefIDInstancesAllocd = 0;
gCurNumRefIDInstances = 0;
}

int
AddChrom(char *pszChrom)	// get unique chromosome identifier
{
tsRefIDChrom *pNewChrom;
unsigned short aHash = CUtility::GenHash16(pszChrom);

if(gpChroms == NULL)		// NULL if first time...
	{
	gpChroms = new tsRefIDChrom [cChromsInitalAllocNum];
	if(gpChroms == NULL)
		return(eBSFerrMem);
	gNumChromsAllocd = cChromsInitalAllocNum;
	gCurNumChroms = 0;
	gMRAChromID = 0;
	}
if(gCurNumChroms == gNumChromsAllocd)
	{
	pNewChrom = new tsRefIDChrom [gNumChromsAllocd + cChromsGrowAllocNum];
	if(pNewChrom == NULL)
		return(eBSFerrMem);
	memcpy(pNewChrom,gpChroms,sizeof(tsRefIDChrom) * gCurNumChroms);
	delete gpChroms;
	gpChroms = pNewChrom;
	gNumChromsAllocd += cChromsGrowAllocNum;
	}

if(gCurNumChroms != 0)
	{
	// high probability that chromosome will be one which was last accessed
	if(gMRAChromID > 0)
		{
		pNewChrom = &gpChroms[gMRAChromID-1];
		if(aHash == pNewChrom->Hash && !stricmp(pszChrom,pNewChrom->szChromName))
			return(pNewChrom->ChromID);
		}
	// not most recently accessed, need to do a linear search
	pNewChrom = gpChroms;
	for(gMRAChromID = 1; gMRAChromID <= gCurNumChroms; gMRAChromID++,pNewChrom++)
		if(aHash == pNewChrom->Hash && !stricmp(pszChrom,pNewChrom->szChromName))
			return(pNewChrom->ChromID);
	}
// new chromosome entry
pNewChrom = &gpChroms[gCurNumChroms++];
pNewChrom->Hash = aHash;
strcpy(pNewChrom->szChromName,pszChrom);
pNewChrom->ChromID = gCurNumChroms;
gMRAChromID = gCurNumChroms;
return(gCurNumChroms);
}

int AddLoci(int SrcID,char *pszChrom,char *pszStrand,int StartLoci,int EndLoci)
{
int ChromID;
tsRefIDInstances *pNewRefIDInst;
if((ChromID = AddChrom(pszChrom)) < 0)
	return(ChromID);

if(gpRefIDInstances == NULL)		// NULL if first time...
	{
	gpRefIDInstances = new tsRefIDInstances [cRefIDInitalAllocNum];
	if(gpRefIDInstances == NULL)
		return(eBSFerrMem);
	gNumRefIDInstancesAllocd = cRefIDInitalAllocNum;
	gCurNumRefIDInstances = 0;
	}
if(gCurNumRefIDInstances == gNumRefIDInstancesAllocd)
	{
	pNewRefIDInst = new tsRefIDInstances [gNumRefIDInstancesAllocd + cRefIDGrowAllocNum];
	if(pNewRefIDInst == NULL)
		return(eBSFerrMem);
	memcpy(pNewRefIDInst,gpRefIDInstances,sizeof(tsRefIDInstances) * gCurNumRefIDInstances);
	delete gpRefIDInstances;
	gpRefIDInstances = pNewRefIDInst;
	gNumRefIDInstancesAllocd += cRefIDGrowAllocNum;
	}

pNewRefIDInst = &gpRefIDInstances[gCurNumRefIDInstances++];
pNewRefIDInst->SrcID = SrcID;
pNewRefIDInst->ChromID = ChromID;
pNewRefIDInst->StartLoci = StartLoci;
pNewRefIDInst->EndLoci = EndLoci;
pNewRefIDInst->Strand = *pszStrand;
pNewRefIDInst->InstCnt = 1;
return(gCurNumRefIDInstances);
}

int
ProcessPSL2CSV(int MinLen,			// filter out elements of less than this length
		int MaxLen,					// or longer than this length
		int MinInstances,			// must be at least this number of instances
		int MaxInstances,			// must be no more than this number of instances
		char *pszInfile,			// csv file containing element loci
		char *pszFilterRefIDFile,	// RefIDs to filter
		char *pszOutfile)			// BED file to create
{
int NumFields;
int Rslt;
int SrcID;
char *pszSrcID;
char *pszChrom;
char *pszStrand;

int StartLoci;
int EndLoci;
int Len;
bool bHeader;

int LineLen;
int hFile;
char szLineBuff[cLineBuffSize];
int NumElsRead;		// number of elements before filtering
int NumElsAccepted;	// number of elements accepted after filtering
int NumFiltRefIDs;	// number of elements filtered out because of matching RefID
int NumFiltLen;		// number of elements filtered out because of length

CFilterRefIDs *pFilterRefIDs = NULL;
ResetRefIDInstances();

if(pszFilterRefIDFile != NULL && pszFilterRefIDFile[0] != '\0')
	{
	pFilterRefIDs = new CFilterRefIDs;
	if((Rslt=pFilterRefIDs->Open(pszFilterRefIDFile)) < 0)
		{
		while(pFilterRefIDs->NumErrMsgs())
			gDiagnostics.DiagOut(eDLFatal,gszProcName,pFilterRefIDs->GetErrMsg());
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open/parse input filter RefID file %s",pszFilterRefIDFile);
		delete pFilterRefIDs;
		return(Rslt);
		}
	}

CCSVFile *pCSV = new CCSVFile;
if(pCSV == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CCSVfile");
	if(pFilterRefIDs != NULL)
		delete pFilterRefIDs;
	return(eBSFerrObj);
	}

if((Rslt=pCSV->Open(pszInfile))!=eBSFSuccess)
	{
	while(pCSV->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,pCSV->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open file: %s",pszInfile);
	delete pCSV;
	if(pFilterRefIDs != NULL)
		delete pFilterRefIDs;
	return(Rslt);
	}

#ifdef _WIN32
if((hFile = open(pszOutfile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
if((hFile = open(pszOutfile, O_RDWR | O_CREAT | O_TRUNC, S_IREAD | S_IWRITE ))==-1)
#endif
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create output BED file: %s - %s",pszOutfile,strerror(errno));
	delete pCSV;
	if(pFilterRefIDs != NULL)
		delete pFilterRefIDs;
	return(eBSFerrCreateFile);
	}

NumElsRead=0;		// number of elements before filtering
NumElsAccepted=0;	// number of elements accepted after filtering
NumFiltRefIDs=0;	// number of elements filtered out because of matching RefID
NumFiltLen=0;		// number of elements filtered out because of length
LineLen = 0;
bHeader = true;		// file has header line as 1st line

while((Rslt=pCSV->NextLine()) > 0)	// onto next line containing fields
	{
	NumFields = pCSV->GetCurFields();
	if(NumFields < 11)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expected at least 11 fields in '%s', GetCurFields() returned '%d'",pszInfile,NumFields);
		delete pCSV;
		if(pFilterRefIDs != NULL)
			delete pFilterRefIDs;
		ResetRefIDInstances();
		return(eBSFerrFieldCnt);
		}

	if(bHeader || (!NumElsRead && pCSV->IsLikelyHeaderLine()))
		{
		bHeader = false;
		continue;
		}

	NumElsRead += 1;
	pCSV->GetInt(2,&Len);
	if(Len < MinLen || Len > MaxLen)
		{
		NumFiltLen += 1;
		continue;
		}

	pCSV->GetText(1,&pszSrcID);
	sscanf(pszSrcID,"lcl|%d",&SrcID);
	if(pFilterRefIDs != NULL && pFilterRefIDs->Locate(SrcID))
		{	
		NumFiltRefIDs += 1;
		continue;
		}

	pCSV->GetText(6,&pszChrom);
	pCSV->GetText(7,&pszStrand);
	pCSV->GetInt(10,&StartLoci);
	pCSV->GetInt(11,&EndLoci);

	if((Rslt=AddLoci(SrcID,pszChrom,pszStrand,StartLoci,EndLoci))<0)
		{
		ResetRefIDInstances();
		return(Rslt);
		}
	NumElsAccepted += 1;
	}

delete pCSV;
if(pFilterRefIDs != NULL)
	delete pFilterRefIDs;

if(gCurNumRefIDInstances > 1)
	qsort(gpRefIDInstances,gCurNumRefIDInstances,sizeof(tsRefIDInstances),SortRefIDInstances);

tsRefIDInstances *pCurRefIDInstance = gpRefIDInstances;
tsRefIDInstances *pStartRefIDInstance = gpRefIDInstances;
int InstCnt = 0;
int CurSrcID = -1;
int NumFiltInstances = 0;
int NumAccepted = 0;
int Idx;
for(Idx=0; Idx < gCurNumRefIDInstances; Idx++, pCurRefIDInstance++)
	{
	if(!InstCnt || pCurRefIDInstance->SrcID == CurSrcID)
		{
		CurSrcID = pCurRefIDInstance->SrcID;
		InstCnt++;
		}
	else
		{
		if(InstCnt >= MinInstances && InstCnt <= MaxInstances)
			{
			if(InstCnt > 999)
				InstCnt = 999;
			while(pStartRefIDInstance != pCurRefIDInstance)
				{
				LineLen += sprintf(&szLineBuff[LineLen],"%s\t%d\t%d\tEl%d\t%d\t%c\n",
						gpChroms[pStartRefIDInstance->ChromID-1].szChromName,pStartRefIDInstance->StartLoci,pStartRefIDInstance->EndLoci,pStartRefIDInstance->SrcID,InstCnt,pStartRefIDInstance->Strand);
				if(LineLen > (cLineBuffSize-200))
					{
					CUtility::SafeWrite(hFile,szLineBuff,LineLen);
					LineLen = 0;
					}
				pStartRefIDInstance+=1;
				NumAccepted += 1;
				}
			}
		else
			NumFiltInstances += InstCnt;
		pStartRefIDInstance = pCurRefIDInstance;
		CurSrcID = pCurRefIDInstance->SrcID;
		InstCnt = 1;
		}
	}
if(InstCnt >= MinInstances && InstCnt <= MaxInstances)
	{
	while(pStartRefIDInstance != pCurRefIDInstance)
		{
		LineLen += sprintf(&szLineBuff[LineLen],"%s\t%d\t%d\tEl%d\t%d\t%c\n",
					gpChroms[pStartRefIDInstance->ChromID-1].szChromName,pStartRefIDInstance->StartLoci,pStartRefIDInstance->EndLoci,SrcID,InstCnt,pStartRefIDInstance->Strand);
		if(LineLen > (cLineBuffSize-200))
			{
			CUtility::SafeWrite(hFile,szLineBuff,LineLen);
			LineLen = 0;
			}
		NumAccepted += 1;
		pStartRefIDInstance+=1;
		}
	}
else
	NumFiltInstances += InstCnt;
if(LineLen)
	CUtility::SafeWrite(hFile,szLineBuff,LineLen);
close(hFile);
if(Rslt >= 0)
gDiagnostics.DiagOut(eDLInfo,gszProcName,"%s: Read: %d, Accepted: %d, FiltRefIDs: %d, FiltLen: %d, FiltInstances: %d",
				pszInfile,NumElsRead,NumAccepted,NumFiltRefIDs,NumFiltLen,NumFiltInstances);
return(Rslt < 0 ? NumElsAccepted : Rslt);
}

int
ProcessBLAST2CSV(int MinLen,			// filter out elements of less than this length
		int MaxLen,					// or longer than this length
		int MinInstances,			// must be at least this number of instances
		int MaxInstances,			// must be no more than this number of instances
		char *pszInfile,			// csv file containing element loci
		char *pszFilterRefIDFile,	// RefIDs to filter
		char *pszOutfile)			// BED file to create
{
int NumFields;
int Rslt;
int SrcID;
char *pszSrcID;
char *pszChrom;
char *pszStrand;

int StartLoci;
int EndLoci;
int Len;
bool bHeader;

int LineLen;
int hFile;
char szLineBuff[cLineBuffSize];
int NumElsRead;		// number of elements before filtering
int NumElsAccepted;	// number of elements accepted after filtering
int NumFiltRefIDs;	// number of elements filtered out because of matching RefID
int NumFiltLen;		// number of elements filtered out because of length

CFilterRefIDs *pFilterRefIDs = NULL;
ResetRefIDInstances();

if(pszFilterRefIDFile != NULL && pszFilterRefIDFile[0] != '\0')
	{
	pFilterRefIDs = new CFilterRefIDs;
	if((Rslt=pFilterRefIDs->Open(pszFilterRefIDFile)) < 0)
		{
		while(pFilterRefIDs->NumErrMsgs())
			gDiagnostics.DiagOut(eDLFatal,gszProcName,pFilterRefIDs->GetErrMsg());
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open/parse input filter RefID file %s",pszFilterRefIDFile);
		delete pFilterRefIDs;
		return(Rslt);
		}
	}

CCSVFile *pCSV = new CCSVFile;
if(pCSV == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CCSVfile");
	if(pFilterRefIDs != NULL)
		delete pFilterRefIDs;
	return(eBSFerrObj);
	}

if((Rslt=pCSV->Open(pszInfile))!=eBSFSuccess)
	{
	while(pCSV->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,pCSV->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open file: %s",pszInfile);
	delete pCSV;
	if(pFilterRefIDs != NULL)
		delete pFilterRefIDs;
	return(Rslt);
	}

#ifdef _WIN32
if((hFile = open(pszOutfile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
if((hFile = open(pszOutfile, O_RDWR | O_CREAT | O_TRUNC, S_IREAD | S_IWRITE ))==-1)
#endif
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create output BED file: %s - %s",pszOutfile,strerror(errno));
	delete pCSV;
	if(pFilterRefIDs != NULL)
		delete pFilterRefIDs;
	return(eBSFerrCreateFile);
	}

NumElsRead=0;		// number of elements before filtering
NumElsAccepted=0;	// number of elements accepted after filtering
NumFiltRefIDs=0;	// number of elements filtered out because of matching RefID
NumFiltLen=0;		// number of elements filtered out because of length
LineLen = 0;
bHeader = true;		// file has header line as 1st line

while((Rslt=pCSV->NextLine()) > 0)	// onto next line containing fields
	{
	NumFields = pCSV->GetCurFields();
	if(NumFields < 12)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expected at least 12 fields in '%s', GetCurFields() returned '%d'",pszInfile,NumFields);
		delete pCSV;
		if(pFilterRefIDs != NULL)
			delete pFilterRefIDs;
		ResetRefIDInstances();
		return(eBSFerrFieldCnt);
		}

	if(bHeader || (!NumElsRead && pCSV->IsLikelyHeaderLine()))
		{
		bHeader = false;
		continue;
		}

	NumElsRead += 1;

	pCSV->GetInt(5,&Len);
	if(Len < MinLen || Len > MaxLen)
		{
		NumFiltLen += 1;
		continue;
		}

	pCSV->GetText(1,&pszSrcID);
	sscanf(pszSrcID,"lcl|%d",&SrcID);
	if(pFilterRefIDs != NULL && pFilterRefIDs->Locate(SrcID))
		{	
		NumFiltRefIDs += 1;
		continue;
		}

	pCSV->GetText(2,&pszChrom);
	pCSV->GetText(3,&pszStrand);
	pCSV->GetInt(10,&StartLoci);
	pCSV->GetInt(11,&EndLoci);

	if((Rslt=AddLoci(SrcID,pszChrom,pszStrand,StartLoci,EndLoci))<0)
		{
		ResetRefIDInstances();
		return(Rslt);
		}
	NumElsAccepted += 1;
	}

delete pCSV;
if(pFilterRefIDs != NULL)
	delete pFilterRefIDs;

if(gCurNumRefIDInstances > 1)
	qsort(gpRefIDInstances,gCurNumRefIDInstances,sizeof(tsRefIDInstances),SortRefIDInstances);

tsRefIDInstances *pCurRefIDInstance = gpRefIDInstances;
tsRefIDInstances *pStartRefIDInstance = gpRefIDInstances;
int InstCnt = 0;
int CurSrcID = -1;
int NumFiltInstances = 0;
int NumAccepted = 0;
int Idx;
for(Idx=0; Idx < gCurNumRefIDInstances; Idx++, pCurRefIDInstance++)
	{
	if(!InstCnt || pCurRefIDInstance->SrcID == CurSrcID)
		{
		CurSrcID = pCurRefIDInstance->SrcID;
		InstCnt++;
		}
	else
		{
		if(InstCnt >= MinInstances && InstCnt <= MaxInstances)
			{
			if(InstCnt > 999)
				InstCnt = 999;
			while(pStartRefIDInstance != pCurRefIDInstance)
				{
				LineLen += sprintf(&szLineBuff[LineLen],"%s\t%d\t%d\tEl%d\t%d\t%c\n",
						gpChroms[pStartRefIDInstance->ChromID-1].szChromName,pStartRefIDInstance->StartLoci,pStartRefIDInstance->EndLoci,pStartRefIDInstance->SrcID,InstCnt,pStartRefIDInstance->Strand);
				if(LineLen > (cLineBuffSize-200))
					{
					CUtility::SafeWrite(hFile,szLineBuff,LineLen);
					LineLen = 0;
					}
				pStartRefIDInstance+=1;
				NumAccepted += 1;
				}
			}
		else
			NumFiltInstances += InstCnt;
		pStartRefIDInstance = pCurRefIDInstance;
		CurSrcID = pCurRefIDInstance->SrcID;
		InstCnt = 1;
		}
	}
if(InstCnt >= MinInstances && InstCnt <= MaxInstances)
	{
	while(pStartRefIDInstance != pCurRefIDInstance)
		{
		LineLen += sprintf(&szLineBuff[LineLen],"%s\t%d\t%d\tEl%d\t%d\t%c\n",
					gpChroms[pStartRefIDInstance->ChromID-1].szChromName,pStartRefIDInstance->StartLoci,pStartRefIDInstance->EndLoci,SrcID,InstCnt,pStartRefIDInstance->Strand);
		if(LineLen > (cLineBuffSize-200))
			{
			CUtility::SafeWrite(hFile,szLineBuff,LineLen);
			LineLen = 0;
			}
		NumAccepted += 1;
		pStartRefIDInstance+=1;
		}
	}
else
	NumFiltInstances += InstCnt;
if(LineLen)
	CUtility::SafeWrite(hFile,szLineBuff,LineLen);
close(hFile);
if(Rslt >= 0)
gDiagnostics.DiagOut(eDLInfo,gszProcName,"%s: Read: %d, Accepted: %d, FiltRefIDs: %d, FiltLen: %d, FiltInstances: %d",
				pszInfile,NumElsRead,NumAccepted,NumFiltRefIDs,NumFiltLen,NumFiltInstances);
return(Rslt < 0 ? NumElsAccepted : Rslt);
}





