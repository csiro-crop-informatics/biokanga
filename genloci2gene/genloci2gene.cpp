// genloci2gene.cpp : Defines the entry point for the console application.
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

const char *cpszProgVer = "1.0.0";		// increment with each release


const int cMinWeighting = 0;	  // can have 0 as weighting if that region not to be analysed
const int cMaxWeighting = 100;	  // max allowed weighting
const int cDfltWIntergenic = 1;	  // default intergenic weighting
const int cDfltWUpstream = 4;	  // default upstream regulatory region weighting
const int cDfltWIntragenic = 5;   // default intragenic weighting
const int cDfltWDnstream = 3;	  // default downstream regulatory region weighting

const int cMinDensityWindow  = 10000;	// minimum density window size
const int cDfltDensityWindow = 1000000; // default density window size
const int cMaxChromLen = 512000000;		// can handle upto this sized maximum length chromosome
const int cMaxNumDensityCnts = 1 + (cMaxChromLen / cMinDensityWindow);	// worst case max number of windows in any chromosome we can handle
const int cMaxGOTerms = 50;		  // maximum number of GO terms of interest which can be handled

// processing modes
typedef enum eProcMode {
	ePMStandard = 0,			  // standard processing
	ePMDistribution,			  // distribution characteristics
	ePMDensity					  // distribution density
} etProcMode;

// format of file containing loci information
typedef enum eLociFileFormat {
	eLFormat9,				// 9 field CSV file, field 4: chromosome, field 5: start, field 6: end
	eLFormat3,				// 3 field CSV file, field 1: chromosome, field 2: start, field 3: end
	eLFUCSCBed				// UCSC BED file
} etLociFileFormat;

// generated CSV gene file will contain three fields, field 1: gene name, field 2: weighting, field3: length 
typedef struct TAG_sLociGene { 
	int Length;				// length
	int Weighting;			// weighting associated with this gene
	char szGeneName[50];	// gene name
} tsLociGene;

typedef struct TAG_sDensityCnts {
	int NumGenes;			// number of genes which start in this window
	int NumGeneBases;		// number of bases within genes that start in this window
	int NumFeats;			// number of features e.g loci which start in this window
	int NumFeatBases;		// number of bases within features which start in this window
	int GOCnts[cMaxGOTerms];
	} tsDensityCnts;


typedef struct TAG_sGOTermMap {
	char szTerm[cMaxOBOGOname];		// GO term
	tsGOTerm *pTerm;				// term details
} tsGOTermMap;


const int cPopGenes2Alloc = 10000;	// allocate in this many increments
CStopWatch gStopWatch;
CDiagnostics gDiagnostics;				// for writing diagnostics messages to log file
char gszProcName[_MAX_FNAME];			// process name

char *Mode2Text(etProcMode Mode);
int TrimQuotes(char *pszTxt);
char *LociFormat2Txt(eLociFileFormat LociFormat);

char *Ontologies2Txt(etOntologies Ontology);

int ParseGOTerms(int MaxGOTerms,char *pszTrackGOTerms, tsGOTermMap *pGOTerms);

int
Add2GO(char *pszGene,
	etOntologies Ontology, // which root ontology to process
	bool bProp,				// propagate counts from GO:Term into parent terms
	CBEDfile *pBED,			// gene BED file
	CGOAssocs *pAssocs,		// gene to ID association file
	CGOTerms *pGOTerms);	// GO ontology file

int
Process(etProcMode iMode,		// 
		int Intergenic,			// weighting to use if loci is intergenic
		int Upstream,			// weighting to use if loci is in upstream regulatory region
		int Intragenic,			// weighting to use if loci is intragenic
		int Dnstream,			// weighting to use if loci is in downstream regulatory region
		int RegLen,				// up/down stream regulatory region length 
		int AssocDist,			// maximum association distance
		int ClustDist,			// treat as clustered if loci is less than this distance from end to start of next loci
		char *pszBED,			// gene BED file
		char Strand,			// loci are on this strand
		etLociFileFormat LociFormat, // expected LociFile format
		char *pszLociFile,		 // file containing genomic loci
		char *pszGeneFile,		 // gene + weighting + length file to generate
		etOntologies Ontology, // which ontology
		char *pszGoAssoc,		// gene to ID association file
		char *pszGOTerms,		// GO ontology file
		char *pszTrackGOTerms);		// GO terms to track

int GenDensity(int WindowSize,CBEDfile *pBED,CFeatLoci *pFeatLoci,FILE * pGeneFile,	
			   etOntologies Ontology, // which ontology
		char *pszGoAssoc,		// gene to ID association file
		char *pszGOTerms,		// GO ontology file
		char *pszTrackGOTerms); // GO terms to track

int								// eBSFSuccess, eBSFerrFeature or eBSFerrChrom
Loci2Gene(etProcMode iMode,		// processing mode 
		int Intergenic,			// weighting to use if loci is intergenic
		int Upstream,			// weighting to use if loci is in upstream regulatory region
		int Intragenic,			// weighting to use if loci is intragenic
		int Dnstream,			// weighting to use if loci is in downstream regulatory region
		int RegLen,				// up/down stream regulatory region length 
		char Strand,			// loci are on this strand
		int AssocDist,			// maximum association distance
		CBEDfile *pBED,			// BED file containing genes to map loci onto
		char *pszChrom,			// chromosome
		int LociStart,			// loci start
		int LociEnd,			// loci end
		char *pszGene,			// where to return gene loci maps onto
		int *pWeighting);		// where to return weighting

int								// eBSFSuccess, eBSFerrFeature or eBSFerrChrom
Loci2GeneA(etProcMode Mode,		// processing mode 
		char Strand,			// loci are on this strand
		CBEDfile *pBED,			// BED file containing genes to map loci onto
		char *pszChrom,			// chromosome
		int LociStart,			// loci start
		int LociEnd,			// loci end
		char *pszGene,			// where to return gene loci maps onto
		int *pGeneLen,			// where to return associated gene length
		int *pDistance);			// where to return distance away from gene if loci not intragenic

int SetBits(int Start,int End,unsigned int *pBitVector);
int CountBitsSet(int Start,int End,unsigned int *pBitVector);



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
	return _T("genloci2gene");
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
int iMode;						// processing mode 0:standard, 1:distribution, 3:density
int	iIntergenic;				// weighting to use if loci is intergenic
int	iUpstream;					// weighting to use if loci is in upstream regulatory region
int	iIntragenic;				// weighting to use if loci is intragenic
int	iDnstream;					// weighting to use if loci is in downstream regulatory region
int iRegLen;					// up/down stream regulatory length
int iAssocDist;					// Maximum intergenic gene association distance
int iClustDist;					// Cluster into single if less than this distance between end and start of next loci
int iLociStrand;				// process on this strand
char cStrand;
etOntologies iOntologies;		// ontoligies to process: 1 Cellular, 2 Biological, 3 Molecular

int iLociFormat;				// szLociFile format
char szBED[_MAX_PATH];			// gene BED file
char szLociFile[_MAX_PATH];		// file containing hits
char szGeneFile[_MAX_PATH];		// results file to generate
char szGoAssoc[_MAX_PATH];		// gene to ID association file
char szGOTerms[_MAX_PATH];		// GO ontology file
char szTrackGOTerms[2048];		// all GO terms of interest for which counts are to be returned

// command line args
struct arg_lit  *help    = arg_lit0("hH","help",                "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_int *ScreenLogLevel=arg_int0("S", "ScreenLogLevel",	"<int>","Level of diagnostics written to logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"diagnostics log file");

struct arg_int  *Mode = arg_int0("m","procmode","<int>",			"processing mode 0:standard, 1:distribution, 3:density");
struct arg_int  *RegLen = arg_int0("L","updnstream","<int>",		"length of 5'up or 3'down  stream regulatory region length (default = 2000) 0..1000000");
struct arg_int  *AssocDist = arg_int0("a","assocdist","<int>",		"Maximum intergenic gene association distance");

struct arg_int  *Intergenic = arg_int0("w","intergenic","<int>",	"intergenic weighting (default 1)");
struct arg_int  *Upstream = arg_int0("x","upstream","<int>",	    "regulatory region upstream weighting (default 4)");
struct arg_int  *Intragenic = arg_int0("y","intragenic","<int>",	"intragenic region weighting (default 5)");
struct arg_int  *Dnstream = arg_int0("z","downstream","<int>",		"regulatory region downstream weighting (default 3)");
struct arg_int  *ClustDist = arg_int0("c","clustdist","<int>",		"Cluster distance between end and start of next loci (default 0)");

struct arg_int  *LociStrand = arg_int0("s","strand","<int>",		"loci strand default 0==both,1=='+',2=='-'");
struct arg_file *BED = arg_file1("b","locibed","<file>",			"input loci to gene BED file");
struct arg_int  *LociFormat = arg_int0("l","lociformat","<int>",	"input genome loci file format: 0: 9field loci, 1:3field loci, 2: BED");
struct arg_file *LociFile = arg_file1("i","loci","<file>",			"input genome loci+size .csv or BED file");
struct arg_file *GeneFile = arg_file1("o","output","<file>",		"output results containing genes + weighting file");

struct arg_int  *Ontology = arg_int0("r","ontologies","<int>",		"Ontology: (default)1:Cellular 2:Biological 3:Molecular");
struct arg_file *GoAssoc = arg_file0("g",NULL,"<file>",				"input GO associations file");
struct arg_file *GOTerms = arg_file0("G",NULL,"<file>",				"input GO terms file");
struct arg_str *TrackGOTerms = arg_str0("t",NULL,"<str>",			"space delimited GO terms of interest");

struct arg_end *end = arg_end(20);

void *argtable[] = {help,version,FileLogLevel,ScreenLogLevel,LogFile,
					Mode,RegLen,AssocDist,Intergenic,Upstream,Intragenic,Dnstream,ClustDist,
					LociStrand,BED,LociFormat,LociFile,GeneFile,Ontology,GoAssoc,GOTerms,TrackGOTerms,
					end};

char **pAllArgs;
int argerrors;
argerrors = CUtility::arg_parsefromfile(argc,(char **)argv,&pAllArgs);
if(argerrors >= 0)
	argerrors = arg_parse(argerrors,pAllArgs,argtable);

/* special case: '--help' takes precedence over error reporting */
if (help->count > 0)
        {
		printf("\n%s Map loci to genes, Version %s\nOptions ---\n", gszProcName,cpszProgVer);
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

	iMode = Mode->count ? Mode->ival[0] : ePMStandard;
	if(iMode < ePMStandard || iMode > ePMDensity)
		{
		printf("\nError: Requested processing mode '-m%d' not supported",iMode);
		exit(1);
		}

	iLociFormat = LociFormat->count ? LociFormat->ival[0] : eLFormat9;
	if(iLociFormat < eLFormat9 || iLociFormat > eLFUCSCBed)
		{
		printf("\nError: Requested input file loci format '-F%d' not supported",iLociFormat);
		exit(1);
		}

	iLociStrand = LociStrand->count ? LociStrand->ival[0] : 0;
	if(iLociStrand < 0 || iLociStrand > 2)
		{
		printf("\nError: Requested strand '-s%d' not supported",iLociStrand);
		exit(1);
		}
	switch(iLociStrand) {
		case 0:
			cStrand = '*';
			break;
		case 1:
			cStrand = '+';
			break;
		case 2:
			cStrand = '-';
			break;
		}

	if(iMode == ePMStandard)
		{
		iIntergenic = Intergenic->count ? Intergenic->ival[0] : cDfltWIntergenic;
		if(iIntergenic < cMinWeighting || iIntergenic > cMaxWeighting)
			{
			printf("\nError: Intergenic weighting '-w%d' not supported",iIntergenic);
			exit(1);
			}

		iUpstream = Upstream->count ? Upstream->ival[0] : cDfltWUpstream;
		if(iUpstream < cMinWeighting || iUpstream > cMaxWeighting)
			{
			printf("\nError: Upstream weighting '-x%d' not supported",iUpstream);
			exit(1);
			}

		iIntragenic = Intragenic->count ? Intragenic->ival[0] : cDfltWIntragenic;
		if(iIntragenic < cMinWeighting || iIntragenic > cMaxWeighting)
			{
			printf("\nError: Intragenic weighting '-y%d' not supported",iIntragenic);
			exit(1);
			}

		iDnstream = Dnstream->count ? Dnstream->ival[0] : cDfltWDnstream;
		if(iDnstream < cMinWeighting || iDnstream > cMaxWeighting)
			{
			printf("\nError: Dnstream weighting '-z%d' not supported",iDnstream);
			exit(1);
			}

		iRegLen = RegLen->count ? RegLen->ival[0] : cDfltRegLen;
		if(iRegLen < cMinRegLen || iRegLen > cMaxRegLen)
			{
			printf("\nError: RegLen '-L%d' not supported",iRegLen);
			exit(1);
			}
		if(iRegLen == cMinRegLen)
			{
			iDnstream = 0;
			iUpstream = 0;
			}

		iAssocDist = AssocDist->count ? AssocDist->ival[0] : 0;
		if(iAssocDist < 0)
			{
			printf("\nError: maximum association distance '-a%d' not supported",iAssocDist);
			exit(1);
			}
		if(iRegLen > iAssocDist)
			{
			printf("\nWarning: maximum association distance '-a%d is less than regulatory region length %d, resetting association distance to 0",iAssocDist,iRegLen);
			iAssocDist = 0;
			}

		if(iAssocDist == 0)
			iIntergenic = 0;
		iOntologies = eONTUndefined;
		}
	else
		{
	
		if(iMode == ePMDensity)
			{
			
			iAssocDist = AssocDist->count ? AssocDist->ival[0] : cDfltDensityWindow;
			if(iAssocDist < cMinDensityWindow)
				{
				printf("\nError: density window size less than %d not supported",cMinDensityWindow);
				exit(1);
				}
			
			if(Ontology->count)
				{
				iOntologies = (etOntologies)(Ontology->count ? Ontology->ival[0] : eONTCellular);
				if(iOntologies < eONTCellular || iOntologies > eONTMolecular)
					{
					printf("\nError: Requested ontology '-O%d' not supported",iOntologies);
					exit(1);
					}

				if(!GoAssoc->count)
					{
					printf("\nError: Ontology specified with '-r%d' but no term association file specified with '-g<file>",Ontology->ival[0]);
					exit(1);
					}
				
				if(!GOTerms->count)
					{						
					printf("\nError: Ontology specified with '-r%d' but no GOterm file specified with '-G<file>",Ontology->ival[0]);
					exit(1);
					}

				if(!TrackGOTerms->count)
					{						
					printf("\nError: Ontology specified with '-r%d' but no GO terms of interest specified with '-t<term list>'",Ontology->ival[0]);
					exit(1);
					}

				strcpy(szTrackGOTerms,TrackGOTerms->sval[0]);
				int NumTerms;
				if((NumTerms = ParseGOTerms(cMaxGOTerms + 1,szTrackGOTerms,NULL)) > cMaxGOTerms)
					{
					printf("\nError: Too many GO terms of interest specified with '-t%s' maximum of %d allowed",szTrackGOTerms,cMaxGOTerms);
					exit(1);
					}
				if(NumTerms < 1)
					{
					printf("\nError: At least 1 GO term must be specified - '-t%s'",szTrackGOTerms);
					exit(1);
					}
				strcpy(szGoAssoc,GoAssoc->filename[0]);
				strcpy(szGOTerms,GOTerms->filename[0]);
				}
			else
				{
				if(GoAssoc->count || GOTerms->count)
					{
					printf("\nError: Ontology not specified with '-r<ontology>'");
					exit(1);
					}
				szGoAssoc[0] = '\0';
				szGOTerms[0] = '\0';
				szTrackGOTerms[0] = '\0';
				iOntologies = eONTUndefined;
				}
			}
		else
			{
			iAssocDist = 0;
			iOntologies = eONTUndefined;
			szGoAssoc[0] = '\0';
			szGOTerms[0] = '\0';
			}
		iIntergenic = 0;
		iRegLen = 0;
		iDnstream = 0;
		}

	iClustDist = ClustDist->count ? ClustDist->ival[0] : 0;
	if(iClustDist < 0)
		{
		printf("\nError: Cluster distance '-c%d' less than 0",iClustDist);
		exit(1);
		}
		
	strcpy(szBED,BED->filename[0]);
	strcpy(szLociFile,LociFile->filename[0]);
	strcpy(szGeneFile,GeneFile->filename[0]);

			// now that command parameters have been parsed then initialise diagnostics log system
	if(!gDiagnostics.Open(szLogFile,(etDiagLevel)iScreenLogLevel,(etDiagLevel)iFileLogLevel,true))
		{
		printf("\nError: Unable to start diagnostics subsystem.");
		if(szLogFile[0] != '\0')
			printf(" Most likely cause is that logfile '%s' can't be opened/created",szLogFile);
		exit(1);
		}

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Version: %s",cpszProgVer);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Processing mode: %s",Mode2Text((etProcMode)iMode));
	
	switch(iMode) {
		case ePMStandard:
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Intergenic weighting: %d",iIntergenic);
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Upstream regulatory region weighting: %d",iUpstream);
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Intragenic weighting: %d",iIntragenic);
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Downstream regulatory weighting: %d",iDnstream);
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Regulatory region length: %d",iRegLen);
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Maximum intergenic gene association distance: %d",iAssocDist);
			break;
		case ePMDensity:
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Density window size: %d",iAssocDist);
			if(iOntologies != eONTUndefined)
				{
				gDiagnostics.DiagOutMsgOnly(eDLInfo,"Ontology: %s",Ontologies2Txt(iOntologies));
				gDiagnostics.DiagOutMsgOnly(eDLInfo,"gene to ID association file: '%s'",szGoAssoc);
				gDiagnostics.DiagOutMsgOnly(eDLInfo,"GO ontology file: '%s'",szGOTerms);
				gDiagnostics.DiagOutMsgOnly(eDLInfo,"GO terms to track: '%s'",szTrackGOTerms);
				}
			break;
		case ePMDistribution:
			break;
		}

	if(iClustDist > 0)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Cluster into single if less than this distance between end and start of next loci: %d",iClustDist);
	else
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"No clustering of loci into single");
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Process loci on strand: %c",cStrand);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"input loci to gene mapping BED file: '%s'",szBED);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"expected CSV file loci format: '%s'",LociFormat2Txt((eLociFileFormat)iLociFormat));
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"input CSV file containing loci: '%s'",szLociFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"gene+weighting CSV file to generate: '%s'",szGeneFile);

	gStopWatch.Start();
#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	Rslt = Process((etProcMode)iMode,iIntergenic,iUpstream,iIntragenic,iDnstream,iRegLen,iAssocDist,iClustDist,szBED,cStrand,(etLociFileFormat)iLociFormat,szLociFile,szGeneFile,
									iOntologies,szGoAssoc,szGOTerms,szTrackGOTerms);
	gStopWatch.Stop();
	Rslt = Rslt < 0 ? 1 : 0;
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Exit Code: %d Total processing time: %s",Rslt,gStopWatch.Read());
	exit(Rslt);
	}
else
	{
	arg_print_errors(stdout,end,gszProcName);
	arg_print_syntax(stdout,argtable,"\nend of help\n");
	exit(1);
	}
}

char *
Mode2Text(etProcMode Mode)
{
switch(Mode) {
	case ePMStandard:		  // standard processing
		return((char *)"Standard loci mapping to gene/feature");
	case ePMDistribution:	  // distribution characteristics
		return((char *)"Loci location distribution characteristics");
	case ePMDensity:		  // distribution density
		return((char *)"Loci gene/feature mapping density");
	default:
		break;
	}
return((char *)"Unsupported processing mode");
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

char *
LociFormat2Txt(eLociFileFormat LociFormat)
{
switch(LociFormat) {
	case eLFormat9: return((char *)"9 field CSV file, field 4: chromosome, field 5: start, field 6: end");
	case eLFormat3: return((char *)"3 field CSV file, field 1: chromosome, field 2: start, field 3: end");
	case eLFUCSCBed: return((char *)"UCSC BED file");
	default: return((char *)"Undefined");
	}
}

char *
TrimWhitespace(char *pTxt)
{
char *pStart;
char Chr;
	// strip leading whitespace
while(Chr = *pTxt++)
	if(!isspace(Chr))
			break;
if(Chr == '\0')					// empty line?
	return(pTxt-1);
pStart = pTxt-1;
while(Chr = *pTxt)			// fast forward to line terminator
	pTxt++;
pTxt-=1;
while(Chr = *pTxt--)
	if(!isspace(Chr))
		break;
pTxt[2] = '\0';
return(pStart);
}

char *
Ontologies2Txt(etOntologies Ontology)
{
char *pszOntology;
switch(Ontology) {
	case eONTCellular:
		pszOntology = (char *)"Cellular component";
		break;
	case eONTBiological:
		pszOntology = (char *)"Biological process";
		break;
	case eONTMolecular:
		pszOntology = (char *)"Molecular function";
		break;
	default:
		pszOntology = (char *)"Unrecognised";
		break;
	}
return(pszOntology);
}



int
ParseGOTerms(int MaxGOTerms,char *pszTrackGOTerms, tsGOTermMap *pGOTerms)
{
// parse out species list
char Chr;
char *pGOTerm;
int NumGOTerms = 0;
bool InToken = false;
if(MaxGOTerms < 1 || pszTrackGOTerms == NULL || *pszTrackGOTerms == '\0')
	return(0);

while(Chr = *pszTrackGOTerms++)
	{
	if(Chr == '"' || Chr == '\'') // change any single or double quotes into spaces
		Chr = ' ';
	if(isspace(Chr) || Chr==',')
		{
		if(!InToken)			// slough whitespace or ',' if not inside a token parse
			continue;
		InToken = false;
		pszTrackGOTerms[-1] = '\0';
		if(pGOTerms != NULL)
			{
			strncpy(pGOTerms[NumGOTerms].szTerm,pGOTerm,cMaxOBOGOname);
			pGOTerms[NumGOTerms].szTerm[cMaxOBOGOname-1] = '\0';
			pGOTerms[NumGOTerms].pTerm = NULL;
			}
		pszTrackGOTerms[-1] = Chr;
		NumGOTerms++;
		if(NumGOTerms >= MaxGOTerms)
			break;
		continue;
		}
	if(!InToken)			// if not already inside token then start token 
		{
		pGOTerm = pszTrackGOTerms-1;
		InToken = true;
		}
	}
if(InToken)
	{
	pszTrackGOTerms[-1] = '\0';
	if(pGOTerms != NULL)
		{
		strncpy(pGOTerms[NumGOTerms].szTerm,pGOTerm,cMaxOBOGOname);
		pGOTerms[NumGOTerms].szTerm[cMaxOBOGOname-1] = '\0';
		pGOTerms[NumGOTerms].pTerm = NULL;
		}
	pszTrackGOTerms[-1] = Chr;
	NumGOTerms++;
	}
return(NumGOTerms);
}


int
GenDensity(int WindowSize,CBEDfile *pBED,CFeatLoci *pFeatLoci,FILE * pGeneFile,
   		etOntologies Ontology,  // which ontology
		char *pszGoAssoc,		// gene to ID association file
		char *pszGOTerms,		// GO ontology file
		char *pszTrackGOTerms) // GO terms to track
{
int Rslt;
CGOAssocs *pAssocs = NULL;		// gene to ID association file
CGOTerms *pGOTerms = NULL;		// GO ontology file
tsGOTerm *pTerm;
int NumWindows;
int WindowStart;
int LociIdx;
tsDensityCnts *pDensityCnts;
tsDensityCnts *pCurCnts;
tsGOTermMap GOTerms[cMaxGOTerms];
unsigned int *pDensityBits;
int NumDensityBitsInts;
int NumTerms;
int TermIdx;
int DensityIdx;
int PrvDensityIdx;
int CurGeneID;
int PrvGeneID;
int GeneStart;
int GeneEnd;
char Strand;
int PrvGeneEnd;
int LastFeatEnd;
int LociStart;
int LociEnd;
char szChrom[128];
char szCurChrom[128];
char szFeatName[128];


if(Ontology != eONTUndefined)
	{
	NumTerms = ParseGOTerms(cMaxGOTerms,pszTrackGOTerms,GOTerms);
	if(NumTerms < 1)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error parsing GO terms of interest '%s'",pszTrackGOTerms);
		return(NumTerms);
		}
	pAssocs = new CGOAssocs();
	if((Rslt = pAssocs->Open(pszGoAssoc))!=eBSFSuccess)
		{
		while(pAssocs->NumErrMsgs())
			gDiagnostics.DiagOut(eDLFatal,gszProcName,pAssocs->GetErrMsg());
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open gene to GO:ID association file '%s'",pszGoAssoc);
		delete pAssocs;
		return(Rslt);
		}

	pGOTerms = new CGOTerms();
	if((Rslt = pGOTerms->Open(pszGOTerms))!=eBSFSuccess)
		{
		while(pGOTerms->NumErrMsgs())
			gDiagnostics.DiagOut(eDLFatal,gszProcName,pGOTerms->GetErrMsg());
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open GO:Term ontology file '%s'",pszGOTerms);
		delete pAssocs;
		delete pGOTerms;
		return(Rslt);
		}

	if((Rslt=pGOTerms->SetBkgndCnts(Ontology,// which class of ontologies to set counts for
		true,						// propagate counts from GO:Term into parent terms
		false,						// background counts proportional to gene lengths
		0,							// background counts are for which strand 0==both,1=='+' and 2=='-'
		pBED,						// gene BED file
		pAssocs))<eBSFSuccess)		// gene to ID association file
		{
		while(pGOTerms->NumErrMsgs())
			gDiagnostics.DiagOut(eDLFatal,gszProcName,pGOTerms->GetErrMsg());
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to set backgrond counts");
		delete pAssocs;
		delete pGOTerms;
		return(Rslt);
		}

	for(TermIdx = 0; TermIdx < NumTerms; TermIdx++)
		{
		pTerm = pGOTerms->LocateGOID(GOTerms[TermIdx].szTerm);
		if(pTerm == NULL)	// error if term not known to GO!
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unrecognised GO:term '%s'",GOTerms[TermIdx].szTerm);
			delete pAssocs;
			delete pGOTerms;
			return(eBSFerrGOID);
			}
		GOTerms[TermIdx].pTerm = pTerm;
		// ensure term root is as expected
		if(pTerm->RootOntology != Ontology)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"GO:term '%s' is not in %s ontology class",GOTerms[TermIdx].szTerm,Ontologies2Txt(Ontology));
			delete pAssocs;
			delete pGOTerms;
			return(eBSFerrGOID);
			}
		}
	}

if((pDensityCnts = new tsDensityCnts [cMaxNumDensityCnts])==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to set allocate memory to hold density counts");
	delete pAssocs;
	delete pGOTerms;
	return(Rslt);
	}

NumDensityBitsInts = 1 + (cMaxChromLen / (8 * sizeof(unsigned int)));
if((pDensityBits = new unsigned int [NumDensityBitsInts])==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to set allocate memory to hold density bits");
	delete pAssocs;
	delete pGOTerms;
	delete pDensityCnts;
	return(Rslt);
	}

fprintf(pGeneFile,"\"Chromosome\",\"WindowID\",\"Start\",\"NumGenes\",\"NumGeneBases\",\"NumFeats\",\"NumFeatBases\"");
for(TermIdx = 0; TermIdx < NumTerms; TermIdx++)
	fprintf(pGeneFile,",\"%s\"",GOTerms[TermIdx].szTerm);
fprintf(pGeneFile,"\n");	

PrvGeneID = 0;
do {
	szCurChrom[0] = '\0';	// processing new chromosome
	LastFeatEnd = 0;
	PrvDensityIdx = -1;
	PrvGeneEnd = -1;
	memset(pDensityCnts,0,sizeof(tsDensityCnts) * cMaxNumDensityCnts);
	memset(pDensityBits,0,sizeof(unsigned int) * NumDensityBitsInts);
	pGOTerms->ClearSampleCounts();
	while((CurGeneID = pBED->GetNextFeatureID(PrvGeneID)) > 0)
		{
		pBED->GetFeature(CurGeneID,szFeatName,szChrom,&GeneStart,&GeneEnd,NULL,&Strand);
		if(szCurChrom[0] != '\0' &&			// finished with current, on to a new chromosome?
			stricmp(szChrom,szCurChrom))
			break;							// time to process loci mapped on to the current chromosome
	
		if(szCurChrom[0] == '\0')
			strcpy(szCurChrom,szChrom);

		if(Strand != '-')					// strand specifies gene transcriptional direction 
			DensityIdx = GeneStart / WindowSize;
		else
			DensityIdx = GeneEnd / WindowSize;
		

		pCurCnts = &pDensityCnts[DensityIdx];
		pCurCnts->NumGenes += 1;

		SetBits(GeneStart,GeneEnd,pDensityBits);

		if(Ontology != eONTUndefined)
			{
			if(DensityIdx != 0 && DensityIdx != PrvDensityIdx)
				{
				// get GO counts for GO:Terms of interest here...
				for(TermIdx = 0; TermIdx < NumTerms; TermIdx++)
					{
					tsGOTermCnts *pTermCnts;

					pTermCnts = pGOTerms->GetExistingTermCnts(GOTerms[TermIdx].pTerm);
					if(pTermCnts != NULL)
						pCurCnts->GOCnts[TermIdx] = pTermCnts->SampleCnt;
					else
						pCurCnts->GOCnts[TermIdx] = 0;
					}
				pGOTerms->ClearSampleCounts();
				PrvDensityIdx = DensityIdx;
				}
			Add2GO(szFeatName,Ontology,true,pBED,pAssocs,pGOTerms);
			}
		PrvGeneID = CurGeneID;
		if(GeneEnd > LastFeatEnd)
			LastFeatEnd = GeneEnd;
		}

	// for each window determine number of gene bases
	NumWindows = 1 + (LastFeatEnd / WindowSize);
	WindowStart = 0;
	pCurCnts = pDensityCnts;
	for(DensityIdx = 0; DensityIdx < NumWindows; DensityIdx++,pCurCnts++,WindowStart += WindowSize)
		pCurCnts->NumGeneBases = CountBitsSet(WindowStart,WindowStart + WindowSize - 1,pDensityBits);

	// now find all loci on same chromosome just processed for genes
	if((LociIdx = pFeatLoci->Locate1stLoci(szCurChrom)) > 0)
		{
		memset(pDensityBits,0,sizeof(unsigned int) * NumDensityBitsInts);
		do {
			pFeatLoci->GetLoci(LociIdx,NULL,NULL,&LociStart,&LociEnd);
			DensityIdx = LociStart / WindowSize;
			pCurCnts = &pDensityCnts[DensityIdx];
			pCurCnts->NumFeats += 1;
			SetBits(LociStart,LociEnd,pDensityBits);
			if(LociEnd > LastFeatEnd)
				LastFeatEnd = LociEnd;
			}
		while((LociIdx = pFeatLoci->LocateNxtLociChrom(LociIdx))>0);

			// for each window determine number of feature bases
		NumWindows = 1 + (LastFeatEnd / WindowSize);
		WindowStart = 0;
		pCurCnts = pDensityCnts;
		for(DensityIdx = 0; DensityIdx < NumWindows; DensityIdx++,pCurCnts++,WindowStart += WindowSize)
			pCurCnts->NumFeatBases = CountBitsSet(WindowStart,WindowStart + WindowSize - 1,pDensityBits);
		}
		
	NumWindows = 1 + (LastFeatEnd / WindowSize);
	pCurCnts = pDensityCnts;
	for(DensityIdx = 0; DensityIdx < NumWindows; DensityIdx++,pCurCnts++)
		{
		fprintf(pGeneFile,"\"%s\",%d,%d,%d,%d,%d,%d",szCurChrom,DensityIdx,WindowSize * DensityIdx,
				pCurCnts->NumGenes,pCurCnts->NumGeneBases,pCurCnts->NumFeats,pCurCnts->NumFeatBases);
		for(TermIdx = 0; TermIdx < NumTerms; TermIdx++)
			fprintf(pGeneFile,",%d",pCurCnts->GOCnts[TermIdx]);
		fprintf(pGeneFile,"\n");	
		}
	}
while(CurGeneID > 0);

if(pDensityBits != NULL)
	delete pDensityBits;
if(pDensityCnts != NULL)
	delete pDensityCnts;
if(pAssocs != NULL)
	delete pAssocs;
if(pGOTerms != NULL)
	delete pGOTerms;
return(0);
}


int
Process(etProcMode Mode,		 // standard, distribution or density 
		int Intergenic,			// weighting to use if loci is intergenic
		int Upstream,			// weighting to use if loci is in upstream regulatory region
		int Intragenic,			// weighting to use if loci is intragenic
		int Dnstream,			// weighting to use if loci is in downstream regulatory region
		int RegLen,				// up/down stream regulatory region length
		int AssocDist,			// maximum association distance
		int ClustDist,			// treat as clustered if loci is less than this distance from end to start of next loci
		char *pszBED,			// gene BED file
		char Strand,			// loci are on this strand
		etLociFileFormat LociFormat, // expected LociFile format
		char *pszLociFile,		 // file containing genomic loci
		char *pszGeneFile,		 // gene + weighting file to generate
		etOntologies Ontology, // which ontology
		char *pszGoAssoc,		// gene to ID association file
		char *pszGOTerms,		// GO ontology file
		char *pszTrackGOTerms)	// GO terms to track
{
int Rslt;
FILE *pGeneFile;
CBEDfile *pBED = NULL;
CCSVFile *pCSVLoci = NULL;
CFeatLoci *pFeatLoci = NULL;
int RefField;
int ChromField;
int StartField;
int EndField;

int NumFields = 0;
int GeneCount = 0;
int NumLoci;

int RefID;
char *pszChrom;
int LociStart;
int LociEnd;
int NumUnclusteredLoci;

int NumAssoc = 0;
int NumWeighted = 0;
char szGene[cMaxFeatNameLen];
int Weighting;

switch(LociFormat) {
	case eLFormat9:		// 9 field CSV file, field 1: RefID, field 4: chromosome, field 5: start, field 6: end
		RefField = 1;
		ChromField = 4;
		StartField = 5;
		EndField = 6;
		break;
	case eLFormat3:		// 3 field CSV file, field 1: chromosome, field 2: start, field 3: end
		RefField = 0;
		ChromField = 1;
		StartField = 2;
		EndField = 3;
		break;
	case eLFUCSCBed:
		break;
	}

pBED = new CBEDfile();
if((Rslt = pBED->Open(pszBED))!=eBSFSuccess)
	{
	while(pBED->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,pBED->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open BED gene file '%s'",pszBED);
	delete pBED;
	return(Rslt);
	}


if((pGeneFile = fopen(pszGeneFile,"w+"))==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open/create results file %s error: %s",pszGeneFile,strerror(errno));
	return(eBSFerrOpnFile);
	}

pFeatLoci = new CFeatLoci();
NumUnclusteredLoci = 0;
if(LociFormat == eLFUCSCBed)
	{
	FILE *pBEDStream;
	int LineNum;
	int Cnt;
	char szLineBuff[cLineBuffLen];
	char szChrom[51];
	char *pTxt;
	Strand = '*';
	if((pBEDStream = fopen(pszLociFile,"r"))==NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open loci file %s error: %s",pszLociFile,strerror(errno));
		fclose(pGeneFile);
		delete pBED;
		delete pFeatLoci;
		return(Rslt);
		}
	LineNum = 0;
	RefID = 0;
	Rslt = eBSFerrParse; // assume the worst..., no features parsed
	while(fgets(szLineBuff,sizeof(szLineBuff),pBEDStream)!= NULL)
		{
		pTxt = TrimWhitespace(szLineBuff);
		if(*pTxt=='\0')	// simply slough lines which were just whitespace
			continue;
		Cnt = sscanf(pTxt," %50s %d %d",
			szChrom,&LociStart,&LociEnd);
		if(Cnt < 3)		// must be rubbish on this line or a comment line of some type
			{
			if(!NumUnclusteredLoci)	
				continue;	// no features processed, hope for the best and keep trying!
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to parse BED line '%s' from %s",szLineBuff,pszLociFile);
			Rslt = eBSFerrParse; // if at least one feature processed then treat rubbish as error
			break;
			}
		if((Rslt = NumUnclusteredLoci = pFeatLoci->AddLoci(++RefID,szChrom,LociStart,LociEnd)) < 0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error adding loci");
			break;
			}
		}
	fclose(pBEDStream);
	}
else
	{
	pCSVLoci = new CCSVFile();
	if((Rslt = pCSVLoci->Open(pszLociFile))!=eBSFSuccess)
		{
		while(pCSVLoci->NumErrMsgs())
			gDiagnostics.DiagOut(eDLFatal,gszProcName,pCSVLoci->GetErrMsg());
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open loci file '%s'",pszLociFile);
		fclose(pGeneFile);
		delete pBED;
		delete pCSVLoci;
		delete pFeatLoci;
		return(Rslt);
		}

	while((Rslt=pCSVLoci->NextLine()) > eBSFSuccess)
		{		
		NumFields = Rslt;
		if(RefField > 0)
			{
			Rslt = pCSVLoci->GetInt(RefField,&RefID);	
			if(Rslt < eBSFSuccess)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to read RefID start from field %d in line %d from %s",RefField,pCSVLoci->GetLineNumber(),pszLociFile);
				break;
				}
			}
		else
			RefID += 1;

		Rslt = pCSVLoci->GetText(ChromField,&pszChrom);	// field is expected to be the chromosome
		if(Rslt < eBSFSuccess)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to read chromosome from field %d in line %d from %s",ChromField,pCSVLoci->GetLineNumber(),pszLociFile);
			break;
			}
		Rslt = pCSVLoci->GetInt(StartField,&LociStart);	
		if(Rslt < eBSFSuccess)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to read loci start from field %d in line %d from %s",StartField,pCSVLoci->GetLineNumber(),pszLociFile);
			break;
			}

		Rslt = pCSVLoci->GetInt(EndField,&LociEnd);	
		if(Rslt < eBSFSuccess)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to read loci end from field %d in line %d from %s",EndField,pCSVLoci->GetLineNumber(),pszLociFile);
			break;
			}

		if((Rslt = NumUnclusteredLoci = pFeatLoci->AddLoci(RefID,pszChrom,LociStart,LociEnd)) < 0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error adding loci");
			break;
			}
		}
	pCSVLoci->Close();
	}
if(Rslt >= eBSFSuccess)
	{
	switch(Mode) {
		case ePMStandard:
			NumLoci = pFeatLoci->ClusterLoci(ClustDist);
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"Loci parsed: %d After clustering: %d",NumUnclusteredLoci,NumLoci);

			for(int Idx = 1; Idx <= NumLoci; Idx++)
				{
				pFeatLoci->GetLoci(Idx,NULL,&pszChrom,&LociStart,&LociEnd);
				szGene[0] = '\0';
				Weighting = 0;
				Rslt=Loci2Gene(Mode,Intergenic,Upstream,Intragenic,Dnstream,RegLen,Strand,
									AssocDist,pBED,pszChrom,LociStart,LociEnd,szGene,&Weighting);
				switch(Rslt) {
					case eBSFSuccess:
						NumAssoc += 1;
						if(Weighting != 0)
							{
							NumWeighted += 1;
							fprintf(pGeneFile,"%s,%d,%d\n",szGene,Weighting,1+LociEnd-LociStart);
							}
						continue;
					case eBSFerrFeature:	// unable to associate with gene/feature
					case eBSFerrChrom:
						Rslt = eBSFSuccess; // not unexpected so treat as success
						continue;
					default:				// must be some other error
						break;
					}
				break;
				}
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"Num associated: %d, Num weighted: %d",NumAssoc,NumWeighted);
			break;

		case ePMDensity:
			NumLoci = pFeatLoci->ClusterLoci(ClustDist);
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"Loci parsed: %d After clustering: %d",NumUnclusteredLoci,NumLoci);
			Rslt = GenDensity(AssocDist,pBED,pFeatLoci,pGeneFile,Ontology,pszGoAssoc,pszGOTerms,pszTrackGOTerms);
			break;

		case ePMDistribution:
			int GeneLen;
			int GeneDistance;
			int PrvLociDist;
			int NxtLociDist;
			int NxtLociStart;
			int NxtLociEnd;
			int PrvLociEnd;
			char *pszNxtChrom = NULL;
			char *pszPrvChrom = NULL;
			NumLoci = pFeatLoci->ClusterLoci(ClustDist);
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"Loci parsed: %d After clustering: %d",NumUnclusteredLoci,NumLoci);
			for(int Idx = 1; Idx <= NumLoci; Idx++)
				{
				pFeatLoci->GetLoci(Idx,&RefID,&pszChrom,&LociStart,&LociEnd);
				szGene[0] = '\0';
				Weighting = 0;
				Rslt=Loci2GeneA(Mode,Strand,pBED,pszChrom,LociStart,LociEnd,szGene,&GeneLen,&GeneDistance);
				
				if(Idx != NumLoci)
					{
					pFeatLoci->GetLoci(Idx+1,NULL,&pszNxtChrom,&NxtLociStart,&NxtLociEnd);
					if(!stricmp(pszChrom,pszNxtChrom))
						NxtLociDist = NxtLociStart - LociEnd;
					else
						NxtLociDist = -1;
					}
				else
					NxtLociDist = -1;

				if(Idx > 1)
					{
					if(!stricmp(pszChrom,pszPrvChrom))
						PrvLociDist = LociStart - PrvLociEnd;
					else
						PrvLociDist = -1;
					}
				else
					PrvLociDist = -1;

				if(Rslt >= eBSFSuccess)
					fprintf(pGeneFile,"%d,\"%s\",%d,%d,%d,\"%s\",%d,%d,%d,%d\n",RefID,pszChrom,LociStart,LociEnd,1 + LociEnd - LociStart,szGene,GeneLen,GeneDistance,PrvLociDist,NxtLociDist);

				PrvLociEnd = LociEnd;	
				pszPrvChrom = pszChrom;
				}
			break;


		}
	}

fclose(pGeneFile);
delete pCSVLoci;
pBED->Close();
delete pBED;
delete pFeatLoci;
return(Rslt);
}

int								// eBSFSuccess, eBSFerrFeature or eBSFerrChrom
Loci2Gene(etProcMode iMode,		 // processing mode
		int Intergenic,			// weighting to use if loci is intergenic
		int Upstream,			// weighting to use if loci is in upstream regulatory region
		int Intragenic,			// weighting to use if loci is intragenic
		int Dnstream,			// weighting to use if loci is in downstream regulatory region
		int RegLen,				// up/down stream regulatory region length 
		char Strand,			// loci are on this strand
		int AssocDist,			// maximum association distance
		CBEDfile *pBED,			// BED file containing genes to map loci onto
		char *pszChrom,			// chromosome
		int LociStart,			// loci start
		int LociEnd,			// loci end
		char *pszGene,			// where to return gene loci maps onto
		int *pWeighting)		// where to return weighting
{
int Ith;
int ChromID;
int FeatID;
int InGeneFeatID;
int UpstreamFeatID;
int DnstreamFeatID;

int RegLociStart;
int RegLociEnd;

int GeneStart;
int GeneEnd;
int InGeneStart;
int InGeneEnd;

int UpstreamDist;
int DnstreamDist;

*pWeighting = 0;
*pszGene = '\0';

// is the chromosome known?
if((ChromID = pBED->LocateChromIDbyName(pszChrom)) < 1)
	return(ChromID);	// most likely eBSFerrChrom
// processing may be strand specific...
pBED->SetStrand(Strand);

RegLociStart = LociStart < RegLen ? 0 : LociStart - RegLen;
RegLociEnd = LociEnd + RegLen;
FeatID = 0;
InGeneFeatID = 0;
UpstreamFeatID = 0;
DnstreamFeatID = 0;


// if loci +/-  RegLen is within at least one gene then...
if(pBED->InAnyFeature(ChromID,RegLociStart,RegLociEnd))
	{
	// iterate over all genes loci maps onto and use shortest gene completely containing loci
	Ith = 1;
	do {
		FeatID = pBED->LocateFeatureIDinRangeOnChrom(ChromID,RegLociStart,RegLociEnd,Ith++);
		if(FeatID > 0)
			{
			if(pBED->InFeature(FeatID,LociStart,LociEnd))
				{
				// could be multiple features hit so choose samllest feature which completely contains loci
				pBED->GetFeature(FeatID,NULL,NULL,&GeneStart,&GeneEnd);
				if(InGeneFeatID == 0)
					{
					InGeneFeatID = FeatID;
					InGeneStart = GeneStart;
					InGeneEnd = GeneEnd;
					}
				else
					{
					if(GeneStart <= LociStart && GeneEnd >= LociEnd)
						{
						if(GeneStart >= InGeneStart && GeneEnd <= InGeneEnd)
							{
							InGeneFeatID = FeatID;
							InGeneStart = GeneStart;
							InGeneEnd = GeneEnd;
							}
						}
					}
				UpstreamFeatID = 0;
				DnstreamFeatID = 0;
				}
			else			// else not in any gene, check up/down stream regulatory
				if(InGeneFeatID == 0 && pBED->In5Upstream(FeatID,LociStart,LociEnd,RegLen))
					{
					UpstreamFeatID = FeatID;
					DnstreamFeatID = 0;
					}
				else
					if(UpstreamFeatID == 0 && pBED->In3Dnstream(FeatID,LociStart,LociEnd,RegLen))
						DnstreamFeatID = FeatID;
			}	
		}
	while(FeatID > 0);
	}

if(InGeneFeatID)	// will be >= 1 if within a gene
	{
	FeatID = InGeneFeatID;
	*pWeighting = Intragenic;
	}
else				// not within gene so need to look at regulatory up/dn
	{
	if(UpstreamFeatID)
		{
		FeatID = UpstreamFeatID;
		*pWeighting = Upstream; 
		}
	else
		{
		if(DnstreamFeatID)
			{	
			FeatID = DnstreamFeatID;
			*pWeighting = Dnstream; 
			}
		}
	}

if(FeatID <= 0)		// if loci not up/dnstream or intragenic then look at intergenic
	{
	if(Intergenic > 0 && AssocDist > 0)
		{
		// loci must be intergenic, use closest gene as a proxy
		UpstreamFeatID = pBED->LocateFeatureBefore(ChromID,LociStart);
		DnstreamFeatID = pBED->LocateFeatureAfter(ChromID,LociEnd);
		if(UpstreamFeatID)
			{
			pBED->GetFeature(UpstreamFeatID,NULL,NULL,NULL,&GeneEnd);
			UpstreamDist = LociStart - GeneEnd;
			if(UpstreamDist > AssocDist)
				GeneEnd = -1;
			}
		else
			GeneEnd = -1;
		
		if(DnstreamFeatID)
			{
			pBED->GetFeature(DnstreamFeatID,NULL,NULL,&GeneStart,NULL);
			DnstreamDist = GeneStart - LociEnd;
			if(DnstreamDist > AssocDist)
				GeneStart = -1;
			}
		else
			GeneStart = -1;

		if(GeneEnd == -1 && GeneStart == -1)	// both < 0 if no feature on chrom or nearest gene more than AssocDist away
			return(eBSFerrFeature);

		// choose nearest gene
		if(GeneEnd == -1)
			FeatID = DnstreamFeatID;
		else
			if(GeneStart == -1)
				FeatID = UpstreamFeatID;
			else
				{
				if(UpstreamDist < DnstreamDist)
					FeatID = UpstreamFeatID;
				else
					FeatID = DnstreamFeatID;
				}
		*pWeighting = Intergenic;
		}
	}

if(FeatID <	1)
	return(eBSFerrFeature);	// unable to associate	

return(pBED->GetFeature(FeatID,pszGene));
}




int								// eBSFSuccess, eBSFerrFeature or eBSFerrChrom
Loci2GeneA(etProcMode Mode,		// processing mode
		char Strand,			// loci are on this strand
		CBEDfile *pBED,			// BED file containing genes to map loci onto
		char *pszChrom,			// chromosome
		int LociStart,			// loci start
		int LociEnd,			// loci end
		char *pszGene,			// where to return gene loci maps onto
		int *pGeneLen,			// where to return associated gene length
		int *pDistance)			// where to return distance away from gene if loci not intragenic
{
int Ith;
int ChromID;
int FeatID;
int InGeneFeatID;

int GeneStart;
int GeneEnd;
int InGeneStart;
int InGeneEnd;

int UpstreamDist;
int DnstreamDist;

*pGeneLen = 0;
*pDistance = 0;
*pszGene = '\0';

// is the chromosome known?
if((ChromID = pBED->LocateChromIDbyName(pszChrom)) < 1)
	return(ChromID);
// processing may be strand specific...
pBED->SetStrand(Strand);

FeatID = 0;
InGeneFeatID = 0;

// if loci is intragenic then...
if(pBED->InAnyFeature(ChromID,LociStart,LociEnd))
	{
	// iterate over all genes loci maps onto and use shortest gene completely containing loci

	Ith = 1;
	do {
		FeatID = pBED->LocateFeatureIDinRangeOnChrom(ChromID,LociStart,LociEnd,Ith++);
		if(FeatID < 1)
			continue;
				// could be multiple features hit so choose smallest feature which completely contains loci
		pBED->GetFeature(FeatID,NULL,NULL,&GeneStart,&GeneEnd);

		if(InGeneFeatID == 0)
			{
			InGeneFeatID = FeatID;
			InGeneStart = GeneStart;
			InGeneEnd = GeneEnd;
			}
		else
			{
			if(GeneStart <= LociStart && GeneEnd >= LociEnd)
				{
				if(GeneStart >= InGeneStart && GeneEnd <= InGeneEnd)
					{
					InGeneFeatID = FeatID;
					InGeneStart = GeneStart;
					InGeneEnd = GeneEnd;
					}
				}
			}
		}
	while(FeatID > 0);
	}

if(InGeneFeatID)
	FeatID = InGeneFeatID;

if(FeatID <= 0)		// if loci not intragenic 
	{
	// loci must be intergenic, use closest gene as a proxy
	int UpstreamFeatID = pBED->LocateFeatureBefore(ChromID,LociStart);
	int DnstreamFeatID = pBED->LocateFeatureAfter(ChromID,LociEnd);
	if(UpstreamFeatID)
		{
		pBED->GetFeature(UpstreamFeatID,NULL,NULL,&GeneStart,&GeneEnd);
		UpstreamDist = LociStart - GeneEnd;
		}
	else
		GeneEnd = -1;
		
	if(DnstreamFeatID)
		{
		pBED->GetFeature(DnstreamFeatID,NULL,NULL,&GeneStart,&GeneEnd);
		DnstreamDist = GeneStart - LociEnd;
		}
	else
		GeneStart = -1;

	if(GeneEnd == -1 && GeneStart == -1)	// both < 0 if no feature on chrom
		return(eBSFerrFeature);

		// choose nearest gene
	if(GeneEnd == -1)
		{
		FeatID = DnstreamFeatID;
		*pDistance = DnstreamDist;
		}
	else
		if(GeneStart == -1)
			{
			FeatID = UpstreamFeatID;
			*pDistance = UpstreamDist;
			}
		else
			{
			if(UpstreamDist < DnstreamDist)
				{
				FeatID = UpstreamFeatID;
				*pDistance = UpstreamDist;
				}
			else
				{
				FeatID = DnstreamFeatID;
				*pDistance = DnstreamDist;
				}
			}
	
	}
if(FeatID < 1)
	return(eBSFerrFeature);	// unable to associate	

pBED->GetFeature(FeatID,pszGene,NULL,&GeneStart,&GeneEnd);
*pGeneLen = 1 + GeneEnd - GeneStart;
return(eBSFSuccess);
}

int
Add2GO(char *pszGene,
	etOntologies Ontology, // which root ontologies to process
	bool bProp,				// propagate counts from GO:Term into parent terms
	CBEDfile *pBED,			// gene BED file
	CGOAssocs *pAssocs,		// gene to ID association file
	CGOTerms *pGOTerms)		// GO ontology file
{
int Rslt;
int NumGOIDs;
char *pszGOID;
char *pszGOIDs[1000];
char szGOIDs[10000];
char *pszConcatGOIDs;
int CntGOIDs;
int Cnt;

NumGOIDs = pAssocs->GetNumGOIDs(pszGene);
if(NumGOIDs > 0)
	{
	pszConcatGOIDs = szGOIDs;
	*pszConcatGOIDs = '\0';
	CntGOIDs = 0;
	for(Cnt = 0; Cnt < NumGOIDs; Cnt++)
		{
		pszGOID = pAssocs->GetGOID(pszGene,Cnt+1);
		if(pszGOID != NULL)
			{
			strcpy(pszConcatGOIDs,pszGOID);
			pszGOIDs[CntGOIDs++] = pszConcatGOIDs;
			pszConcatGOIDs += strlen(pszGOID) + 1;
			}
		}
	if(CntGOIDs)
		Rslt=pGOTerms->AddCount(Ontology,bProp,true,1,CntGOIDs,pszGOIDs);
	}
return(Rslt);
}

const int cBitsPerInt = sizeof(unsigned int) * 8;
const unsigned int cMSBitInt = 0x01 << (cBitsPerInt - 1); 

// CountBitsSet
// Counts number of bits set between Start and End inclusive in bit vector ptd at by pBitVector
int
CountBitsSet(int Start,int End,unsigned int *pBitVector)
{
int NumBitsSet;
int StartBit;
unsigned int Bits;
unsigned int LSBit;
int NumBits;
int TstBits;
int StartOfs;				
unsigned int *pBits;

if(Start < 0 || Start > End || pBitVector == NULL)
	return(-1);

NumBits = 1 + End - Start;
StartOfs = Start/cBitsPerInt;				
pBits = &pBitVector[StartOfs];
NumBitsSet = 0;
if(StartBit = (Start % cBitsPerInt))		// initial bit to test if set, 0 if on an int boundary
	{
	LSBit = 0x01 << StartBit;	// StartBit will be 1..31
	TstBits = cBitsPerInt - StartBit;
	while(--TstBits && --NumBits) 
		{
		if(*pBits & LSBit)
			NumBitsSet += 1;
		LSBit <<= 1;
		}
	pBits += 1;
	}

while(NumBits >= cBitsPerInt)
	{
	if((Bits = *pBits++) == (unsigned int)-1)
		NumBitsSet += cBitsPerInt;
	else
		if(Bits > 0)
			{
			LSBit = 0x00001;
			TstBits = cBitsPerInt;
			while(TstBits--)
				{
				if(Bits & LSBit)
					NumBitsSet += 1;
				LSBit <<= 1;
				}
			}
	NumBits -= cBitsPerInt;
	}

if(NumBits > 0)
	{
	LSBit = 0x00001;
	while(NumBits--)
		{
		if(*pBits & LSBit)
			NumBitsSet += 1;
		LSBit <<= 1;
		}
	}
return(NumBitsSet);
}

// SetBits
// Sets all bits between Start and End inclusive in bit vector ptd at by pBitVector
// Returns number of bits set or < 0 if parameter errors
int
SetBits(int Start,int End,unsigned int *pBitVector)
{
int StartBit;
unsigned int Bits;
unsigned int LSBit;
int NumBits;
int NumSet;
int StartOfs;				
unsigned int *pBits;

if(Start < 0 || Start > End || pBitVector == NULL)
	return(-1);

NumSet = NumBits = 1 + End - Start;
StartOfs = Start/cBitsPerInt;				
pBits = &pBitVector[StartOfs];

if(StartBit = (Start % cBitsPerInt))		// initial bit to set, 0 if on an int boundary
	{
	Bits = LSBit = 0x01 << StartBit;	// StartBit will be 1..BitsPerInt-1
	while(--NumBits && !(Bits & cMSBitInt)) 
		{
		Bits <<= 1;
		Bits |= LSBit;
		}
	*pBits++ |= Bits;
	}

while(NumBits >= cBitsPerInt)
	{
	*pBits++ |= (unsigned int)-1;
	NumBits -= cBitsPerInt;
	}

if(NumBits > 0)
	{
	Bits = 0x01;
	while(--NumBits)
		{
		Bits <<= 1;
		Bits |= 1;
		}
	*pBits |= Bits;
	}
return(NumSet);
}

