// bed2csv.cpp : Defines the entry point for the console application.
#include "stdafx.h"
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#if _WIN32

#include "../libbiokanga/commhdrs.h"
#else
#include "../libbiokanga/commhdrs.h"
#endif

const unsigned int cProgVer = 302;			// increment with each release


// processing mode
typedef enum eProcMode {
	ePMdefault,						// Process biobed (.bsb) or ascii BED (.bed) into csv loci file
	ePMUCSCGG						// Process biobed (.bed) or ascii BED (.bed) into UCSC genome graph csv
} etProcMode;

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


char *ProcMode2Txt(etProcMode ProcMode);
char *Region2Txt(etBEDRegion Region);
int 
Process(etProcMode ProcMode,		// processing mode
			etBEDRegion Region,		// which regions
			char *pszType,			// this as the element type in generated file
			char *pszSpecies,		// this species as the reference species in generated file
			char *pszRelSpecies,		// this as the relative species in generated file
			char *pszInputFile,		// core loci Biobed (.bsb) file to process
			char *pszOutputFile);	// where to write CSV loci

int
Output(etProcMode ProcMode,
	   int hOutFile,				// file handle to use for output
	   int RefID,					// unique reference identifier
	   char *pszType,				// element type
	   char *pszSpecies,			// reference species
	   char *pszChrom,				// chromosome
	   int StartLoci,				// start loci on chromsome 
	   int EndLoci,					// end loci on chromosome
	   char *pszRelSpecies,			// relative species
	   int Region,					// region
	   char Strand,					// strand
	   int Score);					// associated score 0..1000

CStopWatch gStopWatch;
CDiagnostics gDiagnostics;				// for writing diagnostics messages to log file
char gszProcName[_MAX_FNAME];			// process name

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

int iProcMode;
char szOutputFile[_MAX_PATH];
char szInputFile[_MAX_PATH];
char szRefSpecies[cMaxGeneNameLen];		// reference species
char szType[cMaxGeneNameLen];			// element type
char szRelSpecies[cMaxGeneNameLen];		// relative species


int iRegion;			// genomic regions to include: 0:ANY,1:Exons,2:Introns,3:CDS,4:UTRs,5:5'UTR,6:3'UTR (default = ANY)

// command line args
struct arg_lit  *help    = arg_lit0("hH","help",                "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_int *ScreenLogLevel=arg_int0("S", "ScreenLogLevel",	"<int>","Level of diagnostics written to screen 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",	"diagnostics log file");

struct arg_int  *ProcMode = arg_int0("m","procmode","<int>",		"Process Mode, 0: biobed (.bsb) or bed (.bed) into csv loci file, 1: into csv UCSC Genome Graph format");
struct arg_file *InFile   = arg_file1("i",NULL,"<file>",			"input from biobed or raw bed file");
struct arg_file *OutFile  = arg_file1("o",NULL,"<file>",			"output to CSV loci file");

struct arg_str *Type = arg_str1("t","elementtype","<string>",		"text to use as the element type in generated CSV loci file");
struct arg_str *RefSpecies = arg_str1("r","refspecies","<string>",	"reference species to use in generated CSV loci file");
struct arg_str *RelSpecies= arg_str1("R","relspecies","<string>",	"relative species to use in generated CSV loci file");
struct arg_int *Region    = arg_int0("g","genomicregion","<int>",	"Process regions 0:ALL,1:Intergenic,2:Exons,3:Introns,4:CDS,5:UTRs,6:5'UTR,7:3'UTR (default = ALL)");

struct arg_end *end = arg_end(20);

void *argtable[] = {help,version,FileLogLevel,ScreenLogLevel,LogFile,
					ProcMode,InFile,OutFile,Region,Type,RefSpecies,RelSpecies,
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

	iProcMode = ProcMode->count ? ProcMode->ival[0] : ePMdefault;
	if(iProcMode < ePMdefault || iProcMode > ePMUCSCGG)
		{
		printf("\nError: Requested processing mode '-m%d' not supported",iProcMode);
		exit(1);
		}


	iRegion = Region->count ? Region->ival[0] : eMEGRAny;	// default as being any region
	if(iRegion < eMEGRAny)
		{
		printf("\nSpecified region '-g%d' < 0, assuming you meant ANY region",iRegion);
		iRegion = eMEGRAny;
		}
	else
		{
		if(iRegion > eMEG3UTR)
			{
			printf("\nSpecified region '-n%d' > %d, assuming you meant 3'DS",iRegion,eMEG3UTR);
			iRegion = eMEG3UTR;
			}
		}
	

	strncpy(szInputFile,InFile->filename[0],sizeof(szInputFile));
	szInputFile[sizeof(szInputFile)-1] = '\0';
	strncpy(szOutputFile,OutFile->filename[0],sizeof(szOutputFile));
	szOutputFile[sizeof(szOutputFile)-1] = '\0';

	strncpy(szRefSpecies,RefSpecies->sval[0],sizeof(szRefSpecies));
	szRefSpecies[sizeof(szRefSpecies)-1] = '\0';

	strncpy(szRelSpecies,RelSpecies->sval[0],sizeof(szRelSpecies));
	szRelSpecies[sizeof(szRelSpecies)-1] = '\0';

	strncpy(szType,Type->sval[0],sizeof(szType));
	szType[sizeof(szType)-1] = '\0';


		// now that command parameters have been parsed then initialise diagnostics log system
	if(!gDiagnostics.Open(szLogFile,(etDiagLevel)iScreenLogLevel,(etDiagLevel)iFileLogLevel,true))
		{
		printf("\nError: Unable to start diagnostics subsystem.");
		if(szLogFile[0] != '\0')
			printf(" Most likely cause is that logfile '%s' can't be opened/created",szLogFile);
		exit(1);
		}

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Version: %d.%2.2d Processing parameters:",cProgVer/100,cProgVer%100);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Processing mode: %s",ProcMode2Txt((etProcMode)iProcMode));
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Process biobed input file %s",szInputFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Output file: '%s'",szOutputFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Include Region: %s",Region2Txt((etBEDRegion)iRegion));
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Type: '%s'",szType);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Ref Species: '%s'",szRefSpecies);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Rel Species: '%s'",szRelSpecies);

	gStopWatch.Start();
#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	Rslt = Process((etProcMode)iProcMode,		// processing mode
					(etBEDRegion)iRegion,		// which regions
					szType,				// element type
					szRefSpecies,		// reference species
					szRelSpecies,		// relative species
					szInputFile,		// core loci Biobed (.bsb) file to process
					szOutputFile);		// where to write CSV loci

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
ProcMode2Txt(etProcMode ProcMode)
{
switch(ProcMode) {
	case ePMdefault:
		return((char *)"Process biobed (.bsb) or ascii BED (.bed) into csv loci file");
	case ePMUCSCGG:
		return((char *)"Process biobed (.bed) or ascii BED (.bed) into UCSC genome graph csv file");
	default:
		break;
	}
return((char *)"Unsupported");
}

char *
Region2Txt(etBEDRegion Region)
{
switch(Region) {
	case eMEGRAny:		// process any region
		return((char *)"All");

	case eMEGRIntergenic:	// only process intergenic
		return((char *)"Intergenic");

	case eMEGRExons:	// only process exons
		return((char *)"EXONS");

	case eMEGRIntrons:	// only process introns
		return((char *)"INTRONS");

	case eMEGRCDS:		// only process CDSs
		return((char *)"CDS");

	case eMEGUTR:		// only process UTRs
		return((char *)"UTR");

	case eMEG5UTR:		// only process 5'UTRs
		return((char *)"5'UTR");

	case eMEG3UTR:		// only process 3'UTRs
		return((char *)"3'UTR");

	default:
		break;
	}
return((char *)"Unsupported");
}

int 
Process(etProcMode ProcMode,	// processing mode
			etBEDRegion Region,		// which regions
			char *pszType,			// this as the element type in generated file
			char *pszSpecies,		// this species as the reference species in generated file
			char *pszRelSpecies,		// this as the relative species in generated file
			char *pszInputFile,		// core loci Biobed (.bsb) file to process
			char *pszOutputFile)	// where to write CSV loci
{
CBEDfile BEDfile;
CAlignValidate AlignValidate;
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
char szRelComp[128];
char Strand;
int Rslt;
int RefID;
int IntergenicStart;
int Idx;
int hOutFile;

#ifdef _WIN32
if((hOutFile = open(pszOutputFile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
if((hOutFile = open(pszOutputFile, O_RDWR | O_CREAT | O_TRUNC, S_IREAD | S_IWRITE ))==-1)
#endif
	{
	gDiagnostics.DiagOutMsgOnly(eDLFatal,"Unable to create or truncate  output file %s error: %s",pszOutputFile,strerror(errno));
	return(-1);
	}

if((Rslt=BEDfile.Open(pszInputFile,eBTAnyBed)) !=eBSFSuccess)
	{
	gDiagnostics.DiagOutMsgOnly(eDLFatal,"Unable to open '%s' for processing",pszInputFile);
	close(hOutFile);
	return(-1);
	}

if(!BEDfile.ContainsGeneDetail() && Region > eMEGRIntergenic)			// returns true if file contains gene detail (utr/cds/intron etc)
	{
	gDiagnostics.DiagOutMsgOnly(eDLFatal,"biobed file %s does not contain gene regions",pszInputFile);
	close(hOutFile);
	return(-1);
	}

szPrevChrom[0] = '\0';
CurFeatureID = 0;
IntergenicStart = 0;
RefID = 0;
while((CurFeatureID = BEDfile.GetNextFeatureID(CurFeatureID)) > 0)
	{
	BEDfile.GetFeature(CurFeatureID,	// feature instance identifier
				szFeatName,				// where to return feature name
				szChrom,				// where to return chromosome name
				&StartLoci,				// where to return feature start on chromosome (0..n) 
				&EndLoci,				// where to return feature end on chromosome
 				&Score,					// where to return score
 				&Strand);				// where to return strand

	sprintf(szRelComp,"%s_%s",pszRelSpecies,szFeatName);

	if(CurFeatureID == 1 || stricmp(szChrom,szPrevChrom))	// if new chromosome then reset IntergenicStart
		{
		strcpy(szPrevChrom,szChrom);
		IntergenicStart = 0;
		}

	if(Region != eMEGRAny)
		{
		NumExons = BEDfile.GetNumExons(CurFeatureID);		// returns number of exons - includes UTRs + CDS
		NumIntrons = BEDfile.GetNumIntrons(CurFeatureID);
		CDSstart = StartLoci + BEDfile.GetCDSStart(CurFeatureID);		// returns relative start offset of CDS - NOTE add to '+' strand gene start, subtract on '-' strand gene start
		CDSend = StartLoci + BEDfile.GetCDSEnd(CurFeatureID);			// returns relative end offset of CDS - NOTE add to '+' strand gene start, subtract on '-' strand gene start
		}

	switch(Region) {
		case eMEGRAny:					// process any region
			Output(ProcMode,hOutFile,++RefID,pszType,pszSpecies,szChrom,StartLoci,EndLoci,szRelComp,0,Strand,Score);
			continue;

		case eMEGRIntergenic:	// only process intergenic
			if(IntergenicStart < StartLoci)
				Output(ProcMode,hOutFile,++RefID,pszType,pszSpecies,szChrom,IntergenicStart,StartLoci-1,szRelComp,0,Strand,Score);
			if(IntergenicStart <= EndLoci)
				IntergenicStart = EndLoci+1;
			continue;

		case eMEGRExons:		// only process exons
			for(Idx = 1; Idx <= NumExons; Idx++)
				{
				StartLoci = BEDfile.GetExonStart(CurFeatureID,Idx);
				EndLoci   = BEDfile.GetExonEnd(CurFeatureID,Idx);
				if(StartLoci <= EndLoci)
					Output(ProcMode,hOutFile,++RefID,pszType,pszSpecies,szChrom,StartLoci,EndLoci,szRelComp,0,Strand,Score);
				}
			continue;

		case eMEGRIntrons:		// only process introns
			for(Idx = 1; Idx <= NumIntrons; Idx++)
				{
				StartLoci = BEDfile.GetIntronStart(CurFeatureID,Idx);
				EndLoci = BEDfile.GetIntronEnd(CurFeatureID,Idx);
				if(StartLoci <= EndLoci)
					Output(ProcMode,hOutFile,++RefID,pszType,pszSpecies,szChrom,StartLoci,EndLoci,szRelComp,0,Strand,Score);
				}
			continue;

		case eMEGRCDS:			// only process CDSs
			for(Idx = 1; Idx <= NumExons; Idx++)
				{
				StartLoci = BEDfile.GetExonStart(CurFeatureID,Idx);
				EndLoci   = BEDfile.GetExonEnd(CurFeatureID,Idx);
				if(EndLoci < CDSstart || StartLoci > CDSend)
					continue;
				if(StartLoci < CDSstart)
					StartLoci = CDSstart;
				if(EndLoci > CDSend)
					EndLoci = CDSend;
				if(StartLoci <= EndLoci)
					Output(ProcMode,hOutFile,++RefID,pszType,pszSpecies,szChrom,StartLoci,EndLoci,szRelComp,0,Strand,Score);
				}
			continue;

		case eMEGUTR:			// only process UTRs - single exon may have both 5' and 3' UTRs
			for(Idx = 1; Idx <= NumExons; Idx++)
				{
				StartLoci = BEDfile.GetExonStart(CurFeatureID,Idx);
				EndLoci   = BEDfile.GetExonEnd(CurFeatureID,Idx);
				if(EndLoci <= CDSend && StartLoci >= CDSstart) // is exon CDS only?
					continue;

				// check if 5' UTR 
				if(StartLoci < CDSstart)
					{
					if(EndLoci >= CDSstart)
						Output(ProcMode,hOutFile,++RefID,pszType,pszSpecies,szChrom,StartLoci,CDSstart-1,szRelComp,0,Strand,Score);
					else
						Output(ProcMode,hOutFile,++RefID,pszType,pszSpecies,szChrom,StartLoci,EndLoci,szRelComp,0,Strand,Score);
					}
					
				// check if 3'UTR
				if(EndLoci > CDSend)
					{
					if(StartLoci <= CDSend)
						StartLoci = CDSend+1;
					Output(ProcMode,hOutFile,++RefID,pszType,pszSpecies,szChrom,StartLoci,EndLoci,szRelComp,0,Strand,Score);
					}
				}
			continue;

		case eMEG5UTR:			// only process 5'UTRs - strand sensitive
			for(Idx = 1; Idx <= NumExons; Idx++)
				{
				StartLoci = BEDfile.GetExonStart(CurFeatureID,Idx);
				EndLoci   = BEDfile.GetExonEnd(CurFeatureID,Idx);
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
				Output(ProcMode,hOutFile,++RefID,pszType,pszSpecies,szChrom,StartLoci,EndLoci,szRelComp,0,Strand,Score);
				}
			continue;

		case eMEG3UTR:			// only process 3'UTRs  - strand sensitive
			for(Idx = 1; Idx <= NumExons; Idx++)
				{
				StartLoci = BEDfile.GetExonStart(CurFeatureID,Idx);
				EndLoci   = BEDfile.GetExonEnd(CurFeatureID,Idx);
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
				Output(ProcMode,hOutFile,++RefID,pszType,pszSpecies,szChrom,StartLoci,EndLoci,szRelComp,0,Strand,Score);
				}
			continue;
		}
	}

BEDfile.Close(false);
close(hOutFile);
return(0);
}

int
Output(etProcMode ProcMode,
	   int hOutFile,				// file handle to use for output
	   int RefID,					// unique reference identifier
	   char *pszType,				// element type
	   char *pszSpecies,			// reference species
	   char *pszChrom,				// chromosome
	   int StartLoci,				// start loci on chromsome 
	   int EndLoci,					// end loci on chromosome
	   char *pszRelSpecies,			// relative species
	   int Region,					// region
	   char Strand,					// strand
	   int Score)					// associated score 0..1000
{
int Rslt;
char szBuff[4096];
int Len;

switch(ProcMode) {
	case ePMdefault:
		Len = sprintf(szBuff,"%d,\"%s\",\"%s\",\"%s\",%d,%d,%d,\"%s\",%d,\"%c\"\n",
								RefID,pszType,pszSpecies,pszChrom,StartLoci,EndLoci,EndLoci-StartLoci+1,pszRelSpecies,Region,Strand);
		if((Rslt = write(hOutFile,szBuff,Len))!=Len)
			{
			gDiagnostics.DiagOutMsgOnly(eDLFatal,"Write to file failed: %s",strerror(errno));
			return(-1);
			}
		break;

	case ePMUCSCGG:
		Len = sprintf(szBuff,"%s,%d,%1.4f\n",
								pszChrom,StartLoci,((double)Score/1000.0));
		if((Rslt = write(hOutFile,szBuff,Len))!=Len)
			{
			gDiagnostics.DiagOutMsgOnly(eDLFatal,"Write to file failed: %s",strerror(errno));
			return(-1);
			}
		break;
	}
return(Rslt);
}

