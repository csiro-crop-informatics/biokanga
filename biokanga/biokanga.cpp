/*
 * CSIRO Open Source Software License Agreement (GPLv3)
 * Copyright (c) 2017, Commonwealth Scientific and Industrial Research Organisation (CSIRO) ABN 41 687 119 230.
 * See LICENSE for the complete license information (https://github.com/csiro-crop-informatics/biokanga/LICENSE)
 * Contact: Alex Whan <alex.whan@csiro.au>
 */

// biokanga.cpp : Defines the entry point for the console application.
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

#include "biokanga.h"


// when required to static compile
#ifdef _LINKSTATIC_
const char *cpszProgVer = "4.3.6S";		// increment with each release
const char *cpszProcOverview = "BioKanga NGS Processing Toolset - Static Linked";
#else
const char *cpszProgVer = "4.3.6";		// increment with each release
const char *cpszProcOverview = "BioKanga NGS Processing Toolset";
#endif

// Subprocesses 
extern int Blitz(int argc, char* argv[]);
extern int LocateROI(int argc, char* argv[]);
extern int RemapLoci(int argc, char* argv[]);
extern int FilterSAMAlignments(int argc, char* argv[]);
extern int SimReads(int argc, char* argv[]);
extern int fasta2nxx(int argc, char* argv[]);
extern int ArtefactReduce(int argc, char* argv[]);
extern int ReadsetDists(int argc, char* argv[]);
extern int Assemble(int argc, char* argv[]);
extern int ScaffoldContigs(int argc, char* argv[]);
extern int kanga(int argc, char* argv[]);
extern int kangax(int argc, char* argv[]);
extern int kangade(int argc, char* argv[]);
extern int maploci2features(int argc, char* argv[]);
extern int genpseudogenome(int argc, char* argv[]);
extern int gensnpmarkers(int argc, char* argv[]);
extern int genkmarkers(int argc, char* argv[]);
extern int GenMarkerSeqs(int argc, char* argv[]);
extern int gendeseq(int argc, char* argv[]);
extern int Markers2SQLite(int argc, char* argv[]);
extern int SNPs2SQLite(int argc, char* argv[]);
extern int DE2SQLite(int argc, char* argv[]);
extern int PSL2SQLite(int argc, char* argv[]);
extern int kmermarkers(int argc, char* argv[]);
extern int fastaextract(int argc, char* argv[]);
extern int mergeoverlaps(int argc, char* argv[]);
extern int pescaffold(int argc, char* argv[]);
extern int SSRdiscovery(int argc, char* argv[]);
extern int AlignsBootstrap(int argc, char* argv[]);


// inplace text cleaning; any leading/trailing or internal quote characters are removed; excessive whitespace is reduced to single
char *
RemoveQuotes(char *pszRawText)
{
char *pSrcChr;
char *pDstChr;
bool bInSpace;
char Chr;
CUtility::TrimQuotedWhitespcExtd(pszRawText);
pSrcChr = pszRawText;
pDstChr = pSrcChr;
bInSpace = false;
while((Chr = *pSrcChr++)!= '\0')
	{
	if(Chr == '\'' || Chr == '"')
		continue;
	if(Chr == ' ' || Chr == '\t')
		{
		if(bInSpace)
			continue;
		bInSpace = true;
		}
	else
		bInSpace = false;
	*pDstChr++ = Chr;
	}
*pDstChr = '\0';
return(pszRawText);
}

tsSubProcess SubProcesses[] = {
	{ "simreads", "Simulate NGS Reads", "Generate simulated NGS readsets", SimReads },
	{ "ngsqc", "NGS Reads QC", "\tProcess NGS reads and report quality scores with compositional distributions", ReadsetDists },
	{ "fasta2nxx", "Fasta Nxx", "Generate N10..N90 over Fasta sequences", fasta2nxx },
	{"filter","Filter NGS Reads","Filter NGS reads for sequencer errors and/or exact duplicates",ArtefactReduce },
	{"assemb","de Novo Assemble","de Novo assemble filtered reads into contigs",Assemble },
	{"scaffold","Scaffold Contigs","Scaffold de Novo assembled contigs",ScaffoldContigs },
	{"index","Index Assembly","\tGenerate index over genome assembly or sequences",kangax },
	{"kmarkers","K-Mer Markers","NGS reads alignment-less K-mer derived marker sequences generation",genkmarkers},
	{"prekmarkers","K-Mer Prefix Markers","NGS reads alignment-less prefix K-mer derived marker sequences generation",kmermarkers},
	{"pseudogenome","Generate Pseudo-Genome","Concatenate sequences to create pseudo-genome assembly",genpseudogenome},
	{"align","Align NGS reads","\tAlign NGS reads to indexed genome assembly or sequences",kanga},
	{"pescaffold","PE Scaffold","Scaffold assembly contigs using PE read alignments",pescaffold},
	{"ssr","SSR Discovery","\tIdentify SSRs in multifasta sequences",SSRdiscovery},
	{"maploci","Map Loci to Features","Map aligned reads loci to known features",maploci2features},
	{"rnade","RNA-seq Differential Expression","\tRNA-seq differential expression analyser with optional Pearsons generation",kangade},
	{"gendeseq","DESeq","Generate tab delimited counts file for input to DESeq or EdgeR",gendeseq},
	{"xfasta","Extract Fasta","Extract fasta sequences from multifasta file",fastaextract},
	{"mergeoverlaps","Merge PE Overlaps","Merge PE short insert overlap reads",mergeoverlaps},
	{"snpmarkers","SNP Markers","SNP alignment derived marker sequences identification",gensnpmarkers},
	{"markerseqs","Marker Seqs","Generate marker sequences from SNP loci",GenMarkerSeqs},
	{"blitz", "Blat like Local align", "\tBlat like local align genomic sequences", Blitz },
	{"remaploci","Remap Alignment Loci", "Remap alignment loci", RemapLoci },
	{"filtchrom","Filter SAM/BAM by chrom", "Filter SAM/BAM alignments by chromosome", FilterSAMAlignments },
	{"locateroi","Locate Regions of Interest", "Locate and report regions of interest", LocateROI },
	{"alignsbs","Alignment Bootstraps", "Alignments bootstrapper", AlignsBootstrap },
	{"psl2sqlite","SQLite Blat Alignments","Generate SQLite Blat alignment Database from Blat generated PSL alignments",PSL2SQLite},
	{"snpm2sqlite","SQLite SNP Markers","Generate SQLite Marker Database from SNP markers  ",Markers2SQLite},
	{"snps2sqlite","SQLite SNPs","Generate SQLite SNP Database from aligner identified SNPs",SNPs2SQLite},
	{"de2sqlite","SQLite DE","Generate SQLite DE Database from RNA-seq DE",DE2SQLite}
};
const int cNumSubProcesses = (sizeof(SubProcesses) / sizeof(tsSubProcess));

CStopWatch gStopWatch;
CDiagnostics gDiagnostics;				// for writing diagnostics messages to log file
CSQLiteSummaries gSQLiteSummaries;		// for writing processing result summaries to SQLite database
int	gExperimentID = 0;					// SQLite experiment identifier
int gProcessID = 0;						// SQLite process identifier
int	gProcessingID = 0;					// SQLite processing identifier

char gszProcName[_MAX_FNAME];			// this processes name
tsSubProcess *gpszSubProcess;			// selected subprocess


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
	return _T("biokanga");
}
// end of str library required code
#endif


void
GiveHelpSubProcesses(char *pszProcOverview)
{
int Idx;
printf("Help for %s Version %s\n",gszProcName,cpszProgVer);
printf("  %s\n",pszProcOverview);
printf("  Please specify one of the following subprocesses:\n");
for(Idx = 0; Idx < cNumSubProcesses; Idx++)
	printf("  %s\t\t%s\n",SubProcesses[Idx].pszName,SubProcesses[Idx].pszFullDescr);
printf("To obtain parameter help on any subprocess then enter that subprocess name e.g:\n%s %s -h\n",gszProcName,SubProcesses[0].pszName);
}

int
IsValidSubprocess(char *pszSubProcess)
{
int Idx;
tsSubProcess *pSubProcess;
pSubProcess = SubProcesses;
for(Idx = 0; Idx < cNumSubProcesses; Idx++,pSubProcess++)
	if(!stricmp(pSubProcess->pszName,pszSubProcess))
		return(Idx+1);
return(-1);
}

int
ExecSubProcess(int SubProcID,	// which subprocess as returned by IsValidSubprocess()
		int argc,char *argv[])	// parameters for this subprocess
{
int Idx;
char *pArg;
char *pszSubProc;
char szParam[1024];
if(SubProcID < 0 || SubProcID > cNumSubProcesses)
	{
	printf("\nInvalid SubProcID %d",SubProcID);
	return(-1);
	}
// process parameters and remove subprocess specifier
pArg = argv[1];
for(Idx = 1; Idx < argc; Idx++,pArg++)
	{
	if(pArg[0] == '-')
		{
		if(pArg[1] != 'p')
			continue;
		pszSubProc = &pArg[2];
		}
	else
		pszSubProc = pArg;
	strncpy(szParam,pszSubProc,sizeof(szParam) - 1);
	szParam[sizeof(szParam) - 1] = '\0';
	CUtility::TrimQuotedWhitespcExtd(szParam);
	if(!stricmp(szParam,SubProcesses[SubProcID-1].pszName))
		break;
	}
if(Idx != argc) // if located the subprocess specifier then remove from parameters
	{
	while(Idx < argc)
		{
		pArg = argv[Idx];
		argv[Idx] = argv[Idx+1];
		argv[++Idx] = pArg;
		}
	argc -= 1;
	}
gpszSubProcess = &SubProcesses[SubProcID-1];
return((*SubProcesses[SubProcID-1].SubFunct)(argc,argv));
}

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

// check if user has specified the subprocess required, each subprocess has it's own set of parameters
// the subprocess requested is either the first word immediately following the process name or the '-p<subprocess' option
// iterate all command line parameters checking for subprocess
int Rslt;
int Idx;
int SubProcID;
char szSubProc[1024];
char *pszSubProc;
char *pArg = (char *)argv[1];
SubProcID = 0;
for(Idx = 1; Idx < argc; Idx++)
	{
	if(pArg[0] == '-')
		{
		if(pArg[1] != 'p')
			continue;
		pszSubProc = &pArg[2];
		}
	else
		pszSubProc = pArg;
	strncpy(szSubProc,pszSubProc,sizeof(szSubProc) - 1);
	szSubProc[sizeof(szSubProc) - 1] = '\0';
	CUtility::TrimQuotedWhitespcExtd(szSubProc);
	if((SubProcID = IsValidSubprocess(szSubProc)) > 0)
		break;
	}
if(SubProcID > 0)
	Rslt = ExecSubProcess(SubProcID,argc,(char **)argv);
else
	{
	GiveHelpSubProcesses((char *)cpszProcOverview);
	Rslt = -1;
	}

return(Rslt);
}
