// Copyright 2013 CSIRO  ( http://www.csiro.au/ ) 
//   Licensed under the Apache License, Version 2.0 (the "License");
//   you may not use this file except in compliance with the License.
//   You may obtain a copy of the License at
//       http://www.apache.org/licenses/LICENSE-2.0
//
//   Unless required by applicable law or agreed to in writing, software
//   distributed under the License is distributed on an "AS IS" BASIS,
//   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//   See the License for the specific language governing permissions and
//   limitations under the License
//   Please contact stuart.stephen@csiro.au for support or 
//   to submit modifications to this source

// Aligner.cpp : contains the CAligner class implementation
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

#define _DISNPS_ 1    // experimental, exploring the utility of generating DiSNPs

#include "biokanga.h"
#include "Aligner.h"

#include "../libbiokanga/bgzf.h"
#include "../libbiokanga/SAMfile.h"

sNAR CAligner::m_NARdesc[] = {{eNARUnaligned,(char *)"NA",(char *)"Not processed for alignment"},
	{eNARAccepted,(char *)"AA",(char *)"Alignment accepted"},
	{eNARNs,(char *)"EN",(char *)"Excessive indeterminate (Ns) bases"},
	{eNARNoHit,(char *)"NL",(char *)"No potential alignment loci"},
	{eNARMMDelta,(char *)"MH",(char *)"Mismatch delta (minimum Hamming) criteria not met"},
	{eNARMultiAlign,(char *)"ML",(char *)"Aligned to multiloci"},
	{eNARTrim,(char *)"ET",(char *)"Excessively end trimmed"},
	{eNARSpliceJctn,(char *)"OJ",(char *)"Aligned as orphaned splice junction"},
	{eNARmicroInDel,(char *)"OM",(char *)"Aligned as orphaned microInDel"},
	{eNARPCRdup,(char *)"DP",(char *)"Duplicate PCR"},
	{eNARNonUnique,(char *)"DS",(char *)"Duplicate read sequence"},
	{eNARChromFilt,(char *)"FC",(char *)"Aligned to filtered target sequence"},
	{eNARRegionFilt,(char *)"PR",(char *)"Aligned to a priority region"},
	{eNARPEInsertMin,(char *)"UI",(char *)"PE under minimum insert size"},
	{eNARPEInsertMax,(char *)"OI",(char *)"PE over maximum insert size"},
	{eNARPENoHit,(char *)"UP",(char *)"PE partner not aligned"},
	{eNARPEStrand,(char *)"IS",(char *)"PE partner aligned to inconsistent strand"},
	{eNARPEChrom,(char *)"IT",(char *)"PE partner aligned to different target sequence"},
	{eNARPEUnalign,(char *)"NP",(char *)"PE alignment not accepted"},
{eNARLociConstrained,(char *)"LC",(char *)"Alignment violated loci base constraints"}};

CAligner::CAligner(void)
{
Init();
}


CAligner::~CAligner(void)
{
Reset(false);
}


int
CAligner::Align(etPMode PMode,			// processing mode
		etFQMethod Quality,				// quality scoring for fastq sequence files
		bool bSOLiD,					// if true then processing in colorspace
		bool bBisulfite,				// if true then process for bisulfite methylation patterning
		etPEproc PEproc,				// paired reads alignment processing mode
		int PairMinLen,					// accept paired end alignments with apparent length of at least this (default = 100)
		int PairMaxLen,					// accept paired end alignments with apparent length of at most this (default = 300)
		bool bPairStrand,				// accept paired ends if on same strand
		bool bPEcircularised,			// experimental - true if processing for PE spaning circularised fragments
		bool bPEInsertLenDist,			// experimental - true if stats file to include PE insert length distributions for each transcript
		eALStrand AlignStrand,			// align on to watson, crick or both strands of target
		int MinChimericLen,				// minimum chimeric length as a percentage (0 to disable, otherwise 50..99) of probe sequence length: negative if chimeric diagnostics to be reported
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
		int MinFlankExacts,				// trim matched reads on 5' and 3' flanks until at least this number of exactly matching bases in flanks
		int PCRPrimerCorrect,			// initially align with MaxSubs+PCRPrimerCorrect subs allowed but then correct substitutions in 5' 12bp until overall sub rate within MaxSubs
		int MaxRptSAMSeqsThres,			// report all SAM chroms or sequences if number of reference chroms <= this limit (defaults to 10000)
		etFMode FMode,					// output format mode
		teSAMFormat SAMFormat,			// if SAM output format then could be SAM or BAM compressed dependent on the file extension used
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
		char *pszSNPFile,				// Output SNPs (CSV or VCF format) to this file (default is to output file name with '.snp' appended)
		char *pszMarkerFile,			// Output markers to this file
		char *pszSNPCentroidFile,		// Output SNP centroids (CSV format) to this file (default is for no centroid processing)
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
int Idx;
int BuffLen = 0;
int BuffOfs = 0;
int SeqIdx;
char szPEInsertDistFile[_MAX_PATH];
char szOutBAIFile[_MAX_PATH];
Init();


m_pszTrackTitle = pszTrackTitle;

m_pszOutFile = pszOutFile;

if(SAMFormat == etSAMFBAM)
	{
	strcpy(szOutBAIFile,pszOutFile);
	strcat(szOutBAIFile,".bai");
	}
else
	szOutBAIFile[0] = '\0';
m_pszOutBAIFile = szOutBAIFile;

m_pszSNPRsltsFile = pszSNPFile;
if(m_FMode != eFMbed && pszSNPFile != NULL && pszSNPFile[0] != '\0')
	{
	// check if file extension is VCF, if so then output SNPs as VCF instead of the default CSV
	int SNPFileNameLen;
	SNPFileNameLen = (int)strlen(m_pszSNPRsltsFile);

	if(!stricmp(".vcf",&pszSNPFile[SNPFileNameLen-4]))
		m_bSNPsVCF = true;
	}
m_pszSNPCentroidFile = pszSNPCentroidFile;
m_NumPE2InputFiles = NumPE2InputFiles;

m_ppszPE2InputFiles = pszPE2InputFiles;
m_pszMultiAlignFile = pszMultiAlignFile;
m_pszNoneAlignFile = pszNoneAlignFile;
m_pszStatsFile = pszStatsFile;
m_pszSfxFile = pszSfxFile;

if(bPEInsertLenDist && m_pszStatsFile != NULL && m_pszStatsFile[0] != '\0')
	{
	strcpy(szPEInsertDistFile,pszStatsFile);
	strcat(szPEInsertDistFile,".peinserts.csv");
	m_pszPEInsertDistFile = szPEInsertDistFile;
	}
else
	{
	bPEInsertLenDist = false;
	m_pszPEInsertDistFile = NULL;
	}

m_pszSitePrefsFile = pszSitePrefsFile;

m_NumPE1InputFiles = NumPE1InputFiles;
m_ppszPE1InputFiles = pszPE1InputFiles;
m_bFiltPriorityRegions = bFiltPriorityRegions;

m_bPEcircularised = bPEcircularised;		// experimental - true if processing for PE spaning circularised fragments
m_bPEInsertLenDist = bPEInsertLenDist;		// true if stats file to include PE insert length distributions for each transcript

m_PMode = PMode;
m_FMode = FMode;
m_SAMFormat = SAMFormat;
m_PEproc = PEproc;
m_QMethod = Quality;
m_NumThreads = NumThreads;
m_bBisulfite = bBisulfite;
m_MaxMLmatches = MaxMLmatches;
m_bClampMaxMLmatches = bClampMaxMLmatches;
m_bLocateBestMatches = bLocateBestMatches; 
m_MLMode = MLMode;
m_MaxNs = MaxNs;
m_MaxSubs = MaxSubs;
if(PCRPrimerCorrect > 0)
	m_InitalAlignSubs = min(MaxSubs + PCRPrimerCorrect,cMaxAllowedSubs);
else
	m_InitalAlignSubs = MaxSubs;

m_AlignStrand = AlignStrand;
m_MinChimericLen = abs(MinChimericLen);
m_bReportChimerics = MinChimericLen >= 0 ? false : true;
m_microInDelLen = microInDelLen;
m_SpliceJunctLen = SpliceJunctLen;
m_MinSNPreads = MinSNPreads;
m_SNPNonRefPcnt = SNPNonRefPcnt/100.0;

m_QValue = QValue;
m_Marker5Len = MarkerLen/2;
m_Marker3Len = MarkerLen - 1 - m_Marker5Len;
m_MarkerPolyThres = MarkerPolyThres;
m_pszMarkerFile = pszMarkerFile;

m_Trim5 = Trim5;						// trim this number of bases from 5' end of reads when loading the reads
m_Trim3 = Trim3;						// trim this number of bases from 3' end of reads when loading the reads

m_SitePrefsOfs = SitePrefsOfs;			// offset read start sites when processing  octamer preferencing, range -100..100


m_NumIncludeChroms = NumIncludeChroms;
m_NumExcludeChroms = NumExcludeChroms;

m_MaxRptSAMSeqsThres = MaxRptSAMSeqsThres;

if(CreateMutexes()!=eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Failed to create thread synchronisation mutexes");
	Reset(false);
	return(cBSFSyncObjErr);
	}

m_mtqsort.SetMaxThreads(NumThreads);

// load contaminants if user has specified a contaminant sequence file
if(pszContamFile != NULL && pszContamFile[0] != '\0')
	{
	if((m_pContaminants = new CContaminants)==NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CContaminants");
		Reset(false);
		return(eBSFerrObj);
		}
	if((Rslt = m_pContaminants->LoadContaminantsFile(pszContamFile)) < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to load contaminate sequences file '%s'",pszContamFile);
		Reset(false);
		return(Rslt);
		}
	if(Rslt = 0)
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"No contaminant sequences loaded from '%s'",pszContamFile);
	}

// compile include/exclude chromosome regexpr if user has specified alignments to be filtered by chrom
if(NumIncludeChroms > 0 || m_NumExcludeChroms > 0)
	if((Rslt = CompileChromRegExprs(NumIncludeChroms,ppszIncludeChroms,NumExcludeChroms,ppszExcludeChroms)) != eBSFSuccess)
		{
		Reset(false);
		return(Rslt);
		}

// if specified then load high priority regions from BED file
m_pPriorityRegionBED = NULL;
if(pszPriorityRegionFile != NULL && pszPriorityRegionFile[0] != '\0')
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Loading high priority regions BED file '%s'", pszPriorityRegionFile);
	if((m_pPriorityRegionBED = new CBEDfile()) == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CBEDfile");
		Reset(false);
		return(eBSFerrObj);
		}
	if((Rslt=m_pPriorityRegionBED->Open(pszPriorityRegionFile))!=eBSFSuccess)
		{
		while(m_pSfxArray->NumErrMsgs())
			gDiagnostics.DiagOut(eDLFatal,gszProcName,m_pPriorityRegionBED->GetErrMsg());
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open high priority regions BED file '%s'",pszPriorityRegionFile);
		Reset(false);
		return(Rslt);
		}
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"High priority regions BED file '%s' loaded", pszPriorityRegionFile);
	}


// open bioseq file containing suffix array for targeted assembly to align reads against
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Loading suffix array file '%s'", pszSfxFile);
if((m_pSfxArray = new CSfxArrayV3()) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CSfxArrayV3");
	Reset(false);
	return(eBSFerrObj);
	}
if((Rslt=m_pSfxArray->Open(pszSfxFile,false,bBisulfite,bSOLiD))!=eBSFSuccess)
	{
	while(m_pSfxArray->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,m_pSfxArray->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open input bioseq suffix array file '%s'",pszSfxFile);
	Reset(false);
	return(Rslt);
	}

// report to user some sfx array metadata as conformation the targeted assembly is correct
strcpy(m_szTargSpecies,m_pSfxArray->GetDatasetName());
m_bIsSOLiD = m_pSfxArray->IsSOLiD();
tsSfxHeaderV3 SfxHeader;
m_pSfxArray->GetSfxHeader(&SfxHeader);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Genome Assembly Name: '%s' Descr: '%s' Title: '%s' Version: %d",
					 m_szTargSpecies,SfxHeader.szDescription,SfxHeader.szTitle,SfxHeader.Version);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Assembly has blocks: %d, max block size: %llu",SfxHeader.NumSfxBlocks,SfxHeader.SfxBlockSize);

// if user has specified that there are constraints on loci bases then need to load the loci constraints file
if(pszLociConstraintsFile != NULL && pszLociConstraintsFile[0] != '\0')
	{
	if((Rslt=LoadLociConstraints(pszLociConstraintsFile)) < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to load loci constraints");
		Reset(false);
		return(Rslt);
		}
	}


// restrict the max core iterations according to the requested sensitivity
int MaxIter;
switch(PMode) {
	case ePMdefault:			// default processing mode
		MaxIter = cDfltSensCoreIters;
		break;
	case ePMMoreSens:			// more sensitive - slower
		MaxIter = cMoreSensCoreIters;
		break;
	case ePMUltraSens:			// ultra sensitive - much slower
		MaxIter = cUltraSensCoreIters;
		break;
	case ePMLessSens:			// less sensitive - quicker
	default:
		MaxIter = cMinSensCoreIters;
	}
m_pSfxArray->SetMaxIter(MaxIter);

// reads are loaded asynchronously to the alignment processing
if((Rslt=InitiateLoadingReads()) < eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Failed to load reads");
	Reset(false);
	return(Rslt);
	}

// if user wants to utilise reads aligning to multiple loci then need to make an initial alloc for holding these..
if(m_MLMode > eMLrand && m_MLMode != eMLall)
	{
	size_t memreq = (size_t)sizeof(tsReadHit) * cAllocMultihits;

#ifdef _WIN32
	m_pMultiHits = (tsReadHit *) malloc(memreq);	// initial and perhaps the only allocation

	if(m_pMultiHits == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadReads: Memory allocation of %lld bytes - %s",(INT64)memreq,strerror(errno));
		Reset(false);
		return(eBSFerrMem);
		}
#else
	// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
	m_pMultiHits = (tsReadHit *)mmap(NULL,memreq, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
	if(m_pMultiHits == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadReads: Memory allocation of %lld bytes through mmap()  failed - %s",(INT64)memreq,strerror(errno));
		m_pMultiHits = NULL;
		Reset(false);
		return(eBSFerrMem);
		}
#endif
	m_AllocdMultiHits = cAllocMultihits;
	m_AllocdMultiHitsMem = memreq;
	m_NumMultiHits = 0;
	m_NumUniqueMultiHits = 0;
	m_NumProvMultiAligned = 0;
	}

// if outputing multiloci all then need to allocate memory for these
if(m_MLMode >= eMLall || m_FMode == eFMsamAll)
	{
	size_t memreq = (size_t)(sizeof(tsReadHit) + 150) * 5 * cAllocMultihits;	// read sizes not known yet so assume 100bp reads plus long descriptors plus many multiloci - realloc'd as may be required

#ifdef _WIN32
	m_pMultiAll = (tsReadHit *) malloc(memreq);	// initial and perhaps the only allocation

	if(m_pMultiAll == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadReads: Memory allocation of %lld bytes for multiAll - %s",(INT64)memreq,strerror(errno));
		Reset(false);
		return(eBSFerrMem);
		}
#else
	// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
	m_pMultiAll = (tsReadHit *)mmap(NULL,memreq, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
	if(m_pMultiAll == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadReads: Memory allocation of %lld bytes for multiAll through mmap()  failed - %s",(INT64)memreq,strerror(errno));
		m_pMultiAll = NULL;
		Reset(false);
		return(eBSFerrMem);
		}
#endif
	m_AllocMultiAllMem = memreq;
	m_NumMultiAll = 0;
	m_NxtMultiAllOfs = 0;
	}

// reasonably confident that there will be some results to report so create/trunc all required result files
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Creating/truncating result files..");
if((Rslt=CreateOrTruncResultFiles())!=eBSFSuccess)
	{
	// problems.. need to ensure all background threads (at this stage should only be the reads loading or sfx loading thread) are cleanly terminated
	m_TermBackgoundThreads = 1;	// need to immediately self-terminate?
	m_pSfxArray->Reset(false);
#ifdef _WIN32
	if(m_hThreadLoadReads != NULL)
		{
		while(WAIT_TIMEOUT == WaitForSingleObject(m_hThreadLoadReads, 5000))
			{
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"Progress: waiting for reads load thread to terminate");
			}
		CloseHandle(m_hThreadLoadReads);
		}
#else
	if(m_ThreadLoadReadsID != 0)
		{
		struct timespec ts;
		int JoinRlt;
		clock_gettime(CLOCK_REALTIME, &ts);
		ts.tv_sec += 5;
		while((JoinRlt = pthread_timedjoin_np(m_ThreadLoadReadsID, NULL, &ts)) != 0)
			{
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"Progress: waiting for reads load thread to terminate");
			ts.tv_sec += 60;
			}
		}
#endif
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Progress: reads load thread terminated");
	Reset(false);
	return(Rslt);
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Creating/truncating result files completed");
m_CurReadsSortMode = eRSMReadID;			// reads were loaded and assigned ascending read identifiers so that is their initial sort order

// locate all read matches
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Aligning in %s...",bSOLiD ? "colorspace" : "basespace");
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Aligning for %s cored matches...",bBisulfite ? "bisulfite" : "normal");

// heavy lifting now starts!
Rslt = LocateCoredApprox(MinEditDist,m_InitalAlignSubs);

if(Rslt < eBSFSuccess)
	{
	Reset(false);
	return(Rslt);
	}



if(m_MinChimericLen > 0) // any aligned reads accepted as being chimeric are flank trimmed and subsequently processed as if full length aligned
	{
	char szChimericsFile[_MAX_PATH];
	if(m_bReportChimerics)
		{
		strcpy(szChimericsFile,m_pszOutFile);
		strcat(szChimericsFile,".chimericseqs.csv");
		}
	else
		szChimericsFile[0] = '\0';
	TrimChimeric(szChimericsFile);
	}

// if autodetermining max subs that were allowed from actual reads then let user know what the average read length was
size_t TotReadsLen = 0;
tsReadHit *pReadHit;
int AvReadsLen;
int MinReadsLen = -1;
int MaxReadsLen = 0;
pReadHit = m_pReadHits;
for(Idx = 0; Idx < (int)m_NumReadsLoaded; Idx++)
	{
	TotReadsLen += pReadHit->ReadLen;
	if(MinReadsLen > pReadHit->ReadLen || MinReadsLen == -1)
		MinReadsLen = pReadHit->ReadLen;
	if(MaxReadsLen < pReadHit->ReadLen)
		MaxReadsLen = pReadHit->ReadLen;
	pReadHit = (tsReadHit *)((UINT8 *)pReadHit + sizeof(tsReadHit) + pReadHit->ReadLen + pReadHit->DescrLen);
	}
AvReadsLen = (int)(TotReadsLen/m_NumReadsLoaded);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Average length of all reads was: %d (min: %d, max: %d)",AvReadsLen,MinReadsLen,MaxReadsLen);
m_MaxReadsLen = MaxReadsLen;
m_MinReadsLen = MinReadsLen;
m_AvReadsLen = AvReadsLen;
if(gProcessingID > 0)
	{
	gSQLiteSummaries.AddResult(gProcessingID,(char *)"ReadLen",ePTInt32,sizeof(AvReadsLen),"MeanLen",&AvReadsLen);
	gSQLiteSummaries.AddResult(gProcessingID,(char *)"ReadLen",ePTInt32,sizeof(MinReadsLen),"MinLen",&MinReadsLen);
	gSQLiteSummaries.AddResult(gProcessingID,(char *)"ReadLen",ePTInt32,sizeof(MaxReadsLen),"MaxLen",&MaxReadsLen);
	}

if(m_InitalAlignSubs != 0)
	{
	int MeanSubs;
	int MinSubs;
	int MaxSubs;
	MeanSubs = max(1,(AvReadsLen * m_InitalAlignSubs)/100);
	MinSubs = max(1,(MinReadsLen * m_InitalAlignSubs)/100);
	MaxSubs = max(1,(MaxReadsLen * m_InitalAlignSubs)/100);

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Typical allowed aligner induced substitutions was: %d (min: %d, max: %d)", MeanSubs, MinSubs, MaxSubs);
	if(gProcessingID > 0)
		{
		gSQLiteSummaries.AddResult(gProcessingID,(char *)"AllowedSubs",ePTInt32,sizeof(MeanSubs),"MeanSubs",&MeanSubs);
		gSQLiteSummaries.AddResult(gProcessingID,(char *)"AllowedSubs",ePTInt32,sizeof(MinSubs),"MinSubs",&MinSubs);
		gSQLiteSummaries.AddResult(gProcessingID,(char *)"AllowedSubs",ePTInt32,sizeof(MaxSubs),"MaxSubs",&MaxSubs);
		}
	}
else
	if(gProcessingID > 0)
		gSQLiteSummaries.AddResult(gProcessingID,(char *)"AllowedSubs",ePTInt32,sizeof(MaxSubs),"MeanSubs",&m_InitalAlignSubs);
 
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Provisionally accepted %d aligned reads (%d uniquely, %d aligning to multiloci) aligning to a total of %d loci", m_TotAcceptedAsAligned,m_TotAcceptedAsUniqueAligned,m_TotAcceptedAsMultiAligned,m_TotLociAligned);

m_OrigNumReadsLoaded = m_NumReadsLoaded;		// make a copy of actual number of reads loaded as if m_MLMode >= eMLall then will be overwritten with number of multialigned loci 
if(m_MLMode >= eMLall)		// a little involved as need to reuse m_pReadHits ptrs so sorting and reporting of multi-SAM hits will be same as if normal processing...
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Treating accepted %d multialigned reads as uniquely aligned %d source reads in subsequent processing",m_TotAcceptedAsMultiAligned,m_TotLociAligned - m_TotAcceptedAsUniqueAligned);
	if(m_pReadHits != NULL)
		{
#ifdef _WIN32
		free(m_pReadHits);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
		if(m_pReadHits != MAP_FAILED)
			munmap(m_pReadHits,m_AllocdReadHitsMem);
#endif
		m_pReadHits = NULL;
		}
	m_pReadHits = m_pMultiAll;
	m_pMultiAll = NULL;
	m_AllocdReadHitsMem = m_AllocMultiAllMem;
	m_UsedReadHitsMem = m_NxtMultiAllOfs;
	m_NumReadsLoaded = m_NumMultiAll;
	m_FinalReadID = m_NumMultiAll;
	m_AllocMultiAllMem = 0;
	m_NxtMultiAllOfs = 0;
	m_NumMultiAll = 0;
	if(m_ppReadHitsIdx != NULL)
		{
		delete m_ppReadHitsIdx;
		m_ppReadHitsIdx = NULL;
		m_AllocdReadHitsIdx = 0;
		}
	m_ppReadHitsIdx = 0;
	}

if(PEproc != ePEdefault)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Paired end association and partner alignment processing started..");
	if((Rslt=ProcessPairedEnds(PEproc,MinEditDist,PairMinLen,PairMaxLen,bPairStrand,m_InitalAlignSubs)) < eBSFSuccess)
		{
		Reset(false);
		return(Rslt);
		}
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Paired end association and partner alignment processing completed..");
	}

// try and assign multimatch read loci?
if(PEproc == ePEdefault && m_MLMode > eMLrand && m_MLMode != eMLall)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Multialignment processing started..");
	if((Rslt = AssignMultiMatches()) < eBSFSuccess)
		{
		Reset(false);
		return(Rslt);
		}
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Multialignment processing completed");
	}

IdentifyConstraintViolations(PEproc != ePEdefault);

// if requested then attempt to reduce the number of  PCR differential amplification artefacts (reads stacking to same loci)
if(PEproc == ePEdefault && PCRartefactWinLen >= 0)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing to reduce PCR differential amplification artefacts processing started..");
	if((Rslt=ReducePCRduplicates(PCRartefactWinLen)) < eBSFSuccess)
		{
		Reset(false);
		return(Rslt);
		}
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"PCR differential amplification artefacts processing completed");
	}

if(PCRPrimerCorrect > 0) 
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"PCR 5' Primer correction processing started..");
	if((Rslt=PCR5PrimerCorrect(m_MaxSubs)) < eBSFSuccess)
		{
		Reset(false);
		return(Rslt);
		}
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"PCR 5' Primer correction processing completed");
	}

if(MinFlankExacts > 0) // autotrim aligned reads flanks?
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Autotrim aligned read flank processing started..");
	if((Rslt=AutoTrimFlanks(MinFlankExacts)) < eBSFSuccess)
		{
		Reset(false);
		return(Rslt);
		}
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Autotrim aligned read flank processing completed");
	}

// if splice junctions being processed then check for orphans and remove these
if(PEproc == ePEdefault && SpliceJunctLen > 0)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Removal of orphan splice junction processing started..");
	if((Rslt=RemoveOrphanSpliceJuncts(SpliceJunctLen)) < eBSFSuccess)
		{
		Reset(false);
		return(Rslt);
		}
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Removal of orphan splice junction processing completed");
	}

// if processing for microInDels then need to check for orphans and remove these
if(PEproc == ePEdefault && microInDelLen > 0)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Removal of orphan microInDels processing started..");
	if((Rslt=RemoveOrphanMicroInDels(microInDelLen)) < eBSFSuccess)
		{
		Reset(false);
		return(Rslt);
		}
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Removal of orphan microInDels processing completed");
	}

BuffLen = 0;
BuffOfs = 0;
SeqIdx;

// now apply any chromosome filtering that user has specified
if(NumExcludeChroms || NumIncludeChroms)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Filtering aligned reads by chromosome started..");
	if((Rslt=FiltByChroms())<eBSFSuccess)
		{
		Reset(false);
		return(Rslt);
		}
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Filtering aligned reads by chromosome completed");
	}

// apply priority region filtering if requested
if(m_pPriorityRegionBED != NULL && m_bFiltPriorityRegions)
	FiltByPriorityRegions();

// user interested in the nonaligned?
if(m_hNoneAlignFile != -1 || m_gzNoneAlignFile != NULL)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Reporting of non-aligned reads started..");
	if((Rslt=ReportNoneAligned())<eBSFSuccess)
		{
		Reset(false);
		return(Rslt);
		}
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Reporting of non-aligned reads completed");
	}

// user interested in the multialigned?
// these only include those reads which would otherwise have been accepted but aligned to multiple loci
if(PEproc == ePEdefault && (m_hMultiAlignFile != -1 || m_gzMultiAlignFile != NULL))
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Reporting of multialigned reads started..");
	if((Rslt=ReportMultiAlign()) < eBSFSuccess)
		{
		Reset(false);
		return(Rslt);
		}
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Reporting of multialigned reads completed");
	}

// all processing to accept alignments completed, can now report basic stats
if((Rslt=ReportAlignStats()) < eBSFSuccess)
	{
	Reset(false);
	return(Rslt);
	}

if(m_NARAccepted && m_hSitePrefsFile != -1)
	ProcessSiteProbabilites(m_SitePrefsOfs);

// now time to write out the read hits
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Reporting of aligned result set started...");
if(FMode >= eFMsam)
	Rslt = WriteBAMReadHits(FMode,SAMFormat,PEproc == ePEdefault ? false : true,6);	// default to compression level 6
else
	Rslt = WriteReadHits(PEproc == ePEdefault ? false : true);

if(Rslt < eBSFSuccess)
	{
	Reset(false);
	return(Rslt);
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Reporting of aligned result set completed");

if(m_bPEInsertLenDist && m_NARAccepted)
	ReportPEInsertLenDist();

if(Rslt >= eBSFSuccess && m_NARAccepted && m_hStatsFile > 0 && m_MaxAlignLen > 0)
	Rslt = WriteBasicCountStats();

if(Rslt >= eBSFSuccess  && m_NARAccepted && m_hSitePrefsFile > 0)
	Rslt = WriteSitePrefs();

m_TotNumSNPs = 0;
if(Rslt >= eBSFSuccess && m_NARAccepted && MinSNPreads > 0 && m_hSNPfile != -1)
	{
	bool bMarkers;
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing for SNPs and writing out SNPs to file '%s",m_pszSNPRsltsFile);
	if(m_hMarkerFile != -1)
		{
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing for Markers and writing out marker sequences to file '%s",m_pszMarkerFile);
		bMarkers = true;
		}
	else
		bMarkers = false;
	Rslt = ProcessSNPs();			// track title if output format is to be UCSC BED, will have '_SNPs' appended
	if(Rslt >= eBSFSuccess)
		{
		if(bMarkers)
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"Marker processing completed with %d marker sequences writtten to file '%s",m_MarkerID,m_pszMarkerFile);
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"SNP processing completed with %d putative SNPs discovered",m_TotNumSNPs);
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"There are %lld aligned loci bases which are covered by %lld read bases with mean coverage of %1.2f",m_LociBasesCovered,m_LociBasesCoverage,m_LociBasesCoverage/(double)m_LociBasesCovered);
		}
	}
if(gProcessingID != 0)
	gSQLiteSummaries.AddResult(gProcessingID,(char *)"SNPs",ePTInt32,sizeof(AvReadsLen),"Cnt",&m_TotNumSNPs);
Reset(Rslt >= eBSFSuccess ? true : false);
return(Rslt);
}

void
CAligner::Init(void)
{
m_hInFile = -1;
m_hOutFile = -1;

m_hBAIFile = -1;

m_hIndOutFile = -1;
m_hJctOutFile = -1;
m_hStatsFile = -1;
m_hInsertLensFile = -1;
m_hSitePrefsFile = -1;
m_hNoneAlignFile = -1;
m_hMultiAlignFile = -1;
m_hSNPfile = -1;
m_hDiSNPfile = -1;
m_hTriSNPfile = -1;
m_hMarkerFile = -1;	
m_hSNPCentsfile = -1;

m_gzOutFile = NULL;
m_gzIndOutFile = NULL;
m_gzJctOutFile = NULL;
m_gzNoneAlignFile = NULL;
m_gzMultiAlignFile = NULL;
m_gzSNPfile = NULL;
m_gzSNPCentsfile = NULL;

m_bgzOutFile = false;
m_bgzNoneAlignFile = false;
m_bgzMultiAlignFile = false;

m_pReadHits = NULL;
m_ppReadHitsIdx = NULL;
m_ppReadHitsIdx = NULL;
m_pMultiHits = NULL;
m_pMultiAll = NULL;
m_pLociPValues = NULL;
m_pSfxArray = NULL;
m_pPriorityRegionBED = NULL;
m_pAllocsIdentNodes = NULL;
m_pAllocsMultiHitLoci = NULL;
m_pAllocsMultiHitBuff = NULL;
m_pChromSNPs = NULL;
m_pszLineBuff = NULL;
m_pLenDist = NULL;
m_pSNPCentroids = NULL;
m_pConstraintLoci = NULL; 
m_pContaminants = NULL;
m_NumConstraintLoci = 0;
m_NumConstrainedChroms = 0;
m_ConstrainedChromIDs[0] = 0;
m_bSNPsVCF = false;
m_LociBasesCovered = 0;
m_LociBasesCoverage = 0;
m_PrevSizeOf = 0;
m_NumMultiAll = 0;
m_NxtMultiAllOfs = 0;
m_AllocMultiAllMem = 0;
m_szLineBuffIdx = 0;
m_MaxMLmatches = 0;
m_bClampMaxMLmatches = false;
m_bLocateBestMatches = false;
m_TotAllocdIdentNodes = 0;
m_PerThreadAllocdIdentNodes = 0;
m_AllocdReadHitsMem = 0;
m_UsedReadHitsMem = 0;
m_NumReadsLoaded = 0;
m_UsedReadHitsMem = 0;
m_FinalReadID = 0;
m_NumDescrReads = 0;
m_FinalReadID = 0;
m_AllocdReadHitsIdx = 0;
m_AllocdMultiHits = 0;
m_AllocdMultiHitsMem = 0;
m_NumMultiHits = 0;
m_NumUniqueMultiHits = 0;
m_NumProvMultiAligned = 0;
m_MaxAlignLen = 0;
m_MaxNs = 0;
m_MaxRptSAMSeqsThres = 0;
m_bPEcircularised = false;
m_bPEInsertLenDist = false;	
m_bIsSOLiD = false;
m_bBisulfite = false;
m_AlignStrand = eALSnone;
m_MinChimericLen = 0;
m_microInDelLen = 0;
m_SpliceJunctLen = 0;
m_AllocLociPValuesMem = 0;
m_NumLociPValues = 0;
m_QValue = 0.0;
m_MinSNPreads = 0;
m_SNPNonRefPcnt = 0.0; 
m_MaxDiSNPSep = cDfltMaxDiSNPSep;
m_MarkerID = 0;	
m_Marker5Len = 0;
m_Marker3Len = 0;
m_MarkerPolyThres = 0;
m_NumReadsProc = 0;
m_NxtReadProcOfs = 0;
m_ElimPlusTrimed = 0;
m_ElimMinusTrimed = 0;
m_PEproc = ePEdefault;
m_TotNonAligned = 0;
m_NumSloughedNs = 0;
m_TotAcceptedAsAligned = 0;
m_TotLociAligned = 0;
m_TotAcceptedAsUniqueAligned = 0;
m_TotAcceptedAsMultiAligned = 0;
m_TotNotAcceptedDelta = 0;
m_TotAcceptedHitInsts = 0;
m_SitePrefsOfs = cDfltRelSiteStartOfs;
m_pszOutFile = NULL;
m_szIndRsltsFile[0] = '\0';
m_szJctRsltsFile[0] = '\0';
m_szDiSNPFile[0] = '\0';
m_pszSNPRsltsFile = NULL;
m_pszSNPCentroidFile = NULL;
m_pszSitePrefsFile = NULL;
m_MaxReadsLen = 0;
m_MinReadsLen = 0;
m_AvReadsLen = 0;
m_NARAccepted = 0;
m_TermBackgoundThreads = 0;
m_CurClusterFrom = 0;
m_bFiltPriorityRegions = false;
m_SAMFormat = etSAMFformat;
m_CurReadsSortMode = eRSMunsorted;
m_ThreadCoredApproxRslt = 0;

#ifdef _WIN32
m_hThreadLoadReads = NULL;
#endif
memset(m_AlignQSubDist,0,sizeof(m_AlignQSubDist));
memset(m_AlignMSubDist,0,sizeof(m_AlignMSubDist));
memset(m_MultiHitDist,0,sizeof(m_MultiHitDist));
memset(&m_FileHdr,0,sizeof(m_FileHdr));
m_bMutexesCreated = false;
}

void
CAligner::Reset(bool bSync)			// if bSync true then fsync before closing output file handles
{
m_TermBackgoundThreads = 0x01;	// need to require any background threads to self-terminate
if(m_hInFile != -1)
	{
	close(m_hInFile);
	m_hInFile = -1;
	}
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

if(m_hBAIFile != -1)
	{
	if(bSync)
#ifdef _WIN32
		_commit(m_hBAIFile);
#else
		fsync(m_hBAIFile);
#endif
	close(m_hBAIFile);
	m_hBAIFile = -1;
	}

if(m_hJctOutFile != -1)
	{
	if(bSync)
#ifdef _WIN32
		_commit(m_hJctOutFile);
#else
		fsync(m_hJctOutFile);
#endif
	close(m_hJctOutFile);
	m_hJctOutFile = -1;
	}
if(m_hIndOutFile != -1)
	{
	if(bSync)
#ifdef _WIN32
		_commit(m_hIndOutFile);
#else
		fsync(m_hIndOutFile);
#endif
	close(m_hIndOutFile);
	m_hIndOutFile = -1;
	}

if(m_hStatsFile != -1)
	{
	if(bSync)
#ifdef _WIN32
		_commit(m_hStatsFile);
#else
		fsync(m_hStatsFile);
#endif
	close(m_hStatsFile);
	m_hStatsFile = -1;
	}

if(m_hInsertLensFile != -1)
	{
	if(bSync)
#ifdef _WIN32
		_commit(m_hInsertLensFile);
#else
		fsync(m_hInsertLensFile);
#endif
	close(m_hInsertLensFile);
	m_hInsertLensFile = -1;
	}

if(m_hNoneAlignFile != -1)
	{
	if(bSync)
#ifdef _WIN32
		_commit(m_hNoneAlignFile);
#else
		fsync(m_hNoneAlignFile);
#endif
	close(m_hNoneAlignFile);
	m_hNoneAlignFile = -1;
	}
if(m_hMultiAlignFile != -1)
	{
	if(bSync)
#ifdef _WIN32
		_commit(m_hMultiAlignFile);
#else
		fsync(m_hMultiAlignFile);
#endif
	close(m_hMultiAlignFile);
	m_hMultiAlignFile = -1;
	}
if(m_hSitePrefsFile != -1)
	{
	if(bSync)
#ifdef _WIN32
		_commit(m_hSitePrefsFile);
#else
		fsync(m_hSitePrefsFile);
#endif
	close(m_hSitePrefsFile);
	m_hSitePrefsFile = -1;
	}

if(m_hSNPfile != -1)
	{
	if(bSync)
#ifdef _WIN32
		_commit(m_hSNPfile);
#else
		fsync(m_hSNPfile);
#endif
	close(m_hSNPfile);
	m_hSNPfile = -1;
	}

if(m_hDiSNPfile != -1)
	{
	if(bSync)
#ifdef _WIN32
		_commit(m_hDiSNPfile);
#else
		fsync(m_hDiSNPfile);
#endif
	close(m_hDiSNPfile);
	m_hDiSNPfile = -1;
	}
if(m_hTriSNPfile != -1)
	{
	if(bSync)
#ifdef _WIN32
		_commit(m_hTriSNPfile);
#else
		fsync(m_hTriSNPfile);
#endif
	close(m_hTriSNPfile);
	m_hTriSNPfile = -1;
	}

if(m_hMarkerFile != -1)
	{
	if(bSync)
#ifdef _WIN32
		_commit(m_hMarkerFile);
#else
		fsync(m_hMarkerFile);
#endif
	close(m_hMarkerFile);
	m_hMarkerFile = -1;
	}

if(m_hSNPCentsfile != -1)
	{
	if(bSync)
#ifdef _WIN32
		_commit(m_hSNPCentsfile);
#else
		fsync(m_hSNPCentsfile);
#endif
	close(m_hSNPCentsfile);
	m_hSNPCentsfile = -1;
	}

if(m_pSNPCentroids != NULL)
	{
	delete(m_pSNPCentroids);
	m_pSNPCentroids = NULL;
	}


if(m_pConstraintLoci != NULL)
	{
	delete(m_pConstraintLoci);
	m_pConstraintLoci = NULL;
	}

if(m_gzOutFile != NULL)
	{
	gzclose(m_gzOutFile);
	m_gzOutFile = NULL;
	}

if(m_gzIndOutFile != NULL)
	{
	gzclose(m_gzIndOutFile);
	m_gzIndOutFile = NULL;
	}

if(m_gzJctOutFile != NULL)
	{
	gzclose(m_gzJctOutFile);
	m_gzJctOutFile = NULL;
	}

if(m_gzNoneAlignFile != NULL)
	{
	gzclose(m_gzNoneAlignFile);
	m_gzNoneAlignFile = NULL;
	}

if(m_gzMultiAlignFile != NULL)
	{
	gzclose(m_gzMultiAlignFile);
	m_gzMultiAlignFile = NULL;
	}

if(m_gzSNPfile != NULL)
	{
	gzclose(m_gzSNPfile);
	m_gzSNPfile = NULL;
	}

if(m_gzSNPCentsfile != NULL)
	{
	gzclose(m_gzSNPCentsfile);
	m_gzSNPCentsfile = NULL;
	}

if(m_pszLineBuff != NULL)
	{
	delete m_pszLineBuff;
	m_pszLineBuff = NULL;
	}
if(m_pAllocsIdentNodes != NULL)
	{
	delete m_pAllocsIdentNodes;
	m_pAllocsIdentNodes = NULL;
	}

if(m_pAllocsMultiHitBuff != NULL)
	{
	delete m_pAllocsMultiHitBuff;
	m_pAllocsMultiHitBuff = NULL;
	}

if(m_pAllocsMultiHitLoci != NULL)
	{
	delete m_pAllocsMultiHitLoci;
	m_pAllocsMultiHitLoci = NULL;
	}
if(m_pSfxArray != NULL)
	{
	delete m_pSfxArray;
	m_pSfxArray = NULL;
	}
if(m_pPriorityRegionBED != NULL)
	{
	delete m_pPriorityRegionBED;
	m_pPriorityRegionBED = NULL;
	}

if(m_pReadHits != NULL)
	{
#ifdef _WIN32
	free(m_pReadHits);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pReadHits != MAP_FAILED)
		munmap(m_pReadHits,m_AllocdReadHitsMem);
#endif
	m_pReadHits = NULL;
	}

if(m_pMultiAll != NULL)
	{
#ifdef _WIN32
	free(m_pMultiAll);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pMultiAll != MAP_FAILED)
		munmap(m_pMultiAll,m_AllocMultiAllMem);
#endif
	m_pMultiAll = NULL;
	}

if(m_ppReadHitsIdx != NULL)
	{
	delete m_ppReadHitsIdx;
	m_ppReadHitsIdx = NULL;
	}
if(m_pMultiHits != NULL)
	{
#ifdef _WIN32
	free(m_pMultiHits);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pMultiHits != MAP_FAILED)
		munmap(m_pMultiHits,m_AllocdMultiHitsMem);
#endif
	m_pMultiHits = NULL;
	}

if(m_pLociPValues != NULL)
	{
#ifdef _WIN32
	free(m_pLociPValues);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pLociPValues != MAP_FAILED)
		munmap(m_pLociPValues,m_AllocLociPValuesMem);
#endif
	m_pLociPValues = NULL;
	}

if(m_pChromSNPs != NULL)
	{
	delete m_pChromSNPs;
	m_pChromSNPs = NULL;
	}

if(m_pLenDist != NULL)
	{
	delete m_pLenDist;
	m_pLenDist = NULL;
	}

if(m_pContaminants != NULL)
	{
	delete m_pContaminants;
	m_pContaminants = NULL;
	}

DeleteMutexes();

Init();
m_TermBackgoundThreads = 0x0;	// can startup any background thread processing
}


// Load alignment loci base constraints from CSV file
// Expected file format is -
// <chrom>,<startloci>,<endloci>,<baselist>
// Whereby -
// chrom  must the name of a targeted chrom/sequence
// startloci is the start loci on the chrom
// endloci is the end loci on the chrom
// baselist names one or more bases which a read must align with for that aligned read to be accepted
//          if the baselist is '.' then no substitutions allowed at that loci
// examples -
// Chr1A,100,100,RA       reads aligning over Chr1A.100 will only be accepted if the read base is A or the target base
// Chr1A,100,100,C        reads aligning over Chr1A.100 will only be accepted if the read base is C
// Chr1A,100,100,CT       reads aligning over Chr1A.100 will only be accepted if the read base is C or T
// Chr1A,100,107,R        reads aligning over Chr1A.100 to Chr1A.107 will only be accepted if the read bases matches the target bases
// 

teBSFrsltCodes
CAligner::LoadLociConstraints(char *pszLociConstraints)	// load loci constraints from file
{
int Rslt;
int NumLines;
int NumFields;
int CSVLineNum;

int Idx;

CCSVFile InFile;

char *pszTargChrom;
char szPrevTargChrom[cMaxGeneNameLen];

int ChromID;
int StartLoci;
int EndLoci;
UINT8 Constraint;

char *pszBases;
char *pBase;
char Base;

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Loading loci base constraints from CSV file '%s' ...",pszLociConstraints);

if((Rslt=InFile.Open(pszLociConstraints)) !=eBSFSuccess)
	{
	while(InFile.NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,InFile.GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open '%s' for processing",pszLociConstraints);
	return(eBSFerrOpnFile);
	}

if(m_pConstraintLoci == NULL)
	m_pConstraintLoci = new tsConstraintLoci [cMaxConstrainedLoci];


NumLines = 0;
NumFields = 0;
szPrevTargChrom[0] = 0;
ChromID = 0;
while((Rslt=InFile.NextLine()) > 0)	// onto next line containing fields
	{
	NumLines += 1;
	NumFields = InFile.GetCurFields();
	if(NumFields < 4)
		{
		CSVLineNum = InFile.GetLineNumber();
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expected at least 4 fields at line %d in '%s', GetCurFields() returned '%d'",CSVLineNum,pszLociConstraints,NumFields);
		InFile.Close();
		return(eBSFerrFieldCnt);
		}

	if(NumLines == 1)		// slough 1st line if it is a title line, assuming title line contains only text values ... 
		{
		if(InFile.IsLikelyHeaderLine())
			continue;
		}

	InFile.GetText(1,&pszTargChrom);
	InFile.GetInt(2,&StartLoci);
	InFile.GetInt(3,&EndLoci);
	InFile.GetText(4,&pszBases);

	// get chrom identifier
	if(ChromID == 0 || szPrevTargChrom[0] == 0 || stricmp(pszTargChrom,szPrevTargChrom))
		{
		strncpy(szPrevTargChrom,pszTargChrom,sizeof(szPrevTargChrom)-1);
		szPrevTargChrom[sizeof(szPrevTargChrom)-1] = '\0';
		if((ChromID = m_pSfxArray->GetIdent(szPrevTargChrom)) <= 0)
			{
			CSVLineNum = InFile.GetLineNumber();
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to find matching indexed identifier for '%s' at line %d in '%s'",szPrevTargChrom,CSVLineNum,pszLociConstraints);
			InFile.Close();
			return(eBSFerrFieldCnt);
			}
		
		}

	if(StartLoci < 0 || StartLoci > EndLoci)
		{
		CSVLineNum = InFile.GetLineNumber();
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Start loci must be >= 0 and <= end loci for '%s' at line %d in '%s'",szPrevTargChrom,CSVLineNum,pszLociConstraints);
		InFile.Close();
		return(eBSFerrFieldCnt);
		}

	if((UINT32)EndLoci >= m_pSfxArray->GetSeqLen(ChromID))
		{
		CSVLineNum = InFile.GetLineNumber();
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"End loci must be > targeted sequence length for '%s' at line %d in '%s'",szPrevTargChrom,CSVLineNum,pszLociConstraints);
		InFile.Close();
		return(eBSFerrFieldCnt);
		}

	Constraint = 0;
	pBase = pszBases;
	while((Base = *pBase++) != '\0')
		{
		switch(Base) {
			case 'a': case 'A':
				Constraint |= 0x01;
				break;
			case 'c': case 'C':
				Constraint |= 0x02;
				break;
			case 'g': case 'G':
				Constraint |= 0x04;
				break;
			case 't': case 'T':
				Constraint |= 0x08;
				break;
			case 'r': case 'R':
				Constraint |= 0x10;
				break;
			case ' ': case '\t':
				continue;
			default:
				CSVLineNum = InFile.GetLineNumber();
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Illegal base specifiers for '%s' at line %d in '%s'",szPrevTargChrom,CSVLineNum,pszLociConstraints);
				InFile.Close();
				return(eBSFerrFieldCnt);
			}
		}
	if(Constraint == 0)
		{
		CSVLineNum = InFile.GetLineNumber();
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Illegal base specifiers for '%s' at line %d in '%s'",szPrevTargChrom,CSVLineNum,pszLociConstraints);
		InFile.Close();
		return(eBSFerrFieldCnt);
		}

	// chrom, start, end and bases parsed
	for(Idx = 0; Idx < m_NumConstrainedChroms; Idx++)
		{
		if(m_ConstrainedChromIDs[Idx] == ChromID)
			break;
		}
	if(Idx == m_NumConstrainedChroms)
		{
		if(m_NumConstrainedChroms == cMaxConstrainedChroms)
			{
			CSVLineNum = InFile.GetLineNumber();
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Number of constrained chroms would be more than max (%d) allowed for '%s' at line %d in '%s'",cMaxConstrainedChroms,szPrevTargChrom,CSVLineNum,pszLociConstraints);
			InFile.Close();
			return(eBSFerrFieldCnt);
			}
		m_ConstrainedChromIDs[Idx] = ChromID;
		m_NumConstrainedChroms += 1;
		}

	if(m_NumConstraintLoci == cMaxConstrainedLoci)
		{
		CSVLineNum = InFile.GetLineNumber();
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Number of constrained loci would be more than max (%d) allowed for '%s' at line %d in '%s'",cMaxConstrainedLoci,szPrevTargChrom,CSVLineNum,pszLociConstraints);
		InFile.Close();
		return(eBSFerrFieldCnt);
		}
	m_pConstraintLoci[m_NumConstraintLoci].ChromID = ChromID;
	m_pConstraintLoci[m_NumConstraintLoci].StartLoci = StartLoci;
	m_pConstraintLoci[m_NumConstraintLoci].EndLoci = EndLoci;
	m_pConstraintLoci[m_NumConstraintLoci++].Constraint = Constraint;
	}

InFile.Close();

// sort the constraints by chromid.start.end ascending
if(m_NumConstraintLoci > 1)
	m_mtqsort.qsort(m_pConstraintLoci,m_NumConstraintLoci,sizeof(tsConstraintLoci),SortConstraintLoci);

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Completed loading %d loci base constraints for %d target sequences from CSV file '%s'",m_NumConstraintLoci,m_NumConstrainedChroms, pszLociConstraints);
return((teBSFrsltCodes)m_NumConstraintLoci);
}

teBSFrsltCodes
CAligner::Disk2Hdr(char *pszRdsFile)			// read from disk and validate header
{
if(_lseeki64(m_hInFile,0,SEEK_SET)!=0)			// read in header..
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Seek failed to offset 0 - %s",pszRdsFile,strerror(errno));
	Reset(false);			// closes opened files..
	return(eBSFerrFileAccess);
	}

if(sizeof(tsBSFRdsHdr) != read(m_hInFile,&m_FileHdr,sizeof(tsBSFRdsHdr)))
	return(eBSFerrNotBioseq);

// header read, validate it as being a reads file header
if(tolower(m_FileHdr.Magic[0]) != 'b' ||
	tolower(m_FileHdr.Magic[1]) != 'i' ||
	tolower(m_FileHdr.Magic[2]) != 'o' ||
	tolower(m_FileHdr.Magic[3]) != 'r')
	return(eBSFerrNotBioseq);

	// can we handle this version?
if(m_FileHdr.Version < cBSFRdsVersionBack || m_FileHdr.Version > cBSFRdsVersion)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"%s opened as a preprocessed reads file - expected between version %d and %d, file version is %d",
					pszRdsFile,cBSFRdsVersionBack,cBSFRdsVersion,m_FileHdr.Version);
	Reset(false);			// closes opened files..
	return(eBSFerrFileVer);
	}

return(eBSFSuccess);
}

// NOTE: will only return SNP bases, e.g. aligned read sequence must not contain microInDel or span splice junctions
inline eSeqBase
CAligner::AdjAlignSNPBase(tsReadHit *pReadHit,	// aligned read
		   UINT32 ChromID,			// read expected to have aligned to this chromosome
			UINT32 Loci)            // base to be returned is at this alignment loci, base will be complemented if antisense alignment
{
UINT32 AdjStartLoc;
UINT32 AdjEndLoc;
etSeqBase Base;
UINT8 *pBases;
tsSegLoci *pSeg;
tsHitLoci *pHit;

pHit = &pReadHit->HitLoci.Hit;
if(pHit->FlgInDel || pHit->FlgSplice)
	return(eBaseEOS);

pSeg = &pReadHit->HitLoci.Hit.Seg[0];
if(pSeg->ChromID != ChromID)
	return(eBaseEOS);

if(pSeg->Strand == '+')
	{
	AdjStartLoc=((UINT32)pSeg->MatchLoci + pSeg->TrimLeft);
	AdjEndLoc=((UINT32)pSeg->MatchLoci + (pSeg->MatchLen - pSeg->TrimRight - 1));
	}
else
	{
	AdjStartLoc=((UINT32)pSeg->MatchLoci + pSeg->TrimRight);
	AdjEndLoc=((UINT32)pSeg->MatchLoci + (pSeg->MatchLen - pSeg->TrimLeft - 1));
	}

if(AdjStartLoc > Loci || AdjEndLoc < Loci)
	return(eBaseEOS);

pBases = &pReadHit->Read[pReadHit->DescrLen+1];

if(pSeg->Strand == '+')
	Base = pBases[Loci - pSeg->MatchLoci]  & 0x07;
else
	{
	Base = pBases[pSeg->MatchLoci + pSeg->MatchLen - Loci - 1] & 0x07;
	switch(Base) {
		case eBaseA: Base = eBaseT; break;
		case eBaseC: Base = eBaseG; break;
		case eBaseG: Base = eBaseC; break;
		case eBaseT: Base = eBaseA; break;
		default:
			break;
		}
	}
return((eSeqBase)Base);
}

inline UINT32
CAligner::AdjStartLoci(tsSegLoci *pSeg)
{
if(pSeg->Strand == '+')
	return((UINT32)pSeg->MatchLoci + pSeg->TrimLeft);
else
	return((UINT32)pSeg->MatchLoci + pSeg->TrimRight);
}

inline UINT32
CAligner::AdjEndLoci(tsSegLoci *pSeg)
{
if(pSeg->Strand == '+')
	return((UINT32)pSeg->MatchLoci + pSeg->MatchLen + pSeg->TrimRight - 1);
else
	return((UINT32)pSeg->MatchLoci  + pSeg->MatchLen +  pSeg->TrimLeft - 1);
}

inline UINT32
CAligner::AdjHitLen(tsSegLoci *pSeg)
{
return((UINT32)pSeg->MatchLen - pSeg->TrimLeft - pSeg->TrimRight);
}

inline UINT32
CAligner::AdjAlignStartLoci(tsHitLoci *pHit)
{
tsSegLoci *pSeg = &pHit->Seg[0];
if(pSeg->Strand == '+')
	return((UINT32)pSeg->MatchLoci + pSeg->TrimLeft);
else
	return((UINT32)pSeg->MatchLoci + pSeg->TrimRight);
}

inline UINT32
CAligner::AdjAlignEndLoci(tsHitLoci *pHit)
{
tsSegLoci *pSeg;
if(pHit->FlgInDel || pHit->FlgSplice)
	pSeg = &pHit->Seg[1];
else
	pSeg = &pHit->Seg[0];
if(pSeg->Strand == '+')
	return((UINT32)pSeg->MatchLoci + pSeg->MatchLen + pSeg->TrimRight - 1);
else
	return((UINT32)pSeg->MatchLoci  + pSeg->MatchLen +  pSeg->TrimLeft - 1);
}

inline UINT32
CAligner::AdjAlignHitLen(tsHitLoci *pHit)
{
UINT32 HitLen;
tsSegLoci *pSeg = &pHit->Seg[0];
HitLen = (UINT32)pSeg->MatchLen - pSeg->TrimLeft - pSeg->TrimRight;
if(pHit->FlgInDel || pHit->FlgSplice)
	{
	pSeg = &pHit->Seg[1];
	HitLen += (UINT32)pSeg->MatchLen - pSeg->TrimLeft - pSeg->TrimRight;
	}
return(HitLen);
}

inline int
CAligner::AdjAlignMismatches(tsHitLoci *pHit)
{
int Mismatches;
tsSegLoci *pSeg = &pHit->Seg[0];
Mismatches = pSeg->TrimMismatches;
if(pHit->FlgInDel || pHit->FlgSplice)
	{
	pSeg = &pHit->Seg[1];
	Mismatches += pSeg->TrimMismatches;
	}
return(Mismatches);
}

// AutoTrimFlanks
// Intent is that this will be useful for delimiting RNAseq reads covering exon boundaries
// Autotrimmed aligned reads must be at least 50% of their untrimmed length or they will be discarded; exception is that if paired end processing then
// triming of these reads is limited so that at least 1/3rd of the read is retained as a central core
int
CAligner::AutoTrimFlanks(int MinFlankExacts) // Autotrim back aligned read flanks until there are at least MinFlankExacts exactly matching bases in the flanks
{
int Rslt;
tsReadHit *pReadHit;
UINT32 SeqIdx;
UINT32 Idx;
int MinTrimmedLen;
UINT8 *pSeq;

if((Rslt=SortReadHits(eRSMHitMatch,false)) < eBSFSuccess)
	{
	Reset(false);
	return(Rslt);
	}

m_ElimPlusTrimed = 0;
m_ElimMinusTrimed = 0;
if(MinFlankExacts > 0)	// do  3' and 5' autotrim? Note that currently can't trim multiseg hits
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Starting 5' and 3' flank sequence autotrim processing...");
	etSeqBase *pHitSeq;
	pReadHit = NULL;
	int ExactLen;
	int LeftOfs;
	int RightOfs;
	UINT32 MatchLen;
	int TrimMismatches;
	UINT8 ReadSeq[cMaxFastQSeqLen+1];
	UINT8 TargSeq[cMaxFastQSeqLen+1];

	while((pReadHit = IterReads(pReadHit))!=NULL)
		{
		pReadHit->HitLoci.FlagTR = 0;
		if(pReadHit->NAR != eNARAccepted || pReadHit->HitLoci.FlagSegs == 1)
			continue;
		MatchLen = pReadHit->HitLoci.Hit.Seg[0].MatchLen;
		if(MatchLen != pReadHit->ReadLen)
			{
			pReadHit->NumHits = 0;
			pReadHit->NAR = eNARTrim;
			if(pReadHit->HitLoci.Hit.Seg[0].Strand == '+')
				m_ElimPlusTrimed += 1;
			else
				m_ElimMinusTrimed += 1;
			continue;
			}

		MinTrimmedLen = (MatchLen+1)/2;
		if(MinTrimmedLen < 15)
			MinTrimmedLen = 15;

				// copy read into ReadSeq masking any hi-order phred scores (in bits 4..7) and repeat mask (bit 3) out
		pHitSeq = &pReadHit->Read[pReadHit->DescrLen+1];
		pSeq = ReadSeq;
		for(SeqIdx = 0; SeqIdx < MatchLen; SeqIdx++)
			*pSeq++ = *pHitSeq++ & 0x07;

		if(m_bIsSOLiD)
			{
			UINT32 Loci = AdjStartLoci(&pReadHit->HitLoci.Hit.Seg[0]);
			if(pReadHit->HitLoci.Hit.Seg[0].Strand == '+')
				Loci += 1;
			if(m_pSfxArray->GetColorspaceSeq(pReadHit->HitLoci.Hit.Seg[0].ChromID,
									Loci,
									TargSeq,MatchLen) == 0) // get colorspace sequence
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: GetColorspaceSeq() for chrom %d, Loci %d, matchlen %d failed",
										pReadHit->HitLoci.Hit.Seg[0].ChromID,Loci,MatchLen);
				Reset(false);
				return(eBSFerrMem);
				}

			if(pReadHit->HitLoci.Hit.Seg[0].Strand == '-')
				CSeqTrans::ReverseSeq(MatchLen,TargSeq);

			// convert target assembly sequence read into colorspace
			pSeq = ReadSeq;
			UINT8 PrvBase = *pSeq;
			for(SeqIdx = 1; SeqIdx <= MatchLen; SeqIdx++,pSeq++)
				{
				*pSeq = SOLiDmap[PrvBase][pSeq[1]];
				PrvBase = pSeq[1];
				}
			MatchLen -= 1;
			}
		else
			{
			// get basespace sequence
			m_pSfxArray->GetSeq(pReadHit->HitLoci.Hit.Seg[0].ChromID,(UINT32)(UINT64)pReadHit->HitLoci.Hit.Seg[0].MatchLoci,TargSeq,(UINT32)MatchLen);
			if(pReadHit->HitLoci.Hit.Seg[0].Strand == '-')
				CSeqTrans::ReverseComplement(MatchLen,TargSeq);
			}

		// trim from 5' towards 3'
		pSeq = ReadSeq;
		pHitSeq = TargSeq;
		ExactLen = 0;
		TrimMismatches = 0;
		int PEmincore;
		if(m_PEproc == ePEdefault)
			PEmincore = MatchLen;
		else
			PEmincore = MatchLen / 3;
		for(Idx = 0; Idx <= (MatchLen-(UINT32)MinTrimmedLen) && Idx < (UINT32)PEmincore; Idx++,pSeq++,pHitSeq++)
			{
			if(*pSeq != *pHitSeq)
				{
				ExactLen = 0;
				TrimMismatches += 1;
				continue;
				}
			ExactLen += 1;
			if(ExactLen == MinFlankExacts)
				break;
			}

		if(m_PEproc == ePEdefault)
			{
			// if can't trim to a 5' flank of at least MinFlankExacts then this read is to be sloughed
			if((Idx + MinTrimmedLen) > MatchLen || ExactLen < MinFlankExacts)
				{
				pReadHit->NumHits = 0;
				pReadHit->NAR = eNARTrim;
				if(pReadHit->HitLoci.Hit.Seg[0].Strand == '+')
					m_ElimPlusTrimed += 1;
				else
					m_ElimMinusTrimed += 1;
				continue;
				}
			}
		LeftOfs = Idx - (MinFlankExacts - 1);

		// trim from 3' towards 5'
		pHitSeq = &TargSeq[MatchLen-1];
		pSeq = &ReadSeq[MatchLen-1];
		ExactLen = 0;
		if(m_PEproc == ePEdefault)
			PEmincore = 0;
		else
			PEmincore = (MatchLen * 2) / 3;

		for(Idx = MatchLen-1; Idx >= (UINT32)(LeftOfs+MinTrimmedLen) && Idx > (UINT32)PEmincore; Idx--,pSeq--,pHitSeq--)
			{
			if(*pSeq != *pHitSeq)
				{
				ExactLen = 0;
				TrimMismatches += 1;
				continue;
				}
			ExactLen += 1;
			if(ExactLen == MinFlankExacts)
				break;
			}

		if(m_PEproc == ePEdefault)
			{
			// if can't trim to a 3' flank of at least MinFlankExacts then this read is to be sloughed
			if(ExactLen != MinFlankExacts || Idx < (UINT32)(LeftOfs + MinTrimmedLen))
				{
				pReadHit->NumHits = 0;
				pReadHit->NAR = eNARTrim;
				if(pReadHit->HitLoci.Hit.Seg[0].Strand == '+')
					m_ElimPlusTrimed += 1;
				else
					m_ElimMinusTrimed += 1;
				continue;
				}
			}
		RightOfs = Idx + MinFlankExacts;

		// left and right offset in the read at which the exact matching flanks start is now known
		if(m_bIsSOLiD)
			{
			MatchLen += 1;
			RightOfs += 1;
			}

		pReadHit->HitLoci.Hit.Seg[0].TrimLeft = LeftOfs;
		pReadHit->HitLoci.Hit.Seg[0].TrimRight = MatchLen - RightOfs;
		if(LeftOfs || (MatchLen - RightOfs))
			{
			pReadHit->HitLoci.Hit.Seg[0].TrimMismatches = pReadHit->HitLoci.Hit.Seg[0].Mismatches - TrimMismatches;
			pReadHit->HitLoci.FlagTR = 1;
			}
		else
			pReadHit->HitLoci.FlagTR = 0;
		}

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Finished 5' and 3' flank sequence autotriming, %d plus strand and %d minus strand aligned reads removed",m_ElimPlusTrimed,m_ElimMinusTrimed);

	if((Rslt=SortReadHits(eRSMHitMatch,false,true)) < eBSFSuccess)
		{
		Reset(false);
		return(Rslt);
		}

	if(gProcessingID)
		{
		gSQLiteSummaries.AddResult(gProcessingID,(char *)"Filtering",ePTInt32,sizeof(m_ElimPlusTrimed),"5' trimmed",&m_ElimPlusTrimed);
		gSQLiteSummaries.AddResult(gProcessingID,(char *)"Filtering",ePTInt32,sizeof(m_ElimMinusTrimed),"3' trimmed",&m_ElimMinusTrimed);
		}
	}
return(eBSFSuccess);
}

// TrimChimeric
// Aligned reads marked as chimeric aligned are left and right flank trimmed so that in subsequent processing these reads are treated as matching over their full length
const int cChimericSeqBuffLen = 0x07fffff;  // use this sized buffer when reporting the chimeric sequences

int
CAligner::TrimChimeric(char *pszChimericSeqs)	// trim back aligned chimeric read flanks with option to write chimeric sequences to file pszChimericSeqs 
{
int hChimerics;
char Strand;
char szChromName[128];
UINT32 Ofs;
UINT32 PrevChromID;

char *pszLineBuff;
int BuffIdx;
tsReadHit *pCurReadHit;
tsReadHit *pNxtReadHit;
UINT32 TrimLeft;
UINT32 TrimRight;
UINT32 NumLeftRightTrimmed;
UINT32 SeqIdx;
UINT32 NumTrimmed;
UINT32 NumLeftTrimmed;
UINT32 NumRightTrimmed;
UINT32 MatchLen;
UINT8 *pSeq;
UINT8 *pChimericSeq;
UINT8 *pTo;
int CopyLen;
UINT32 NumReads;


gDiagnostics.DiagOut(eDLInfo,gszProcName,"Starting chimeric flank sequence trim processing...");
PrevChromID = 0;
hChimerics = -1;
pszLineBuff = NULL;
if(pszChimericSeqs != NULL && pszChimericSeqs[0] != '\0')
	{
#ifdef _WIN32
	hChimerics = open(pszChimericSeqs,( O_WRONLY | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC),(_S_IREAD | _S_IWRITE) );
#else
	if((hChimerics = open(pszChimericSeqs,O_WRONLY | O_CREAT,S_IREAD | S_IWRITE))!=-1)
		if(ftruncate(hChimerics,0)!=0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate chimeric sequences file %s - %s",pszChimericSeqs,strerror(errno));
			return(eBSFerrCreateFile);
			}
#endif

	if(hChimerics < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate chimeric sequences file '%s'",pszChimericSeqs);
		return(eBSFerrCreateFile);
		}	
	if((pszLineBuff = new char [cChimericSeqBuffLen]) == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to allocate memory %d for chimeric sequences buffering",cChimericSeqBuffLen);
		return(eBSFerrMem);
		}	
	BuffIdx = sprintf(pszLineBuff,"\"Chrom\",\"Loci\",\"Strand\",\"ReadDescr\",\"5'TrimLen\",\"5'TrimSeq\",\"AlignLen\",\"AlignSeq\",\"3'TrimLen\",\"3'TrimSeq\"\n");
	CUtility::SafeWrite(hChimerics,pszLineBuff,BuffIdx);
	BuffIdx = 0;
	}
NumTrimmed = 0;
NumLeftTrimmed = 0;
NumRightTrimmed = 0;
NumLeftRightTrimmed = 0;

pCurReadHit = m_pReadHits;
pTo = (UINT8 *)pCurReadHit;
NumReads = 0;
while(pCurReadHit != NULL) {
	NumReads += 1;
	CopyLen = sizeof(tsReadHit) + pCurReadHit->DescrLen + pCurReadHit->ReadLen;
	if(pCurReadHit->ReadID != m_FinalReadID)
		pNxtReadHit = (tsReadHit *)((UINT8 *)pCurReadHit + CopyLen);
	else
		pNxtReadHit = NULL;
	pCurReadHit->HitLoci.FlagTR = 0;
	if(pCurReadHit->NAR != eNARAccepted || pCurReadHit->HitLoci.Hit.FlgChimeric == 0)
		{
		if(NumTrimmed > 0)
			memmove(pTo,pCurReadHit,CopyLen);
		pTo += CopyLen;
		pCurReadHit = pNxtReadHit;
		if(pCurReadHit != NULL)
			pCurReadHit->PrevSizeOf = CopyLen;
		continue;
		}

	TrimLeft = pCurReadHit->HitLoci.Hit.Seg[0].TrimLeft;
	TrimRight = pCurReadHit->HitLoci.Hit.Seg[0].TrimRight;
	if(TrimLeft == 0 && TrimRight == 0)	// if accepted as a chimeric then should have had at least one flank to be trimmed, treat as full length match ...
		{
		pCurReadHit->HitLoci.Hit.FlgChimeric = 0;
		if(NumTrimmed > 0)
			memmove(pTo,pCurReadHit,CopyLen);
		pTo += CopyLen;
		pCurReadHit = pNxtReadHit;
		if(pCurReadHit != NULL)
			pCurReadHit->PrevSizeOf = CopyLen;
		continue;
		}
	if(TrimLeft > 0)
		NumLeftTrimmed += 1;
	if(TrimRight > 0)
		NumRightTrimmed += 1;
	if(TrimLeft > 0 && TrimRight > 0)
		NumLeftRightTrimmed += 1;

	MatchLen = pCurReadHit->ReadLen - (TrimLeft + TrimRight);
	pSeq = &pCurReadHit->Read[pCurReadHit->DescrLen+1];
	pChimericSeq = pSeq + TrimLeft;
	
	if(hChimerics != -1)
		{
		Strand = pCurReadHit->HitLoci.Hit.Seg[0].Strand;
		Ofs = (UINT32)pCurReadHit->HitLoci.Hit.Seg[0].MatchLoci;
		if(PrevChromID == 0 || pCurReadHit->HitLoci.Hit.Seg[0].ChromID != PrevChromID)
			{
			m_pSfxArray->GetIdentName(pCurReadHit->HitLoci.Hit.Seg[0].ChromID,sizeof(szChromName),szChromName);
			PrevChromID = pCurReadHit->HitLoci.Hit.Seg[0].ChromID;
			}
		BuffIdx += sprintf(&pszLineBuff[BuffIdx],"\"%s\",%d,\"%c\",",szChromName,Ofs,Strand);
		}
	if(TrimLeft > 0)
		{
		if(hChimerics != -1)
			{
			BuffIdx += sprintf(&pszLineBuff[BuffIdx],"\"%s\",%d,\"",(char *)pCurReadHit->Read,TrimLeft);
			CSeqTrans::MapSeq2Ascii(pSeq,TrimLeft,&pszLineBuff[BuffIdx]);
			BuffIdx += TrimLeft;
			BuffIdx += sprintf(&pszLineBuff[BuffIdx],"\",%d,\"",MatchLen);
			CSeqTrans::MapSeq2Ascii(pChimericSeq,MatchLen,&pszLineBuff[BuffIdx]);
			BuffIdx += MatchLen;
			BuffIdx += sprintf(&pszLineBuff[BuffIdx],"\",%d,\"",TrimRight);
			if(TrimRight)
				CSeqTrans::MapSeq2Ascii(&pChimericSeq[MatchLen],TrimRight,&pszLineBuff[BuffIdx]);
			BuffIdx += TrimRight;
			BuffIdx += sprintf(&pszLineBuff[BuffIdx],"\"\n");
			}
		for(SeqIdx = 0; SeqIdx < MatchLen; SeqIdx++)
			*pSeq++ = *pChimericSeq++;
		}
	else
		if(hChimerics != -1)
			{
			BuffIdx += sprintf(&pszLineBuff[BuffIdx],"\"%s\",0,\"\",%d,\"",(char *)pCurReadHit->Read,MatchLen);
			CSeqTrans::MapSeq2Ascii(pSeq,MatchLen,&pszLineBuff[BuffIdx]);
			BuffIdx += MatchLen;
			BuffIdx += sprintf(&pszLineBuff[BuffIdx],"\",%d,\"",TrimRight);
			if(TrimRight)
				CSeqTrans::MapSeq2Ascii(&pSeq[MatchLen],TrimRight,&pszLineBuff[BuffIdx]);
			BuffIdx += TrimRight;
			BuffIdx += sprintf(&pszLineBuff[BuffIdx],"\"\n");
			}

	if(hChimerics != -1 && BuffIdx > (cChimericSeqBuffLen / 2))
		{
		CUtility::SafeWrite(hChimerics,pszLineBuff,BuffIdx);
		BuffIdx = 0;
		}

	pCurReadHit->ReadLen = MatchLen;
	pCurReadHit->HitLoci.Hit.Seg[0].MatchLen = MatchLen;
	pCurReadHit->HitLoci.Hit.Seg[0].TrimLeft = 0;
	pCurReadHit->HitLoci.Hit.Seg[0].TrimRight = 0;

	pCurReadHit->HitLoci.Hit.Seg[0].TrimMismatches = pCurReadHit->HitLoci.Hit.Seg[0].Mismatches;

	CopyLen = sizeof(tsReadHit) + pCurReadHit->DescrLen + pCurReadHit->ReadLen;
	if(NumTrimmed > 0)
		memmove(pTo,pCurReadHit,CopyLen);
	pTo += CopyLen;
	pCurReadHit = pNxtReadHit;
	if(pCurReadHit != NULL)
		pCurReadHit->PrevSizeOf = CopyLen;
	NumTrimmed += 1;
	}
if(hChimerics != -1)
	{
	if(BuffIdx)
		CUtility::SafeWrite(hChimerics,pszLineBuff,BuffIdx);

#ifdef _WIN32
	_commit(hChimerics);
#else
	fsync(hChimerics);
#endif
	close(hChimerics);
	hChimerics = -1;
	}
if(pszLineBuff != NULL)
	delete pszLineBuff;
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Completed trimming %u chimeric flanks, %u 5', %u 3', %u both 5' and 3' flanks trimmed",NumTrimmed,NumLeftTrimmed,NumRightTrimmed, NumLeftRightTrimmed);
return(eBSFSuccess);
}



// PCR5PrimerCorrect
// Intent is that this will be useful for aligning with tight substitution constraints whereby end PCR hexamer artefacts dominate the base substitution profile
// Substitutions in the first 5' KLen bases are progressively corrected, using the targeted sequence as template
int 
CAligner::PCR5PrimerCorrect(int MaxAllowedSubRate,	// after corrections overall sub rate for read must be no more than this 
					int KLen) // progressively correct substitutions in first Klen 5' read bases - assumed to be PCR random primer artefacts - until overall read meets MaxAllowedSubRate
{
int Rslt;
tsReadHit *pReadHit;
UINT8 *pSeq;
int MaxMMs;
int Ofs;
int CurMMs;
etSeqBase Base;
int NumCorrectedReads;
int NumCorrectedBases;
int NumUnacceptedReads;
etSeqBase *pHitSeq;
pReadHit = NULL;
UINT32 MatchLen;
UINT8 TargSeq[cMaxFastQSeqLen+1];


if(m_bIsSOLiD || KLen < 1)
	return(0);

if((Rslt=SortReadHits(eRSMHitMatch,false)) < eBSFSuccess)
	{
	Reset(false);
	return(Rslt);
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Starting PCR 5' %dbp primer correction processing targeting substitution rate of %d ...",KLen, MaxAllowedSubRate);

NumCorrectedReads = 0;
NumCorrectedBases = 0;
NumUnacceptedReads = 0;
while((pReadHit = IterReads(pReadHit))!=NULL)
	{
	if(pReadHit->NAR != eNARAccepted || pReadHit->HitLoci.FlagSegs == 1)
		continue;

	MaxMMs = ((MaxAllowedSubRate * pReadHit->ReadLen)+50) / 100;

	if(pReadHit->LowMMCnt <= MaxMMs)				// if already meeting targeted sub rate then no correction required
		continue;

	MatchLen = pReadHit->HitLoci.Hit.Seg[0].MatchLen;
	if(MatchLen != pReadHit->ReadLen)     // should usually be the same but some earlier processing function may have modified the length so don't attempt correction
		continue;

	// get target sequence
	m_pSfxArray->GetSeq(pReadHit->HitLoci.Hit.Seg[0].ChromID,(UINT32)(UINT64)pReadHit->HitLoci.Hit.Seg[0].MatchLoci,TargSeq,(UINT32)MatchLen);
	if(pReadHit->HitLoci.Hit.Seg[0].Strand == '-')
		CSeqTrans::ReverseComplement(MatchLen,TargSeq);

	// progressively look for mismatches over first KLen of read and correct until meeting  MaxAllowedSubRate
	pSeq =  &pReadHit->Read[pReadHit->DescrLen+1];
	pHitSeq = TargSeq;
	CurMMs = pReadHit->LowMMCnt;
	for(Ofs = 0; Ofs < (int)KLen; Ofs++,pSeq++,pHitSeq++)
		{
		if((*pSeq & 0x07) != *pHitSeq)					
			{
			if(--CurMMs <= MaxMMs)
				break;
			}
		}
	if(CurMMs <= MaxMMs)
		{
		pSeq =  &pReadHit->Read[pReadHit->DescrLen+1];
		pHitSeq = TargSeq;
		CurMMs = pReadHit->LowMMCnt;
		for(Ofs = 0; Ofs < (int)KLen; Ofs++,pSeq++,pHitSeq++)
			{
			if(((Base = *pSeq) & 0x07) != *pHitSeq)					// required an aligner induced sub?
				{
				*pSeq = (Base & 0xf8) | *pHitSeq;
				NumCorrectedBases += 1;
				if(--CurMMs <= MaxMMs)
					break;
				}
			}
		pReadHit->LowMMCnt = CurMMs;
		pReadHit->HitLoci.Hit.Seg[0].Mismatches = CurMMs;
		NumCorrectedReads += 1;
		}
	else
		{
		pReadHit->NumHits = 0;
		pReadHit->NAR = eNARNoHit;
		NumUnacceptedReads += 1;
		}
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Completed PCR 5' primer correction, %d reads with %d bases corrected, %d reads with excessive substitutions rejected",NumCorrectedReads,NumCorrectedBases, NumUnacceptedReads);

if((Rslt=SortReadHits(eRSMHitMatch,false,true)) < eBSFSuccess)
	{
	Reset(false);
	return(Rslt);
	}

if(gProcessingID)
	{
	gSQLiteSummaries.AddResult(gProcessingID,(char *)"PCR5PrimerCorrect",ePTInt32,sizeof(NumCorrectedReads),"NumCorrectedReads",&NumCorrectedReads);
	gSQLiteSummaries.AddResult(gProcessingID,(char *)"PCR5PrimerCorrect",ePTInt32,sizeof(NumCorrectedBases),"NumCorrectedBases",&NumCorrectedBases);
	gSQLiteSummaries.AddResult(gProcessingID,(char *)"PCR5PrimerCorrect",ePTInt32,sizeof(NumUnacceptedReads),"NumUnacceptedReads",&NumUnacceptedReads);
	}

return(eBSFSuccess);
}


//
int
CAligner::NumUniqueAlignedReads(void)		// return count of accepted uniquely aligned reads
{
tsReadHit *pReadHit = NULL;
int NumUniques = 0;
while((pReadHit = IterReads(pReadHit))!=NULL)
	{
	if(pReadHit->NAR == eNARAccepted)
		NumUniques += 1;
	}
return(NumUniques);
}

int
CAligner::DedupeNonuniqueReads(void)	// dedupe reads
{
int NumSubDups;

tsReadHit *pReadHit1;
tsReadHit *pReadHitMark;
tsReadHit *pReadHit = NULL;

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Deduping aligned reads which after substitutions are no longer unique");
SortReadHits(eRSMHitMatch,false);
NumSubDups = 0;
while((pReadHit = IterSortedReads(pReadHit))!=NULL)
	{
	if(pReadHit->NAR != eNARAccepted)		// only interested in reads accepted as aligned
		continue;

	pReadHit1 = pReadHit;			// start iterating forward checking for multiple probe hits onto same locus
	pReadHitMark = pReadHit;
	while((pReadHit1 = IterSortedReads(pReadHit1)) != NULL)
		{
		if(pReadHit1->NAR != eNARAccepted)	// only interested in reads with a single hits
			break;

		if(	pReadHit->HitLoci.Hit.Seg[0].ChromID == pReadHit1->HitLoci.Hit.Seg[0].ChromID && // iterated probe hitting same locus as original probe?
			AdjStartLoci(&pReadHit->HitLoci.Hit.Seg[0]) == AdjStartLoci(&pReadHit1->HitLoci.Hit.Seg[0]) &&
			pReadHit->HitLoci.Hit.Seg[0].Strand == pReadHit1->HitLoci.Hit.Seg[0].Strand)
			{
			// the iterated probe will be of equal or lower level than original probe
			// simply treat iterated probe as though it had multiple matches
			pReadHit1->NumHits = 0;
			pReadHit1->LowHitInstances = 0;
			pReadHit1->NAR = eNARNonUnique;
			pReadHitMark = pReadHit1;
			NumSubDups += 1;
			}
		else
			break;
		}
	pReadHit = pReadHitMark;
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Deduping completed, removed %d non-unique reads",NumSubDups);

if(NumSubDups)
	SortReadHits(eRSMHitMatch,false,true);

if(gProcessingID)
	gSQLiteSummaries.AddResult(gProcessingID,(char *)"Filtering",ePTInt32,sizeof(NumSubDups),"Nonunique",&NumSubDups);

return(eBSFSuccess);
}


//
// PCR results in differential read duplications with some sequence templates being many fold duplicated relative to others
// Differential artifact reduction is very crude!
// Each unique loci to which a read is aligned is iterated, then the number of unique read alignment sites up/down stream within WinLen bases is counted and this count is then used
// to limit the number of reads attributed to the loci currently being iterated.
// The limit is set to be a function of the maximum of the up/downstream unique aligned loci within WinLen...
int
CAligner::ReducePCRduplicates(int WinLen)		// remove potential PCR artefacts
{
int Rslt;
int NumSubDups;
int CurChromID;
int CurLen;
int CurStart;
UINT8 CurStrand;
int UpUniques;
int DnUniques;
int LimitDups;
int PropWin;

tsReadHit *pReadHit1;
tsReadHit *pReadHitMark;
tsReadHit *pReadHit = NULL;

// sort reads by match loci
if((Rslt=SortReadHits(eRSMHitMatch,false)) < eBSFSuccess)
	{
	Reset(false);
	return(Rslt);
	}

NumSubDups = 0;
while((pReadHit = IterSortedReads(pReadHit))!=NULL)
	{
	if(pReadHit->NAR != eNARAccepted)		// only interested in reads with unique hits
		continue;

	if(WinLen > 0)
		{
		UpUniques = NumUpUniques(pReadHit,WinLen,true);
		DnUniques = NumDnUniques(pReadHit,WinLen,true);
		LimitDups = max(UpUniques,DnUniques);
		PropWin = (int)(((double)LimitDups/WinLen) * 100.0);
		if(PropWin < 5)
			LimitDups = 1;
		else
			if(PropWin <= 10)
				LimitDups = 2;
			else
				if(PropWin <= 20)
					LimitDups = 3;
				else
					if(PropWin <= 40)
						LimitDups = 4;
					else
						if(PropWin <= 60)
							LimitDups = 5;
						else
							if(PropWin <= 80)
								LimitDups = 10;
							else
								LimitDups = 50;
		}
	else
		LimitDups = 0;

	pReadHit1 = pReadHit;			// start iterating forward checking for multiple probe hits onto same locus
	pReadHitMark = pReadHit;
	CurChromID = pReadHit->HitLoci.Hit.Seg[0].ChromID;
	CurStart = AdjStartLoci(&pReadHit->HitLoci.Hit.Seg[0]);
	CurLen = AdjHitLen(&pReadHit->HitLoci.Hit.Seg[0]);
	CurStrand = pReadHit->HitLoci.Hit.Seg[0].Strand;
	while((pReadHit1 = IterSortedReads(pReadHit1)) != NULL)
		{
		if(pReadHit1->NAR != eNARAccepted)	// only interested in reads with a unique hits
			continue;

		if(CurChromID == pReadHit1->HitLoci.Hit.Seg[0].ChromID && // iterated probe hitting same locus as original probe?
			CurStart == AdjStartLoci(&pReadHit1->HitLoci.Hit.Seg[0]) &&
			CurStrand == pReadHit1->HitLoci.Hit.Seg[0].Strand)
			{
			if(CurLen != AdjHitLen(&pReadHit1->HitLoci.Hit.Seg[0]))
				continue;
			if(LimitDups > 0)
				{
				LimitDups -= 1;
				continue;
				}

			pReadHit1->NumHits = 0;
			pReadHit1->LowHitInstances = 0;
			pReadHit1->NAR = eNARPCRdup;
			pReadHitMark = pReadHit1;
			NumSubDups += 1;
			}
		else
			break;
		}
	pReadHit = pReadHitMark;
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Removed %d potential PCR artefact reads",NumSubDups);
if(NumSubDups)
	SortReadHits(eRSMHitMatch,false,true);
if(gProcessingID)
	gSQLiteSummaries.AddResult(gProcessingID,(char *)"Filtering",ePTInt32,sizeof(NumSubDups),"PCRartefacts",&NumSubDups);
return(eBSFSuccess);
}


int
CAligner::RemoveOrphanSpliceJuncts(int SpliceJunctLen)	// remove unsupported orphan splice junctions
{
int Idx;
if(SpliceJunctLen > 0)
	{
	SortReadHits(eRSMHitMatch,false);
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Filtering out any orphan splice junction reads");
	int NumSpliceJuncs = 0;
	int	NumSpliceAccepted = 0;
	int	NumSpliceNotAccepted = 0;
	tsSegJuncts *pSegJuncts = NULL;
	tsSegJuncts *pJunct;
	tsSegJuncts *pNxtJunct;
	tsReadHit *pCurHit = NULL;
	tsReadHit *pReadHit = NULL;
	while((pReadHit = IterReads(pReadHit))!=NULL)
		{
		if(pReadHit->NAR != eNARAccepted || pReadHit->HitLoci.Hit.FlgSplice == 0)
			continue;
		pCurHit = pReadHit;
		NumSpliceJuncs += 1;
		}
	if(NumSpliceJuncs > 1)
		{
		// allocate for this number of splice junctions, then sort and dedupe
		pSegJuncts = new tsSegJuncts[NumSpliceJuncs];
		if(pSegJuncts == NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to allocate memory (%d bytes) for splice junction processing",NumSpliceJuncs * sizeof(tsSegJuncts));
			Reset(false);
			return(eBSFerrMem);
			}
		pJunct = pSegJuncts;
		pReadHit = NULL;
		while((pReadHit = IterReads(pReadHit))!=NULL)
			{
			if(pReadHit->NAR != eNARAccepted || pReadHit->HitLoci.Hit.FlgSplice == 0)
				continue;
			pJunct->pRead = pReadHit;
			pJunct->ChromID = pReadHit->HitLoci.Hit.Seg[0].ChromID;
			pJunct->Starts = AdjEndLoci(&pReadHit->HitLoci.Hit.Seg[0]);
			pJunct->Ends = AdjStartLoci(&pReadHit->HitLoci.Hit.Seg[1]);
			pJunct->Cnt = 1;
			pJunct += 1;
			}
		// now sort ascending by chrom, start, end
		m_mtqsort.qsort(pSegJuncts,NumSpliceJuncs,sizeof(tsSegJuncts),SortSegJuncts);
		// iterate and determine multiplicity of junctions
		pJunct = pSegJuncts;
		pNxtJunct = &pSegJuncts[1];
		for(Idx = 0; Idx < (NumSpliceJuncs-1); Idx++,pJunct++,pNxtJunct++)
			{
			if(pJunct->ChromID == pNxtJunct->ChromID &&
			   (pJunct->Starts <= (pNxtJunct->Starts + 3) && pJunct->Starts >= (pNxtJunct->Starts - 3)) &&
			   (pJunct->Ends <= (pNxtJunct->Ends + 3) && pJunct->Ends >= (pNxtJunct->Ends - 3)))
				{
				pJunct->pRead->HitLoci.Hit.FlgNonOrphan = 1;
				pNxtJunct->pRead->HitLoci.Hit.FlgNonOrphan = 1;
				}
			}
		// junctions not supported by at least two reads are to be treated as simply unaligned
		pJunct = pSegJuncts;
		for(Idx = 0; Idx < NumSpliceJuncs; Idx++,pJunct++)
			{
			if(pJunct->pRead->HitLoci.Hit.FlgNonOrphan != 1)
				{
				pJunct->pRead->NAR = eNARSpliceJctn;
				pJunct->pRead->NumHits = 0;	// treat as unaligned
				pJunct->pRead->LowHitInstances = 0;
				NumSpliceNotAccepted += 1;
				}
			else
				NumSpliceAccepted += 1;
			}
		delete pSegJuncts;
		}
	else
		if(NumSpliceJuncs > 0 && pCurHit != NULL)
			{
			pCurHit->NAR = eNARSpliceJctn;
			pCurHit->NumHits = 0;	// treat as unaligned
			pCurHit->LowHitInstances = 0;
			NumSpliceNotAccepted += 1;
			}
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"From %d reads with putative splice sites %d orphans were removed", NumSpliceJuncs,NumSpliceNotAccepted);
	if(gProcessingID)
		gSQLiteSummaries.AddResult(gProcessingID,(char *)"Filtering",ePTInt32,sizeof(NumSpliceNotAccepted),"OrphanSpliceJuncts",&NumSpliceNotAccepted);
	if(NumSpliceNotAccepted)
		SortReadHits(eRSMHitMatch,false,true);
	}

return(eBSFSuccess);
}

int
CAligner::RemoveOrphanMicroInDels(int microInDelLen) // remove any unsupported orphan microInDels
{
int Idx;
int NumInDelJuncs = 0;
int	NumInDelsAccepted = 0;
int	NumInDelsNotAccepted = 0;

tsSegJuncts *pSegJuncts = NULL;
tsSegJuncts *pJunct;
tsSegJuncts *pNxtJunct;
tsReadHit *pCurHit = NULL;
tsReadHit *pReadHit = NULL;

if(microInDelLen > 0)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Filtering out any orphan microIndel reads");
	SortReadHits(eRSMHitMatch,false);
	while((pReadHit = IterReads(pReadHit))!=NULL)
		{
		if(pReadHit->NAR != eNARAccepted || pReadHit->HitLoci.Hit.FlgInDel == 0)
			continue;
		pCurHit = pReadHit;
		NumInDelJuncs += 1;
		}
	if(NumInDelJuncs > 1)
		{
		// allocate for this number of InDel junctions, then sort and dedupe
		pSegJuncts = new tsSegJuncts[NumInDelJuncs];
		if(pSegJuncts == NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to allocate memory (%d bytes) for microInDel processing",NumInDelJuncs * sizeof(tsSegJuncts));
			Reset(false);
			return(eBSFerrMem);
			}
		pJunct = pSegJuncts;
		pReadHit = NULL;
		while((pReadHit = IterReads(pReadHit))!=NULL)
			{
			if(pReadHit->NAR != eNARAccepted || pReadHit->HitLoci.Hit.FlgInDel == 0)
				continue;
			pJunct->pRead = pReadHit;
			pJunct->ChromID = pReadHit->HitLoci.Hit.Seg[0].ChromID;
			pJunct->Starts = AdjEndLoci(&pReadHit->HitLoci.Hit.Seg[0]);
			pJunct->Ends = AdjStartLoci(&pReadHit->HitLoci.Hit.Seg[1]);
			pJunct->Cnt = 1;
			pJunct += 1;
			}
		// now sort ascending by chrom, start, end
		m_mtqsort.qsort(pSegJuncts,NumInDelJuncs,sizeof(tsSegJuncts),SortSegJuncts);
		// iterate and determine multiplicity of junctions
		pJunct = pSegJuncts;
		pNxtJunct = &pSegJuncts[1];
		for(Idx = 0; Idx < (NumInDelJuncs-1); Idx++,pJunct++,pNxtJunct++)
			{
			if(pJunct->ChromID == pNxtJunct->ChromID &&
			   (pJunct->Starts <= (pNxtJunct->Starts + 3) && pJunct->Starts >= (pNxtJunct->Starts - 3)) &&
			   (pJunct->Ends <= (pNxtJunct->Ends + 3) && pJunct->Ends >= (pNxtJunct->Ends - 3)))				{
				pJunct->pRead->HitLoci.Hit.FlgNonOrphan = 1;
				pNxtJunct->pRead->HitLoci.Hit.FlgNonOrphan = 1;
				}
			}
		// InDels not supported by at least two reads are to be treated as simply unaligned
		pJunct = pSegJuncts;
		for(Idx = 0; Idx < NumInDelJuncs; Idx++,pJunct++)
			{
			if(pJunct->pRead->HitLoci.Hit.FlgNonOrphan != 1)
				{
				pJunct->pRead->NAR = eNARmicroInDel;
				pJunct->pRead->NumHits = 0;	// treat as unaligned
				pJunct->pRead->LowHitInstances = 0;
				NumInDelsNotAccepted += 1;
				}
			else
				NumInDelsAccepted += 1;
			}
		delete pSegJuncts;
		}
	else
		if(NumInDelJuncs > 0 && pCurHit != NULL)
			{
			pCurHit->NAR = eNARmicroInDel;
			pCurHit->NumHits = 0;	// treat as unaligned
			pCurHit->LowHitInstances = 0;
			NumInDelsNotAccepted += 1;
			}
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"From %d reads with putative microIndels %d orphans were removed", NumInDelJuncs,NumInDelsNotAccepted);
	if(gProcessingID)
		gSQLiteSummaries.AddResult(gProcessingID,(char *)"Filtering",ePTInt32,sizeof(NumInDelsNotAccepted),"OrphanInDels",&NumInDelsNotAccepted);

	if(NumInDelsNotAccepted)
		SortReadHits(eRSMHitMatch,false,true);
	}
return(eBSFSuccess);
}



int								// -1: base not meeting constraints, 0: chromID has no constraints, 1: ChromID constrained but base accepted
CAligner::AcceptBaseConstraint(UINT32 ChromID,			// base aligned to this chrom/sequence
                                  UINT32 Loci,			// aligned to this loci
								  etSeqBase Base)		// base in read
{
int Idx;
etSeqBase TargBase;
int Rslt;
tsConstraintLoci *pConstraintLoci;

if(m_NumConstrainedChroms == 0 || m_NumConstraintLoci == 0 || m_pConstraintLoci == NULL) // there may be no constraints!
	return(0);

Rslt = 0;		// assuming no constraints on the targeted ChromID

// current implementation not expecting too many constraints, and only a few constrained chroms, so just a simple linear search for matching constraints
pConstraintLoci = m_pConstraintLoci;
for(Idx = 0; Idx < m_NumConstraintLoci; Idx++,pConstraintLoci++)
	{
	if(pConstraintLoci->ChromID < ChromID)
		continue;
	if(pConstraintLoci->ChromID > ChromID)
		break;

	// there is at least one constraint on targeted chrom
	Rslt = 1;
	if(pConstraintLoci->StartLoci <= Loci && pConstraintLoci->EndLoci >= Loci)
		{
		if(pConstraintLoci->Constraint & 0x10)	// accept if reads base same as targets base
			{
			TargBase = m_pSfxArray->GetBase(ChromID,Loci);
			if(TargBase == Base)
				continue;
			}

		if(pConstraintLoci->Constraint & 0x01 && Base == eBaseA)
			continue;
		if(pConstraintLoci->Constraint & 0x02 && Base == eBaseC)
			continue;
		if(pConstraintLoci->Constraint & 0x04 && Base == eBaseG)
			continue;
		if(pConstraintLoci->Constraint & 0x08 && Base == eBaseT)
			continue;
		return(-1);
		}
	}
return(Rslt);
}

bool												  // true if read alignment meets any loci base constraints
CAligner::AcceptLociConstraints(tsReadHit *pReadHit)   // read alignment to check
{
int Rslt;
int Idx;
int SeqIdx;
UINT32 ChromID;
UINT8 *pHitSeq;
UINT8 *pSeq;
UINT8 ReadSeq[cMaxSeqLen+1];
etSeqBase Base;

UINT32 StartLoci;
UINT32 EndLoci; 
UINT32 CurLoci;

if(pReadHit == NULL || pReadHit->NAR != eNARAccepted)	  // only interested in reads which have, thus far, been accepted as being aligned
	return(true);

if(m_NumConstrainedChroms == 0 || m_NumConstraintLoci == 0 || m_pConstraintLoci == NULL) // there may be no constraints!
	return(true);

ChromID = pReadHit->HitLoci.Hit.Seg[0].ChromID;

// check if aligned chrom has any constraint loci
for(Idx = 0; Idx < m_NumConstrainedChroms; Idx++)
	if(ChromID == m_ConstrainedChromIDs[Idx])
		break;
if(Idx == m_NumConstrainedChroms)
	return(true);


// copy read into ReadSeq masking any hi-order phred scores (in bits 4..7) and repeat mask (bit 3) out
pHitSeq = &pReadHit->Read[pReadHit->DescrLen+1];
pSeq = ReadSeq;
for(SeqIdx = 0; SeqIdx < pReadHit->ReadLen; SeqIdx++)
	*pSeq++ = *pHitSeq++ & 0x07;

// if read was aligned antisense to target then reverse complement
if(pReadHit->HitLoci.Hit.Seg[0].Strand == '-')
	CSeqTrans::ReverseComplement(pReadHit->ReadLen,ReadSeq);

StartLoci = AdjStartLoci(&pReadHit->HitLoci.Hit.Seg[0]);
EndLoci = AdjEndLoci(&pReadHit->HitLoci.Hit.Seg[0]);
SeqIdx = pReadHit->HitLoci.Hit.Seg[0].ReadOfs + pReadHit->HitLoci.Hit.Seg[0].TrimLeft;
for(CurLoci = StartLoci; CurLoci <= EndLoci; CurLoci++, SeqIdx++)
	{
	Base = ReadSeq[SeqIdx];
	Rslt = AcceptBaseConstraint(ChromID,CurLoci,Base);
	if(Rslt != 1)
		return(Rslt == -1 ? false : true);
	}

if(pReadHit->HitLoci.FlagSegs)
	{
	StartLoci = AdjStartLoci(&pReadHit->HitLoci.Hit.Seg[1]);
	EndLoci = AdjEndLoci(&pReadHit->HitLoci.Hit.Seg[1]);
	SeqIdx = pReadHit->HitLoci.Hit.Seg[1].ReadOfs + pReadHit->HitLoci.Hit.Seg[1].TrimLeft;
	for(CurLoci = StartLoci; CurLoci <= EndLoci; CurLoci++, SeqIdx++)
		{
		Base = ReadSeq[SeqIdx];
		Rslt = AcceptBaseConstraint(ChromID,CurLoci,Base);
		if(Rslt != 1)
			return(Rslt == -1 ? false : true);
		}
	}
return(true);
}

// Identify and mark with eNARLociConstrained any currently accepted alignments which violate a loci base constraint
int						// number of new alignments identified as violating a loci base constraint
CAligner::IdentifyConstraintViolations(bool bPEread) // true if processing for PE's, false if processing for SE
{
int NumIdentified;
tsReadHit *pReadHit;
tsReadHit *pPEReadHit;
bool bIsPE2;
if(m_NumConstrainedChroms == 0 || m_NumConstraintLoci == 0 || m_pConstraintLoci == NULL)
	return(0);

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Identifying %s loci base constraint violations ...",bPEread ? "PE" : "SE");
SortReadHits(eRSMHitMatch,false);

NumIdentified = 0;
pReadHit = NULL;
while((pReadHit = IterSortedReads(pReadHit))!=NULL)
	{
	if(pReadHit->NAR == eNARAccepted && !AcceptLociConstraints(pReadHit))
		{
		pReadHit->NAR = eNARLociConstrained;
		pReadHit->NumHits = 0;
		pReadHit->LowHitInstances = 0;
		NumIdentified += 1;
		}

	if(bPEread && pReadHit->NAR == eNARLociConstrained)		// if PE processsing then ensure that the partner read also marked 
		{
		// locate partner read for current read
		// if current read is PE1 then PE2 will be next read, if PE2 then PE1 will be previous read
		bIsPE2 = pReadHit->PairReadID & 0x80000000 ? true : false;	// top bit set if PE2
		if(bIsPE2)	
			pPEReadHit = (tsReadHit *)((UINT8 *)pReadHit - pReadHit->PrevSizeOf);   
		else
			pPEReadHit = (tsReadHit *)((UINT8 *)pReadHit + sizeof(tsReadHit) + pReadHit->ReadLen + pReadHit->DescrLen);
		if(pPEReadHit->NAR != eNARLociConstrained)
			{
			pPEReadHit->NAR = eNARLociConstrained;
			pPEReadHit->NumHits = 0;
			pPEReadHit->LowHitInstances = 0;
			NumIdentified += 1;
			}
		}
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Identified %d %s loci base constraint violations",NumIdentified,bPEread ? "PE" : "SE");
if(NumIdentified > 0)
	SortReadHits(eRSMHitMatch,false,true);
return(NumIdentified);
}



bool					// true if chrom is accepted, false if chrom not accepted
CAligner::AcceptThisChromID(UINT32 ChromID)
{
int IncChromIdx;
int ExclChromIdx;
char szChromName[128];
bool bProcChrom = false;
int MatchesFiltOut = 0;

if(!(m_NumExcludeChroms || m_NumIncludeChroms))
	return(true);

#ifdef _WIN32
RegexpMatch mc;
#else
regmatch_t mc;
int RegErr;					// regular expression parsing error
char szRegErr[128];			// to hold RegErr as textual representation ==> regerror();
#endif

m_pSfxArray->GetIdentName(ChromID,sizeof(szChromName),szChromName);

AcquireSerialise();

	// check if to be excluded
bProcChrom = true;
for(ExclChromIdx = 0; ExclChromIdx < m_NumExcludeChroms; ExclChromIdx++)
	{
#ifdef _WIN32
	if(m_ExcludeChromsRE[ExclChromIdx]->Match(szChromName,&mc))
#else
	if(!regexec(&m_ExcludeChromsRE[ExclChromIdx],szChromName,1,&mc,0))
#endif
		{
		bProcChrom = false;
		break;
		}
	}

	// to be included?
if(bProcChrom && m_NumIncludeChroms > 0)
	{
	bProcChrom = false;
	for(IncChromIdx = 0; IncChromIdx < m_NumIncludeChroms; IncChromIdx++)
		{
#ifdef _WIN32
		if(m_IncludeChromsRE[IncChromIdx]->Match(szChromName,&mc))
#else
		if(!regexec(&m_IncludeChromsRE[IncChromIdx],szChromName,1,&mc,0))
#endif
			{
			bProcChrom = true;
			break;
			}
		}
	}

ReleaseSerialise();

if(!bProcChrom)
	return(false);
return(true);
}


// returns
// > 0: accepted with return value the fragment (insert) size
//   0: not paired end alignment, both ends not uniquely aligned with NumHits == 1
// -1: alignment strands not consistent
// -2: aligned to different chromosomes
// -3: both ends on filtered chrom
// -4: 5' end on filtered chrom
// -5: 3' end on filtered chrom
// -6: fragment < PairMinLen
// -7: fragment > PairMaxLen
int				// > 0: accepted with return value the fragment (insert) size
CAligner::AcceptProvPE(int PairMinLen,		// only accept paired reads with a combined sequence length of at least this
			 int PairMaxLen,				// only accept paired reads with a combined sequence length of no more than this
			 bool bPairStrand,				// accept paired ends if on same strand
			 tsReadHit *pFwdReadHit,
			 tsReadHit *pRevReadHit)
{
UINT32 FwdChrom;
UINT8 FwdStrand;
UINT32 FwdProbeHitLen;
UINT32 FwdStartLoci;
UINT32 FwdEndLoci;

UINT32 RevChrom;
UINT8 RevStrand;
UINT32 RevProbeHitLen;
UINT32 RevStartLoci;
UINT32 RevEndLoci;

// to be a pair then expecting both ends to have been aligned
if(!(pFwdReadHit->NumHits == 1 && pRevReadHit->NumHits == 1))
	return(0);

FwdChrom = pFwdReadHit->HitLoci.Hit.Seg[0].ChromID;
FwdStrand = pFwdReadHit->HitLoci.Hit.Seg[0].Strand;
FwdProbeHitLen = AdjHitLen(&pFwdReadHit->HitLoci.Hit.Seg[0]);
FwdStartLoci = AdjStartLoci(&pFwdReadHit->HitLoci.Hit.Seg[0]);
if(pFwdReadHit->HitLoci.Hit.FlgInDel || pFwdReadHit->HitLoci.Hit.FlgSplice)
	{
	FwdEndLoci = AdjEndLoci(&pFwdReadHit->HitLoci.Hit.Seg[1]);
	FwdProbeHitLen += AdjHitLen(&pFwdReadHit->HitLoci.Hit.Seg[1]);
	}
else
	FwdEndLoci = AdjEndLoci(&pFwdReadHit->HitLoci.Hit.Seg[0]);
RevChrom = pRevReadHit->HitLoci.Hit.Seg[0].ChromID;
RevStrand = pRevReadHit->HitLoci.Hit.Seg[0].Strand;
RevProbeHitLen = AdjHitLen(&pRevReadHit->HitLoci.Hit.Seg[0]);
RevStartLoci = AdjStartLoci(&pRevReadHit->HitLoci.Hit.Seg[0]);
if(pRevReadHit->HitLoci.Hit.FlgInDel || pRevReadHit->HitLoci.Hit.FlgSplice)
	{
	RevEndLoci = AdjEndLoci(&pRevReadHit->HitLoci.Hit.Seg[1]);
	RevProbeHitLen += AdjHitLen(&pRevReadHit->HitLoci.Hit.Seg[1]);
	}
else
	RevEndLoci = AdjEndLoci(&pRevReadHit->HitLoci.Hit.Seg[0]);

bool bFwdChrom = AcceptThisChromID(FwdChrom);
if(FwdChrom != RevChrom)									
	{
	bool bRevChrom = AcceptThisChromID(RevChrom);
	if(bRevChrom == true && bFwdChrom == true)
		return(-2);
	if(bRevChrom == false && bFwdChrom == false)
		return(-3);
	if(bFwdChrom == false)
		return(-4);
	return(-5);
	}
else
	if(bFwdChrom == false)
		return(-3);

return(PEInsertSize(PairMinLen,PairMaxLen,bPairStrand,FwdStrand,FwdStartLoci,FwdEndLoci,RevStrand,RevStartLoci,RevEndLoci));
}

// calculates PE alignment insert size
// returns
// > 0: accepted with return value the fragment (insert) size
//   0: not paired end alignment, both ends not uniquely aligned with NumHits == 1
// -1: alignment strands not consistent
// -2: aligned to different chromosomes
// -3: both ends on filtered chrom
// -4: 5' end on filtered chrom
// -5: 3' end on filtered chrom
// -6: fragment < PairMinLen
// -7: fragment > PairMaxLen
int										// returned PE insert size, <= 0 if errors
CAligner::PEInsertSize(int PairMinLen,	// only accept paired reads with a combined sequence length of at least this
			 int PairMaxLen,			// only accept paired reads with a combined sequence length of no more than this
			 bool bPairStrand,			// accept paired ends if on same strand		
			 UINT8 PE1Strand,			// PE1 aligned on to this strand
		     UINT32 PE1StartLoci,		// PE read starts at this loci
		     UINT32 PE1EndLoci,			// PE1 read ends at this loci
             UINT8 PE2Strand,			// PE2 aligned on to this strand
		     UINT32 PE2StartLoci,		// PE2 read starts at this loci
		     UINT32 PE2EndLoci)			// PE2 read ends at this loci
{
int SeqFragLen;

if((bPairStrand && PE1Strand != PE2Strand) ||				
    (!bPairStrand && PE1Strand == PE2Strand))
		return(-1);

// if processing for circulars then 
//		if fwd is on '+' then distance = FwdEndLoci - RevStartLoci
//		if fwd is on '-' then distance = RevStartLoci - FwdEndLoci
// else
//		if fwd is on '+' then distance = RevEndLoci - FwdStartLoci
//		if fwd is on '-' then distance = FwdStartLoci - RevStartLoci
if(m_bPEcircularised)
	{
	if(PE1Strand == '+')
		SeqFragLen = 1 + (int)PE1EndLoci - (int)PE2StartLoci;
	else
		SeqFragLen = 1 + (int)PE2StartLoci - (int)PE1EndLoci;
	}
else
	{
	if(PE1Strand == '+')
		SeqFragLen = 1 + (int)PE2EndLoci - (int)PE1StartLoci;
	else
		SeqFragLen = 1 + (int)PE1EndLoci - (int)PE2StartLoci;
	}

if(SeqFragLen < 0)			// treat as inconsistent strand if fragment length is negative 
	return(-1);

if(SeqFragLen < PairMinLen)
	return(-6);

if(SeqFragLen > PairMaxLen)
	return(-7);

// can accept
return(SeqFragLen);
}

//ProcessPairedEnds
// Matches read pairs and if one read is aligned but partner is not due to muiltihit loci then
// attempts to find a unique loci for that multihit read within the expected approximate pair distance range
// If user has optionally requested then will accept reads which are ophaned (other PE not aligned) or PEs where PE1 and PE2 align to separate chroms/contigs 
int
CAligner::ProcessPairedEnds(etPEproc PEproc, // paired reads alignment processing mode
				  int MinEditDist, // accepted alignments must be at least this Hamming away from other putative alignments
				  int PairMinLen,  // only accept paired reads with a combined sequence length of at least this
				  int PairMaxLen,  // only accept paired reads with a combined sequence length of no more than this
				  bool bPairStrand,	// accept paired ends if on same strand
				  int MaxSubs)	   // aligned reads can have at most this many subs per 100bp
{
int Rslt;
UINT32 PrevChromID;
int SeqFragLen;
UINT32 PairReadIdx;
tsReadHit *pFwdReadHit;
tsReadHit *pRevReadHit;

UINT32 OrphStartLoci;
UINT32 OrphEndLoci;

int NumPaired = 0;
int NegPairs = 0;
int UnderLenPairs = 0;
int OverLenPairs = 0;
int LongestSeqFragLen = 0;
int OverPairMaxLen = 0;

bool bPartnerPaired;
int PartnerUnpaired = 0;
int PartnerPaired = 0;
int UnalignedPairs = 0;
int UniquePairCnt = 0;
int NumFilteredByChrom = 0;

int AcceptedNumPaired = 0;
int AcceptedNumSE = 0;

UINT8 ReadSeq[cMaxSeqLen+1];

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Generating paired reads index over %d paired reads", m_NumReadsLoaded/2);
SortReadHits(eRSMPairReadID,false,true);

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Iterating putative paired reads ...");

if(m_pLenDist != NULL)
	{
	delete m_pLenDist;
	m_pLenDist = NULL;
	}
if((m_pLenDist = new int[cPairMaxLen+1])==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate %d bytes memory for paired read sequence length distributions",
								sizeof(int) * (cPairMaxLen+1));
	Reset(false);
	return(eBSFerrMem);
	}
memset(m_pLenDist,0,sizeof(int) * (cPairMaxLen+1));

tsHitLoci HitLoci;
UINT8 *pHitSeq;
UINT8 *pSeq;
UINT32 SeqIdx;
bool b3primeExtend;
bool bAntisense;

m_szLineBuffIdx = 0;
PrevChromID = -1;
m_PrevSAMTargEntry = 0;
time_t Started = time(0);
UINT32 PrevPairReadIdx = 0;
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processed putative 0 pairs, accepted 0");
for(PairReadIdx=0; PairReadIdx < m_NumReadsLoaded;PairReadIdx+=2)
	{
	// is it time to let user know progress
	if(PairReadIdx > (PrevPairReadIdx + 50000))
		{
		time_t Now = time(0);
		unsigned long ElapsedSecs = (unsigned long) (Now - Started);
		if(ElapsedSecs >= (60 * 10))
			{
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processed putative %u pairs, accepted %u",PairReadIdx/2,AcceptedNumPaired);
			Started = Now;
			}
		PrevPairReadIdx = PairReadIdx;
		}

	pFwdReadHit = m_ppReadHitsIdx[PairReadIdx];
	pRevReadHit = m_ppReadHitsIdx[PairReadIdx+1];
	pFwdReadHit->FlgPEAligned = 0;
	pRevReadHit->FlgPEAligned = 0;

	// at least one of the ends must have been accepted as being aligned in order to process as PE
	if(!(pFwdReadHit->NAR <= eNARAccepted || pRevReadHit->NAR == eNARAccepted))
		{
		UnalignedPairs += 1;
		continue;
		}

	// even if both ends have been accepted as aligned it could be that these alignments can't be accepted as a PE - may be to different chroms, or overlength etc
	if(pFwdReadHit->NAR == eNARAccepted && pRevReadHit->NAR == eNARAccepted)
		{
		// both ends were accepted as being aligned but can these be accepted as PE within allowed insert size range ...
		SeqFragLen = AcceptProvPE(PairMinLen,PairMaxLen,bPairStrand,pFwdReadHit,pRevReadHit);	
		if(SeqFragLen > 0)
			{
			// this paired end alignment has been accepted
			pFwdReadHit->FlgPEAligned = 1;
			pRevReadHit->FlgPEAligned = 1;
			LongestSeqFragLen = max(LongestSeqFragLen,SeqFragLen);
			m_pLenDist[SeqFragLen] += 1;
			AcceptedNumPaired += 1;
			continue;	// onto next pair
			}

		// although ends aligned, am unable to accept as a PE alignment
		switch(SeqFragLen) {
			case -1:			// alignment strands not consistent
				pFwdReadHit->NAR = eNARPEStrand;
				pRevReadHit->NAR = eNARPEStrand;
				break;

			case -2:			// aligned to different chromosomes
				pFwdReadHit->NAR = eNARPEChrom;
				pRevReadHit->NAR = eNARPEChrom;
				break;

			case -3:			// both ends were to a filtered chrom so can't use either end as an anchor
				NumFilteredByChrom += 1;
				pFwdReadHit->NumHits = 0;
				pRevReadHit->NumHits = 0;
				pFwdReadHit->LowHitInstances = 0;
				pRevReadHit->LowHitInstances = 0;
				pFwdReadHit->NAR = eNARChromFilt;	
				pRevReadHit->NAR = eNARChromFilt;	
				continue;	// try next pair

			case -4:		// 5' end to filtered chrom so can't be used as an anchor
				pFwdReadHit->NAR = eNARChromFilt;
				pFwdReadHit->LowHitInstances = 0;
				pFwdReadHit->NumHits = 0;
				break;

			case -5:		// 3' end to filtered chrom so can't be used as an anchor
				pRevReadHit->NAR = eNARChromFilt;
				pRevReadHit->NumHits = 0;
				pRevReadHit->LowHitInstances = 0;
				break;

			case -6:		// under min insert size
				pFwdReadHit->NAR = eNARPEInsertMin;
				pRevReadHit->NAR = eNARPEInsertMin;
				break;

			case -7:		// over max insert size
				pFwdReadHit->NAR = eNARPEInsertMax;
				pRevReadHit->NAR = eNARPEInsertMax;
				break;
			}

		// if not allowed to orphan recover or treat as SE alignments then try next pair of reads
		if(PEproc == ePEunique)
			{
			pFwdReadHit->NumHits = 0;
			pRevReadHit->NumHits = 0;
			pFwdReadHit->LowHitInstances = 0;
			pRevReadHit->LowHitInstances = 0;
			if(pFwdReadHit->NAR == eNARAccepted)
				pFwdReadHit->NAR = eNARPENoHit;
			if(pRevReadHit->NAR == eNARAccepted)
				pRevReadHit->NAR = eNARPENoHit;
			PartnerUnpaired += 1;
			continue;
			}
		}

	// at least one end was uniquely aligning although not accepted as a PE
	PartnerUnpaired += 1;
	if(PEproc == ePEorphan || PEproc == ePEorphanSE)		// allowed to try and recover?
		{
		if(pFwdReadHit->NumHits == 1 && pRevReadHit->NAR != eNARNs)
			{
			if(AcceptThisChromID(pFwdReadHit->HitLoci.Hit.Seg[0].ChromID))
				{
				// first try using the 5' alignment as an anchor, if that doesn't provide a pair then will later try using the 3' as the anchor
				bPartnerPaired = false;
				pHitSeq = &pRevReadHit->Read[pRevReadHit->DescrLen+1];
				pSeq = ReadSeq;
				for(SeqIdx = 0; SeqIdx < pRevReadHit->ReadLen; SeqIdx++)
					*pSeq++ = *pHitSeq++ & 0x07;

				// AlignPartnerRead
				// Have been able to unquely align one read out of a pair, now need to align the other read
				// if not bPairStrand
				//		if PE1 was to sense strand then expect PE2 on the antisense strand downstream towards the 3' end of targeted chrom:		b3primeExtend=true,bAntisense=true
				//		if PE1 was to antisense strand then expect PE2 on the sense strand upstream towards the 5' end of targeted chrom:		b3primeExtend=false,bAntisense=false
				// if bPairStrand
				//		if PE1 was to sense strand then expect PE2 on the sense strand downstream towards the 3' end of targeted chrom:			b3primeExtend=true,bAntisense=false
				//		if PE1 was to antisense strand then expect PE2 on the antisense strand upstream towards the 5' end of targeted chrom:	b3primeExtend=false,bAntisense=true
				b3primeExtend = pFwdReadHit->HitLoci.Hit.Seg[0].Strand == '+' ? true : false;
				if(bPairStrand)
					bAntisense = pFwdReadHit->HitLoci.Hit.Seg[0].Strand == '+' ? false : true;
				else
					bAntisense = pFwdReadHit->HitLoci.Hit.Seg[0].Strand == '+' ? true : false;
				if(m_bPEcircularised)
					b3primeExtend = !b3primeExtend;

				OrphStartLoci = AdjStartLoci(&pFwdReadHit->HitLoci.Hit.Seg[0]);
				OrphEndLoci = AdjEndLoci(&pFwdReadHit->HitLoci.Hit.Seg[0]);

				Rslt = m_pSfxArray->AlignPairedRead(b3primeExtend,bAntisense,
							pFwdReadHit->HitLoci.Hit.Seg[0].ChromID,	  // accepted aligned read was on this chromosome
							OrphStartLoci,			// accepted aligned read started at this loci
							OrphEndLoci,			// and ending at this loci
							PairMinLen,			// expecting partner to align at least this distance away from accepted aligned read
							PairMaxLen,			// but no more than this distance away
							MaxSubs == 0 ? 0 : max(1,(pRevReadHit->ReadLen * MaxSubs)/100),		// any accepted alignment can have at most this many mismatches
							MinEditDist,					// and must be at least this Hamming away from the next best putative alignment
							pRevReadHit->ReadLen,		  // length of read excluding any eBaseEOS
							ReadSeq,	  // pts to 5' start of read sequence
							&HitLoci);	  // where to return any paired read alignment loci

				if(Rslt == 1)
					{
					SeqFragLen = PEInsertSize(PairMinLen,PairMaxLen,bPairStrand,pFwdReadHit->HitLoci.Hit.Seg[0].Strand,OrphStartLoci,OrphEndLoci,HitLoci.Seg[0].Strand,AdjStartLoci(&HitLoci.Seg[0]),AdjEndLoci(&HitLoci.Seg[0]));
					if(SeqFragLen <= 0)
						Rslt = 0;
					}

				if(Rslt == 1)	
					{
					// with 5' anchor was able to find an alignment within the min/max insert size and it is known chrom accepted
					pRevReadHit->HitLoci.Hit = HitLoci;
					pRevReadHit->NumHits = 1;
					pRevReadHit->LowMMCnt = HitLoci.Seg[0].Mismatches;
					pRevReadHit->LowHitInstances = 1;
					// this paired end alignment has been accepted
					pFwdReadHit->FlgPEAligned = 1;
					pRevReadHit->FlgPEAligned = 1;
					pFwdReadHit->NAR = eNARAccepted;
					pRevReadHit->NAR = eNARAccepted;
					LongestSeqFragLen = max(LongestSeqFragLen,SeqFragLen);
					m_pLenDist[SeqFragLen] += 1;
					AcceptedNumPaired += 1;
					PartnerPaired += 1;
					continue;	// try next pair
					}
				}
			else
				if(pFwdReadHit->NAR == eNARAccepted)
					{
					pFwdReadHit->NumHits = 0;
					pFwdReadHit->LowHitInstances = 0;
					pFwdReadHit->NAR = eNARChromFilt;
					}
			}


		if(pRevReadHit->NumHits == 1  && pFwdReadHit->NAR != eNARNs)
			{
			if(AcceptThisChromID(pRevReadHit->HitLoci.Hit.Seg[0].ChromID))
				{
				pHitSeq = &pFwdReadHit->Read[pFwdReadHit->DescrLen+1];
				pSeq = ReadSeq;
				for(SeqIdx = 0; SeqIdx < pFwdReadHit->ReadLen; SeqIdx++)
					*pSeq++ = *pHitSeq++ & 0x07;
				// AlignPartnerRead
				// Have been able to unquely align one read out of a pair, now need to align the other read
				// if not bPairStrand
				//		if PE2 was to sense strand then expect PE1 on the antisense strand downstream towards the 3' end of targeted chrom:		b3primeExtend=true,bAntisense=true
				//		if PE2 was to antisense strand then expect PE1 on the sense strand upstream towards the 5' end of targeted chrom:		b3primeExtend=false,bAntisense=false
				// if bPairStrand
				//		if PE2 was to sense strand then expect PE1 on the sense strand upstream towards the 5' end of targeted chrom:			b3primeExtend=false,bAntisense=false
				//		if PE2 was to antisense strand then expect PE2 on the antisense strand downstream towards the 3' end of targeted chrom:	b3primeExtend=true,bAntisense=true
				b3primeExtend = pRevReadHit->HitLoci.Hit.Seg[0].Strand == '+' ? true : false;
				bAntisense = pRevReadHit->HitLoci.Hit.Seg[0].Strand == '+' ? true : false;
				if(bPairStrand)
					{
					b3primeExtend = !b3primeExtend;
					bAntisense = !bAntisense;
					}
				if(m_bPEcircularised)
					b3primeExtend = !b3primeExtend;

				OrphStartLoci = AdjStartLoci(&pRevReadHit->HitLoci.Hit.Seg[0]);
				OrphEndLoci = AdjEndLoci(&pRevReadHit->HitLoci.Hit.Seg[0]);

				Rslt = m_pSfxArray->AlignPairedRead(b3primeExtend,bAntisense,
							pRevReadHit->HitLoci.Hit.Seg[0].ChromID,	  // accepted aligned read was on this chromosome
							OrphStartLoci,			// accepted aligned read started at this loci
							OrphEndLoci,		  // and ending at this loci
							PairMinLen,			// expecting partner to align at least this distance away from accepted aligned read
							PairMaxLen,		// but no more than this distance away
							MaxSubs == 0 ? 0 : max(1,(pFwdReadHit->ReadLen * MaxSubs)/100), // any accepted alignment can have at most this many mismatches
							MinEditDist,					// and must be at least this Hamming away from the next best putative alignment
							pFwdReadHit->ReadLen,			// length of read excluding any eBaseEOS
							ReadSeq,	  // pts to 5' start of read sequence
							&HitLoci);	  // where to return any paired read alignment loci

				if(Rslt == 1)
					{
					SeqFragLen = PEInsertSize(PairMinLen,PairMaxLen,bPairStrand,HitLoci.Seg[0].Strand,AdjStartLoci(&HitLoci.Seg[0]),AdjEndLoci(&HitLoci.Seg[0]),pRevReadHit->HitLoci.Hit.Seg[0].Strand,OrphStartLoci,OrphEndLoci);
					if(SeqFragLen <= 0)
						Rslt = 0;
					}

				if(Rslt == 1)
					{
					// with 3' anchor was able to find an alignment within the min/max insert size
					pFwdReadHit->HitLoci.Hit = HitLoci;
					pFwdReadHit->LowMMCnt = HitLoci.Seg[0].Mismatches;
					pFwdReadHit->NumHits = 1;
					pFwdReadHit->LowHitInstances = 1;
					// this paired end alignment has been accepted
					pFwdReadHit->FlgPEAligned = 1;
					pRevReadHit->FlgPEAligned = 1;
					pFwdReadHit->NAR = eNARAccepted;
					pRevReadHit->NAR = eNARAccepted;
					LongestSeqFragLen = max(LongestSeqFragLen,SeqFragLen);
					m_pLenDist[SeqFragLen] += 1;
					AcceptedNumPaired += 1;
					PartnerPaired += 1;
					continue;	// try next pair
					}
				}
			else
				if(pRevReadHit->NAR == eNARAccepted)
					{
					pFwdReadHit->NumHits = 0;
					pFwdReadHit->LowHitInstances = 0;
					pFwdReadHit->NAR = eNARChromFilt;
					}
			}
		}

	// unable to accept as PE
	if(pFwdReadHit->NAR == eNARChromFilt || pRevReadHit->NAR == eNARChromFilt)
		NumFilteredByChrom += 1;
	if(pFwdReadHit->NAR == eNARPEInsertMin || pRevReadHit->NAR == eNARPEInsertMin)
		UnderLenPairs += 1;
	if(pFwdReadHit->NAR == eNARPEInsertMax || pRevReadHit->NAR == eNARPEInsertMax)
		OverLenPairs += 1;

	if(!(PEproc == ePEorphanSE || PEproc == ePEuniqueSE))		
		{
		pFwdReadHit->NumHits = 0;
		pFwdReadHit->LowHitInstances = 0;
		pRevReadHit->NumHits = 0;
		pRevReadHit->LowHitInstances = 0;
		if(pFwdReadHit->NAR == eNARAccepted)
			pFwdReadHit->NAR = eNARPENoHit;
		if(pRevReadHit->NAR == eNARAccepted)
			pRevReadHit->NAR = eNARPENoHit;
		continue;
		}

	// allowed to accept as being SE if was able to uniquely align
	bool bChromFilt;
	if(pFwdReadHit->NumHits == 1)
		bChromFilt = AcceptThisChromID(pFwdReadHit->HitLoci.Hit.Seg[0].ChromID); 
	else 
		bChromFilt = false;
	if(pFwdReadHit->NumHits != 1 || !bChromFilt)
		{
		pFwdReadHit->NumHits = 0;
		pFwdReadHit->LowHitInstances = 0;
		if(pFwdReadHit->NAR == eNARAccepted)
			pFwdReadHit->NAR = bChromFilt ? eNARChromFilt : eNARPEUnalign;
		}
	else
		{
		pFwdReadHit->NAR = eNARAccepted;
		AcceptedNumSE += 1;
		}

	if(pRevReadHit->NumHits == 1)
		bChromFilt = AcceptThisChromID(pRevReadHit->HitLoci.Hit.Seg[0].ChromID);
	else 
		bChromFilt = false;

	if(pRevReadHit->NumHits != 1 || !bChromFilt)
		{
		pRevReadHit->NumHits = 0;
		pRevReadHit->LowHitInstances = 0;
		if(pRevReadHit->NAR == eNARAccepted)
			pRevReadHit->NAR = bChromFilt ? eNARChromFilt : eNARPEUnalign;
		}
	else
		{
		pRevReadHit->NAR = eNARAccepted;
		AcceptedNumSE += 1;
		}
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processed putative %u pairs, accepted %u",PairReadIdx/2,AcceptedNumPaired);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"From %d putative pairs there were %d accepted (%d from recovered orphans)",
					 PairReadIdx/2,AcceptedNumPaired,PartnerPaired);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"%d putative pairs unrecoverable as still orphan partnered",PartnerUnpaired - PartnerPaired);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"%d under length, %d over length pairs not accepted",UnderLenPairs,OverLenPairs);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"%d aligned pairs filtered out by chromosome",NumFilteredByChrom);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"%d putative pairs have neither end uniquely aligned",UnalignedPairs);
if(PEproc == ePEuniqueSE || PEproc == ePEorphanSE)
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"%d reads accepted as SE aligned",AcceptedNumSE);

if(gProcessingID > 0)
	{
	UINT32 NumPairsLoaded;
	NumPairsLoaded = m_NumReadsLoaded/2;
	gSQLiteSummaries.AddResult(gProcessingID,(char *)"PEProcessing",ePTUint32,sizeof(NumPairsLoaded),"Loaded",&NumPairsLoaded);
	gSQLiteSummaries.AddResult(gProcessingID,(char *)"PEProcessing",ePTUint32,sizeof(AcceptedNumPaired),"Accepted",&AcceptedNumPaired);
	gSQLiteSummaries.AddResult(gProcessingID,(char *)"PEProcessing",ePTUint32,sizeof(PartnerPaired),"RecoveredOrphans",&PartnerPaired);
	gSQLiteSummaries.AddResult(gProcessingID,(char *)"PEProcessing",ePTUint32,sizeof(UnderLenPairs),"UnderInsertSize",&UnderLenPairs);
	gSQLiteSummaries.AddResult(gProcessingID,(char *)"PEProcessing",ePTUint32,sizeof(OverLenPairs),"OverInsertSize",&OverLenPairs);
	gSQLiteSummaries.AddResult(gProcessingID,(char *)"PEProcessing",ePTUint32,sizeof(NumFilteredByChrom),"FilteredByChrom",&NumFilteredByChrom);
	gSQLiteSummaries.AddResult(gProcessingID,(char *)"PEProcessing",ePTUint32,sizeof(UnalignedPairs),"NonUniqueEnds",&UnalignedPairs);
	gSQLiteSummaries.AddResult(gProcessingID,(char *)"PEProcessing",ePTUint32,sizeof(AcceptedNumSE),"AcceptedNumSE",&AcceptedNumSE);
	}

if(m_hStatsFile != -1)
	{
	char szLineBuff[cMaxReadLen*2];
	int BuffIdx = 0;
	for(PairReadIdx = 0; PairReadIdx <= (UINT32)LongestSeqFragLen; PairReadIdx++)
		{
		BuffIdx += sprintf(&szLineBuff[BuffIdx],"%d,%d\n",PairReadIdx,m_pLenDist[PairReadIdx]);
		if(BuffIdx + 100 > sizeof(szLineBuff))
			{
			CUtility::SafeWrite(m_hStatsFile,szLineBuff,BuffIdx);
			BuffIdx = 0;
			}
		}
	if(BuffIdx)
		{
		CUtility::SafeWrite(m_hStatsFile,szLineBuff,BuffIdx);
		BuffIdx = 0;
		}
	}
if((Rslt=SortReadHits(eRSMHitMatch,false,true)) < eBSFSuccess)
	{
	Reset(false);
	return(Rslt);
	}
return(eBSFSuccess);
}


int
CAligner::ReportAlignStats(void)		// report basic alignment statistics
{
int Rslt;
char szChromName[128];
int NumUniques = 0;
int NumPlusHits = 0;
int NumNoMatches = 0;
int NumChimeric = 0;
int NumMultiMatches = 0;
int NumHamming = 0;
tsReadHit *pReadHit = NULL;
bool bSimReads = false;
UINT32 NumReads1EdgeAligned = 0;
UINT32 NumReads2EdgeAligned = 0;
UINT32 NumReadsMisaligned = 0;
UINT32 NumIndels = 0;
UINT32 NumSpliced = 0;
UINT32 PrevTargEntry = 0;
UINT32 NumTrimmed = 0;



UINT32 NARUnaligned = 0;				// read has yet to be aligned
UINT32 NARAccepted = 0;					// read has been accepted as being aligned
UINT32 NARNs = 0;						// not accepted because contains too many indeterminate bases
UINT32 NARNoHit = 0;					// not accepted as aligned because unable to find any potential hits
UINT32 NARMMDelta = 0;					// not accepted as aligned because MMDelta criteria not met
UINT32 NARMultiAlign = 0;				// not accepted as aligned because was multiloci aligned
UINT32 NARTrim = 0;						// not accepted as aligned because aligned read excessively trimmed
UINT32 NARSpliceJctn = 0;				// not accepted as aligned because aligned read orphan splice junction
UINT32 NARmicroInDel = 0;				// not accepted as aligned because aligned read orphan microInDel
UINT32 NARPCRdup = 0;					// not accepted as aligned because aligned read was a PCR duplicate
UINT32 NARNonUnique = 0;				// not accepted as aligned because aligned read was a duplicate sequence
UINT32 NARChromFilt = 0;				// not accepted as aligned because aligned to a filtered chrom
UINT32 NARRegionFilt = 0;				// not accepted as aligned because not aligned to a priority region
UINT32 NARPEInsertMin = 0;				// not accepted as aligned because both PE ends align but less than minimum insert size
UINT32 NARPEInsertMax = 0;				// not accepted as aligned because both PE ends align but more than maximum insert size
UINT32 NARPENoHit = 0;					// not accepted as aligned because PE aligning and although able to align this read was unable to align partner end
UINT32 NARPEStrand = 0;					// not accepted as aligned because PE aligning and although able to align this read other read was aligned to inconsistent strand
UINT32 NARPEChrom = 0;					// not accepted as aligned because PE aligning and although able to align this read other read was aligned to different chromosome
UINT32 NARPEUnalign = 0;				// not accepted as aligned because PE aligning and unable to accept this alignment
UINT32 NARLociConstrained = 0;			// not accepted as aligned because alignment violated loci base constraints
if((Rslt=SortReadHits(eRSMHitMatch,false)) < eBSFSuccess)
	{
	Reset(false);
	return(Rslt);
	}
while((pReadHit = IterReads(pReadHit))!=NULL)
	{
	switch(pReadHit->NAR) {
		case eNARUnaligned:				// read has yet to be aligned
			NARUnaligned += 1;
			break;

		case eNARAccepted:									// accepted hit, accumulate strand count
			NARAccepted += 1;
			NumUniques += 1;
			if(pReadHit->HitLoci.Hit.Seg[0].Strand == '+')
				NumPlusHits += 1;
			if(pReadHit->HitLoci.Hit.FlgInDel)
				NumIndels += 1;
			if(pReadHit->HitLoci.Hit.FlgChimeric)
				NumChimeric += 1;
			if(pReadHit->HitLoci.Hit.FlgSplice)
				NumSpliced += 1;
			if(pReadHit->HitLoci.FlagTR)
				NumTrimmed += 1;
			pReadHit->HitLoci.FlagIA = 0;
			pReadHit->HitLoci.FlagCA = 0;
			if(bSimReads || NumUniques == 1)				  // check if simulated reads containing loci
				{
				bool bRandReadSeq;
				int x1;
				char szxChrom[100];
				char szx1Chrom[100];
				char szx2Chrom[100];
				char szSimType[100];
				int xStart;
				int xEnd;
				int xLen;
				char xStrand;
				int xErrs;
				int Its;
				bRandReadSeq = false;
				// >lcl|usimreads|00000001|gnl|UG|Ta#S58887126|330|429|100|+|0|0|0
				// >lcl|usimreads|00000001  
				// |gnl|UG|Ta#S58887126
				// |330|429|100|+|0|0|0
				Its = sscanf((char *)pReadHit->Read,"%[^|]|usimreads|%d|%[^|]|%d|%d|%d|%c|%d",szSimType,&x1,szxChrom,&xStart,&xEnd,&xLen,&xStrand,&xErrs);
				if(Its >= 6)
					{
					bSimReads = true;
					if(!stricmp(szSimType,"lcl"))
						bRandReadSeq = false;
					else
						bRandReadSeq = true;
					}
				else
					{
					Its = sscanf((char *)pReadHit->Read,"%[^|]|usimreads|%d|%[^|]|%[^|]|%[^|]|%d|%d|%d|%c|%d",szSimType,&x1,szxChrom,szx1Chrom,szx2Chrom,&xStart,&xEnd,&xLen,&xStrand,&xErrs);
					if(Its < 8) 
						bSimReads = false;
					else
						{
						strcat(szxChrom,"|");
						strcat(szxChrom,szx1Chrom);
						strcat(szxChrom,"|");
						strcat(szxChrom,szx2Chrom);
						bSimReads = true;
						if(!stricmp(szSimType,"lcl"))
							bRandReadSeq = false;
						else
							bRandReadSeq = true;
						}
					}

				if(pReadHit->HitLoci.Hit.Seg[0].ChromID != (UINT32)PrevTargEntry)
					{
					m_pSfxArray->GetIdentName(pReadHit->HitLoci.Hit.Seg[0].ChromID,sizeof(szChromName),szChromName);
					PrevTargEntry = pReadHit->HitLoci.Hit.Seg[0].ChromID;
					}

				if(!bSimReads)
					break;

				if(stricmp(szChromName,szxChrom))
					{
					NumReadsMisaligned += 1;
					pReadHit->HitLoci.FlagIA = 1;
					break;
					}

				int LeftOfs;
				int RightOfs;

				// the determination as to if a read has been correctly aligned is a little crude but will generally be correct
				// the crudeness arises because a read may not have been segmented during the alignment process although it originated from a microInDel or
				// splice junction so only one of the edges will be correct. This could have been because there were insufficient mismatches in a flank
				// to have been worth the cost of exploring as a potential segmented alignment.
				// Another complication is auto-edge trimming resulting in alignments with edges not coincident with the source simulated read loci.
				// The flow here is to classify the alignments into three groups:
				// a)	both edges are coincident with the source simulated reads edge loci
				// b)   one edge only is coincident
				// c)   neither edge is coincident
				LeftOfs = (int32)pReadHit->HitLoci.Hit.Seg[0].MatchLoci;
				if(pReadHit->HitLoci.FlagSegs)
					RightOfs = (int32)pReadHit->HitLoci.Hit.Seg[1].MatchLoci + pReadHit->HitLoci.Hit.Seg[1].MatchLen - 1;
				else
					RightOfs = LeftOfs + pReadHit->HitLoci.Hit.Seg[0].MatchLen - 1;

				if(LeftOfs == xStart || RightOfs == xEnd)
					{
					if(LeftOfs != xStart || RightOfs != xEnd)
						NumReads1EdgeAligned += 1;
					else
						NumReads2EdgeAligned += 1;
					pReadHit->HitLoci.FlagCA = 1;
					}
				else
					{
					NumReadsMisaligned += 1;
					pReadHit->HitLoci.FlagIA = 1;
					}
				}
			break;

		case eNARNs:					// not accepted because contains too many indeterminate bases
//			NARNs += 1;
//			NumNoMatches += 1;
			break;
		case eNARNoHit:					// not accepted as aligned because unable to find any potential hits
			NARNoHit += 1;
			NumNoMatches += 1;
			break;
		case eNARMMDelta:				// not accepted as aligned because MMDelta (Hamming) criteria not met
			NARMMDelta += 1;
			NumHamming += 1;
			break;
		case eNARMultiAlign:			// not accepted as aligned because was multiloci aligned
			NARMultiAlign += 1;
			NumMultiMatches+= 1;
			break;
		case eNARTrim:					// not accepted as aligned because aligned read excessively trimmed
			NARTrim += 1;
			break;
		case eNARSpliceJctn:			// not accepted as aligned because aligned read orphan splice junction
			NARSpliceJctn += 1;
			break;
		case eNARmicroInDel:			// not accepted as aligned because aligned read orphan microInDel
			NARmicroInDel += 1;
			break;
		case eNARPCRdup:				// not accepted as aligned because aligned read was a PCR duplicate
			NARPCRdup += 1;
			break;
		case eNARNonUnique:				// not accepted as aligned because aligned read was a duplicate sequence
			NARNonUnique += 1;
			break;
		case eNARChromFilt:				// not accepted as aligned because aligned to a filtered chrom
			NARChromFilt += 1;
			break;
		case eNARRegionFilt:			// not accepted as aligned because not aligned to a priority region
			NARRegionFilt += 1;
			break;
		case eNARPEInsertMin:			// not accepted as aligned because both PE ends align but less than minimum insert size
			NARPEInsertMin += 1;
			break;
		case eNARPEInsertMax:			// not accepted as aligned because both PE ends align but more than maximum insert size
			NARPEInsertMax += 1;
			break;
		case eNARPENoHit:				// not accepted as aligned because PE aligning and although able to align this read was unable to align partner end
			NARPENoHit += 1;
			break;
		case eNARPEStrand:				// not accepted as aligned because PE aligning and although able to align this read other read was aligned to inconsistent strand
			NARPEStrand += 1;
			break;
		case eNARPEChrom:				// not accepted as aligned because PE aligning and although able to align this read other read was aligned to different chromosome
			NARPEChrom += 1;
			break;
		case eNARPEUnalign:				// not accepted as aligned because PE aligning and unable to accept this alignment
			NARPEUnalign += 1;
			break;
		case eNARLociConstrained:		// not accepted as aligned because alignment violated loci base constraints
			NARLociConstrained += 1;
			break;
		default:							
			break;
		}
	}

if(m_MLMode >= eMLall)
	NumNoMatches = m_TotNonAligned;
NumNoMatches += m_NumSloughedNs;
NARNs = m_NumSloughedNs;
gDiagnostics.DiagOut(eDLInfo,gszProcName,"From %d source reads there are %d accepted alignments, %d on '+' strand, %d on '-' strand", m_OrigNumReadsLoaded,NumUniques,NumPlusHits,NumUniques-NumPlusHits);
if(bSimReads)
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"There are %d (%d 2 edge, %d 1 edge) high confidence aligned simulated reads with %d misaligned",NumReads2EdgeAligned + NumReads1EdgeAligned,NumReads2EdgeAligned,NumReads1EdgeAligned,NumReadsMisaligned);
if(NumIndels)
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Of the accepted aligned reads, %d contained microInDels", NumIndels);
if(NumSpliced)
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Of the accepted aligned reads, %d contained splice junctions", NumSpliced);
if(NumChimeric)
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Of the accepted aligned reads, %d were chimeric", NumChimeric);

gDiagnostics.DiagOut(eDLInfo,gszProcName,"A further %d multiloci aligned reads could not accepted as hits because they were unresolvable",NumMultiMatches);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"A further %d aligned reads were not accepted as hits because of insufficient Hamming edit distance",NumHamming);

gDiagnostics.DiagOut(eDLInfo,gszProcName,"A further %d '+' and %d '-' strand aligned reads not accepted because of flank trimming (%d were trimmed) requirements",m_ElimPlusTrimed,m_ElimMinusTrimed,NumTrimmed);

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Unable to align %d source reads of which %d were not aligned as they contained excessive number of indeterminate 'N' bases",NumNoMatches,NARNs);


gDiagnostics.DiagOut(eDLInfo,gszProcName,"Read nonalignment reason summary:");
gDiagnostics.DiagOut(eDLInfo,gszProcName,"   %u (%s) %s",NARUnaligned,m_NARdesc[eNARUnaligned].pszNAR, m_NARdesc[eNARUnaligned].pszNARdescr);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"   %u (%s) %s",NARAccepted,m_NARdesc[eNARAccepted].pszNAR, m_NARdesc[eNARAccepted].pszNARdescr);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"   %u (%s) %s",NARNs,m_NARdesc[eNARNs].pszNAR, m_NARdesc[eNARNs].pszNARdescr);

if(m_MLMode >= eMLall)
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"   %u (%s) %s",NumNoMatches,m_NARdesc[eNARNoHit].pszNAR, m_NARdesc[eNARNoHit].pszNARdescr);
else
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"   %u (%s) %s",NARNoHit,m_NARdesc[eNARNoHit].pszNAR, m_NARdesc[eNARNoHit].pszNARdescr);

gDiagnostics.DiagOut(eDLInfo,gszProcName,"   %u (%s) %s",NARMMDelta,m_NARdesc[eNARMMDelta].pszNAR, m_NARdesc[eNARMMDelta].pszNARdescr);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"   %u (%s) %s",NARMultiAlign,m_NARdesc[eNARMultiAlign].pszNAR, m_NARdesc[eNARMultiAlign].pszNARdescr);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"   %u (%s) %s",NARTrim,m_NARdesc[eNARTrim].pszNAR, m_NARdesc[eNARTrim].pszNARdescr);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"   %u (%s) %s",NARSpliceJctn,m_NARdesc[eNARSpliceJctn].pszNAR, m_NARdesc[eNARSpliceJctn].pszNARdescr);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"   %u (%s) %s",NARmicroInDel,m_NARdesc[eNARmicroInDel].pszNAR, m_NARdesc[eNARmicroInDel].pszNARdescr);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"   %u (%s) %s",NARPCRdup,m_NARdesc[eNARPCRdup].pszNAR, m_NARdesc[eNARPCRdup].pszNARdescr);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"   %u (%s) %s",NARNonUnique,m_NARdesc[eNARNonUnique].pszNAR, m_NARdesc[eNARNonUnique].pszNARdescr);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"   %u (%s) %s",NARChromFilt,m_NARdesc[eNARChromFilt].pszNAR, m_NARdesc[eNARChromFilt].pszNARdescr);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"   %u (%s) %s",NARRegionFilt,m_NARdesc[eNARRegionFilt].pszNAR, m_NARdesc[eNARRegionFilt].pszNARdescr);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"   %u (%s) %s",NARPEInsertMin,m_NARdesc[eNARPEInsertMin].pszNAR, m_NARdesc[eNARPEInsertMin].pszNARdescr);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"   %u (%s) %s",NARPEInsertMax,m_NARdesc[eNARPEInsertMax].pszNAR, m_NARdesc[eNARPEInsertMax].pszNARdescr);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"   %u (%s) %s",NARPENoHit,m_NARdesc[eNARPENoHit].pszNAR, m_NARdesc[eNARPENoHit].pszNARdescr);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"   %u (%s) %s",NARPEStrand,m_NARdesc[eNARPEStrand].pszNAR, m_NARdesc[eNARPEStrand].pszNARdescr);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"   %u (%s) %s",NARPEChrom,m_NARdesc[eNARPEChrom].pszNAR, m_NARdesc[eNARPEChrom].pszNARdescr);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"   %u (%s) %s",NARPEUnalign,m_NARdesc[eNARPEUnalign].pszNAR, m_NARdesc[eNARPEUnalign].pszNARdescr);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"   %u (%s) %s",NARLociConstrained,m_NARdesc[eNARLociConstrained].pszNAR, m_NARdesc[eNARLociConstrained].pszNARdescr);

if(gProcessingID > 0)
	{
	UINT32 Cricks;
	Cricks = NumUniques-NumPlusHits;
	if(m_MLMode >= eMLall)
		gSQLiteSummaries.AddResult(gProcessingID,(char *)"Alignments",ePTUint32,sizeof(m_OrigNumReadsLoaded),"NumLoaded",&m_OrigNumReadsLoaded);
	else
		gSQLiteSummaries.AddResult(gProcessingID,(char *)"Alignments",ePTUint32,sizeof(m_NumReadsLoaded),"NumLoaded",&m_NumReadsLoaded);
	gSQLiteSummaries.AddResult(gProcessingID,(char *)"Alignments",ePTUint32,sizeof(NumUniques),"AcceptedAligned",&NumUniques);
	gSQLiteSummaries.AddResult(gProcessingID,(char *)"Alignments",ePTUint32,sizeof(NumPlusHits),"AcceptedSense",&NumPlusHits);
	gSQLiteSummaries.AddResult(gProcessingID,(char *)"Alignments",ePTUint32,sizeof(Cricks),"AcceptedAntisense",&Cricks);
	gSQLiteSummaries.AddResult(gProcessingID,(char *)"Alignments",ePTUint32,sizeof(NumChimeric),"AcceptedChimeric",&NumChimeric);
	gSQLiteSummaries.AddResult(gProcessingID,(char *)"Alignments",ePTUint32,sizeof(NumIndels),"AcceptedIndels",&NumIndels);
	gSQLiteSummaries.AddResult(gProcessingID,(char *)"Alignments",ePTUint32,sizeof(NumSpliced),"AcceptedSpliced",&NumSpliced);
	gSQLiteSummaries.AddResult(gProcessingID,(char *)"Alignments",ePTUint32,sizeof(NumMultiMatches),"RjctMultiMatches",&NumMultiMatches);
	gSQLiteSummaries.AddResult(gProcessingID,(char *)"Alignments",ePTUint32,sizeof(NumHamming),"RjctNumHamming",&NumHamming);
	gSQLiteSummaries.AddResult(gProcessingID,(char *)"Alignments",ePTUint32,sizeof(m_ElimPlusTrimed),"RjctPlusTrim",&m_ElimPlusTrimed);
	gSQLiteSummaries.AddResult(gProcessingID,(char *)"Alignments",ePTUint32,sizeof(m_ElimMinusTrimed),"RjctMinusTrim",&m_ElimMinusTrimed);
	gSQLiteSummaries.AddResult(gProcessingID,(char *)"Alignments",ePTUint32,sizeof(NumTrimmed),"NumTrimmed",&NumTrimmed);
	gSQLiteSummaries.AddResult(gProcessingID,(char *)"Alignments",ePTUint32,sizeof(NumNoMatches),"Unalignable",&NumNoMatches);
	gSQLiteSummaries.AddResult(gProcessingID,(char *)"Alignments",ePTUint32,sizeof(m_NumSloughedNs),"RjctExcessNs",&m_NumSloughedNs);


	gSQLiteSummaries.AddResult(gProcessingID,(char *)"Alignments",ePTUint32,sizeof(NARUnaligned),"NARUnaligned",&NARUnaligned);
	gSQLiteSummaries.AddResult(gProcessingID,(char *)"Alignments",ePTUint32,sizeof(NARAccepted),"NARAccepted",&NARAccepted);
	gSQLiteSummaries.AddResult(gProcessingID,(char *)"Alignments",ePTUint32,sizeof(NARNs),"NARNs",&NARNs);

	if(m_MLMode >= eMLall)
		gSQLiteSummaries.AddResult(gProcessingID,(char *)"Alignments",ePTUint32,sizeof(NumNoMatches),"NARNoHit",&NumNoMatches);
	else
		gSQLiteSummaries.AddResult(gProcessingID,(char *)"Alignments",ePTUint32,sizeof(NARNoHit),"NARNoHit",&NARNoHit);
	
	gSQLiteSummaries.AddResult(gProcessingID,(char *)"Alignments",ePTUint32,sizeof(NARMMDelta),"NARMMDelta",&NARMMDelta);
	gSQLiteSummaries.AddResult(gProcessingID,(char *)"Alignments",ePTUint32,sizeof(NARMultiAlign),"NARMultiAlign",&NARMultiAlign);
	gSQLiteSummaries.AddResult(gProcessingID,(char *)"Alignments",ePTUint32,sizeof(NARTrim),"NARTrim",&NARTrim);
	gSQLiteSummaries.AddResult(gProcessingID,(char *)"Alignments",ePTUint32,sizeof(NARSpliceJctn),"NARSpliceJctn",&NARSpliceJctn);
	gSQLiteSummaries.AddResult(gProcessingID,(char *)"Alignments",ePTUint32,sizeof(NARmicroInDel),"NARmicroInDel",&NARmicroInDel);
	gSQLiteSummaries.AddResult(gProcessingID,(char *)"Alignments",ePTUint32,sizeof(NARPCRdup),"NARPCRdup",&NARPCRdup);
	gSQLiteSummaries.AddResult(gProcessingID,(char *)"Alignments",ePTUint32,sizeof(NARNonUnique),"NARNonUnique",&NARNonUnique);
	gSQLiteSummaries.AddResult(gProcessingID,(char *)"Alignments",ePTUint32,sizeof(NARChromFilt),"NARChromFilt",&NARChromFilt);
	gSQLiteSummaries.AddResult(gProcessingID,(char *)"Alignments",ePTUint32,sizeof(NARRegionFilt),"NARRegionFilt",&NARRegionFilt);
	gSQLiteSummaries.AddResult(gProcessingID,(char *)"Alignments",ePTUint32,sizeof(NARPEInsertMin),"NARPEInsertMin",&NARPEInsertMin);
	gSQLiteSummaries.AddResult(gProcessingID,(char *)"Alignments",ePTUint32,sizeof(NARPEInsertMax),"NARPEInsertMax",&NARPEInsertMax);
	gSQLiteSummaries.AddResult(gProcessingID,(char *)"Alignments",ePTUint32,sizeof(NARPENoHit),"NARPENoHit",&NARPENoHit);
	gSQLiteSummaries.AddResult(gProcessingID,(char *)"Alignments",ePTUint32,sizeof(NARPEChrom),"NARPEChrom",&NARPEChrom);
	gSQLiteSummaries.AddResult(gProcessingID,(char *)"Alignments",ePTUint32,sizeof(NARPEUnalign),"NARPEUnalign",&NARPEUnalign);
	gSQLiteSummaries.AddResult(gProcessingID,(char *)"Alignments",ePTUint32,sizeof(NARLociConstrained),"NARLociConstrained",&NARLociConstrained);
	}

m_NARAccepted = NARAccepted;		// further downstream processing may be interested in number of reads being reported as accepted
return(eBSFSuccess);
}


int
CAligner::ReportNoneAligned(void)
{
int ReadLen;
int SeqOfs;
int	NxtFastaCol;
int NumCols;
tsReadHit *pReadHit;
UINT32 SeqIdx;
etSeqBase Sequence[cMaxFastQSeqLen+1];	// to hold sequence (sans quality scores) for current read
UINT8 *pSeqVal;
etSeqBase *pSeq;

int NumNoMatches = 0;
int NumMultiMatches = 0;
int NumHamming = 0;

int LineLen = 0;
// user interested in the non-alignable reads?
// these only include those which had no alignment at all
if(m_hNoneAlignFile != -1 || m_gzNoneAlignFile != NULL)
	{
	SortReadHits(eRSMHitMatch,false);
	pReadHit = NULL;
	while((pReadHit = IterSortedReads(pReadHit))!=NULL)
		{
		if(!(pReadHit->NAR == eNARNs || pReadHit->NAR == eNARNoHit))
			continue;
		ReadLen = pReadHit->ReadLen;
		pSeqVal = &pReadHit->Read[pReadHit->DescrLen+1];
		pSeq = Sequence;
		for(SeqIdx = 0; SeqIdx < pReadHit->ReadLen; SeqIdx++,pSeq++,pSeqVal++)
			*pSeq = (*pSeqVal & 0x07);

		LineLen += sprintf(&m_pszLineBuff[LineLen],">lcl|nonealign %s %d|%d|%d\n",pReadHit->Read,pReadHit->ReadID,pReadHit->NumReads,pReadHit->ReadLen);
		if((LineLen + (2 * cMaxFastQSeqLen)) > cAllocLineBuffSize)
			{
			if(m_hNoneAlignFile != -1)
				CUtility::SafeWrite(m_hNoneAlignFile,m_pszLineBuff,LineLen);
			else
				CUtility::SafeWrite_gz(m_gzNoneAlignFile,m_pszLineBuff,LineLen);
			LineLen = 0;
			}
		SeqOfs = 0;
		NxtFastaCol = 0;
		while(ReadLen)
			{
			NumCols = ReadLen > 70 ? 70 : ReadLen;
			if((NumCols + NxtFastaCol) > 70)
				NumCols = 70 - NxtFastaCol;
			CSeqTrans::MapSeq2Ascii(&Sequence[SeqOfs],NumCols,&m_pszLineBuff[LineLen]);
			LineLen += NumCols;
			NxtFastaCol += NumCols;
			SeqOfs += NumCols;
			ReadLen -= NumCols;
			if(!ReadLen || NxtFastaCol >= 70)
				{
				LineLen += sprintf(&m_pszLineBuff[LineLen],"\n");
				NxtFastaCol = 0;
				}

			if((LineLen + (2 * cMaxFastQSeqLen)) > cAllocLineBuffSize)
				{
				if(m_hNoneAlignFile != -1)
					CUtility::SafeWrite(m_hNoneAlignFile,m_pszLineBuff,LineLen);
				else
					CUtility::SafeWrite_gz(m_gzNoneAlignFile,m_pszLineBuff,LineLen);
				LineLen = 0;
				}
			}
		}
	if(LineLen)
		{
		if(m_hNoneAlignFile != -1)
			CUtility::SafeWrite(m_hNoneAlignFile,m_pszLineBuff,LineLen);
		else
			CUtility::SafeWrite_gz(m_gzNoneAlignFile,m_pszLineBuff,LineLen);
		}

	if(m_hNoneAlignFile != -1)
		{
#ifdef _WIN32
		_commit(m_hNoneAlignFile);
#else
		fsync(m_hNoneAlignFile);
#endif
		close(m_hNoneAlignFile);
		m_hNoneAlignFile = -1;
		}
	else
		{
		gzclose(m_gzNoneAlignFile);
		m_gzNoneAlignFile = NULL;
		}
	}
return(eBSFSuccess);
}

int
CAligner::ReportMultiAlign(void)
{
int ReadLen;
int SeqOfs;
int	NxtFastaCol;
int NumCols;
tsReadHit *pReadHit;
UINT32 SeqIdx;
etSeqBase Sequence[cMaxFastQSeqLen+1];	// to hold sequence (sans quality scores) for current read
UINT8 *pSeqVal;
etSeqBase *pSeq;
int LineLen = 0;

if(m_hMultiAlignFile != -1 || m_gzMultiAlignFile != NULL)
	{
	SortReadHits(eRSMHitMatch,false);
	pReadHit = NULL;
	while((pReadHit = IterSortedReads(pReadHit))!=NULL)
		{
		if(pReadHit->NAR != eNARMultiAlign)
			continue;

		ReadLen = pReadHit->ReadLen;
		pSeqVal = &pReadHit->Read[pReadHit->DescrLen+1];
		pSeq = Sequence;
		for(SeqIdx = 0; SeqIdx < pReadHit->ReadLen; SeqIdx++,pSeq++,pSeqVal++)
			*pSeq = (*pSeqVal & 0x07);

		LineLen += sprintf(&m_pszLineBuff[LineLen],">lcl|multialign %s %d|%d|%d\n",pReadHit->Read,pReadHit->ReadID,pReadHit->NumReads,pReadHit->ReadLen);
		if((LineLen + (2 * cMaxFastQSeqLen)) > cAllocLineBuffSize)
			{
			if(m_hMultiAlignFile != -1)
				CUtility::SafeWrite(m_hMultiAlignFile,m_pszLineBuff,LineLen);
			else
				CUtility::SafeWrite_gz(m_gzMultiAlignFile,m_pszLineBuff,LineLen);
			LineLen = 0;
			}
		SeqOfs = 0;
		NxtFastaCol = 0;
		while(ReadLen)
			{
			NumCols = ReadLen > 70 ? 70 : ReadLen;
			if((NumCols + NxtFastaCol) > 70)
				NumCols = 70 - NxtFastaCol;
			CSeqTrans::MapSeq2Ascii(&Sequence[SeqOfs],NumCols,&m_pszLineBuff[LineLen]);
			LineLen += NumCols;
			NxtFastaCol += NumCols;
			SeqOfs += NumCols;
			ReadLen -= NumCols;
			if(!ReadLen || NxtFastaCol >= 70)
				{
				LineLen += sprintf(&m_pszLineBuff[LineLen],"\n");
				NxtFastaCol = 0;
				}

			if((LineLen + (2 * cMaxFastQSeqLen)) > cAllocLineBuffSize)
				{
				if(m_hMultiAlignFile != -1)
					CUtility::SafeWrite(m_hMultiAlignFile,m_pszLineBuff,LineLen);
				else
					CUtility::SafeWrite_gz(m_gzMultiAlignFile,m_pszLineBuff,LineLen);
				LineLen = 0;
				}
			}
		}
	if(LineLen)
		{
		if(m_hMultiAlignFile != -1)
			CUtility::SafeWrite(m_hMultiAlignFile,m_pszLineBuff,LineLen);
		else
			CUtility::SafeWrite_gz(m_gzMultiAlignFile,m_pszLineBuff,LineLen);
		}

	if(m_hMultiAlignFile != -1)
		{
#ifdef _WIN32
		_commit(m_hMultiAlignFile);
#else
		fsync(m_hMultiAlignFile);
#endif
		close(m_hMultiAlignFile);
		m_hMultiAlignFile = -1;
		}
	else
		{
		gzclose(m_gzMultiAlignFile);
		m_gzMultiAlignFile = NULL;
		}
	}
return(eBSFSuccess);
}



int
CAligner::FiltByChroms(void)
{
tBSFEntryID PrevTargEntry;
tBSFEntryID ExcludeTargEntryID;
int IncChromIdx;
int ExclChromIdx;
char szChromName[128];
bool bProcChrom = false;
int MatchesFiltOut = 0;
tsReadHit *pReadHit;

#ifdef _WIN32
RegexpMatch mc;
#else
regmatch_t mc;
int RegErr;					// regular expression parsing error
char szRegErr[128];			// to hold RegErr as textual representation ==> regerror();
#endif

SortReadHits(eRSMHitMatch,false);

if(m_NumExcludeChroms || m_NumIncludeChroms)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Now filtering matches by chromosome");
	PrevTargEntry = 0;
	pReadHit = NULL;
	while((pReadHit = IterSortedReads(pReadHit))!=NULL)
		{
		if(pReadHit->NAR == eNARAccepted)
			{
			if(pReadHit->HitLoci.Hit.Seg[0].ChromID != (UINT32)PrevTargEntry)
				{
				m_pSfxArray->GetIdentName(pReadHit->HitLoci.Hit.Seg[0].ChromID,sizeof(szChromName),szChromName);
				PrevTargEntry = pReadHit->HitLoci.Hit.Seg[0].ChromID;
				ExcludeTargEntryID = (tBSFEntryID)-1;
				}
			else
				{
				if(pReadHit->HitLoci.Hit.Seg[0].ChromID == (UINT32)ExcludeTargEntryID)
					{
					if(pReadHit->NAR == eNARAccepted)
						MatchesFiltOut += 1;
					pReadHit->NAR = eNARChromFilt;
					pReadHit->NumHits = 0;
					pReadHit->LowHitInstances = 0;
					}
				continue;
				}

			// to be included?
			bProcChrom = false;
			for(IncChromIdx = 0; IncChromIdx < m_NumIncludeChroms; IncChromIdx++)
				{
#ifdef _WIN32
				if(m_IncludeChromsRE[IncChromIdx]->Match(szChromName,&mc))
#else
				if(!regexec(&m_IncludeChromsRE[IncChromIdx],szChromName,1,&mc,0))
#endif
					{
					bProcChrom = true;
					break;
					}
				}

			// if not explicitly included then check if to be excluded?
			if(!bProcChrom && !m_NumIncludeChroms)
				{
				bProcChrom = true;
				for(ExclChromIdx = 0; ExclChromIdx < m_NumExcludeChroms; ExclChromIdx++)
					{
#ifdef _WIN32
					if(m_ExcludeChromsRE[ExclChromIdx]->Match(szChromName,&mc))
#else
					if(!regexec(&m_ExcludeChromsRE[ExclChromIdx],szChromName,1,&mc,0))
#endif
						{
						bProcChrom = false;
						break;
						}
					}
				}

			if(!bProcChrom)
				{
				ExcludeTargEntryID = pReadHit->HitLoci.Hit.Seg[0].ChromID;
				if(pReadHit->NAR  == eNARAccepted)
					MatchesFiltOut += 1;
				pReadHit->NAR = eNARChromFilt;
				pReadHit->NumHits = 0;
				pReadHit->LowHitInstances = 0;
				}
			else
				ExcludeTargEntryID = -1;
			}
		}
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Filtering by chromosome completed - removed %d  matches",MatchesFiltOut);

	if(MatchesFiltOut)
		SortReadHits(eRSMHitMatch,false,true);
	if(gProcessingID > 0)
		gSQLiteSummaries.AddResult(gProcessingID,(char *)"Filtered",ePTInt32,sizeof(MatchesFiltOut),"Chroms",&MatchesFiltOut);
		
	}
return(eBSFSuccess);
}

// remove hits not aligned into a priority region
int
CAligner::FiltByPriorityRegions(void) // remove hits not aligned into a priority region
{
char szPriorityChromName[100];
int PriorityChromID;

UINT32 MatchesFiltOut = 0;
UINT32 MatchesFiltIn = 0;

tsReadHit *pReadHit;
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Now filtering matches by prioritorised regions");
SortReadHits(eRSMHitMatch,false);
pReadHit = NULL;
while((pReadHit = IterSortedReads(pReadHit))!=NULL)
	{
	if(pReadHit->NAR == eNARAccepted)
		{
		// check if hit loci within region designated as being a priority exact matching region
		if(m_pSfxArray->GetIdentName(pReadHit->HitLoci.Hit.Seg[0].ChromID,sizeof(szPriorityChromName),szPriorityChromName)!=eBSFSuccess)
			{
			MatchesFiltOut += 1;
			pReadHit->NAR = eNARRegionFilt;
			pReadHit->NumHits = 0;
			pReadHit->HitLoci.Hit.Seg[0].ChromID = 0;
			pReadHit->HitLoci.Hit.Seg[1].ChromID = 0;
			continue;
			}
		if((PriorityChromID = m_pPriorityRegionBED->LocateChromIDbyName(szPriorityChromName)) < 1)
			{
			MatchesFiltOut += 1;
			pReadHit->NAR = eNARRegionFilt;
			pReadHit->NumHits = 0;
			pReadHit->HitLoci.Hit.Seg[0].ChromID = 0;
			pReadHit->HitLoci.Hit.Seg[1].ChromID = 0;
			continue;
			}

		if(!m_pPriorityRegionBED->InAnyFeature(PriorityChromID,(int)pReadHit->HitLoci.Hit.Seg[0].MatchLoci,
										(int)(pReadHit->HitLoci.Hit.Seg[0].MatchLoci + pReadHit->HitLoci.Hit.Seg[0].MatchLen-1)))
			{
			MatchesFiltOut += 1;
			pReadHit->NAR = eNARRegionFilt;
			pReadHit->NumHits = 0;
			pReadHit->HitLoci.Hit.Seg[0].ChromID = 0;
			pReadHit->HitLoci.Hit.Seg[1].ChromID = 0;
			}
		else
			MatchesFiltIn += 1;
		}
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Filtering by prioritorised regions completed - retained %u, removed %u matches",MatchesFiltIn,MatchesFiltOut);

if(MatchesFiltOut)
	SortReadHits(eRSMHitMatch,false,true);

if(gProcessingID > 0)
	{
	gSQLiteSummaries.AddResult(gProcessingID,(char *)"Filtered",ePTInt32,sizeof(MatchesFiltIn),"RegionRetained",&MatchesFiltIn);
	gSQLiteSummaries.AddResult(gProcessingID,(char *)"Filtered",ePTInt32,sizeof(MatchesFiltOut),"RegionRemoved",&MatchesFiltOut);
	}
return(eBSFSuccess);
}


int
CAligner::WriteBasicCountStats(void)
{
int Idx;
int SeqOfs;
char szLineBuff[cMaxReadLen * 8];
int BuffIdx;
int PhredBand;

if(m_hStatsFile > 0 && m_MaxAlignLen > 0)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Writing out basic count stats to file");

	// report multiple hit distribution if MLMode not simply the default
	BuffIdx = 0;
	if(m_MLMode > eMLdefault)
		{
		BuffIdx = sprintf(szLineBuff,"\"Multihit distribution\",");
		for(Idx=0; Idx < m_MaxMLmatches; Idx++)
			{
			BuffIdx += sprintf(&szLineBuff[BuffIdx],",%d",Idx+1);
			if((BuffIdx + 512) > sizeof(szLineBuff))
				{
				CUtility::SafeWrite(m_hStatsFile,szLineBuff,BuffIdx);
				BuffIdx = 0;
				}
			}
		BuffIdx += sprintf(&szLineBuff[BuffIdx],"\n,\"Instances\"");
		for(Idx=0; Idx < m_MaxMLmatches; Idx++)
			{
			BuffIdx += sprintf(&szLineBuff[BuffIdx],",%d",m_MultiHitDist[Idx]);
			if((BuffIdx + 512) > sizeof(szLineBuff))
				{
				CUtility::SafeWrite(m_hStatsFile,szLineBuff,BuffIdx);
				BuffIdx = 0;
				}
			}
		BuffIdx += sprintf(&szLineBuff[BuffIdx],"\n");
		}

	BuffIdx += sprintf(&szLineBuff[BuffIdx],"\"Phred Score Instances\",");
	for(SeqOfs = 0; SeqOfs < m_MaxAlignLen; SeqOfs++)
		{
		BuffIdx += sprintf(&szLineBuff[BuffIdx],",%d",SeqOfs+1);
		if((BuffIdx + 512) > sizeof(szLineBuff))
			{
			CUtility::SafeWrite(m_hStatsFile,szLineBuff,BuffIdx);
			BuffIdx = 0;
			}
		}

	for(PhredBand = 0; PhredBand <= 3; PhredBand++)
		{
		switch(PhredBand) {
			case 0:
				BuffIdx += sprintf(&szLineBuff[BuffIdx],"\n,\"Phred 0..9\"");
				break;
			case 1:
				BuffIdx += sprintf(&szLineBuff[BuffIdx],"\n,\"Phred 10..19\"");
				break;
			case 2:
				BuffIdx += sprintf(&szLineBuff[BuffIdx],"\n,\"Phred 20..29\"");
				break;
			case 3:
				BuffIdx += sprintf(&szLineBuff[BuffIdx],"\n,\"Phred 30+\"");
				break;
			}

		for(SeqOfs = 0; SeqOfs < m_MaxAlignLen; SeqOfs++)
			{
			BuffIdx += sprintf(&szLineBuff[BuffIdx],",%d",m_AlignQSubDist[PhredBand][SeqOfs].QInsts);
			if((BuffIdx + 512) > sizeof(szLineBuff))
				{
				CUtility::SafeWrite(m_hStatsFile,szLineBuff,BuffIdx);
				BuffIdx = 0;
				}
			}
		}

	BuffIdx += sprintf(&szLineBuff[BuffIdx],"\n\n\"Aligner Induced Subs\",");
	for(SeqOfs = 0; SeqOfs < m_MaxAlignLen; SeqOfs++)
		{
		BuffIdx += sprintf(&szLineBuff[BuffIdx],",%d",SeqOfs+1);
		if((BuffIdx + 512) > sizeof(szLineBuff))
			{
			CUtility::SafeWrite(m_hStatsFile,szLineBuff,BuffIdx);
			BuffIdx = 0;
			}
		}

	for(PhredBand = 0; PhredBand <= 3; PhredBand++)
		{
		switch(PhredBand) {
			case 0:
				BuffIdx += sprintf(&szLineBuff[BuffIdx],"\n,\"Phred 0..8\"");
				break;
			case 1:
				BuffIdx += sprintf(&szLineBuff[BuffIdx],"\n,\"Phred 9..19\"");
				break;
			case 2:
				BuffIdx += sprintf(&szLineBuff[BuffIdx],"\n,\"Phred 20..29\"");
				break;
			case 3:
				BuffIdx += sprintf(&szLineBuff[BuffIdx],"\n,\"Phred 30+\"");
				break;
			}

		for(SeqOfs = 0; SeqOfs < m_MaxAlignLen; SeqOfs++)
			{
			BuffIdx += sprintf(&szLineBuff[BuffIdx],",%d",m_AlignQSubDist[PhredBand][SeqOfs].Subs);
			if((BuffIdx + 512) > sizeof(szLineBuff))
				{
				CUtility::SafeWrite(m_hStatsFile,szLineBuff,BuffIdx);
				BuffIdx = 0;
				}
			}
		}

	BuffIdx += sprintf(&szLineBuff[BuffIdx],"\n\n\"Multiple substitutions\",");
	for(SeqOfs = 0; SeqOfs < m_MaxAlignLen; SeqOfs++)
		{
		BuffIdx += sprintf(&szLineBuff[BuffIdx],",%d",SeqOfs);
		if((BuffIdx + 512) > sizeof(szLineBuff))
			{
			CUtility::SafeWrite(m_hStatsFile,szLineBuff,BuffIdx);
			BuffIdx = 0;
			}
		}
	BuffIdx += sprintf(&szLineBuff[BuffIdx],"\n,\"Instances\"");
	for(SeqOfs = 0; SeqOfs < m_MaxAlignLen; SeqOfs++)
		{
		BuffIdx += sprintf(&szLineBuff[BuffIdx],",%d",m_AlignMSubDist[SeqOfs]);
		if((BuffIdx + 512) > sizeof(szLineBuff))
			{
			CUtility::SafeWrite(m_hStatsFile,szLineBuff,BuffIdx);
			BuffIdx = 0;
			}
		}
	BuffIdx += sprintf(&szLineBuff[BuffIdx],"\n");
	CUtility::SafeWrite(m_hStatsFile,szLineBuff,BuffIdx);
	BuffIdx = 0;
	}
return(eBSFSuccess);
}


bool								// true if file to be generated compressed with gzopen/gzwrite/gzclose
CAligner::FileReqWriteCompr(char *pszFile) // If last 3 chars of file name is ".gz" then this file is assumed to require compression
{
int Len;
if(pszFile == NULL || pszFile[0] == '\0')
	return(false);
if((Len = (int)strlen(pszFile)) < 4)
	return(false);
return(stricmp(".gz",&pszFile[Len-3]) == 0 ? true : false);
}

// Create, or if file already exists then truncate, the multitude of results files which user may have requested
// If a primary output (m_pszOutFile, m_pszNoneAlignFile or m_pszMultiAlignFile) file name ends with '.gz' then will create file with gzopen ready for output compression
// If output is to SAM or BAM then CSAMfile will handle the file processing
int
CAligner::CreateOrTruncResultFiles(void)
{
int LineLen;

// determine which primary output files need to be generated as compressed
m_bgzOutFile = FileReqWriteCompr(m_pszOutFile);
m_bgzNoneAlignFile = FileReqWriteCompr(m_pszNoneAlignFile);
m_bgzMultiAlignFile = FileReqWriteCompr(m_pszMultiAlignFile);


// create/truncate all output reporting files
if(m_FMode < eFMsam)
	{
	if(!m_bgzOutFile)
		{
#ifdef _WIN32
		m_hOutFile = open(m_pszOutFile,( O_WRONLY | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC),(_S_IREAD | _S_IWRITE));
#else
		if((m_hOutFile = open(m_pszOutFile,O_WRONLY | O_CREAT,S_IREAD | S_IWRITE))!=-1)
			if(ftruncate(m_hOutFile,0)!=0)
					{
					gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate %s - %s",m_pszOutFile,strerror(errno));
					return(eBSFerrCreateFile);
					}
#endif
		if(m_hOutFile < 0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate output file '%s'",m_pszOutFile);
			return(eBSFerrCreateFile);
			}
		}
	else
		{
		m_gzOutFile = gzopen(m_pszOutFile,"wb");
		if(m_gzOutFile == NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate output file '%s'",m_pszOutFile);
			return(eBSFerrCreateFile);
			}
		gzbuffer(m_gzOutFile,cAllocLineBuffSize);		// large buffer to reduce number of writes required
		}
	}
else
	{
	m_hOutFile = -1;
	m_gzOutFile = NULL;
	}

if((m_pszLineBuff = new char [cAllocLineBuffSize])==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to allocate memory (%d bytes) for output line buffer",cAllocLineBuffSize);
	return(eBSFerrMem);
	}

if(m_FMode == eFMbed && m_MLMode == eMLall)
	{
	LineLen = sprintf(m_pszLineBuff,"track type=bed name=\"%s\" description=\"%s\"\n",m_pszTrackTitle,m_pszTrackTitle);
	if(!m_bgzOutFile)
		CUtility::SafeWrite(m_hOutFile,m_pszLineBuff,LineLen);
	else
		CUtility::SafeWrite_gz(m_gzOutFile,m_pszLineBuff,LineLen);
	LineLen = 0;
	}
LineLen = 0;
m_szLineBuffIdx = 0;

// If BED output requested then separate BEDs for InDels and splice junctions
if(m_microInDelLen && m_FMode == eFMbed)
	{
	strcpy(m_szIndRsltsFile,m_pszOutFile);
	if(m_bgzOutFile)								// if compressing the primary alignment results file then remove the ".gz' file suffix 
		m_szJctRsltsFile[strlen(m_pszOutFile)-3] = '\0';
	strcat(m_szIndRsltsFile,".ind");
#ifdef _WIN32
	m_hIndOutFile = open(m_szIndRsltsFile,( O_WRONLY | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC),(_S_IREAD | _S_IWRITE) );
#else
	if((m_hIndOutFile = open(m_szIndRsltsFile,O_WRONLY | O_CREAT,S_IREAD | S_IWRITE))!=-1)
		if(ftruncate(m_hIndOutFile,0)!=0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate microInDel file %s - %s",m_szIndRsltsFile,strerror(errno));
			return(eBSFerrCreateFile);
			}
#endif
	if(m_hIndOutFile < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate microInDel output file '%s'",m_szIndRsltsFile);
		return(eBSFerrCreateFile);
		}
	}

if(m_SpliceJunctLen && m_FMode == eFMbed)
	{
	strcpy(m_szJctRsltsFile,m_pszOutFile);
	if(m_bgzOutFile)								// if compressing the primary alignment results file then remove the ".gz' file suffix 
		m_szJctRsltsFile[strlen(m_pszOutFile)-3] = '\0';
	strcat(m_szJctRsltsFile,".jct");
#ifdef _WIN32
	m_hJctOutFile = open(m_szJctRsltsFile,( O_WRONLY | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC),(_S_IREAD | _S_IWRITE) );
#else
	if((m_hJctOutFile = open(m_szJctRsltsFile,O_WRONLY | O_CREAT,S_IREAD | S_IWRITE))!=-1)
		if(ftruncate(m_hJctOutFile,0)!=0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate  splice junct  %s - %s",m_szJctRsltsFile,strerror(errno));
			return(eBSFerrCreateFile);
			}
#endif
	if(m_hJctOutFile < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate splice junct output file '%s'",m_szJctRsltsFile);
		return(eBSFerrCreateFile);
		}
	}


if(m_pszSitePrefsFile != NULL && m_pszSitePrefsFile[0] != '\0')
	{
#ifdef _WIN32
	m_hSitePrefsFile = open(m_pszSitePrefsFile,( O_WRONLY | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC),(_S_IREAD | _S_IWRITE) );
#else
	if((m_hSitePrefsFile = open(m_pszSitePrefsFile,O_WRONLY | O_CREAT,S_IREAD | S_IWRITE))!=-1)
	   if(ftruncate(m_hSitePrefsFile,0)!=0)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate %s - %s",m_pszSitePrefsFile,strerror(errno));
				return(eBSFerrCreateFile);
				}
#endif

	if(m_hSitePrefsFile < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate site preferencing file '%s'",m_pszSitePrefsFile);
		return(eBSFerrCreateFile);
		}
	}
else
	m_hSitePrefsFile = -1;


if(m_pszSNPRsltsFile != NULL && m_pszSNPRsltsFile[0] != '\0' && m_MinSNPreads > 0 && m_FMode <= eFMsam)
	{
#ifdef _WIN32
	m_hSNPfile = open(m_pszSNPRsltsFile,( O_WRONLY | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC),(_S_IREAD | _S_IWRITE) );
#else
	if((m_hSNPfile = open(m_pszSNPRsltsFile,O_WRONLY | O_CREAT,S_IREAD | S_IWRITE))!=-1)
		if(ftruncate(m_hSNPfile,0)!=0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate SNP file  %s - %s",m_pszSNPRsltsFile,strerror(errno));
			return(eBSFerrCreateFile);
			}
#endif
	if(m_hSNPfile < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate SNP output file '%s'",m_pszSNPRsltsFile);
		return(eBSFerrCreateFile);
		}

	if(m_pszSNPCentroidFile != NULL && m_pszSNPCentroidFile[0] != '\0')
		{
#ifdef _WIN32
		m_hSNPCentsfile = open(m_pszSNPCentroidFile,( O_WRONLY | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC),(_S_IREAD | _S_IWRITE) );
#else
		if((m_hSNPCentsfile = open(m_pszSNPCentroidFile,O_WRONLY | O_CREAT,S_IREAD | S_IWRITE))!=-1)
			if(ftruncate(m_hSNPCentsfile,0)!=0)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate SNP file  %s - %s",m_pszSNPCentroidFile,strerror(errno));
				return(eBSFerrCreateFile);
				}
#endif
		if(m_hSNPCentsfile < 0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate SNP output file '%s'",m_pszSNPCentroidFile);
			return(eBSFerrCreateFile);
			}
		}
	else
		m_hSNPCentsfile = -1;

#ifdef _DISNPS_
	strcpy(m_szDiSNPFile,m_pszSNPRsltsFile);	// the DiSNP file name should be parameterised!!!
	strcat(m_szDiSNPFile,(char *)".disnp.csv");
	strcpy(m_szTriSNPFile,m_pszSNPRsltsFile);	// the DiSNP file name should be parameterised!!!
	strcat(m_szTriSNPFile,(char *)".trisnp.csv");


	if(m_szDiSNPFile[0] != '\0')
		{
#ifdef _WIN32
		m_hDiSNPfile = open(m_szDiSNPFile,( O_WRONLY | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC),(_S_IREAD | _S_IWRITE) );
#else
		if((m_hDiSNPfile = open(m_szDiSNPFile,O_WRONLY | O_CREAT,S_IREAD | S_IWRITE))!=-1)
			if(ftruncate(m_hDiSNPfile,0)!=0)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate DiSNP file  %s - %s",m_szDiSNPFile,strerror(errno));
				return(eBSFerrCreateFile);
				}
#endif
		if(m_hDiSNPfile < 0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate DiSNP output file '%s'",m_szDiSNPFile);
			return(eBSFerrCreateFile);
			}
		}
	else
		m_hDiSNPfile = -1;

	if(m_szTriSNPFile[0] != '\0')
		{
#ifdef _WIN32
		m_hTriSNPfile = open(m_szTriSNPFile,( O_WRONLY | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC),(_S_IREAD | _S_IWRITE) );
#else
		if((m_hTriSNPfile = open(m_szTriSNPFile,O_WRONLY | O_CREAT,S_IREAD | S_IWRITE))!=-1)
			if(ftruncate(m_hTriSNPfile,0)!=0)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate TriSNP file  %s - %s",m_szTriSNPFile,strerror(errno));
				return(eBSFerrCreateFile);
				}
#endif
		if(m_hTriSNPfile < 0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate TriSNP output file '%s'",m_szTriSNPFile);
			return(eBSFerrCreateFile);
			}
		}
	else
		m_hTriSNPfile = -1;
#endif

	if(m_pszMarkerFile != NULL && m_pszMarkerFile[0] != '\0')
		{
#ifdef _WIN32
		m_hMarkerFile = open(m_pszMarkerFile,( O_WRONLY | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC),(_S_IREAD | _S_IWRITE) );
#else
		if((m_hMarkerFile = open(m_pszMarkerFile,O_WRONLY | O_CREAT,S_IREAD | S_IWRITE))!=-1)
			if(ftruncate(m_hMarkerFile,0)!=0)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate SNP file  %s - %s",m_pszMarkerFile,strerror(errno));
				return(eBSFerrCreateFile);
				}
#endif
		if(m_hMarkerFile < 0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate SNP output file '%s'",m_pszMarkerFile);
			return(eBSFerrCreateFile);
			}
		}
	else
		m_hMarkerFile = -1;
	}
else
	{
	m_hSNPfile = -1;
	m_hSNPCentsfile = -1;
	m_hMarkerFile = -1;
	m_hDiSNPfile = -1;
	m_hTriSNPfile = -1;
	};

// none-aligned fasta reads file to be also generated?
if(m_pszNoneAlignFile != NULL && m_pszNoneAlignFile[0] != '\0')
	{
	if(!m_bgzNoneAlignFile)
		{
#ifdef _WIN32
		m_hNoneAlignFile = open(m_pszNoneAlignFile,( O_WRONLY | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC),(_S_IREAD | _S_IWRITE) );
#else
		if((m_hNoneAlignFile = open(m_pszNoneAlignFile,O_WRONLY | O_CREAT,S_IREAD | S_IWRITE))!=-1)
			if(ftruncate(m_hNoneAlignFile,0)!=0)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate nonaligned %s - %s",m_pszNoneAlignFile,strerror(errno));
				return(eBSFerrCreateFile);
				}
#endif

		if(m_hNoneAlignFile < 0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate output nonaligned file '%s'",m_pszNoneAlignFile);
			return(eBSFerrCreateFile);
			}
		}
	else
		{
		m_gzNoneAlignFile = gzopen(m_pszNoneAlignFile,"wb");
		if(m_gzNoneAlignFile == NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate nonaligned file '%s'",m_pszNoneAlignFile);
			return(eBSFerrCreateFile);
			}
		gzbuffer(m_gzNoneAlignFile,cAllocLineBuffSize);		// large buffer to reduce number of writes required
		}
	}
else
	{
	m_hNoneAlignFile = -1;
	m_gzNoneAlignFile = NULL;
	}

// multialigned fasta reads file to be also generated?
if(m_pszMultiAlignFile != NULL && m_pszMultiAlignFile[0] != '\0')
	{
	if(!m_bgzMultiAlignFile)
		{
#ifdef _WIN32
		m_hMultiAlignFile = open(m_pszMultiAlignFile,( O_WRONLY | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC),(_S_IREAD | _S_IWRITE) );
#else
		if((m_hMultiAlignFile = open(m_pszMultiAlignFile,O_WRONLY | O_CREAT,S_IREAD | S_IWRITE))!=-1)
			if(ftruncate(m_hMultiAlignFile,0)!=0)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate multiple alignment file %s - %s",m_pszMultiAlignFile,strerror(errno));
				return(eBSFerrCreateFile);
				}
#endif

		if(m_hMultiAlignFile < 0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate multiple alignment file '%s'",m_pszMultiAlignFile);
			return(eBSFerrCreateFile);
			}
		}
	else
		{
		m_gzMultiAlignFile = gzopen(m_pszMultiAlignFile,"wb");
		if(m_gzMultiAlignFile == NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate multialigned file '%s'",m_pszMultiAlignFile);
			return(eBSFerrCreateFile);
			}
		gzbuffer(m_gzMultiAlignFile,cAllocLineBuffSize);		// large buffer to reduce number of writes required
		}
	}
else
	{
	m_hMultiAlignFile = -1;
	m_bgzMultiAlignFile = NULL;
	}

// substitution stats to be also generated? Note that these can't be generated for InDels or splice junctions
if(m_pszStatsFile != NULL && m_pszStatsFile[0] != '\0')
	{
#ifdef _WIN32
	m_hStatsFile = open(m_pszStatsFile,( O_WRONLY | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC),(_S_IREAD | _S_IWRITE) );
#else
	if((m_hStatsFile = open(m_pszStatsFile,O_WRONLY | O_CREAT,S_IREAD | S_IWRITE))!=-1)
		if(ftruncate(m_hStatsFile,0)!=0)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate stats %s - %s",m_pszStatsFile,strerror(errno));
				return(eBSFerrCreateFile);
				}
#endif

	if(m_hStatsFile < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate output stats file '%s'",m_pszStatsFile);
		return(eBSFerrCreateFile);
		}
	}
else
	m_hStatsFile = -1;


// PE insert lengths to be also generated? Note that these can't be generated for InDels or splice junctions
if(m_pszPEInsertDistFile != NULL && m_pszPEInsertDistFile[0] != '\0')
	{
#ifdef _WIN32
	m_hInsertLensFile = open(m_pszPEInsertDistFile,( O_WRONLY | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC),(_S_IREAD | _S_IWRITE) );
#else
	if((m_hInsertLensFile = open(m_pszPEInsertDistFile,O_WRONLY | O_CREAT,S_IREAD | S_IWRITE))!=-1)
		if(ftruncate(m_hInsertLensFile,0)!=0)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate stats %s - %s",m_pszPEInsertDistFile,strerror(errno));
				return(eBSFerrCreateFile);
				}
#endif

	if(m_hInsertLensFile < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate output stats file '%s'",m_pszPEInsertDistFile);
		return(eBSFerrCreateFile);
		}
	}
else
	m_hInsertLensFile = -1;

return(eBSFSuccess);
}

int
CAligner::CompileChromRegExprs(int	NumIncludeChroms,	// number of chromosome regular expressions to include
		char **ppszIncludeChroms,		// array of include chromosome regular expressions
		int	NumExcludeChroms,			// number of chromosome expressions to exclude
		char **ppszExcludeChroms)		// array of exclude chromosome regular expressions
{
int Idx;

#ifndef _WIN32
int RegErr;				// regular expression parsing error
char szRegErr[128];			// to hold RegErr as textual representation ==> regerror();
#endif

#ifdef _WIN32
try {
	for(Idx=0;Idx < NumIncludeChroms;Idx++)
		{
		m_IncludeChromsRE[Idx] = new Regexp();
		m_IncludeChromsRE[Idx]->Parse(ppszIncludeChroms[Idx],false);	// note case insensitive
		}
	}
catch(...)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to process include regexpr chrom '%s'",ppszIncludeChroms[Idx]);
	return(eBSFerrMem);
	}
try {
	for(Idx=0;Idx < NumExcludeChroms;Idx++)
		{
		m_ExcludeChromsRE[Idx] = new Regexp();
		m_ExcludeChromsRE[Idx]->Parse(ppszExcludeChroms[Idx],false);	// note case insensitive
		}
	}
catch(...)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to process exclude regexpr chrom '%s'",ppszExcludeChroms[Idx]);
	return(eBSFerrMem);
	}

#else
for(Idx=0;Idx < NumIncludeChroms;Idx++)
	{

	RegErr=regcomp(&m_IncludeChromsRE[Idx],ppszIncludeChroms[Idx],REG_EXTENDED | REG_ICASE);	// note case insensitive
	if(RegErr)
		{
		regerror(RegErr,&m_IncludeChromsRE[Idx],szRegErr,sizeof(szRegErr));
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to process include chrom '%s' error: %s",ppszIncludeChroms[Idx],szRegErr);
		return(eBSFerrMem);
		}
	}
for(Idx=0;Idx < NumExcludeChroms;Idx++)
	{
	RegErr = regcomp(&m_ExcludeChromsRE[Idx],ppszExcludeChroms[Idx],REG_EXTENDED | REG_ICASE);	// note case insensitive
	if(RegErr)
		{
		regerror(RegErr,&m_ExcludeChromsRE[Idx],szRegErr,sizeof(szRegErr));
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to process exclude chrom '%s' error: %s",ppszExcludeChroms[Idx],szRegErr);
		return(eBSFerrMem);
		}
	}
#endif
return(eBSFSuccess);
}

#ifdef _WIN32
unsigned __stdcall LoadReadFilesThread(void * pThreadPars)
#else
void *LoadReadFilesThread(void * pThreadPars)
#endif
{
int Rslt;
tsLoadReadsThreadPars *pPars = (tsLoadReadsThreadPars *)pThreadPars;			// makes it easier not having to deal with casts!
CAligner *pAligner = (CAligner *)pPars->pThis;
Rslt = pAligner->ProcLoadReadFiles(pPars);
pPars->Rslt = Rslt;
#ifdef _WIN32
_endthreadex(0);
return(eBSFSuccess);
#else
pthread_exit(NULL);
#endif
}

int
CAligner::InitiateLoadingReads(void)
{
tsLoadReadsThreadPars ThreadPars;

// initiate loading the reads
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Loading reads from file...");
m_ThreadLoadReadsRslt = -1;


ThreadPars.pRslt = &m_ThreadLoadReadsRslt;
ThreadPars.pThis = this;
ThreadPars.Rslt = 0;

#ifdef _WIN32
m_hThreadLoadReads = ThreadPars.threadHandle = (HANDLE)_beginthreadex(NULL,0x0fffff,LoadReadFilesThread,&ThreadPars,0,&m_ThreadLoadReadsID);
#else
int ThreadRslt = ThreadPars.threadRslt = pthread_create (&m_ThreadLoadReadsID , NULL , LoadReadFilesThread , &ThreadPars );
#endif

// wait a few seconds, if major problems with loading reads then should show very quickly
#ifdef _WIN32
if(WAIT_TIMEOUT != WaitForSingleObject(m_hThreadLoadReads, 3000))
	{
	CloseHandle(m_hThreadLoadReads);
	m_hThreadLoadReads = NULL;
	return(m_ThreadLoadReadsRslt);
	}
#else
struct timespec ts;
int JoinRlt;
clock_gettime(CLOCK_REALTIME, &ts);
ts.tv_sec += 3;
if((JoinRlt = pthread_timedjoin_np(m_ThreadLoadReadsID, NULL, &ts)) == 0)
	{
	m_ThreadLoadReadsID = 0;
	return(m_ThreadLoadReadsRslt);
	}
#endif
return(eBSFSuccess);
}


#ifdef _WIN32
unsigned __stdcall AssignMultiMatchesThread(void * pThreadPars)
#else
void *AssignMultiMatchesThread(void * pThreadPars)
#endif
{
int Rslt;
tsClusterThreadPars *pPars = (tsClusterThreadPars *)pThreadPars; // makes it easier not having to deal with casts!
CAligner *pThis = (CAligner *)pPars->pThis;
Rslt = pThis->ProcAssignMultiMatches(pPars);
pPars->Rslt = Rslt;
#ifdef _WIN32
_endthreadex(0);
return(eBSFSuccess);
#else
pthread_exit(NULL);
#endif
}

int
CAligner::RunClusteringThreads(int NumThreads)
{
int ThreadIdx;
UINT32 ClusterStartIdx;
tsClusterThreadPars ClusterThreads[cMaxWorkerThreads];

ClusterStartIdx = 0;
memset(ClusterThreads,0,sizeof(ClusterThreads));
for(ThreadIdx = 0; ThreadIdx < NumThreads; ThreadIdx++)
	{
	ClusterThreads[ThreadIdx].ThreadIdx = ThreadIdx + 1;
	ClusterThreads[ThreadIdx].pThis = this;
#ifdef _WIN32
	ClusterThreads[ThreadIdx].threadHandle = (HANDLE)_beginthreadex(NULL,0x0fffff,AssignMultiMatchesThread,&ClusterThreads[ThreadIdx],0,&ClusterThreads[ThreadIdx].threadID);
#else
	ClusterThreads[ThreadIdx].threadRslt =	pthread_create (&ClusterThreads[ThreadIdx].threadID , NULL , AssignMultiMatchesThread , &ClusterThreads[ThreadIdx] );
#endif
	}

for(ThreadIdx = 0; ThreadIdx < NumThreads; ThreadIdx++)
	{
#ifdef _WIN32
	while(WAIT_TIMEOUT == WaitForSingleObject( ClusterThreads[ThreadIdx].threadHandle, 60000))
		{
		}
	CloseHandle( ClusterThreads[ThreadIdx].threadHandle);
#else
	struct timespec ts;
	int JoinRlt;
	clock_gettime(CLOCK_REALTIME, &ts);
	ts.tv_sec += 60;
	while((JoinRlt = pthread_timedjoin_np(ClusterThreads[ThreadIdx].threadID, NULL, &ts)) != 0)
		{
		ts.tv_sec += 60;
		}
#endif
	}
return(eBSFSuccess);
}

// GetClusterStartEnd
// Clustering threads call this function to obtain the next from, until inclusive indexes into m_pMultiHits[] to be clustered
// Returns number of
int												// returns 0 if finished clustering or cnt of multihits to be processed by this thread
CAligner::GetClusterStartEnd(UINT32 *pMatchFrom,			// cluster from this inclusive index
					UINT32 *pMatchUntil)		// until this inclusive index
{
int Rslt;
UINT32 NumLeft2Cluster;
UINT32 Num4Thread;
AcquireSerialise();
NumLeft2Cluster = m_NumMultiHits - m_CurClusterFrom;
if(NumLeft2Cluster > 0)
	{
	if(NumLeft2Cluster < 100)	// if < 100 yet to be processed then give it all to the one thread
		Num4Thread = NumLeft2Cluster;
	else
		{
		Num4Thread = min(2000,m_NumThreads + (NumLeft2Cluster / (UINT32)m_NumThreads));
		Num4Thread = min(Num4Thread,NumLeft2Cluster);
		}

	*pMatchFrom = m_CurClusterFrom;
	*pMatchUntil = m_CurClusterFrom + Num4Thread - 1;
	m_CurClusterFrom = min(1 + *pMatchUntil,m_NumMultiHits);
	Rslt = Num4Thread;
	}
else
	Rslt = 0;
ReleaseSerialise();
return(Rslt);
}


int
CAligner::ProcAssignMultiMatches(tsClusterThreadPars *pPars)
{
int Rslt;
UINT32 MatchFrom;
UINT32 MatchUntil;
UINT32 HitIdx;
UINT32 Score;
UINT32 ClustHitIdx;
tsReadHit *pCurHit;
tsReadHit *pClustHit;
tsReadHit *pPrevProcCurHit;
int Overlap;
UINT32 ClustEndLoci;

while((Rslt = GetClusterStartEnd(&MatchFrom,&MatchUntil)) > 0)
	{
	pCurHit = &m_pMultiHits[MatchFrom];
	pPrevProcCurHit = NULL;
	for(HitIdx=MatchFrom;HitIdx <= MatchUntil; HitIdx++,pCurHit++)
		{
		if(!pCurHit->HitLoci.FlagMH)			// only interested in assigning reads which align to multiple hit loci
			continue;							// reads which map to a single hit loci do not require processing

		if(pPrevProcCurHit != NULL &&
					AdjStartLoci(&pPrevProcCurHit->HitLoci.Hit.Seg[0]) == AdjStartLoci(&pCurHit->HitLoci.Hit.Seg[0]) &&
					AdjHitLen(&pPrevProcCurHit->HitLoci.Hit.Seg[0]) == AdjHitLen(&pCurHit->HitLoci.Hit.Seg[0]) &&
					pPrevProcCurHit->HitLoci.Hit.Seg[0].Strand == pCurHit->HitLoci.Hit.Seg[0].Strand &&
					pPrevProcCurHit->HitLoci.Hit.Seg[0].ChromID == pCurHit->HitLoci.Hit.Seg[0].ChromID)
			{
			pCurHit->HitLoci.Hit.Score = pPrevProcCurHit->HitLoci.Hit.Score;
			continue;
			}

		pPrevProcCurHit = NULL;
		pCurHit->HitLoci.Hit.Score = 0;
		pClustHit = pCurHit;
		ClustHitIdx = HitIdx;
		while(ClustHitIdx-- > 0)				// checking for clustering upstream of current hit loci
			{
			pClustHit -= 1;

			// can't cluster with reads on a different chrom!
			if(pClustHit->HitLoci.Hit.Seg[0].ChromID != pCurHit->HitLoci.Hit.Seg[0].ChromID)
				break;

			// can't cluster if much too far away
			if((AdjStartLoci(&pCurHit->HitLoci.Hit.Seg[0]) - AdjStartLoci(&pClustHit->HitLoci.Hit.Seg[0])) >= (int)m_MaxReadsLen)
				break;
			// only cluster if >= cClustMultiOverLap
			ClustEndLoci = AdjEndLoci(&pClustHit->HitLoci.Hit.Seg[0]);
			if((int)ClustEndLoci < (AdjStartLoci(&pCurHit->HitLoci.Hit.Seg[0]) + cClustMultiOverLap))
				continue;
			Overlap = min(AdjHitLen(&pCurHit->HitLoci.Hit.Seg[0]),(int)ClustEndLoci - AdjStartLoci(&pCurHit->HitLoci.Hit.Seg[0]));

			if((m_MLMode == eMLuniq && pClustHit->HitLoci.FlagMH) ||
				(pCurHit->HitLoci.Hit.Score & cUniqueClustFlg && (pCurHit->HitLoci.Hit.Score & ~cUniqueClustFlg) >= 0x01fff))
				continue;

			if(pClustHit->HitLoci.Hit.Seg[0].Strand == pCurHit->HitLoci.Hit.Seg[0].Strand && pClustHit->ReadID != pCurHit->ReadID)
				{
				if(!pClustHit->HitLoci.FlagMH)	// clustering to a unique aligned reads has much higher priority than to other multialigned reads
					{
					Score = 1 + (Overlap * cClustUniqueScore)/cClustScaleFact;
					if(pCurHit->HitLoci.Hit.Score & cUniqueClustFlg)
						Score += pCurHit->HitLoci.Hit.Score & ~cUniqueClustFlg;
					if(Score > 0x01fff)			// clamp upstream scores to be no more than 0x01fff so still room for dnstream scores
						Score = 0x01fff;
					pCurHit->HitLoci.Hit.Score = (UINT16)(Score | cUniqueClustFlg);	// flag that this score is because now clustering to uniquely aligned reads
					if(Score == 0x01fff)
						break;					// this is a good candidate loci
					}
				else							// never seen a uniquely aligned so still clustering to other non-unique aligned reads
					{
					if(!(pCurHit->HitLoci.Hit.Score & cUniqueClustFlg))
						{
						Score = 1 + (Overlap * cClustMultiScore)/cClustScaleFact;
						Score += pCurHit->HitLoci.Hit.Score & ~cUniqueClustFlg;
						if(Score > 0x01fff)			// clamp upstream scores to be no more than 0x01fff so still room for dnstream scores
							Score = 0x01fff;
						pCurHit->HitLoci.Hit.Score = Score;
						}
					}
				}
			}

		// now cluster downstream
		pClustHit = pCurHit;
		ClustHitIdx = HitIdx;
		while(++ClustHitIdx < m_NumMultiHits)				// checking for clustering downstream of current hit loci
			{
			pClustHit += 1;
			// can't cluster with reads on a different chrom!
			if(pClustHit->HitLoci.Hit.Seg[0].ChromID != pCurHit->HitLoci.Hit.Seg[0].ChromID)
				break;

			// can't score if no overlap of at least cClustMultiOverLap
			ClustEndLoci = AdjEndLoci(&pCurHit->HitLoci.Hit.Seg[0]);
			if(AdjStartLoci(&pClustHit->HitLoci.Hit.Seg[0]) > (ClustEndLoci - (UINT32)cClustMultiOverLap))
				break;

			Overlap = min(AdjHitLen(&pClustHit->HitLoci.Hit.Seg[0]),(int)ClustEndLoci - AdjStartLoci(&pClustHit->HitLoci.Hit.Seg[0]));

			if((m_MLMode == eMLuniq && pClustHit->HitLoci.FlagMH) ||
				(pCurHit->HitLoci.Hit.Score & cUniqueClustFlg && (pCurHit->HitLoci.Hit.Score & ~cUniqueClustFlg) >= 0x03fff))
				continue;

			if(pClustHit->HitLoci.Hit.Seg[0].Strand == pCurHit->HitLoci.Hit.Seg[0].Strand && pClustHit->ReadID != pCurHit->ReadID)
				{
				if(!pClustHit->HitLoci.FlagMH)	// clustering to a unique aligned reads has much higher priority than to other multialigned reads
					{
					Score = 1 + (Overlap * cClustUniqueScore)/cClustScaleFact;
					if(pCurHit->HitLoci.Hit.Score & cUniqueClustFlg)
						Score += pCurHit->HitLoci.Hit.Score & ~cUniqueClustFlg;
					if(Score > 0x3fff)
						Score = 0x3fff;
					pCurHit->HitLoci.Hit.Score = (UINT16)(Score | cUniqueClustFlg);
					if(Score == 0x3fff)
						break;
					}
				else
					{
					if(!(pCurHit->HitLoci.Hit.Score & cUniqueClustFlg))
						{
						Score = 1 + (Overlap * cClustMultiScore)/cClustScaleFact;
						Score += pCurHit->HitLoci.Hit.Score & ~cUniqueClustFlg;
						if(Score > 0x3fff)
							Score = 0x3fff;
						pCurHit->HitLoci.Hit.Score = (UINT16)Score;
						}
					}
				}
			}
		pPrevProcCurHit = pCurHit;
		}
	}
pPars->Rslt = 1;
#ifdef _WIN32
_endthreadex(0);
return(eBSFSuccess);
#else
pthread_exit(NULL);
#endif
}

// AssignMultiMatches
// Use clustering (within a sliding window) to determine which matching loci should be assigned to those reads with multiple hits
// This is multipass clustering process as need to cluster multimatches with uniquely matched before clustering with other multimatched.
int
CAligner::AssignMultiMatches(void) // false to cluster with uniques, true to cluster with multimatches
{
UINT32 HitIdx;
UINT32 ClustHitIdx;
UINT32 Distance;
int NumAssigned;
int NumUnlocated;
tsReadHit *pCurHit;
tsReadHit *pClustHit;

int ClusterUniqueAssigned;
int ClusterAllAssigned;

if(m_MLMode <= eMLrand)		// nothing to do?
	return(0);

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Assigning %d reads which aligned to multiple loci to a single loci",m_NumProvMultiAligned);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Sorting...");
SortReadHits(eRSMReadID,false);
m_mtqsort.qsort(m_pMultiHits,m_NumMultiHits,sizeof(tsReadHit),SortMultiHits);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Sorting completed, now clustering...");

RunClusteringThreads(m_NumThreads);

// sort now by ascending ReadID and descending scores
// and assign the read match with the highest score to that read
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Assigning...");
m_mtqsort.qsort(m_pMultiHits,m_NumMultiHits,sizeof(tsReadHit),SortMultiHitReadIDs);
UINT32 CurReadID = 0;
UINT32 BestScore;
UINT32 NxtBestScore;

tsReadHit *pAssign2Read;
NumAssigned = 0;
NumUnlocated = 0;
ClusterUniqueAssigned = 0;
ClusterAllAssigned = 0;
pCurHit = m_pMultiHits;
for(HitIdx=0;HitIdx < m_NumMultiHits; HitIdx++,pCurHit++)
	{
	if(!pCurHit->HitLoci.FlagMH)	 // only interested in reads with multiple hits
		continue;
	if(CurReadID == pCurHit->ReadID) // only interested in the first for each read
		continue;
	CurReadID = pCurHit->ReadID;
	NumAssigned += 1;

	// only interested if best score for read is at least cMHminScore and that score is at least 2x next best score for read
	BestScore = (UINT32)(pCurHit->HitLoci.Hit.Score & ~cUniqueClustFlg);
	if(BestScore < cMHminScore)
		continue;

	if((pCurHit->HitLoci.Hit.Score & cUniqueClustFlg) == (pCurHit[1].HitLoci.Hit.Score & cUniqueClustFlg))
		{
		NxtBestScore = (UINT32)(pCurHit[1].HitLoci.Hit.Score & ~cUniqueClustFlg);
		if(BestScore < (NxtBestScore * 2))
			continue;
		}

	// accept this loci as the loci as being the putative read hit loci for this multihit
	pCurHit->HitLoci.FlagMHA = 1;
	if(pCurHit->HitLoci.Hit.Score & cUniqueClustFlg)
		{
		ClusterUniqueAssigned += 1;
		pCurHit->HitLoci.FlagHL = eHLclustunique;
		}
	else
		{
		ClusterAllAssigned += 1;
		pCurHit->HitLoci.FlagHL = eHLclustany;
		}
	}

// need to ensure that as a result of the assignments no assigned multialigned read is now actually an orphan
// orphans are those assigned as eHLclustany if within the window there are no other assigned multireads...
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Checking for orphans (unclustered) from %d putative assignments..",NumAssigned);
m_mtqsort.qsort(m_pMultiHits,m_NumMultiHits,sizeof(tsReadHit),SortMultiHits);
NumUnlocated = 0;
ClusterUniqueAssigned = 0;
ClusterAllAssigned = 0;
int PutativeAssignments = NumAssigned;
bool bAcceptMulti;
NumAssigned = 0;
pCurHit = m_pMultiHits;
for(HitIdx=0;HitIdx < m_NumMultiHits; HitIdx++,pCurHit++)
	{
	if(!pCurHit->HitLoci.FlagMHA)			// only interested in multihit reads which have been assigned to a single loci
		continue;

	bAcceptMulti = false;
	if(pCurHit->HitLoci.FlagHL == eHLclustany) // if was clustered with another multihit then ensure that the other multihit still in play..
		{
		pClustHit = pCurHit;
		ClustHitIdx = HitIdx;
		while(ClustHitIdx-- > 0)				// checking for clustering upstream of current hit loci
			{
			pClustHit -= 1;
			Distance = AdjStartLoci(&pCurHit->HitLoci.Hit.Seg[0]) - AdjStartLoci(&pClustHit->HitLoci.Hit.Seg[0]);
			if(Distance > (UINT32)(cClustMultiOverLap + (int)AdjHitLen(&pClustHit->HitLoci.Hit.Seg[0]))) // finish if much too far away
				break;
			if(pClustHit->HitLoci.Hit.Seg[0].ChromID != pCurHit->HitLoci.Hit.Seg[0].ChromID)
				break;							// can't cluster with reads on a different chrom!
			if(!pClustHit->HitLoci.FlagMH || pClustHit->HitLoci.FlagMHA == 1)
				{
				bAcceptMulti = true;
				break;
				}
			}

		if(!bAcceptMulti)
			{
			// now cluster downstream
			pClustHit = pCurHit;
			ClustHitIdx = HitIdx;
			while(++ClustHitIdx < m_NumMultiHits)				// checking for clustering downstream of current hit loci
				{
				pClustHit += 1;
				Distance = AdjStartLoci(&pClustHit->HitLoci.Hit.Seg[0]) - AdjStartLoci(&pCurHit->HitLoci.Hit.Seg[0]);
				if(Distance > (UINT32)(cClustMultiOverLap + (int)AdjHitLen(&pCurHit->HitLoci.Hit.Seg[0])))				// can't score if much too far away
					break;

				if(pClustHit->HitLoci.Hit.Seg[0].ChromID != pCurHit->HitLoci.Hit.Seg[0].ChromID)
					break;
				if(!pClustHit->HitLoci.FlagMH || pClustHit->HitLoci.FlagMHA == 1)
					{
					bAcceptMulti = true;
					break;
					}
				}
			}
		if(!bAcceptMulti)
			pCurHit->HitLoci.FlagMHA = 0;
		}
	else
		bAcceptMulti = true;

	if(bAcceptMulti == true)
		{
		if((pAssign2Read = LocateRead(pCurHit->ReadID))==NULL)
			{
			if(NumUnlocated++ < 10)
				gDiagnostics.DiagOut(eDLWarn,gszProcName,"A read (ReadID = %d) was not located in call to LocateRead()",pCurHit->ReadID);
			continue;
			}
		pAssign2Read->HitLoci = pCurHit->HitLoci;
		pAssign2Read->NumHits = 1;
		pAssign2Read->NAR = eNARAccepted;
		pAssign2Read->LowHitInstances = 1;
		if(pCurHit->HitLoci.FlagHL == eHLclustunique)
			ClusterUniqueAssigned += 1;
		else
			if(pCurHit->HitLoci.FlagHL == eHLclustany)
				ClusterAllAssigned += 1;
		NumAssigned += 1;
		}
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Clustering completed, removed %d unclustered orphans from %d putative resulting in %d (%d clustered near unique, %d clustered near other multiloci reads) multihit reads accepted as assigned",
					 PutativeAssignments-NumAssigned,PutativeAssignments,NumAssigned,ClusterUniqueAssigned,ClusterAllAssigned);
SortReadHits(eRSMHitMatch,false,true);
return(0);
}

// Median calculation
UINT32
CAligner::MedianInsertLen(UINT32 NumInserts,				// number of insert lengths in pInsertLens
				UINT32 *pInsertLens)		// insert lengths
{
UINT32 Idx;
UINT32 Idy;
UINT32 NumUnder;
UINT32 NumOver;
UINT32 *pInsertLen;
UINT32 *pInsertLen1;
UINT32 Median;
UINT32 CurMedian;
UINT32 LoMedian;
UINT32 HiMedian;

if(NumInserts == 0 || pInsertLens == NULL)
	return(0);
if(NumInserts == 1)
	return(*pInsertLens);
if(NumInserts == 2)
	return((pInsertLens[0] + pInsertLens[1])/2);

LoMedian = 0;
HiMedian = 0;
pInsertLen = pInsertLens;
for(Idx = 0; Idx < NumInserts; Idx++,pInsertLen++)
	{
	if(*pInsertLen < LoMedian || LoMedian == 0)
		LoMedian = *pInsertLen;
	if(*pInsertLen > HiMedian || HiMedian == 0)
		HiMedian = *pInsertLen;
	}
if(HiMedian == LoMedian)
	return(HiMedian);

pInsertLen = pInsertLens;
Median = *pInsertLen++;
for(Idx = 1; Idx < NumInserts; Idx++,pInsertLen++)
	{
	CurMedian = *pInsertLen;
	if(CurMedian < LoMedian || CurMedian > HiMedian)
		continue;
	
	NumUnder = 0;
	NumOver = 0;
	pInsertLen1 = pInsertLens;
	for(Idy = 0; Idy < NumInserts; Idy++,pInsertLen1++)
		{
		if(Idy == Idx || *pInsertLen1 == CurMedian)
			continue;
		if(*pInsertLen1 < CurMedian)
			NumUnder += 1;
		else
			NumOver += 1;
		}
	if(NumUnder == NumOver)
		return(CurMedian);
	if(NumUnder < NumOver)
		LoMedian = CurMedian;
	else
		HiMedian = CurMedian;
	}
return((LoMedian + HiMedian)/2);
}

// Report PE insert length distributions for each transcript or assembly contig/sequence
int
CAligner::ReportPEInsertLenDist(void)
{
UINT32 PE1Start;
UINT32 PE2Start;
UINT32 PE1Len;
UINT32 PE2Len;
UINT32 TLen;
UINT32 LenDist[53];				// holds counts for insert lengths from 100 to 600 in 10bp increments plus under 100bp and over 600bp
char szDistBuff[0x7fff];
int BuffOfs;
int Idx;
UINT32 NumPEs;
UINT64 SumInsertLens;
UINT32 MedianInsert;
UINT32 CurTargID;
UINT32 NumTargIDs;
char szTargChromName[128];
UINT32 TargChromLen;
tsReadHit *pPE1ReadHit;
tsReadHit *pPE2ReadHit;
UINT32 *pInsertLens;

if(m_hInsertLensFile == -1 || !m_bPEInsertLenDist)		// only process PE insert distributions if actually aligning PEs and user requested the insert distributions
	return(0);

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Reporting PE insert lengths for each transcript or assembly sequence");
SortReadHits(eRSMPEHitMatch,false);

pInsertLens = new UINT32 [5000000];	// big assumption that transcript will have no more than 1M paired end hits - check is made and only 1st 5M inserts are processed!

BuffOfs = sprintf(szDistBuff,"\"TargSeq\",\"TargLen\",\"TotPEs\",\"MedianInsertLen\",\"MeanInsertLen\",\"Insert <100bp\"");
for(Idx = 0; Idx < 50; Idx++)
	BuffOfs += sprintf(&szDistBuff[BuffOfs],",\"Insert %dbp\"",100 + (10 * Idx));
BuffOfs += sprintf(&szDistBuff[BuffOfs],",\"Insert >600bp\"\n");
CUtility::SafeWrite(m_hInsertLensFile,szDistBuff,BuffOfs);
BuffOfs = 0;

pPE1ReadHit = NULL;
CurTargID = 0;
NumTargIDs = 0;
NumPEs = 0;
SumInsertLens = 0;
memset(LenDist,0,sizeof(LenDist));
while((pPE1ReadHit = IterSortedReads(pPE1ReadHit))!=NULL)
	{
	if((pPE2ReadHit = IterSortedReads(pPE1ReadHit)) == NULL)
		break;

	// both ends of pair must be accepted as aligned and both must be aligned as PE
	if(!(pPE1ReadHit->NAR == eNARAccepted && pPE1ReadHit->FlgPEAligned &&  pPE2ReadHit->NAR == eNARAccepted && pPE2ReadHit->FlgPEAligned))
		break;

	// both pPE1ReadHit and pPE2ReadHit have been accepted as being PE aligned
	if(pPE1ReadHit->HitLoci.Hit.Seg[0].ChromID != (UINT32)CurTargID)  // processing a different transcript or assembly sequence?
		{
		if(NumTargIDs && NumPEs)
			{
			MedianInsert = MedianInsertLen(min(5000000,NumPEs),pInsertLens);
			BuffOfs += sprintf(&szDistBuff[BuffOfs],"\"%s\",%u,%u,%u,%u",szTargChromName,TargChromLen,NumPEs,MedianInsert,(UINT32)(SumInsertLens/NumPEs));
			for(Idx = 0; Idx < 52; Idx++)
				BuffOfs += sprintf(&szDistBuff[BuffOfs],",%d",LenDist[Idx]);
			BuffOfs += sprintf(&szDistBuff[BuffOfs],"\n");
			if((BuffOfs + 4096) > sizeof(szDistBuff))
				{
				CUtility::SafeWrite(m_hInsertLensFile,szDistBuff,BuffOfs);
				BuffOfs = 0;
				}
			}

		CurTargID = pPE1ReadHit->HitLoci.Hit.Seg[0].ChromID;
		m_pSfxArray->GetIdentName(CurTargID,sizeof(szTargChromName),szTargChromName);
		TargChromLen = m_pSfxArray->GetSeqLen(CurTargID);
		NumTargIDs += 1;
		NumPEs = 0;
		SumInsertLens = 0;
		memset(LenDist,0,sizeof(LenDist));
		}

	PE1Start = AdjAlignStartLoci(&pPE1ReadHit->HitLoci.Hit) + 1;
	PE2Start = AdjAlignStartLoci(&pPE2ReadHit->HitLoci.Hit) + 1;
	PE1Len = AdjAlignHitLen(&pPE1ReadHit->HitLoci.Hit);
	PE2Len = AdjAlignHitLen(&pPE2ReadHit->HitLoci.Hit);
	if(PE1Start <= PE2Start)
		TLen = (PE2Start - PE1Start) + PE2Len;
	else
		TLen = (PE1Start - PE2Start) + PE1Len;

	SumInsertLens += TLen;
	if(TLen < 100)
		LenDist[0] += 1;
	else
		{
		if(TLen > 600)
			LenDist[51] += 1;
		else
			LenDist[1 + ((TLen - 99)/10)] += 1;
		}
	if(NumPEs < 1000000)
		pInsertLens[NumPEs] = TLen;
	NumPEs += 1;
	pPE1ReadHit = pPE2ReadHit;
	}

if(NumTargIDs)
	{
	if(NumPEs)
		{
		MedianInsert = MedianInsertLen(min(5000000,NumPEs),pInsertLens);
		BuffOfs += sprintf(&szDistBuff[BuffOfs],"\"%s\",%u,%u,%u,%u",szTargChromName,TargChromLen,NumPEs,MedianInsert,(UINT32)(SumInsertLens/NumPEs));
		for(Idx = 0; Idx < 52; Idx++)
			BuffOfs += sprintf(&szDistBuff[BuffOfs],",%u",LenDist[Idx]);
		BuffOfs += sprintf(&szDistBuff[BuffOfs],"\n");
		}

	if(BuffOfs)
		CUtility::SafeWrite(m_hInsertLensFile,szDistBuff,BuffOfs);
	}

delete pInsertLens;
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Reporting PE insert lengths for %u transcripts or assembly sequences completed",NumTargIDs);
SortReadHits(eRSMHitMatch,false);
return(NumTargIDs);
}


// Write results as BAM or SAM format
int
CAligner::WriteBAMReadHits(etFMode ProcMode,	   // eFMsam or eFMsamAll
							teSAMFormat SAMFormat, // if SAM output format then could be SAM or BAM compressed dependent on the file extension used
							bool bPEProc,		   // true if processing paired ends
							int ComprLev)		   // BAM to be BGZF compressed at the requested level (0..9)
{
int Rslt;
tsBAMalign BAMalign;			// to hold each SAM or BAM alignment as it is constructed
eSAMFileType FileType;
CSAMfile *pSAMfile;

char szChromName[128];
int NumAlignedToSeqs;
int NumSeqsInHdr;
UINT32 NumReportedBAMreads;
UINT32 PrevNumReportedBAMreads;
tsReadHit *pReadHit;
tBSFEntryID PrevTargEntry;
int NumChroms;
bool bRptAllChroms;
UINT32 ChromSeqLen;

int ChromID;
UINT32 CurChromID;
UINT16 EntryFlags;

if((pSAMfile = new CSAMfile) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"WriteBAMReadHits: Unable to instantiate class CSAMfile");
	return(eBSFerrInternal);
	}

switch(SAMFormat) {
	case etSAMFformat:			// output SAM
		if(m_bgzOutFile)
			FileType = eSFTSAMgz;
		else
			FileType = eSFTSAM;
		break;
	case etSAMFBAM:				// output as BAM compressed with bgzf
		FileType = eSFTBAM_BAI;
		break;
	}

if((Rslt = pSAMfile->Create(FileType,m_pszOutFile,ComprLev,(char *)cpszProgVer)) < eBSFSuccess)
	{
	delete pSAMfile;
	return(Rslt);
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Sorting alignments by ascending chrom.loci");
SortReadHits(eRSMHitMatch,false);

// mark entries (chroms/contigs/sequences) for which there is at least one alignment so only these marked entries
// will be written to the SAM or BAM header
pReadHit = NULL;
PrevTargEntry = 0;
m_PrevSAMTargEntry = 0;
CurChromID = 0;
NumChroms = m_pSfxArray->GetNumEntries();
bRptAllChroms = m_MaxRptSAMSeqsThres >= NumChroms ? true : false;

// identify chroms to be reported
// only reporting those which have accepted alignments unless reporting all chromosomes even if not all have alignments
while((pReadHit = IterSortedReads(pReadHit))!=NULL)
	{
	if(pReadHit->NAR == eNARAccepted)
		{
		if(CurChromID == 0 || CurChromID != pReadHit->HitLoci.Hit.Seg[0].ChromID)
			{
			CurChromID = pReadHit->HitLoci.Hit.Seg[0].ChromID;
			if(CurChromID > 0)
				m_pSfxArray->SetResetIdentFlags(CurChromID,0x01,0x00);
			}
		}
	}

NumAlignedToSeqs = 0;
NumSeqsInHdr = 0;
for(ChromID = 1; ChromID <= NumChroms; ChromID++)
	{
	EntryFlags = m_pSfxArray->GetIdentFlags(ChromID);
	if(EntryFlags & 0x01 || bRptAllChroms)
		{
		m_pSfxArray->GetIdentName(ChromID,sizeof(szChromName),szChromName);
		ChromSeqLen=m_pSfxArray->GetSeqLen(ChromID);
		if((Rslt = pSAMfile->AddRefSeq(m_szTargSpecies,szChromName,ChromSeqLen)) < 1)
			{
			delete pSAMfile;
			return(Rslt);
			}
		NumSeqsInHdr += 1;
		if(EntryFlags & 0x01)
			NumAlignedToSeqs += 1;
		}
	}
pSAMfile->StartAlignments();

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Header written with references to %d sequences of which %d have at least 1 alignments",NumSeqsInHdr, NumAlignedToSeqs);
NumAlignedToSeqs = 0;
NumSeqsInHdr = 0;

// now write out each alignment in BAM format
pReadHit = NULL;

PrevTargEntry = 0;
m_PrevSAMTargEntry = 0;
PrevNumReportedBAMreads = 0;
NumReportedBAMreads = 0;
int RefID;
int BAMRefID;
tsReadHit *pNxtRead; 
bool bLastAligned;
bool bAcceptPE;

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Reported %s %u read alignments",SAMFormat == etSAMFformat ? "SAM" : "BAM",NumReportedBAMreads);
RefID = 0;
time_t Started = time(0);
while((pReadHit = IterSortedReads(pReadHit))!=NULL)
	{
	if(pReadHit->NAR == eNARAccepted || ProcMode == eFMsamAll)
		{
		if(pReadHit->NAR == eNARAccepted)
			{
			if(pReadHit->HitLoci.Hit.Seg[0].ChromID != (UINT32)m_PrevSAMTargEntry)
				{
				m_pSfxArray->GetIdentName(pReadHit->HitLoci.Hit.Seg[0].ChromID,sizeof(m_szSAMTargChromName),m_szSAMTargChromName);
				m_PrevSAMTargEntry = pReadHit->HitLoci.Hit.Seg[0].ChromID;
				}
			BAMRefID = 0;
			bAcceptPE = bPEProc && pReadHit->FlgPEAligned ? true : false;
			}
		else   // else also reporting reads not accepted as being aligned
			{
			BAMRefID = -1;
			m_szSAMTargChromName[0] = '*';
			m_szSAMTargChromName[1] = '\0';
			bAcceptPE = false;
			}

		if((Rslt = ReportBAMread(pReadHit,BAMRefID,bAcceptPE,&BAMalign)) < eBSFSuccess)
			return(Rslt);
		strcpy(BAMalign.szRefSeqName,m_szSAMTargChromName);

		// look ahead to check if current read is the last accepted aligned read
		bLastAligned = false;
		if(pReadHit->NAR == eNARAccepted)
			{
			pNxtRead = IterSortedReads(pReadHit);
			if(pNxtRead == NULL || pNxtRead->NAR != eNARAccepted)
				bLastAligned = true;
			}

		if((Rslt = pSAMfile->AddAlignment(&BAMalign,bLastAligned)) < eBSFSuccess)
			return(Rslt);

		NumReportedBAMreads += 1;

		if(NumReportedBAMreads > (PrevNumReportedBAMreads + 50000))
			{
			time_t Now = time(0);
			unsigned long ElapsedSecs = (unsigned long) (Now - Started);
			if(ElapsedSecs >= (60 * 10))
				{
				gDiagnostics.DiagOut(eDLInfo,gszProcName,"Reported %s %u read alignments",SAMFormat == etSAMFformat ? "SAM" : "BAM",NumReportedBAMreads);
				Started = Now;
				}
			PrevNumReportedBAMreads = NumReportedBAMreads;
			}
		// user may be interested in the distribution of the aligner induced substitutions, after any auto-trimming of flanks,
		// along the length of the reads and how this distribution relates to the quality scores
		if(m_hStatsFile != -1)
			WriteSubDist(pReadHit);
		}
	}

pSAMfile->Close();
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Completed reporting %s %u read alignments",SAMFormat == etSAMFformat ? "SAM" : "BAM",NumReportedBAMreads);
return(0);
}

// BAM index
// UINT8 magic[4];    // "BAI\1"
// UINT32 n_rf;       // number of reference sequences following
//    UINT32 n_bin;   // number of distinct bins for current reference sequence
//        UINT32 bin; // distinct bin
//        UINT32 chunks; // number of chunks following
//            UINT64 chumk_beg;		// virtual file offset at which chunk starts
//            UINT64 chumk_end;		// virtual file offset at which chunk ends
//    UINT32 n_intv;  // number of 16kb intervals for linear index
//        UINT64 ioffset;   // virtual file offset of first alignment in interval


// following BAM bin functions are copied from the specification at http://samtools.sourceforge.net/SAMv1.pdf
/* calculate bin given an alignment covering [beg,end) (zero-based, half-close-half-open) */
int CAligner::BAMreg2bin(int beg, int end)
{
--end;
if (beg>>14 == end>>14) return ((1<<15)-1)/7 + (beg>>14);
if (beg>>17 == end>>17) return ((1<<12)-1)/7 + (beg>>17);
if (beg>>20 == end>>20) return ((1<<9)-1)/7 + (beg>>20);
if (beg>>23 == end>>23) return ((1<<6)-1)/7 + (beg>>23);
if (beg>>26 == end>>26) return ((1<<3)-1)/7 + (beg>>26);
return 0;
}

/* calculate the list of bins that may overlap with region [beg,end) (zero-based) */
#define MAX_BIN (((1<<18)-1)/7)
int CAligner::BAMreg2bins(int beg, int end, UINT16 *plist)
{
int i = 0, k;
--end;
plist[i++] = 0;
for (k = 1 + (beg>>26); k <= 1 + (end>>26); ++k) plist[i++] = k;
for (k = 9 + (beg>>23); k <= 9 + (end>>23); ++k) plist[i++] = k;
for (k = 73 + (beg>>20); k <= 73 + (end>>20); ++k) plist[i++] = k;
for (k = 585 + (beg>>17); k <= 585 + (end>>17); ++k) plist[i++] = k;
for (k = 4681 + (beg>>14); k <= 4681 + (end>>14); ++k) plist[i++] = k;
return i;
}

int
CAligner::ReportBAMread(tsReadHit *pReadHit,	// read to report
			int RefID,							// read aligns to this BAM refID (-1) if unaligned
 			bool bPEread,						// true if read is from a PE
			tsBAMalign *pBAMalign)				// BAM alignment to return
{
tsReadHit *pPEReadHit;
char szSEQName[128];
char *pszSEQName;
char *pszPEQName;
char *pszQName;
char *pszRNext;
bool bIsPE2;

int SeqIdx;
etSeqBase Sequence[cMaxFastQSeqLen+1];	// to hold sequence (sans quality scores) for current read

int SumScores;
UINT8 *pSeqVal;
etSeqBase *pSeq;
char *pQScore;
char ExchScore;

int Flags;
int MAPQ;
int GapLen;
int QNameLen;

int CigarIdx;
int Seg0RightTrimLen;
int Seg0LeftTrimLen;
int Seg0Hitlen;
int Seg1Hitlen;
int Seg1RightTrimLen;
int Seg1LeftTrimLen;

int SEStart;
int PEStart;
int SELen;
int PELen;
int PNext;
int TLen;

int LimitGapErrs;
if(pReadHit == NULL)
	return(eBSFerrInternal);

memset(pBAMalign,0,sizeof(tsBAMalign));

strcpy(szSEQName,ReplaceTabs((char *)pReadHit->Read));
pszSEQName = szSEQName;
pszQName = pszSEQName;

if(pReadHit->NAR == eNARAccepted)
	{
	SEStart = AdjAlignStartLoci(&pReadHit->HitLoci.Hit);
	SELen = AdjAlignHitLen(&pReadHit->HitLoci.Hit);
	}
else
	{
	SEStart = 0;
	SELen = 0;
	}

if(bPEread)
	{
	// locate partner read for current read
	// if current read is PE1 then PE2 will be next read, if PE2 then PE1 will be previous read
	bIsPE2 = pReadHit->PairReadID & 0x80000000 ? true : false;	// top bit set if PE2
	if(bIsPE2)	
		pPEReadHit = (tsReadHit *)((UINT8 *)pReadHit - pReadHit->PrevSizeOf);   
	else
		pPEReadHit = (tsReadHit *)((UINT8 *)pReadHit + sizeof(tsReadHit) + pReadHit->ReadLen + pReadHit->DescrLen);

	pszRNext = (char *)"=";
	PNext = 0;
	PEStart = AdjAlignStartLoci(&pPEReadHit->HitLoci.Hit);
	PELen = AdjAlignHitLen(&pPEReadHit->HitLoci.Hit);

	if(SEStart <= PEStart)
		TLen = (PEStart - SEStart) + PELen;
	else
		TLen = (SEStart - PEStart) + SELen;
	}
else    // else SE or if PE then treating as if SE with no partner to be reported
	{
	pszRNext = (char *)"*";
	pszPEQName = NULL;
	bIsPE2 = false;
	PNext = 0;
	TLen = 0;
	}

LimitGapErrs = 10;

pSeqVal = &pReadHit->Read[pReadHit->DescrLen+1];
pSeq = Sequence;
pQScore = (char *)pBAMalign->qual;
SumScores = 0;
for(SeqIdx = 0; SeqIdx < pReadHit->ReadLen; SeqIdx++,pSeq++,pQScore++,pSeqVal++)
	{
	*pSeq = (*pSeqVal & 0x07);		// not interested in any soft masking within a SAM sequence
	SumScores += (*pSeqVal >> 4) & 0x0f;	// any quality scores would be in bits 4..7
	*pQScore = (char)(33 + ((((*pSeqVal >> 4) & 0x0f) * 40))/15);
	}

if(SumScores == 0)		// if there were no associated quality scores
	memset(pBAMalign->qual,0x0ff,pReadHit->ReadLen);
else
	{
	*pQScore = '\0';
	if(pReadHit->NAR == eNARAccepted && pReadHit->HitLoci.Hit.Seg[0].Strand != '+')	// sequence is reversed so quality scores also should be
		{
		pQScore -= 1;
		for(SeqIdx = 0; SeqIdx < (pReadHit->ReadLen/2); SeqIdx++,pQScore--)
			{
			ExchScore = (char)pBAMalign->qual[SeqIdx];
			pBAMalign->qual[SeqIdx] = *pQScore;
			*pQScore = ExchScore;
			}
		}
	}

if(pReadHit->NAR == eNARAccepted)
	{
	Flags = pReadHit->HitLoci.Hit.Seg[0].Strand == '+' ? 0x00 : 0x010;

	if(bPEread)						// if PE
		{
		Flags |= 0x03;				// assumes if PE then both have been mapped
		Flags |= bIsPE2 ? 0x080 : 0x040;
// although the SAM specification treats templates as linear some toolsets are assuming that the PE2 read references the sense of the PE1 read as being the next in the template in the flags field!!!! 
//		if(!bIsPE2) // originally pre-3.7.4 was folllowing specification
		Flags |= pPEReadHit->HitLoci.Hit.Seg[0].Strand == '+' ? 0x00 : 0x020; // post-3.7.4 to enable samtools stats to process sense distributions
		}
	}
else
	Flags = 0x04;			// flags as being unmapped

MAPQ = 255;

CigarIdx = 0;
if(pReadHit->NAR == eNARAccepted)
	{
	Seg0RightTrimLen = pReadHit->HitLoci.Hit.Seg[0].TrimRight;
	Seg0LeftTrimLen = pReadHit->HitLoci.Hit.Seg[0].TrimLeft;
	Seg0Hitlen = AdjHitLen(&pReadHit->HitLoci.Hit.Seg[0]);
	Seg1Hitlen = 0;
	Seg1RightTrimLen = 0;
	Seg1LeftTrimLen = 0;

	if(pReadHit->HitLoci.Hit.Seg[0].Strand == '+')
		{
		if(Seg0LeftTrimLen != 0)
			pBAMalign->cigar[CigarIdx++] = Seg0LeftTrimLen << 4 | 4;		// 'S'
		pBAMalign->cigar[CigarIdx++] = Seg0Hitlen << 4 | 0;					// 'M'
		if(Seg0RightTrimLen > 0)
			pBAMalign->cigar[CigarIdx++] = Seg0RightTrimLen << 4 | 4;		// 'S'
		}
	else
		{
		if(Seg0RightTrimLen != 0)
			pBAMalign->cigar[CigarIdx++] = Seg0RightTrimLen << 4 | 4;		// 'S'
		pBAMalign->cigar[CigarIdx++] = Seg0Hitlen << 4 | 0;					// 'M'	
		if(Seg0LeftTrimLen > 0)
			pBAMalign->cigar[CigarIdx++] = Seg0LeftTrimLen << 4 | 4;		// 'S'
		}

	if(pReadHit->HitLoci.FlagSegs != 0)
		{
		Seg1Hitlen = AdjHitLen(&pReadHit->HitLoci.Hit.Seg[1]);
		Seg1LeftTrimLen = pReadHit->HitLoci.Hit.Seg[1].TrimLeft;
		Seg1RightTrimLen = pReadHit->HitLoci.Hit.Seg[1].TrimRight;
		if(pReadHit->HitLoci.Hit.FlgSplice == 1)  // splice
			{
			GapLen = (int)pReadHit->HitLoci.Hit.Seg[1].MatchLoci - (int)(pReadHit->HitLoci.Hit.Seg[0].MatchLoci + pReadHit->HitLoci.Hit.Seg[0].MatchLen);
			pBAMalign->cigar[CigarIdx++] = GapLen << 4 | 3;		// 'N'
			}
		else   // else if an InDel
			{
			if(pReadHit->HitLoci.Hit.FlgInsert)
				{
				GapLen = pReadHit->ReadLen -
							((pReadHit->HitLoci.Hit.Seg[0].MatchLen - pReadHit->HitLoci.Hit.Seg[0].TrimRight) +
								(pReadHit->HitLoci.Hit.Seg[1].MatchLen  - pReadHit->HitLoci.Hit.Seg[1].TrimLeft));
				if(LimitGapErrs-- > 0 && GapLen <= 0)
					gDiagnostics.DiagOut(eDLInfo,gszProcName,"Check - an apparent insertion GapLen of %d",GapLen);
				pBAMalign->cigar[CigarIdx++] = GapLen << 4 | 1;		// 'I'
				}
			else
				{
				GapLen = (int)pReadHit->HitLoci.Hit.Seg[1].MatchLoci -
					(int)(pReadHit->HitLoci.Hit.Seg[0].MatchLoci + pReadHit->HitLoci.Hit.Seg[0].MatchLen);
				if(LimitGapErrs-- > 0 && GapLen <= 0)
					gDiagnostics.DiagOut(eDLInfo,gszProcName,"Check - an apparent deletion GapLen of %d, Seg[1].MatchLoci:%d, Seg[0].MatchLoci: %d,Seg[0].MatchLen: %d",GapLen,
										pReadHit->HitLoci.Hit.Seg[1].MatchLoci,pReadHit->HitLoci.Hit.Seg[0].MatchLoci,pReadHit->HitLoci.Hit.Seg[0].MatchLen);
				pBAMalign->cigar[CigarIdx++] = abs(GapLen) << 4 | 2;		// 'D'
				}
			}

		if(pReadHit->HitLoci.Hit.Seg[0].Strand == '+')
			{
			if(Seg1LeftTrimLen != 0)
				pBAMalign->cigar[CigarIdx++] = Seg1LeftTrimLen << 4 | 4;		// 'S'
			pBAMalign->cigar[CigarIdx++] = Seg1Hitlen << 4 | 0;		// 'M'
			if(Seg1RightTrimLen > 0)
				pBAMalign->cigar[CigarIdx++] = Seg1RightTrimLen << 4 | 4;		// 'S'
			}
		else
			{
			if(Seg1RightTrimLen != 0)
				pBAMalign->cigar[CigarIdx++] = Seg1RightTrimLen << 4 | 4;		// 'S'
			pBAMalign->cigar[CigarIdx++] = Seg1Hitlen << 4 | 0;		// 'M'
			if(Seg1LeftTrimLen > 0)
				pBAMalign->cigar[CigarIdx++] = Seg1LeftTrimLen << 4 | 4;		// 'S'
			}
		}

	if(bPEread)
		PNext = AdjStartLoci(&pPEReadHit->HitLoci.Hit.Seg[0]);
	else
		PNext = -1;
	QNameLen = 1 + (int)strlen(pszQName);
	pBAMalign->NumReadNameBytes = QNameLen;
	strcpy(pBAMalign->read_name,pszQName);
	pBAMalign->NumCigarBytes = CigarIdx * sizeof(UINT32);
	pBAMalign->flag_nc = Flags << 16 | CigarIdx;
	pBAMalign->refID = RefID;
	pBAMalign->next_refID = PNext == -1 ? -1 : RefID;
	pBAMalign->pos = SEStart;
	pBAMalign->end = SEStart+SELen-1;
	pBAMalign->bin_mq_nl = BAMreg2bin(pBAMalign->pos,SEStart+SELen) << 16 | MAPQ << 8 | QNameLen;
	pBAMalign->next_pos = PNext;
	if(TLen >= 0)
		pBAMalign->tlen = TLen;
	else
		pBAMalign->tlen = -1 * abs(TLen);
	pBAMalign->l_seq = pReadHit->ReadLen;
	pBAMalign->NumSeqBytes = (pReadHit->ReadLen + 1)/2;
	}
else   // treating as being unaligned
	{
    QNameLen = 1 + (int)strlen(pszQName);
	pBAMalign->NumReadNameBytes = QNameLen;
	strcpy(pBAMalign->read_name,pszQName);
	pBAMalign->NumCigarBytes = 4;
	pBAMalign->cigar[0] = pReadHit->ReadLen << 4 | 0;
	pBAMalign->flag_nc = (Flags << 16) | 0x01;
	pBAMalign->refID = -1;
	pBAMalign->next_refID = -1;
	pBAMalign->pos = -1;
	pBAMalign->end = 0;
	pBAMalign->bin_mq_nl = MAPQ << 8 | QNameLen;
	pBAMalign->next_pos = -1;
	pBAMalign->tlen = 0;
	pBAMalign->l_seq = pReadHit->ReadLen;
	pBAMalign->NumSeqBytes = (pReadHit->ReadLen + 1)/2;
	pBAMalign->NumAux = 1;
	pBAMalign->auxData[0].tag[0] = 'Y';
	pBAMalign->auxData[0].tag[1] = 'U';
	pBAMalign->auxData[0].val_type = 'Z';
	pBAMalign->auxData[0].NumVals = 1;
	pBAMalign->auxData[0].array_type = 'A';		// not actually required because value type is a string
	strcpy((char *)pBAMalign->auxData[0].value,m_NARdesc[pReadHit->NAR].pszNAR);
	}

if(pReadHit->NAR == eNARAccepted && pReadHit->HitLoci.Hit.Seg[0].Strand != '+')   // 1.1.6 seems that downstream applications are expecting the read sequences to be reverse complemented if sequence was mapped to Crick strand????
	CSeqTrans::ReverseComplement(pReadHit->ReadLen,Sequence);

UINT8 Byte;
int Ofs;

etSeqBase Base;

Byte = 0;
Ofs = 0;
do
	{
	Base = Sequence[Ofs];
	switch(Base) {       // `=ACMGRSVTWYHKDBN' --> [0, 15];
		case eBaseA:
			Byte |= 1;
			break;
		case eBaseC:
			Byte |= 2;
			break;
		case eBaseG:
			Byte |= 4;
			break;
		case eBaseT:
			Byte |= 8;
			break;
		default:
			Byte |= 15;
		}
	if(!(Ofs & 0x01))
		Byte <<= 4;

	if((Ofs & 0x01) || Ofs == pReadHit->ReadLen-1)
		{
		pBAMalign->seq[Ofs/2] = Byte;
		Byte = 0;
		}
	}
while(++Ofs < pReadHit->ReadLen);

return(eBSFSuccess);
}




// ReplaceTabs
char *
CAligner::ReplaceTabs(char *pszTabTxt) // Inplace replacement of any tabs with a single space char
{
char Chr;
char *pTxt = pszTabTxt;
while(Chr = *pTxt++)
	if(Chr == '\t')
		pTxt[-1]=' ';
return(pszTabTxt);
}

int
CAligner::AppendStr(char *pszBuff,	// write to this buffer
		  char LeadSep,		// if != '\0' then prefix with this separator (usually ',' or '\t')
		  char Quote,		// if != '\0' then quote string with this char (usually single or double quote char)
		  char *pStr,		// '\0' terminated string
		  char TrailSep)	// if != '\0' then suffix with this separator (usually ',' or '\t' or '\n')
{
int Len = 0;
if(LeadSep != '\0')
	{
	*pszBuff++ = LeadSep;
	Len += 1;
	}
if(Quote != '\0')
	{
	*pszBuff++ = Quote;
	Len += 1;
	}
while(*pStr)
	{
	*pszBuff++ = *pStr++;
	Len += 1;
	}
if(Quote != '\0')
	{
	*pszBuff++ = Quote;
	Len += 1;
	}
if(TrailSep != '\0')
	{
	*pszBuff++ = TrailSep;
	Len += 1;
	}
*pszBuff = '\0';
return(Len);
}

int
CAligner::AppendChrs(char *pszBuff,	// write to this buffer
		  char LeadSep,		// if != '\0' then prefix with this separator (usually ',' or '\t')
		  char Quote,		// if != '\0' then quote chars with this char (usually single or double quote char)
		  int NumChrs,		// number of chars to append
		  char *Chrs,		// pts to chars to append
		  char TrailSep)	// if != '\0' then suffix with this separator (usually ',' or '\t' or '\n')
{
int Len = 0;
if(LeadSep != '\0')
	{
	*pszBuff++ = LeadSep;
	Len += 1;
	}
if(Quote != '\0')
	{
	*pszBuff++ = Quote;
	Len += 1;
	}
if(NumChrs)
	{
	while(NumChrs--)
		{
		*pszBuff++ = *Chrs++;
		Len += 1;
		}
	}
if(Quote != '\0')
	{
	*pszBuff++ = Quote;
	Len += 1;
	}
if(TrailSep != '\0')
	{
	*pszBuff++ = TrailSep;
	Len += 1;
	}
*pszBuff = '\0';
return(Len);
}


// very fast version of uitoa
int							// length written
CAligner::AppendUInt(char *pszBuff,	// write to this buffer
		  char LeadSep,		// if > '\0' then prefix with this separator (usually ',' or '\t')
		  UINT32 Value,
		  char TrailSep)	// if > '\0' then suffix with this separator (usually ',' or '\t' or '\n')
{
int Len = 0;
char *pChr;
char *pMark;
char Tmp;
if(LeadSep != '\0')
	{
	*pszBuff++ = LeadSep;
	Len += 1;
	}

if(Value)
	{
	pChr = pszBuff;
	while(Value)
		{
		*pChr++ = '0' + (char)(Value%10);
		Value/=10;
		Len += 1;
		}
	pMark = pChr;
	*pChr-- = '\0';
	while(pszBuff < pChr)
		{
		Tmp = *pChr;
		*pChr-- = *pszBuff;
		*pszBuff++ = Tmp;
		}
	pszBuff = pMark;
	}
else
	{
	Len += 1;
	*pszBuff++ = '0';
	}
if(TrailSep)
	{
	Len += 1;
	*pszBuff++ = TrailSep;
	}
*pszBuff = '\0';
return(Len);
}

// user may be interested in the distribution of the aligner induced substitutions, after any auto-trimming of flanks,
// along the length of the reads and how this distribution relates to the quality scores
int
CAligner::WriteSubDist(tsReadHit *pReadHit)
{
int QScoreIdx;
int NumMSubs;
int SeqIdx;
UINT8 *pSeqVal;
etSeqBase *pAssembSeq;
tsSegLoci *pSeg;

etSeqBase AssembSeq[cMaxFastQSeqLen+1];	// to hold targeted genome assembly sequence

if(m_hStatsFile == -1 || pReadHit == NULL || pReadHit->NAR != eNARAccepted || pReadHit->HitLoci.FlagSegs !=0) // if read was segmented because of InDel or splice junction then can't handle, simply slough
	return(eBSFSuccess);

pSeg = &pReadHit->HitLoci.Hit.Seg[0];
if(pSeg->ChromID == 0)
	return(eBSFSuccess);
if(pSeg->Strand == '\0')	// default strand to be sense if not specified
	pSeg->Strand = '+';

m_pSfxArray->GetSeq(pSeg->ChromID,AdjStartLoci(pSeg),AssembSeq,AdjHitLen(pSeg));	// get sequence for entry starting at offset and of length len
if(pSeg->Strand == '-')
	CSeqTrans::ReverseComplement(AdjHitLen(pSeg),AssembSeq);
		
if(m_MaxAlignLen < pReadHit->ReadLen)
	m_MaxAlignLen = pReadHit->ReadLen;

pSeqVal = &pReadHit->Read[pReadHit->DescrLen+1];
pSeqVal += pSeg->ReadOfs + pSeg->TrimLeft;
pAssembSeq = AssembSeq;
NumMSubs = 0;
for(SeqIdx = pSeg->ReadOfs + pSeg->TrimLeft; SeqIdx < (pReadHit->ReadLen - pSeg->TrimRight); SeqIdx++, pSeqVal++, pAssembSeq++)
	{
	QScoreIdx = (*pSeqVal >> 4) & 0x0f;
	// note that Phred scores were scaled to fit within 4bits (((Phred + 2)* 15)/40) by genreads
	if(QScoreIdx <= 3)		//Phred 0..8?
		QScoreIdx = 0;
	else
		if(QScoreIdx <= 7)	// Phred 9..19?
			QScoreIdx = 1;
		else
			if(QScoreIdx <= 11) // Phred 20..29?
				QScoreIdx = 2;
			else
				QScoreIdx = 3;	// Phred 30+
	m_AlignQSubDist[QScoreIdx][SeqIdx].QInsts += 1;
	if((*pSeqVal & 0x07) != (*pAssembSeq & 0x07))
		{
		m_AlignQSubDist[QScoreIdx][SeqIdx].Subs += 1;
		NumMSubs += 1;
		}
	}
if(m_MaxMSubDist < NumMSubs)
	m_MaxMSubDist = NumMSubs;
m_AlignMSubDist[NumMSubs] += 1;

return(eBSFSuccess);
}


int
CAligner::WriteReadHits(bool bPEProc)		   // true if processing paired ends	
{
const char *pszBsMap;
int LineLen;
char szChromName[128];
int SeqIdx;
etSeqBase ReadSeq[cMaxFastQSeqLen+1];	// to hold sequence (sans quality scores) for current read
etSeqBase AssembSeq[cMaxFastQSeqLen+1];	// to hold targeted genome assembly sequence

UINT8 *pSeqVal;
etSeqBase *pReadSeq;

tsReadHit *pReadHit;
tBSFEntryID PrevTargEntry;

m_MaxAlignLen = 0;

if(m_FMode == eFMbed)
	{
	LineLen = sprintf(m_pszLineBuff,"track type=bed name=\"%s\" description=\"%s\"\n",m_pszTrackTitle,m_pszTrackTitle);
	if(!m_bgzOutFile)
		CUtility::SafeWrite(m_hOutFile,m_pszLineBuff,LineLen);
	else
		CUtility::SafeWrite_gz(m_gzOutFile,m_pszLineBuff,LineLen);

	if(m_hJctOutFile != -1)
		{
		LineLen = sprintf(m_pszLineBuff,"track type=bed name=\"JCT_%s\" description=\"%s\"\n",m_pszTrackTitle,m_pszTrackTitle);
		CUtility::SafeWrite(m_hJctOutFile,m_pszLineBuff,LineLen);
		}

	if(m_hIndOutFile != -1)
		{
		LineLen = sprintf(m_pszLineBuff,"track type=bed name=\"IND_%s\" description=\"%s\"\n",m_pszTrackTitle,m_pszTrackTitle);
		CUtility::SafeWrite(m_hIndOutFile,m_pszLineBuff,LineLen);
		}
	LineLen = 0;
	}

pReadHit = NULL;
LineLen = 0;
PrevTargEntry = 0;
const char *pszAlignType;

bool bPrevInDelSeg = false;
bool bPrevJunctSeg = false;
bool bPrevAlignSeg = false;

while((pReadHit = IterSortedReads(pReadHit))!=NULL)
	{
	if(pReadHit->NAR == eNARAccepted)
		{
		if(pReadHit->HitLoci.Hit.FlgInDel)
			pszAlignType = "ari";
		else
			if(pReadHit->HitLoci.Hit.FlgSplice)
				pszAlignType = "arj";
			else
				pszAlignType = "ar";

		if(pReadHit->HitLoci.FlagIA == 1)
			{
			if(pReadHit->HitLoci.Hit.FlgInDel)
				pszAlignType = "iari";
			else
				if(pReadHit->HitLoci.Hit.FlgSplice)
					pszAlignType = "iarj";
				else
					pszAlignType = "iar";
			}
		else
			{
			if(pReadHit->HitLoci.Hit.FlgInDel)
				pszAlignType = "ari";
			else
				if(pReadHit->HitLoci.Hit.FlgSplice)
					pszAlignType = "arj";
				else
					pszAlignType = "ar";
			}

		tsSegLoci *pSeg;
		int SegIdx;

		pSeg = &pReadHit->HitLoci.Hit.Seg[0];
		if(pSeg->ChromID != (UINT32)PrevTargEntry)
			{
			m_pSfxArray->GetIdentName(pSeg->ChromID,sizeof(szChromName),szChromName);
			PrevTargEntry = pSeg->ChromID;
			}

		if(pSeg->Strand == '\0')	// default strand to be sense if not specified
			pSeg->Strand = '+';

		int Score = (int)min(1000.0,(999 * m_OctSitePrefs[pSeg->Strand == '+' ? 0 : 1][pReadHit->SiteIdx].RelScale));

		if(m_FMode == eFMbed)
			{
			if(pReadHit->HitLoci.FlagSegs==0)
				{
				if(bPrevInDelSeg || bPrevJunctSeg)
					{
					if(LineLen > 0)
						{
						CUtility::SafeWrite(bPrevInDelSeg ? m_hIndOutFile : m_hJctOutFile,m_pszLineBuff,LineLen);
						LineLen = 0;
						}
					bPrevInDelSeg = false;
					bPrevJunctSeg = false;
					}
				pSeg = &pReadHit->HitLoci.Hit.Seg[0];

				LineLen+=sprintf(&m_pszLineBuff[LineLen],"%s\t%d\t%d\t%s\t%d\t%c\n",
					szChromName,AdjStartLoci(pSeg),AdjEndLoci(pSeg) + 1,pszAlignType,Score,pSeg->Strand);
				bPrevAlignSeg = true;
				}
			else // segmented alignment
				{
				if(bPrevAlignSeg)
					{
					if(LineLen > 0)
						{
						if(!m_bgzOutFile)
							CUtility::SafeWrite(m_hOutFile,m_pszLineBuff,LineLen);
						else
							CUtility::SafeWrite_gz(m_gzOutFile,m_pszLineBuff,LineLen);
						LineLen = 0;
						}
					bPrevAlignSeg = false;
					}
				if(pReadHit->HitLoci.Hit.FlgInDel)
					{
					if(bPrevJunctSeg)
						{
						if(LineLen > 0)
							{
							CUtility::SafeWrite(m_hJctOutFile,m_pszLineBuff,LineLen);
							LineLen = 0;
							}
						bPrevJunctSeg = false;
						}
					bPrevInDelSeg = true;
					}
				else
					{
					if(bPrevInDelSeg)
						{
						if(LineLen > 0)
							{
							CUtility::SafeWrite(m_hIndOutFile,m_pszLineBuff,LineLen);
							LineLen = 0;
							}
						bPrevInDelSeg = false;
						}
					bPrevJunctSeg = true;
					}

			    int AjAlignStartLoci;
				int AjAlignEndLoci;

				AjAlignStartLoci = AdjAlignStartLoci(&pReadHit->HitLoci.Hit);
				AjAlignEndLoci = AdjAlignEndLoci(&pReadHit->HitLoci.Hit);
				LineLen+=sprintf(&m_pszLineBuff[LineLen],"%s\t%d\t%d\t%s\t%d\t%c\t%d\t%d\t0\t2\t%d,%d\t0,%d\n",
					szChromName,AjAlignStartLoci,AjAlignEndLoci+1,pszAlignType,Score,pSeg->Strand,AjAlignStartLoci,AjAlignEndLoci+1,
					       AdjHitLen(pSeg),AdjHitLen(&pSeg[1]),AdjStartLoci(&pSeg[1])-AdjStartLoci(pSeg));
				}

			if((cAllocLineBuffSize - LineLen) < 1000)
				{
				if(bPrevAlignSeg)
					{
					if(!m_bgzOutFile)
						CUtility::SafeWrite(m_hOutFile,m_pszLineBuff,LineLen);
					else
						CUtility::SafeWrite_gz(m_gzOutFile,m_pszLineBuff,LineLen);
					LineLen = 0;
					bPrevAlignSeg = false;
					}
				else
					{
					CUtility::SafeWrite(bPrevInDelSeg ? m_hIndOutFile : m_hJctOutFile,m_pszLineBuff,LineLen);
					LineLen = 0;
					bPrevInDelSeg = false;
					bPrevJunctSeg = false;
					}
				}
			continue;
			}

		if(!m_bIsSOLiD && m_FMode >= eFMread)
			{
			pSeqVal = &pReadHit->Read[pReadHit->DescrLen+1];
			pReadSeq = ReadSeq;
			for(SeqIdx = 0; SeqIdx < pReadHit->ReadLen; SeqIdx++,pReadSeq++,pSeqVal++)
				*pReadSeq = (*pSeqVal & 0x07);
			}
		bPrevAlignSeg = true;
		for(SegIdx = 0; SegIdx < 2; SegIdx++)
			{
			pSeg = &pReadHit->HitLoci.Hit.Seg[SegIdx];
			if(pSeg->ChromID == 0)
				continue;
			if(pSeg->Strand == '\0')	// default strand to be sense if not specified
				pSeg->Strand = '+';
			if(pSeg->ChromID != (UINT32)PrevTargEntry)
				{
				m_pSfxArray->GetIdentName(pSeg->ChromID,sizeof(szChromName),szChromName);
				PrevTargEntry = pSeg->ChromID;
				}

			LineLen += AppendUInt(&m_pszLineBuff[LineLen],0,pReadHit->ReadID,',');
			LineLen += AppendStr(&m_pszLineBuff[LineLen],0,'"',(char *)pszAlignType,',');
			LineLen += AppendStr(&m_pszLineBuff[LineLen],0,'"',m_szTargSpecies,',');
			LineLen += AppendStr(&m_pszLineBuff[LineLen],0,'"',szChromName,',');

			LineLen += AppendUInt(&m_pszLineBuff[LineLen],0,AdjStartLoci(pSeg),',');
			LineLen += AppendUInt(&m_pszLineBuff[LineLen],0,AdjEndLoci(pSeg),',');
			LineLen += AppendUInt(&m_pszLineBuff[LineLen],0,AdjHitLen(pSeg),',');
			LineLen += AppendChrs(&m_pszLineBuff[LineLen],0,'"',1,(char *)&pSeg->Strand,',');
			LineLen += AppendUInt(&m_pszLineBuff[LineLen],0,Score,',');
			LineLen += AppendUInt(&m_pszLineBuff[LineLen],0,0,'\0');

			if(m_bBisulfite) {
				switch(pReadHit->HitLoci.Hit.BisBase) {
					case eBaseT:
						pszBsMap = "TC:C";
						break;
					case eBaseA:
						pszBsMap = "AG:T";
						break;
					default:
						pszBsMap = "?:?";
					}
				}
			else
				pszBsMap = "N/A";

			LineLen += AppendUInt(&m_pszLineBuff[LineLen],',',pReadHit->NumReads,',');
			LineLen += AppendUInt(&m_pszLineBuff[LineLen],0,pSeg->TrimMismatches,',');
			LineLen += AppendStr(&m_pszLineBuff[LineLen],0,'"',(char *)pszBsMap,',');
			LineLen += AppendStr(&m_pszLineBuff[LineLen],0,'"',(char *)pReadHit->Read,'\0');

			if(m_FMode >= eFMread)
				LineLen += AppendStr(&m_pszLineBuff[LineLen],',','"',CSeqTrans::MapSeq2Ascii(&ReadSeq[pSeg->ReadOfs+pSeg->TrimLeft],AdjHitLen(pSeg)),0);
			if(m_FMode == eFMmatch || m_FMode == eFMreadmatch)
				{
				m_pSfxArray->GetSeq(pSeg->ChromID,AdjStartLoci(pSeg),AssembSeq,AdjHitLen(pSeg));	// get sequence for entry starting at offset and of length len
				if(pSeg->Strand == '-')
					CSeqTrans::ReverseComplement(AdjHitLen(pSeg),AssembSeq);
				LineLen += AppendStr(&m_pszLineBuff[LineLen],',','"',CSeqTrans::MapSeq2Ascii(AssembSeq,AdjHitLen(pSeg)),0);
				}

			LineLen += AppendStr(&m_pszLineBuff[LineLen],0,0,(char *)"\n",0);
			if(LineLen + ((cMaxFastQSeqLen * 2) + 1024) > cAllocLineBuffSize)
				{
				if(!m_bgzOutFile)
					CUtility::SafeWrite(m_hOutFile,m_pszLineBuff,LineLen);
				else
					CUtility::SafeWrite_gz(m_gzOutFile,m_pszLineBuff,LineLen);
				LineLen = 0;
				}
			}

		// user may be interested in the distribution of the aligner induced substitutions, after any auto-trimming of flanks,
		// along the length of the reads and how this distribution relates to the quality scores
		if(m_hStatsFile != -1)
			WriteSubDist(pReadHit);
		}
	}
if(LineLen)
	{
	if(bPrevAlignSeg)
		{
		if(!m_bgzOutFile)
			CUtility::SafeWrite(m_hOutFile,m_pszLineBuff,LineLen);
		else
			CUtility::SafeWrite_gz(m_gzOutFile,m_pszLineBuff,LineLen);

		LineLen = 0;
		bPrevAlignSeg = false;
		}
	else
		{
		CUtility::SafeWrite(bPrevInDelSeg ? m_hIndOutFile : m_hJctOutFile,m_pszLineBuff,LineLen);
		LineLen = 0;
		bPrevInDelSeg = false;
		bPrevJunctSeg = false;
		}
	}
return(eBSFSuccess);
}


int
CAligner::AddMultiHit(tsReadHit *pReadHit)
{
int NumMultiHits;
size_t HitLen;
tsReadHit *pMultiHit;

HitLen = sizeof(tsReadHit) + pReadHit->ReadLen + pReadHit->DescrLen; 
#ifdef _WIN32
DWORD WaitRslt = WaitForSingleObject(m_hMtxMultiMatches,INFINITE);
#else
pthread_mutex_lock(&m_hMtxMultiMatches);
#endif
if(m_NxtMultiAllOfs + (2 * HitLen) >= m_AllocMultiAllMem)
	{
	size_t memreq = m_AllocMultiAllMem + (cAllocMultihits * HitLen);
#ifdef _WIN32
	pMultiHit = (tsReadHit *) realloc(m_pMultiAll,memreq);
#else
	pMultiHit = (tsReadHit *)mremap(m_pMultiAll,m_AllocMultiAllMem,memreq,MREMAP_MAYMOVE);
	if(pMultiHit == MAP_FAILED)
		pMultiHit = NULL;
#endif
	if(pMultiHit == NULL)
		{
#ifdef _WIN32
ReleaseMutex(m_hMtxMultiMatches);
#else
pthread_mutex_unlock(&m_hMtxMultiMatches);
#endif
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddMultiHit: Memory re-allocation to %lld bytes - %s",memreq,strerror(errno));
		return(eBSFerrMem);
		}
	m_pMultiAll = pMultiHit;
	m_AllocMultiAllMem = memreq;
	}
pMultiHit = (tsReadHit *)((UINT8 *)m_pMultiAll + m_NxtMultiAllOfs);
memmove(pMultiHit,pReadHit,HitLen);
m_NxtMultiAllOfs += HitLen;
m_NumMultiAll += 1;
pMultiHit->ReadID = m_NumMultiAll;
NumMultiHits = m_NumMultiAll;
#ifdef _WIN32
ReleaseMutex(m_hMtxMultiMatches);
#else
pthread_mutex_unlock(&m_hMtxMultiMatches);
#endif
return((int)NumMultiHits);
}

int				// normally NumHits, but will be actual number of hits if unable to accept any of the loci hit because of chromosome filtering
CAligner::WriteHitLoci(tsThreadMatchPars *pThreadPars,tsReadHit *pReadHit,int NumHits,tsHitLoci *pHits)
{
int Rslt;
tsHitLoci *pHit;
tsSegLoci *pSeg;

int ReadHitBuffIdx;						// index into output szReadHits

tBSFEntryID PrevTargEntry;
UINT8 ReadHit[sizeof(tsReadHit) + cMaxDescrLen + cMaxFastQSeqLen + 10];
tsReadHit *pMultiHit;

m_MaxAlignLen = 0;
PrevTargEntry = 0;
ReadHitBuffIdx = 0;

// check if hits are to chromosomes which are to be retained
if(NumHits > 0 && (m_NumExcludeChroms || m_NumIncludeChroms))
	{
	tsHitLoci *pAcceptHit;
	int AcceptHits;
	AcceptHits = 0;
	pHit = pHits;
	pAcceptHit = pHit;
	for(int HitIdx = 0; HitIdx < NumHits; HitIdx++,pHit++)
		{
		pSeg = &pHit->Seg[0];
		if(AcceptThisChromID(pSeg->ChromID))
			{
			AcceptHits += 1;
			if(pAcceptHit != pHit)
				*pAcceptHit++ = *pHit;
			}
		}
	NumHits = AcceptHits;
	if(m_FMode != eFMsamAll && NumHits == 0)
		return(0);
	}

PrevTargEntry = 0;
pHit = pHits;

pMultiHit = (tsReadHit *)&ReadHit;
memmove(&ReadHit,pReadHit,sizeof(tsReadHit) + pReadHit->ReadLen + pReadHit->DescrLen);
if(NumHits == 0)
	{
	pMultiHit->NumHits = 0;
	NumHits = 1;
	}
else
	pMultiHit->NumHits = 1;


for(int HitIdx = 0; HitIdx < NumHits; HitIdx++,pHit++)
	{
	if(pMultiHit->NumHits > 0)
		{
		memmove(&pMultiHit->HitLoci.Hit,pHit,sizeof(tsHitLoci));
		pMultiHit->HitLoci.FlagSegs = (pHit->FlgInDel == 1 || pHit->FlgSplice == 1) ? 1 : 0;
		}
	else
		pMultiHit->HitLoci.FlagSegs = 0;

	if((Rslt = AddMultiHit(pMultiHit)) < eBSFSuccess)
		return(Rslt);
	}

return(NumHits);
}


// OutputSNPs
// Currently can't process for SNPs in InDels or splice junctions
// FDR: Benjamini–Hochberg
// QValue == acceptable FDR e.g. 0.05% or 0.01%
// PValue == 1.0 - Stats.Binomial(TotBasesInColumn,NumBasesInColMismatching+1,GlobalSeqErrRate);
// PValueIdx == sorted index of PValue and locus pairs, 1..k
// Generate PValues for all alignment columns meeting minimum constraints into an array of structures containing column loci and associated PValues
// Sort array of structures ascending on PValues
// Iterate array 1 to k and accept as SNPs those elements with PValues < (PValueIdx/k) * QValue
int
CAligner::OutputSNPs(void)
{
double PValue;
double GlobalSeqErrRate;
double LocalSeqErrRate;
tsSNPcnts *pSNP;
UINT32 Loci;
int Idx;
int NumSNPs;
int TotBases;
char szChromName[cMaxDatasetSpeciesChrom+1];
int LineLen;
double Proportion;
double AdjPValue;
int RelRank;
tsLociPValues *pLociPValues;
size_t memreq;
tsSNPcnts *pSNPWinL;
tsSNPcnts *pSNPWinR;
UINT32 LocalBkgndRateWindow;
UINT32 LocalBkgndRateWinFlank;
UINT32 LocalTotMismatches;
UINT32 LocalTotMatches;
UINT32 LocTMM;
UINT32 LocTM;
CStats Stats;

int CurDiSNPLoci;
int PrevDiSNPLoci;
int DiSNPBuffIdx;
char szDiSNPs[4000];

int CurTriSNPLoci;
int PrevTriSNPLoci;
int FirstTriSNPLoci;
int TriSNPBuffIdx;
char szTriSNPs[4000];

UINT8 SNPFlanks[9];
UINT8 *pSNPFlank;
int SNPFlankIdx;
int SNPCentroidIdx;
UINT8 Base;
tsSNPCentroid *pCentroid;

if(m_pLociPValues == NULL)					// will be NULL first time in
	{
	memreq = cAllocLociPValues * sizeof(tsLociPValues);
#ifdef _WIN32
	m_pLociPValues = (tsLociPValues *) malloc((size_t)memreq);
	if(m_pLociPValues == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"OutputSNPs: Memory allocation of %lld bytes failed",(INT64)memreq);
		return(eBSFerrMem);
		}
#else
	// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
	m_pLociPValues = (tsLociPValues *)mmap(NULL,(size_t)memreq, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
	if(m_pLociPValues == MAP_FAILED)
		{
		m_pLociPValues = NULL;
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"OutputSNPs: Memory allocation of %lld bytes through mmap()  failed",(INT64)memreq,strerror(errno));
		return(eBSFerrMem);
		}
#endif
	m_AllocLociPValuesMem = memreq;
	m_NumLociPValues = 0;
	}

// NOTE: set a floor on the global (whole chromosome) sequencing error rate
GlobalSeqErrRate = max(cMinSeqErrRate,(double)m_pChromSNPs->TotMismatch / (double)(1 + m_pChromSNPs->TotMatch + m_pChromSNPs->TotMismatch));

m_pSfxArray->GetIdentName(m_pChromSNPs->ChromID,sizeof(szChromName),szChromName);

pSNPWinR = &m_pChromSNPs->Cnts[0];
LocalBkgndRateWinFlank = cSNPBkgndRateWindow / 2;
LocalBkgndRateWindow = (LocalBkgndRateWinFlank * 2) + 1;
LocalTotMismatches = 0;
LocalTotMatches = 0;
for(Loci = 0; Loci < min(LocalBkgndRateWindow,m_pChromSNPs->ChromLen); Loci++,pSNPWinR++)
	{
	LocalTotMismatches += pSNPWinR->NumNonRefBases;
	LocalTotMatches += pSNPWinR->NumRefBases;
	}

pLociPValues = m_pLociPValues;
m_NumLociPValues = 0;
pSNP = &m_pChromSNPs->Cnts[0];
pSNPWinL = pSNP;
LineLen = 0;
DiSNPBuffIdx = 0;
TriSNPBuffIdx = 0;
CurDiSNPLoci = 0;
PrevDiSNPLoci = -1;
CurTriSNPLoci = 0;
PrevTriSNPLoci = -1;
FirstTriSNPLoci = -1;
m_MaxDiSNPSep = m_pChromSNPs->MeanReadLen;
for(Loci = 0; Loci < m_pChromSNPs->ChromLen;Loci++, pSNP++)
	{
	// determine background expected error rate from window surrounding the current loci
	if(Loci > LocalBkgndRateWinFlank && (Loci + LocalBkgndRateWinFlank) < m_pChromSNPs->ChromLen)
		{
		// need to ensure that LocalTotMismatches and LocalTotMismatches will never underflow
		if(LocalTotMismatches >= pSNPWinL->NumNonRefBases)
			LocalTotMismatches -= pSNPWinL->NumNonRefBases;
		else
			LocalTotMismatches = 0;
		if(LocalTotMatches >= pSNPWinL->NumRefBases)
			LocalTotMatches -= pSNPWinL->NumRefBases;
		else
			LocalTotMatches = 0;
		LocalTotMismatches += pSNPWinR->NumNonRefBases;
		LocalTotMatches += pSNPWinR->NumRefBases;
		pSNPWinL += 1;
		pSNPWinR += 1;
		}

	TotBases = pSNP->NumNonRefBases + pSNP->NumRefBases;
	if(TotBases > 0)
		{
		m_LociBasesCovered += 1;
		m_LociBasesCoverage += TotBases;
		}

	if(TotBases < m_MinSNPreads)
		continue;

	if(m_hSNPCentsfile != -1)
		{
		// get 4bases up/dn stream from loci with SNP and use these to inc centroid counts of from/to counts
		if(Loci >= cSNPCentfFlankLen && Loci < (m_pChromSNPs->ChromLen - cSNPCentfFlankLen))
			{
			m_pSfxArray->GetSeq(m_pChromSNPs->ChromID,Loci-(UINT32)cSNPCentfFlankLen,SNPFlanks,cSNPCentroidLen);
			pSNPFlank = &SNPFlanks[cSNPCentroidLen-1];
			SNPCentroidIdx = 0;
			for(SNPFlankIdx = 0; SNPFlankIdx < cSNPCentroidLen; SNPFlankIdx++,pSNPFlank--)
				{
				Base = *pSNPFlank & 0x07;
				if(Base > eBaseT)
					break;
				SNPCentroidIdx |= (Base << (SNPFlankIdx * 2));
				}
			if(SNPFlankIdx == cSNPCentroidLen)
				m_pSNPCentroids[SNPCentroidIdx].NumInsts += 1;
			}
		}


	if(pSNP->NumNonRefBases < cMinSNPreads)
		continue; 
	Proportion = (double)pSNP->NumNonRefBases/TotBases;
	if(Proportion < m_SNPNonRefPcnt)	// needs to be at least this proportion of non-ref bases to be worth exploring as being SNP
		continue;

	// needing to allocate more memory? NOTE: allowing small safety margin of 10 tsLociPValues
	if(m_AllocLociPValuesMem  < (sizeof(tsLociPValues) * (m_NumLociPValues + 10)))
		{
		size_t memreq = m_AllocLociPValuesMem + (cAllocLociPValues * sizeof(tsLociPValues));
#ifdef _WIN32
		pLociPValues = (tsLociPValues *) realloc(m_pLociPValues,memreq);
		if(pLociPValues == NULL)
			{
#else
		pLociPValues = (tsLociPValues *)mremap(m_pLociPValues,m_AllocLociPValuesMem,memreq,MREMAP_MAYMOVE);
		if(pLociPValues == MAP_FAILED)
			{
#endif
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"OutputSNPs: Memory reallocation to %lld bytes failed - %s",memreq,strerror(errno));
			return(eBSFerrMem);
			}
		m_pLociPValues = pLociPValues;
		m_AllocLociPValuesMem = memreq;
		pLociPValues = &m_pLociPValues[m_NumLociPValues];
		}



	if(pSNP->NumNonRefBases <= LocalTotMismatches)
		LocTMM = LocalTotMismatches - pSNP->NumNonRefBases;
	else
		LocTMM = 0;

	if(pSNP->NumRefBases < LocalTotMatches)
		LocTM = LocalTotMatches - pSNP->NumRefBases;
	else
		LocTM = 0;


	if((LocTMM + LocTM) == 0)
		LocalSeqErrRate = GlobalSeqErrRate;
	else
		{
		LocalSeqErrRate = (double)LocTMM / (double)(LocTMM + LocTM);
		if(LocalSeqErrRate < GlobalSeqErrRate)
			LocalSeqErrRate = GlobalSeqErrRate;
		}
	if(LocalSeqErrRate > cMaxBkgdNoiseThres)	// don't bother attempting to call if the background is too noisy
		continue;

	// accepting as being a SNP

#ifdef _DISNPS_
	if(!m_bIsSOLiD && m_hDiSNPfile != -1)
		{
		// try to find all reads which are overlapping this SNP plus the prev within 300bp SNP
		tsReadHit *pCurOverlappingRead;
		UINT8 PrevDiSNPBase;
		UINT8 CurDiSNPBase;
		UINT8 FirstTriSNPBase;
		UINT8 PrevTriSNPBase;
		UINT8 CurTriSNPBase;

		int NumHaplotypes;
		int HaplotypeCntThres;
		int NumReadsOverlapping;
		int NumReadsAntisense;
		int DiSNPIdx;
		int DiSNPCnts[64];

		CurDiSNPLoci = Loci;
		if(PrevDiSNPLoci != -1 && CurDiSNPLoci > 0 && ((CurDiSNPLoci - PrevDiSNPLoci) <= m_MaxDiSNPSep))
			{
			NumReadsOverlapping = 0;
			NumReadsAntisense = 0;
			memset(DiSNPCnts,0,sizeof(DiSNPCnts));
			while((pCurOverlappingRead = IterateReadsOverlapping(false,m_pChromSNPs,PrevDiSNPLoci, CurDiSNPLoci)) != NULL)
				{
				// get bases at both SNP loci
				PrevDiSNPBase = AdjAlignSNPBase(pCurOverlappingRead,m_pChromSNPs->ChromID,PrevDiSNPLoci);
				if(PrevDiSNPBase > eBaseT)
					continue;
				CurDiSNPBase = AdjAlignSNPBase(pCurOverlappingRead,m_pChromSNPs->ChromID,CurDiSNPLoci);
				if(CurDiSNPBase > eBaseT)
					continue;
				NumReadsOverlapping += 1;
				if(pCurOverlappingRead->HitLoci.Hit.Seg[0].Strand == '-')
					NumReadsAntisense += 1;
				DiSNPIdx = ((PrevDiSNPBase & 0x03) << 2) | (CurDiSNPBase & 0x03);
				DiSNPCnts[DiSNPIdx] += 1;
				}
			if(NumReadsOverlapping >= m_MinSNPreads)
				{
				NumHaplotypes = 0;       // very simplistic - calling haplotype if at least 3 for any putative DiSNP and more than 5% of read depth
				HaplotypeCntThres = max(3,NumReadsOverlapping / 20);
				for(DiSNPIdx = 0; DiSNPIdx < 16; DiSNPIdx++) 
					if(DiSNPCnts[DiSNPIdx] >= HaplotypeCntThres)
						NumHaplotypes += 1;

				DiSNPBuffIdx += sprintf(&szDiSNPs[DiSNPBuffIdx],"\"%s\",%d,%d,%d,%d,%d",szChromName,PrevDiSNPLoci,CurDiSNPLoci,NumReadsOverlapping,NumReadsAntisense,NumHaplotypes);

				for(DiSNPIdx = 0; DiSNPIdx < 16; DiSNPIdx++)
					DiSNPBuffIdx += sprintf(&szDiSNPs[DiSNPBuffIdx],",%d",DiSNPCnts[DiSNPIdx]);	
				DiSNPBuffIdx += sprintf(&szDiSNPs[DiSNPBuffIdx],"\n");

				if((DiSNPBuffIdx + 200) > sizeof(szDiSNPs))
					{
					CUtility::SafeWrite(m_hDiSNPfile,szDiSNPs,DiSNPBuffIdx);
					DiSNPBuffIdx  = 0;
					}
				}
			}
		
		CurTriSNPLoci = Loci;
		if(FirstTriSNPLoci != -1 && PrevTriSNPLoci > 0 && CurTriSNPLoci > 0 && ((CurTriSNPLoci - FirstTriSNPLoci) <= m_MaxDiSNPSep))
			{
			NumReadsOverlapping = 0;
			NumReadsAntisense = 0;
			memset(DiSNPCnts,0,sizeof(DiSNPCnts));
			while((pCurOverlappingRead = IterateReadsOverlapping(true,m_pChromSNPs,FirstTriSNPLoci, CurTriSNPLoci)) != NULL)
				{
				// get bases at all three SNP loci
				FirstTriSNPBase = AdjAlignSNPBase(pCurOverlappingRead,m_pChromSNPs->ChromID,FirstTriSNPLoci);
				if(FirstTriSNPBase > eBaseT)
					continue;
				PrevTriSNPBase = AdjAlignSNPBase(pCurOverlappingRead,m_pChromSNPs->ChromID,PrevTriSNPLoci);
				if(PrevTriSNPBase > eBaseT)
					continue;
				CurTriSNPBase = AdjAlignSNPBase(pCurOverlappingRead,m_pChromSNPs->ChromID,CurTriSNPLoci);
				if(CurTriSNPBase > eBaseT)
					continue;
				NumReadsOverlapping += 1;
				if(pCurOverlappingRead->HitLoci.Hit.Seg[0].Strand == '-')
					NumReadsAntisense += 1;
				DiSNPIdx = ((FirstTriSNPBase & 0x03) << 4) | ((PrevTriSNPBase & 0x03) << 2) | (CurTriSNPBase & 0x03);
				DiSNPCnts[DiSNPIdx] += 1;
				}
			if(NumReadsOverlapping >= m_MinSNPreads)
				{
				NumHaplotypes = 0;       // very simplistic - calling haplotype if at least 3 for any putative DiSNP and more than 5% of read depth
				HaplotypeCntThres = max(3,NumReadsOverlapping / 20);
				for(DiSNPIdx = 0; DiSNPIdx < 64; DiSNPIdx++) 
					if(DiSNPCnts[DiSNPIdx] >= HaplotypeCntThres)
						NumHaplotypes += 1;

				TriSNPBuffIdx += sprintf(&szTriSNPs[TriSNPBuffIdx],"\"%s\",%d,%d,%d,%d,%d,%d",szChromName,FirstTriSNPLoci,PrevTriSNPLoci,CurTriSNPLoci,NumReadsOverlapping,NumReadsAntisense,NumHaplotypes);

				for(DiSNPIdx = 0; DiSNPIdx < 64; DiSNPIdx++)
					TriSNPBuffIdx += sprintf(&szTriSNPs[TriSNPBuffIdx],",%d",DiSNPCnts[DiSNPIdx]);	
				TriSNPBuffIdx += sprintf(&szTriSNPs[TriSNPBuffIdx],"\n");

				if((TriSNPBuffIdx + 500) > sizeof(szTriSNPs))
					{
					CUtility::SafeWrite(m_hTriSNPfile,szTriSNPs,TriSNPBuffIdx);
					TriSNPBuffIdx  = 0;
					}
				}
			}
		PrevDiSNPLoci = CurDiSNPLoci;
		FirstTriSNPLoci = PrevTriSNPLoci;
		PrevTriSNPLoci = CurTriSNPLoci;
		}
#endif

	// if outputting as marker sequence then get SNP up/dnstream sequence and report
	int MarkerStartLoci;
	int MarkerSeqIdx;
	int AllelicIdx;
	tsSNPcnts *pMarkerBase;
	etSeqBase MarkerSequence[2000];
	etSeqBase *pMarkerSeq;
	int TotMarkerLociBases;
	double MarkerLociBaseProportion;
	int MarkerLen;
	int NumPolymorphicSites;
	if(m_hMarkerFile != -1)			// output marker sequences? 
		{
		// ensure putative marker sequence would be completely contained within the chromosome
		if(Loci < (UINT32)m_Marker5Len)
			continue;
		if((Loci + m_Marker3Len) >= m_pChromSNPs->ChromLen)
			continue;
		TotMarkerLociBases = pSNP->NumNonRefBases + pSNP->NumRefBases;
		MarkerLociBaseProportion = (double)pSNP->NumNonRefBases/TotMarkerLociBases;
		if(MarkerLociBaseProportion < 0.5)
			continue;

		NumPolymorphicSites = 0;
		MarkerLen = 1 + m_Marker5Len + m_Marker3Len;
		MarkerStartLoci = Loci - m_Marker5Len;
		pMarkerSeq = MarkerSequence;
		pMarkerBase = &m_pChromSNPs->Cnts[MarkerStartLoci];
		// check there are alignments covering the complete putative marker sequence
		// and that at any loci covered by the marker has a significant allelic base
		for(MarkerSeqIdx = 0; MarkerSeqIdx < MarkerLen; MarkerSeqIdx++,pMarkerBase++,pMarkerSeq++)
			{
			if((TotMarkerLociBases = pMarkerBase->NumNonRefBases + pMarkerBase->NumRefBases) < m_MinSNPreads)	// must be at least enough reads covering to have confidence in base call
				break;
			MarkerLociBaseProportion = (double)pMarkerBase->NumNonRefBases/TotMarkerLociBases;
			if(MarkerLociBaseProportion <= m_MarkerPolyThres)													// if no more than polymorphic threshold then can simply accept RefBase
				{
				if(MarkerLociBaseProportion > 0.1)
					NumPolymorphicSites += 1;
				*pMarkerSeq = CSeqTrans::MapBase2Ascii(pMarkerBase->RefBase); 
				continue;
				}
			// need to find a major allelic base - base must account for very high proportion of counts
			for(AllelicIdx = 0; AllelicIdx < 5; AllelicIdx++)
				if(pMarkerBase->NonRefBaseCnts[AllelicIdx] > 0 && (MarkerLociBaseProportion = ((double)pMarkerBase->NonRefBaseCnts[AllelicIdx]/TotMarkerLociBases)) >= (1.0 - m_MarkerPolyThres))
					{
					if(MarkerLociBaseProportion < 0.9)
						NumPolymorphicSites += 1;
					*pMarkerSeq = CSeqTrans::MapBase2Ascii(AllelicIdx); 
					break;
					}
			if(AllelicIdx == 5)
				break;
			}
		if(MarkerSeqIdx != MarkerLen)			// only reporting SNPs which are consistent with reported markers
			continue;
		char SNPbase;
		char RefBase;
		RefBase = CSeqTrans::MapBase2Ascii(pSNP->RefBase);
		SNPbase = MarkerSequence[m_Marker5Len];
		if(RefBase == SNPbase)					// double check that the reference base is not being called as being the SNP base
			continue;
		MarkerSequence[MarkerLen] = '\0';

		// accepted marker sequence
		m_MarkerID += 1;
		pLociPValues->MarkerID = m_MarkerID;
		pLociPValues->NumPolymorphicSites = NumPolymorphicSites;
		// >MarkerNNN  Chrom StartLoci|MarkerLen|SNPLoci|Marker5Len,SNPbase|RefBase|NumPolymorphicSites
		LineLen+=sprintf(&m_pszLineBuff[LineLen],">Marker%d %s %d|%d|%d|%d|%c|%c|%d\n%s\n",
										m_MarkerID,szChromName,MarkerStartLoci,MarkerLen,Loci,m_Marker5Len,SNPbase,RefBase,NumPolymorphicSites,MarkerSequence);

		if((LineLen + cMaxSeqLen) > cAllocLineBuffSize)
			{
			CUtility::SafeWrite(m_hMarkerFile,m_pszLineBuff,LineLen);
			LineLen = 0;
			}
		}

	if(m_hMarkerFile == -1)
		{
		pLociPValues->MarkerID = 0;
		pLociPValues->NumPolymorphicSites = 0;
		}
	PValue = 1.0 - Stats.Binomial(TotBases,pSNP->NumNonRefBases,LocalSeqErrRate);
	pLociPValues->PValue = PValue;
	pLociPValues->Loci = Loci;
	pLociPValues->Rank = 0;
	pLociPValues->LocalBkGndSubRate = LocalSeqErrRate;
	pLociPValues->LocalReads = LocTMM + LocTM;
	pLociPValues->LocalSubs = LocTMM;
	pLociPValues->NumReads = TotBases;
	pLociPValues->SNPcnts = *pSNP;
	pLociPValues->NumSubs = pSNP->NumNonRefBases;
	pLociPValues += 1;
	m_NumLociPValues += 1;
	}

if(m_hMarkerFile != -1 && LineLen)
	{
	CUtility::SafeWrite(m_hMarkerFile,m_pszLineBuff,LineLen);
	LineLen = 0;
	}

if(m_hDiSNPfile != -1 && DiSNPBuffIdx > 0)
	{
	CUtility::SafeWrite(m_hDiSNPfile,szDiSNPs,DiSNPBuffIdx);
	DiSNPBuffIdx  = 0;
	}
if(m_hTriSNPfile != -1 && TriSNPBuffIdx > 0)
	{
	CUtility::SafeWrite(m_hTriSNPfile,szTriSNPs,TriSNPBuffIdx);
	DiSNPBuffIdx  = 0;
	}

if(m_NumLociPValues == 0)
	return(eBSFSuccess);
if(m_NumLociPValues > 1)
	m_mtqsort.qsort(m_pLociPValues,m_NumLociPValues,sizeof(tsLociPValues),SortLociPValues);
pLociPValues = m_pLociPValues;
NumSNPs = 0;
for(Idx = 0; Idx < (int)m_NumLociPValues; Idx++,pLociPValues++)
	{
	AdjPValue = ((Idx+1)/(double)m_NumLociPValues) * m_QValue;
	if(pLociPValues->PValue >= AdjPValue)
		break;
	NumSNPs += 1;
	pLociPValues->Rank = Idx + 1;
	}
m_NumLociPValues = NumSNPs;
if(m_NumLociPValues > 1)
	m_mtqsort.qsort(m_pLociPValues,m_NumLociPValues,sizeof(tsLociPValues),SortPValuesLoci);

LineLen = 0;
pLociPValues = m_pLociPValues;
for(Idx = 0; Idx < (int)m_NumLociPValues; Idx++,pLociPValues++)
	{
	m_TotNumSNPs += 1;
	RelRank = max(1,999 - ((999 * pLociPValues->Rank) / m_NumLociPValues));
	if(m_FMode == eFMbed)
		{
		LineLen+=sprintf(&m_pszLineBuff[LineLen],"%s\t%d\t%d\tSNP_%d\t%d\t+\n",
				szChromName,pLociPValues->Loci,pLociPValues->Loci+1,m_TotNumSNPs,RelRank);
		}
	else   // else could be either CSV or VCF
		{
		if(m_bSNPsVCF)
			{
			char szALTs[100];
			char szAltFreq[100];
			int AltOfs;
			int SNPPhred;
			int AltFreqOfs;
			int AltIdx;
			UINT32 CntsThres;		// only reporting cnts which are at least 10% of the highest non-ref base counts. 
                                    // otherwise too many noise cnt bases are reported 

			CntsThres = 0;
			for(AltIdx = 0; AltIdx < eBaseN; AltIdx++)
				{
				if(AltIdx == pLociPValues->SNPcnts.RefBase)
					continue;
				if(pLociPValues->SNPcnts.NonRefBaseCnts[AltIdx] > CntsThres)
					CntsThres = pLociPValues->SNPcnts.NonRefBaseCnts[AltIdx];
				}
			CntsThres = max((CntsThres+5)/10,1);
			AltOfs = 0;
			AltFreqOfs = 0;
			for(AltIdx = 0; AltIdx < eBaseN; AltIdx++)
				{
				if(AltIdx == pLociPValues->SNPcnts.RefBase)
					continue;
				if(pLociPValues->SNPcnts.NonRefBaseCnts[AltIdx] >= CntsThres)
					{
					if(AltOfs > 0)
						{
						szALTs[AltOfs++] = ',';
						szAltFreq[AltFreqOfs++] = ',';
						}
					szALTs[AltOfs++] = CSeqTrans::MapBase2Ascii(AltIdx);
					szALTs[AltOfs] = '\0';
					AltFreqOfs += sprintf(&szAltFreq[AltFreqOfs], "%1.4f",(double)pLociPValues->SNPcnts.NonRefBaseCnts[AltIdx]/pLociPValues->NumReads);
					}
				}
			if(pLociPValues->PValue < 0.0000000001)
				SNPPhred = 100;
			else
				SNPPhred = (int)(0.5 + (10.0*log10(1.0/pLociPValues->PValue)));
			LineLen+=sprintf(&m_pszLineBuff[LineLen],"%s\t%u\tSNP%d\t%c\t%s\t%d\tPASS\tAF=%s;DP=%d\n",
													szChromName,pLociPValues->Loci+1,m_TotNumSNPs,CSeqTrans::MapBase2Ascii(pLociPValues->SNPcnts.RefBase),
													szALTs,SNPPhred,szAltFreq,pLociPValues->NumReads);
			}
		else
			LineLen+=sprintf(&m_pszLineBuff[LineLen],"%d,\"SNP\",\"%s\",\"%s\",%d,%d,1,\"+\",%d,%f,%d,%d,\"%c\",%d,%d,%d,%d,%d,%f,%d,%d,%d,%d\n",
					m_TotNumSNPs,m_szTargSpecies,szChromName,pLociPValues->Loci,pLociPValues->Loci,RelRank,pLociPValues->PValue,
								pLociPValues->NumReads,pLociPValues->NumSubs,
								CSeqTrans::MapBase2Ascii(pLociPValues->SNPcnts.RefBase),
								pLociPValues->SNPcnts.NonRefBaseCnts[0],pLociPValues->SNPcnts.NonRefBaseCnts[1],pLociPValues->SNPcnts.NonRefBaseCnts[2],pLociPValues->SNPcnts.NonRefBaseCnts[3],pLociPValues->SNPcnts.NonRefBaseCnts[4],
								pLociPValues->LocalBkGndSubRate,pLociPValues->LocalReads,pLociPValues->LocalSubs,pLociPValues->MarkerID,pLociPValues->NumPolymorphicSites);
		}
	if((LineLen + cMaxSeqLen + 1) > cAllocLineBuffSize)
		{
		CUtility::SafeWrite(m_hSNPfile,m_pszLineBuff,LineLen);
		LineLen = 0;
		}

	if(m_hSNPCentsfile != -1)
		{
		// get 4bases up/dn stream from loci with SNP and use these to inc centroid counts of from/to counts
		if(pLociPValues->Loci >= cSNPCentfFlankLen && pLociPValues->Loci < (m_pChromSNPs->ChromLen - cSNPCentfFlankLen))
			{
			m_pSfxArray->GetSeq(m_pChromSNPs->ChromID,pLociPValues->Loci-(UINT32)cSNPCentfFlankLen,SNPFlanks,cSNPCentroidLen);
			pSNPFlank = &SNPFlanks[cSNPCentroidLen-1];
			SNPCentroidIdx = 0;
			for(SNPFlankIdx = 0; SNPFlankIdx < cSNPCentroidLen; SNPFlankIdx++,pSNPFlank--)
				{
				Base = *pSNPFlank & 0x07;
				if(Base > eBaseT)
					break;
				SNPCentroidIdx |= (Base << (SNPFlankIdx * 2));
				}
			if(SNPFlankIdx != cSNPCentroidLen)
				continue;

			pSNP = &m_pChromSNPs->Cnts[pLociPValues->Loci];
			pCentroid = &m_pSNPCentroids[SNPCentroidIdx];
			pCentroid->CentroidID = SNPCentroidIdx;
			pCentroid->RefBaseCnt += pSNP->NumRefBases;
			pCentroid->NonRefBaseCnts[0] += pSNP->NonRefBaseCnts[0];
			pCentroid->NonRefBaseCnts[1] += pSNP->NonRefBaseCnts[1];
			pCentroid->NonRefBaseCnts[2] += pSNP->NonRefBaseCnts[2];
			pCentroid->NonRefBaseCnts[3] += pSNP->NonRefBaseCnts[3];
			pCentroid->NonRefBaseCnts[4] += pSNP->NonRefBaseCnts[4];
			pCentroid->NumSNPs += 1;
			}
		}

	}
if(LineLen)
	CUtility::SafeWrite(m_hSNPfile,m_pszLineBuff,LineLen);
return(eBSFSuccess);
}


int
CAligner::ProcessSNPs(void)
{
int Rslt;
int LineLen;

UINT32 SeqIdx;
etSeqBase ReadSeq[cMaxFastQSeqLen+1];	// to hold sequence (sans quality scores) for current read
etSeqBase AssembSeq[cMaxFastQSeqLen+1];	// to hold targeted genome assembly sequence

etSeqBase TargBases[3];
etSeqBase ReadBase;
tsSNPcnts *pSNP;
UINT8 *pSeqVal;
etSeqBase *pReadSeq;
etSeqBase *pAssembSeq;
tsSegLoci *pSeg;
tsReadHit *pReadHit;
tBSFEntryID PrevTargEntry;
UINT32 ChromLen;
UINT32 HitLoci;
UINT32 MatchLen;
UINT32 PrevMMChromID;
UINT32 PrevMMLoci;

if(m_FMode == eFMbed)
	{
	LineLen = sprintf(m_pszLineBuff,"track type=bed name=\"%s_SNPs\" description=\"%s SNPs\"\n",m_pszTrackTitle,m_pszTrackTitle);
	CUtility::SafeWrite(m_hSNPfile,m_pszLineBuff,LineLen);
	LineLen = 0;
	}
else			// else must be either CSV or VCF
	{
	if(m_bSNPsVCF)
		{
		LineLen = sprintf(m_pszLineBuff,"##fileformat=VCFv4.2\n##source=biokangaV%s\n##reference=%s\n##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">\n##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n",
													cpszProgVer,m_pszSfxFile);
		LineLen += sprintf(&m_pszLineBuff[LineLen],"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
		}
	else
		LineLen = sprintf(m_pszLineBuff,"\"SNP_ID\",\"ElType\",\"Species\",\"Chrom\",\"StartLoci\",\"EndLoci\",\"Len\",\"Strand\",\"Rank\",\"PValue\",\"Bases\",\"Mismatches\",\"RefBase\",\"MMBaseA\",\"MMBaseC\",\"MMBaseG\",\"MMBaseT\",\"MMBaseN\",\"BackgroundSubRate\",\"TotWinBases\",\"TotWinMismatches\",\"MarkerID\",\"NumPolymorphicSites\"\n");
	CUtility::SafeWrite(m_hSNPfile,m_pszLineBuff,LineLen);
	LineLen = 0;
	}

if(m_hDiSNPfile != -1)
	{
	int Idx;
	char szDiSNPs[3];
	LineLen = sprintf(m_pszLineBuff,"\"Chrom\",\"SNP1Loci\",\"SNP2Loci\",\"Depth\",\"Antisense\",\"Haplotypes\"");
	for(Idx = 0; Idx < 16; Idx++)
		{
		switch(Idx & 0x03) {
			case 0: szDiSNPs[1] = 'a'; break;
			case 1: szDiSNPs[1] = 'c'; break;
			case 2: szDiSNPs[1] = 'g'; break;
			case 3: szDiSNPs[1] = 't'; break;
			}
		switch((Idx >> 2) & 0x03) {
			case 0: szDiSNPs[0] = 'a'; break;
			case 1: szDiSNPs[0] = 'c'; break;
			case 2: szDiSNPs[0] = 'g'; break;
			case 3: szDiSNPs[0] = 't'; break;
			}
		szDiSNPs[2] = '\0';
		LineLen += sprintf(&m_pszLineBuff[LineLen],",\"%s\"",szDiSNPs);
		}
	LineLen += sprintf(&m_pszLineBuff[LineLen],"\n");
	CUtility::SafeWrite(m_hDiSNPfile,m_pszLineBuff,LineLen);
	LineLen = 0;
	}

if(m_hTriSNPfile != -1)
	{
	int Idx;
	char szTriSNPs[4];
	LineLen = sprintf(m_pszLineBuff,"\"Chrom\",\"SNP1Loci\",\"SNP2Loci\",\"SNP3Loci\",\"Depth\",\"Antisense\",\"Haplotypes\"");
	for(Idx = 0; Idx < 64; Idx++)
		{
		switch(Idx & 0x03) {
			case 0: szTriSNPs[2] = 'a'; break;
			case 1: szTriSNPs[2] = 'c'; break;
			case 2: szTriSNPs[2] = 'g'; break;
			case 3: szTriSNPs[2] = 't'; break;
			}
		switch((Idx >> 2) & 0x03) {
			case 0: szTriSNPs[1] = 'a'; break;
			case 1: szTriSNPs[1] = 'c'; break;
			case 2: szTriSNPs[1] = 'g'; break;
			case 3: szTriSNPs[1] = 't'; break;
			}
		switch((Idx >> 4) & 0x03) {
			case 0: szTriSNPs[0] = 'a'; break;
			case 1: szTriSNPs[0] = 'c'; break;
			case 2: szTriSNPs[0] = 'g'; break;
			case 3: szTriSNPs[0] = 't'; break;
			}
		szTriSNPs[3] = '\0';
		LineLen += sprintf(&m_pszLineBuff[LineLen],",\"%s\"",szTriSNPs);
		}
	LineLen += sprintf(&m_pszLineBuff[LineLen],"\n");
	CUtility::SafeWrite(m_hTriSNPfile,m_pszLineBuff,LineLen);
	LineLen = 0;
	}

// need to check that there are accepted aligned reads to process for SNPs!!!
if(m_hSNPCentsfile != -1)
	{
	if(m_pSNPCentroids == NULL)
		{
		int CentroidIdx;
		tsSNPCentroid *pCentroid;
		if((m_pSNPCentroids = new tsSNPCentroid[cSNPCentroidEls + 16])==NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessSNPs: Memory allocation of %d SNP centroid elements failed",cSNPCentroidEls + 16);
			Reset(false);
			return(eBSFerrMem);
			}
		memset(m_pSNPCentroids,0,sizeof(tsSNPCentroid) * (cSNPCentroidEls+16));
		pCentroid = m_pSNPCentroids;
		for(CentroidIdx = 1; CentroidIdx <= cSNPCentroidEls; CentroidIdx++,pCentroid++)
			pCentroid->CentroidID = CentroidIdx;
		}
	}

pReadHit = NULL;
LineLen = 0;
PrevTargEntry = 0;

m_TotNumSNPs = 0;
while((pReadHit = IterSortedReads(pReadHit))!=NULL)
	{
	if(pReadHit->NAR == eNARAccepted)
		{
		if(pReadHit->HitLoci.Hit.FlgInDel || pReadHit->HitLoci.Hit.FlgSplice)
			continue;

		pSeg = &pReadHit->HitLoci.Hit.Seg[0];
		if(pSeg->ChromID != (UINT32)PrevTargEntry)
			{
			if(m_pChromSNPs != NULL)
				{
				// this is where the SNPs for the previously processed chrom need to saved off as new chrom is about to be processed
				m_pChromSNPs->MeanReadLen = (UINT32)(((m_pChromSNPs->TotReadLen + m_pChromSNPs->NumReads - 1) / m_pChromSNPs->NumReads));
				if((Rslt=OutputSNPs())!=eBSFSuccess)
					{
					Reset(false);
					return(Rslt);
					}
				}

    		PrevTargEntry = pSeg->ChromID;
			ChromLen = m_pSfxArray->GetSeqLen(PrevTargEntry);
			if(m_pChromSNPs == NULL || (m_pChromSNPs != NULL && (ChromLen + 16) > m_pChromSNPs->AllocChromLen))
				{
				if(m_pChromSNPs != NULL)
					{
					delete m_pChromSNPs;
					m_pChromSNPs = NULL;
					}
				size_t AllocSize = sizeof(tsChromSNPs) + ((ChromLen + 16) * sizeof(tsSNPcnts));
				if((m_pChromSNPs = (tsChromSNPs *)new UINT8[AllocSize])==NULL)
					{
					gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessSNPs: Memory allocation of %lld bytes - %s",(INT64)AllocSize,strerror(errno));
					Reset(false);
					return(eBSFerrMem);
					}
				m_pChromSNPs->AllocChromLen = ChromLen + 16;
				}
			memset(&m_pChromSNPs->Cnts,0,((ChromLen + 16) * sizeof(tsSNPcnts)));
			m_pChromSNPs->ChromLen = (UINT32)ChromLen;
			m_pChromSNPs->ChromID = pSeg->ChromID;
			m_pChromSNPs->TotMatch = 0;
			m_pChromSNPs->TotMismatch = 0;
			m_pChromSNPs->MeanReadLen = 0;
			m_pChromSNPs->NumReads = 0;
			m_pChromSNPs->TotReadLen = 0;
			m_pChromSNPs->AdjacentSNPs[0].StartLoci = 0;
			m_pChromSNPs->AdjacentSNPs[0].EndLoci = 0;
			m_pChromSNPs->AdjacentSNPs[0].pFirstIterReadHit = NULL;
			m_pChromSNPs->AdjacentSNPs[0].pPrevIterReadHit = 0;
			m_pChromSNPs->AdjacentSNPs[1].StartLoci = 0;
			m_pChromSNPs->AdjacentSNPs[1].EndLoci = 0;
			m_pChromSNPs->AdjacentSNPs[1].pFirstIterReadHit = NULL;
			m_pChromSNPs->AdjacentSNPs[1].pPrevIterReadHit = 0;
			m_pChromSNPs->pFirstReadHit = NULL;
			m_pChromSNPs->pLastReadHit = NULL;
			PrevTargEntry = m_pChromSNPs->ChromID;
			PrevMMChromID = 0;
			PrevMMLoci = -1;
			}

		// get target genome sequence
		if(m_bIsSOLiD)
			{
			MatchLen = AdjHitLen(pSeg);
			HitLoci = AdjStartLoci(pSeg);
			HitLoci += 1;
			MatchLen -= 1;
			m_pSfxArray->GetColorspaceSeq(pSeg->ChromID,
									HitLoci,
									AssembSeq,MatchLen);	// get colorspace sequence


			}
		else
			{
			// get target assembly sequence for entry starting at offset and of length len
			MatchLen = AdjHitLen(pSeg);
			HitLoci = AdjStartLoci(pSeg);
			m_pSfxArray->GetSeq(pSeg->ChromID,HitLoci,AssembSeq,MatchLen);
			pAssembSeq = AssembSeq;
			for(SeqIdx = 0; SeqIdx < MatchLen; SeqIdx++,pAssembSeq++)
				*pAssembSeq = *pAssembSeq & 0x07;
			}

			// get accepted aligned read sequence
		pSeqVal = &pReadHit->Read[pReadHit->DescrLen+1];
		pSeqVal += pSeg->ReadOfs + pSeg->TrimLeft;
		pReadSeq = ReadSeq;

		if(m_bIsSOLiD)
			{
			// convert read sequence into colorspace
			UINT8 PrvBase = *pSeqVal & 0x07;
			for(SeqIdx = 1; SeqIdx <= MatchLen; SeqIdx++,pReadSeq++,pSeqVal++)
				{
				*pReadSeq = SOLiDmap[PrvBase][pSeqVal[1] & 0x07];
				PrvBase = pSeqVal[1] & 0x07;
				}
			// reverse, not complement, sequence if hit was onto '-' strand
			if(pSeg->Strand == '-')
				CSeqTrans::ReverseSeq(MatchLen,ReadSeq);
			}
		else
			{
			for(SeqIdx = 0; SeqIdx < MatchLen; SeqIdx++,pReadSeq++,pSeqVal++)
				*pReadSeq = *pSeqVal & 0x07;
			if(pSeg->Strand == '-')
				CSeqTrans::ReverseComplement(MatchLen,ReadSeq);
			}

		// double check not about to update snp counts past the expected chrom length
		if((HitLoci + MatchLen) > ChromLen)
			{
			if((MatchLen = (int)ChromLen - HitLoci) < 10)
				continue;
			}

		if(m_pChromSNPs->pFirstReadHit == NULL)
			m_pChromSNPs->pFirstReadHit = pReadHit;
		m_pChromSNPs->pLastReadHit = pReadHit;
		m_pChromSNPs->TotReadLen += MatchLen;
		m_pChromSNPs->NumReads += 1;

		// now iterate read bases and if mismatch then update appropriate counts
		pSNP = &m_pChromSNPs->Cnts[HitLoci];
		pAssembSeq = &AssembSeq[0];
		pReadSeq = &ReadSeq[0];
		UINT32 Loci = HitLoci;
		bool bPairMM = false;
		int SeqMM = 0;
		for(SeqIdx = 0; SeqIdx < MatchLen; SeqIdx++, Loci++,pReadSeq++, pAssembSeq++,pSNP++)
			{
			if(*pAssembSeq >= eBaseN || (m_bIsSOLiD && *pReadSeq >= eBaseN) || *pReadSeq > eBaseN)
				{
				SeqMM += 1;
				continue;
				}

			if(m_bIsSOLiD)		// in colorspace, unpaired mismatches assumed to be sequencer errors and simply sloughed when identifying SNPs
				{
				if(Loci == 0)	// too problematic with SNPs in colorspace at the start of the target sequence, simply slough
					continue;

				if(!bPairMM && *pAssembSeq != *pReadSeq)
					{
					if(SeqIdx < (1+MatchLen))
						{
						if(pReadSeq[1] == pAssembSeq[1])
							{
							pSNP->NumRefBases += 1;
							m_pChromSNPs->TotMatch += 1;
							SeqMM += 1;
							continue;
							}
						}

					// get the previous target sequence base and use this + read colorspace space to derive the mismatch in basespace
					if(pSeg->ChromID != PrevMMChromID || Loci != PrevMMLoci)
						{
						PrevMMChromID = pSeg->ChromID;
						PrevMMLoci = Loci;
						m_pSfxArray->GetSeq(pSeg->ChromID,Loci-1,&TargBases[0],2);
						if(TargBases[0] > eBaseN)
							TargBases[0] = eBaseN;
						if(TargBases[1] > eBaseN)
							TargBases[1] = eBaseN;
						pSNP->RefBase = TargBases[1];
						}

					ReadBase = *pReadSeq;
					if(ReadBase > eBaseT)
						ReadBase = eBaseN;

					if(SeqMM == 0)
						ReadBase = SOLiDmap[TargBases[0]][ReadBase];
					else
						ReadBase = eBaseN;

					// sometimes it seems that a colorspace read may have had a sequencing error earlier in the read
					// or some mismatch such that the current loci mismatches in colorspace but matches in basespace
					// these strange bases are treated as though they are undefined and accrue counts as being eBaseN's
					if(ReadBase == pSNP->RefBase)
						ReadBase = eBaseN;
					pSNP->NonRefBaseCnts[ReadBase] += 1;
					pSNP->NumNonRefBases += 1;
					m_pChromSNPs->TotMismatch += 1;
					bPairMM = true;
					SeqMM += 1;
					}
				else
					{
					pSNP->NumRefBases += 1;
					m_pChromSNPs->TotMatch += 1;
					bPairMM = false;
					SeqMM = 0;
					}
				}
			else				// in basespace any mismatch is counted as a NonRefCnt
				{
				ReadBase = *pReadSeq & 0x07;
				TargBases[0] = *pAssembSeq & 0x07;

				pSNP->RefBase = TargBases[0];
				if(TargBases[0] == ReadBase)
					{
					pSNP->NumRefBases += 1;
					m_pChromSNPs->TotMatch += 1;
					}
				else
					{
					if(ReadBase > eBaseT)
						ReadBase = eBaseN;
					pSNP->NonRefBaseCnts[ReadBase] += 1;
					pSNP->NumNonRefBases += 1;
					m_pChromSNPs->TotMismatch += 1;
					}
				}
			}
		}
	}
if((Rslt=OutputSNPs())!=eBSFSuccess)
	{
	Reset(false);
	return(Rslt);
	}
if(m_hSNPfile != -1)
	{
#ifdef _WIN32
	_commit(m_hSNPfile);
#else
	fsync(m_hSNPfile);
#endif
	close(m_hSNPfile);
	m_hSNPfile = -1;
	}

if(m_hDiSNPfile != -1)
	{
#ifdef _WIN32
	_commit(m_hDiSNPfile);
#else
	fsync(m_hDiSNPfile);
#endif
	close(m_hDiSNPfile);
	m_hDiSNPfile = -1;
	}

if(m_hTriSNPfile != -1)
	{
#ifdef _WIN32
	_commit(m_hTriSNPfile);
#else
	fsync(m_hTriSNPfile);
#endif
	close(m_hTriSNPfile);
	m_hTriSNPfile = -1;
	}

if(m_hSNPCentsfile != -1)
	{
	// report on the SNP centroid distributions
	int SNPCentroidIdx;
	int CentroidSeq;
	UINT8 Bases[cSNPCentroidLen];
	int BaseIdx;
	tsSNPCentroid *pCentroid;
	char szCentroids[4096];
	int BuffIdx;

	BuffIdx = sprintf(szCentroids,"\"CentroidID\",\"Seq\",\"NumInsts\",\"NumSNPs\",\"RefBase\",\"RefBaseCnt\",\"BaseA\",\"BaseC\",\"BaseG\",\"BaseT\",\"BaseN\"\n");
	pCentroid = m_pSNPCentroids;
	for(SNPCentroidIdx = 0; SNPCentroidIdx < cSNPCentroidEls; SNPCentroidIdx++, pCentroid++)
		{
		CentroidSeq = SNPCentroidIdx;
		for(BaseIdx = cSNPCentroidLen-1; BaseIdx >= 0; BaseIdx--)
			{
			Bases[BaseIdx] = CentroidSeq & 0x03;
			CentroidSeq >>= 2;
			}

		BuffIdx += sprintf(&szCentroids[BuffIdx],"%d,\"%s\",%d,%d,\"%c\",%d,%d,%d,%d,%d,%d\n",
							SNPCentroidIdx+1,CSeqTrans::MapSeq2Ascii(Bases,cSNPCentroidLen),pCentroid->NumInsts,pCentroid->NumSNPs,CSeqTrans::MapBase2Ascii(Bases[cSNPCentfFlankLen]),
							pCentroid->RefBaseCnt,pCentroid->NonRefBaseCnts[0],pCentroid->NonRefBaseCnts[1],pCentroid->NonRefBaseCnts[2],pCentroid->NonRefBaseCnts[3],pCentroid->NonRefBaseCnts[4]);


		if(BuffIdx + 200 > sizeof(szCentroids))
			{
			CUtility::SafeWrite(m_hSNPCentsfile,szCentroids,BuffIdx);
			BuffIdx = 0;
			}
		}
	if(BuffIdx)
		CUtility::SafeWrite(m_hSNPCentsfile,szCentroids,BuffIdx);
	}

if(m_hSNPCentsfile != -1)
	{
#ifdef _WIN32
	_commit(m_hSNPCentsfile);
#else
	fsync(m_hSNPCentsfile);
#endif
	close(m_hSNPCentsfile);
	m_hSNPCentsfile = -1;
	}

if(m_hMarkerFile != -1)
	{
#ifdef _WIN32
	_commit(m_hMarkerFile);
#else
	fsync(m_hMarkerFile);
#endif
	close(m_hMarkerFile);
	m_hMarkerFile = -1;
	}

if(m_pChromSNPs != NULL)
	{
	delete m_pChromSNPs;
	m_pChromSNPs = NULL;
	}
return(eBSFSuccess);
}

//--- the following function ProcessSiteProbabilites() is targeted for use in RNA-seq processing
// which should help in identifying differentially expressed transcripts
int
CAligner::ProcessSiteProbabilites(int RelSiteStartOfs)	// offset the site octamer by this relative start offset (read start base == 0)
{
int LineLen;

int SeqIdx;
int SiteIdx;
etSeqBase AssembSeq[cMaxFastQSeqLen+1];	// to hold genome assembly sequence for current read loci

tsOctSitePrefs *pSitePrefs;
tsOctSitePrefs *pSitePref;
etSeqBase *pAssembSeq;
tsReadHit *pReadHit;
tsSegLoci *pSeg;
tBSFEntryID PrevTargEntry;
UINT32 HitLoci;
UINT32 PrevLoci;
UINT32 MatchLen;
UINT32 CurChromLen;

int TotOccs;
int TotSites;

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing for alignment site probabilities...");

TotOccs = 0;
TotSites = 0;
memset(m_OctSitePrefs,0,sizeof(m_OctSitePrefs));
pSitePrefs = m_OctSitePrefs[0];
pSitePref = m_OctSitePrefs[1];
for(SiteIdx = 0; SiteIdx < cNumOctamers; SiteIdx++,pSitePref++,pSitePrefs++)
	{
	pSitePref->Octamer = SiteIdx;
	pSitePrefs->Octamer = SiteIdx;
	}
pReadHit = NULL;
LineLen = 0;
PrevTargEntry = 0;
PrevLoci = -1;
CurChromLen = -1;
// iterate all accepted as aligned reads
// reads are assumed to have been sorted by chrom and start loci
while((pReadHit = IterSortedReads(pReadHit))!=NULL)
	{
	if(pReadHit->NAR == eNARAccepted)
		{
		if(pReadHit->HitLoci.Hit.FlgInDel || pReadHit->HitLoci.Hit.FlgSplice)
			continue;

		pSeg = &pReadHit->HitLoci.Hit.Seg[0];
		if(pSeg->ChromID != (UINT32)PrevTargEntry)
			{
    		PrevTargEntry = pSeg->ChromID;
			PrevLoci = -1;
			CurChromLen = m_pSfxArray->GetSeqLen(pSeg->ChromID);
			}

		// get target assembly sequence for entry starting at MatchLoci, offset by RelSiteStartOfs, and of octamer length
		MatchLen = (UINT32)pSeg->MatchLen;
		HitLoci = (UINT32)pSeg->MatchLoci;

		if(pSeg->Strand == '+')
			HitLoci += RelSiteStartOfs;
		else
			{
			HitLoci += pSeg->MatchLen - 1;
			HitLoci -= RelSiteStartOfs;
			HitLoci -= 7;
			}

		if(HitLoci < 0)						// force octamer to start/end within the chrom sequence
			HitLoci = 0;
		else
			if((HitLoci + 8) >= CurChromLen)
				HitLoci = CurChromLen - 9;

		m_pSfxArray->GetSeq(pSeg->ChromID,HitLoci,AssembSeq,8);
		pAssembSeq = AssembSeq;
		for(SeqIdx = 0; SeqIdx < 8; SeqIdx++,pAssembSeq++)
			*pAssembSeq = *pAssembSeq & 0x07;

		if(pSeg->Strand == '-')
			{
			pSitePrefs = m_OctSitePrefs[1];
			CSeqTrans::ReverseComplement(8,AssembSeq);
			}
		else
			pSitePrefs = m_OctSitePrefs[0];
		pAssembSeq = &AssembSeq[0];

        // generate site index
		SiteIdx = 0;
		for(SeqIdx = 0; SeqIdx < 8; SeqIdx++)
			{
			if(*pAssembSeq > eBaseT)
				break;
			SiteIdx <<= 2;
			SiteIdx |= *pAssembSeq++;
			}
		if(SeqIdx != 8)
			continue;
		pReadHit->SiteIdx = SiteIdx;
		pSitePref = &pSitePrefs[SiteIdx];

		pSitePref->NumOccs += 1;
		TotOccs += 1;
		if(HitLoci != PrevLoci)
			{
			pSitePref->NumSites += 1;
			PrevLoci = HitLoci;
			TotSites += 1;
			}
		}
	}

// now to generate the relative abundance scores
// firstly determine the mean number of read alignment occurrences for each octamer
// find the top 64 (~ 0.1%) with highest mean occurrences
// all octamers are then scaled to this mean of the top 0.1%
pSitePrefs = m_OctSitePrefs[0];
pSitePref = m_OctSitePrefs[1];
for(SiteIdx = 0; SiteIdx < cNumOctamers; SiteIdx++,pSitePrefs++,pSitePref++)
	{
	if(pSitePrefs->NumSites >= 1)
		pSitePrefs->RelScale = (double)pSitePrefs->NumOccs/pSitePrefs->NumSites;
	else
		pSitePrefs->RelScale = 0.0;
	if(pSitePref->NumSites >= 1)
		pSitePref->RelScale = (double)pSitePref->NumOccs/pSitePref->NumSites;
	else
		pSitePref->RelScale = 0;
	}

// now sort ascending  by RelScale so top 0.1% can be determined and their mean determined
m_mtqsort.qsort(m_OctSitePrefs[0],cNumOctamers,sizeof(tsOctSitePrefs),SortSiteRelScale);
m_mtqsort.qsort(m_OctSitePrefs[1],cNumOctamers,sizeof(tsOctSitePrefs),SortSiteRelScale);

double TopWatsonMean = 0.0;
double TopCrickMean = 0.0;
pSitePrefs = &m_OctSitePrefs[0][0x0ffc0];
pSitePref = &m_OctSitePrefs[1][0x0ffc0];
for(SiteIdx = 0x0ffc0; SiteIdx < cNumOctamers; SiteIdx++,pSitePrefs++,pSitePref++)
	{
	TopWatsonMean += pSitePrefs->RelScale;
	pSitePrefs->RelScale = 1.0;
	TopCrickMean +=	pSitePref->RelScale;
	pSitePref->RelScale = 1.0;
	}
TopWatsonMean /= 64;
TopCrickMean /= 64;

// top means known now set normalisation scale factors
pSitePrefs = m_OctSitePrefs[0];
pSitePref = m_OctSitePrefs[1];
for(SiteIdx = 0; SiteIdx < 0x0ffc0; SiteIdx++,pSitePrefs++,pSitePref++)
	{
	if(pSitePrefs->RelScale > 0.0)
		pSitePrefs->RelScale = max(0.0001,pSitePrefs->RelScale / TopWatsonMean);
	if(pSitePref->RelScale)
		pSitePref->RelScale = max(0.0001,pSitePref->RelScale / TopCrickMean);
	}

// restore m_OctSitePrefs to original octamer ascending order
m_mtqsort.qsort(m_OctSitePrefs[0],cNumOctamers,sizeof(tsOctSitePrefs),SortSiteRelOctamer);
m_mtqsort.qsort(m_OctSitePrefs[1],cNumOctamers,sizeof(tsOctSitePrefs),SortSiteRelOctamer);

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Completed alignment site probabilities");
return(eBSFSuccess);
}



char *
CAligner::Octamer2Txt(int Octamer) // Report on site octamer site preferencing distribution
{
static char szOctamer[0x09];	// to contain '\0' terminated octamer as nucleotide bases 'a'..'t'
char *pChr;
int Idx;
pChr = &szOctamer[8];
*pChr-- = '\0';
for(Idx = 0; Idx < 8; Idx++,pChr--)
	{
	switch(Octamer & 0x03) {
		case 0:
			*pChr = 'a';
			break;
		case 1:
			*pChr = 'c';
			break;
		case 2:
			*pChr = 'g';
			break;
		case 3:
			*pChr = 't';
			break;
		}
	Octamer >>= 2;
	}
return(szOctamer);
return(NULL);
}

int
CAligner::WriteSitePrefs(void)
{
char szBuff[0x03fff];
int BuffIdx;
int SiteIdx;
tsOctSitePrefs *pSite;
BuffIdx = sprintf(szBuff,"\"Id\",\"Strand\",\"Octamer\",\"TotalHits\",\"UniqueLoci\",\"RelScale\"\n");
pSite = m_OctSitePrefs[0];
for(SiteIdx = 0; SiteIdx < 0x0ffff; SiteIdx++,pSite++)
	{
	BuffIdx += sprintf(&szBuff[BuffIdx],"%d,\"+\",\"%s\",%d,%d,%1.3f\n",SiteIdx+1,Octamer2Txt(pSite->Octamer),pSite->NumOccs,pSite->NumSites,pSite->RelScale);
	if(BuffIdx + 200 > sizeof(szBuff))
		{
		CUtility::SafeWrite(m_hSitePrefsFile,szBuff,BuffIdx);
		BuffIdx = 0;
		}
	}
pSite = m_OctSitePrefs[1];
for(SiteIdx = 0; SiteIdx < 0x0ffff; SiteIdx++,pSite++)
	{
	BuffIdx += sprintf(&szBuff[BuffIdx],"%d,\"-\",\"%s\",%d,%d,%1.3f\n",SiteIdx+1,Octamer2Txt(pSite->Octamer),pSite->NumOccs,pSite->NumSites,pSite->RelScale);
	if(BuffIdx + 200 > sizeof(szBuff))
		{
		CUtility::SafeWrite(m_hSitePrefsFile,szBuff,BuffIdx);
		BuffIdx = 0;
		}
	}

if(BuffIdx > 0)
	{
	CUtility::SafeWrite(m_hSitePrefsFile,szBuff,BuffIdx);
	BuffIdx = 0;
	}
return(eBSFSuccess);
}

int
CAligner::LoadReads(char *pszRdsFile)	// file containing preprocessed reads (genreads output)
{
int Rslt;
size_t memreq;
int RdLen;
tsRawReadV5 *pReadV5;					// current preprocessed read being processed if V5
tsRawReadV6 *pReadV6;					// current preprocessed read being processed if V6

UINT8 *pReadBuff;						// alloc'd to buffer the preprocessed reads from disk
tsReadHit *pReadHit;					// current read hit
int BuffLen;
int BuffOfs;
char *pChr;
char Chr;
// check if file name contains any wildcard chars
pChr = pszRdsFile;
while(Chr = *pChr++)
	if(Chr == '*' || Chr == '?' || Chr == '[' || Chr == ']')
		return(eBSFerrNotBioseq);

#ifdef _WIN32
m_hInFile = open(pszRdsFile, O_READSEQ ); // file access is normally sequential..
#else
m_hInFile = open64(pszRdsFile, O_READSEQ ); // file access is normally sequential..
#endif

if(m_hInFile == -1)					// check if file open succeeded
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open %s - %s",pszRdsFile,strerror(errno));
	return(eBSFerrOpnFile);
	}

// expecting a preprocessed .rds file as input, header processing will confirm!
if((Rslt=Disk2Hdr(pszRdsFile))!=eBSFSuccess)
	{
	close(m_hInFile);
	m_hInFile = -1;
	return((teBSFrsltCodes)Rslt);
	}

if((pReadBuff = new UINT8 [cRdsBuffAlloc])==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadReads: Memory allocation of %d bytes - %s",cRdsBuffAlloc,strerror(errno));
	close(m_hInFile);
	m_hInFile = -1;
	return(eBSFerrMem);
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Reads file '%s' generator version: %d",pszRdsFile,m_FileHdr.Version);

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Contains total of %d reads from %d source reads with duplicate sequences %s, PMode was %d, Quality was %d, %d bases 5' trimmed, %d bases 3' trimmed",
			m_FileHdr.NumRds,m_FileHdr.OrigNumReads,m_FileHdr.FlagsK ? "retained":"removed", m_FileHdr.PMode,m_FileHdr.QMode,m_FileHdr.Trim5,m_FileHdr.Trim3);
if(m_FileHdr.FlagsPR)
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Reads were processed as %s",m_FileHdr.FlagsPR ? "paired end reads" : "single end reads");

gDiagnostics.DiagOutMsgOnly(eDLInfo,"Reads were processed from %d files",m_FileHdr.NumFiles);
char *pszSrcFile = (char *)m_FileHdr.FileNames;
for(BuffOfs=0;BuffOfs<m_FileHdr.NumFiles;BuffOfs++)
	{
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Source file: '%s'",pszSrcFile);
	pszSrcFile += strlen(pszSrcFile) + 1;
	}

// ensure there is at least one read...
if(m_FileHdr.NumRds == 0 || m_FileHdr.TotReadsLen == 0)
	{
	gDiagnostics.DiagOut(eDLWarn,gszProcName,"Nothing to do, '%s' contains no reads...",pszRdsFile);
	delete pReadBuff;
	close(m_hInFile);
	m_hInFile = -1;
	return(eBSFerrOpnFile);
	}

// initial allocation of memory to hold all pre-processed reads plus a little safety margin (10000) bytes)
memreq = ((size_t)m_FileHdr.NumRds * sizeof(tsReadHit)) + (size_t)m_FileHdr.TotReadsLen + 10000;
AcquireSerialise();
AcquireLock(true);
#ifdef _WIN32
m_pReadHits = (tsReadHit *) malloc(memreq);	// initial and perhaps the only allocation

if(m_pReadHits == NULL)
	{
	ReleaseLock(true);
	ReleaseSerialise();
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadReads: Memory allocation of %lld bytes - %s",(INT64)memreq,strerror(errno));
	delete pReadBuff;
	close(m_hInFile);
	m_hInFile = -1;
	return(eBSFerrMem);
	}
#else
// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
m_pReadHits = (tsReadHit *)mmap(NULL,memreq, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
if(m_pReadHits == MAP_FAILED)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadReads: Memory allocation of %lld bytes through mmap()  failed - %s",(INT64)memreq,strerror(errno));
	m_pReadHits = NULL;
	ReleaseLock(true);
	ReleaseSerialise();
	delete pReadBuff;
	close(m_hInFile);
	m_hInFile = -1;
	return(eBSFerrMem);
	}
#endif

m_AllocdReadHitsMem = memreq;
m_UsedReadHitsMem = 0;
m_FinalReadID = 0;
m_NumReadsLoaded = 0;
m_NumDescrReads = 0;
ReleaseLock(true);
ReleaseSerialise();

// iterate each read sequence starting from the first
lseek(m_hInFile,(long)m_FileHdr.RdsOfs,SEEK_SET);
BuffLen = 0;
BuffOfs = 0;

int SizeOfRawRead = m_FileHdr.Version == 5 ? sizeof(tsRawReadV5) : sizeof(tsRawReadV6);
int CurReadLen;
int CurDescrLen;

while((RdLen = read(m_hInFile,&pReadBuff[BuffLen],cRdsBuffAlloc - BuffLen)) > 0)
	{
	BuffLen += RdLen;
	BuffOfs = 0;
	while((BuffLen - BuffOfs) >=  SizeOfRawRead)
		{
		if(m_FileHdr.Version == 5)
			{
			pReadV5 = (tsRawReadV5 *)&pReadBuff[BuffOfs];
			CurDescrLen = (int)pReadV5->DescrLen;
			CurReadLen = (int)pReadV5->ReadLen;
			}
		else
			{
			pReadV6 = (tsRawReadV6 *)&pReadBuff[BuffOfs];
			CurDescrLen = (int)pReadV6->DescrLen;
			CurReadLen = (int)pReadV6->ReadLen;
			}

		if((int)(CurDescrLen + CurReadLen + SizeOfRawRead) > (BuffLen - BuffOfs))
			break;
		BuffOfs += CurDescrLen + CurReadLen + SizeOfRawRead;
			// pReadV5/V6 now pts at a read
			// shouldn't but is there a need to allocate more memory?
		if(m_UsedReadHitsMem + (sizeof(tsReadHit) + CurReadLen + CurDescrLen) >= (m_AllocdReadHitsMem - 5000))
			{
			AcquireSerialise();
			AcquireLock(true);
			memreq = m_AllocdReadHitsMem + ((sizeof(tsReadHit) + (size_t)cDfltReadLen) * cReadsHitReAlloc);
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"LoadReads: Needing memory re-allocation to %lld bytes from %lld",(INT64)m_AllocdReadHitsMem,(INT64)memreq);

#ifdef _WIN32
			pReadHit = (tsReadHit *) realloc(m_pReadHits,memreq);
#else
			pReadHit = (tsReadHit *)mremap(m_pReadHits,m_AllocdReadHitsMem,memreq,MREMAP_MAYMOVE);
			if(pReadHit == MAP_FAILED)
				pReadHit = NULL;
#endif
			if(pReadHit == NULL)
				{
				ReleaseLock(true);
				ReleaseSerialise();
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadReads: Memory re-allocation to %lld bytes - %s",(INT64)(m_AllocdReadHitsMem + ((sizeof(tsReadHit) + cDfltReadLen)*cReadsHitReAlloc)),strerror(errno));
				delete pReadBuff;
				close(m_hInFile);
				m_hInFile = -1;
				return(eBSFerrMem);
				}
			m_pReadHits = pReadHit;
			m_AllocdReadHitsMem = memreq;
			ReleaseLock(true);
			ReleaseSerialise();
			}

		pReadHit = (tsReadHit *)((UINT8 *)m_pReadHits + m_UsedReadHitsMem);
		m_UsedReadHitsMem += sizeof(tsReadHit) + CurReadLen + CurDescrLen;
		memset(pReadHit,0,sizeof(tsReadHit));
		pReadHit->HitLoci.Hit.Seg[0].Strand = '?';
		pReadHit->ReadLen = CurReadLen;
		pReadHit->DescrLen = CurDescrLen;

		if(m_FileHdr.Version == 5)
			{
			pReadHit->ReadID = pReadV5->ReadID;
			pReadHit->PairReadID = pReadV5->PairReadID;
			pReadHit->NumReads = pReadV5->NumReads;

			memmove(pReadHit->Read,pReadV5->Read,CurDescrLen + 1 + CurReadLen);
			m_FinalReadID = pReadV5->ReadID;
			}
		else
			{
			pReadHit->ReadID = pReadV6->ReadID;
			pReadHit->PairReadID = pReadV6->PairReadID;
			pReadHit->NumReads = pReadV6->NumReads;
			memmove(pReadHit->Read,pReadV6->Read,CurDescrLen + 1 + CurReadLen);
			m_FinalReadID = pReadV6->ReadID;
			}

		m_NumDescrReads += 1;

		// processing threads are only updated with actual number of loaded reads every 100000 reads so as
		// to minimise disruption to the actual aligner threads which will also be serialised through m_hMtxIterReads
		if(m_NumDescrReads > 0 && !(m_NumDescrReads % 100000))
			{
			AcquireSerialise();
			m_FinalReadID = m_NumDescrReads;
			m_NumReadsLoaded = m_NumDescrReads;
			ReleaseSerialise();
			}
		}
	BuffLen -= BuffOfs;
	if(BuffLen)
		memmove(pReadBuff,&pReadBuff[BuffOfs],BuffLen);
	}
delete pReadBuff;
close(m_hInFile);
m_hInFile = -1;
if(m_NumDescrReads != m_NumReadsLoaded)
	{
	AcquireSerialise();
	m_FinalReadID = m_NumDescrReads;
	m_NumReadsLoaded = m_NumDescrReads;
	ReleaseSerialise();
	}
return(m_NumDescrReads);
}

int
CAligner::CreateMutexes(void)
{
if(m_bMutexesCreated)
	return(eBSFSuccess);

#ifdef _WIN32
InitializeSRWLock(&m_hRwLock);
#else
if(pthread_rwlock_init (&m_hRwLock,NULL)!=0)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: unable to create rwlock");
	return(eBSFerrInternal);
	}
#endif

#ifdef _WIN32
if((m_hMtxIterReads = CreateMutex(NULL,false,NULL))==NULL)
	{
#else
if(pthread_mutex_init (&m_hMtxIterReads,NULL)!=0)
	{
	pthread_rwlock_destroy(&m_hRwLock);
#endif
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: unable to create mutex");
	return(eBSFerrInternal);
	}

#ifdef _WIN32
if((m_hMtxMHReads = CreateMutex(NULL,false,NULL))==NULL)
	{
#else
if(pthread_mutex_init (&m_hMtxMHReads,NULL)!=0)
	{
#endif
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: unable to create mutex");
#ifdef _WIN32
	CloseHandle(m_hMtxIterReads);
#else
	pthread_rwlock_destroy(&m_hRwLock);
	pthread_mutex_destroy(&m_hMtxIterReads);
#endif
	return(eBSFerrInternal);
	}
if(m_MLMode != eMLdefault)
	{
#ifdef _WIN32
	if((m_hMtxMultiMatches = CreateMutex(NULL,false,NULL))==NULL)
		{
		CloseHandle(m_hMtxIterReads);
		CloseHandle(m_hMtxMHReads);
#else
	if(pthread_mutex_init (&m_hMtxMultiMatches,NULL)!=0)
		{
		pthread_mutex_destroy(&m_hMtxIterReads);
		pthread_mutex_destroy(&m_hMtxMHReads);
		pthread_rwlock_destroy(&m_hRwLock);
#endif
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: unable to create mutex");
		return(eBSFerrInternal);
		}
	}
m_bMutexesCreated = true;
return(eBSFSuccess);
}

void
CAligner::DeleteMutexes(void)
{
if(!m_bMutexesCreated)
	return;
#ifdef _WIN32
CloseHandle(m_hMtxIterReads);
CloseHandle(m_hMtxMHReads);
if(m_MLMode != eMLdefault)
	CloseHandle(m_hMtxMultiMatches);
#else
pthread_mutex_destroy(&m_hMtxIterReads);
pthread_mutex_destroy(&m_hMtxMHReads);
pthread_rwlock_destroy(&m_hRwLock);
if(m_MLMode != eMLdefault)
	pthread_mutex_destroy(&m_hMtxMultiMatches);
#endif
m_bMutexesCreated = false;
}

#ifdef _WIN32
unsigned __stdcall CoredApproxThread(void * pThreadPars)
#else
void *CoredApproxThread(void * pThreadPars)
#endif
{
int Rslt;
tsThreadMatchPars *pPars = (tsThreadMatchPars *)pThreadPars; // makes it easier not having to deal with casts!
CAligner *pAligner = (CAligner *)pPars->pThis;
Rslt = pAligner->ProcCoredApprox(pPars);
pPars->Rslt = Rslt;
#ifdef _WIN32
_endthreadex(Rslt < 0 ? 1 : 0);
return(Rslt < 0 ? 1 : 0);
#else
pthread_exit(NULL);
#endif
}

// LocateCoredApprox
// Locates all cored approximates
int
CAligner::LocateCoredApprox(int MinEditDist,	// any matches must have at least this edit distance to the next best match
							int MaxSubs)		// maximum number of substitutions allowed per 100bp of accepted aligned read length
{
int Rslt;
int CurBlockID;							// current suffix block being processed
tBSFEntryID CurChromID;				    // current suffix array entry being processed
UINT32 TotNumReadsProc;					// total number of reads processed
UINT32 PlusHits;
UINT32 MinusHits;
UINT32 ChimericHits;

UINT32 CurReadsAligned;
UINT32 PrevReadsAligned;
UINT32 CurReadsLoaded;
UINT32 PrevReadsLoaded;
int MaxNumSlides;

int ThreadIdx;
tsThreadMatchPars WorkerThreads[cMaxWorkerThreads];

m_PerThreadAllocdIdentNodes = cMaxNumIdentNodes;
m_TotAllocdIdentNodes = m_PerThreadAllocdIdentNodes * m_NumThreads;
if((m_pAllocsIdentNodes = new tsIdentNode [m_TotAllocdIdentNodes])==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: unable to allocate memory for %d tsIdentNodes",m_TotAllocdIdentNodes);
	Reset(false);
	return(eBSFerrMem);
	}

if((m_pAllocsMultiHitLoci = new tsHitLoci [m_NumThreads * (m_MaxMLmatches + cPriorityExacts)])==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: unable to allocate memory for %d tsHitLoci",m_NumThreads * m_MaxMLmatches);
	Reset(false);
	return(eBSFerrMem);
	}

if(m_MLMode == eMLall)
	{
	if((m_pAllocsMultiHitBuff = new UINT8 [m_NumThreads * cReadHitBuffLen])==NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: unable to allocate memory for %d multihit record buffering",m_NumThreads * cReadHitBuffLen);
		Reset(false);
		return(eBSFerrMem);
		}
	}
else
	m_pAllocsMultiHitBuff = NULL;




// load single SfxBlock, expected to contain all chromosomes, and process all reads against that block
PlusHits = 0;
MinusHits = 0;
ChimericHits = 0;
TotNumReadsProc = 0;
PrevReadsAligned = 0;
PrevReadsLoaded = 0;

CurChromID = 0;
CurBlockID = 1;
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Loading genome assembly suffix array...");
if((Rslt=m_pSfxArray->SetTargBlock(CurBlockID))<0)
	{
	while(m_pSfxArray->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,m_pSfxArray->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: unable to load genome assembly suffix array");
	return(Rslt);
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Genome assembly suffix array loaded");

	// determine minimum core length from targeted sequence length
	// core length is a balance between sensitivity and throughput
	// reducing core size has a relatively minor effect on sensitivity but significantly reduces throughput
	// large genomes require larger cores, more sensitive alignments require smaller cores
m_BlockTotSeqLen = m_pSfxArray->GetTotSeqsLen();

if(m_BlockTotSeqLen < 20000000)				    // covers yeast
	m_MinCoreLen = cMinCoreLen;
else
	if(m_BlockTotSeqLen < 250000000)		    // covers arabidopsis and fly
		m_MinCoreLen = cMinCoreLen+3; 
	else
		m_MinCoreLen = cMinCoreLen+6; 			// covers the big guys...

switch(m_PMode) {
	case ePMUltraSens:				// ultra sensitive - much slower
		MaxNumSlides = 8;			// leave m_MinCoreLen at it's minimum
		break;
	case ePMMoreSens:				// more sensitive - slower
		m_MinCoreLen += 1;
		MaxNumSlides = 7;
		break;
	case ePMdefault:				// default processing mode
		m_MinCoreLen += 3;
		MaxNumSlides = 6;
		break;
	case ePMLessSens:				// less sensitive - quicker
		m_MinCoreLen += 5;
		MaxNumSlides = 5;
		break;
	default:
		MaxNumSlides = 4;
		break;
	}


gDiagnostics.DiagOut(eDLInfo,gszProcName,"Now aligning with minimum core size of %d...\n",m_MinCoreLen);
m_ThreadCoredApproxRslt = 0;
ResetThreadedIterReads();
memset(WorkerThreads,0,sizeof(WorkerThreads));
for(ThreadIdx = 0; ThreadIdx < m_NumThreads; ThreadIdx++)
	{
	WorkerThreads[ThreadIdx].ThreadIdx = ThreadIdx + 1;
	WorkerThreads[ThreadIdx].pThis = this;
	WorkerThreads[ThreadIdx].NumIdentNodes = m_PerThreadAllocdIdentNodes;
	WorkerThreads[ThreadIdx].pIdentNodes = &m_pAllocsIdentNodes[m_PerThreadAllocdIdentNodes * ThreadIdx];
	WorkerThreads[ThreadIdx].pMultiHits = &m_pAllocsMultiHitLoci[(m_MaxMLmatches + cPriorityExacts) * ThreadIdx];
	WorkerThreads[ThreadIdx].CurBlockID = CurBlockID;
	WorkerThreads[ThreadIdx].MinEditDist = MinEditDist;
	WorkerThreads[ThreadIdx].MaxSubs = MaxSubs;
	WorkerThreads[ThreadIdx].AlignStrand = m_AlignStrand;
	WorkerThreads[ThreadIdx].microInDelLen = m_microInDelLen;
	WorkerThreads[ThreadIdx].SpliceJunctLen = m_SpliceJunctLen;
	WorkerThreads[ThreadIdx].MaxNumSlides = MaxNumSlides;
	WorkerThreads[ThreadIdx].MinChimericLen = m_MinChimericLen;			
	if(m_MLMode == eMLall)
		WorkerThreads[ThreadIdx].pszOutBuff = &m_pAllocsMultiHitBuff[cReadHitBuffLen * ThreadIdx];
	else
		WorkerThreads[ThreadIdx].pszOutBuff = NULL;
	WorkerThreads[ThreadIdx].OutBuffIdx = 0;
#ifdef _WIN32
	WorkerThreads[ThreadIdx].threadHandle = (HANDLE)_beginthreadex(NULL,0x0fffff,CoredApproxThread,&WorkerThreads[ThreadIdx],0,&WorkerThreads[ThreadIdx].threadID);
#else
	WorkerThreads[ThreadIdx].threadRslt =	pthread_create (&WorkerThreads[ThreadIdx].threadID , NULL , CoredApproxThread , &WorkerThreads[ThreadIdx] );
#endif
	}

// allow threads a few seconds to startup
#ifdef _WIN32
	Sleep(5000);
#else
	sleep(5);
#endif

// let user know that Kanga is working hard...
ApproxNumReadsAligned(&PrevReadsAligned,&PrevReadsLoaded);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Progress: %u reads aligned from %d loaded",PrevReadsAligned,PrevReadsLoaded);

// wait for all threads to have completed
for(ThreadIdx = 0; ThreadIdx < m_NumThreads; ThreadIdx++)
	{
#ifdef _WIN32
	while(WAIT_TIMEOUT == WaitForSingleObject( WorkerThreads[ThreadIdx].threadHandle, 60000 * 10))
		{
		ApproxNumReadsAligned(&CurReadsAligned,&CurReadsLoaded);
		if(CurReadsAligned > PrevReadsAligned || CurReadsLoaded > PrevReadsLoaded)
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"Progress: %u reads aligned from %d loaded",CurReadsAligned,CurReadsLoaded);
		PrevReadsAligned = CurReadsAligned;
		PrevReadsLoaded = CurReadsLoaded;
		}
	CloseHandle( WorkerThreads[ThreadIdx].threadHandle);
#else
	struct timespec ts;
	int JoinRlt;
	clock_gettime(CLOCK_REALTIME, &ts);
	ts.tv_sec += 60 * 10;
	while((JoinRlt = pthread_timedjoin_np(WorkerThreads[ThreadIdx].threadID, NULL, &ts)) != 0)
		{
		ApproxNumReadsAligned(&CurReadsAligned,&CurReadsLoaded);
		if(CurReadsAligned > PrevReadsAligned || CurReadsLoaded > PrevReadsLoaded)
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"Progress: %u reads aligned from %d loaded",CurReadsAligned,CurReadsLoaded);
		PrevReadsAligned = CurReadsAligned;
		PrevReadsLoaded = CurReadsLoaded;
		ts.tv_sec += 60 * 10;
		}

#endif
	PlusHits += WorkerThreads[ThreadIdx].PlusHits;
	MinusHits += WorkerThreads[ThreadIdx].MinusHits;
	ChimericHits += WorkerThreads[ThreadIdx].ChimericHits;
	TotNumReadsProc += WorkerThreads[ThreadIdx].NumReadsProc;

	if(WorkerThreads[ThreadIdx].OutBuffIdx != 0)
		{
#ifdef _WIN32
		WaitForSingleObject(m_hMtxMultiMatches,INFINITE);
#else
		pthread_mutex_lock(&m_hMtxMultiMatches);
#endif
		if((cAllocLineBuffSize - m_szLineBuffIdx) < (int)(WorkerThreads[ThreadIdx].OutBuffIdx + ((cMaxFastQSeqLen * 2) + 1024)))
			{
			if(!m_bgzOutFile)
				CUtility::SafeWrite(m_hOutFile,m_pszLineBuff,m_szLineBuffIdx);
			else
				CUtility::SafeWrite_gz(m_gzOutFile,m_pszLineBuff,m_szLineBuffIdx);
			m_szLineBuffIdx = 0;
			}
		memmove(&m_pszLineBuff[m_szLineBuffIdx],WorkerThreads[ThreadIdx].pszOutBuff,WorkerThreads[ThreadIdx].OutBuffIdx);
		m_szLineBuffIdx += WorkerThreads[ThreadIdx].OutBuffIdx;
#ifdef _WIN32
		ReleaseMutex(m_hMtxMultiMatches);
#else
		pthread_mutex_unlock(&m_hMtxMultiMatches);
#endif
		WorkerThreads[ThreadIdx].OutBuffIdx = 0;
		}
	}

// pickup the read loader thread, if the reads processing threads all finished then the loader thread should also have finished
#ifdef _WIN32
if(m_hThreadLoadReads != NULL)
	{
	while(WAIT_TIMEOUT == WaitForSingleObject(m_hThreadLoadReads, 5000))
		{
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Progress: waiting");
		}
	CloseHandle(m_hThreadLoadReads);
	}
#else
if(m_ThreadLoadReadsID != 0)
	{
	struct timespec ts;
	int JoinRlt;
	clock_gettime(CLOCK_REALTIME, &ts);
	ts.tv_sec += 5;
	while((JoinRlt = pthread_timedjoin_np(m_ThreadLoadReadsID, NULL, &ts)) != 0)
		{
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Progress: waiting");
		ts.tv_sec += 60;
		}
	}
#endif

// Checking here that the reads were all loaded w/o any major dramas!
if(m_ThreadLoadReadsRslt < 0 || m_ThreadCoredApproxRslt < 0)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Progress: Early terminated");
	Reset(false);
	return(m_ThreadLoadReadsRslt < 0 ? m_ThreadLoadReadsRslt : m_ThreadCoredApproxRslt);
	}
ApproxNumReadsAligned(&CurReadsAligned,&CurReadsLoaded);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Alignment of %u from %u loaded completed",CurReadsAligned,CurReadsLoaded);

m_PerThreadAllocdIdentNodes = 0;
m_TotAllocdIdentNodes = 0;
if(m_pAllocsIdentNodes != NULL)
	{
	delete m_pAllocsIdentNodes;
	m_pAllocsIdentNodes = NULL;
	}
if(m_pAllocsMultiHitLoci != NULL)
	{
	delete m_pAllocsMultiHitLoci;
	m_pAllocsMultiHitLoci = NULL;
	}

if(m_pAllocsMultiHitBuff != NULL)
	{
	delete m_pAllocsMultiHitBuff;
	m_pAllocsMultiHitBuff = NULL;
	}

if((m_FMode == eFMsam || m_FMode == eFMsamAll) && m_MLMode == eMLall)
	m_pszLineBuff[m_szLineBuffIdx] = '\0';

if(m_szLineBuffIdx > 0)
	{
	if(!m_bgzOutFile)
		CUtility::SafeWrite(m_hOutFile,m_pszLineBuff,m_szLineBuffIdx);
	else
		CUtility::SafeWrite_gz(m_gzOutFile,m_pszLineBuff,m_szLineBuffIdx);
	m_szLineBuffIdx = 0;
	}

return(eBSFSuccess);
}



int
CAligner::ProcCoredApprox(tsThreadMatchPars *pPars)
{
etSeqBase Sequence[cMaxFastQSeqLen+1];	// to hold sequence (sans quality scores) for current read
etSeqBase PrevSequence[cMaxFastQSeqLen+1];	// to hold sequence (sans quality scores) for previous read

int MatchLen;							// match length
int PrevMatchLen;						// previous read's match length

char szPriorityChromName[100];
int PriorityChromID;

int SeqIdx;
int NumNs;
UINT8 *pSeqVal;
etSeqBase *pSeq;

tsReadsHitBlock ReadsHitBlock;			// block of reads for this thread to process
int ReadsHitIdx;						// index of current read in ReadsHitBlock.pReadsHits[]
tsReadHit *pReadHit;					// current read being processed
tsReadHit *pPrevReadHit;				// previous to current read

UINT32 NumNonAligned = 0;				// number of reads for which no alignment was discovered
UINT32 NumAcceptedAsAligned = 0;		// number of reads accepted as being aligned
UINT32 NumLociAligned = 0;				// number of loci aligned which have been reported
UINT32 NumNotAcceptedDelta = 0;			// number of reads aligned but not accepted because of insufficient hamming
UINT32 NumAcceptedHitInsts = 0;			// number of reads aligned which were accepted even though there were too many instances
UINT32 TotAcceptedAsUniqueAligned = 0;  // number of reads aligned which were uniquely aligned
UINT32 TotAcceptedAsMultiAligned = 0;   // number of reads accepted as aligned which aligned to multiple loci

int Rslt;
int HitRslt;
int PrevHitRslt;
int ProbeLen;
int RandIdx;
int HitIdx;
tsHitLoci *pHit;
int MaxTotMM;
int CoreLen;
int CoreDelta;
int MaxIter;

bool bForceNewAlignment;
int LowHitInstances;
int LowMMCnt;
int NxtLowMMCnt;

int	PrevLowHitInstances;
int	PrevLowMMCnt;
int	PrevNxtLowMMCnt;

int NumSloughedNs;

UINT32 ExtdProcFlags;

tsReadHit HitReads[cMaxMultiHits];
int MultiHitDist[cMaxMultiHits];		// used to record the multihit distribution

// iterate each read sequence starting from the first
// assumes that reads will have been sorted by ReadID
pPars->PlusHits = 0;
pPars->MinusHits = 0;
pPars->ChimericHits = 0;
pPars->NumReadsProc = 0;
pReadHit = NULL;
pPrevReadHit = NULL;
PrevMatchLen = 0;
NumSloughedNs = 0;
ReadsHitBlock.MaxReads = cMaxReadsPerBlock;
ReadsHitBlock.NumReads = 0;
Rslt = 0;
memset(MultiHitDist,0,sizeof(MultiHitDist));
ExtdProcFlags = 0;
pPrevReadHit = NULL;
PrevMatchLen = 0;
MaxIter = m_pSfxArray->GetMaxIter();


// m_hRwLock will be released and regained within ThreadedIterReads so need to always have acquired a read lock before calling ThreadedIterReads
AcquireLock(false);
bForceNewAlignment = true;			
while(ThreadedIterReads(&ReadsHitBlock))
	{
	Rslt = 0;
	for(ReadsHitIdx = 0; ReadsHitIdx < ReadsHitBlock.NumReads; ReadsHitIdx++)
		{
		pReadHit = ReadsHitBlock.pReadHits[ReadsHitIdx];
		pReadHit->NAR = eNARNoHit;			// assume unable to align read
		pReadHit->FlgPEAligned = 0;
		pReadHit->LowHitInstances = 0;
		pReadHit->LowMMCnt = 0;
		pReadHit->NumHits = 0;
		pPars->NumReadsProc+=1;

		// get sequence for read and remove any packed quality values
		pSeqVal = &pReadHit->Read[pReadHit->DescrLen+1];
		pSeq = Sequence;
		NumNs = 0;
		int MaxNsSeq = 0;		// maximum allowed for this sequence
		if(m_MaxNs)
			MaxNsSeq = max(((pReadHit->ReadLen * m_MaxNs) / 100),m_MaxNs);

		for(SeqIdx = 0; SeqIdx < pReadHit->ReadLen; SeqIdx++,pSeq++,pSeqVal++)
			{
			if((*pSeq = (*pSeqVal & 0x07)) > eBaseN)
				break;
			*pSeq &= ~cRptMskFlg;
			if(*pSeq == eBaseN)
				{
				if(++NumNs > MaxNsSeq)
					break;
				}
			}
		if(SeqIdx != pReadHit->ReadLen) // if too many 'N's...
			{
			pPrevReadHit = NULL;
			pReadHit->NAR = eNARNs;
			PrevMatchLen = 0;
			NumSloughedNs += 1;
			continue;
			}

		if(m_bIsSOLiD)	// if SOLiD colorspace then need to convert read back into colorspace before attempting to locate
			{
			pSeq = Sequence;
			UINT8 PrvBase = *pSeq;
			for(SeqIdx = 1; SeqIdx < pReadHit->ReadLen; SeqIdx++,pSeq++)
				{
				*pSeq = SOLiDmap[PrvBase][pSeq[1]];
				PrvBase = pSeq[1];
				}
			ProbeLen = pReadHit->ReadLen;
			MatchLen = ProbeLen-1;
			}
		else
			{
			ProbeLen = pReadHit->ReadLen;
			MatchLen = ProbeLen;
			}

		MaxTotMM = pPars->MaxSubs == 0 ? 0 : max(1,(MatchLen * pPars->MaxSubs)/100); // any accepted alignment can have at most this many mismatches
		if(MaxTotMM > cMaxTotAllowedSubs)		// irrespective of length allow at most this many subs
			MaxTotMM = cMaxTotAllowedSubs;

		// The window core length is set to be read length / (subs+1) for minimum Hamming difference of 1, and
		// to be read length / (subs+2) for minimum Hamming difference of 2
		// The window core length is clamped to be at least m_MinCoreLen
		CoreLen = max(m_MinCoreLen,MatchLen/(pPars->MinEditDist == 1 ? MaxTotMM+1 : MaxTotMM+2));
		CoreDelta = CoreLen;

		LowHitInstances = pReadHit->LowHitInstances;
		LowMMCnt = pReadHit->LowMMCnt;
		NxtLowMMCnt = pReadHit->NxtLowMMCnt;

		ExtdProcFlags &= ~0x0001;			// no trace

		// for priority alignments to known reference sequences (example would be cDNA transcripts assembled with a de Novo assembly) then
		// a) align for exact matches allowing say 10 multiloci hits
		// b) iterate these multiloci hits and discard any not aligning to a ref sequence
		// c) process those aligning to reference as if these were uniquely aligning
		int RefExacts = cPriorityExacts;
		bool bProcNorm;

		if(m_pPriorityRegionBED != NULL)
			RefExacts = cPriorityExacts;
		else
			RefExacts = 0;
			
		// Heuristic is that if read sequence is identical to previously processed sequence then
		// simply reuse previously AlignReads() hit results - saves a lot of processing time
		if(pPrevReadHit == NULL || bForceNewAlignment || MatchLen != PrevMatchLen || memcmp(Sequence,PrevSequence,MatchLen))
			{
			memset(pPars->pMultiHits,0,sizeof(tsHitLoci));
			bProcNorm = true;
			if(RefExacts > 0)
				{
				HitRslt = m_pSfxArray->AlignReads(ExtdProcFlags,						// flags indicating if lower levels need to do any form of extended processing with this specific read...
														pReadHit->ReadID,				// identifies this read
														pPars->MinChimericLen,			// minimum chimeric length as a percentage (0 to disable, otherwise 50..99) of probe sequence
														MaxTotMM,						// max number of mismatches allowed
														CoreLen,						// core window length
														CoreDelta,						// core window offset increment (1..n)
														pPars->MaxNumSlides,			// max number of times to slide core on each strand
														pPars->MinEditDist,				// minimum (1..n) mismatch difference between the best and next best core alignment
														pPars->AlignStrand,				// watson, crick or both?
														0,								// microInDel length maximum
														0,								// maximum splice junction length when aligning RNAseq reads
														&LowHitInstances,				// Out number of match instances for lowest number of mismatches thus far for this read
														&LowMMCnt,						// Out lowest number of mismatches thus far for this read
														&NxtLowMMCnt,					// Out next to lowest number of mismatches thus far for this read
														Sequence,						// probe
														MatchLen,						// probe length
														m_MaxMLmatches+RefExacts,		// (IN) process for at most this number of hits
														pPars->pMultiHits,				// where to return the loci for each hit by current read
														pPars->NumIdentNodes,			// memory has been allocated by caller for holding up to this many tsIdentNodes
														pPars->pIdentNodes);			// memory allocated by caller for holding tsIdentNodes

				if(HitRslt == eHRhits)
					{
					int NumInPriorityRegions;
					int NumInNonPriorityRegions;
					tsHitLoci *pPriorityHits;
					NumInPriorityRegions = 0;
					NumInNonPriorityRegions = 0;
					pHit = pPars->pMultiHits;
					pPriorityHits = pHit;
					for(HitIdx=0; HitIdx < min(LowHitInstances,m_MaxMLmatches+RefExacts); HitIdx++,pHit++)
						{
						// check if hit loci within region designated as being a priority exact matching region
						if(m_pSfxArray->GetIdentName(pHit->Seg->ChromID,sizeof(szPriorityChromName),szPriorityChromName)!=eBSFSuccess)
							{
							NumInNonPriorityRegions += 1;
							continue;
							}
						if((PriorityChromID = m_pPriorityRegionBED->LocateChromIDbyName(szPriorityChromName)) < 1)
							{
							NumInNonPriorityRegions += 1;
							continue;
							}
						if(!m_pPriorityRegionBED->InAnyFeature(PriorityChromID,(int)pHit->Seg->MatchLoci,(int)(pHit->Seg->MatchLoci+pHit->Seg->MatchLen-1)))
							{
							NumInNonPriorityRegions += 1;
							continue;
							}
						if(NumInNonPriorityRegions > 0)
							*pPriorityHits++ = *pHit;
						NumInPriorityRegions += 1;
						}
					if(NumInPriorityRegions > 0 && (m_bClampMaxMLmatches || NumInPriorityRegions <= m_MaxMLmatches))
						{
						if(NumInPriorityRegions > m_MaxMLmatches)
							NumInPriorityRegions = m_MaxMLmatches;
						LowHitInstances = NumInPriorityRegions;
						bProcNorm = false;
						}
					}
				if(bProcNorm)
					{
					LowHitInstances = pReadHit->LowHitInstances;
					LowMMCnt = pReadHit->LowMMCnt;
					NxtLowMMCnt = pReadHit->NxtLowMMCnt;
					}
				}
			else
				bProcNorm = true;

			if(bProcNorm)
				{
				if(m_bLocateBestMatches)
					{
					HitRslt =						// < 0 if errors, 0 if no matches, 1..MaxHits, or MaxHits+1 if additional matches have been sloughed
						m_pSfxArray->LocateBestMatches(pReadHit->ReadID,			// identifies this read
											MaxTotMM,			        // return matches having at most this number of mismatches
											CoreLen,					// core window length
											CoreDelta,					// core window offset increment (1..n)
											pPars->MaxNumSlides,		// max number of times to slide core on each strand
											pPars->AlignStrand,			// watson, crick or both?
 											Sequence,MatchLen,
											m_MaxMLmatches,				// process for at most this many hits by current read
											&LowHitInstances,			// returned number of match instances in pHits
											pPars->pMultiHits,			// where to return the loci for each hit by current read
											 MaxIter,					// max allowed iterations per subsegmented sequence when matching that subsegment
											 pPars->NumIdentNodes,		// memory has been allocated by caller for holding up to this many tsIdentNodes
											 pPars->pIdentNodes);		// memory allocated by caller for holding tsIdentNodes
					if(HitRslt == 0)
						HitRslt = eHRnone;
					else
						if(HitRslt >= 1)
							HitRslt = eHRhits;
					}
				else
					HitRslt = m_pSfxArray->AlignReads(ExtdProcFlags,					// flags indicating if lower levels need to do any form of extended processing with this specific read...
														pReadHit->ReadID,				// identifies this read
														pPars->MinChimericLen,			// minimum chimeric length as a percentage (0 to disable, otherwise 50..99) of probe sequence
														MaxTotMM,CoreLen,CoreDelta,
														pPars->MaxNumSlides,
														pPars->MinEditDist,
														pPars->AlignStrand,				// watson, crick or both?
														pPars->microInDelLen,			// microInDel length maximum
														pPars->SpliceJunctLen,			// maximum splice junction length when aligning RNAseq reads
														&LowHitInstances,
														&LowMMCnt,
														&NxtLowMMCnt,
														Sequence,MatchLen,
														m_MaxMLmatches,					// process for at most this many hits by current read
														pPars->pMultiHits,				// where to return the loci for each hit by current read
														pPars->NumIdentNodes,
														pPars->pIdentNodes);
				}

			// user may be interested in multihits upto m_MaxMLmatches limit even if there were actually many more so in this case remap the hit result
			if(LowHitInstances > m_MaxMLmatches)
				LowHitInstances = m_MaxMLmatches + 1;
			if(m_bClampMaxMLmatches && HitRslt == eHRHitInsts)
				{
				LowHitInstances=m_MaxMLmatches;
				HitRslt = eHRhits;
				bForceNewAlignment = true;
				}
			PrevHitRslt = HitRslt;
			PrevLowHitInstances = LowHitInstances;
			PrevLowMMCnt = LowMMCnt;
			PrevNxtLowMMCnt = NxtLowMMCnt;
			PrevMatchLen = MatchLen;
			pPrevReadHit = pReadHit;
			memmove(PrevSequence,Sequence,MatchLen);
			}
		else		// else, identical sequence, reuse previous hit results from AlignReads
			{
			HitRslt = PrevHitRslt;
			LowHitInstances = PrevLowHitInstances;
			LowMMCnt = PrevLowMMCnt;
			NxtLowMMCnt = PrevNxtLowMMCnt;
			}

		// if SOLiD colorspace then need to normalise the hits back to as if aligned in basespace
		if(m_bIsSOLiD)
			{
			switch(HitRslt) {
				case eHRHitInsts:
					if(m_MLMode != eMLall)
						break;
				case eHRhits:
					pHit = pPars->pMultiHits;
					for(HitIdx=0; HitIdx <  min(LowHitInstances,m_MaxMLmatches); HitIdx++,pHit++)
						{
						pHit->Seg[0].MatchLoci -= 1;
						pHit->Seg[0].MatchLen += 1;
						if(pHit->FlgInDel == 1 || pHit->FlgSplice == 1)
							pHit->Seg[1].ReadOfs += 1;
						}
					break;
				default:
					break;
				}
			}

		// ensure that TrimLeft/Right/Mismatches are consistent with the fact that trimming is a post alignment phase
		// if alignment was the result of a chimeric alignment then accept the flank trimming
		if(HitRslt == eHRHitInsts ||  HitRslt == eHRhits)
			{
			pHit = pPars->pMultiHits;
			for(HitIdx=0; HitIdx < min(LowHitInstances,m_MaxMLmatches); HitIdx++,pHit++)
				{
				if(pHit->FlgChimeric != 1)
					{
					pHit->Seg[0].TrimLeft = 0;
					pHit->Seg[0].TrimRight = 0;
					}
				pHit->Seg[0].TrimMismatches = pHit->Seg[0].Mismatches;
				if(pHit->FlgInDel == 1 || pHit->FlgSplice == 1)
					{
					pHit->Seg[1].TrimLeft = 0;
					pHit->Seg[1].TrimRight = 0;
					pHit->Seg[1].TrimMismatches = pHit->Seg[1].Mismatches;
					}
				}
			}


		Rslt = 0;
		switch(HitRslt) {
			case eHRnone:			// no change or no hits
				pReadHit->NAR = eNARNoHit;
				pReadHit->LowMMCnt = 0;
				pReadHit->NumHits = 0;
				pReadHit->LowHitInstances = 0;
				if(m_FMode == eFMsamAll)
					{
					pReadHit->LowMMCnt = (INT8)0;
					if((Rslt = WriteHitLoci(pPars,pReadHit,0,pPars->pMultiHits)) < 0)
						break;
					}
				NumNonAligned += 1;
				bForceNewAlignment = false;
				break;

			case eHRhits:			// MMDelta criteria met and within the max allowed number of hits
				if(LowHitInstances > 1)
					bForceNewAlignment = true;
				else
					bForceNewAlignment = false;

				pReadHit->NAR = eNARAccepted;
				if(m_MLMode == eMLall)
					{
					int NumWriteHitLoci;
					// report each hit here...
					pReadHit->LowMMCnt = (INT8)LowMMCnt;
					if((Rslt = NumWriteHitLoci = WriteHitLoci(pPars,pReadHit,LowHitInstances,pPars->pMultiHits)) < 0)
						break;
					if(NumWriteHitLoci)
						{
						NumAcceptedAsAligned += 1;
						NumLociAligned += NumWriteHitLoci;
						if(NumWriteHitLoci == 1)
							TotAcceptedAsUniqueAligned += 1;
						else
							TotAcceptedAsMultiAligned += 1;
						}
					HitRslt = eHRHitInsts;
					break;
					}

				NumAcceptedAsAligned += 1;
				NumLociAligned += LowHitInstances;

				if(LowHitInstances == 1)
					TotAcceptedAsUniqueAligned += 1;
				else
					TotAcceptedAsMultiAligned += 1;


				MultiHitDist[LowHitInstances-1] += 1;
				if(m_MLMode == eMLrand)
					RandIdx = m_MaxMLmatches == 1 ? 0 : (rand() % LowHitInstances);
				else
					RandIdx = 0;
				if(m_MLMode <= eMLrand || LowHitInstances == 1)
					{
					if((m_MLMode == eMLdist && LowHitInstances == 1) || m_MLMode != eMLdist)
						{
						if(LowHitInstances == 1)	// was this a unique hit?
							pReadHit->HitLoci.FlagHL = (int)eHLunique;
						else
							pReadHit->HitLoci.FlagHL = (int)eHLrandom;
						LowHitInstances = 1;			// currently just accepting one
						pReadHit->NumHits = 1;
						pReadHit->HitLoci.Hit = pPars->pMultiHits[RandIdx];
						pReadHit->HitLoci.FlagSegs = (pReadHit->HitLoci.Hit.FlgInDel == 1 || pReadHit->HitLoci.Hit.FlgSplice == 1) ? 1 : 0;
						}
					else
						{
						pReadHit->NAR = eNARMultiAlign;
						pReadHit->NumHits = 0;
						}
					}
				else
					if(LowHitInstances == 1)
						{
						pReadHit->NAR = eNARAccepted;
						pReadHit->NumHits = 1;
						}
					else
						{
						pReadHit->NAR = eNARMultiAlign;
						pReadHit->NumHits = 0;
						}

				pReadHit->LowHitInstances = (INT16)LowHitInstances;
				pReadHit->LowMMCnt = (INT8)LowMMCnt;
				pReadHit->NxtLowMMCnt = (INT8)NxtLowMMCnt;
// handling multiply aligned reads
				if(m_MLMode > eMLrand)
					{
					int HitIdx;
					tsHitLoci *pHit;
					tsReadHit *pMHit;
					pHit = pPars->pMultiHits;
					pMHit = &HitReads[0];
					for(HitIdx=0; HitIdx < LowHitInstances; HitIdx++,pMHit++,pHit++)
						{
						memcpy(pMHit,pReadHit,sizeof(tsReadHit));
						pMHit->FlgPEAligned = false;
						pMHit->HitLoci.Hit = *pHit;
						pMHit->HitLoci.FlagSegs = (pMHit->HitLoci.Hit.FlgInDel == 1 || pMHit->HitLoci.Hit.FlgSplice == 1) ? 1 : 0;
						pMHit->HitLoci.FlagMH = LowHitInstances > 1 ? 1 : 0;
						pMHit->HitLoci.FlagMHA = 0;
						}
					if((Rslt=AddMHitReads(LowHitInstances,&HitReads[0])) < 0)		// pts to array of hit loci
						break;
					}
// finish handling multiply aligned reads
				break;

			case eHRMMDelta:			// same or a new LowMMCnt unique hit but MMDelta criteria not met
				NumNotAcceptedDelta += 1;
				memset(&pReadHit->HitLoci.Hit.Seg[0],0,sizeof(tsSegLoci));
				pReadHit->NAR = eNARMMDelta;
				pReadHit->NumHits = 0;
				pReadHit->HitLoci.Hit.Seg[0].Strand = '?';
				pReadHit->HitLoci.Hit.BisBase = eBaseN;
				pReadHit->HitLoci.Hit.Seg[0].MatchLen = (UINT16)ProbeLen;
				pReadHit->LowHitInstances = (INT16)LowHitInstances;
				pReadHit->LowMMCnt = (INT8)LowMMCnt;
				pReadHit->NxtLowMMCnt = (INT8)NxtLowMMCnt;
				bForceNewAlignment = false;
				break;

			case eHRHitInsts:			// same or new LowMMCnt but now simply too many multiple hit instances, treat as none-aligned
				pReadHit->NAR = eNARMultiAlign;
				if(m_MLMode == eMLall && m_FMode == eFMsamAll)
					{
					pReadHit->LowMMCnt = (INT8)0;
					if((Rslt = WriteHitLoci(pPars,pReadHit,0,pPars->pMultiHits)) < 0)
						break;
					}
				NumNonAligned += 1;
				bForceNewAlignment = false;
				break;

			case eHRRMMDelta:			// reduced NxtLowMMCnt only
				pReadHit->NxtLowMMCnt = (INT8)NxtLowMMCnt;
				bForceNewAlignment = false;
				break;
			}

		if(Rslt < 0)
			{
			ReleaseLock(false);
			AcquireSerialise();
			m_ThreadCoredApproxRslt = Rslt;
			ReleaseSerialise(); 
			return(-1);
			}
			

		if(HitRslt == eHRHitInsts)
			{
			pReadHit->NumHits = 0;
			pReadHit->NAR = eNARMultiAlign;
			memset(&pReadHit->HitLoci.Hit.Seg[0],0,sizeof(tsSegLoci));
			pReadHit->HitLoci.Hit.Seg[0].Strand = '?';
			pReadHit->HitLoci.Hit.BisBase = eBaseN;
			pReadHit->HitLoci.Hit.Seg[0].MatchLen = (UINT16)ProbeLen;
			pReadHit->LowHitInstances = (INT16)LowHitInstances;
			pReadHit->LowMMCnt = (INT8)LowMMCnt;
			pReadHit->NxtLowMMCnt = (INT8)NxtLowMMCnt;
			}

		// if SOLiD colorspace then need to restore any hits to what they were immediately following return from AlignReads
		// as when SOLiD processing MatchLoci/MatchLen and ReadOfs would have all been normalised as if basespace processing and if
		// a duplicate of this read sequence follows then hits will be simply reused w/o the overhead of a call to AlignReads
		if(m_bIsSOLiD)
			{
			switch(HitRslt) {
				case eHRHitInsts:
					if(m_MLMode != eMLall)
						break;
				case eHRhits:
					pHit = pPars->pMultiHits;
					for(HitIdx=0; HitIdx <  min(LowHitInstances,m_MaxMLmatches); HitIdx++,pHit++)
						{
						pHit->Seg[0].MatchLoci += 1;
						pHit->Seg[0].MatchLen -= 1;
						if(pHit->FlgInDel == 1 || pHit->FlgSplice == 1)
							pHit->Seg[1].ReadOfs -= 1;
						}
					break;
				default:
					break;
				}
			}
		}
	}

ReleaseLock(false);
AcquireSerialise();
m_NumSloughedNs += NumSloughedNs;
m_TotNonAligned += NumNonAligned;
m_TotAcceptedAsUniqueAligned += TotAcceptedAsUniqueAligned;
m_TotAcceptedAsMultiAligned += TotAcceptedAsMultiAligned;
m_TotAcceptedAsAligned += NumAcceptedAsAligned;
m_TotLociAligned += NumLociAligned;
m_TotNotAcceptedDelta += NumNotAcceptedDelta;
m_TotAcceptedHitInsts += NumAcceptedHitInsts;

if(m_MLMode != eMLall)
	{
	int Idx;
	for(Idx = 0; Idx < m_MaxMLmatches; Idx++)
		m_MultiHitDist[Idx] += MultiHitDist[Idx];
	}

ReleaseSerialise();
return(1);
}

// LocateRead
tsReadHit *
CAligner::LocateRead(UINT32 ReadID) // Locate read with requested ReadID
{
int Rslt;
tsReadHit *pProbe;
int Lo,Mid,Hi;	// search limits

if(m_ppReadHitsIdx == NULL || m_AllocdReadHitsIdx < m_NumReadsLoaded)
	return(NULL);

Lo = 0;
Hi = m_NumReadsLoaded-1;
while(Hi >= Lo) {
	Mid = (Hi + Lo)/2;
	pProbe = m_ppReadHitsIdx[Mid];
	Rslt = pProbe->ReadID - ReadID;
	if(Rslt > 0)
		{
		Hi = Mid - 1;
		continue;
		}
	if(Rslt < 0)
		{
		Lo = Mid + 1;
		continue;
		}
	return(pProbe);
	}
return(NULL);
}


int
CAligner::AddMHitReads(UINT32 NumHits,	// number of multimatches loci in pHits
		tsReadHit *pHits)		// pts to array of hit loci
{
size_t memreq;
tsReadHit *pDstHits;
// ensure actually processing multihits
if(m_MLMode <= eMLrand)
	return(0);					// silently slough these hits

AcquireSerialiseMH();

if((m_AllocdMultiHits - m_NumMultiHits) < (NumHits+1000))	// need to realloc? -- added 1000 to provide a little safety margin
	{
	memreq = (m_AllocdMultiHits + cAllocMultihits) * sizeof(tsReadHit);
#ifdef _WIN32
	pDstHits = (tsReadHit *) realloc(m_pMultiHits,memreq);
#else
	pDstHits = (tsReadHit *)mremap(m_pMultiHits,m_AllocdMultiHitsMem,memreq,MREMAP_MAYMOVE);
	if(pDstHits == MAP_FAILED)
		pDstHits = NULL;
#endif
	if(pDstHits == NULL)
		{
		ReleaseSerialiseMH();
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddMHitReads: Memory re-allocation to %lld bytes - %s",memreq,strerror(errno));
		return(eBSFerrMem);
		}
	m_pMultiHits = pDstHits;
	m_AllocdMultiHitsMem = memreq;
	m_AllocdMultiHits += cAllocMultihits;
	}
pDstHits = &m_pMultiHits[m_NumMultiHits];
memmove(pDstHits,pHits,sizeof(tsReadHit) * NumHits);
m_NumMultiHits += NumHits;
if(NumHits == 1)
	m_NumUniqueMultiHits += 1;
else
	m_NumProvMultiAligned += 1;
ReleaseSerialiseMH();
return((int)NumHits);
}


// ThreadedIterReads
// use to return a block of reads reserved for processing by calling thread instance
// Returns false if no more reads available for processing
void
CAligner::ResetThreadedIterReads(void) // must be called by master thread prior to worker threads calling ThreadedIterReads()
{
m_NumReadsProc = 0;
m_NxtReadProcOfs = 0;
m_ProcessingStartSecs = gStopWatch.ReadUSecs();
}


UINT32		// Returns the number of reads thus far loaded and aligned
CAligner::ApproxNumReadsAligned(UINT32 *pNumAligned,UINT32 *pNumLoaded)
{
UINT32 NumAligned;
AcquireSerialise();
NumAligned = m_NumReadsProc;
if(pNumAligned != NULL)
	*pNumAligned = NumAligned;
*pNumLoaded = m_NumReadsLoaded;
ReleaseSerialise();
return(NumAligned);
}

// ThreadedIterReads
// Iterates over all reads
// The body of this function is serialised
bool	// returns false if no more reads availing for processing by calling thread
CAligner::ThreadedIterReads(tsReadsHitBlock *pRetBlock)
{
UINT32 NumReadsLeft;
UINT32 MaxReads2Proc;
tsReadHit *pCurReadHit;
pRetBlock->NumReads = 0;


ReleaseLock(false);
while(1) {
	AcquireSerialise();
	AcquireLock(false);
	if(m_bAllReadsLoaded || m_NumReadsProc != m_NumReadsLoaded || m_ThreadCoredApproxRslt < 0)
    	break;

	ReleaseLock(false);
	ReleaseSerialise();
#ifdef _WIN32
	Sleep(2000);			// must have caught up to the reads loader, allow it some breathing space to parse and load some more reads...
#else
	sleep(2);
#endif
	}

if(m_pReadHits == NULL ||
	m_ThreadCoredApproxRslt < 0 ||
	(m_bAllReadsLoaded && (m_LoadReadsRslt != eBSFSuccess)) ||
	m_bAllReadsLoaded && (m_NumReadsLoaded == 0 || m_NumReadsProc == m_NumReadsLoaded)) // if all reads have been loaded and all processed then time to move onto next processing phase
	{
	pRetBlock->NumReads = 0;
	pRetBlock->pReadHits[0] = NULL;
	ReleaseSerialise();
	return(false);
	}

// adjust pRetBlock->MaxReads according to the number of reads remaining and threads still processing these reads
// idea is to maximise the number of threads still processing when most reads have been processed so that
// the last thread processing doesn't end up with a large block of reads needing lengthly processing
NumReadsLeft = m_NumReadsLoaded - m_NumReadsProc;
if(NumReadsLeft < cMaxReadsPerBlock/4)	// if < cMaxReadsPerBlock/4 yet to be processed then give it all to the one thread
	MaxReads2Proc = NumReadsLeft;
else
	MaxReads2Proc = min((UINT32)pRetBlock->MaxReads,10 + (NumReadsLeft / (UINT32)m_NumThreads));
MaxReads2Proc = min(MaxReads2Proc,NumReadsLeft);
if(!m_NumReadsProc)
	m_NxtReadProcOfs = 0;
pCurReadHit = (tsReadHit *)((UINT8 *)m_pReadHits + m_NxtReadProcOfs);

while(MaxReads2Proc)
	{
	pRetBlock->pReadHits[pRetBlock->NumReads++] = pCurReadHit;
	pCurReadHit = (tsReadHit *)((UINT8 *)pCurReadHit + sizeof(tsReadHit) + pCurReadHit->ReadLen + pCurReadHit->DescrLen);
	MaxReads2Proc -= 1;
	}

m_NumReadsProc += pRetBlock->NumReads;
m_NxtReadProcOfs = (size_t)((UINT8 *)pCurReadHit - (UINT8 *)m_pReadHits);

ReleaseSerialise();
return(true);
}


// IterReads
// use to iterate over reads returning ptr to next read following the current read
// Reads need not be sorted
// To start from first read then pass in NULL as pCurReadHit
// Returns NULL if all read hits have been iterated
tsReadHit *
CAligner::IterReads(tsReadHit *pCurReadHit) // to start from first read then pass in NULL as pCurReadHit
{
tsReadHit *pNxtReadHit = NULL;
if(pCurReadHit == NULL)
	pNxtReadHit = m_pReadHits;
else
	if(pCurReadHit->ReadID != m_FinalReadID)
		pNxtReadHit = (tsReadHit *)((UINT8 *)pCurReadHit + sizeof(tsReadHit) + pCurReadHit->ReadLen + pCurReadHit->DescrLen);
return(pNxtReadHit);
}

// IterSortedReads
// use to iterate over sorted reads returning ptr to next read following the current read
// To start from first read then pass in NULL as pCurReadHit
// Returns NULL if all read hits have been iterated
tsReadHit *
CAligner::IterSortedReads(tsReadHit *pCurReadHit)
{
tsReadHit *pNxtReadHit = NULL;
if(pCurReadHit == NULL)
	pNxtReadHit = m_ppReadHitsIdx[0];
else
	if(pCurReadHit->ReadHitIdx < m_NumReadsLoaded)
		pNxtReadHit = m_ppReadHitsIdx[pCurReadHit->ReadHitIdx];
return(pNxtReadHit);
}

tsReadHit *		// returned read which overlaps StartLoci and EndLoci, NULL if no read located
CAligner::IterateReadsOverlapping(bool bTriSNPs, // false if iterating DiSNPs, true if iterating TriSNPs
						tsChromSNPs *pChromSNPs, // processing SNPs on this chromosome
						int StartLoci,				// returned reads are required to overlap both this starting and
						int EndLoci)				// this ending loci
{
bool bNewStartEnd;
int CurStart;
int CurEnd;
tsAdjacentSNPs *pAdjacentSNPs;
tsReadHit *pNxtReadHit = NULL;

pAdjacentSNPs = bTriSNPs ? &pChromSNPs->AdjacentSNPs[1] : &pChromSNPs->AdjacentSNPs[0];
bNewStartEnd = false;
if(pAdjacentSNPs->pFirstIterReadHit == NULL ||				// only NULL if iterations just starting for SNPs on this chrom
	StartLoci <  pAdjacentSNPs->StartLoci || EndLoci <  pAdjacentSNPs->EndLoci)	// or if StartLoci/EndLoci now 5' to previous so needing to search from start
	{
	pNxtReadHit = pChromSNPs->pFirstReadHit;
	pAdjacentSNPs->pFirstIterReadHit = pChromSNPs->pFirstReadHit;
	pAdjacentSNPs->pPrevIterReadHit = NULL;
	pAdjacentSNPs->StartLoci = StartLoci;
	pAdjacentSNPs->EndLoci = EndLoci;
	bNewStartEnd = true;
	}
else
	{
	if(StartLoci == pAdjacentSNPs->StartLoci && EndLoci == pAdjacentSNPs->EndLoci) // same start/end loci as previous so must be iterating same start/end loci
		{
		if(pAdjacentSNPs->pPrevIterReadHit == pChromSNPs->pLastReadHit)		// check if last returned was last on chrom
			return(NULL);

		pNxtReadHit = m_ppReadHitsIdx[pAdjacentSNPs->pPrevIterReadHit->ReadHitIdx]; // recommence search from next
		}
	else   // else starting new Start/EndLoci which is 3' to previous
		{
		pNxtReadHit = pAdjacentSNPs->pFirstIterReadHit;
		pAdjacentSNPs->pPrevIterReadHit = pNxtReadHit;
		pAdjacentSNPs->StartLoci = StartLoci;
		pAdjacentSNPs->EndLoci = EndLoci;
		bNewStartEnd = true;
		}	
	}

	// have the current read, iterate forward returning reads until reads align past the StartLoci
do {
	if(pNxtReadHit == NULL)
		return(NULL);
	if(pNxtReadHit->NAR == eNARAccepted &&
		!(pNxtReadHit->HitLoci.Hit.FlgInDel || pNxtReadHit->HitLoci.Hit.FlgSplice))
		{
		if(pNxtReadHit->HitLoci.Hit.Seg[0].ChromID !=  pChromSNPs->ChromID) // double check read is aligning on to expected chrom
			return(NULL);

		CurStart = AdjStartLoci(&pNxtReadHit->HitLoci.Hit.Seg[0]); 
		CurEnd = AdjEndLoci(&pNxtReadHit->HitLoci.Hit.Seg[0]);

		if(CurStart <= StartLoci && CurEnd >= EndLoci)   // this read overlaps both start and end?
			{
			pAdjacentSNPs->pPrevIterReadHit = pNxtReadHit;
			if(bNewStartEnd)		
				pAdjacentSNPs->pFirstIterReadHit = pNxtReadHit;
			return(pNxtReadHit);
			}
		else
			if(CurStart > StartLoci)
				return(NULL);
		}

	}
while((pNxtReadHit != pChromSNPs->pLastReadHit) && (pNxtReadHit = m_ppReadHitsIdx[pNxtReadHit->ReadHitIdx])!=NULL);

return(NULL);
}

// NumDnUniques
// Returns number of unique loci to which reads map which are downstream of the specified current read
int											// returned number of unique reads
CAligner::NumDnUniques(tsReadHit *pCurReadHit,		// current read
				int WinLen,					// only interested in unique reads starting within this window
				bool bStrandDep)			// if true then unique loci reads must be on current read stand
{
int NumUniques;
UINT32 CurChromID;
UINT8 CurStrand;
int CurStart;
int NxtStart;
int PrvNxtStart;
UINT32 NxtReadIdx;
tsReadHit *pNxtReadHit = NULL;

if(pCurReadHit == NULL)									// if NULL then start from first
	NxtReadIdx = 0;
else
	if((NxtReadIdx = pCurReadHit->ReadHitIdx) == m_NumReadsLoaded)
		return(0);

// have the next read, iterate forward counting the uniques until loci outside of window
CurChromID = pCurReadHit->HitLoci.Hit.Seg[0].ChromID;
CurStart = AdjStartLoci(&pCurReadHit->HitLoci.Hit.Seg[0]);
PrvNxtStart = CurStart;
CurStrand = pCurReadHit->HitLoci.Hit.Seg[0].Strand;
NumUniques = 0;
do {
	pNxtReadHit = m_ppReadHitsIdx[NxtReadIdx];
	NxtReadIdx = pNxtReadHit->ReadHitIdx;
	if(CurChromID != pNxtReadHit->HitLoci.Hit.Seg[0].ChromID)	// must be on same chromosome
		return(NumUniques);
	if(pNxtReadHit->NAR != eNARAccepted)						// skip over any read not uniquely aligned
		continue;
	NxtStart = AdjStartLoci(&pNxtReadHit->HitLoci.Hit.Seg[0]);

	if((CurStart + WinLen) < NxtStart)					// outside of window?
		return(NumUniques);

	if(bStrandDep && CurStrand != pNxtReadHit->HitLoci.Hit.Seg[0].Strand) // has to be on same strand?
		continue;

	if(NxtStart != PrvNxtStart)
		{
		NumUniques += 1;
		PrvNxtStart = NxtStart;
		}
	}
while(NxtReadIdx != m_NumReadsLoaded);
return(NumUniques);
}


// NumUpUniques
// Returns number of unique loci to which reads map which are upstream of the specified current read
int											// returned number of unique reads
CAligner::NumUpUniques(tsReadHit *pCurReadHit,		// current read
				int WinLen,					// only interested in unique reads starting within this window
				bool bStrandDep)			// if true then unique loci reads must be on current read stand
{
int NumUniques;
UINT32 CurChromID;
UINT8 CurStrand;
int CurStart;
int NxtStart;
int PrvNxtStart;
UINT32 NxtReadIdx;
tsReadHit *pNxtReadHit = NULL;

if(pCurReadHit == NULL || (NxtReadIdx = pCurReadHit->ReadHitIdx)==1) // if NULL or first then
	return(0);								// can't be any upstream from the first!

CurChromID = pCurReadHit->HitLoci.Hit.Seg[0].ChromID;
CurStart = AdjStartLoci(&pCurReadHit->HitLoci.Hit.Seg[0]);
PrvNxtStart = CurStart;
CurStrand = pCurReadHit->HitLoci.Hit.Seg[0].Strand;
NumUniques = 0;
do {
	pNxtReadHit = m_ppReadHitsIdx[NxtReadIdx-1];
	if(CurChromID != pNxtReadHit->HitLoci.Hit.Seg[0].ChromID)		// must be on same chromosome
		return(NumUniques);
	if(pNxtReadHit->NAR != eNARAccepted)						// skip over any read not uniquely aligned
		continue;
	NxtStart = AdjStartLoci(&pNxtReadHit->HitLoci.Hit.Seg[0]);

	if(CurStart > WinLen && (CurStart - WinLen) > NxtStart)					// outside of window?
		return(NumUniques);

	if(bStrandDep && CurStrand != pNxtReadHit->HitLoci.Hit.Seg[0].Strand) // has to be on same strand?
		continue;

	if(NxtStart != PrvNxtStart)
		{
		NumUniques += 1;
		PrvNxtStart = NxtStart;
		}
	}
while(--NxtReadIdx > 0);
return(NumUniques);
}

int
CAligner::SortReadHits(etReadsSortMode SortMode,		// sort mode required
				bool bSeqSorted,			// used to optimise eRSMSeq processing, if it is known that reads are already sorted in sequence order (loaded from pre-processed .rds file)
				bool bForce)				// if true then force sort
{
tsReadHit *pReadHit;
UINT32 Idx;

if(!bForce && SortMode == m_CurReadsSortMode && m_ppReadHitsIdx != NULL && m_AllocdReadHitsIdx >= m_NumReadsLoaded)
	return(eBSFSuccess);					// if already in requested mode

if(m_ppReadHitsIdx == NULL || m_AllocdReadHitsIdx < m_NumReadsLoaded)
	{
	if(m_ppReadHitsIdx != NULL)
		{
		delete m_ppReadHitsIdx;
		m_ppReadHitsIdx = NULL;
		m_AllocdReadHitsIdx = 0;
		}
	if((m_ppReadHitsIdx = (tsReadHit **) new tsReadHit * [m_NumReadsLoaded])==NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"SortReadHits: Memory reads index allocation for %d ptrs - %s",m_NumReadsLoaded,strerror(errno));
		return(eBSFerrMem);
		}
	m_AllocdReadHitsIdx = m_NumReadsLoaded;
	}

if(SortMode != eRSMSeq || bSeqSorted)
	{
	pReadHit = NULL;
	tsReadHit **pIdx = m_ppReadHitsIdx;
	while((pReadHit = IterReads(pReadHit))!=NULL)
		*pIdx++ = pReadHit;
	}
else
	{
	pReadHit = m_pReadHits;
	size_t BuffOfs = 0;
	for(Idx = 0; Idx < m_NumReadsLoaded; Idx++)
		{
		m_ppReadHitsIdx[Idx] = pReadHit;
		BuffOfs += sizeof(tsReadHit) + pReadHit->ReadLen + pReadHit->DescrLen;
		pReadHit = (tsReadHit *)((char *)m_pReadHits+BuffOfs);
		}
	}

switch(SortMode) {
	case eRSMReadID:
		m_mtqsort.qsort(m_ppReadHitsIdx,m_NumReadsLoaded,sizeof(tsReadHit *),SortReadIDs);
		break;
	case eRSMPairReadID:
		m_mtqsort.qsort(m_ppReadHitsIdx,m_NumReadsLoaded,sizeof(tsReadHit *),SortPairReadIDs);
		break;
	case eRSMHitMatch:
		m_mtqsort.qsort(m_ppReadHitsIdx,m_NumReadsLoaded,sizeof(tsReadHit *),SortHitMatch);
		break;

	case eRSMPEHitMatch:
		m_mtqsort.qsort(m_ppReadHitsIdx,m_NumReadsLoaded,sizeof(tsReadHit *),SortPEHitMatch);
		break;		

	case eRSMSeq:
		if(!bSeqSorted)
			m_mtqsort.qsort(m_ppReadHitsIdx,m_NumReadsLoaded,sizeof(tsReadHit *),SortReadSeqs);
		break;
	default:
		break;
	}

// m_ppReadHitsIdx now in requested order, assign sequentially incrementing ReadHitIdx to the reads
for(Idx = 1; Idx <= m_NumReadsLoaded; Idx++)
	m_ppReadHitsIdx[Idx-1]->ReadHitIdx = Idx;

m_CurReadsSortMode = SortMode;
return(eBSFSuccess);
}


// SortReadIDs
// Sort reads by ascending read identifiers
int
CAligner::SortReadIDs(const void *arg1, const void *arg2)
{
tsReadHit *pEl1 = *(tsReadHit **)arg1;
tsReadHit *pEl2 = *(tsReadHit **)arg2;

if(pEl1->ReadID < pEl2->ReadID )
		return(-1);
if(pEl1->ReadID > pEl2->ReadID )
	return(1);
return(0);
}

// SortPairReadIDs
// Sort paired reads reads by ascending PairReadID identifiers
int
CAligner::SortPairReadIDs(const void *arg1, const void *arg2)
{
tsReadHit *pEl1 = *(tsReadHit **)arg1;
tsReadHit *pEl2 = *(tsReadHit **)arg2;
UINT32 El1ID = pEl1->PairReadID & 0x7fffffff;
UINT32 El2ID = pEl2->PairReadID & 0x7fffffff;

if(El1ID < El2ID)
	return(-1);
if(El1ID > El2ID)
	return(1);
if(pEl1->PairReadID & 0x80000000)		// if same PairReadID then the 3' read is after the 5' read
	return(1);
return(0);
}

// SortPEHitmatch
// Sort by NAR eNARAccepted with FlgPEAligned, then ascending chrom,  then PairReadID, then 5' forward read
int
CAligner::SortPEHitMatch(const void *arg1, const void *arg2)
{
tsReadHit *pEl1 = *(tsReadHit **)arg1;
tsReadHit *pEl2 = *(tsReadHit **)arg2;
UINT32 El1ID;
UINT32 El2ID;

// order by eNARAccepted.FlgPEAligned
if((pEl1->NAR == eNARAccepted && pEl1->FlgPEAligned) && !(pEl2->NAR == eNARAccepted && pEl2->FlgPEAligned))
	return(-1);
if(!(pEl1->NAR == eNARAccepted && pEl1->FlgPEAligned) && (pEl2->NAR == eNARAccepted && pEl2->FlgPEAligned))
	return(1);
if(!(pEl1->NAR == eNARAccepted && pEl1->FlgPEAligned &&  pEl2->NAR == eNARAccepted && pEl2->FlgPEAligned))
	return(0);

// sort by Chrom ascending
if(pEl1->HitLoci.Hit.Seg[0].ChromID < pEl2->HitLoci.Hit.Seg[0].ChromID)
	return(-1);
if(pEl1->HitLoci.Hit.Seg[0].ChromID > pEl2->HitLoci.Hit.Seg[0].ChromID)
	return(1);

// sort by PairReadID ascending
El1ID = pEl1->PairReadID & 0x7fffffff;
El2ID = pEl2->PairReadID & 0x7fffffff;
if(El1ID < El2ID)
	return(-1);
if(El1ID > El2ID)
	return(1);
// same PairReadID, order by 5' followed by 3' read
if(pEl1->PairReadID & 0x80000000)		
	return(1);
return(0);
}


// SortHitmatch
// Sort by ascending read NAR,NumHits(1,0,2,3..), chrom, loci, len, strand, LowMMCnt
int
CAligner::SortHitMatch(const void *arg1, const void *arg2)
{
tsReadHit *pEl1 = *(tsReadHit **)arg1;
tsReadHit *pEl2 = *(tsReadHit **)arg2;

if(pEl1->NAR < pEl2->NAR)
	return(-1);
if(pEl1->NAR > pEl2->NAR)
	return(1);
if(pEl1->NumHits == 1 && pEl2->NumHits != 1)	// hit counts are to be sorted 1,0,2,3,...
	return(-1);
if(pEl1->NumHits != 1 && pEl2->NumHits == 1)
	return(1);
if(pEl1->NumHits != 1 && pEl2->NumHits != 1)
	{
	if(pEl1->NumHits < pEl2->NumHits)
		return(-1);
	if(pEl1->NumHits > pEl2->NumHits)
		return(1);
	return(0);
	}

if(pEl1->HitLoci.Hit.Seg[0].ChromID < pEl2->HitLoci.Hit.Seg[0].ChromID)
		return(-1);
if(pEl1->HitLoci.Hit.Seg[0].ChromID > pEl2->HitLoci.Hit.Seg[0].ChromID)
	return(1);
if(AdjStartLoci(&pEl1->HitLoci.Hit.Seg[0]) < AdjStartLoci(&pEl2->HitLoci.Hit.Seg[0]))
		return(-1);
if(AdjStartLoci(&pEl1->HitLoci.Hit.Seg[0]) > AdjStartLoci(&pEl2->HitLoci.Hit.Seg[0]))
	return(1);

if(AdjHitLen(&pEl1->HitLoci.Hit.Seg[0]) < AdjHitLen(&pEl2->HitLoci.Hit.Seg[0]))
		return(-1);
if(AdjHitLen(&pEl1->HitLoci.Hit.Seg[0]) > AdjHitLen(&pEl2->HitLoci.Hit.Seg[0]))
	return(1);

if(pEl1->HitLoci.Hit.Seg[0].Strand < pEl2->HitLoci.Hit.Seg[0].Strand)
		return(-1);
if(pEl1->HitLoci.Hit.Seg[0].Strand > pEl2->HitLoci.Hit.Seg[0].Strand)
	return(1);
if(pEl1->LowMMCnt < pEl2->LowMMCnt)
		return(-1);
if(pEl1->LowMMCnt > pEl2->LowMMCnt)
	return(1);
return(0);
}

// SortMultiHits
// Sort by ascending ChromID, HitLoci, HitLen, Hitmismatches, Strand, then ReadID
int
CAligner::SortMultiHits(const void *arg1, const void *arg2)
{
tsReadHit *pEl1 = (tsReadHit *)arg1;
tsReadHit *pEl2 = (tsReadHit *)arg2;

if(pEl1->HitLoci.Hit.Seg[0].ChromID < pEl2->HitLoci.Hit.Seg[0].ChromID )
		return(-1);
if(pEl1->HitLoci.Hit.Seg[0].ChromID > pEl2->HitLoci.Hit.Seg[0].ChromID )
	return(1);
if(AdjStartLoci(&pEl1->HitLoci.Hit.Seg[0]) < AdjStartLoci(&pEl2->HitLoci.Hit.Seg[0]))
		return(-1);
if(AdjStartLoci(&pEl1->HitLoci.Hit.Seg[0]) > AdjStartLoci(&pEl2->HitLoci.Hit.Seg[0]))
	return(1);
if(AdjHitLen(&pEl1->HitLoci.Hit.Seg[0]) < AdjHitLen(&pEl2->HitLoci.Hit.Seg[0]))
	return(-1);
if(AdjHitLen(&pEl1->HitLoci.Hit.Seg[0]) > AdjHitLen(&pEl2->HitLoci.Hit.Seg[0]))
	return(1);

if(pEl1->HitLoci.Hit.Seg[0].Mismatches < pEl2->HitLoci.Hit.Seg[0].Mismatches )
	return(-1);
if(pEl1->HitLoci.Hit.Seg[0].Mismatches > pEl2->HitLoci.Hit.Seg[0].Mismatches )
	return(1);

if(pEl1->HitLoci.Hit.Seg[0].Strand < pEl2->HitLoci.Hit.Seg[0].Strand )
		return(-1);
if(pEl1->HitLoci.Hit.Seg[0].Strand > pEl2->HitLoci.Hit.Seg[0].Strand )
	return(1);

if(pEl1->ReadID < pEl2->ReadID )
	return(-1);
if(pEl1->ReadID > pEl2->ReadID )
	return(1);
return(0);
}

// SortMultiHitReadIDs
// sort by ascending ReadID then by descending scores
// to break tied scores, further sorted by
int
CAligner::SortMultiHitReadIDs(const void *arg1, const void *arg2)
{
tsReadHit *pEl1 = (tsReadHit *)arg1;
tsReadHit *pEl2 = (tsReadHit *)arg2;

if(pEl1->ReadID < pEl2->ReadID)
		return(-1);
if(pEl1->ReadID > pEl2->ReadID)
	return(1);

if(pEl1->HitLoci.Hit.Score > pEl2->HitLoci.Hit.Score)
		return(-1);
if(pEl1->HitLoci.Hit.Score < pEl2->HitLoci.Hit.Score)
	return(1);

if(pEl1->HitLoci.Hit.Seg[0].ChromID < pEl2->HitLoci.Hit.Seg[0].ChromID )
		return(-1);
if(pEl1->HitLoci.Hit.Seg[0].ChromID > pEl2->HitLoci.Hit.Seg[0].ChromID )
	return(1);

if(AdjHitLen(&pEl1->HitLoci.Hit.Seg[0]) < AdjHitLen(&pEl2->HitLoci.Hit.Seg[0]))
	return(-1);
if(AdjHitLen(&pEl1->HitLoci.Hit.Seg[0]) > AdjHitLen(&pEl2->HitLoci.Hit.Seg[0]))
	return(1);

if(pEl1->HitLoci.Hit.Seg[0].Mismatches < pEl2->HitLoci.Hit.Seg[0].Mismatches )
	return(-1);
if(pEl1->HitLoci.Hit.Seg[0].Mismatches > pEl2->HitLoci.Hit.Seg[0].Mismatches )
	return(1);

if(AdjStartLoci(&pEl1->HitLoci.Hit.Seg[0]) < AdjStartLoci(&pEl2->HitLoci.Hit.Seg[0]))
		return(-1);
if(AdjStartLoci(&pEl1->HitLoci.Hit.Seg[0]) > AdjStartLoci(&pEl2->HitLoci.Hit.Seg[0]))
	return(1);

if(pEl1->HitLoci.Hit.Seg[0].Strand < pEl2->HitLoci.Hit.Seg[0].Strand )
		return(-1);
if(pEl1->HitLoci.Hit.Seg[0].Strand > pEl2->HitLoci.Hit.Seg[0].Strand )
	return(1);

return(0);
}

// SortSegJuncts
// sort by ascending chrom, start, end
int
CAligner::SortSegJuncts(const void *arg1, const void *arg2)
{
tsSegJuncts *pEl1 = (tsSegJuncts *)arg1;
tsSegJuncts *pEl2 = (tsSegJuncts *)arg2;

if(pEl1->ChromID < pEl2->ChromID)
		return(-1);
if(pEl1->ChromID > pEl2->ChromID)
	return(1);

if(pEl1->Starts < pEl2->Starts)
		return(-1);
if(pEl1->Starts > pEl2->Starts)
	return(1);
if(pEl1->Ends < pEl2->Ends)
		return(-1);
if(pEl1->Ends > pEl2->Ends)
	return(1);
return(0);
}

int
CAligner::SortReadSeqs(const void *arg1, const void *arg2)
{
int Idx;
UINT8 *pSeq1;
UINT8 *pSeq2;

tsReadHit *pEl1 = *(tsReadHit **)arg1;
tsReadHit *pEl2 = *(tsReadHit **)arg2;

pSeq1 = &pEl1->Read[pEl1->DescrLen+1];
pSeq2 = &pEl2->Read[pEl2->DescrLen+1];

for(Idx = 0; Idx < pEl1->ReadLen; Idx++, pSeq1++, pSeq2++)
	{
	if((*pSeq1 & 0x07) < (*pSeq2 & 0x07))
		return(-1);
	if((*pSeq1 & 0x07) > (*pSeq2 & 0x07))
		return(1);
	if(Idx >= pEl2->ReadLen)
		return(1);
	}
if(pEl1->ReadLen < pEl2->ReadLen)
	return(-1);
return(0);
}

int
CAligner::SortLociPValues(const void *arg1, const void *arg2)
{
tsLociPValues *pEl1 = (tsLociPValues *)arg1;
tsLociPValues *pEl2 = (tsLociPValues *)arg2;
if(pEl1->PValue < pEl2->PValue)
		return(-1);
if(pEl1->PValue > pEl2->PValue)
	return(1);
return(0);
}

int
CAligner::SortPValuesLoci(const void *arg1, const void *arg2)
{
tsLociPValues *pEl1 = (tsLociPValues *)arg1;
tsLociPValues *pEl2 = (tsLociPValues *)arg2;
if(pEl1->Loci < pEl2->Loci)
		return(-1);
if(pEl1->Loci > pEl2->Loci)
	return(1);
return(0);
}


int
CAligner::SortSiteRelScale(const void *arg1, const void *arg2)
{
tsOctSitePrefs *pEl1 = (tsOctSitePrefs *)arg1;
tsOctSitePrefs *pEl2 = (tsOctSitePrefs *)arg2;
if(pEl1->RelScale < pEl2->RelScale)
		return(-1);
if(pEl1->RelScale > pEl2->RelScale)
	return(1);
return(0);
}

int
CAligner::SortSiteRelOctamer(const void *arg1, const void *arg2)
{
tsOctSitePrefs *pEl1 = (tsOctSitePrefs *)arg1;
tsOctSitePrefs *pEl2 = (tsOctSitePrefs *)arg2;
if(pEl1->Octamer < pEl2->Octamer)
		return(-1);
if(pEl1->Octamer > pEl2->Octamer)
	return(1);
return(0);
}

int
CAligner::SortConstraintLoci(const void *arg1, const void *arg2)
{
tsConstraintLoci *pEl1 = (tsConstraintLoci *)arg1;
tsConstraintLoci *pEl2 = (tsConstraintLoci *)arg2;
if(pEl1->ChromID < pEl2->ChromID)
		return(-1);
if(pEl1->ChromID > pEl2->ChromID)
	return(1);

if(pEl1->StartLoci < pEl2->StartLoci)
		return(-1);
if(pEl1->StartLoci > pEl2->StartLoci)
	return(1);

if(pEl1->EndLoci < pEl2->EndLoci)
		return(-1);
if(pEl1->EndLoci > pEl2->EndLoci)
	return(1);
return(0);
}

void
CAligner::AcquireSerialise(void)
{
#ifdef _WIN32
WaitForSingleObject(m_hMtxIterReads,INFINITE);
#else
pthread_mutex_lock(&m_hMtxIterReads);
#endif
}

void
CAligner::ReleaseSerialise(void)
{
#ifdef _WIN32
ReleaseMutex(m_hMtxIterReads);
#else
pthread_mutex_unlock(&m_hMtxIterReads);
#endif
}

void
CAligner::AcquireSerialiseMH(void)
{
#ifdef _WIN32
WaitForSingleObject(m_hMtxMHReads,INFINITE);
#else
pthread_mutex_lock(&m_hMtxMHReads);
#endif
}

void
CAligner::ReleaseSerialiseMH(void)
{
#ifdef _WIN32
ReleaseMutex(m_hMtxMHReads);
#else
pthread_mutex_unlock(&m_hMtxMHReads);
#endif
}

void
CAligner::AcquireLock(bool bExclusive)
{
#ifdef _WIN32
if(bExclusive)
	AcquireSRWLockExclusive(&m_hRwLock);
else
	AcquireSRWLockShared(&m_hRwLock);
#else
if(bExclusive)
	pthread_rwlock_wrlock(&m_hRwLock);
else
	pthread_rwlock_rdlock(&m_hRwLock);
#endif
}

void
CAligner::ReleaseLock(bool bExclusive)
{

#ifdef _WIN32
if(bExclusive)
	ReleaseSRWLockExclusive(&m_hRwLock);
else
	ReleaseSRWLockShared(&m_hRwLock);
#else
pthread_rwlock_unlock(&m_hRwLock);
#endif
}

int
CAligner::ProcLoadReadFiles(tsLoadReadsThreadPars *pPars)
{
teBSFrsltCodes Rslt;
int *pRslt = pPars->pRslt;
int Idx;
char *pszInfile;

int NumInputFilesProcessed;

AcquireLock(true);
m_bAllReadsLoaded = false;
m_LoadReadsRslt = eBSFSuccess;
ReleaseLock(true);

// first try to load as pre-processed reads '.rds' (as generated by 'kangar' or 'genreads')
// if unable to load as a '.rds' then try to load as raw fasta/fastq

if(((Rslt = (teBSFrsltCodes)LoadReads(m_ppszPE1InputFiles[0])) < eBSFSuccess) && Rslt != eBSFerrNotBioseq )
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to load reads from %s",m_ppszPE1InputFiles[0]);
	AcquireLock(true);
	*pRslt = Rslt;
	m_bAllReadsLoaded = true;
	ReleaseLock(true);
	return(Rslt);
	}

if(Rslt != eBSFerrNotBioseq)
	{
	if(Rslt == 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"There were no reads in %s",m_ppszPE1InputFiles[0]);
		Rslt = eBSFerrNoEntries;
		}
	else
		Rslt = eBSFSuccess;
	AcquireLock(true);
	*pRslt = Rslt;
	m_bAllReadsLoaded = true;
	ReleaseLock(true);
	return(Rslt);
	}

if(m_hInFile != -1)
	{
	close(m_hInFile);
	m_hInFile = -1;
	}

AcquireLock(true);
m_bAllReadsLoaded = false;
memset(&m_FileHdr,0,sizeof(tsBSFRdsHdr));
m_FileHdr.Magic[0] = 'b'; m_FileHdr.Magic[1] = 'i'; m_FileHdr.Magic[2] = 'o'; m_FileHdr.Magic[3] = 'r';
m_FileHdr.Version = cBSFRdsVersion;
m_FileHdr.FlagsK = 1;
m_FileHdr.FlagsCS = m_bIsSOLiD;
m_FileHdr.FlagsPR = m_PEproc == ePEdefault ? 0 : 1;
m_FileHdr.PMode = (UINT8)0;
m_FileHdr.QMode = m_QMethod;
m_FileHdr.Trim5 = m_Trim5;
m_FileHdr.Trim3 = m_Trim3;
m_pReadHits = NULL;
m_AllocdReadHitsMem = 0;
m_NumReadsLoaded = 0;
m_UsedReadHitsMem = 0;
m_FinalReadID = 0;
m_NumDescrReads = 0;
ReleaseLock(true);

if(m_FileHdr.FlagsPR && (m_NumPE1InputFiles < 1 || (m_NumPE1InputFiles != m_NumPE2InputFiles)))
	{
	Rslt = eBSFerrParams;
	AcquireLock(true);
	*pRslt = Rslt;
	m_bAllReadsLoaded = true;
	ReleaseLock(true);
	return(Rslt);
	}

NumInputFilesProcessed = 0;
if(!m_FileHdr.FlagsPR)
	{
	CSimpleGlob glob(SG_GLOB_FULLSORT);
	for(Idx = 0; Idx < m_NumPE1InputFiles; Idx++)
		{
		glob.Init();
		if(glob.Add(m_ppszPE1InputFiles[Idx]) < SG_SUCCESS)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to glob '%s",m_ppszPE1InputFiles[Idx]);
			Rslt = eBSFerrOpnFile;
			AcquireLock(true);
			*pRslt = Rslt;
			m_bAllReadsLoaded = true;
			m_LoadReadsRslt = Rslt;
			ReleaseLock(true);
			return(Rslt);	// treat as though unable to open file
			}

		if(glob.FileCount() <= 0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to locate any source raw reads file matching '%s",m_ppszPE1InputFiles[Idx]);
			continue;
			}

		Rslt = eBSFSuccess;
		for (int FileID = 0; Rslt >= eBSFSuccess &&  FileID < glob.FileCount(); ++FileID)
			{
			pszInfile = glob.File(FileID);

			gDiagnostics.DiagOut(eDLInfo,gszProcName,"Loading and parsing reads from raw sequence file '%s'\n",pszInfile);
			Rslt = LoadRawReads(false,NumInputFilesProcessed+1,pszInfile,NULL);
			if(Rslt != eBSFSuccess)
				{
				if(m_TermBackgoundThreads == 0)
					gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: Load failed for input raw sequence file '%s'\n",pszInfile);
				AcquireLock(true);
				*pRslt = Rslt;
				m_bAllReadsLoaded = true;
				m_LoadReadsRslt = Rslt;
				ReleaseLock(true);
				return(Rslt);
				}
			NumInputFilesProcessed += 1;
			}
		}
	}
else
	{
	Rslt = eBSFSuccess;
	NumInputFilesProcessed = 0;
	for(int FileID = 0; Rslt >= eBSFSuccess &&  FileID < m_NumPE1InputFiles; ++FileID)
		{
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Loading and parsing reads from 5' PE1 raw sequence file '%s'",m_ppszPE1InputFiles[FileID]);
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Loading and parsing reads from 3' PE2 raw sequence file '%s'",m_ppszPE2InputFiles[FileID]);

		Rslt = LoadRawReads(true,NumInputFilesProcessed+1,m_ppszPE1InputFiles[FileID],m_ppszPE2InputFiles[FileID]);
		if(Rslt != eBSFSuccess)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: Load failed for input raw sequence PE files '%s' and '%s'\n",m_ppszPE1InputFiles[FileID],m_ppszPE2InputFiles[FileID]);
			AcquireLock(true);
			*pRslt = Rslt;
			m_bAllReadsLoaded = true;
			m_LoadReadsRslt = Rslt;
			ReleaseLock(true);
			return(Rslt);
			}
		NumInputFilesProcessed += 2;
		}
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"%d raw sequence file%s parsed and reads loaded for aligning", NumInputFilesProcessed, NumInputFilesProcessed == 1 ? " was" : "s were");
if(NumInputFilesProcessed == 0)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Nothing to do, no raw sequence files to be filtered");
	AcquireLock(true);
	*pRslt = Rslt;
	m_bAllReadsLoaded = true;
	m_LoadReadsRslt = eBSFerrNoEntries;
	ReleaseLock(true);
	return(Rslt);
	}

// can now update header with final numbers
AcquireLock(true);
m_FileHdr.NumRds = m_NumDescrReads;
m_FileHdr.OrigNumReads = m_NumDescrReads;
m_FileHdr.TotReadsLen = m_DataBuffOfs;
m_FinalReadID = m_NumDescrReads;
m_NumReadsLoaded = m_NumDescrReads;
m_LoadReadsRslt = m_NumReadsLoaded > 0 ? eBSFSuccess : eBSFerrNoEntries;
m_bAllReadsLoaded = true;
*pRslt = Rslt;
ReleaseLock(true);
return(Rslt);
}



int
CAligner::AddEntry(bool bIsPairRead,	// false if SE or PE1, true if this is the paired read PE2
		 UINT32 PairReadID,		// identifies partner of this read if paired read processing (0 if no partner read)
		 UINT8 FileID,			// identifies file from which this read was parsed
		 int DescrLen,			// length of following descriptor
		 char *pszReadDescr,	// copy of descriptor, used to pair reads with matching descriptors
		 int ReadLen,			// length of following read
		 UINT8 *pszReadBuff)	// packed read + phred score
{
UINT8 *pTmpAlloc;
tsReadHit *pReadHit;
size_t memreq;

if(m_pReadHits == NULL)
	{
	memreq = cDataBuffAlloc;
	AcquireSerialise();
	AcquireLock(true);
#ifdef _WIN32
	m_pReadHits = (tsReadHit *) malloc((size_t)memreq);
	if(m_pReadHits == NULL)
		{
		ReleaseLock(true);
		ReleaseSerialise();
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddEntry: Memory allocation of %lld bytes failed",(INT64)memreq);
		return(eBSFerrMem);
		}
#else
	// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
	m_pReadHits = (tsReadHit *)mmap(NULL,(size_t)memreq, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
	if(m_pReadHits == MAP_FAILED)
		{
		m_pReadHits = NULL;
		ReleaseLock(true);
		ReleaseSerialise();
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddEntry: Memory allocation of %lld bytes through mmap()  failed",(INT64)memreq,strerror(errno));
		return(eBSFerrMem);
		}
#endif
	m_AllocdReadHitsMem = memreq;
	m_DataBuffOfs = 0;
	ReleaseLock(true);
	ReleaseSerialise();
	}

// need to allocate more memory? NOTE: allowing margin of 1K
if((m_AllocdReadHitsMem - m_DataBuffOfs) < (sizeof(tsReadHit) +  ReadLen + DescrLen + 0x03ff))
	{
	memreq = m_AllocdReadHitsMem + cDataBuffAlloc;
	AcquireSerialise();
	AcquireLock(true);
#ifdef _WIN32
	pTmpAlloc = (UINT8 *) realloc(m_pReadHits,memreq);
	if(pTmpAlloc == NULL)
		{
#else
	pTmpAlloc = (UINT8 *)mremap(m_pReadHits,m_AllocdReadHitsMem,memreq,MREMAP_MAYMOVE);
	if(pTmpAlloc == MAP_FAILED)
		{
		pTmpAlloc = NULL;
#endif

		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddEntry: Memory reallocation to %lld bytes failed - %s",memreq,strerror(errno));
		ReleaseLock(true);
		ReleaseSerialise();
		return(eBSFerrMem);
		}
	m_pReadHits = (tsReadHit *)pTmpAlloc;
	m_AllocdReadHitsMem = memreq;
	ReleaseLock(true);
	ReleaseSerialise();
	}

pReadHit = (tsReadHit *)((UINT8 *)m_pReadHits + m_DataBuffOfs);
m_DataBuffOfs += sizeof(tsReadHit) + ReadLen + DescrLen;
memset(pReadHit,0,sizeof(tsReadHit));
pReadHit->PrevSizeOf = m_PrevSizeOf;
pReadHit->ReadID = ++m_NumDescrReads;
pReadHit->PairReadID = PairReadID;
if(bIsPairRead)
	pReadHit->PairReadID |= 0x80000000;
pReadHit->NumReads = 1;
pReadHit->HitLoci.Hit.Seg[0].Strand = '?';
pReadHit->ReadLen = ReadLen;
pReadHit->DescrLen = (UINT8)DescrLen;
if(DescrLen > 0)
	memcpy((char *)&pReadHit->Read[0],pszReadDescr,DescrLen+1);
else
	pReadHit->Read[0] = '\0';
memmove(&pReadHit->Read[DescrLen+1],pszReadBuff,ReadLen);
m_PrevSizeOf = (UINT32)sizeof(tsReadHit) + ReadLen + DescrLen;

// processing threads are only updated with actual number of loaded reads every 100K reads so as
// to minimise disruption to the actual aligner threads which will also be serialised through m_hMtxIterReads
if(m_NumDescrReads > 0 && (m_NumDescrReads - m_NumReadsLoaded) >= 100000)
	{
	AcquireSerialise();
	m_FinalReadID = m_NumDescrReads;
	m_NumReadsLoaded = m_NumDescrReads;
	ReleaseSerialise();
	}
return(eBSFSuccess);
}


double												// returned prob of read being error free
CAligner::GenProbErrFreeRead(int QSSchema,			// guestimated scoring schema - 0: no scoring, 1: Solexa, 2: Illumina 1.3+, 3: Illumina 1.5+, 4: Illumina 1.8+ or could be Sanger 
				   int ReadLen,						// read length
				   UINT8 *pQScores)					// Phred quality scores
{
int SeqOfs;
int Score;
double ProbErr;
double ProbNoReadErr;
if(QSSchema == 0 || ReadLen < 1 || pQScores == NULL || pQScores[0] == '\0')		// if can't score then have to assume the best - read is sequencing base call error free
	return(1.0);

ProbNoReadErr = 1.0;
for (SeqOfs = 0; SeqOfs < ReadLen; SeqOfs++, pQScores++)
	{
	switch(QSSchema) { // quality scoring schema 
		case 4: // Illumina 1.8+ or Sanger. Sanger is '!' (0) to 'J' (41) and Illumina is '#' (2) to 'J' (41)
			Score = (int)(*pQScores - (UINT8)'!');
			break;
		default:		// all other scoring schemas have 40 at 'h'
			if(*pQScores <= '@')
				Score = 0;
			else
				Score = (int)(*pQScores - (UINT8)'@');
			break;
			}
	if(Score < 0)			// force scores to always be in the range 0..41 --- shouldn't be outside 0..41 but some qscores have been observed to be outside the expected range!
		Score = 0;
	else
		if(Score > 41)
			Score = 41;

	ProbErr = 1.0 / pow(10.0,(double)Score/10.0);
	ProbNoReadErr *= (1.0 - ProbErr); 
	}
if(ProbNoReadErr < 0.0)
	ProbNoReadErr = 0.0;
else
	if(ProbNoReadErr > 1.0)
		ProbNoReadErr = 1.0;
return(ProbNoReadErr);
}

teBSFrsltCodes
CAligner::LoadRawReads(bool bIsPairReads,	// true if paired end processing - PE1 reads in pszPE1File and PE2 reads in pszPE2File
		  int FileID,						// uniquely identifies source file for PE1, FileID + 1 uniquely identifies PE2 file
		  char *pszPE1File,					// process PE1 reads from this file
		  char *pszPE2File)					// optionally process PE2 reads from this file
{
static int FileNamesOfs = 0;
teBSFrsltCodes Rslt;
int Idx;
bool bIsFastq;

UINT8 *pReadBuff;
UINT8 *pQualBuff;
UINT8 Qphred;

bool bPE1SimReads;
int PE1NumDescrReads;
int PE1DescrLen;
UINT8 szPE1DescrBuff[cMaxDescrLen];
int PE1ReadLen;
UINT8 szPE1ReadBuff[cMaxReadLen];
int PE1QualLen;
UINT8 szPE1QualBuff[cMaxReadLen];
int PE1NumReadsAccepted;
int PE1NumInvalValues;
int PE1NumUnsupportedBases;
int PE1NumUnderlength;

bool bPE2SimReads;
int PE2NumDescrReads;
int PE2DescrLen;
UINT8 szPE2DescrBuff[cMaxDescrLen];
int PE2ReadLen;
UINT8 szPE2ReadBuff[cMaxReadLen];
int PE2QualLen;
UINT8 szPE2QualBuff[cMaxReadLen];
int PE2NumReadsAccepted;
int PE2NumInvalValues;
int PE2NumUnsupportedBases;
int PE2NumUnderlength;

int ContamLen5PE1;
int ContamLen3PE1;
int ContamLen5PE2;
int ContamLen3PE2;

int NumContamLen5PE1;
int NumContamLen3PE1;
int NumContamLen5PE2;
int NumContamLen3PE2;

CFasta PE1Fasta;
CFasta PE2Fasta;

INT32 EstScoreSchema;				// guestimated scoring schema - 0: no scoring, 1: Solexa, 2: Illumina 1.3+, 3: Illumina 1.5+, 4: Illumina 1.8+ or could be Sanger 
INT32 EstSeqLen;
INT32 EstDescrLen;
UINT32 EstNumSeqs;
INT64 ReqAllocSize;

// try and guestimate memory requirements so these can be allocated upfront - will be realloc'd if guestimate is wrong
if((EstNumSeqs = (teBSFrsltCodes)PE1Fasta.FastaEstSizes(pszPE1File,NULL,NULL,&EstDescrLen,NULL,&EstSeqLen,&EstScoreSchema)) == 0)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to estimate memory requirements for file '%s'",pszPE1File);
	return(eBSFerrOpnFile);
	}
// ReqAllocSize assumes -
// a) EstNumSeqs is not accurate so adds another 100000 to reduce chance of subsequent realloc required
// b) descriptors will be trimmed to 1st whitespace
ReqAllocSize = (INT64)(EstNumSeqs + 100000) * (INT64)(sizeof(tsReadHit) + EstSeqLen + max(20,(EstDescrLen/2)));	
if(bIsPairReads)
	ReqAllocSize *= 2;

UINT32 PairReadID;

PairReadID = (m_NumDescrReads/2) + 1;				// if bIsPairReads then start paired reads identifiers from this value and increment after each read processed
if((Rslt=(teBSFrsltCodes)PE1Fasta.Open(pszPE1File,true))!=eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Load: Unable to open '%s' [%s] %s",pszPE1File,PE1Fasta.ErrText((teBSFrsltCodes)Rslt),PE1Fasta.GetErrMsg());
	return(Rslt);
	}
if(!m_FileHdr.NumFiles)
	FileNamesOfs = 0;

m_FileHdr.NumFiles += 1;
if((FileNamesOfs + strlen(pszPE1File) + 1) < sizeof(m_FileHdr.FileNames))
	{
	strcpy((char *)&m_FileHdr.FileNames[FileNamesOfs],pszPE1File);
	FileNamesOfs += (int)strlen(pszPE1File) + 1;
	}

if(m_bIsSOLiD != PE1Fasta.IsSOLiD())
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Load: Aligning for %s but reads file '%s' is in %s",m_bIsSOLiD ? "Colorspace" : "Basespace",pszPE1File,m_bIsSOLiD ? "Basespace" : "Colorspace");
	PE1Fasta.Close();
	return(eBSFerrCvrtType);
	}

bIsFastq = PE1Fasta.IsFastq();

if(bIsPairReads)	
	{
	if((Rslt=(teBSFrsltCodes)PE2Fasta.Open(pszPE2File,true))!=eBSFSuccess)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Load: Unable to open '%s' [%s] %s",pszPE2File,PE2Fasta.ErrText((teBSFrsltCodes)Rslt),PE2Fasta.GetErrMsg());
		PE1Fasta.Close();
		return(Rslt);
		}

	m_FileHdr.NumFiles += 1;
	if((FileNamesOfs + strlen(pszPE2File) + 1) < sizeof(m_FileHdr.FileNames))
		{
		strcpy((char *)&m_FileHdr.FileNames[FileNamesOfs],pszPE2File);
		FileNamesOfs += (int)strlen(pszPE2File) + 1;
		}

	if(m_bIsSOLiD != PE2Fasta.IsSOLiD())
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Load: Aligning for %s but reads file '%s' is in %s",m_bIsSOLiD ? "Colorspace" : "Basespace",pszPE2File,m_bIsSOLiD ? "Basespace" : "Colorspace");
		PE1Fasta.Close();
		PE2Fasta.Close();
		return(eBSFerrCvrtType);
		}

	if(bIsFastq != PE2Fasta.IsFastq())
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Load: Paired end files not of same type");
		PE1Fasta.Close();
		PE2Fasta.Close();
		return(eBSFerrCvrtType);
		}
	}

AcquireSerialise();
AcquireLock(true);
if(m_pReadHits == NULL)
	{
	m_AllocdReadHitsMem = (size_t)ReqAllocSize;
#ifdef _WIN32
	m_pReadHits = (tsReadHit *) malloc((size_t)m_AllocdReadHitsMem);
	if(m_pReadHits == NULL)
		{
		ReleaseLock(true);
		ReleaseSerialise();
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadReads: Concatenated sequences memory allocation of %lld bytes - %s",(INT64)m_AllocdReadHitsMem,strerror(errno));
		PE1Fasta.Close();
		m_AllocdReadHitsMem = 0;
		return(eBSFerrMem);
		}
#else
	// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
	m_pReadHits = (tsReadHit *)mmap(NULL,m_AllocdReadHitsMem, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
	if(m_pReadHits == MAP_FAILED)
		{
		ReleaseLock(true);
	    ReleaseSerialise();
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadReads: Concatenated sequences memory of %lld bytes through mmap()  failed - %s",(INT64)m_AllocdReadHitsMem,strerror(errno));
		PE1Fasta.Close();
		m_pReadHits = NULL;
		return(eBSFerrMem);
		}
#endif
	m_DataBuffOfs = 0;

	}
else
	{
	UINT8 *pDstSeq;
	size_t memreq;
	if((m_DataBuffOfs + ReqAllocSize + 0x0fffff) >= m_AllocdReadHitsMem)		// 1M as a small safety margin!
		{
		memreq = (size_t)(m_AllocdReadHitsMem + ReqAllocSize);
#ifdef _WIN32
		pDstSeq = (UINT8 *) realloc(m_pReadHits,memreq);
#else
		pDstSeq = (UINT8 *)mremap(m_pReadHits,m_AllocdReadHitsMem,memreq,MREMAP_MAYMOVE);
		if(pDstSeq == MAP_FAILED)
			pDstSeq = NULL;
#endif
		if(pDstSeq == NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddSeq: Memory re-allocation to %lld bytes - %s",memreq,strerror(errno));
			PE1Fasta.Close();
			return(eBSFerrMem);
			}
		m_AllocdReadHitsMem = memreq;
		m_pReadHits = (tsReadHit *)pDstSeq;
		}
	}
ReleaseLock(true);
ReleaseSerialise();

PE1NumUnsupportedBases = 0;
PE1NumDescrReads = 0;
PE1NumReadsAccepted = 0;
PE1NumInvalValues = 0;
PE1NumUnderlength = 0;
PE2NumUnsupportedBases = 0;
PE2NumDescrReads = 0;
PE2NumReadsAccepted = 0;
PE2NumInvalValues = 0;
PE2NumUnderlength = 0;

NumContamLen5PE1 = 0;
NumContamLen3PE1 = 0;
NumContamLen5PE2 = 0;
NumContamLen3PE2 = 0;

bPE1SimReads = false;
bPE2SimReads = false;
while((Rslt = (teBSFrsltCodes)(PE1ReadLen = PE1Fasta.ReadSequence(szPE1ReadBuff,sizeof(szPE1ReadBuff)-1,true,false))) > eBSFSuccess)
	{
	if(m_TermBackgoundThreads != 0)	// need to immediately self-terminate?
		break;
	PE1NumDescrReads += 1;
	if(PE1ReadLen == eBSFFastaDescr)		// just read a descriptor line which would be as expected for multifasta or fastq
		{
		PE1DescrLen = PE1Fasta.ReadDescriptor((char *)szPE1DescrBuff,sizeof(szPE1DescrBuff)-1);
		szPE1DescrBuff[cMaxDescrLen-1] = '\0';
		PE1ReadLen = PE1Fasta.ReadSequence(szPE1ReadBuff,sizeof(szPE1ReadBuff)-1);
		if(PE1ReadLen < 0 || PE1ReadLen == eBSFFastaDescr || PE1ReadLen > cMaxReadLen)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Problem parsing sequence after %d reads parsed",PE1NumDescrReads);
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Last descriptor parsed: %s",szPE1DescrBuff);
			PE1Fasta.Close();
			if(bIsPairReads)
				PE2Fasta.Close();
			return(eBSFerrParse);
			}

			// check if these reads have been simulated, in which case the descriptor contains the loci of where the read was simulated from
		if(PE1NumDescrReads == 1)
			{
			if(!strncmp((char *)szPE1DescrBuff,"lcl|usimreads|",14) || !strncmp((char *)szPE1DescrBuff,"lcr|usimreads|",14))
				bPE1SimReads = true;
			else
				bPE1SimReads = false;
			}

		// if paired end processing then also load PE2 read
		if(bIsPairReads)
			{
			if((Rslt = (teBSFrsltCodes)(PE2ReadLen = PE2Fasta.ReadSequence(szPE2ReadBuff,sizeof(szPE2ReadBuff)-1,true,false))) <= eBSFSuccess)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Problem parsing sequence after %d reads parsed",PE2NumDescrReads);
				if(PE2NumDescrReads)
					gDiagnostics.DiagOut(eDLFatal,gszProcName,"Last descriptor parsed: %s",szPE2DescrBuff);
				PE1Fasta.Close();
				PE2Fasta.Close();
				return(eBSFerrParse);
				}

			if(m_TermBackgoundThreads != 0)	// need to immediately self-terminate?
				break;

			PE2NumDescrReads += 1;
			if(PE2ReadLen == eBSFFastaDescr)		// just read a descriptor line which would be as expected for multifasta or fastq
				{
				PE2DescrLen = PE2Fasta.ReadDescriptor((char *)szPE2DescrBuff,sizeof(szPE2DescrBuff)-1);
				szPE2DescrBuff[cMaxDescrLen-1] = '\0';
				PE2ReadLen = PE2Fasta.ReadSequence(szPE2ReadBuff,sizeof(szPE2ReadBuff)-1);
				if(PE2ReadLen < 0  || PE1ReadLen == eBSFFastaDescr  || PE1ReadLen > cMaxReadLen)
					{
					gDiagnostics.DiagOut(eDLFatal,gszProcName,"Problem parsing sequence after %d reads parsed",PE2NumDescrReads);
					gDiagnostics.DiagOut(eDLFatal,gszProcName,"Last descriptor parsed: %s",szPE2DescrBuff);
					PE1Fasta.Close();
					PE2Fasta.Close();
					return(eBSFerrParse);
					}
				}
			else
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"raw sequence file '%s' processing error: %s ",pszPE2File,
													Rslt > 0 ? "over length read" : "not a multifasta short reads or fastq file");
				PE1Fasta.Close();
				if(bIsPairReads)
					PE2Fasta.Close();
				return(eBSFerrParse);
				}

					// check if these reads have been simulated, in which case the descriptor contains the loci of where the read was simulated from
			if(PE2NumDescrReads == 1)
				{
				if(!strncmp((char *)szPE2DescrBuff,"lcl|usimreads|",14) || !strncmp((char *)szPE2DescrBuff,"lcr|usimreads|",14))
					bPE2SimReads = true;
				else
					bPE2SimReads = false;
				}
			}


		if(m_pContaminants != NULL)
			{
			// currently treating any contaminant match errors as if simply there was no overlap - should really report these!!!
			if((ContamLen5PE1 = m_pContaminants->MatchContaminants(eAOF5PE1Targ,1,m_Trim5+1,PE1ReadLen,szPE1ReadBuff)) <= m_Trim5)
				ContamLen5PE1 = 0;
			else
				{
				if(ContamLen5PE1 > m_Trim5)
					ContamLen5PE1 -= m_Trim5;
				}
			if((ContamLen3PE1 = m_pContaminants->MatchContaminants(eAOF3PE1Targ,1,m_Trim3+1,PE1ReadLen,szPE1ReadBuff)) <= m_Trim3)
				ContamLen3PE1 = 0;
			else
				{
				if(ContamLen3PE1 > m_Trim3)
					ContamLen3PE1 -= m_Trim3;
				}

			if(bIsPairReads)
				{
				if((ContamLen5PE2 = m_pContaminants->MatchContaminants(eAOF5PE2Targ,1,m_Trim5+1,PE2ReadLen,szPE2ReadBuff)) <= m_Trim5)
					ContamLen5PE2 = 0;
				else
					{
					if(ContamLen5PE2 > m_Trim5)
						ContamLen5PE2 -= m_Trim5;
					}
				if((ContamLen3PE2 = m_pContaminants->MatchContaminants(eAOF3PE2Targ,1,m_Trim3+1,PE2ReadLen,szPE2ReadBuff)) <= m_Trim3)
					ContamLen3PE2 = 0;
				else
					{
					if(ContamLen3PE2 > m_Trim3)
						ContamLen3PE2 -= m_Trim3;
					}
				}
			else
				{
				ContamLen5PE2 = 0;
				ContamLen3PE2 = 0;
				}
			}
		else
			{
			ContamLen5PE1 = 0;
			ContamLen3PE1 = 0;
			ContamLen5PE2 = 0;
			ContamLen3PE2 = 0;
			}

		// ensure would still have a sequence of at least cMinSeqLen after any end trims were applied
		if((m_Trim5 + m_Trim3 + ContamLen5PE1 + ContamLen3PE1 + (int)cMinSeqLen) > PE1ReadLen)
			{
			PE1NumUnderlength += 1;
			if(PE1NumUnderlength <= 10)
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Load: under length (%d) sequence in '%s' after end trims has been sloughed..",PE1ReadLen,pszPE1File);
			continue;
			}

		if(bIsPairReads)
			{
			if((m_Trim5 + m_Trim3 + ContamLen5PE2 + ContamLen3PE2  + (int)cMinSeqLen) > PE2ReadLen)
				{
				PE2NumUnderlength += 1;
				if(PE2NumUnderlength <= 10)
					gDiagnostics.DiagOut(eDLFatal,gszProcName,"Load: under length (%d) sequence in '%s' after end trims has been sloughed..",PE2ReadLen,pszPE2File);
				continue;
				}
			}

		if(bIsFastq && m_QMethod != eFQIgnore)
			{
			PE1QualLen = PE1Fasta.ReadQValues((char *)szPE1QualBuff,sizeof(szPE1QualBuff)-1);
			if(PE1QualLen != PE1ReadLen)		// must be same...
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Load: quality length (%d) not same as read length (%d) for '%s' entry in file '%s'",PE1QualLen,PE1ReadLen,szPE1DescrBuff,pszPE1File);
				PE1Fasta.Close();
				if(bIsPairReads)
					PE2Fasta.Close();
				return(eBSFerrParse);
				}
			// normalise the quality score to be in range 0..15 (needs to fit into 4 bits!)
			pQualBuff = szPE1QualBuff;
			for(Idx = 0; Idx < PE1ReadLen; Idx++,pQualBuff++)
				{
				switch(m_QMethod) {
					case eFQIgnore:		// simply treat as the minimum phred
						Qphred = 0;
						break;

					case eFQSanger:		// Qphred = -10 log10(P), where Qphred is in range 0..93; Illumina 1.8+ is essentially the same as Sanger
						if(*pQualBuff < 33 || *pQualBuff >= 126)
							{
							if(!PE1NumInvalValues++)
								gDiagnostics.DiagOut(eDLWarn,gszProcName,"Load: quality value of '%c' outside range expected in read '%s' from file '%s' for Sanger",*(char *)pQualBuff,szPE1DescrBuff,pszPE1File);
							if(*pQualBuff < 33)
								*pQualBuff = 33;
							else
								*pQualBuff = 125;
							}
						Qphred = *pQualBuff - 33;	// Sanger encodes into ascii starting from decimal 33 '!'
						if(Qphred > 40)				// clamp at phred equiv to 0.0001
							Qphred = 40;
						break;

					case eFQIllumia:	//Qphred = -10 log10(P), where Qphred is in range 0..63
						if(*pQualBuff < 64 || *pQualBuff >= 126)
							{
							if(!PE1NumInvalValues++)
								gDiagnostics.DiagOut(eDLWarn,gszProcName,"Load: quality value of '%c' outside range expected in read '%s' from file '%s' for Illumina 1.3+",*(char *)pQualBuff,szPE1DescrBuff,pszPE1File);
							if(*pQualBuff < 64)
								*pQualBuff = 64;
							else
								*pQualBuff = 125;
							}

						Qphred = *pQualBuff - 64;	//Illumia encodes into ascii starting from decimal 64
						if(Qphred > 40)				// clamp at phred equiv to 0.0001
							Qphred = 40;
						break;

					case eFQSolexa:		// SolexaQ = -10 log10(P/(1-P)), where SolexaQ is in range -5 to 62, note the negative value
										// negative values will result if P > 0.5
										// $Q = 10 * log(1 + 10 ** (ord(SolexaQphred) - 64) / 10.0)) / log(10);
										// once Qphred is over about 15 then essentially same as Sanger and Illumina 1.3+ so
										// is it worth doing the full conversion????
						if(*pQualBuff < 59 || *pQualBuff >= 126)
							{
							if(!PE1NumInvalValues++)
								gDiagnostics.DiagOut(eDLWarn,gszProcName,"Load: quality value of '%c' outside range expected in read '%s' from file '%s' for Solexa/Illumina pre 1.3",*(char *)pQualBuff,szPE1DescrBuff,pszPE1File);
							if(*pQualBuff < 64)
								*pQualBuff = 64;
							else
								*pQualBuff = 125;
							}
						Qphred = *pQualBuff - 59;	// Solexa/Illumina encodes into ascii starting from decimal 59
						Qphred = (UINT8)(10 * log(1 + pow(10.0,((double)Qphred/10.0) / log(10.0))));	//
						if(Qphred > 40)				// clamp at phred equiv to 0.0001
							Qphred = 40;
						break;
					}
				*pQualBuff = (UINT8)((((UINT32)Qphred+2)*15)/40);
				}

			// pack the read and quality, read into the low order bits 0..3, quality into bits 4..7
			pQualBuff = szPE1QualBuff;
			pReadBuff = szPE1ReadBuff;
			for(Idx = 0; Idx < PE1ReadLen; Idx++,pQualBuff++,pReadBuff++)
				szPE1ReadBuff[Idx] |= *pQualBuff << 4;

			if(bIsPairReads)
				{
				PE2QualLen = PE2Fasta.ReadQValues((char *)szPE2QualBuff,sizeof(szPE2QualBuff)-1);
				if(PE2QualLen != PE2ReadLen)		// must be same...
					{
					gDiagnostics.DiagOut(eDLFatal,gszProcName,"Load: quality length (%d) not same as read length (%d) for '%s' entry in file '%s'",PE2QualLen,PE2ReadLen,szPE2DescrBuff,pszPE2File);
					PE1Fasta.Close();
					PE2Fasta.Close();
					return(eBSFerrParse);
					}
				// normalise the quality score to be in range 0..15 (needs to fit into 4 bits!)
				pQualBuff = szPE2QualBuff;
				for(Idx = 0; Idx < PE2ReadLen; Idx++,pQualBuff++)
					{
					switch(m_QMethod) {
						case eFQIgnore:		// simply treat as the minimum phred
							Qphred = 0;
							break;

						case eFQSanger:		// Qphred = -10 log10(P), where Qphred is in range 0..93
							if(*pQualBuff < 33 || *pQualBuff >= 126)
								{
								if(!PE2NumInvalValues++)
									gDiagnostics.DiagOut(eDLWarn,gszProcName,"Load: quality value of '%c' outside range expected in read '%s' from file '%s' for Sanger",*(char *)pQualBuff,szPE2DescrBuff,pszPE2File);
								if(*pQualBuff < 33)
									*pQualBuff = 33;
								else
									*pQualBuff = 125;
								}
							Qphred = *pQualBuff - 33;	// Sanger encodes into ascii starting from decimal 33 '!'
							if(Qphred > 40)				// clamp at phred equiv to 0.0001
								Qphred = 40;
							break;

						case eFQIllumia:	//Qphred = -10 log10(P), where Qphred is in range 0..63
							if(*pQualBuff < 64 || *pQualBuff >= 126)
								{
								if(!PE2NumInvalValues++)
									gDiagnostics.DiagOut(eDLWarn,gszProcName,"Load: quality value of '%c' outside range expected in read '%s' from file '%s' for Illumina 1.3+",*(char *)pQualBuff,szPE2DescrBuff,pszPE2File);
								if(*pQualBuff < 64)
									*pQualBuff = 64;
								else
									*pQualBuff = 125;
								}

							Qphred = *pQualBuff - 64;	//Illumia encodes into ascii starting from decimal 64
							if(Qphred > 40)				// clamp at phred equiv to 0.0001
								Qphred = 40;
							break;

						case eFQSolexa:		// SolexaQ = -10 log10(P/(1-P)), where SolexaQ is in range -5 to 62, note the negative value
											// negative values will result if P > 0.5
											// $Q = 10 * log(1 + 10 ** (ord(SolexaQphred) - 64) / 10.0)) / log(10);
											// once Qphred is over about 15 then essentially same as Sanger and Illumina 1.3+ so
											// is it worth doing the full conversion????
							if(*pQualBuff < 59 || *pQualBuff >= 126)
								{
								if(!PE2NumInvalValues++)
									gDiagnostics.DiagOut(eDLWarn,gszProcName,"Load: quality value of '%c' outside range expected in read '%s' from file '%s' for Solexa/Illumina pre 1.3",*(char *)pQualBuff,szPE2DescrBuff,pszPE2File);
								if(*pQualBuff < 64)
									*pQualBuff = 64;
								else
									*pQualBuff = 125;
								}
							Qphred = *pQualBuff - 59;	// Solexa/Illumina encodes into ascii starting from decimal 59
							Qphred = (UINT8)(10 * log(1 + pow(10.0,((double)Qphred/10.0) / log(10.0))));	//
							if(Qphred > 40)				// clamp at phred equiv to 0.0001
								Qphred = 40;
							break;
						}
					*pQualBuff = (UINT8)((((UINT32)Qphred+2)*15)/40);
					}

					// pack the read and quality, read into the low order bits 0..3, quality into bits 4..7
				pQualBuff = szPE2QualBuff;
				pReadBuff = szPE2ReadBuff;
				for(Idx = 0; Idx < PE2ReadLen; Idx++,pQualBuff++,pReadBuff++)
					szPE2ReadBuff[Idx] |= *pQualBuff << 4;
				}
			}

		// apply any end trims, note that because quality scores packed into same byte as the base then end trims also trim the quality scores...
		if(m_Trim5 || ContamLen5PE1)
			{
			PE1ReadLen -= (m_Trim5 + ContamLen5PE1);
			memmove(szPE1ReadBuff,&szPE1ReadBuff[m_Trim5 + ContamLen5PE1],PE1ReadLen);
			}
		if(m_Trim3 || ContamLen3PE1)
			PE1ReadLen -= (m_Trim3 + ContamLen3PE1);
		if(ContamLen5PE1 > 0)
			NumContamLen5PE1 += 1;
		if(ContamLen3PE1 > 0)
			NumContamLen3PE1 += 1;


		// truncate descriptors at 1st whitespace unless fasta was generated by simulating reads in which case
		// the whole descriptor is retained as where the read was simulated from is of interest
		if(!bPE1SimReads)	// if not simulated reads
			{
			for(Idx = 0; Idx < cMaxDescrIDLen-1; Idx++)
				{
				if(szPE1DescrBuff[Idx] == '\0' || isspace(szPE1DescrBuff[Idx]))
					break;
				}
			szPE1DescrBuff[Idx] = '\0';
			PE1DescrLen = Idx;
			}

		if(bIsPairReads)
			{
			// apply any end trims
			if(m_Trim5 || ContamLen5PE2)
				{
				PE2ReadLen -= (m_Trim5 + ContamLen5PE2);
				memmove(szPE2ReadBuff,&szPE2ReadBuff[m_Trim5 + ContamLen5PE2],PE2ReadLen);
				}
			if(m_Trim3 || ContamLen3PE2)
				PE2ReadLen -= (m_Trim3 + ContamLen3PE2);
			if(ContamLen5PE2 > 0)
				NumContamLen5PE2 += 1;
			if(ContamLen3PE2 > 0)
				NumContamLen3PE2 += 1;

			// truncate descriptors at 1st whitespace unless fasta was generated by simulating reads in which case
			// the whole descriptor is retained as where the read was simulated from is of interest
			if(!bPE2SimReads)	// if not simulated reads
				{
				for(Idx = 0; Idx < cMaxDescrIDLen-1; Idx++)
					{
					if(szPE2DescrBuff[Idx] == '\0' || isspace(szPE2DescrBuff[Idx]))
						break;
					}
				szPE2DescrBuff[Idx] = '\0';
				PE2DescrLen = Idx;
				}
			}

		if(m_TermBackgoundThreads != 0)	// need to immediately self-terminate?
			break;

		if(!bIsPairReads)
			{
			if((Rslt=(teBSFrsltCodes)AddEntry(false,0,FileID,PE1DescrLen,(char *)szPE1DescrBuff,PE1ReadLen,szPE1ReadBuff))!=eBSFSuccess)
				break;
			PE1NumReadsAccepted += 1;
			}
		else
			{
			if((Rslt=(teBSFrsltCodes)AddEntry(false,PairReadID,FileID,PE1DescrLen,(char *)szPE1DescrBuff,PE1ReadLen,szPE1ReadBuff))!=eBSFSuccess)
				break;
			PE1NumReadsAccepted += 1;
			if((Rslt=(teBSFrsltCodes)AddEntry(true,PairReadID,FileID+1,PE2DescrLen,(char *)szPE2DescrBuff,PE2ReadLen,szPE2ReadBuff))!=eBSFSuccess)
				break;
			PE2NumReadsAccepted += 1;
			PairReadID += 1;
			}
		continue;
		}
	else
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"raw sequence file '%s' processing error: %s ",pszPE1File,
											Rslt > 0 ? "over length read" : "not a multifasta short reads or fastq file");
		PE1Fasta.Close();
		if(bIsPairReads)
			PE2Fasta.Close();
		return(eBSFerrParse);
		}
	}
if(Rslt != eBSFSuccess)
	{
	if(m_TermBackgoundThreads == 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Errors processing file: %s ",pszPE1File);
		while(PE1Fasta.NumErrMsgs())
			gDiagnostics.DiagOut(eDLFatal,gszProcName,PE1Fasta.GetErrMsg());
		}
	PE1Fasta.Close();
	if(bIsPairReads)
		PE2Fasta.Close();
	return(Rslt);
	}
PE1Fasta.Close();
if(bIsPairReads)
	PE2Fasta.Close();

if(m_TermBackgoundThreads != 0)	// need to immediately self-terminate?
	return(eBSErrSession);

AcquireSerialise();
m_FinalReadID = m_NumDescrReads;
m_NumReadsLoaded = m_NumDescrReads;
ReleaseSerialise();

gDiagnostics.DiagOut(eDLInfo,gszProcName,"LoadReads: Total of %1.9d reads parsed and loaded from %s",PE1NumReadsAccepted,pszPE1File);
if(PE1NumInvalValues > 0)
	gDiagnostics.DiagOut(eDLWarn,gszProcName,"Load: total of %d unexpected quality values read from file '%s'",PE1NumInvalValues,pszPE1File);
if(PE1NumUnsupportedBases > 0)
	gDiagnostics.DiagOut(eDLWarn,gszProcName,"Load: total of %d unsupported bases read from file '%s'",PE1NumUnsupportedBases,pszPE1File);
if(PE1NumUnderlength > 0)
	gDiagnostics.DiagOut(eDLWarn,gszProcName,"Load: total of %d under length sequences sloughed from file '%s'",PE1NumUnderlength,pszPE1File);
if(m_pContaminants != NULL)
	{
	gDiagnostics.DiagOut(eDLWarn,gszProcName,"Load: total of %d sequences PE1 sequences were 5' contaminate trimmed",NumContamLen5PE1);
	gDiagnostics.DiagOut(eDLWarn,gszProcName,"Load: total of %d sequences PE1 sequences were 3' contaminate trimmed",NumContamLen3PE1);
	}

if(bIsPairReads)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"LoadReads: Total of %1.9d reads parsed and loaded from %s",PE2NumReadsAccepted,pszPE2File);
	if(PE2NumInvalValues > 0)
		gDiagnostics.DiagOut(eDLWarn,gszProcName,"Load: total of %d unexpected quality values read from file '%s'",PE2NumInvalValues,pszPE2File);
	if(PE2NumUnsupportedBases > 0)
		gDiagnostics.DiagOut(eDLWarn,gszProcName,"Load: total of %d unsupported bases read from file '%s'",PE2NumUnsupportedBases,pszPE2File);
	if(PE2NumUnderlength > 0)
		gDiagnostics.DiagOut(eDLWarn,gszProcName,"Load: total of %d under length sequences sloughed from file '%s'",PE2NumUnderlength,pszPE2File);
	if(m_pContaminants != NULL)
		{
		gDiagnostics.DiagOut(eDLWarn,gszProcName,"Load: total of %d sequences PE1 sequences were 5' contaminant trimmed",NumContamLen5PE2);
		gDiagnostics.DiagOut(eDLWarn,gszProcName,"Load: total of %d sequences PE1 sequences were 3' contaminant trimmed",NumContamLen3PE2);
		}
	}
return(eBSFSuccess);
}





