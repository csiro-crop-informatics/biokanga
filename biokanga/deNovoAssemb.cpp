// Copyright 2013, 2014 CSIRO  ( http://www.csiro.au/ ) m_bRawDedupe
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

#include "./biokanga.h"
#include "./Kangadna.h"
#include "./deNovoAssemb.h"


// relies on base classes constructors
CdeNovoAssemb::CdeNovoAssemb(void)
{
m_bProcPE = false;
m_bTermPass = false;
m_bSenseStrandOnly = false;
m_bSingleEnded = false;
m_EarlyOverlapTermThres = 0.0;
m_NReduceThresSteps = 0;
m_pAllocdThreadSeqs = NULL;
m_AllocdThreadSeqsSize = 0;
memset(m_ThreadSeqBlocks,0,sizeof(m_ThreadSeqBlocks));
}


// relies on base classes destructors
CdeNovoAssemb::~CdeNovoAssemb(void)
{
if(m_pAllocdThreadSeqs != NULL)
	{
#ifdef _WIN32
	free(m_pAllocdThreadSeqs);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pAllocdThreadSeqs != MAP_FAILED)
		munmap(m_pAllocdThreadSeqs,m_AllocdThreadSeqsSize);
#endif	
	m_pAllocdThreadSeqs = NULL;
	}
}

teBSFrsltCodes
CdeNovoAssemb::LoadSeqsOnly(bool bSenseStrandOnly,	// process sequences as strand specific
					bool bSingleEnded,				// treat all sequences as being single ended even if loaded as paired ends
					int OrientatePE,					// PE end orientations 0: sense/antisense, 1: sense/sense, 2: antisense/sense, 3: antisense/antisense 
					char *pszPE1File,				// input high confidence seed PE1 sequences file
					char *pszPE2File,				// input high confidence seed PE2 sequences file
					char *pszSeedContigsFile,		// input high confidence seed SE contigs file
					char *pszInArtReducfile)		// optional input preprocessed artefact reduced packed reads from this file
{
int Rslt;
m_bSenseStrandOnly = bSenseStrandOnly;		// sequences from sense strand specific
m_bSingleEnded = bSingleEnded;				// treat all sequences as being single ended even if loaded as paired ends

// firstly, if specified then load any high confidence seed contigs or SE fragments
// these are likely to be longer than any SE/PE reads and thus will take longer to process for overlaps
if(pszSeedContigsFile != NULL && pszSeedContigsFile[0] != '\0')
	if((Rslt = LoadSeedContigs(pszSeedContigsFile))  < eBSFSuccess)
		return((teBSFrsltCodes)Rslt);

// if specified then load packed artefact reduced reads
if(pszInArtReducfile != NULL && pszInArtReducfile[0] != '\0')
	if((Rslt = LoadPackedSeqsFromFile(pszInArtReducfile)) < eBSFSuccess)
		return((teBSFrsltCodes)Rslt);

// next, if specified then load filtered PE1 and PE2 reads
if(pszPE1File != NULL && pszPE1File[0] != '\0')
	if((Rslt = LoadSeedPEs(pszPE1File,pszPE2File,OrientatePE))  < eBSFSuccess)
		return((teBSFrsltCodes)Rslt);

// initialise header flags
gDiagnostics.DiagOut(eDLInfo,gszProcName,"AssembReads: Initialising sequence headers ...");
if(bSingleEnded)
	{
	UpdateAllSeqHeaderFlags(0,~(cFlgNonOverlap),false);			// any PEs are now single ended
	m_Sequences.bPESeqs = false;					
	}
else
	UpdateAllSeqHeaderFlags(0,~(cFlgSeqPE2 | cFlgSeqPE | cFlgNonOverlap),false);  // just retain flags for paired end or single end sequences 

if((Rslt=GenRdsSfx(1, 2)) < eBSFSuccess)
	return((teBSFrsltCodes)Rslt);

	// generate array of sequence starts plus array of flags from sequence headers
if((Rslt=GenSeqStarts(true)) < eBSFSuccess)
	return((teBSFrsltCodes)Rslt);
return(eBSFSuccess);
}


teBSFrsltCodes
CdeNovoAssemb::AssembReads(etdeNovoPMode PMode,		// processing mode, currently either eAMEAssemble (default), eAMESAssemble (stringent) or eAMQAssemble (quick)
					int TrimInputEnds,				// trim input sequences, both 5' and 3' ends by this many bases
				    int MinInputSeqLen,				// only accept for assembly sequences which are, after any trimming, of at least this length
					int TrimPE2SE,					// trim PEs both 5' and 3' ends by this many bases before treating as SE
					bool bSenseStrandOnly,			// process sequences as strand specific
					bool bAllowSE2PE,				// if true then if SE overlaps PE1 or PE2 ends singularly but not both then merge the SE with the relevant overlapped PE end and retain as a PE
					bool bSingleEnded,				// treat all sequences as being single ended even if loaded as paired ends
					int MaxPasses,					// limit number of de Novo assembly passes to this maximum (quick mode defaults to 10, exhaustive defaults to 100) set to 0 for no limit
					int OutPass2File,				// output partially assembled starting from this pass - 0 if only final is to be written to file
					int NReduceThresSteps,			// reduce overlap thresholds over this many steps (defaults: 3 quick, 5 standard, 8 stringent assemble)");
					int Subs1Kbp,					// allow this many induced substitutions per Kbp overlapping sequence fragments
					int MaxEnd12Subs,				// allow at the initial 12bp (2 x hexamer primer artefacts) of the 5' or 3' of overlap to have this many base mismatches in addition to the overall allowed MaxSubs 
					int InitSEOvlp,					// initial minimal SE overlap required to merge SEs
					int FinSEOvlp,					// final minimal SE overlap required to merge SEs
					int InitPEOvlp,					// initial minimal PE overlap required to merge PEs
					int FinPEOvlp,					// final minimal PE overlap required to merge PEs
					int MinPE2SEOvlp,				// minimal overlap of PE1 onto PE2 required to merge as SE
					int PE2SESteps,					// when less than this many steps remaining then treat PE1 and PE2 as individual SE sequences if excessive lengths (defaults to 2, set 0 to disable)");
					int OrientatePE,				// PE end orientations 0: sense/antisense, 1: sense/sense, 2: antisense/sense, 3: antisense/antisense 
					char *pszPE1File,				// optional input high confidence seed PE1 sequences file
					char *pszPE2File,				// optional input high confidence seed PE2 sequences file
					char *pszSeedContigsFile,		// optional input high confidence seed SE contigs file
					char *pszInArtReducfile,		// optional input preprocessed artefact reduced packed reads from this file
					char *pszAssembFragsFile)		// where to write assembled sequence fragments as contigs ("SE" appended)
{
int Rslt;
int CurPass;				// incremented every merge pass over sequences
int CurMinReqPEPrimOverlap;	// if primary probe is overlapping onto a PE then the initial primary overlap (onto PE1 or PE2) must be of at least this length
int CurMinReqPESecOverlap;	// if primary probe was overlapping onto a PE then the secondary probe overlap (onto PE1 or PE2) must be of at least this length
int CurMinReqPESumOverlap;	// if primary probe was overlapping onto a PE then the sum of the PE1 and PE2 overlap must be of at least this length
int CurMinReqSEPrimOverlap;	// if primary probe is overlapping onto a SE then the overlap must be of at least this length
int CurMinPEMergeOverlap;	// if probe PE1 and probe PE2 being considered for merging then there must be an overlap of at least this many bases

int	CurMinPETotSeqLen2SE;	// cuurent minimum total of PE1 and PE2 end sequence lengths
int	CurMinPESeqLen2SE;      // current minimum of either PE1 or PE2 end sequence lengths

UINT32 CurTotNumPEs;	// total number of paired ends
int CurPE1MinLen;		// returned PE1 min length
int CurPE1MeanLen;		// returned PE1 mean length
int CurPE1MaxLen;		// returned PE1 max length
int CurPE2MinLen;		// returned PE2 min length
int CurPE2MeanLen;		// returned PE2 mean length
int CurPE2MaxLen;		// returned PE2 max length
UINT32 CurTotNumSEs;	// total number of single ends
int CurSEMinLen;		// returned SE min length
int CurSEMeanLen;		// returned SE mean length
int CurSEMaxLen;		// returned SE max length

int AllowedSubsKbp;		// allowed number of substutions per overlapping Kbp in current pass 
int AllowedEnd12Subs;	// allowed number of end primer artefact substitutions in current pass
bool bSubsIncr;			// set true when AllowedSubsKbp or AllowedEnd12Subs is incremented, false if not incremented
bool bUsingSubs;        // set true when using subs
double MergedPercentage; // percentage of sequences merged in this processing pass
bool bNewThres;			// set true when thresholds have been changed
bool bProcPE;			// set true if any PE processing, false if only SE
bool bAtMinThres;		// set true when processing with minimal thresholds
bool bAtInterThres;		// set true when processing with intermediate thresholds
int RemainingThresSteps;  // remaining threhold reduction steps		

UINT32 PrevNumPartialSeqs2Assemb;
UINT32 PrevNumSeqs2Assemb;

m_bSenseStrandOnly = bSenseStrandOnly;		// sequences from sense strand specific
m_bSingleEnded = bSingleEnded;				// treat all sequences as being single ended even if loaded as paired ends

m_TrimInputEnds = TrimInputEnds;
m_TrimPE2SE = TrimPE2SE;

// firstly, if specified then load any high confidence seed contigs or SE fragments
// these are likely to be longer than any SE/PE reads and thus will take longer to process for overlaps
if(pszSeedContigsFile != NULL && pszSeedContigsFile[0] != '\0')
	if((Rslt = LoadSeedContigs(pszSeedContigsFile,TrimInputEnds,MinInputSeqLen))  < eBSFSuccess)
		return((teBSFrsltCodes)Rslt);

// if specified then load packed artefact reduced reads
if(pszInArtReducfile != NULL && pszInArtReducfile[0] != '\0')
	if((Rslt = LoadPackedSeqsFromFile(pszInArtReducfile)) < eBSFSuccess)
		return((teBSFrsltCodes)Rslt);

// next, if specified then load filtered PE1 and PE2 reads
if(pszPE1File != NULL && pszPE1File[0] != '\0')
	if((Rslt = LoadSeedPEs(pszPE1File,pszPE2File,OrientatePE,TrimInputEnds,MinInputSeqLen))  < eBSFSuccess)
		return((teBSFrsltCodes)Rslt);

// initialise header flags
gDiagnostics.DiagOut(eDLInfo,gszProcName,"AssembReads: Initialising sequence headers ...");
if(bSingleEnded)
	{
	UpdateAllSeqHeaderFlags(0,~(cFlgNonOverlap),false);			// any PEs are now single ended
	m_Sequences.bPESeqs = false;					
	}
else
	UpdateAllSeqHeaderFlags(0,~(cFlgSeqPE2 | cFlgSeqPE | cFlgNonOverlap),false);  // just retain flags for paired end or single end sequences 

if(MaxPasses < 0)	// caller just wanted the sequences loaded?
	{
	if((Rslt=GenRdsSfx(1, 2)) < eBSFSuccess)
		return((teBSFrsltCodes)Rslt);

	// generate array of sequence starts plus array of flags from sequence headers
	if((Rslt=GenSeqStarts(true)) < eBSFSuccess)
		return((teBSFrsltCodes)Rslt);
	return(eBSFSuccess);
	}

// there will be multiple passes until assembly is deemed to have completed
CurPass = 0;
MergedPercentage = 100.0;			// ensure will not terminate on 1st pass without trying to merge! 
bNewThres = false;
bProcPE = false;
bAtInterThres = false;
bAtMinThres = false;
AllowedSubsKbp = 0;
AllowedEnd12Subs = 0;
bUsingSubs = false;
RemainingThresSteps = NReduceThresSteps;

CurMinPETotSeqLen2SE = (cMinPETotSeqLen2SE*3)/2;
CurMinPESeqLen2SE = (cMinPESeqLen2SE*3)/2;

while(m_Sequences.NumSeqs2Assemb) {       
	CurPass += 1;
	bSubsIncr = false;

	if((bUsingSubs && !bSubsIncr) || ((Subs1Kbp || MaxEnd12Subs) &&
		MergedPercentage <= 1.0 &&
		bAtMinThres == true))						
		{
		bUsingSubs = true;
		if(AllowedSubsKbp < Subs1Kbp)
			{
			AllowedSubsKbp = Subs1Kbp;
			bSubsIncr = true;
			}
		if(AllowedEnd12Subs < MaxEnd12Subs)
			{
			AllowedEnd12Subs = MaxEnd12Subs;
			bSubsIncr = true;
			}
		}
	else
		{
		if(!bUsingSubs)
			{
			AllowedSubsKbp = 0;
			AllowedEnd12Subs = 0;
			}
		}

	// generate index - if no subs specified then index is required on just on initial SeqWrd of sequence, otherwise it is over 
	// the 1st 4 SeqWrds as this allows subs in the first 60 bases to be discovered
	// as an memory optimisation don't create index over the last 2 SeqWrds (could be between 16 and 30 bases in these)     
	if((Rslt=GenRdsSfx(AllowedSubsKbp == 0 && AllowedEnd12Subs == 0 ? 1 : 4, 2)) < eBSFSuccess)
		return((teBSFrsltCodes)Rslt);

	// generate array of sequence starts plus array of flags from sequence headers
	if((Rslt=GenSeqStarts(true)) < eBSFSuccess)
		return((teBSFrsltCodes)Rslt);

#ifdef _DEBUG
#ifdef _WIN32
	ValidateSeqs2AssembStarts(); // only bother with validation whist debugging
#endif
#endif

	// get current min, mean, max lengths for PE1, PE2 and SEs
	GetSeqLenDist(&CurTotNumPEs,
			  &CurPE1MinLen,&CurPE1MeanLen,&CurPE1MaxLen,
			  &CurPE2MinLen,&CurPE2MeanLen,&CurPE2MaxLen,
			  &CurTotNumSEs,
			  &CurSEMinLen,&CurSEMeanLen,&CurSEMaxLen);

	// on the first pass establish the assembly thresholds
	// thresholds are proportional to read lengths and as to if a quick assembly has been requested
	// if high stringency requested then the initial and final thresholds are both increased by 10% 
	if(CurPass == 1)
		{
		m_NReduceThresSteps = NReduceThresSteps;
	
		if(CurTotNumPEs)							// if there is at least one PE to be processed
			{
			bProcPE = true;							// PE processing so initialise thresholds
			m_MinReqPESumOverlap = FinPEOvlp;
			m_MinReqPEPrimOverlap = max((FinPEOvlp - 5) / 2,cMinPEOvlp/2);
		
			m_MinReqPESecOverlap = m_MinReqPEPrimOverlap;

			m_InitialReqPESumOverlap = InitPEOvlp;
			m_InitialReqPEPrimOverlap = max((InitPEOvlp - 5) / 2,cMinPEOvlp/2);			
			m_InitialReqPESecOverlap = m_InitialReqPEPrimOverlap;

			m_MinPE2SEOvlp = MinPE2SEOvlp;
			}
		else    // else it's all SE processing
			{
			bProcPE = false;						
			m_MinReqPEPrimOverlap = 0;
			m_MinReqPESecOverlap = 0;
			m_MinReqPESumOverlap = 0;
			m_InitialReqPEPrimOverlap = 0;
			m_InitialReqPESecOverlap = 0;
			m_InitialReqPESumOverlap = 0;
			m_MinPE2SEOvlp = 0;
			}

		m_InitialReqSEPrimOverlap = InitSEOvlp;
		m_MinReqSEPrimOverlap = FinSEOvlp;

		bAtMinThres = false;
		bAtInterThres = false;

		if(bProcPE)
			{
			m_MinReqPESepDist = cMinReqPESepDist;
			m_MaxReqPESepDist = cMaxReqPESepDist;
			m_MinReqMergeLen = 20 + min(CurPE1MinLen,CurPE2MinLen);
			}
		else
			{
			m_MinReqPESepDist = 0;
			m_MaxReqPESepDist = 0;
			m_MinReqMergeLen = 0;
			}


		CurMinReqPEPrimOverlap = m_InitialReqPEPrimOverlap;		
		CurMinReqPESecOverlap = m_InitialReqPESecOverlap;		
		CurMinReqPESumOverlap = m_InitialReqPESumOverlap;		
		CurMinReqSEPrimOverlap = m_InitialReqSEPrimOverlap;		
		CurMinPEMergeOverlap = m_MinPE2SEOvlp;

		bNewThres = true;
		m_bProcPE = bProcPE;
		}

		// let user know the current breakdown on the PE/SE numbers/lengths plus current assembly parameters for this pass
	if(bProcPE)
		{
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"AssembReads: Num PEs: %u PE1 lens - min: %u mean: %u max: %u  ",CurTotNumPEs,CurPE1MinLen,CurPE1MeanLen,CurPE1MaxLen);
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"AssembReads: Num PEs: %u PE2 lens - min: %u mean: %u max: %u  ",CurTotNumPEs,CurPE2MinLen,CurPE2MeanLen,CurPE2MaxLen);
		}
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"AssembReads: Num SEs: %u  SE lens - min: %u mean: %u max: %u  ",CurTotNumSEs,CurSEMinLen,CurSEMeanLen,CurSEMaxLen);



	// if unable to merge few new sequences and merge thresholds already reduced to their minimum then assembly is complete
	if(CurPass > 1  &&								// do at least one pass
		bSubsIncr == false &&						// and AllowedSubsKbp hasn't been incremented for this pass
		AllowedSubsKbp == Subs1Kbp &&				// with allowed subs at their maximum
		AllowedEnd12Subs == MaxEnd12Subs &&			// with allowed end primer artefact subs at their maximum 
		bNewThres == false &&						// always do at least one pass with new thresholds!
		(MergedPercentage < (PMode == eAMQAssemble ?  0.1 : PMode == eAMEAssemble ? 0.05 : 0.025))	&&	// terminate assembly merge thresholds
		bAtMinThres == true)
		break;

	
	if(MaxPasses > 0 && CurPass > MaxPasses)   // number of passes is limited
		break;

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"AssembReads: Starting pass %d allowing %d substitutions per overlapping Kbp, overlap end primer %d substitutions, processing total of %u sequences with total length %llu...",CurPass,AllowedSubsKbp,AllowedEnd12Subs,m_Sequences.NumSeqs2Assemb,m_Sequences.Seqs2AssembLen);
	if(bProcPE)
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"AssembReads: MinReqPEPrimOverlap: %d, MinReqPESecOverlap: %d, MinReqPESumOverlap: %d, MinReqSEPrimOverlap :%d, MinPEMergeOverlap: %d, MinReqPESepDist: %d, MaxReqPESepDist: %d",
								CurMinReqPEPrimOverlap,CurMinReqPESecOverlap,CurMinReqPESumOverlap,CurMinReqSEPrimOverlap,CurMinPEMergeOverlap,m_MinReqPESepDist,m_MaxReqPESepDist);
	else
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"AssembReads: MinReqSEPrimOverlap :%d",CurMinReqSEPrimOverlap);

	// for diagnostics it can be useful to check on the partially assembled sequences...
	// user can request for these to be output to file starting from specified pass
	if(OutPass2File > 0 && CurPass >= OutPass2File)
		Rslt = SaveAssembSeqs(pszAssembFragsFile,CurPass,2000);

	if(bSenseStrandOnly)
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"AssembReads: Processing for Sense overlapping onto sense ...");
	else
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"AssembReads: Processing for Sense and antisense overlapping onto sense ...");

	// do some real work - locate overlaps and extend sequences with current overlap thresholds
	m_NextProcSeqID = 0;
	if(bNewThres == false)
		m_EarlyOverlapTermThres = PMode == eAMQAssemble ? 1000.0 : PMode == eAMEAssemble ? 500.0 : 250.0;  // early terminate pass if gaining few additional merges over 3 minutes
	else
		m_EarlyOverlapTermThres = 0.0;			// no early termination of overlap processing if a new threshold is being used
	bNewThres = false;

	Rslt = BuildOverlapExtensions(CurPass,false,bAllowSE2PE,AllowedSubsKbp,AllowedEnd12Subs,CurPE1MinLen,CurPE2MinLen,CurSEMinLen,CurMinReqPEPrimOverlap,CurMinReqPESecOverlap,CurMinReqPESumOverlap,CurMinReqSEPrimOverlap,CurMinPEMergeOverlap,MinPE2SEOvlp,m_MinReqPESepDist,m_MaxReqPESepDist,TrimPE2SE);
	if(Rslt < eBSFSuccess)
		return((teBSFrsltCodes)Rslt);

	if(!bSenseStrandOnly)
		{
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"AssembReads: Reverse complementing sequences ready for processing sense overlapping onto antisense ...");
		// reverse complement all sequences including PE's
		PackedRevCplAllIncPEs();

		// regenerate the index on the reverse complemented sequences
		if((Rslt=GenRdsSfx(AllowedSubsKbp == 0 && AllowedEnd12Subs == 0 ? 1 : 4,2)) < eBSFSuccess)
			return((teBSFrsltCodes)Rslt);

		// generate array of sequence starts but do not overwrite existing array of existing flags as these will have been updated during the overlap onto sense processing
		if((Rslt=GenSeqStarts(false)) < eBSFSuccess)
			return((teBSFrsltCodes)Rslt);

		// do some real work - locate overlaps and extend sequences with current overlap thresholds
		m_NextProcSeqID = 0;
		Rslt = BuildOverlapExtensions(CurPass,true,bAllowSE2PE,AllowedSubsKbp,AllowedEnd12Subs,CurPE1MinLen,CurPE2MinLen,CurSEMinLen,CurMinReqPEPrimOverlap,CurMinReqPESecOverlap,CurMinReqPESumOverlap,CurMinReqSEPrimOverlap,CurMinPEMergeOverlap,MinPE2SEOvlp,m_MinReqPESepDist,m_MaxReqPESepDist,TrimPE2SE);
		if(Rslt < eBSFSuccess)
			return((teBSFrsltCodes)Rslt);

		gDiagnostics.DiagOut(eDLInfo,gszProcName,"AssembReads: Reverse complementing sequences back to original sense ...");

		// reverse complement all sequences including PE's back to their original sense
		PackedRevCplAllIncPEs();
		// generate array of sequence starts but do not overwrite existing array of existing flags
		if((Rslt=GenSeqStarts(false)) < eBSFSuccess)
			return((teBSFrsltCodes)Rslt);

		gDiagnostics.DiagOut(eDLInfo,gszProcName,"AssembReads: Generated %d partial merges ...",m_NumPartialSeqs2Assemb);	// report number of newly merged sequences
		}

	// if at least one partial then need to iterate over all sequences and add those which have not been flagged as (cFlgAsmbSeed | cFlgAsmbExtn | cFlgAsmbCplt) to
	// the partial sequences ready for the next merge pass
	PrevNumPartialSeqs2Assemb = m_NumPartialSeqs2Assemb;
	PrevNumSeqs2Assemb = m_Sequences.NumSeqs2Assemb;
	if(PE2SESteps > 0 && RemainingThresSteps <= PE2SESteps)
		{
		if(RemainingThresSteps == PE2SESteps)
			Rslt = CombinePartialsWithUnmerged(false, CurMinPETotSeqLen2SE,CurMinPESeqLen2SE,TrimPE2SE);
		else
			Rslt = CombinePartialsWithUnmerged(false, 1,CurMinPESeqLen2SE,TrimPE2SE);	// time to force all PEs to be converted into SE sequences subject to the minimum sequence length requested and any trimming
		}
	else
		Rslt = CombinePartialsWithUnmerged(false,0,0,0);
	if(Rslt < 0)
		return((teBSFrsltCodes)Rslt);

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"AssembReads: Generated %d partial merges resulting in total sequence length %lld ...",m_NumPartialSeqs2Assemb,m_Sequences.Seqs2AssembLen);	// report number of newly merged sequences
	// what percentage of sequences were merged in current pass - this percentage is used to determine if the thresholds should be relaxed
	MergedPercentage = 100.0 * ((double)(PrevNumSeqs2Assemb - m_Sequences.NumSeqs2Assemb) / (double)PrevNumSeqs2Assemb);
	if(MergedPercentage < 0.0)		// possible if too many overlength PEs are being classified as SE relative to sequences merged
		MergedPercentage *= -1.0;

	if(bAtMinThres)					// if already at minimum thresholds then obviously can't reduce thresholds any further
		continue;

	// no change to thresholds whilst percentage of merged remains high
	if(MergedPercentage >=  2.0)			// if less than 2% change in number of sequences then progressively reduce thresholds
		continue;

	// progressively, in max of m_NReduceThresSteps  steps, reduce thresholds down towards the minimum
	CurMinReqPEPrimOverlap -= ((m_NReduceThresSteps-1) + m_InitialReqPEPrimOverlap - m_MinReqPEPrimOverlap) / m_NReduceThresSteps;
	CurMinReqPESecOverlap -= ((m_NReduceThresSteps-1) + m_InitialReqPESecOverlap - m_MinReqPESecOverlap) / m_NReduceThresSteps;
	CurMinReqPESumOverlap -= ((m_NReduceThresSteps-1) + m_InitialReqPESumOverlap - m_MinReqPESumOverlap) / m_NReduceThresSteps;
	CurMinReqSEPrimOverlap -= ((m_NReduceThresSteps-1) + m_InitialReqSEPrimOverlap - m_MinReqSEPrimOverlap) / m_NReduceThresSteps;

	CurMinPETotSeqLen2SE -= ((m_NReduceThresSteps-1) + ((cMinPETotSeqLen2SE*3)/2) - cMinPETotSeqLen2SE) / m_NReduceThresSteps;
	CurMinPESeqLen2SE -= ((m_NReduceThresSteps-1) + ((cMinPESeqLen2SE*3)/2) - cMinPESeqLen2SE) / m_NReduceThresSteps;

	// if within 10 bases of minimum threshold then don't continue steps and prolong the agony - use the minimum!
	if(CurMinReqPEPrimOverlap <= (m_MinReqPEPrimOverlap + 10))
		CurMinReqPEPrimOverlap = m_MinReqPEPrimOverlap;

	if(CurMinReqPESecOverlap <= (m_MinReqPESecOverlap + 10))
		CurMinReqPESecOverlap = m_MinReqPESecOverlap;

	if(CurMinReqPESumOverlap <= (m_MinReqPESumOverlap + 10))
		CurMinReqPESumOverlap = m_MinReqPESumOverlap;

	if(CurMinReqSEPrimOverlap <= (m_MinReqSEPrimOverlap + 10))
		CurMinReqSEPrimOverlap = m_MinReqSEPrimOverlap;

	if(CurMinPETotSeqLen2SE <= cMinPETotSeqLen2SE)
		CurMinPETotSeqLen2SE = cMinPETotSeqLen2SE;
	if(CurMinPESeqLen2SE <= cMinPESeqLen2SE)
		CurMinPESeqLen2SE = cMinPESeqLen2SE;

	// if still above minimum thresholds then set bAtInterThres otherwise set bAtMinThres
	if(CurMinReqPEPrimOverlap == m_MinReqPEPrimOverlap &&
	   CurMinReqPESecOverlap == m_MinReqPESecOverlap &&
	   CurMinReqPESumOverlap == m_MinReqPESumOverlap &&
	   CurMinReqSEPrimOverlap == m_MinReqSEPrimOverlap)
		{
		bAtMinThres = true;
		bAtInterThres = false;
		RemainingThresSteps = 0;
		}
	else
		{
		bAtInterThres = true;
		bAtMinThres = false;
		if(RemainingThresSteps > 0)
			RemainingThresSteps -= 1;
		}

	// thresholds now updated ready for next pass
	MergedPercentage = 100.0;
	bNewThres = true;
	}
return((teBSFrsltCodes)Rslt);
}


// PEs (PE1, PE2) can be converted to be processed as if two separate SE sequences if:
// a) MinPETotSeqLen2SE is > 0 and
//	a1) PE1 or PE2 length > 3/5ths of MinPETotSeqLen2SE
//	a2) PE1 and PE2 length totals to be at least MinPETotSeqLen2SE
// But if the length of the putative SE after any end trimming by TrimPE2SE would be less than MinPESeqLen2SE then that putative SE will be sloughed 
// The function caller can force all PEs to be converted into SEs by simply setting MinPETotSeqLen2SE to 1

// when PE processing initials PE state to 2, if 2 then PE1 is being processed, if 1 then PE2 is being processed, if 0 then a SE is being processed
// if PE1 being processed 

// Use PE1, use PE2      11   ;; if processing PEs, if bPE1 then process sequence and then reset bPE1     
// Use PE1, not PE2      10
// not PE1, use PE2      01
// not PE1, not PE2      00
//
int				// combine partial sequences back with unmerged sequences ready for next pass
CdeNovoAssemb::CombinePartialsWithUnmerged(bool bTrim15bp,						// if true then trim 15bp of the 5' ends to reduce the primer induced errors
									int MinPETotSeqLen2SE,						// combine partial sequences if length PE1 plus length PE2 is at least this length - 0 if no combining of PE1 and PE2s
								    int MinPESeqLen2SE,							// only combine if individual PE1 and PE2 lengths at least this length
									int TrimPE2SE)								// trim overlength non-overlapping PEs both 5' and 3' ends by this many bases before treating as SE
{
tSeqID SeqID;
UINT16 *pSeqFlags;
UINT16 SeqFlags;
tSeqWrd4 *pPE1SeqWrd;
UINT32 PE1SeqLen;
tSeqWrd4 *pPE2SeqWrd;
UINT32 PE2SeqLen;
tSeqWrd4 *pSeqWrd;
tSeqWrd4 *pSeqWrds;
tSeqWrd4 *pPackSeq;
int SeqLen;
int NumSeqWrds;

int NumSEs;
int NumPEs;
int NumPEsCvrt2SE;
int TrimPEEndsBy;
int CvtPEs2SE;
bool bPE1toSE;
bool bPE2toSE;
int TrimSeqLen;
int TrimNumSeqWrds;
tSeqWrd4 *pTrimSeqWrd;

INT64 PartialSeqsLen;

// if at least one partial then need to iterate over all sequences and add those which have not been flagged as (cFlgAsmbSeed | cFlgAsmbExtn | cFlgAsmbCplt) to
// the partial sequences ready for the next merge pass
if(m_NumPartialSeqs2Assemb)	// almost certainly there will be at least one but better check!
	{
	// identify those sequences which were not merged and treat these as if merged so they will be retained for next pass
	NumSEs = 0;
	NumPEs = 0;
	NumPEsCvrt2SE = 0;
	pSeqFlags = &m_Sequences.pSeqFlags[0];
	for(SeqID = 1; SeqID <= m_Sequences.NumSeqs2Assemb; SeqID++,pSeqFlags++)
		{
		if((SeqFlags = *pSeqFlags) & (cFlgAsmbSeed | cFlgAsmbExtn | cFlgAsmbCplt)) // can skip if sequence has been marked as already added to partials
			{
			if(SeqFlags & cFlgSeqPE)
				{
				SeqID += 1;
				pSeqFlags += 1;
				}
			continue;
			}
		// sequence not already added to partials
		pPE1SeqWrd = (tSeqWrd4 *)GetSeqHeader(SeqID,NULL,NULL,&PE1SeqLen,false);
		if(SeqFlags & cFlgSeqPE)
			{
			pPE2SeqWrd = (tSeqWrd4 *)GetSeqHeader(SeqID+1,NULL,NULL,&PE2SeqLen,false);
			SeqID += 1;
			pSeqFlags += 1;
			NumPEs += 1;
			}
		else
			{
			NumSEs += 1;
			pPE2SeqWrd = NULL;
			PE2SeqLen = 0;
			}

		// add to partials ready for next merge pass
		if((PartialSeqsLen = SavePartialSeqs(PE1SeqLen,pPE1SeqWrd,PE2SeqLen,pPE2SeqWrd)) < (INT64)0)
			return((int)PartialSeqsLen);
		}

	// all sequences to be retained for next pass are now in m_pPartialSeqs2Assemb, copy these back into m_Sequences.pSeqs2Assemb
	// at the same time, any PE sequences which have lengths cMinPESeqLen2SE and cMinPETotSeqLen2SE will be classified as being two separate SE sequences
	NumSEs = 0;
	NumPEs = 0;
	pSeqWrd = (tSeqWrd4 *)m_pPartialSeqs2Assemb;
	pSeqWrd+=1;
	
	pSeqWrds = (tSeqWrd4 *)m_Sequences.pSeqs2Assemb;
	*pSeqWrds++ = cSeqWrd4BOS;
	pPackSeq = pSeqWrds;
	m_Sequences.Seqs2AssembOfs = 1;
	m_Sequences.NumSeqs2Assemb = 0;
	m_Sequences.Seqs2AssembLen = 0;

	NumPEsCvrt2SE = 0;
	CvtPEs2SE = 0;
	TrimPEEndsBy = 0;

	for(SeqID = 1; SeqID <= m_NumPartialSeqs2Assemb; SeqID++)
		{
		SeqLen = *pSeqWrd++;
		if((tSeqWrd4)SeqLen == cSeqWrd4EOS)
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"AssembReads: unexpected EOS at SeqID = %d ...",SeqID);

		SeqFlags = SeqLen & 0x03;
		SeqLen >>= 2;
		SeqLen &= 0x3fffffff;

		// PEs (PE1, PE2) can be converted to be processed as if two separate SE sequences if:
		// a) MinPETotSeqLen2SE is > 0 and
		//	a1) PE1 or PE2 length > 3/5ths of MinPETotSeqLen2SE
		//	a2) PE1 and PE2 length totals to be at least MinPETotSeqLen2SE
		// But if the length of the putative SE after any end trimming by TrimPE2SE would be less than MinPESeqLen2SE then that putative SE will be sloughed 
		// The function caller can force all PEs to be converted into SEs by simply setting MinPETotSeqLen2SE to 1
		if((SeqFlags & cFlgSeqPE) && !(SeqFlags & cFlgSeqPE2))	// if PE1 then check if lengths for PE1 and/or PE2 are such that should be reclassified as being SE
			{
			UINT16 SeqFlagsPE2;
			tSeqWrd4 SeqLenPE2;
			NumSeqWrds = (SeqLen + 14) / 15;
			SeqLenPE2 = pSeqWrd[NumSeqWrds];
			SeqFlagsPE2 = SeqLenPE2 & 0x03;
			SeqLenPE2 >>= 2;
			SeqLenPE2 &= 0x3fffffff;

			// if sum of PE lengths sufficently long, or an individual PE more than 3/5ths of MinPETotSeqLen2SE then change to be two SE's
			if(MinPETotSeqLen2SE > 0 && ((SeqLen > (3 * MinPETotSeqLen2SE) / 5) || ((int)SeqLenPE2 > (3 * MinPETotSeqLen2SE) / 5) || ((SeqLen + (int)SeqLenPE2) >= MinPETotSeqLen2SE)))
				{
				TrimPEEndsBy = TrimPE2SE;			// always trim PE ends when accepting as SE sequences
				CvtPEs2SE = 2;
				bPE1toSE = true;
				bPE2toSE = true;
				SeqFlags = 0;						// PE1 is now a SE
				pSeqWrd[NumSeqWrds] &= 0xfffffffc;	// changes PE2 to also being single ended
				NumPEsCvrt2SE += 1;
				if((SeqLen - (2 * TrimPE2SE)) < MinPESeqLen2SE)   // if seqlen is too short then mark as to be discarded - assumes that PE1 end extension stopped because of errors in the sequence
					bPE1toSE = false;
				if(((int)SeqLenPE2 - (2 * TrimPE2SE)) < MinPESeqLen2SE) // if SeqlenPE2 is too short then mark as to be discarded - assumes that PE2 end extension stopped because of errors in the sequence
					bPE2toSE = false;
				}
			else
				{
				TrimPEEndsBy = 0;
				CvtPEs2SE = 0;
				bPE1toSE = false;
				bPE2toSE = false;
				}
			}

		NumSeqWrds = (SeqLen + 14) / 15;

		if(CvtPEs2SE == 0 || (CvtPEs2SE == 1 && bPE2toSE == true) || (CvtPEs2SE == 2 && bPE1toSE == true))
			{
			if(SeqFlags & cFlgSeqPE2)
				NumPEs += 1;
			else
				if(!(SeqFlags & cFlgSeqPE))
					NumSEs += 1;


			if(bTrim15bp)	// if trimming 5' by 15bp
				{
				TrimSeqLen = SeqLen - 15;
				TrimNumSeqWrds = NumSeqWrds-1;
				pTrimSeqWrd = pSeqWrd+1;
				}
			else
				{
				TrimSeqLen = SeqLen;
				TrimNumSeqWrds = NumSeqWrds;
				pTrimSeqWrd = pSeqWrd;
				}

			if(TrimPEEndsBy)
				{
				TrimSeqLen -= 2 * TrimPEEndsBy;
				TrimNumSeqWrds = (TrimSeqLen+14)/15;
				}
			m_Sequences.Seqs2AssembLen += TrimSeqLen;
			m_Sequences.NumSeqs2Assemb += 1;

			if((pPackSeq = (tSeqWrd4 *)SetSeqHeader(pPackSeq,m_Sequences.NumSeqs2Assemb,1,SeqFlags,TrimSeqLen,NULL))==NULL)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"CombinePartialsWithUnmerged: SetSeqHeader() failed");
				return(eBSFerrInternal);
				}

			if(TrimPEEndsBy)
				GetNHSeqWrdSubSeq(TrimPEEndsBy,TrimSeqLen,pSeqWrd,pPackSeq,false);
			else
				memcpy(pPackSeq,pTrimSeqWrd,sizeof(tSeqWrd4) * TrimNumSeqWrds); 
			m_Sequences.Seqs2AssembOfs += TrimNumSeqWrds;
			pPackSeq += TrimNumSeqWrds;
			}
		if(CvtPEs2SE > 0)
			CvtPEs2SE -= 1;
		if(CvtPEs2SE == 0)
			TrimPEEndsBy = 0;
		pSeqWrd += NumSeqWrds;
		}
	*pPackSeq = cSeqWrd4EOS;
	// m_Sequences.pSeqs2Assemb now ready for next pass

	m_NumPartialSeqs2Assemb = 0;
	m_LenPartialSeqs2Assemb = 0;
	*(tSeqWrd4 *)m_pPartialSeqs2Assemb = cSeqWrd4BOS;
	m_PartialSeqs2AssembOfs = 1;
	if(m_bProcPE)
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"CombinePartialsWithUnmerged: Reclassified %d PEs as SEs because of excessive length",NumPEsCvrt2SE);
	}

// if significantly reducing memory requirements then release this surplus memory
size_t memreq;
void *pAllocd;

memreq = ((15 * m_Sequences.Seqs2AssembOfs * sizeof(tSeqWrd4))/10);		// 50% to allow for a subsequent increase, unlikely but if SE converted to PEs then could happen
if((memreq * 2) < m_Sequences.AllocMemSeqs2Assemb)
	{
#ifdef _WIN32
	pAllocd = realloc(m_Sequences.pSeqs2Assemb,memreq);
#else
	pAllocd = mremap(m_Sequences.pSeqs2Assemb,m_Sequences.AllocMemSeqs2Assemb,memreq,MREMAP_MAYMOVE);
	if(pAllocd == MAP_FAILED)
		pAllocd = NULL;
#endif
	if(pAllocd == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"CombinePartialsWithUnmerged: Memory re-allocation to %d bytes - %s",memreq,strerror(errno));
		Reset(false);
		return(eBSFerrMem);
		}

	m_Sequences.pSeqs2Assemb = pAllocd;
	m_Sequences.AllocMemSeqs2Assemb = (UINT64)memreq;
	}

memreq =  m_Sequences.AllocMemSeqs2Assemb;	
if(memreq < m_AllocdPartialSeqs2Assemb)
	{
#ifdef _WIN32
	pAllocd = realloc(m_pPartialSeqs2Assemb,memreq);
#else
	pAllocd = mremap(m_pPartialSeqs2Assemb,m_AllocdPartialSeqs2Assemb,memreq,MREMAP_MAYMOVE);
	if(pAllocd == MAP_FAILED)
		pAllocd = NULL;
#endif
	if(pAllocd == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"CombinePartialsWithUnmerged: Memory re-allocation to %d bytes - %s",memreq,strerror(errno));
		Reset(false);
		return(eBSFerrMem);
		}

	m_pPartialSeqs2Assemb = pAllocd;
	m_AllocdPartialSeqs2Assemb = (UINT64)memreq;
	}
return(eBSFSuccess);
}

#ifdef _WIN32
unsigned __stdcall ThreadedOverlapExtend(void * pThreadPars)
#else
void * ThreadedOverlapExtend(void * pThreadPars)
#endif
{
int Rslt = 0;
tsThreadOverlapExtendPars *pPars = (tsThreadOverlapExtendPars *)pThreadPars; // makes it easier not having to deal with casts!
CdeNovoAssemb *pThis = (CdeNovoAssemb *)pPars->pThis;
Rslt = pThis->ProcOverlapExtend(pPars);
pPars->Rslt = Rslt;
#ifdef _WIN32
_endthreadex(0);
return(eBSFSuccess);
#else
pthread_exit(NULL);
#endif
}


const int cMinSeqsThreadBlock = 2;     // minimum number of sequences to be processed by a single thread - if less then interthread syncronisation cost may be too high

int
CdeNovoAssemb::InitGetSeqProc(int NumThreads)	// number of threads which will be requesting, via GetSeqProc(), sequence identifiers for processing
{
if(NumThreads < 1 || m_Sequences.NumSeqs2Assemb == 0)
	return(eBSFerrParams);

if(NumThreads > cMaxWorkerThreads)					// clamp number of threads
	NumThreads = cMaxWorkerThreads;

m_Sequences.NumProcessed = 0;	
memset(m_ThreadSeqBlocks,0,sizeof(m_ThreadSeqBlocks));
m_NxtSeqs2BlockAlloc = 1;		// allocate next thread sequence block starting with this sequence
m_SeqsPerThreadBlk = max(cMinSeqsThreadBlock,(m_Sequences.NumSeqs2Assemb + m_NumThreads - 1) / ((UINT32)m_NumThreads * 1000));	// nominal number of sequences per thread processing block
return(m_SeqsPerThreadBlk);
}

// GetSeqProc
// Returns identifier and flags for next sequence to be processed
// Identifers returned are for both PE1 and PE2, or for SE/contig 
int				// 0 if all returned, 1 if PE1 only, 2 if both PE1 and PE2
CdeNovoAssemb::GetSeqProc(tSeqID *pPE1SeqID,	// returned SE or PE1 sequence identifier
					  UINT32 *pPE1SeqFlags, // SE or PE1 flags
					  tSeqID *pPE2SeqID,	// if PE then PE2 sequence identifier
					  UINT32 *pPE2SeqFlags, // if PE then PE2 flags
					  int ThreadIdx)		// uniquely identifies calling thread, 1..NumThreads
{
tsSeqBlock *pSeqBlock;
UINT16 *pFlags;
int PE1Flags;
int PE2Flags;
bool bTermPass;
int Rslt;

Rslt = 0;
*pPE1SeqID = 0;			// default to no sequence identifier available
*pPE2SeqID = 0;
*pPE1SeqFlags = 0;
*pPE2SeqFlags = 0;
bTermPass = false;

AcquireSerialiseSeqFlags();
pSeqBlock = &m_ThreadSeqBlocks[ThreadIdx-1];

while(!bTermPass && Rslt == 0)		// whilst unable to locate sequences not already claimed or seeds  
	{
	if(pSeqBlock->NxtID == 0)		// when completed current block then try another block
		{
		// any unprocessed sequences remaining?
		if(m_NxtSeqs2BlockAlloc == 0 || m_NxtSeqs2BlockAlloc > m_Sequences.NumSeqs2Assemb)
			{
			m_NxtSeqs2BlockAlloc = 0;
			break;
			}
		
		// possible that the sequences near the end may require more processing as these are likely to be SE and thus longer than PE's
        // so start using smaller blocks of sequences per thread so as to more evenly spread the processing load and reduce possibilities of a single
        // thread at the end shouldering all the load.... 
		if(m_SeqsPerThreadBlk > cMinSeqsThreadBlock)			// do not to reduce the load for a thread below cMinSeqsThreadBlock, the cost per thread of sync'ing with other threads is not worth it...
			{
			if((m_Sequences.NumSeqs2Assemb - m_NxtSeqs2BlockAlloc) < m_SeqsPerThreadBlk * 50)
				m_SeqsPerThreadBlk = (m_SeqsPerThreadBlk * 2) / 3;
			if(m_SeqsPerThreadBlk < cMinSeqsThreadBlock)
				m_SeqsPerThreadBlk = cMinSeqsThreadBlock;
			}		
		if((m_Sequences.NumSeqs2Assemb - m_NxtSeqs2BlockAlloc) < (m_SeqsPerThreadBlk * 3)/2) // if next block after this would be < 50% of nominal size then make this current block the last
			m_SeqsPerThreadBlk = 1 + m_Sequences.NumSeqs2Assemb - m_NxtSeqs2BlockAlloc;

		pSeqBlock->StartID = m_NxtSeqs2BlockAlloc;
		pSeqBlock->NxtID = pSeqBlock->StartID;
		m_NxtSeqs2BlockAlloc += m_SeqsPerThreadBlk;
		pSeqBlock->NumSeqIDs = m_SeqsPerThreadBlk;
		pSeqBlock->EndID = m_NxtSeqs2BlockAlloc - 1;
		pFlags = &m_Sequences.pSeqFlags[pSeqBlock->EndID-1];

		if(*pFlags & cFlgSeqPE && !(*pFlags & cFlgSeqPE2))		// if would end on a PE 5' then adjust end to be on the PE 3'
			{
			m_NxtSeqs2BlockAlloc += 1;
			pSeqBlock->NumSeqIDs += 1;
			pSeqBlock->EndID += 1;
			}
		}

	AcquireLock(true);
	while(pSeqBlock->NxtID != 0)
		{
		pFlags = &m_Sequences.pSeqFlags[pSeqBlock->NxtID-1];
		PE1Flags = *pFlags;
		if(PE1Flags & cFlgSeqPE)
			{
			if((bTermPass = m_bTermPass) == false)
				m_Sequences.NumProcessed += 2;
			PE2Flags = pFlags[1];
			}
		else
			{
			if((bTermPass = m_bTermPass) == false)
				m_Sequences.NumProcessed += 1;
			PE2Flags = 0;
			}
	
		if(bTermPass)
			break;

		if(!((PE1Flags | PE2Flags) & (cFlgAsmbSeed | cFlgAsmbExtn | cFlgAsmbCplt)))	// if not already claimed or a seed, then return this sequence
			{
			PE1Flags |= cFlgAsmbSeed;				// accepting this (PE1 or a SE) as being a seed
			*pFlags++ = PE1Flags;
			*pPE1SeqID = pSeqBlock->NxtID;
			*pPE1SeqFlags = PE1Flags;
			Rslt = 1;
			pSeqBlock->NxtID += 1;
			if(PE1Flags & cFlgSeqPE)
				{
				PE2Flags |= cFlgAsmbSeed;			// accepting this PE2 as being a seed	
				*pFlags = PE2Flags;
				*pPE2SeqID = pSeqBlock->NxtID;
				*pPE2SeqFlags = PE2Flags;
				Rslt = 2;
				pSeqBlock->NxtID += 1;
				}
			if(pSeqBlock->NxtID > pSeqBlock->EndID) // last in block?
				pSeqBlock->NxtID = 0;
			break;
			}
		pSeqBlock->NxtID += 1;
		if(PE1Flags & cFlgSeqPE)
			pSeqBlock->NxtID += 1;
		if(pSeqBlock->NxtID > pSeqBlock->EndID)
			pSeqBlock->NxtID = 0;
		}
	ReleaseLock(true);
	}

ReleaseSerialiseSeqFlags();
return(Rslt);
}


// BuildOverlapExtensions
// build overlap extended sequences
int
CdeNovoAssemb::BuildOverlapExtensions(int CurPass,	// current overlap extension iteration (1..N)
						bool bAntisenseOnly,		// true if to process with antisense only probes
						bool bAllowSE2PE,			// if true then if SE overlaps PE1 or PE2 ends singularly but not both then merge the SE with the relevant overlapped PE end and retain as a PE
						int SubsKbp,				// allow this many induced substitutions per 1Kbp overlapping sequence fragments (defaults to 0, range 1..20)
						int MaxEnd12Subs,			// allow at the initial 12bp (2 x hexamer primer artefacts) of the 5' or 3' of overlap to have this many base mismatches in addition to the overall allowed MaxSubs		
						int CurMinPE1Len,			// currently the minimum PE1 length for any sequence of a paired end
						int CurMinPE2Len,			// currently the minimum PE2 length for any sequence of a paired end
						int CurMinSELen,			// currently the minimum length for any SE sequence
						int MinReqPEPrimOverlap,	// if primary probe is overlapping onto a PE then the initial primary overlap (onto PE1 or PE2) must be of at least this length
						int MinReqPESecOverlap,		// if primary probe was overlapping onto a PE then the secondary probe overlap (onto PE1 or PE2) must be of at least this length
						int MinReqPESumOverlap,		// if primary probe was overlapping onto a PE then the sum of the PE1 and PE2 overlap must be of at least this length
						int MinReqSEPrimOverlap,	// if primary probe is overlapping onto a SE then the overlap must be of at least this length
						int MinPEMergeOverlap,		// if probe PE1 and probe PE2 being considered for merging then there must be an overlap of at least this many bases
						int MinPE2SEOverlapLen,		// in final few passes can attempt to merge PEs into a SE with far smaller overlap of PE1 onto PE2 required
						int MinReqPESepDist,		// PE start sites expected to be separated by at least this many bases
						int MaxReqPESepDist,		// PE start sites expected be separated by no more than this many bases		
						int TrimPE2SE)				// trim PEs both 5' and 3' ends by this many bases before treating as SE
{
tsThreadOverlapExtendPars *pThreadParams;
tsThreadOverlapExtendPars *pCurThread;
int ThreadIdx;
int NumThreads;
tSeqID CurStartSeqID;
UINT32 CurSeqAllocIdx;

double OverlapRate[3];

UINT32 CurNumProcessed = 0;
UINT32 PrevNumProcessed = 0;
UINT32 PrevNumOverlapped = 0;
UINT32 NumOverlapped = 0;

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Starting overlap sequence extensions ...");

NumThreads = m_NumThreads;

// balance the number threads vs the number of sequences to minimise the thread startup costs
if(m_Sequences.NumSeqs2Assemb < 5)
	NumThreads = 1;
else
	{
	NumThreads = (m_Sequences.NumSeqs2Assemb + 4) / 5;
	if(NumThreads > m_NumThreads)
		NumThreads = m_NumThreads;
	}

// now that number of threads to use is known then initialise the sequence dispenser (GetSeqProc)
InitGetSeqProc(NumThreads);

if((pThreadParams = new tsThreadOverlapExtendPars[NumThreads])==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate memory for threads...");
	Reset(false);
	return(eBSFerrMem);
	}
memset(pThreadParams,0,sizeof(tsThreadOverlapExtendPars) * NumThreads);

m_AllocdThreadSeqsSize = (m_Sequences.bPESeqs ?  5 : 3) * NumThreads * (cMaxOvrlapSeqWrds + 128) * sizeof(tSeqWrd4); // allocating for 5 (3 if SE) blocks per thread plus a 128 word separators between individual thread buffers
if((m_pAllocdThreadSeqs = (tSeqWrd4 *)malloc(m_AllocdThreadSeqsSize)) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate memory for thread sequence buffering..."); 
	Reset(false);
	return(eBSFerrMem);
	}

m_bTermPass = false;
CurStartSeqID = 1;
m_Sequences.NumProcessed = 0;
m_Sequences.NumDuplicates = 0;
m_Sequences.NumOverlapping = 0;
m_Sequences.NumOverlapped = 0;

CurNumProcessed = 0;
PrevNumProcessed = 0;
NumOverlapped = 0;
PrevNumOverlapped = 0;
CurSeqAllocIdx = 0;

m_FinalProcSeqID = m_Sequences.NumSeqs2Assemb;
m_NextProcSeqID = 1;
m_StartProcSeqID = 1;
m_NumProcSeqIDs = cMaxMultiSeqFlags;
m_ThreadsProcessing = NumThreads;
pCurThread = pThreadParams;
#ifndef _WIN32
	// set the default stack 
	size_t defaultStackSize;
	struct timespec ts;
	int JoinRslt;
	pthread_attr_t threadattr; 
	pthread_attr_init(&threadattr);
	pthread_attr_getstacksize(&threadattr, &defaultStackSize);
	if(defaultStackSize != cWorkThreadStackSize)
		pthread_attr_setstacksize(&threadattr, cWorkThreadStackSize);
#endif
for(ThreadIdx = 1; ThreadIdx <= NumThreads; ThreadIdx++,pCurThread++)
	{
	pCurThread->ThreadIdx = ThreadIdx;
	pCurThread->bPESeqs = m_Sequences.bPESeqs ? true : false;
	pCurThread->pThis = this;

	pCurThread->bAntisenseOnly = bAntisenseOnly;	// true if to process with antisense only probes
	pCurThread->CurPass = CurPass;					// current overlap extension iteration (1..N)
	pCurThread->MaxSubs1K = SubsKbp;				// allowing for sub max per 1K overlapping bases
	pCurThread->MaxEnd12Subs = MaxEnd12Subs;		// allow at the initial 12bp (2 x hexamer primer artefacts) of the 5' or 3' of overlap to have this many base mismatches in addition to the overall allowed MaxSubs

	pCurThread->bAllowSE2PE = bAllowSE2PE;			// allowing SE to merge with a PE end resulting in a merged PE

	pCurThread->MinReqMergeLen = m_MinReqMergeLen;   // PE's must be at least this length before they can be merged into a SE
	pCurThread->CurMinPE1Len = CurMinPE1Len;		// currently the minimum PE1 length for any sequence of a paired end
	pCurThread->CurMinPE2Len = CurMinPE2Len;	    // currently the minimum PE2 length for any sequence of a paired end
	pCurThread->CurMinSELen= CurMinSELen;		    // currently the minimum length for any SE sequence

	pCurThread->MinReqPEPrimOverlap = MinReqPEPrimOverlap;		// if primary probe is overlapping onto a PE then the initial primary overlap (onto PE1 or PE2) must be of at least this length
	pCurThread->MinReqPESecOverlap = MinReqPESecOverlap;		// if primary probe was overlapping onto a PE then the secondary probe overlap (onto PE1 or PE2) must be of at least this length
	pCurThread->MinReqPESumOverlap = MinReqPESumOverlap;		// if primary probe was overlapping onto a PE then the sum of the PE1 and PE2 overlap must be of at least this length
	pCurThread->MinReqSEPrimOverlap = MinReqSEPrimOverlap;		// if primary probe is overlapping onto a SE then the overlap must be of at least this length
	pCurThread->MinPEMergeOverlap = MinPEMergeOverlap;			// if probe PE1 and probe PE2 being considered for merging then there must be an overlap of at least this many bases
	pCurThread->TrimPE2SE = TrimPE2SE;							// trim PEs both 5' and 3' ends by this many bases before treating as SE
	pCurThread->ProbeMinReqOverlap = pCurThread->MinReqPEPrimOverlap;
	if(pCurThread->ProbeMinReqOverlap == 0 || pCurThread->ProbeMinReqOverlap > pCurThread->MinReqSEPrimOverlap)
		pCurThread->ProbeMinReqOverlap = pCurThread->MinReqSEPrimOverlap;

	if(MinPE2SEOverlapLen > 16)
		pCurThread->MinPE2SEOverlapLen = MinPE2SEOverlapLen;	// in final few passes can attempt to merge PEs into a SE with far smaller overlap of PE1 onto PE2 required
	else
		pCurThread->MinPE2SEOverlapLen = MinPEMergeOverlap;

	pCurThread->MinReqPESepDist = MinReqPESepDist;			// PE start sites expected to be separated by at least this many bases
	pCurThread->MaxReqPESepDist = MaxReqPESepDist;			// PE start sites expected be separated by no more than this many bases

	pCurThread->MaxPE1SeqWrds = cMaxOvrlapSeqWrds;
	pCurThread->AllocPE1Seq = cMaxOvrlapSeqWrds * sizeof(tSeqWrd4);
	pCurThread->pPE1Seq = &m_pAllocdThreadSeqs[CurSeqAllocIdx];
	CurSeqAllocIdx += cMaxOvrlapSeqWrds + 128;

	pCurThread->MaxOverlapSeqWrds = cMaxOvrlapSeqWrds;
	pCurThread->AllocOverlapSeq = cMaxOvrlapSeqWrds * sizeof(tSeqWrd4);
	pCurThread->pOverlapSeq = &m_pAllocdThreadSeqs[CurSeqAllocIdx];
	CurSeqAllocIdx += cMaxOvrlapSeqWrds + 128;

	pCurThread->MaxTmpPE1SeqWrds = cMaxOvrlapSeqWrds;
	pCurThread->AllocTmpPE1Seq = cMaxOvrlapSeqWrds * sizeof(tSeqWrd4);
	pCurThread->pTmpPE1Seq = &m_pAllocdThreadSeqs[CurSeqAllocIdx];
	CurSeqAllocIdx += cMaxOvrlapSeqWrds + 128;

	if(m_Sequences.bPESeqs)
		{
		pCurThread->MaxPE2SeqWrds = cMaxOvrlapSeqWrds;
		pCurThread->AllocPE2Seq = cMaxOvrlapSeqWrds * sizeof(tSeqWrd4);
		pCurThread->pPE2Seq = &m_pAllocdThreadSeqs[CurSeqAllocIdx];
		CurSeqAllocIdx += cMaxOvrlapSeqWrds + 128;
		pCurThread->MaxTmpPE2SeqWrds = cMaxOvrlapSeqWrds;
		pCurThread->AllocTmpPE2Seq = cMaxOvrlapSeqWrds * sizeof(tSeqWrd4);
		pCurThread->pTmpPE2Seq = &m_pAllocdThreadSeqs[CurSeqAllocIdx];
		CurSeqAllocIdx += cMaxOvrlapSeqWrds + 128;
		}
	else
		{
		pCurThread->AllocPE2Seq = 0;
		pCurThread->pPE2Seq = NULL;
		pCurThread->AllocTmpPE2Seq = 0;
		pCurThread->pTmpPE2Seq = NULL;
		}

#ifdef _WIN32
	pCurThread->threadHandle = (HANDLE)_beginthreadex(NULL,cWorkThreadStackSize,ThreadedOverlapExtend,pCurThread,0,&pCurThread->threadID);
#else
	pCurThread->threadRslt = pthread_create (&pCurThread->threadID , &threadattr , ThreadedOverlapExtend , pCurThread);
#endif
	}
#ifndef _WIN32
pthread_attr_destroy(&threadattr);		// no longer required
#endif

OverlapRate[0] = m_EarlyOverlapTermThres * 10.0;		// will terminate this pass when overlap rate per minute has dropped below this threshold per 1M per minute for at least 3 minutes
OverlapRate[1] = OverlapRate[0];		// so initialise much higher than threshold to ensure at least 3 minutes of processing per pass
OverlapRate[2] = OverlapRate[0];

// allow threads a few seconds to startup
#ifdef _WIN32
	Sleep(2000);
#else
	sleep(2);
#endif
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Progress: 0 sequences processed");

int ReportMins = 0;
bool bErrTerm = false;
int ErrRet = 0;
pCurThread = pThreadParams;
for(ThreadIdx = 0; ThreadIdx < NumThreads; ThreadIdx++,pCurThread++)
	{
#ifdef _WIN32
	while(WAIT_TIMEOUT == WaitForSingleObject(pCurThread->threadHandle, 60000))
#else
	struct timespec ts;
	int JoinRlt;
	clock_gettime(CLOCK_REALTIME, &ts);
	ts.tv_sec += 60;
	while((JoinRlt = pthread_timedjoin_np(pCurThread->threadID, NULL, &ts)) != 0)
#endif
		{
		AcquireLock(false);
		CurNumProcessed = m_Sequences.NumProcessed;
		NumOverlapped = m_NumPartialSeqs2Assemb; 
		ReleaseLock(false);
		if(CurNumProcessed >= PrevNumProcessed)
			{
			if(ReportMins >= 5)		// reporting every 5 mins so as to not end up with excessively long log files - user really just needs assurance that the assembly is being progressed 
				{
				gDiagnostics.DiagOut(eDLInfo,gszProcName,"Progress: %u sequences processed with %u merges",CurNumProcessed,NumOverlapped);
				ReportMins = 0;
				}
			ReportMins += 1;
		
			// when thresholds have just been changed will not be any check for early termination of pass as m_EarlyOverlapTermThres will be set to be 0.0 
			if(m_EarlyOverlapTermThres > 0.0)
				{
				OverlapRate[0] = OverlapRate[1];
				OverlapRate[1] = OverlapRate[2];
				OverlapRate[2] = ((1000000.0 * (NumOverlapped - PrevNumOverlapped)) / (double)(CurNumProcessed - PrevNumProcessed));
				PrevNumProcessed = CurNumProcessed;
				PrevNumOverlapped = NumOverlapped;

				if(OverlapRate[0] < m_EarlyOverlapTermThres &&
					OverlapRate[1] < m_EarlyOverlapTermThres &&
					OverlapRate[2] < m_EarlyOverlapTermThres)
					m_bTermPass = true;
				}
			}

#ifdef _WIN32
		}
	CloseHandle( pCurThread->threadHandle);
#else
		ts.tv_sec += 60;
		}
#endif
	if(pCurThread->Rslt < 0)
		{
		bErrTerm = true;
		ErrRet = pCurThread->Rslt;
		}
	}

if(m_pAllocdThreadSeqs != NULL)
	{
#ifdef _WIN32
	free(m_pAllocdThreadSeqs);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pAllocdThreadSeqs != MAP_FAILED)
		munmap(m_pAllocdThreadSeqs,m_AllocdThreadSeqsSize);
#endif	
	m_pAllocdThreadSeqs = NULL;
	m_AllocdThreadSeqsSize = 0;
	}

if(pThreadParams)
	delete pThreadParams;

if(bErrTerm)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"One or more threads early terminated due to unrecoverable errors");
	return(ErrRet);
	}
if(m_bTermPass)
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Early terminated: merge rate less than %1.1f per million sequences over 3 minutes",m_EarlyOverlapTermThres);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Completed: %u sequences processed with %u merges",m_Sequences.NumSeqs2Assemb,m_NumPartialSeqs2Assemb);

return(eBSFSuccess);
}

const int cMaxMergeInDelLen = 3;			// microInDels of at most this length processed for sequence merging 
const int cMaxMergeMMLen = 3;			    // mismatches of at most this number processed for sequence merging
const int cMaxMergeWithInDelSeqLen = 1000;	// only attempt to merge sequences containing microInDels of upto cMaxMergeInDelLen for sequences no longer than this length

int					// if > 0 then length of merged sequence returned in pRetMerged
CdeNovoAssemb::MergeWithInDel(int MaxInDelLen,		// can accept an InDel of at most this length (clamped to be no more than cMaxMergeInDelLen)
			int MaxMMs,				// can accept overlaps with no more than this number of mismatches (clamped to be no more than 4)
			int OverlapLen,			// must be at least this many bases overlapping between Seq1 and Seq2 excluding any microInDel with overlap containing at most MaxMMs 
			int Seq1Len,			// number of packed bases in sequence1 
			tSeqWrd4 *pSeq1,        // pts to sequence1
			int Seq2Len,			// number of packed bases sequence2
			tSeqWrd4 *pSeq2,		// pts to sequence2
			etSeqBase *pMerged)     // if sequence1 and sequence2 can be merged then copy merged sequence into this buffer
{
UINT8 Seq1Bases[cMaxMergeWithInDelSeqLen+1];				// seq1 is unpacked into this buffer
UINT8 Seq2Bases[cMaxMergeWithInDelSeqLen+1];				// seq2 is unpacked into this buffer
UINT8 *pSeq1Base;
UINT8 *pSeq2Base;
if(Seq1Len > cMaxMergeWithInDelSeqLen || Seq2Len > cMaxMergeWithInDelSeqLen)
	return(-1);
if(MaxInDelLen < 0 || MaxInDelLen > cMaxMergeInDelLen)
	MaxInDelLen = cMaxMergeInDelLen;
if(MaxMMs < 0 || MaxMMs > cMaxMergeMMLen)
	MaxMMs = cMaxMergeMMLen;
GetSeq(pSeq1,Seq1Bases,Seq1Len);
Seq1Bases[Seq1Len] = eBaseEOS;

GetSeq(pSeq2,Seq2Bases,Seq2Len);
Seq2Bases[Seq2Len] = eBaseEOS;

// try for Seq1 5' overlapping onto Seq2 3'; if that fails then try for Seq2 5' overlapping onto Seq1 3' 
pSeq1Base = Seq1Bases;
etSeqBase Base1;
etSeqBase Base2;
int Seq1Idx3;
int Seq2Idx3;
int Seq1Idx5;
int Seq2Idx5;

int Matching5;
int Matching3;
int Mismatches5;
int Mismatches3;

int PutOverlapLen;

// attempt to locate a 5' overlap
for(Seq1Idx5 = 0; Seq1Idx5 < Seq1Len - 100; Seq1Idx5++)
	{
	pSeq1Base = &Seq1Bases[Seq1Idx5];
	pSeq2Base = Seq2Bases; 
	Matching5 = 0;
	Matching3 = 0;
	Mismatches5 = 0;
	Mismatches3 = 0;
	for(Seq2Idx5 = 0; Seq2Idx5 < min(Seq2Len,Seq1Len-Seq1Idx5); Seq2Idx5++,pSeq1Base++,pSeq2Base++)
		{
		Base1 = *pSeq1Base;
		Base2 = *pSeq2Base;
		if(Base1 == Base2)
			{
			Matching5+=1;
			continue;
			}
		if(Matching5 < 20)
			break;
		Mismatches5+=1;
		if(Mismatches5 < MaxMMs)
			continue;
		   
		PutOverlapLen = min(Seq2Len,Seq1Len-Seq1Idx5);         
		Seq1Idx3 = Seq1Idx5 + PutOverlapLen;
		PutOverlapLen += MaxInDelLen;
		if(PutOverlapLen > Seq2Len)
			PutOverlapLen = Seq2Len;
		int Idx;
        for(Idx = 0; Idx <= MaxInDelLen; Idx++)
			{
			pSeq1Base = &Seq1Bases[Seq1Idx3-1];
			pSeq2Base = &Seq2Bases[Seq2Len-1]; 
			Matching3 = 0;
			for(Seq2Idx3 = 0; Seq2Idx3 < Seq2Len; Seq2Idx3++,pSeq1Base-=1,pSeq2Base-=1)
				{
				Base1 = *pSeq1Base;
				Base2 = *pSeq2Base;
				if(Base1 == Base2)
					{
					Matching3+=1;
					continue;
					}
				Mismatches3+=1;
				if(Mismatches3 > 3 || Matching3 < 20)
					break;
				}
			}
		}
	OverlapLen = PutOverlapLen;

	if((MaxInDelLen + Matching5 + Matching3) > OverlapLen )	// accepting overlap even though it contains a microInDel and/or a couple of mismatches
		{

		}		
	}

return(0);
}

// ProcOverlapExtend
// Identify single read/contig or paired end overlaps for extensions
int
CdeNovoAssemb::ProcOverlapExtend(tsThreadOverlapExtendPars *pPars)
{
int TooMAnyWarnings;

int NumProbes;

int NumAcceptedOvl;
int NumExtdProbes;

tSeqID PE1ProbeSeqID;
UINT32 PE1ProbeSeqFlags;
UINT32 PE1ProbeSeqLen;
tSeqWrd4 *pPE1ProbeStartSeqWrd;

tSeqID PE2ProbeSeqID;
UINT32 PE2ProbeSeqFlags;	
UINT32 PE2ProbeSeqLen;
tSeqWrd4 *pPE2ProbeStartSeqWrd;

int OverlapRslt;

UINT32 NumProcessed;

UINT32 NumOverlapped;
UINT32 NumContained;
UINT32 TotNumPEOverlaps;
UINT32 TotNumAccepted;

int DiagLevel = gDiagnostics.GetFileDiagLevel();

gDiagnostics.DiagOut(eDLDebug,gszProcName,"Thread %d startup for %s extensions...",pPars->ThreadIdx, 
										pPars->bPESeqs ? "paired end sequences" : "sequences");

int NumSeqsToProc;
int NumProblems = 0;
int NumThisProbe;
int NumProbesWithOverlap;

int UnderProbeMinReqOverlap;


TooMAnyWarnings = 0;
NumOverlapped = 0;
NumContained = 0;
TotNumPEOverlaps = 0;
TotNumAccepted = 0;
UnderProbeMinReqOverlap = 0;

time_t Started = time(0);

// iterate over sequences (could be either PE (PE1 and PE2 returned), or SE (PE1 only returned) to be processed for overlap extensions
NumProbes = 0;
PE1ProbeSeqID = 0;
NumProcessed = 0;
NumProbesWithOverlap = 0;

// iterating over all sequences via GetSeqProc() with each call returning either 1 (SE/contig) or 2 (5' and 3' paired ends) sequences
while((NumSeqsToProc = GetSeqProc(&PE1ProbeSeqID,&PE1ProbeSeqFlags,&PE2ProbeSeqID,&PE2ProbeSeqFlags,pPars->ThreadIdx)) > 0)	// get next sequence(s)
	{
	NumProbes += 1;
	NumThisProbe = 0;
	NumProcessed += NumSeqsToProc;   // not yet actually processed so maybe a little presumptuous :-)
	if(!(NumProbes % 5))
		{
		time_t Now = time(0);
		unsigned long ElapsedSecs = (unsigned long) (Now - Started);
		if(ElapsedSecs >= 30)
			{
			pPars->NumProcessed += NumProcessed;
			AcquireLock(true);
			m_Sequences.NumOverlapped += NumProbesWithOverlap;
			ReleaseLock(true);
			NumProcessed = 0;
			NumProbesWithOverlap = 0;
			NumProbes = 0;
			Started = Now;
			}
		}

	if(PE1ProbeSeqFlags & cFlgSeqPE2) // probe PE2 overlapping is not currently supported
		continue;

	// take local copy of PE1 and PE2 (if applicable) sequences as sequences may require to be inplace merged with any overlaid target sequences 
	if((pPE1ProbeStartSeqWrd = GetSeqHeader(PE1ProbeSeqID,NULL,NULL,&PE1ProbeSeqLen,false)) == NULL) // serious problem if can't get header, currently just writes to log
		{
		if((TooMAnyWarnings+=1) < 10)
			gDiagnostics.DiagOut(eDLWarn,gszProcName,"Thread %d Couldn't find sequence header for known sequence %d...",pPars->ThreadIdx,PE1ProbeSeqID);
		continue;
		}
	if(PE2ProbeSeqID)
		{
		if((pPE2ProbeStartSeqWrd = GetSeqHeader(PE2ProbeSeqID,NULL,NULL,&PE2ProbeSeqLen,false))==NULL) // serious problem if can't get header, currently just writes to log
			{
			if((TooMAnyWarnings+=1) < 10)
				gDiagnostics.DiagOut(eDLWarn,gszProcName,"Thread %d Couldn't find PE2 sequence header for SeqID %u...",pPars->ThreadIdx,PE2ProbeSeqID);
			continue;
			}
		}
	else
		{
		pPE2ProbeStartSeqWrd = NULL;
		PE2ProbeSeqLen = 0;
		}


	// don't bother with probes which are smaller than the minimum required overlap; with luck other longer seed sequences will overlap these
    // also if PE probes have been significantly extended but still won't overlap PE1 onto PE2 then don't let these extend further
	if(PE1ProbeSeqLen < (UINT32)pPars->ProbeMinReqOverlap ||		
		(PE2ProbeSeqID != 0 && 
			((PE2ProbeSeqLen < (UINT32)pPars->ProbeMinReqOverlap) || ((PE1ProbeSeqLen + PE2ProbeSeqLen) > cMaxPEExtndLen))))
		{
		AcquireSerialiseSeqFlags();
		m_Sequences.pSeqFlags[PE1ProbeSeqID-1] &= ~(cFlgAsmbSeed | cFlgAsmbExtn | cFlgAsmbCplt);
		if(PE2ProbeSeqID)
			m_Sequences.pSeqFlags[PE2ProbeSeqID-1] &= ~(cFlgAsmbSeed | cFlgAsmbExtn | cFlgAsmbCplt);
		ReleaseSerialiseSeqFlags();
		UnderProbeMinReqOverlap +=1;
		continue;
		}


	// probes accepted as being potential seeds so make the copy
	GetSeqWrdSubSeq(0,PE1ProbeSeqLen, pPE1ProbeStartSeqWrd, (tSeqWrd4 *)pPars->pPE1Seq);
	if(PE2ProbeSeqID)
		GetSeqWrdSubSeq(0,PE2ProbeSeqLen, pPE2ProbeStartSeqWrd, (tSeqWrd4 *)pPars->pPE2Seq);

	pPars->PE1ProbeSeqID = PE1ProbeSeqID;
	pPars->PE1ProbeSeqFlags = PE1ProbeSeqFlags;
	pPars->PE1ProbeSeqLen = PE1ProbeSeqLen;

	pPars->PE2ProbeSeqID = PE2ProbeSeqID;
	pPars->PE2ProbeSeqFlags = PE2ProbeSeqFlags;
	pPars->PE2ProbeSeqLen = PE2ProbeSeqLen;

	NumAcceptedOvl = 0;
	NumExtdProbes = 0;

		// if PEs then check if possible to have an overlap at least MinPEMergeOverlap between PE1 and PE2; if so then merge and treat as if SE
	if(PE2ProbeSeqID)
		{
		int ProbeTargOverlap;
		int MergeLen;
		int LeftFlankLen;
		int MinPEMergeOverlap;
		if(pPars->MinPE2SEOverlapLen >= 16)
			MinPEMergeOverlap = min(pPars->MinPE2SEOverlapLen,pPars->MinPEMergeOverlap);
		else
			MinPEMergeOverlap = 16;

		ProbeTargOverlap = GetOverlapAB(MinPEMergeOverlap,&LeftFlankLen,pPars->PE1ProbeSeqLen,(tSeqWrd4 *)pPars->pPE1Seq,pPars->PE2ProbeSeqLen,(tSeqWrd4 *)pPars->pPE2Seq,true,pPars->MaxSubs1K,pPars->MaxEnd12Subs);
		if(ProbeTargOverlap)
			MergeLen = pPars->PE1ProbeSeqLen + pPars->PE2ProbeSeqLen - ProbeTargOverlap;
		else
			MergeLen = 0;
		if(ProbeTargOverlap && MergeLen >= pPars->MinReqMergeLen && MergeLen >= min(pPars->PE1ProbeSeqLen,pPars->PE2ProbeSeqLen))
			{
			pPars->PE1ProbeSeqLen = SeqMerge(LeftFlankLen,pPars->PE2ProbeSeqLen,(tSeqWrd4 *)pPars->pPE2Seq,pPars->PE1ProbeSeqLen,(tSeqWrd4 *)pPars->pPE1Seq,(tSeqWrd4 *)pPars->pTmpPE1Seq,(char *)"Called from ProcOverlapExtend");
			if(pPars->TrimPE2SE)
				pPars->PE1ProbeSeqLen = NHSeqTrim(pPars->TrimPE2SE, pPars->TrimPE2SE, pPars->PE1ProbeSeqLen,(tSeqWrd4 *)pPars->pPE1Seq);
			pPars->PE2ProbeSeqID = 0;					// was merged into the probe PE1 sequence so probe PE2 no longer relevant
			pPars->PE2ProbeSeqLen = 0;
			pPars->PE2ProbeSeqFlags = 0;
			NumProbesWithOverlap += 1;
			NumAcceptedOvl += 1;
			}
		}

	// try to locate all overlaps sense/sense with maximal extensions
	pPars->bPESeqsRevCpl = false;
	do {
		NumExtdProbes = 0;						// incremented if any extensions of SE/PE1 or PE2

		if(!pPars->bAntisenseOnly)				// if not processing for antisense probes only  
			{
			if((OverlapRslt = MergeOverlaps(pPars)) > 0)
				{
				NumAcceptedOvl += 1;		// accepted at least one overlap and/or extension for these probes
				NumProbesWithOverlap += 1;
				if(OverlapRslt == 2)
					{		
					NumExtdProbes += 1;			// at least one probe extension
					// if PE probes have been significantly extended but still won't overlap PE1 onto PE2 then don't let these extend further
					if(pPars->PE2ProbeSeqLen != 0 && (pPars->PE1ProbeSeqLen + pPars->PE2ProbeSeqLen) > cMaxPEExtndLen)
						break;
					}
				}
			if(OverlapRslt == 3)			// major problem???
				break;
			}

		if(!m_bSenseStrandOnly)
			{
			// try merging with probe antisense sequences
			PackedRevCpl((tSeqWrd4 *)pPars->pPE1Seq);				// will always be a PE1 sequence
			if(pPars->PE2ProbeSeqLen)					// check as PE2 have been previously merged into PE1
				{
				PackedRevCpl((tSeqWrd4 *)pPars->pPE2Seq);

				void *pPETmp = pPars->pPE1Seq;
				pPars->pPE1Seq = pPars->pPE2Seq;
				pPars->pPE2Seq = pPETmp;

				PE2ProbeSeqID = pPars->PE1ProbeSeqID;
				pPars->PE1ProbeSeqID = pPars->PE2ProbeSeqID;
				pPars->PE2ProbeSeqID = PE2ProbeSeqID;

				PE2ProbeSeqFlags = pPars->PE1ProbeSeqFlags;
				pPars->PE1ProbeSeqFlags = pPars->PE2ProbeSeqFlags;
				pPars->PE2ProbeSeqFlags = PE2ProbeSeqFlags;

				PE2ProbeSeqLen = pPars->PE1ProbeSeqLen;
				pPars->PE1ProbeSeqLen = pPars->PE2ProbeSeqLen;
				pPars->PE2ProbeSeqLen = PE2ProbeSeqLen;
				}

			// try to locate all overlaps antisense/sense with maximal extensions
			pPars->bPESeqsRevCpl = true;
			if((OverlapRslt = MergeOverlaps(pPars)) > 0)
				{
				NumAcceptedOvl += 1;				// accepted at least one overlap  and/or extension for these probes
				NumProbesWithOverlap += 1;
				if(OverlapRslt == 2)
					{		
					NumExtdProbes += 1;			// at least one probe extension
					// if PE probes have been significantly extended but still won't overlap PE1 onto PE2 then don't let these extend further
					if(pPars->PE2ProbeSeqLen != 0 && (pPars->PE1ProbeSeqLen + pPars->PE2ProbeSeqLen) > cMaxPEExtndLen)
						break;
					}
				}
			if(OverlapRslt == 3)					// bail out if major problems
					break;	
			}


		// if probes not overlaying any target sequences then probe no longer used as an assembly seed
		// and can now be potentially overlaid by some other probe assembly seed ...
		if(NumAcceptedOvl == 0)
			{
			AcquireSerialiseSeqFlags();
			if(pPars->PE1ProbeSeqID)
				m_Sequences.pSeqFlags[pPars->PE1ProbeSeqID-1] &=  ~(cFlgAsmbSeed | cFlgAsmbExtn | cFlgAsmbCplt);
			if(pPars->PE2ProbeSeqID)
				m_Sequences.pSeqFlags[pPars->PE2ProbeSeqID-1] &= ~(cFlgAsmbSeed | cFlgAsmbExtn | cFlgAsmbCplt);
			ReleaseSerialiseSeqFlags();
			break;
			}

		if(!m_bSenseStrandOnly)
			{
			PackedRevCpl((tSeqWrd4 *)pPars->pPE1Seq);			// will always be a PE1 sequence
	
			if(pPars->PE2ProbeSeqID)				// need to check as PE2 may have been merged in with PE1
				{
				PackedRevCpl((tSeqWrd4 *)pPars->pPE2Seq);

				void *pPETmp = pPars->pPE1Seq;
				pPars->pPE1Seq = pPars->pPE2Seq;
				pPars->pPE2Seq = pPETmp;
				PE2ProbeSeqID = pPars->PE1ProbeSeqID;
				pPars->PE1ProbeSeqID = pPars->PE2ProbeSeqID;
				pPars->PE2ProbeSeqID = PE2ProbeSeqID;

				PE2ProbeSeqFlags = pPars->PE1ProbeSeqFlags;
				pPars->PE1ProbeSeqFlags = pPars->PE2ProbeSeqFlags;
				pPars->PE2ProbeSeqFlags = PE2ProbeSeqFlags;

				PE2ProbeSeqLen = pPars->PE1ProbeSeqLen;
				pPars->PE1ProbeSeqLen = pPars->PE2ProbeSeqLen;
				pPars->PE2ProbeSeqLen = PE2ProbeSeqLen;
				}
			pPars->bPESeqsRevCpl = false;
			}
		}
	while(NumExtdProbes > 0);

			// save probe sequences ready for next merge phase
	if(NumAcceptedOvl > 0)
		{
		INT64 PartialSeqsLen;
		PartialSeqsLen = SavePartialSeqs(pPars->PE1ProbeSeqLen,(tSeqWrd4 *)pPars->pPE1Seq,pPars->PE2ProbeSeqLen,pPars->PE2ProbeSeqLen == 0 ? NULL : (tSeqWrd4 *)pPars->pPE2Seq);
		if(PartialSeqsLen < (INT64)0)
			return((int)PartialSeqsLen);		// errors
		}
	}

pPars->NumProcessed += NumProcessed;
pPars->NumOverlapped += NumProbesWithOverlap;
AcquireLock(true);
m_Sequences.NumOverlapped += NumProbesWithOverlap;
ReleaseLock(true);
return(1);		// success
}

#pragma pack(1)
const int cSeqIDCacheEntries = 0x07fff;

typedef struct TAG_sPrevPutTarg {
	tSeqID TargID;		// putative target sequence identifer
	int ProbeStartWrdOfs;    // original overlap would have started at this probe ofs
} tsPrevPutTarg;
#pragma pack()


// Now for the hardwork - find and accept targets which are maximally overlaid by SE/PE1 and PE2
int				// returns 0: no merges, 1: merge but no extension, 2: merge with extension
CdeNovoAssemb::MergeOverlaps(tsThreadOverlapExtendPars *pPars)
{
UINT64 SfxWrdIdx;
tSeqWrd4 *pHit;
tSeqWrd4 *pCurProbeSeq;
tSeqWrd4 *pCurTargStartSeqWrd;
int CurProbeSeqLen;
UINT16 CurProbeSeqFlags;
tSeqID CurProbeSeqID;

tSeqID ProcSeqIDs[cSeqIDCacheEntries + 1];
tsPrevPutTarg PrevPutTags[cSeqIDCacheEntries + 1];

int MinOvrlp;
int MinOvrlpSubs;
int ActSubs;

int RelSfxWrdOfs;

int PrimOverlapLen;		
int SecOverlapLen;
int CurTargSeqLen;

tSeqID PE1TargSeqID;
UINT16 PE1TargSeqFlags;
int PE1TargSeqLen;
int PE1TargSeq5SubOfs;
void *pPE1TargStartSeqWrd;

tSeqID PE2TargSeqID;
UINT16 PE2TargSeqFlags;	
int PE2TargSeqLen;
int PE2TargSeq5SubOfs;
void *pPE2TargStartSeqWrd;

int ProbeStartWrdOfs;		// index into probe tSeqWrds of 1st word which could be overlapping onto the target

bool bProbeContainsTarg;
bool bProbeOverlapsTarg;

int LeftFlankLen;
int ProbeTargOverlap;

tSeqWrd4 CurTargSeqID;
UINT16 CurTargSeqFlags;

int CmpRslt;
int SubOfs;
bool bPE1accepted;
bool bAcceptedOvl = false;	  // will be set true only if probe overlay (or PE1 and PE2 if a paired end) accepted
bool bMergedExtn = false;	  // will be set true only if probe overlay (or PE1 and PE2 if a paired end) accepted resulted in probe extension
bool bPE1ContainsTarg = false;
bool bPE1OverlapsTarg = false;

int MergedLen;
int DeltaPE1PE2;
UINT32 NumAcceptedOvl;

CurProbeSeqLen = pPars->PE1ProbeSeqLen;
CurProbeSeqFlags = pPars->PE1ProbeSeqFlags;
CurProbeSeqID = pPars->PE1ProbeSeqID;
pCurProbeSeq = (tSeqWrd4 *)pPars->pPE1Seq;
bPE1accepted = false;
ActSubs = 0;
NumAcceptedOvl = 0;

int	NonPEOverlapped = 0;
int PEOverlapped = 0;
int NumMergedExtns = 0;
int NumThisProbe = 0;
GetNHSeqWrdSubSeq(0,CurProbeSeqLen, pCurProbeSeq, (tSeqWrd4 *)pPars->pOverlapSeq);

memset(ProcSeqIDs,0,sizeof(ProcSeqIDs));
memset(PrevPutTags,0,sizeof(PrevPutTags));

tSeqWrd4 PrefixSeq[1000];
int CurMinOvrlp;
for(SubOfs = 0; SubOfs <= (int)(CurProbeSeqLen - pPars->ProbeMinReqOverlap); SubOfs++) // if only the 1st SeqWrd has been indexed
	{
	if(pPars->MaxSubs1K || pPars->MaxEnd12Subs)		// if allowing subs then ensure that any seed overlap will be at least cMinErrSeedLen
		{
		if(CurProbeSeqLen - SubOfs < cMinErrSeedLen)
			break;
		}

	// if no subs allowed then use maximally sized minimum overlap otherwise cMinErrSeedLen as a seed so subs can be handled reasonably efficiently
	if(!(pPars->MaxSubs1K || pPars->MaxEnd12Subs))
		MinOvrlp = min((int)m_Sequences.MinSeqLen,(int)(CurProbeSeqLen - SubOfs));
	else
		{
		MinOvrlp = cMinErrSeedLen;									// subs allowed so need to use exactly matching SeqWrds as a seed; one SeqWrd contains 15 bases
		MinOvrlpSubs = min(pPars->ProbeMinReqOverlap,(int)(CurProbeSeqLen - SubOfs));
		}

	if(!(SubOfs % 90))
		CurMinOvrlp = GetNHSeqWrdSubSeq(SubOfs,min(MinOvrlp+150,CurProbeSeqLen - SubOfs), pCurProbeSeq, PrefixSeq);
	else
		CurMinOvrlp = ShfLeftPackedSeq(CurMinOvrlp,PrefixSeq);


	SfxWrdIdx =	LocateFirstExact(m_Sequences.SfxElSize,	// sizeof elements in pSfxArray - currently will be either 4 or 5 bytes
					PrefixSeq,							// pts to probes flank subsequence
					MinOvrlp,							// probe length (in bases, not tSeqWrd4's) required to minimally exactly match over
					m_Sequences.pSeqs2Assemb,			// target sequence
					(UINT8 *)m_Sequences.pSuffixArray,	// target sequence suffix array
					0,									// low index in pSfxArray
					m_Sequences.NumSuffixEls-1);		// high index in pSfxArray

	if(!SfxWrdIdx)		// will be 0 if unable to find any targets prefix (ProbeMinReqOverlap) sequences exactly matching that of the probes subsequence (pPars->pOverlapSeq)
		continue;
	GetNHSeqWrdSubSeq(SubOfs,CurProbeSeqLen - SubOfs, pCurProbeSeq, (tSeqWrd4 *)pPars->pOverlapSeq);

	// iterating over all overlapped target sequences
	do {
		CmpRslt = 0;
		ActSubs = 0;
		if((pHit = (tSeqWrd4 *)SfxIdxToFirstSeqWrd(SfxWrdIdx++,&RelSfxWrdOfs))==NULL)  // get ptr to targets starting (immediately following header) seqword containing 5' bases 
			break;									   // should only be a NULL if exhusted all potential overlaps

		// which target was hit and what length is it?
		pCurTargStartSeqWrd = (tSeqWrd4 *)GetSeqHeader(pHit,&CurTargSeqID,NULL,NULL,(UINT32 *)&CurTargSeqLen,false);  
		if(CurTargSeqID == pPars->PE1ProbeSeqID || CurTargSeqID == pPars->PE2ProbeSeqID)	// if self hit to probe sequence, or effectively a self hit from probes PE1 to PE2, then not interested, try for another target sequence
			continue;

		if(ProcSeqIDs[CurTargSeqID & cSeqIDCacheEntries] == CurTargSeqID)	    // has thread already processed and accepted this target sequence as being overlaid?
			continue;

		// has thread already tried this putative target?
		ProbeStartWrdOfs = (SubOfs - (RelSfxWrdOfs * 15))/15;
		if(PrevPutTags[CurTargSeqID & cSeqIDCacheEntries].TargID == CurTargSeqID &&   
			PrevPutTags[CurTargSeqID & cSeqIDCacheEntries].ProbeStartWrdOfs == ProbeStartWrdOfs)
			continue;
		PrevPutTags[CurTargSeqID & cSeqIDCacheEntries].TargID = CurTargSeqID;
		PrevPutTags[CurTargSeqID & cSeqIDCacheEntries].ProbeStartWrdOfs = ProbeStartWrdOfs;

		// if overlap is less than required minimum overlap then must have exhusted all potential overlaps for current probe
		if((PrimOverlapLen = GetExactMatchLen((tSeqWrd4 *)pPars->pOverlapSeq,&pHit[RelSfxWrdOfs],(pPars->MaxSubs1K || pPars->MaxEnd12Subs) ? MinOvrlp : 0)) < MinOvrlp)
			break;

		if((RelSfxWrdOfs * 15) > SubOfs)		// only interested in extending seed if sure that probe will overlap on to target 
			continue;

			// if allowing for substitutions then get full overlap including any subs ensuring it is of at least MinOvrlpSubs long
		if(pPars->MaxSubs1K || pPars->MaxEnd12Subs)
			{
			// as sequences are extended then it gets inefficent to compare complete sequences... 
			// perhaps should make use of seed and try extending left and right from seed?????
			// SubOfs (0..n) is base offset in pCurProbeSeq corresponding to pPars->pOverlapSeq at which the exact match starts
			// pCurTargStartSeqWrd pts to start of target header
			// pHit pts to 1st tSeqWrd immediately following target header
			// pHit[RelSfxWrdOfs] is the tSeqWrd which exactly matched pPars->pOverlapSeq

			int ReqMatchLen;
			int AcceptedLeftFlankLen;
			int AcceptedActSubs;
			bool bAccepted;
#ifndef USETHISCODE
			AcceptedLeftFlankLen = SubOfs - (RelSfxWrdOfs * 15);
			ReqMatchLen = min(CurTargSeqLen,CurProbeSeqLen - AcceptedLeftFlankLen);
			if(!(bAccepted = IsMatched(ReqMatchLen,AcceptedLeftFlankLen ,CurProbeSeqLen,pCurProbeSeq,0,CurTargSeqLen,pHit,pPars->MaxSubs1K,pPars->MaxEnd12Subs,&AcceptedActSubs)))
				continue;
			LeftFlankLen = AcceptedLeftFlankLen;
			PrimOverlapLen = ReqMatchLen;
			ActSubs = AcceptedActSubs;
#else
			int MaxProbeOvrlapLen = CurProbeSeqLen - (ProbeStartWrdOfs * 15);
			if(MaxProbeOvrlapLen > CurTargSeqLen + 16)
				MaxProbeOvrlapLen = CurTargSeqLen + 16;

			PrimOverlapLen = GetOverlapAB(MinOvrlpSubs,&LeftFlankLen,MaxProbeOvrlapLen,&pCurProbeSeq[ProbeStartWrdOfs],CurTargSeqLen,(tSeqWrd4 *)pHit,true,pPars->MaxSubs1K,pPars->MaxEnd12Subs,&ActSubs);
			if(PrimOverlapLen < MinOvrlpSubs)
				continue;
			LeftFlankLen += ProbeStartWrdOfs * 15;
#endif
			}

		// if a potential overlapped target is either another seed or already part of some extended fragment (could even be part of current fragment!)
		// then don't incorporate into current fragment extensions
		// not serialising read access to the flags with AcquireSerialiseSeqFlags() and  ReleaseSerialiseSeqFlags() there is a significant throughput hit with the serialisation and
		// later will be checking flags with serialisation before updating flags within AtomicSeqMerge()
		UINT16 *pFlag;
		pFlag = &m_Sequences.pSeqFlags[CurTargSeqID-1];
		CurTargSeqFlags = *pFlag;

		if(CurProbeSeqFlags & cFlgSeqPE && CurTargSeqFlags & cFlgSeqPE2)		// only accept probe SE overlaping onto PE2's, not PE1 overlapping onto PE2
			{
			ProcSeqIDs[CurTargSeqID & cSeqIDCacheEntries] = CurTargSeqID;					// mark so target sequence not reprocessed by thread for this probe 
			continue;
			}

		CurTargSeqFlags = *pFlag;
		if(CurTargSeqFlags & cFlgSeqPE2)										// will be true only if probe is SE
			{
			PE1TargSeqFlags = pFlag[-1];
			PE2TargSeqFlags = CurTargSeqFlags;
			}
		else
			{
			if(CurTargSeqFlags & cFlgSeqPE)
				PE2TargSeqFlags = pFlag[1];
			if(!(CurTargSeqFlags & cFlgSeqPE))
				PE2TargSeqFlags = 0;
			PE1TargSeqFlags = CurTargSeqFlags;
			}

			// slough this potential hit if target still being used as a seed or has already been claimed as extended by this or another thread
		if((PE1TargSeqFlags | PE2TargSeqFlags) & (cFlgAsmbSeed | cFlgAsmbExtn | cFlgAsmbCplt))
			{
			pPars->NumAlreadyClaimed += 1;
			ProcSeqIDs[CurTargSeqID & cSeqIDCacheEntries] = CurTargSeqID;		    // mark so target sequence not reprocessed by this thread for this probe
			continue;
			}
		
		// classify overlapped target as either being completely contained within probe or partially overlapped by probe
		bProbeContainsTarg = false;
		bProbeOverlapsTarg = false;
		if(PrimOverlapLen >= (int)CurTargSeqLen)
			bProbeContainsTarg = true;
		else
			if(PrimOverlapLen >= (int)CurProbeSeqLen - SubOfs)
				bProbeOverlapsTarg = true;

		// if probe neither fully overlaps or contains target then try for another target sequence 
		if(!(bProbeContainsTarg || bProbeOverlapsTarg))
			continue;

		ProcSeqIDs[CurTargSeqID & cSeqIDCacheEntries] = CurTargSeqID;		    // mark so target sequence not reprocessed by this thread for this probe 


		// provisionally accepting the overlap of probe SE/PE1 onto a target; determine min required secondary overlap 
		// which will be dependent on if probe is SE or PE and if target was SE or PE
		if(CurTargSeqFlags & cFlgSeqPE)						// if overlapped target is either the PE1 or PE2 of a paired end
			{
			if(PrimOverlapLen < pPars->MinReqPEPrimOverlap) // check for overlap meeting PE primary 
				continue;
			SecOverlapLen = max(pPars->MinReqPESumOverlap - PrimOverlapLen,pPars->MinReqPESecOverlap);
			}
		else												// else overlaid target was a SE
			{
			if(!(CurProbeSeqFlags & cFlgSeqPE))				// if probe was a SE overlaping another SE then need at least MinReqSEPrimOverlap
				{
				if(PrimOverlapLen < pPars->MinReqSEPrimOverlap) // check for overlap meeting SE primary
					continue;
				SecOverlapLen = 0;							// no secondary overlap will be required
				}
			else											// else a PE onto a target SE so need the PE's mate end to also be overlapping the same target SE 
				{
				if(PrimOverlapLen < pPars->MinReqPEPrimOverlap)
					continue;
				SecOverlapLen = max(pPars->MinReqPESumOverlap - PrimOverlapLen,pPars->MinReqPESecOverlap);
				}
			}


		// PE1 probe contains or overlaps onto a putative target PE1 or PE2
		// initalise target PE1/PE1Targ's SeqIDs, SeqFlags, SeqLen, StartSeqWrd accordingly
		if(!(CurTargSeqFlags & cFlgSeqPE2))		// if a PE1 (could be a SE) was the overlaid target
			{
			PE1TargSeqID = CurTargSeqID;
			PE1TargSeqLen = CurTargSeqLen;
			PE1TargSeq5SubOfs = SubOfs;
			pPE1TargStartSeqWrd = pCurTargStartSeqWrd;
			if(CurTargSeqFlags & cFlgSeqPE)      // if target PE1 was paired end (not single end) then also need to initialise target PE2 
				{
				PE2TargSeqID = PE1TargSeqID + 1;
				pPE2TargStartSeqWrd = GetSeqHeader(PE2TargSeqID,NULL,NULL,(UINT32 *)&PE2TargSeqLen,false);
				PE2TargSeq5SubOfs = 0;			// actual is determined in subsequent overlay processing
				}
			else    // else was single ended or a contig so no putative overlap onto target PE2 to be processed 
				{
				PE2TargSeqID = 0;	
				PE2TargSeqLen = 0;
				PE2TargSeq5SubOfs = 0;	
				pPE2TargStartSeqWrd = 0;
				}
			}
		else    // else a PE2 of a paired end was overlapped so need to process target PE1
			{
			PE2TargSeqID = CurTargSeqID;
			PE2TargSeqLen = CurTargSeqLen;
			pPE2TargStartSeqWrd = pCurTargStartSeqWrd;
			PE2TargSeq5SubOfs = SubOfs;
			PE1TargSeqID = PE2TargSeqID - 1;
			pPE1TargStartSeqWrd = GetSeqHeader(PE1TargSeqID,NULL,NULL,(UINT32 *)&PE1TargSeqLen,false);
			PE1TargSeq5SubOfs = 0;			// actual is determined in subsequent processing
			}

		// specific processing dependent on if probe is a SE/contig, or the PE1 or PE2 from a paired end
		bAcceptedOvl = false;    // will be set true only if probe overlay (or PE1 and PE2 if a paired end) has been accepted 
		bMergedExtn = false;	  // will be set true only if probe overlay (or PE1 and PE2 if a paired end) accepted resulted in probe extension
			
		if(!(CurProbeSeqFlags & cFlgSeqPE))		// if probe was single ended 
			{
			if(!(CurTargSeqFlags & cFlgSeqPE))	// can immediately accept if overlaying or containing another single ended providing not already claimed
				{
				// try the merge, if already claimed then will fail with MergedLen == 0
				MergedLen = AtomicSeqMerge(PE1TargSeq5SubOfs,PE1TargSeqLen,(tSeqWrd4 *)pPE1TargStartSeqWrd,
												pPars->PE1ProbeSeqLen,(tSeqWrd4 *)pPars->pPE1Seq,(tSeqWrd4 *)pPars->pTmpPE1Seq,CurTargSeqID);
				if(!MergedLen)					// 0 if already claimed
					{
					pPars->NumAlreadyClaimed += 1;
					continue;
					}
				if(pPars->PE1ProbeSeqLen < MergedLen)	// was there a change in length; if so then flag as being extended 
					{
					pPars->PE1ProbeSeqLen = MergedLen;
					bMergedExtn = true;
					}
				bAcceptedOvl = true;
				NumAcceptedOvl += 1;
				}
			else                // else a SE probe is overlaying either the PE1 or PE2 of a paired end
								// normally the PE's mate required to also be overlaid
				{
				if(!(CurTargSeqFlags & cFlgSeqPE2))  // if overlaying a PE1 then also need an overlap onto it's mate PE2 before accepting
					{
					// try for the probe SE overlap onto PE2; will require overlap to be in the range MinReqPESepDist .. MaxReqPESepDist downstream from the PE1 overlap
					ProbeTargOverlap = GetOverlapAB(SecOverlapLen,&LeftFlankLen,pPars->PE1ProbeSeqLen,(tSeqWrd4 *)pPars->pPE1Seq,PE2TargSeqLen,(tSeqWrd4 *)pPE2TargStartSeqWrd,true,pPars->MaxSubs1K,pPars->MaxEnd12Subs);
					if(ProbeTargOverlap < SecOverlapLen  || LeftFlankLen < PE1TargSeq5SubOfs) // need at least the secondary overlap and PE1 can't start 5' to PE2
						{
						if(!pPars->bAllowSE2PE)		 
							continue;

						// Experimental: TODO: if requested then merge the SE with PE1 to form an extended PE  ....
						int MergedPE1ProbeLen;
						int MergedPE2ProbeLen;

						PE2TargSeq5SubOfs = 0;

						memcpy(pPars->pPE2Seq,pPE2TargStartSeqWrd, 4 * ((PE2TargSeqLen + 14) / 15));
						pPars->PE2ProbeSeqLen = PE2TargSeqLen;

						MergedLen = AtomicSeqMerge(PE1TargSeq5SubOfs,(int)PE1TargSeqLen,(tSeqWrd4 *)pPE1TargStartSeqWrd,
											PE2TargSeq5SubOfs,(int)PE2TargSeqLen,(tSeqWrd4 *)pPE2TargStartSeqWrd,
											pPars->PE1ProbeSeqLen,&MergedPE1ProbeLen,(tSeqWrd4 *)pPars->pPE1Seq,(tSeqWrd4 *)pPars->pTmpPE1Seq,
											pPars->PE2ProbeSeqLen,&MergedPE2ProbeLen,(tSeqWrd4 *)pPars->pPE2Seq,(tSeqWrd4 *)pPars->pTmpPE2Seq,
											PE1TargSeqID,PE2TargSeqID);
						if(!MergedLen)		// 0 if targets already claimed
							{
							pPars->PE2ProbeSeqLen = 0;		// revert probe back to being SE
							continue;
							}

						if(MergedLen != (pPars->PE1ProbeSeqLen + pPars->PE2ProbeSeqLen))
							{
							pPars->PE1ProbeSeqLen = MergedPE1ProbeLen;
							pPars->PE2ProbeSeqLen = MergedPE2ProbeLen;
							bMergedExtn = true;
							}

						// always the possibilty that the merged PE1 and PE2 themselves can be merged into a single SE or contig sequence!
						// check if a putative merge would result in an overlap of at least MinPEMergeOverlap
						ProbeTargOverlap = GetOverlapAB(pPars->MinPEMergeOverlap,&LeftFlankLen,pPars->PE1ProbeSeqLen,(tSeqWrd4 *)pPars->pPE1Seq,pPars->PE2ProbeSeqLen,(tSeqWrd4 *)pPars->pPE2Seq,true,pPars->MaxSubs1K,pPars->MaxEnd12Subs);
						int MergeLen;
						MergeLen = abs(LeftFlankLen) + ProbeTargOverlap;

						if(ProbeTargOverlap && MergeLen >= pPars->MinReqMergeLen && MergeLen >= min(pPars->PE1ProbeSeqLen,pPars->PE2ProbeSeqLen))
							{
							pPars->PE1ProbeSeqLen = SeqMerge(LeftFlankLen,pPars->PE2ProbeSeqLen,(tSeqWrd4 *)pPars->pPE2Seq,pPars->PE1ProbeSeqLen,(tSeqWrd4 *)pPars->pPE1Seq,(tSeqWrd4 *)pPars->pTmpPE1Seq,(char *)"Called from MergeOverlaps A");
							if(pPars->TrimPE2SE)
								pPars->PE1ProbeSeqLen = NHSeqTrim(pPars->TrimPE2SE, pPars->TrimPE2SE, pPars->PE1ProbeSeqLen,(tSeqWrd4 *)pPars->pPE1Seq);

							pPars->PE2ProbeSeqID = 0;					// was merged into the probe PE1 sequence so probe PE2 no longer relevant
							pPars->PE2ProbeSeqLen = 0;
							pPars->PE2ProbeSeqFlags = 0;
							bMergedExtn = true;
							PEOverlapped += 1;
 							}
						else
							NonPEOverlapped += 1;

						bAcceptedOvl = true;
						NumAcceptedOvl += 1;
						return(2);
						}
					else
						{
						//  ensure that PE2 starts within the max allowed range relative to PE1 start
						DeltaPE1PE2 = LeftFlankLen - PE1TargSeq5SubOfs;
						if(DeltaPE1PE2 > pPars->MaxReqPESepDist)
							continue;					
						
						// accepting the probe SE overlay onto target PE1 and PE2, try to merge into the probe SE 
						PE2TargSeq5SubOfs = LeftFlankLen;
						MergedLen = AtomicSeqMerge(PE1TargSeq5SubOfs,PE1TargSeqLen,(tSeqWrd4 *)pPE1TargStartSeqWrd,
														PE2TargSeq5SubOfs,PE2TargSeqLen,(tSeqWrd4 *)pPE2TargStartSeqWrd,
														pPars->PE1ProbeSeqLen,(tSeqWrd4 *)pPars->pPE1Seq,(tSeqWrd4 *)pPars->pTmpPE1Seq,PE1TargSeqID,PE2TargSeqID);
						if(!MergedLen)
							continue;
						if(pPars->PE1ProbeSeqLen != MergedLen)
							{
							pPars->PE1ProbeSeqLen = MergedLen;
							bMergedExtn = true;
							}
						bAcceptedOvl = true;
						NumAcceptedOvl += 1;
						}
					}
				else						// else SE contig overlays or contains PE2
											// normally require the PE1 to also be contained or overlaid by the SE
					{
					ProbeTargOverlap = GetOverlapAB(SecOverlapLen,&LeftFlankLen,pPars->PE1ProbeSeqLen,(tSeqWrd4 *)pPars->pPE1Seq,PE1TargSeqLen,(tSeqWrd4 *)pPE1TargStartSeqWrd,true,pPars->MaxSubs1K,pPars->MaxEnd12Subs);
					if(ProbeTargOverlap < SecOverlapLen || LeftFlankLen > PE2TargSeq5SubOfs)
						{
						// Experimental: TODO: if requested then merge the SE with PE2 to form an extended PE ....
						if(!pPars->bAllowSE2PE)		 
							continue;
						int MergedPE1ProbeLen;
						int MergedPE2ProbeLen;

						PE1TargSeq5SubOfs = 0;

						memcpy(pPars->pPE2Seq,pPars->pPE1Seq, 4 * ((pPars->PE1ProbeSeqLen + 14) / 15)); 
						pPars->PE2ProbeSeqLen = pPars->PE1ProbeSeqLen;

						memcpy(pPars->pPE1Seq,pPE1TargStartSeqWrd, 4 * ((PE1TargSeqLen + 14) / 15));
						pPars->PE1ProbeSeqLen = PE1TargSeqLen;

						MergedLen = AtomicSeqMerge(PE1TargSeq5SubOfs,(int)PE1TargSeqLen,(tSeqWrd4 *)pPE1TargStartSeqWrd,
											PE2TargSeq5SubOfs,(int)PE2TargSeqLen,(tSeqWrd4 *)pPE2TargStartSeqWrd,
											pPars->PE1ProbeSeqLen,&MergedPE1ProbeLen,(tSeqWrd4 *)pPars->pPE1Seq,(tSeqWrd4 *)pPars->pTmpPE1Seq,
											pPars->PE2ProbeSeqLen,&MergedPE2ProbeLen,(tSeqWrd4 *)pPars->pPE2Seq,(tSeqWrd4 *)pPars->pTmpPE2Seq,
											PE1TargSeqID,PE2TargSeqID);
						if(!MergedLen)		// 0 if targets already claimed
							{
							memcpy(pPars->pPE1Seq,pPars->pPE2Seq, 4 * ((pPars->PE2ProbeSeqLen + 14) / 15));			// restore probe
							pPars->PE1ProbeSeqLen = pPars->PE2ProbeSeqLen;
							pPars->PE2ProbeSeqLen = 0;
							continue;
							}

						if(MergedLen != (pPars->PE1ProbeSeqLen + pPars->PE2ProbeSeqLen))
							{
							pPars->PE1ProbeSeqLen = MergedPE1ProbeLen;
							pPars->PE2ProbeSeqLen = MergedPE2ProbeLen;
							bMergedExtn = true;
							}

						// always the possibilty that the merged PE1 and PE2 themselves can be merged into a single SE or contig sequence!
						// check if a putative merge would result in an overlap of at least MinPEMergeOverlap
						ProbeTargOverlap = GetOverlapAB(pPars->MinPEMergeOverlap,&LeftFlankLen,pPars->PE1ProbeSeqLen,(tSeqWrd4 *)pPars->pPE1Seq,pPars->PE2ProbeSeqLen,(tSeqWrd4 *)pPars->pPE2Seq,true,pPars->MaxSubs1K,pPars->MaxEnd12Subs);
						int MergeLen;
						MergeLen = abs(LeftFlankLen) + ProbeTargOverlap;

						if(ProbeTargOverlap && MergeLen >= pPars->MinReqMergeLen && MergeLen >= min(pPars->PE1ProbeSeqLen,pPars->PE2ProbeSeqLen))
							{
							pPars->PE1ProbeSeqLen = SeqMerge(LeftFlankLen,pPars->PE2ProbeSeqLen,(tSeqWrd4 *)pPars->pPE2Seq,pPars->PE1ProbeSeqLen,(tSeqWrd4 *)pPars->pPE1Seq,(tSeqWrd4 *)pPars->pTmpPE1Seq,(char *)"Called from MergeOverlaps B");
							if(pPars->TrimPE2SE)
								pPars->PE1ProbeSeqLen = NHSeqTrim(pPars->TrimPE2SE, pPars->TrimPE2SE, pPars->PE1ProbeSeqLen,(tSeqWrd4 *)pPars->pPE1Seq);

							pPars->PE2ProbeSeqID = 0;					// was merged into the probe PE1 sequence so probe PE2 no longer relevant
							pPars->PE2ProbeSeqLen = 0;
							pPars->PE2ProbeSeqFlags = 0;
							bMergedExtn = true;
							PEOverlapped += 1;
 							}
						else
							NonPEOverlapped += 1;

						bAcceptedOvl = true;
						NumAcceptedOvl += 1;
						return(2);
						}

					// length of overlap is acceptable but check for PE1 starting 5' upstream of PE2 and within the required separation range
					if(LeftFlankLen < 0)		// will be negative if PE1Targ is overlaying onto probe
						DeltaPE1PE2 = abs(LeftFlankLen) + PE2TargSeq5SubOfs;
					else
						DeltaPE1PE2 = PE2TargSeq5SubOfs - LeftFlankLen;

					if(DeltaPE1PE2 > pPars->MaxReqPESepDist)
						continue;					

					PE1TargSeq5SubOfs = LeftFlankLen;
					MergedLen = AtomicSeqMerge(PE1TargSeq5SubOfs,PE1TargSeqLen,(tSeqWrd4 *)pPE1TargStartSeqWrd,
												PE2TargSeq5SubOfs,PE2TargSeqLen,(tSeqWrd4 *)pPE2TargStartSeqWrd,
												pPars->PE1ProbeSeqLen,(tSeqWrd4 *)pPars->pPE1Seq,(tSeqWrd4 *)pPars->pTmpPE1Seq,PE1TargSeqID,PE2TargSeqID);

					if(!MergedLen)
						continue;
					if(pPars->PE1ProbeSeqLen != MergedLen)
						{
						pPars->PE1ProbeSeqLen = MergedLen;
						bMergedExtn = true;
						}
					bAcceptedOvl = true;
					NumAcceptedOvl += 1;
					}
				}
			}
	
		// if bAcceptedOvl not already set then overlapping probe must be PE1 of a paired ended and putative overlapped target could be SE or another PE
		if(!bAcceptedOvl && !(PE1TargSeqFlags & cFlgSeqPE)) // if overlapping a SE sequence then mate PE2 also needs to be overlapping
			{
			// require probe PE2 to also be overlapping
			ProbeTargOverlap = GetOverlapAB(SecOverlapLen,&LeftFlankLen,pPars->PE2ProbeSeqLen,(tSeqWrd4 *)pPars->pPE2Seq,PE1TargSeqLen,(tSeqWrd4 *)pPE1TargStartSeqWrd,true,pPars->MaxSubs1K,pPars->MaxEnd12Subs);
			if(ProbeTargOverlap < SecOverlapLen || LeftFlankLen < PE1TargSeq5SubOfs) // if couldn't meet secondary overlap requirments then try for another overlap
				continue;

			// There is a possibility that the PE2 sequence has been 3' extended, but not merged with the PE1 sequence because of some artefact in either the PE1 or PE2 sequences.
			// This PE2 sequence 3' extension is such that now the PE2 sequence actually starts 5' to the PE1 sequence but the PE1 and PE2 can't be merged into a SE because of the artefact in either PE1 or PE2
		    // But contained within the PE1 and PE2 sequence is a suffiently long subsequence which exactly overlaps or contains a SE sequence
			// Need to explicitly check for PE2 starting 5' of the PE1 using the overlap onto the SE as the anchor!
			// there is an overlap by both PE1 and PE2 so check for PE1 overlaying 5' to PE2 within expected distance range
			DeltaPE1PE2 = PE1TargSeq5SubOfs - LeftFlankLen;
			if(DeltaPE1PE2 < (int)pPars->MinReqPESepDist ||
				DeltaPE1PE2 > (int)pPars->MaxReqPESepDist)
				continue;					

			PE2TargSeq5SubOfs = LeftFlankLen;
	
				// merge the target contig plus PE2 onto the PE1 probe
			MergedLen = AtomicSeqMerge(PE1TargSeq5SubOfs * -1,PE1TargSeqLen,(tSeqWrd4 *)pPE1TargStartSeqWrd,
												PE2TargSeq5SubOfs,pPars->PE2ProbeSeqLen,(tSeqWrd4 *)pPars->pPE2Seq,
												pPars->PE1ProbeSeqLen,(tSeqWrd4 *)pPars->pPE1Seq,(tSeqWrd4 *)pPars->pTmpPE1Seq,PE1TargSeqID,0);
			if(!MergedLen)
				continue;


			pPars->PE2ProbeSeqID = 0;					// was merged into the probe PE1 sequence so probe PE2 no longer relevant
			pPars->PE2ProbeSeqLen = 0;
			pPars->PE2ProbeSeqFlags = 0;
			if(pPars->PE1ProbeSeqLen != MergedLen)
				pPars->PE1ProbeSeqLen = MergedLen;
			bMergedExtn = true;
			bAcceptedOvl = true;
			NumAcceptedOvl += 1;
			}


		// if bAcceptedOvl not already set then overlapping probe must be paired ended and putative overlapped target must also be paired ended
		if(!bAcceptedOvl)			
			{
			if(!(CurTargSeqFlags & cFlgSeqPE2)) // probe PE1 was overlapping onto a target PE1 then need PE2 to be either overlapping onto the target PE2 or alternatively onto the target PE1
				{
					// first try the more likely probe PE2 onto target PE2
				ProbeTargOverlap = GetOverlapAB(SecOverlapLen,&LeftFlankLen,pPars->PE2ProbeSeqLen,(tSeqWrd4 *)pPars->pPE2Seq,PE2TargSeqLen,(tSeqWrd4 *)pPE2TargStartSeqWrd,true,pPars->MaxSubs1K,pPars->MaxEnd12Subs);
				if(ProbeTargOverlap >= SecOverlapLen)
					{
					// have a probe PE2 onto target PE2 so try to merge probe PE1 with target PE1 and probe PE2 with target PE2
					int MergedPE1ProbeLen;
					int MergedPE2ProbeLen;

					PE2TargSeq5SubOfs = LeftFlankLen;

					MergedLen = AtomicSeqMerge(PE1TargSeq5SubOfs,(int)PE1TargSeqLen,(tSeqWrd4 *)pPE1TargStartSeqWrd,
											PE2TargSeq5SubOfs,(int)PE2TargSeqLen,(tSeqWrd4 *)pPE2TargStartSeqWrd,
											pPars->PE1ProbeSeqLen,&MergedPE1ProbeLen,(tSeqWrd4 *)pPars->pPE1Seq,(tSeqWrd4 *)pPars->pTmpPE1Seq,
											pPars->PE2ProbeSeqLen,&MergedPE2ProbeLen,(tSeqWrd4 *)pPars->pPE2Seq,(tSeqWrd4 *)pPars->pTmpPE2Seq,
											PE1TargSeqID,PE2TargSeqID);
					if(!MergedLen)		// 0 if targets already claimed
						continue;

					if(MergedLen != (pPars->PE1ProbeSeqLen + pPars->PE2ProbeSeqLen))
						{
						pPars->PE1ProbeSeqLen = MergedPE1ProbeLen;
						pPars->PE2ProbeSeqLen = MergedPE2ProbeLen;
						bMergedExtn = true;
						}


					// always the possibilty that the merged PE1 and PE2 themselves can be merged into a single SE or contig sequence!
					// check if a putative merge would result in an overlap of at least MinPEMergeOverlap
					ProbeTargOverlap = GetOverlapAB(pPars->MinPEMergeOverlap,&LeftFlankLen,pPars->PE1ProbeSeqLen,(tSeqWrd4 *)pPars->pPE1Seq,pPars->PE2ProbeSeqLen,(tSeqWrd4 *)pPars->pPE2Seq,true,pPars->MaxSubs1K,pPars->MaxEnd12Subs);
					int MergeLen;
					if(ProbeTargOverlap)
						MergeLen = pPars->PE1ProbeSeqLen + pPars->PE2ProbeSeqLen - ProbeTargOverlap;
					else
						MergeLen = 0;
  
					if(ProbeTargOverlap && MergeLen >= pPars->MinReqMergeLen && MergeLen >= min(pPars->PE1ProbeSeqLen,pPars->PE2ProbeSeqLen))
						{
						pPars->PE1ProbeSeqLen = SeqMerge(LeftFlankLen,pPars->PE2ProbeSeqLen,(tSeqWrd4 *)pPars->pPE2Seq,pPars->PE1ProbeSeqLen,(tSeqWrd4 *)pPars->pPE1Seq,(tSeqWrd4 *)pPars->pTmpPE1Seq,(char *)"Called from MergeOverlaps C");
						if(pPars->TrimPE2SE)
								pPars->PE1ProbeSeqLen = NHSeqTrim(pPars->TrimPE2SE, pPars->TrimPE2SE, pPars->PE1ProbeSeqLen,(tSeqWrd4 *)pPars->pPE1Seq);

						pPars->PE2ProbeSeqID = 0;					// was merged into the probe PE1 sequence so probe PE2 no longer relevant
						pPars->PE2ProbeSeqLen = 0;
						pPars->PE2ProbeSeqFlags = 0;
						bMergedExtn = true;
						PEOverlapped += 1;
 						}
					else
						NonPEOverlapped += 1;

					bAcceptedOvl = true;
					NumAcceptedOvl += 1;
					}
				else
					continue;
				}
			else
				continue;
			}
		}
	while(!ActSubs && !bMergedExtn && SfxWrdIdx <= m_Sequences.NumSuffixEls);
	if(bMergedExtn)
		{
		NumMergedExtns += 1;

		// check extensions were not PEs which are now excessively etended
		if(pPars->PE2ProbeSeqID && ((pPars->PE1ProbeSeqLen + pPars->PE2ProbeSeqLen) > cMaxPEExtndLen))
			break;

		bMergedExtn = false;
		CurProbeSeqFlags = pPars->PE1ProbeSeqFlags;
		CurProbeSeqLen = pPars->PE1ProbeSeqLen;
		pCurProbeSeq = (tSeqWrd4 *)pPars->pPE1Seq;
		CurMinOvrlp = GetNHSeqWrdSubSeq(SubOfs,CurProbeSeqLen - SubOfs, pCurProbeSeq, (tSeqWrd4 *)pPars->pOverlapSeq);
		CurMinOvrlp = GetNHSeqWrdSubSeq(0,min(MinOvrlp+150,CurMinOvrlp), (tSeqWrd4 *)pPars->pOverlapSeq, PrefixSeq);
		continue;
		}
	}
if(!NumMergedExtns && !bMergedExtn && !NumAcceptedOvl)
	return(0);

return(NumMergedExtns > 0 ? 2 : 1);
}


