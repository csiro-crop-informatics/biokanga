// Copyright 2013 CSIRO  ( http://www.csiro.au/ ) m_bRawDedupe
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

#pragma once

const size_t cWorkThreadStackSize = (1024*1024*2);	// working threads (can be multiple) stack size
const int cMaxPEExtndLen = ((cMinPETotSeqLen2SE * 3)/2);	// only allow PE extended length (sum of PE1 and PE2 lengths) to grow to this limit if no overlaps of PE1 onto PE2
const int cMinErrSeedLen = 60;						// when allowing for substitutions then use this minimum seed length

typedef enum TAG_edeNovoPMode {
	eAMEAssemble,				// standard de Novo assemble
	eAMESAssemble,				// more stringent de Novo assemble
	eAMQAssemble				// quick assemble packed reads with low stringency
} etdeNovoPMode;

#pragma pack(1)

typedef struct TAG_sSeqBlock {
	UINT32 NumSeqIDs;			    // number of sequence identifiers in this block
	tSeqID StartID;					// starting sequence identifier in this block
	tSeqID NxtID;					// next sequence identifier to be processed from this block
	tSeqID EndID;					// ending sequence identifier in this block
	} tsSeqBlock;

typedef struct TAG_sThreadOverlapExtendPars {
	int ThreadIdx;					// index of this thread (1..m_NumThreads)
	void *pThis;					// will be initialised to pt to CKangadna instance

#ifdef _WIN32
	HANDLE threadHandle;			// handle as returned by _beginthreadex()
	unsigned int threadID;			// identifier as set by _beginthreadex()
#else
	int threadRslt;					// result as returned by pthread_create ()
	pthread_t threadID;				// identifier as set by pthread_create ()
#endif
	int Rslt;						// returned result code

	bool bAntisenseOnly;			// true if to process with antisense only probes
	int CurPass;					// current overlap extension iteration

	int MinReqMergeLen;				// only merge PE ends if would result in SE of at least this length

	int MinReqPEPrimOverlap;		// if primary probe is overlapping onto a PE then the initial primary overlap (onto PE1 or PE2) must be of at least this length
	int MinReqPESecOverlap;			// if primary probe was overlapping onto a PE then the secondary probe overlap (onto PE1 or PE2) must be of at least this length
	int MinReqPESumOverlap;			// if primary probe was overlapping onto a PE then the sum of the PE1 and PE2 overlap must be of at least this length
	int MinReqSEPrimOverlap;		// if primary probe is overlapping onto a SE then the overlap must be of at least this length
	int ProbeMinReqOverlap;			// look for overlaps of at least this length when locating all potential overlaps of a probe onto some other sequence 
	int MinPEMergeOverlap;			// if probe PE1 and probe PE2 being considered for merging then there must be an overlap of at least this many bases
	int TrimPE2SE;					// trim PEs both 5' and 3' ends by this many bases before treating as SE

	int MinPE2SEOverlapLen;			// can override normal MinPEMergeOverlap in last few passes when PE1 and PE2 sequence extensions have been established

	int MinReqPESepDist;			// PE start sites expected to be separated by at least this many bases (minimum insert size - PE2 len)
	int MaxReqPESepDist;			// PE start sites expected to be separated by no more than this many bases (maximum insert size - PE2 len)

	int MaxSubs1K;					// allow at most this many substitutions per 1K overlapping bases
	int MaxEnd12Subs;				// allow at the initial 12bp (2 x hexamer primer artefacts) of the 5' or 3' of overlap to have this many base mismatches in addition to the overall allowed MaxSubs 

	bool bAllowSE2PE;				// if true then if SE overlaps PE1 or PE2 ends singularly but not both then merge the SE with the relevant overlapped PE end and retain as a PE

	int CurMinPE1Len;				// currently the minimum PE1 length for any sequence of a paired end
	int CurMinPE2Len;			    // currently the minimum PE2 length for any sequence of a paired end
	int CurMinSELen;			    // currently the minimum length for any SE sequence		

	bool bPESeqs;					// sequences include paired end reads (could also include SE/contigs); if reset then only SE/contig sequences
	
	tSeqID PE1ProbeSeqID;			// probe sequence identifier for PE1
	UINT16 PE1ProbeSeqFlags;		// probe sequence flags for PE1
	int PE1ProbeSeqLen;				// initial probe sequence PE1 length
	UINT32 AllocPE1Seq;				// memory allocated to pPE1Seq
	UINT32 MaxPE1SeqWrds;			// pPE1Seq allocated to hold at most this many tSeqWrds
	void *pPE1Seq;					// to hold an extended packed PE1 sequence
	 
	tSeqID PE2ProbeSeqID;			// probe sequence identifier for PE2
	UINT16 PE2ProbeSeqFlags;		// probe sequence flags for PE2
	int PE2ProbeSeqLen;				// initial probe sequence PE2 length
	UINT32 AllocPE2Seq;				// memory allocated to pPE2Seq
	UINT32 MaxPE2SeqWrds;			// pPE2Seq is allocated to hold at most this many tSeqWrds
	void *pPE2Seq;					// to hold an extended packed PE2 sequence

	bool bPESeqsRevCpl;				// true if PE1 and PE2 sequences are currently revcpl'd
	
	UINT32 AllocOverlapSeq;			// memory allocated to pOverlapSeq
	UINT32 MaxOverlapSeqWrds;		// pOverlapSeq allocated to hold at most this many tSeqWrds
	void *pOverlapSeq;				// to hold sequence to be overlapped
	
	UINT32 AllocTmpPE1Seq;			// memory allocated to pTmpPE1Seq
	UINT32 MaxTmpPE1SeqWrds;		// pTmpPE1Seq allocated to hold at most this many tSeqWrds
	void *pTmpPE1Seq;				// used for temp holding sequences whilst processing

	UINT32 AllocTmpPE2Seq;			// memory allocated to pTmpPE2Seq
	UINT32 MaxTmpPE2SeqWrds;		// pTmpPE2Seq allocated to hold at most this many tSeqWrds
	void *pTmpPE2Seq;				// used for temp holding sequences whilst processing

	UINT32 NumProcessed;			// number processed
	UINT64 NumAlreadyClaimed;		// number of overlaps which failed because already claimed by some other thread
	UINT64 NumOverlapped;			// number of sequences determined as being overlapped
	UINT32 NumPE1Overlapping;		// number of PE1 sequences which overlapped other sequences
	UINT32 NumPE2Overlapping;		// number of PE2 sequences which overlapped other sequences
} tsThreadOverlapExtendPars;

#pragma pack()

class CdeNovoAssemb : public CKangadna
{
	int m_TrimInputEnds;				// trim input sequences, both 5' and 3' ends by this many bases
	int m_TrimPE2SE;					// trim PEs both 5' and 3' ends by this many bases before treating as SE

	bool m_bSenseStrandOnly;			// sequences from sense strand specific
	bool m_bSingleEnded;				// treat all sequences as being single ended even if loaded as paired ends
	bool m_bTermPass;					// true if current overlap processing pass is to be early terminated 
	double m_EarlyOverlapTermThres;		// terminate current overlap processing pass if overlap rate drops below this threshold for 3 minutes

	UINT32 m_NxtSeqs2BlockAlloc;		// allocate next thread sequence block starting with this sequence
	UINT32 m_SeqsPerThreadBlk;			// nominal number of sequences per thread processing block
	tsSeqBlock m_ThreadSeqBlocks[cMaxWorkerThreads];	// blocks of sequence identifiers to be processed by each thread

	size_t m_AllocdThreadSeqsSize;		// m_pAllocdThreadSeqs was allocated to hold this many bytes
	tSeqWrd4 *m_pAllocdThreadSeqs;		// allocated to hold all sequence buffering as required by threads, each thead buffer (pTmpPE1Seq etc) is sub blocked from this buffer

	bool m_bProcPE;						// true if assembly processing includes PE sequences, false if for SE or contigs only

	int	// returns 0: no merges, 1: merge but no extension, 2: merge with extension
		MergeOverlaps(tsThreadOverlapExtendPars *pPars);

	int	InitGetSeqProc(int NumThreads);	// number of threads which will be requesting, via GetSeqProc(), sequence identifiers for processing
	
	int				// 0 if all returned, 1 if PE1 only, 2 if both PE1 and PE2
			GetSeqProc(tSeqID *pPE1SeqID,	// returned SE or PE1 sequence identifier
					  UINT32 *pPE1SeqFlags, // SE or PE1 flags
					  tSeqID *pPE2SeqID,	// if PE then PE2 sequence identifier
					  UINT32 *pPE2SeqFlags, // if PE then PE2 flags
					  int ThreadIdx);		// uniquely identifies calling thread, 1..NumThreads

		int					// if > 0 then length of merged sequence returned in pRetMerged
			MergeWithInDel(int MaxInDelLen,		// can accept an InDel of at most this length (clamped to be no more than cMaxMergeInDelLen)
					int MaxMMs,				// can accept overlaps with no more than this number of mismatches (clamped to be no more than 4)
					int OverlapLen,			// must be at least this many bases overlapping between Seq1 and Seq2 excluding any microInDel with overlap containing at most MaxMMs 
					int Seq1Len,			// number of packed bases in sequence1 
					tSeqWrd4 *pSeq1,        // pts to sequence1
					int Seq2Len,			// number of packed bases sequence2
					tSeqWrd4 *pSeq2,		// pts to sequence2
					etSeqBase *pMerged);     // if sequence1 and sequence2 can be merged then copy merged sequence into this buffer

	int				// combine partial sequences back with unmerged sequences ready for next pass
				CombinePartialsWithUnmerged(bool bTrim15bp = false,					// if true then trim 15bp of the 5' ends to reduce the primer induced errors
									int MinPETotSeqLen2SE = cMinPETotSeqLen2SE,		// combine partial sequences if length PE1 plus length PE2 is at least this length - 0 if no combining of PE1 and PE2s
								    int MinPESeqLen2SE = cMinPESeqLen2SE,			// only combine if individual PE1 and PE2 lengths at least this length
									int TrimPE2SE = 0);								// trim overlength non-overlapping PEs both 5' and 3' ends by this many bases before treating as SE


public:
	CdeNovoAssemb(void);
	~CdeNovoAssemb(void);

	teBSFrsltCodes AssembReads(	etdeNovoPMode PMode,	  // processing mode, currently either eAMEAssemble (default), eAMESAssemble (stringent) or eAMQAssemble (quick)
								int TrimInputEnds,		  // trim input sequences, both 5' and 3' ends by this many bases
							    int MinInputSeqLen,		  // only accept for assembly sequences which are, after any trimming, of at least this length
								int TrimPE2SE,			  // trim PEs both 5' and 3' ends by this many bases before treating as SE
								bool bSenseStrandOnly,	  // sequences from sense strand specific
								bool bAllowSE2PE,		  // if true then if SE overlaps PE1 or PE2 ends singularly but not both then merge the SE with the relevant overlapped PE end and retain as a PE
								bool bSingleEnded,		  // treat all sequences as being single ended even if loaded as paired ends
								int MaxPasses,			  // limit number of de Novo assembly passes to this maximum (quick mode defaults to 10, exhaustive defaults to 100) set to 0 for no limit
								int OutPass2File,		  // output partially assembled starting from this pass - 0 if only final is to be written to file
								int NReduceThresSteps,	  // reduce overlap thresholds over this many steps (defaults: 3 quick, 5 standard, 8 stringent assemble)");
								int Subs1Kbp,	 		  // allow this many induced substitutions per Kbp overlapping sequence fragments
								int MaxEnd12Subs,		  // allow at the initial 12bp (2 x hexamer primer artefacts) of the 5' or 3' of overlap to have this many base mismatches in addition to the overall allowed MaxSubs 
								int InitSEOvlp,			  // initial minimal SE overlap required to merge SEs
								int FinSEOvlp,			  // final minimal SE overlap required to merge SEs
								int InitPEOvlp,			  // initial minimal PE overlap required to merge PEs
								int FinPEOvlp,			  // final minimal PE overlap required to merge PEs
								int MinPE2SEOvlp,		  // minimal overlap of PE1 onto PE2 required to merge as SE
								int PE2SESteps,			  // when less than this many steps remaining then treat PE1 and PE2 as individual SE sequences if excessive lengths (defaults to 2, set 0 to disable)");
								int OrientatePE,		  // PE end orientations 0: sense/antisense, 1: sense/sense, 2: antisense/sense, 3: antisense/antisense 
								char *pszPE1File,		  // optional input high confidence seed PE1 sequences file
								char *pszPE2File,		  // optional input high confidence seed PE2 sequences file
								char *pszSeedContigsFile, // optional input high confidence seed SE contigs file
								char *pszInArtReducfile,  // optional input preprocessed artefact reduced packed reads from this file
								char *pszAssembFragsFile);  // where to write assembled sequence fragments as contigs ("SE" appended)

	teBSFrsltCodes LoadSeqsOnly(bool bSenseStrandOnly,			// process sequences as strand specific
								bool bSingleEnded,				// treat all sequences as being single ended even if loaded as paired ends
								int OrientatePE,					// PE end orientations 0: sense/antisense, 1: sense/sense, 2: antisense/sense, 3: antisense/antisense 
								char *pszPE1File,				// input high confidence seed PE1 sequences file
								char *pszPE2File,				// input high confidence seed PE2 sequences file
								char *pszSeedContigsFile,		// optional input high confidence seed SE contigs file
								char *pszInArtReducfile = NULL);	// optional input preprocessed artefact reduced packed reads from this file


		// build overlap extended sequences	
	int BuildOverlapExtensions(int CurPass,			// current overlap extension iteration (1..N)
						bool bAntisenseOnly,		// true if to process with antisense only probes
					    bool bAllowSE2PE,			// if true then if SE overlaps PE1 or PE2 ends singularly but not both then merge the SE with the relevant overlapped PE end and retain as a PE
						int Subs100bp,				// allow this many induced substitutions per 100bp overlapping sequence fragments (defaults to 1, range 0..5)
						int MaxEnd12Subs,			// allow at the initial 12bp (2 x hexamer primer artefacts) of the 5' or 3' of overlap to have this many base mismatches in addition to the overall allowed MaxSubs 
						int CurMinPE1Len,			// currently the minimum PE1 length for any sequence of a paired end
						int CurMinPE2Len,			// currently the minimum PE2 length for any sequence of a paired end
						int CurMinSELen,			// currently the minimum length for any SE sequence
						int MinReqPEPrimOverlap,	// if primary probe is overlapping onto a PE then the initial primary overlap (onto PE1 or PE2) must be of at least this length
						int MinReqPESecOverlap,		// if primary probe was overlapping onto a PE then the secondary probe overlap (onto PE1 or PE2) must be of at least this length
						int MinReqPESumOverlap,		// if primary probe was overlapping onto a PE then the sum of the PE1 and PE2 overlap must be of at least this length
						int MinReqSEPrimOverlap,	// if primary probe is overlapping onto a SE then the overlap must be of at least this length
						int MinPEMergeOverlap,		// if probe PE1 and probe PE2 being considered for merging then there must be an overlap of at least this many bases
						int MinPE2SEOverlapLen,		// in final few passes can attempt to merge PEs into a SE with far smaller (< 16 to disable) overlap of PE1 onto PE2 required
						int MinReqPESepDist,		// PE start sites expected to be separated by at least this many bases
						int MaxReqPESepDist,		// PE start sites expected be separated by no more than this many bases
						int TrimPE2SE);				// trim PEs both 5' and 3' ends by this many bases before treating as SE


	int ProcOverlapExtend(tsThreadOverlapExtendPars *pPars);	// threaded processing for overlaps
};

