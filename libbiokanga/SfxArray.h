#pragma once
#include "./commdefs.h"


typedef enum eEMRMode {
	eEMRAddMatch = 0,	// add new match
	eEMRIsAlreadyMatched, // check if already matched
	eEMRNearlyExtend	// nearly extended
	} etEMRMode;

typedef enum eRPTMasking {
	eRPTHignor = 0,		// treat all bases as not being a repeat (ignore any cRptMskFlg)
	eRPTHmasked,		// treat any base cRptMskFlg as a repeat masked base 
	eRPTHunmasked		// treat any base without cRptMskFlg as a repeat masked base 
	} etRPTMasking;

#pragma pack(1)

// Probes
typedef struct TAG_sAMProbe {
	INT32 ProbID;							// uniquely identifies this probe
	INT32 TotPlusHits;						// total hits on '+' strand
	INT32 TotMinHits;						// total hits on '-' strand
	INT32 ProbeLen;							// probe length
	char szDescr[80];						// as parsed from fasta descriptor or bioseq file
	char Strand;							// assumed orginal sequence strand when parsed, usually assume as being the '+' strand
	char CurStrand;						    // which strand '+','-' is currently being processed, if '-' and Strand == '+' then Bases holds reverse complement
	INT64 SeqRefLoc;						// from where descriptor was parsed/read from fasta or bioseqfile
	int CacheSeqOfs;						// if >=0 then offset in segment cache
} sAMProbe;


const int cAllocPutCores = 100000;		// allocate putative core instances in this many sized chunks
typedef struct TAG_sPutCore {
	UINT32 ProbeOfs;		// probe offset at which hit starts 
	UINT32 HitLen;		// number of bases hit
	UINT32 TargOfs;		// target loci hit	
} tsPutCore;

// entries contain the following header by which file type and construction parameters can be identified
typedef struct TAG_sSfxHeader {
	UINT32 ChromID;						// identifies chromosome (1..n) within this suffix file
	UINT32 ElSize;						// size of each suffix array element (in bytes)
	UINT32 SeqLen;						// sequence length (also suffix array size in # of elements)
} tsSfxHeader;

#pragma pack()

typedef struct TAG_sNMerFreq {
	int Freq;			// Frequency of this NMer sequence
	etSeqBase *pSeq;	// pts to this seq of m_NMerLen length
} tsNMerFreq;

class CSfxNMerFreq {
friend class CSfxArray;
	bool m_bOverRep;		// true if processing for over represented instances
	int m_NMerLen;			// freqs for NMers of this length
	int m_MinFreqCutoff;	// minimum freq cutoff
	int m_MaxFreqCutoff;	// maximum freq cutoff
	int m_MaxInstances;		// maximum number of NMer instances
	int m_CurNumInstances;	// current number of NMer instances
	tsNMerFreq *m_pInstances; // pts to array of NMer instances
	etSeqBase *m_pSeq;		// holds all NMer sequences concatenated

protected:
	int Add(int Freq,etSeqBase *pEl1);

public:
	CSfxNMerFreq(void);
	CSfxNMerFreq(bool m_bOverRep,int MinFreqCutoff,int MaxFreqCutoff, int MaxInstances, int NMerLen);
	~CSfxNMerFreq(void);
	void Reset(void);
	int Init(bool m_bOverRep,int MinFreqCutoff,int MaxFreqCutoff, int MaxInstances, int NMerLen);
	int GetNumInstances(void);
	int GetNMerLen(void);
};


class CSfxArray : public CBioSeqFile
{
	tsSfxHeader m_CurSfxHeader;					// currently loaded suffix header
	UINT32 *m_pSfxArray;					// currently loaded suffix array
	UINT32 m_CurMaxSfxSize;					// number of bytes currently allocated to hold sfx array
	etSeqBase *m_pSfxSeq;						// currently loaded sequence
	UINT32 m_CurMaxSeqSize;				// how much memory is allocated to hold sequence
	UINT32 m_CurSfxSeqChromID;
	UINT32 m_CurSfxArrayChromID;

    bool m_ChkSelfMatches;
	UINT32 m_ProbeEntryID;
	UINT32 m_TargEntryID;

	bool m_bProbeRepeatSense;				// false if bases marked as masked to be treated as repeats, otherwise unmarked bases treated as repeats
	bool m_bTargRepeatSense;				// false if bases marked as masked to be treated as repeats, otherwise unmarked bases treated as repeats
	UINT32 m_ProbeRepeatPC;			// percentage (0..100) repeat masked bases allowed in probe sequence
	UINT32 m_TargRepeatPC;			// percentage (0..100) repeat masked bases allowed in suffix sequence

	UINT32 m_Modes;
	UINT32 m_CurMatchMode;
	UINT32 m_CurSampleSize;			// current window probe sample size
    UINT32 m_CurMinMatchLen;			// minimal required match length
	UINT32 m_CurMaxMatchLen;			// maximum required match length
    UINT32 m_CurMaxNumMatches;		// only look for this many matches

	UINT32 m_ProbeLen;				// total length of probe

	UINT32 m_MinExactSeedLen;			// minimum seed length for identity extension 
	UINT32 m_MaxMutsCnt;				// maximum number of mutations allowed when extending
	UINT32 m_MaxNsCnt;				// max number of indeterminate bases allowed when extending
	UINT32 m_MaxMutBlockLen;			// maximum number of mutations alllowed in any sequence block

	bool m_bAllowContained;					// true if contained sequences are to be reported
	bool m_bDontExtend;						// true if matches are not to be maximally extended
    bool m_bSeqOut;							// true if matching sequences to be output in results

	int m_NumPutCores;						// number of putative cores
	int m_MaxPutCores;						// maximum number of cores which can be held in m_pPutCores
	tsPutCore *m_pPutCores;					// pts to memory allocated to hold putative matching cores


	int (* m_pRsltsFileHandler)(int Mode,UINT32 ProbePsn,UINT32 TargPsn,UINT32 MatchLen,UINT32 IdentCnt,UINT32 FHMode);

	tsSfxHeader *LoadSfxHeader(UINT32 EntryID);
	UINT32 *LoadSfxArray(UINT32 EntryID);
	etSeqBase *LoadSfxSeq(UINT32 EntryID,bool bClearSeqMsk = false);
	bool IsAlreadyMatched(UINT32 ProbePsn,UINT32 TargPsn,UINT32 MatchLen);
	UINT32 ExtendLeft(etSeqBase *pProbe,	// sequence containing probe
		  UINT32 ProbePsn,		// start position in probe sequence at which to extend left 
			  etSeqBase *pTarg,			// target sequence
			  UINT32 TargPsn,	// start position in target sequence at which to extend left
			  UINT32 Max2Extend);

	UINT32								// #bases extended 
			ExtendRight(etSeqBase *pProbe,	// sequence containing probe
			  etSeqBase *pTarg,				// target sequence
			  UINT32 Max2Extend);		// extend right for at most this many bases

	bool AcceptMasked(etSeqBase *pTarg,etSeqBase *pProbe,int MatchLen); // determine if repeat masked within limits

	int CmpProbeTarg(etSeqBase *pEl1,etSeqBase *pEl2,int Len);

		int			// index in pSfxArray of first exactly matching probe or -1 if no match					
		LocateFirstExact(etSeqBase *pProbe,  // sequence containing probe
				  UINT32 ProbeLen,	// probe length
				  etSeqBase *pTarg,			// target sequence
				  UINT32 *pSfxArray,	// target sequence suffix array
				  UINT32 TargStart,   // position in pTarg (0..n) corresponding to start of suffix array
				  int SfxLo,				// low index in pSfxArray
				  int SfxHi,				// high index in pSfxArray
	 			  int MaxNumExacts = 0);		// only interested in at most this number of exact matches (0 if all exacts are of interest)
	
	int AddPutCore(UINT32 ProbeOfs,	// probe offset at which hit starts 
				UINT32 HitLen,		// number of bases in core
				UINT32 TargOfs);		// target loci hit	


	int FiltPutCores(bool bRemoveOverlaps = true);		// filter putative cores - remove contained cores and optionally any overlap cores
	int FiltContainedPutCores(int ProbeIdx, int TargIdx); // remove contained cores
	int FiltOverlapPutCores(int ProbeIdx, int TargIdx);	// remove overlap cores

	static int SortPutCores(const void *arg1, const void *arg2);

public:
	CSfxArray(void);
	~CSfxArray(void);
	int Reset(bool bFlush = true);			// reset state back to that immediately following instantiation
	int Open(char *pszSeqFile,				// specifies file to open or create
			   bool bCreate = false);		// create file if it doesn't already exist, truncate if it does
	int Close(bool bFlush = true);			// closes opened file
	int	AddEntry(char *pszSource,		    // from where was this sequence originated
				char *pszDescription,		// descriptive text about this sequence
				etSeqBase   *pSeq,			// sequence to generate suffix indexes over
				UINT32 SeqLen);		// sequence length

	bool ProcessMatch(etSeqBase *pProbe,	// pts to start of matching probe sequence
						UINT32 ProbePsn, // offset in probe ( pProbe = &Probe[ProbePsn] )
						etSeqBase *pTarg,	// pts to start of matching target sequence
						UINT32 TargPsn,	// offset in target ( pTarg = &Target[TargPsn] )
						UINT32 MatchLen,	// number of bases known matching
						UINT32 MaximalMatchLen = 20000);	// extend matches left/right to this maximum

	int										// returns number of matching sequences
		Locate(etSeqBase *pProbe,			// sequence containing probe
				  UINT32 ProbePsn,	// start position in probe sequence at which to search 
				  etSeqBase *pTarg,			// target sequence
				  UINT32 *pSfxArray,	// target sequence suffix array
				  UINT32 TargStart,   // position in pTarg (0..n) corresponding to start of suffix array
				  UINT32 SfxLen);		// number of suffixs in pSfxArray


	bool IsLocated(etSeqBase *pProbe,	// sequence containing probe
				  UINT32 ProbePsn,	// start position in probe sequence at which to search
				  UINT32 MinNoMatchLen,	// minimal probe length 
				  etSeqBase *pTarg,			// target sequence
				  UINT32 *pSfxArray,	// target sequence suffix array
				  UINT32 TargStart,  // position in pTarg (0..n) corresponding to start of suffix array
				  UINT32 SfxLen);		// number of suffixs in pSfxArray

				  
	int													// returns number of matches or -1 if errors
		LocateMatches(UINT32 EntryID,				// target entry identifier
					     UINT32 ProbeID,			// probe entry identifier, uniquely identify this probe when ProcessThisMatch is called
					     UINT32 MinMatchLen,		// matches located must be at least this length
 						 UINT32 MinExactLen,		// and contain a subsequence of at least this 100% identity length
						 UINT32 MaxNs,			// maximum number of indeterminate bases allowed
						 UINT32 MaxMutations,		// maximum number of mutations allowed
						 UINT32 MaxMutBlockLen,	// max length of any mutated sequence block
   						 UINT32 MaxNumMatchesPPS,	// max number of matches per MinMatchLen probe subsequence
						 etSeqBase *pProbe,				// probe sequence (ProbeLen bases)
						UINT32 ProbeLen,			// number of bases ptd to by pProbe
						UINT32 Modes,				// bit 0 fwd (sense), bit 1 fwd rev, bit 2 cpl, bit 3 cpl rev (antisense)
						int RsltsFileHandler(int Mode,UINT32 ProbePsn,UINT32 TargPsn,UINT32 MatchLen,UINT32 IdentCnt,UINT32 FHMode),
						bool bChkSelfMatches,			// true if self matches to be filtered out
						bool bFiltComplexity,			// true if low complexity sequences to be filtered
						bool bNotSampled,				// true if sampling not to be used
						bool bAllowContained,			// true if contained sequences are to be reported
						bool bDontExtend,				// true if matches are not to be maximally extended
					    bool bSeqOut,					// true if matching sequences to be output in results
						bool bProbeRepeatSense=false,	// false if bases marked as masked to be treated as repeats, otherwise unmarked bases treated as repeats
						bool bTargRepeatSense=false,	// false if bases marked as masked to be treated as repeats, otherwise unmarked bases treated as repeats
						UINT32 ProbeRepeatPC=100,	// percentage of bases in any probe matches allowed to have been repeat masked (0..100)
    					UINT32 TargRepeatPC=100);	// percentage of bases in any (Suffix) target matches allowed to have been repeat masked (0..100)

	bool LocateNoMatches(UINT32 EntryID,			// target entry identifier
					     UINT32 ProbeID,			// probe entry identifier, uniquely identify this probe when ProcessThisMatch is called
					     UINT32 MinMatchLen,		// non-matches located must be at least this length
   						 UINT32 MaxNumMatchesPPS,	// max number of matches per MinMatchLen probe subsequence
						 etSeqBase *pProbe,				// probe sequence (ProbeLen bases)
						UINT32 ProbeLen,			// number of bases ptd to by pProbe
						UINT32 Modes,				// bit 0 fwd (sense), bit 1 fwd rev, bit 2 cpl, bit 3 cpl rev (antisense)
						int RsltsFileHandler(int Mode,UINT32 ProbePsn,UINT32 TargPsn,UINT32 MatchLen,UINT32 IdentCnt,UINT32 FHMode),
						bool bChkSelfMatches);			// true if self matches to be filtered out
	
	// Returns true if sequence is to be filtered out either because it contains a non-canonical base or the percentage of repeat masked bases is too high
	bool FiltOutRepeats(etRPTMasking RPTMasking,	// how to interpret cRptMskFlg'd bases
						 UINT32 MaskPercent,	// filter out if percentage of repeats is above this percentage (0-100)
						 etSeqBase *pSeq,			// pts to sequence
						 int SeqLen);				// sequence length

	// LocateExacts returns count of all exact matching sequences
	int	LocateExacts(sAMProbe *pProbe,			// probe
				  UINT32 ProbeLen,	    // probe length
				  etSeqBase *pProbeSeq,			// probe sequence
				  UINT32 MinHitLen,	    // minimum target hit length required
  				  // following is call back on all putative hits of probe into suffix entry, allows caller to filter putative hits
				  int	// < 0 to terminate, 0 to stop processing current probe, > 0 to continue 
						(* FilterCore)(sAMProbe *pProbe,				// probe
										UINT32 ProbeOfs,		// probe offset at which exact core starts 
										UINT32 ProbeCoreLen,	// length of exact core
										UINT32 TargetEntryID,	// target suffix array entry 
										UINT32 TargOfs,		// target loci hit	
										etSeqBase *pSeq),			// target sequence, hit starts at &pSeq[TargOfs]

 				  // following is call back on all accepted hits of probe into suffix entry
				  int	// < 0 to terminate, 0 to stop processing current probe, > 0 to continue 
					(* ProcessCore)(sAMProbe *pProbe,			// probe
										UINT32 ProbeOfs,		// probe offset at which hit starts 
										UINT32 HitLen,		// number of bases hit
										UINT32 TargEntryID,	// target suffix array entry 
										UINT32 TargOfs,		// target loci hit	
										etSeqBase *pSeq),			// target sequence, hit starts at &pSeq[TargOfs]
				UINT32 TargEntryID,							// which suffix entry to match probe on
				etRPTMasking RPTMasking = eRPTHmasked,				// how to interpret probe bases with cRptMskFlg set
				UINT32 MaskPercent=80);						// if matches have more than this percent repeats then slough

	// LocateNearExacts
	// Returns count of all nearly exact maximally extended matching sequences
	int LocateNearExacts(sAMProbe *pProbe,		// probe
 				  UINT32 ProbeLen,	    // probe length
				  etSeqBase *pProbeSeq,			// probe sequence
				  UINT32 MinCoreLen,	    // minimal core length from which to extend left + right
				  UINT32 LeftMaxExtend,   // maximum length to extend left
				  UINT32 RightMaxExtend,  // maximum length to extend right
  				  UINT32 LeftMaxMismatches,  // total mismatches allowed in left flank
   				  UINT32 RightMaxMismatches, // total mismatches allowed in right flank
				  UINT32 MinHitLen,	       // putative hits must be of at least this length
				  // following is call back on all putative hits of probe into suffix entry, allows caller to filter putative hits
				  int	// < 0 to terminate, 0 to stop processing current probe, > 0 to continue 
				  (* FilterCore)(sAMProbe *pProbe,					// probe
										UINT32 ProbeOfs,		// probe offset at which probe hit starts 
										UINT32 HitLen,		// number of bases hit
										UINT32 ProbeCoreOfs,	// probe offset at which exact core starts
										UINT32 ProbeCoreLen,	// length of exact core
										UINT32 NumLeftMismatches, // number of mismatches in left flank
										UINT32 NumRightMismatches, // number of mismatches in right flank
										UINT32 TargetEntryID,	// target suffix array entry 
										UINT32 TargOfs,		// target loci hit	
										etSeqBase *pSeq),			// target sequence, hit starts at &pSeq[TargOfs]
				  
				 // following is call back on all accepted hits of probe into suffix entry
				  int		// < 0 to terminate, 0 to stop processing current probe, > 0 to continue 
				  (* ProcessCore)(sAMProbe *pProbe,			// probe
										UINT32 ProbeOfs,		// probe offset at which hit starts 
										UINT32 HitLen,		// number of bases hit
										UINT32 TargEntryID,	// target suffix array entry 
										UINT32 TargOfs,		// target loci hit	
										etSeqBase *pSeq),			// target sequence, hit starts at &pSeq[TargOfs]
				UINT32 TargEntryID,							// which suffix entry to match probe on
				etRPTMasking RPTMasking = eRPTHmasked,				// how to interpret probe bases with cRptMskFlg set
				UINT32 MaskPercent=80);						// if matches have more than this percent repeats then slough

	// GetNMerFreq
				int	GetNMerFreq(UINT32 TargEntryID,			// which suffix entry
						  CSfxNMerFreq *pNMerFreq,			// returned
						  				 // following is call back on all putative hits of probe into suffix entry
						  bool		// false if instance hit to be sloughed 
				  (* InstanceHit)(UINT32 NMerLen,		// number of bases hit
								UINT32 TargEntryID,	// target suffix array entry 
								UINT32 TargOfs,		// target loci hit	
								etSeqBase *pSeq));			// target sequence, hit starts at &pSeq[TargOfs]

			int SetTargEntry(UINT32 TargEntryID);		// which suffix entry (chromosome) to process with IterateExacts

			 int IterateExacts( etSeqBase *pProbeSeq,			// probe
 						 UINT32 ProbeLen,			// probe length
						 UINT32 PrevHitIdx,		// 0 if starting new sequence, otherwise set to return value of previous successful iteration return
						 UINT32 *pHitLoci,		// if match then where to return loci
 	 	 			     int MaxNumExacts = 0);		// only interested in at most this number of exact matches (0 if all exacts are of interest)

			 int						// < 0 if errors, 0 if no matches, 1 if a unique match, 2 if multiple matches
				LocateApproxUniques(int AccumReadHits,				// how many reads have already been matched 
   						etSeqBase *pProbeSeq,			// probe
 						 etSeqBase *pPatternSeq,		// contains pattern to match with, etBaseN represents wildcard bases to match against any in the target
														// will be updated on return with etBaseN's changed to actual subsitution bases
						 UINT32 ProbeLen,			// probe, and also pattern, length
 						 int MinSubCnt,					// minimum allowed pattern wildcard matches
						 int MaxSubCnt,					// maximum allowed pattern wildcard matches
						 UINT32 *pHitLoci,		// if unique match then where to return loci
						 int *pSubCnt);					// used to return number of pattern wildcard substitutions actually required
	
    int GetSfxSeqLen(UINT32 EntryID);				// get length of sequence for requested entry identifier
	int GetSfxSeq(UINT32 EntryID,UINT32 Offset,etSeqBase *pRetSeq,UINT32 Len,bool bClearSeqMsk = false);	// get sequence for entry starting at offset and of length len
	static int CompareExtdMatches( const void *arg1, const void *arg2 );
};
