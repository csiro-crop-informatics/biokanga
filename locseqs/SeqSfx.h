#pragma once

const int cAllocCPutCoreEls = 100000;	// allocate putative core instances in this many sized chunks

#pragma pack(1)

// Sequence header (one per target sequence)
typedef struct TAG_sSeqHdr {
	INT32 SeqIdx;					// unique Seq index (1..N) - NOTE: not same as SeqID as user allocated identifiers need not be unique or contiguous
	UINT32 TargetRef;			// user assigned target reference 
	INT32 SeqLen;					// Seq length
	etSeqBase *pStart;				// where in concatenated Seqs this Seq sequence starts
} tsSeqHdr;


typedef struct TAG_sCSfxSeqs {
	INT32	SeqIdx;				// unique Seq index - NOTE: not same as SeqID as user allocated identifiers need not be unique or contiguous
	INT32   SeqOfs;				// offset in Seq at which this suffix starts
	etSeqBase *pBase;			// where in concatenated sequences does this suffix starts
} sCSfxSeqs;


typedef struct TAG_sCPutCore {
	UINT32 SeqIdx;				// identifies which target sequence header  
	UINT32 SeqOfs;				// target sequence loci hit	
	UINT32 HitLen;				// number of bases hit
	UINT32 ProbeOfs;			// probe offset at which hit starts 
} tsCPutCore;
#pragma pack()

class CSeqSfx
{
	int m_NumSeqs;				// current number of Seqs
	int m_AllocSeqEls;			// number tsCSeqHdr allocated in m_pSeqs;
	tsSeqHdr *m_pSeqHdrs;			// to hold array of sequence headers

	int m_ConcatSeqLens;			// concatenated Seq length
	int m_AllocSeqLen;				// number bases allocated in m_pConcatSeqs 
	etSeqBase *m_pConcatSeqs;		// to hold concatenated Seqs
	
	int m_NumSfxEls;				// number of elements in m_pSfxEls
	int m_AllocSfxEls;				// number of sCSfxSeqs allocated in m_pSfxEls
	sCSfxSeqs *m_pSfxEls;			// pts to generated sfx array

	int m_NumPutCores;				// number of putative cores
	int m_AllocPutCoreEls;			// number currently allocd for m_pPutCores
	tsCPutCore *m_pPutCores;		// putative cores

	int	AddPutCore(unsigned int ProbeOfs,	// probe offset at which hit starts 
		unsigned int HitLen,			// number of bases hit
		UINT32 SeqIdx,					// identifies which tsSeqHdr 
		unsigned int TargOfs);			// target loci hit


	int								// returns number of cores to be deleted
		FiltContainedPutCores(int ProbeIdx, int TargIdx);
	int							// returns number of cores marked for deletion
		FiltOverlapPutCores(int SeqIdx, int TargIdx);
	int FiltPutCores(bool bRemoveOverlaps = true);		// filter putative cores - remove contained cores and optionally any overlap cores

		// FiltOutRepeats
// Returns true if sequence is to be filtered out either because it contains a non-canonical base or the percentage of repeat masked bases is too high
bool	FiltOutRepeats(etRPTMasking RPTMasking,	// how to interpret cRptMskFlg'd bases
						 unsigned int MaskPercent,	// filter out if percentage of repeats is above this percentage (0-100)
						 etSeqBase *pSeq,			// pts to sequence
						 int SeqLen);				// sequence length

	static int SortPutCores( const void *arg1, const void *arg2);
	static int SortSeqSfxs( const void *arg1, const void *arg2);

public:
	CSeqSfx(void);
	~CSeqSfx(void);
	void Reset(void);

	int Add(unsigned int TargetRef,	// user assigned target reference to associate with any hit to this sequence
			int SeqLen,			// sequence length
			etSeqBase *pProbe);	// Seq sequence
	int GenSfx(void);			// generate suffixes over all concatenated Seqs

	// LocateFirstExact
	// Locates first instance (lowest pSfxArray[] positional index) of exactly matching pTarg 
	int			// index in pSfxArray of exactly matching probe or -1 if no match					
		LocateFirstExact(etSeqBase *pProbe,  // pts to probe sequence
				  unsigned int ProbeLen,		// probe length to exactly match over
				  int SfxLo,					// low index in pSfxArray
				  int SfxHi);					// high index in pSfxArray

	int						
		LocateNearExacts(sAMProbe *pProbe,			// probe
 				  unsigned int ProbeLen,	    // probe length
				  etSeqBase *pProbeSeq,						// ptr to probe sequence
				  unsigned int MinCoreLen,	    // minimal core length from which to extend left + right
				  unsigned int LeftMaxExtend,   // maximum length to extend left
				  unsigned int RightMaxExtend,  // maximum length to extend right
  				  unsigned int LeftMaxMismatches,  // total mismatches allowed in left flank
   				  unsigned int RightMaxMismatches, // total mismatches allowed in right flank
				  unsigned int MinHitLen,	       // putative hits must be of at least this length
				  
				  // following is call back on all putative hits of probe into suffix entry, allows caller to filter putative hits
				  int (* FilterCore)(sAMProbe *pProbe,				// probe
										unsigned int ProbeOfs,		// probe offset at which probe hit starts 
										unsigned int HitLen,		// number of bases hit
										unsigned int ProbeCoreOfs,	// probe offset at which exact core starts
										unsigned int ProbeCoreLen,	// length of exact core
										unsigned int NumLeftMismatches, // number of mismatches in left flank
										unsigned int NumRightMismatches, // number of mismatches in right flank
										unsigned int TargetRef,			// user assigned target reference associated with any hit to this sequence
										unsigned int TargOfs,		// target loci hit	
										etSeqBase *pSeq),			// target sequence, hit starts at &pSeq[TargOfs]
				  
				 // following is call back on all accepted hits of probe into suffix entry
				  int (* ProcessCore)(sAMProbe *pProbe,				// probe
										unsigned int ProbeOfs,		// probe offset at which hit starts 
										unsigned int HitLen,		// number of bases hit
										unsigned int TargetRef,			// user assigned target reference associated with any hit to this sequence
										unsigned int TargOfs,		// target loci hit	
										etSeqBase *pSeq),			// target sequence, hit starts at &pSeq[TargOfs]
				etRPTMasking RPTMasking = eRPTHmasked,				// how to interpret probe bases with cRptMskFlg set
				unsigned int MaskPercent = 80);						// if matches have more than this percent repeats then slough

	// LocateExacts
// Returns count of all exact matching sequences
// These matching sequences are exact (100% identity) and contain no mismatches
			int	LocateExacts(sAMProbe *pProbe,						// probe
				  unsigned int ProbeLen,							// probe length
				  etSeqBase *pProbeSeq,								// ptr to probe sequence
				  unsigned int MinHitLen,							// minimum target hit length required

  				  // following is call back on all putative hits of probe into suffix entry, allows caller to filter putative hits
				  int (* FilterCore)(sAMProbe *pProbe,				// probe
										unsigned int ProbeOfs,		// probe offset at which exact core starts 
										unsigned int ProbeCoreLen,	// length of exact core
										unsigned int TargetRef,				// user assigned target reference associated with any hit to this sequence
										unsigned int TargOfs,		// target loci hit	
										etSeqBase *pSeq),			// target sequence, hit starts at &pSeq[TargOfs]

 				  // following is call back on all hits of probe into suffix entry
				  int (* ProcessCore)(sAMProbe *pProbe,			// probe
										unsigned int ProbeOfs,		// probe offset at which hit starts 
										unsigned int HitLen,		// number of bases hit
										unsigned int TargetRef,				// user assigned target reference associated with any hit to this sequence 
										unsigned int TargOfs,		// target loci hit	
										etSeqBase *pSeq),			// target sequence, hit starts at &pSeq[TargOfs]
				etRPTMasking RPTMasking = eRPTHmasked,				// how to interpret probe bases with cRptMskFlg set
				unsigned int MaskPercent = 80);						// if matches have more than this percent repeats then slough

};
