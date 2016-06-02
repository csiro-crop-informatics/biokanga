#pragma once

// following group are limits are reasonably arbitrary - could be increased with little risk but memory requirements could be an issue
const UINT32 cMaxNumRefSeqs = 0x0ffffff;		// can process up to this many reference sequences (16M)
const UINT64 cMaxTotRefSeqLens = 0x07ffffffff;	// which total in length to no more than this many bases (28Gbp)
const UINT32 cMinRefSeqLen = 100;				// only accepting individual reference sequences of at least this length (100bp)
const UINT32 cMaxRefSeqLen = 0x00fffffff;		// only accepting individual reference sequences no longer than this length (256Mbp)

#pragma pack(1)

// multiple alignment columns
typedef struct TAG_sMAlignConCol {
	UINT32 RefID;			// column is in this reference sequence
	UINT64 ColIdx;			// index + 1 in m_pMACols[] at which this alignment column starts
    UINT64 NxtColIdx;		// index + 1 in m_pMACols[] at which next alignment column starts, 0 if this is the last alignment column in RefID reference sequence
    UINT64 PrevColIdx;		// index + 1 in m_pMACols[] at which previous alignment column starts, 0 if this is the first alignment column in RefID reference sequence
	UINT32 Ofs;				// alignment col is for this reference sequence relative offset (1..n)
	UINT16  Extn;			// if > 0 then relative offset to RefLoci of an inserted column 
	UINT32  Depth;			// currently with this total coverage 
	UINT8   ConsBase;       // initialised with the reference base, later updated when all alignments known with the consensus base
	UINT32  BaseCnts[eBaseInDel+1];	// counts for number of bases, indexed by eBaseA..eBaseInDel  	
} tsMAlignConCol;


typedef struct TAG_sMARefSeq {
	UINT32 RefID;					// uniquely identifies this reference sequence
	UINT32 SeqLen;					// reference sequence is this length
	UINT64 StartColIdx;				// sequence starts, inclusive, from this this alignment column in m_pMACols[]
	UINT64 EndColIdx;				// sequence ends, inclusive, at this this alignment column in m_pMACols[]
	} tsMARefSeq;


#pragma pack()

class CMAConsensus
{

	bool m_bStartedMultiAlignments;     // set true after 1st call to StartMultiAlignments(), reset after call to GenMultialignConcensus() 

	UINT32 m_NumRefSeqs;			// consensus is over this many reference sequences
	UINT32 m_AllocdRefSeqs;			// m_pRefSeqs currently allocated to hold this many tsRefSeqs
	tsMARefSeq *m_pRefSeqs;			// allocated to hold reference sequence detail

	UINT64 m_MATotRefSeqLen;		// can accept up to total reference sequence length 
	UINT64 m_MACurTotRefSeqLen;		// actual current total reference sequence length
	UINT64 m_MACurCols;				// actual current number of multialignment columns used
	UINT64 m_AllocMACols;			// this number of columns have been allocated
	size_t m_AllocMAColsSize;		// allocated memory size for multialignment cols
	tsMAlignConCol *m_pMACols;				// pts to allocated memory for multialignment cols 

	tsMAlignConCol *						// inserted column or NULL if errors inserting
		InsertCol(UINT64 PrevCol);	// allocate and insert new column after PrevCol

public:
	CMAConsensus();
	~CMAConsensus();

	void Reset(void);							// reset state back to that immediately following instantiation

	int Init(UINT32 NumRefSeqs,					// max number of reference sequences which will be added
			 UINT64 TotRefSeqLen);				// max total bases of all reference sequences which will be added

	UINT32								// returned reference sequence indentifier to be used when adding alignments to this reference sequence
		AddRefSeq(UINT32 SeqLen,		// reference sequence to which other sequences are aligned is this length
					etSeqBase *pRefSeq); // reference sequence 

	int												// eBSFSuccess or error
		AddMultiAlignment(UINT32 RefSeqID,			// alignment is against this sequence
					  UINT32 RefStartOfs,			// alignment starts at this reference sequence offset (1..n)
					  UINT32 RefEndOfs,				// alignment ends at this reference sequence offset inclusive
					  UINT32 ProbeStartOfs,			// alignment starts at this probe sequence offset (1..n)
					  UINT32 ProbeEndOfs,			// alignment ends at this probe sequence offset inclusive
					  etSeqBase *pProbeSeq,			// alignment probe sequence
					  UINT32 NumMAAlignOps,			// number of alignment operators
					   tMAOp *pMAAlignOps);			// alignment operators

	int	GenMultialignConcensus(void);          // generate multiple alignment consensus over all multialignment columns at m_pMACols

	int											// number of bases in consensus 
		GetConsensus(UINT32 RefSeqID);			// for this reference sequence identifier

	int											// number of bases returned in pRetBases 
		GetConsensus(UINT32 RefSeqID,			// consensus sequence for this reference identifier
					    UINT32 RefStartOfs,		// starting offset 1..N
						UINT32 RefLen,			// return at most this many bases
						etSeqBase *pRetBases);  // where to return bases

};

