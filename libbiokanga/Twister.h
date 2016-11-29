#pragma once
#include "./commdefs.h"

const int cNumHistBins = 200;		// number of bins used to hold histogram counts
const int cMaxNumDPsInSeg = cDSFMaxDataPts; // max number of data points in any one segment

#pragma pack(1)

typedef struct TAG_sStructStats {
	int Mean;				// the mean (average) over specified window
	int Min;				// the minimum
	int Max;				// the maximum 
	int StdDev;				// standard deviation * 100
} tsStructStats;

#pragma pack()

class CTwister : public CConformation
{
	
	int GenStructParamStats(void);				// generate basic structure parameter stats
	int Interpolate(unsigned int Step,			// which step to interpolate (1..7)
			etSeqBase *pSeq,					// known octamer sequence left (step 5..7) or right (step 1..3) filled with eBaseA
			tsOctStructParam *pRetStructParam);	// where to returned interpolated structural parameters

	int											// returned index (-1 if any base indeterminate) 
	StructParamIdx(etSeqBase *pSeq);			// octamer sequence
	int *MapStructParam2Ptr(teOctStructStats Param,// which parameter
					 tsOctStructParam *pStruct);   // which parameters instance

public:
	CTwister(void);
	~CTwister(void);

	int LoadStructParams(char *pszStructParamsFile); // load structural parameters from file


	int  GetStructParams(unsigned int Step,	 // which step (1..SeqLen-1)
				unsigned int SeqLen,			 // sequence length (8..n)
				etSeqBase *pSeq,				 // sequence
				tsOctStructParam *pRetStructParam); // where to return structural parameters

	int
		GetStructParam(unsigned int Step,	// which step (1..SeqLen-1)
				unsigned int SeqLen,		// sequence length (8..n)
				etSeqBase *pSeq,			// sequence
				teOctStructStats Param);		// which structural parameter value to return


	int
	ProcessSequence(int hRslts,					// file to write sequence structure into
				  unsigned int RefID,			// non-unique identifier by which related sequences can be later associated 
				  char *pszSpecies,				// species description
				  unsigned int EntryID,			// sequence entry identifier (could be chromosome) for this sequence
				  char *pszEntry,				// entry or chromosome description
				  unsigned int EntryOfs,		// entry or chromosome offset (0..n) corresponding to left base of 1st step
				  unsigned int iStartOfs,		// initial starting offset (0..n) in pSeq
				  unsigned int iNumSteps,		// number of steps to process starting at pSeq[iStartPsn]|pSeq[iStartPsn+1]
  				  unsigned int SeqLen,			// total length of sequence
				  etSeqBase *pSeq,				// sequence to be processed
				bool bXML);						// results file type: false == CSV, true == XML

	int
	GetSequenceConformation(teOctStructStats confparam,	// get values for this conformational parameter
				  int iStartOfs,			// initial starting offset (0..n) in pSeq
				  int iNumSteps,			// number of steps (0 for all) to process starting at pSeq[iStartPsn]|pSeq[iStartPsn+1]
  				  int SeqLen,				// total length of sequence
				  etSeqBase *pSeq,			// sequence to be processed
				  int *pRetValues);			// where to return step conformational values

	char *
	Fmt2FixedDec(int Value,char *pszRet = NULL); // format supplied value into a string a d.dddd where Value is assumed to have 4 decimal places

	int	CalcTwistStats(unsigned int iStartOfs,		// initial starting offset (0..n) in pSeq
					unsigned int iNumSteps,		// number of steps to process starting at pSeq[iStartPsn]|pSeq[iStartPsn+1]
				    unsigned int SeqLen,		// total length of sequence
					etSeqBase *pSeq,			// start of sequence
					teOctStructStats StructParam,	// which structural parameter to gen stats for
					tsStructStats *pStats);		// returned stats

    int	CalcDiffStats(unsigned int iRefStartOfs,	// initial starting offset (0..n) in pRefSeq
					unsigned int iNumSteps,		// number of steps to process starting at pSeq[iStartPsn]|pSeq[iStartPsn+1]
				    unsigned int RefSeqLen,		// total length of reference sequence
					etSeqBase *pRefSeq,			// start of reference sequence
					unsigned int iRelStartOfs,	// initial starting offset (0..n) in pRelSeq
					unsigned int RelSeqLen,		// total length of reference sequence
					etSeqBase *pRelSeq,			// start of relative sequence
					teOctStructStats StructParam,	// which structural parameter to gen stats for
					tsStructStats *pStats); 	// returned stats
	
	int
	CalcStats(char *pszSpecies,					// species - e.g 'hg17'
				  unsigned int EntryID,			// sequence entry identifier for this sequence
				  char *pszChrom,				// from this chromosome
				  unsigned int iStartOfs,		// initial starting offset (0..n) in pSeq 
				  unsigned int iNumSteps,		// number of steps to process starting at pSeq[iStartPsn]|pSeq[iStartPsn+1]
				  unsigned int SeqLen,			// total length of sequence
				  etSeqBase *pSeq,				// sequence 
  				  teOctStructStats StructParam,	// which structural parameter to gen stats for
  				  unsigned int WinSteps, 		// sliding window size in dinucleotide steps (5..100000)
				  int hRslts,					// results file to write into
				  bool bXML);					// results file type: false == CSV, true == XML
				  

	void OutResult(bool bXML,				// true==XML, false==CSV
					int hRslts,				// file to output into
					char *pszSpecies,		// species name e.g hs17
					char *pszChrom,			// chromosome e.g chr1
					char *pszStat,			// statistics type e.g mean
					int EntryID,			// entry ident corresponding to pszChrom
					int Stat,				// statistics type: 0 == Min, 1== Max, 2==Mean, 4==StdDev
					int Steps,				// number of dinucleotide steps stats generated over (1..n)
					int Value,				// stats value
					int Freq);				// frequency or number of instances of Value

	int CalcStruct(char *pszBioSeqFile,		// sequence file
				char *pszSpecies,			// reference species name
				int iSeqEntry,				// which entry (1..n) or all entries (0) in sequence file
				int iStartOfs,				// process from start position (0..n)
				int  iNumSteps,				// process this many steps starting at Seq[StartOfs]|Seq[StartOfs+1]
				char *pszResultsFile,		// where to write results
				bool bXML,					// results file type: false == CSV, true == XML
				char *pszStructParamsFile=NULL);// structural parameters file, NULL if to use existing loaded

	int GenRefStructDist(char *pszResultsFile,				  // file to contain results
				bool bXML,					// results file type: false == CSV, true == XML
				char *pszStructParamsFile); // structural parameters file, NULL if to use existing loaded
	void SetMissing(tsOctStructParam *pRetStructParam);
	
	tsStructStats m_StructParamStats[eSSNumStatParams];	// structure parameter basic stats

	static const char *MapStructParam2Txt(teOctStructStats Param);

};
