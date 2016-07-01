#pragma once
#include "./commdefs.h"

const int cMaxProb100Int = 10000;	   // scaled int value corresponding to prob of 1.00 (int)(prob * 10000.0)
const int cMaxStatTransPeriods = 1000; // max number of periods over which stationary probabilities can be calculated
#pragma pack(1)

typedef struct TAG_sCentroidParam {
		union {
			struct {
				int IGFixProb;				// intergenic fixation
				int USFixProb;				// upstream fixation
				int UTR5FixProb;			// UTR5 fixation
				int CDSFixProb;				// CDS fixation
				int IntronFixProb;			// intronic fixation
				int UTR3FixProb;			// UTR3 fixation
				int DSFixProb;				// downstream fixation
			} Param;
			int Params[eSSNumStatParams];
		};
} tsCentroidParam;

typedef struct TAG_sTransProbMatrix {
	double Els[7][4][4];						// 7x4x4 matrix 
												// [n][x][x] ==region 
												// [x][n][x] == a,c,g,t current state
												// [x][x][n] == a,c,g,t state to transited into
												
	} tsTransProbMatrix;

const int cNumCentroidParams = 16384;
const int cCentroidParamAllocSize = sizeof(tsCentroidParam) * cNumCentroidParams;
const int cTransMatixAllocSize = sizeof(tsTransProbMatrix) * cNumCentroidParams;



#pragma pack()

class CCentroid: public CErrorCodes
{
	int m_NumCentroids;							// how many centroids have been loaded - will be 0 (none), 4 (1-mer),64 (3-mer),1024 (5-mer) or 16384 (7-mer)
	int m_CentroidNMer;							// m_NumCentroids mapped into N-Mer
	tsCentroidParam *m_pCentroidParams;			// holds all centroid parameters

	int m_NumProbMatrices;						// how many transistional probabilities matrices have been loaded
	int m_TransMatricesNMer;					// m_NumProbMatrices mapped into N-Mer
	tsTransProbMatrix *m_pTransMatrices;		// holds all transistional probabilities matrices

	char m_szCentroidParamFile[_MAX_PATH];		// file from which 7-mer centroid parameters have been loaded
	char m_szTransMatricesFile[_MAX_PATH];		// file from which transitional matrices have been loaded

public:
	CCentroid(void);
	~CCentroid(void);
	bool CentroidParamsLoaded(void);
	teBSFrsltCodes LoadCentroidParams(char *pszCentroidParamsFile); // load centroid parameters from file
	teBSFrsltCodes LoadTransMatrices(char *pszTransMatricesFile); // load transitional probabilities matrix from file

	int OligoIdx(etSeqBase *pOligo);	// determine oligo index to use from sequence
	int OligoIdx(char *pszOligo);		// determine oligo index to use from sequence

	teBSFrsltCodes GetSequenceCentroids(teFuncRegion Param,		// which centroid parameter value to return
				 unsigned int iStartOfs,		  // initial starting offset (0..n) in pSeq
				  unsigned int iNumSteps,		  // number of steps (0 for all) to process starting at pSeq[iStartPsn]|pSeq[iStartPsn+1]
  				  unsigned int SeqLen,			  // total length of sequence
				  etSeqBase *pSeq,				  // sequence to be processed
  				  int *pRetCentValue,			  // where to return centroid values (0...10000 where 10000 == 1.0)
				  int UndefBaseValue = -1);			  // value to return for undefined or indeterminate ('N') bases


	int						// centroid value (0...10000 where 10000 == Pfixed = 1.0)
		CentroidValue(teFuncRegion Param,	// which centroid value to return
			unsigned int Step,					// which step in sequence to return centroid value for
			unsigned int SeqLen,				// total length of sequence
			etSeqBase *pSeq,					// sequence to be processed
			int UndefBaseValue = -1);				// value to return for undefined or indeterminate ('N') bases 


	int EvolveSeq(teFuncRegion Region, // sequence is in this region
				   etSeqBase *pSeq,			// sequence to be evolved (input and output)
				   int SeqLen,				// sequence length
				   int RandLociSeed = -1,		// -1 == use rand() as seed, >= 0 then use this as random loci seed
				   int RandBaseSeed = -1);		// -1 == use rand() as seed, >= 0 then use this as random base seed


	teBSFrsltCodes StationaryCentroidProbs(teFuncRegion Region, // region to calc stationary probs for
			    int OligoIdx,		// uniquely identifies oligo (0..n) see CentroidParamIdx()
			    int NumPeriods,			// number of time periods
				double *pProbA,		// where to return probabilities in each time period for A
				double *pProbC,		// where to return probabilities in each time period for C
				double *pProbG,		// where to return probabilities in each time period for G
				double *pProbT);		// where to return probabilities in each time period for T


	teBSFrsltCodes StationarySeqProbs(teFuncRegion Region, // region to calc stationary probs for
				   char *pszSeq,	// sequence to calc stationary probs over
				   int SeqLen,		// sequence length
				   int Period,		// which period is of interest (1..n)
				   double *pToRetA,	// where to return probabilities of A at Period N
				   double *pToRetC, // where to return probabilities of C at Period N
				   double *pToRetG, // where to return probabilities of G at Period N
				   double *pToRetT); // where to return probabilities of T at Period N

	teBSFrsltCodes StationarySeqProbs(teFuncRegion Region, // region to calc stationary probs for
				   etSeqBase *pSeq,	// sequence to calc stationary probs over
				   int SeqLen,		// sequence length
				   int Period,		// which period is of interest (1..n)
				   double *pToRetA,	// where to return probabilities of A at Period N
				   double *pToRetC, // where to return probabilities of C at Period N
				   double *pToRetG, // where to return probabilities of G at Period N
				   double *pToRetT); // where to return probabilities of T at Period N
};
