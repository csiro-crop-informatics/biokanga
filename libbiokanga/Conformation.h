#pragma once
#include "./commdefs.h"

#pragma pack(1)

typedef enum eStructStats {
	eSSenergy = 0,				// minimal energy int(energy * 10000) e.g. -408.2632 ==> -4082632
	eSSminorgroove,				// minor groove int (dimensions * 10000) e.g 10.784 ==> 107840
	eSStwist,					// twist int(angle * 10000) e.g 	37.6262 ==> 376262
	eSSroll,					// roll int(angle * 10000) e.g 	2.4139 ==> 24139
	eSStilt,					// tilt int(tilt * 10000) e.g 	-0.022 ==> -220
	eSSrise,					// rise int(rise * 10000) e.g 	3.1409 ==> 31409
	eSSslide,					// slide int(slide * 10000) e.g 	-0.0968 ==> -968	
	eSSshift,					// shift int(shift * 10000) e.g 	0.0645 ==> 645
	eSSrmsd,					// rmsd int(rmsd * 10000) e.g 	0.3078 ==> 3078
	eSSORChidVal,				// hydroxyl radical cleavage value from ORChid dataset
	eSSNumStatParams			// placeholder which also defines the enumeration range
} teStructStats;

typedef struct TAG_sStructParam {
	union {
		struct {
		int energy;					// minimal energy int(energy * 10000) e.g. -408.2632 ==> -4082632
		int minorgroove;			// groove int(dimensions * 10000) e.g 10.784 ==> 107840
		int	twist;					// twist int(angle * 10000) e.g 	37.6262 ==> 376262
		int roll;					// roll int(angle * 10000) e.g 	2.4139 ==> 24139
		int tilt;					// tilt int(tilt * 10000) e.g 	-0.022 ==> -220
		int rise;					// rise int(rise * 10000) e.g 	3.1409 ==> 31409
		int slide;					// slide int(slide * 10000) e.g 	-0.0968 ==> -968	
		int shift;					// shift int(shift * 10000) e.g 	0.0645 ==> 645
		int rmsd;					// rmsd int(rmsd * 10000) e.g 	0.3078 ==> 3078
		int orchid;				// hydroxyl radical cleavage value from ORChid dataset
		} Param;
		int Params[eSSNumStatParams];
		};
} tsStructParam;

const int cNumParamOctamers = 4 * 4 * 4 * 4 * 4 * 4 * 4 * 4;	// number of structural parameter entries for octamer (4^8)
const int cStructParamAllocSize = sizeof(tsStructParam) * (1 + cNumParamOctamers);



#pragma pack()

class CConformation : public CErrorCodes
{
protected:
	tsStructParam *m_pStructParams;				// holds all structural base stacking parameters
	char m_szStructParamFile[_MAX_PATH];		// file from which octamer structural parameters have been loaded

public:
	CConformation(void);
	~CConformation(void);
	bool StructParamsLoaded(void);
	teBSFrsltCodes LoadStructParams(char *pszStructParamsFile); // load structural parameters + ORChid from file

	int StructParamIdx(etSeqBase *pOctamer);		// sequence
	teBSFrsltCodes GetSequenceConformation(teStructStats Param,	// which structural parameter value to return
				 unsigned int iStartOfs,			// initial starting offset (0..n) in pSeq
				  unsigned int iNumSteps,		  // number of steps (0 for all) to process starting at pSeq[iStartPsn]|pSeq[iStartPsn+1]
  				  unsigned int SeqLen,			  // total length of sequence
				  etSeqBase *pSeq,				  // sequence to be processed
 				  int *pRetConfValue,			  // where to return conformation
				  int UndefBaseValue);

	int StructValue(teStructStats Param,	// which structural parameter value to return
			unsigned int Step,				// which step in sequence to return structural value for
			unsigned int SeqLen,			// total length of sequence
			etSeqBase *pSeq,				// sequence to be processed
			int UndefBaseValue);			// value to return for undefined or indeterminate ('N') bases 

	
	int PsudeoRandomise();			// reproducably psudeo-randomise the conformation characteristics

};
