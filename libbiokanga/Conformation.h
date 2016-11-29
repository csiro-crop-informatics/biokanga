#pragma once
#include "./commdefs.h"

#pragma pack(1)

const int cMaxDiNameLen = 80;   // truncating dinucleotide characteristic names to be no longer than this length
const int cMaxDiSupFieldsLen = 1023;   // truncating dinucleotide characteristic concatenated (comma separated) supplementary fields to be no longer than this length

typedef struct TAG_stsDiStructParam {
	UINT32 ID;					// uniquely identifies this structural characteristic, currently ranges from 1 to 125
	char szName[cMaxDiNameLen+1];	// null terminated characteristic name
	double DiValue[16];			// characteristic value for each dinucleotide from AA through  to TT
	char szSupFields[cMaxDiSupFieldsLen+1];		// comma separated original supplementary fields in order of 'nucleicacid', 'strand', 'author', 'pubyear', 'reference', 'PMID', 'howcreated','type','scaleunit', 'description'
} tsDiStructParam;

// each DiStructID will be one of these enumerations
typedef enum eDiStructStats {
eDiUndefined = 0,			// identifiers range from 1 through to 125 inclusive
eDiTwist = 1,				// Twist - B-DNA
eDiStackingEnergy,			// Stacking Energy - B-DNA
eDiRise,					// Rise - B-DNA
eDiBend,					// Bend - B-DNA
eDiTip,						// Tip - B-DNA
eDiInclination,				// Inclination - DNA
eDiMajorGrooveWidth,		// Major Groove Width - B-DNA
eDiMajorGrooveDepth,		// Major Groove Depth - B-DNA
eDiMajorGrooveSize,			// Major Groove Size - B-DNA
eDiMajorGrooveDistance,		// Major Groove Distance - B-DNA
eDiMinorGrooveWidth,		// Minor Groove Width - B-DNA
eDiMinorGrooveDepth,		// Minor Groove Depth - B-DNA
eDiMinorGrooveSize,			// Minor Groove Size - B-DNA
eDiMinorGrooveDistance,		// Minor Groove Distance - B-DNA
eDiPersistanceLength,		// Persistance Length - B-DNA - Describes the rigidity of DNA. Is defined as P = EI / kT. Where E = related to the stress which develops when the long axis of a rod is strained, I = surface moment of inertia for a right cylinder, k = Boltzmann const., T = abs. temp.si
eDiMeltingTemperature,		// Melting Temperature
eDiProbabilityContactingNucleosomeCore, // Probability contacting nucleosome core - B-DNA
eDiMobilityBendMajorGroove,	// Mobility to bend towards major groove - DNA
eDiMobilityBendMinorGroove,	// Mobility to bend towards minor groove - DNA
eDiPropellerTwist,			// Propeller Twist - B-DNA
eDiClashStrength,			// Clash Strength - B-DNA
eDiEnthalpy,				// Enthalpy - B-DNA
eDiEntropy,					// Entropy - B-DNA
eDiShiftRNA,				// Shift (RNA) - A-RNA
eDiRollDNAproteincomplex,	// Roll (DNA-protein complex) - B_DNA
eDiTwistDNAproteincomplex,	// Twist (DNA-protein complex)
eDiTiltDNAproteincomplex,	// Tilt (DNA-protein complex)
eDiSlideDNAproteincomplex,	// Slide (DNA-protein complex)
eDiHydrophilicityRNA,		// Hydrophilicity (RNA) - RNA Relative hydrophilicities Rf Values for 16 Dinucleoside Monophosphates in 10/90 v/v 1.0 M ammonium acetate/saturated ammonium sulfate Original values are in 3'->5' direction. 5'->3' direction can be found in Lacey et al., Org.L.Ev.Bios.(1983)13,3-42
eDiShiftDNAproteincomplex,	// Shift (DNA-protein complex)
eDiHydrophilicityRNA1,		// Hydrophilicity (RNA) - RNA Rf from paper chromatography (80/18/2:V/V/V). Original values are in 3'->5' direction. 5'->3' direction can be found in Lacey and Mullins, Orig. Life Evol. Biosph.(1983) 13, 3-42.
eDiRiseDNAproteincomplex,	// Rise (DNA-protein complex) - B-DNA
eDiStackingEnergy1,			// Stacking energy - B-DNA
eDiFreeEnergy,				// Free energy - B-DNA
eDiFreeEnergy1,				// Free energy - B-DNA
eDiFreeEnergy3,				// Free energy - B-DNA
eDiTwistDNAproteincomplex1,	// Twist (DNA-protein complex) - B-DNA
eDiFreeEnergy4,				// Free energy - B-DNA
eDiTwist_twist,				// Twist_twist - B-DNA
eDiTilt_tilt,				// Tilt_tilt - B-DNA
eDiRoll_roll,				// Roll_roll - B-DNA
eDiTwist_tilt,				// Twist_tilt - B-DNA
eDiTwist_roll,				// Twist_roll - B-DNA
eDiTilt_roll,				// Tilt_roll - B-DNA
eDiShift_shift,				// Shift_shift - B-DNA
eDiSlide_slide,				// Slide_slide - B-DNA
eDiRise_rise,				// Rise_rise - B-DNA
eDiShift_slide,				// Shift_slide - B-DNA
eDiShift_rise,				// Shift_rise - B-DNA
eDiSlide_rise,				// Slide_rise - B-DNA
eDiTwist_shift,				// Twist_shift - B-DNA
eDiTwist_slide,				// Twist_slide - B-DNA
eDiTwist_rise,				// Twist_rise - B-DNA
eDiTilt_shift,				// Tilt_shift - B-DNA
eDiTilt_slide,				// Tilt_slide - B-DNA
eDiTilt_rise,				// Tilt_rise - B-DNA
eDiRoll_shift,				// Roll_shift - B-DNA
eDiRoll_slide,				// Roll_slide - B-DNA
eDiRoll_rise,				// Roll_rise - B-DNA
eDiStackingEnergy2,			//
eDiTwist1,					//
eDiTilt,					//
eDiRoll,					//
eDiShift,					//
eDiSlide,					//
eDiRise1,					//
eDiSlideStiffness,			//
eDiShiftStiffness,			//
eDiRollStiffness,			//
eDiTiltStiffness,			//
eDiTwistStiffness,			//
eDiFreeEnergy5,				//
eDiFreeEnergy6,				//
eDiFreeEnergy7,				//
eDiGCcontent,				//
eDiPurineAGcontent,			//
eDiKetoGTcontent,			//
eDiAdeninecontent,			//
eDiGuaninecontent,			//
eDiCytosinecontent,			//
eDiThyminecontent,			//
eDiTiltDNAproteincomplex1,	//
eDiRollDNAproteincomplex1,	//
eDiShiftDNAproteincomplex1,	//
eDiSlideDNAproteincomplex1,	//
eDiRiseDNAproteincomplex1,	//
eDiTwist2,					//
eDiTilt1,					//
eDiRoll1,					//
eDiSlide1,					//
eDiTwist3,					//
eDiTilt2,					//
eDiRoll2,					//
eDiShift2,					//
eDiSlide2,					//
eDiRise2,					//
eDiTwist4,					//
eDiWedge,					//
eDiDirection,				//
eDiSlideRNA,				//
eDiRiseRNA,					//
eDiTiltRNA,					//
eDiRollRNA,					//
eDiTwistRNA,				//
eDiStackingEnergyRNA,		//
eDiRiseStiffness,			//
eDiMeltingTemperature1,		//
eDiStackingEnergy3,			//
eDiEnthalpyRNA,				//
eDiEntropyRNA,				//
eDiFreeEnergyRNA,			//
eDiFreeEnergyRNA1,			//
eDiEnthalpyRNA1,			//
eDiEntropyRNA1,				//
eDiRoll3,					//
eDiTilt3,					//
eDiTwist5,					//
eDiRoll4,					//
eDiTwist6,					//
eDiFlexibility_slide,		//
eDiFlexibility_shift,		//
eDiEnthalpy1,				//
eDiEntropy2,				//
eDiFreeEnergy8,			//
eDiNumStatParams			// placeholder which also defines the enumeration range
} teDiStructStats;

typedef enum eOctStructStats {
	eSSenergy = 0,				// minimal energy int(energy * 10000) e.g. -408.2632 ==> -4082632
	eSSminorgroove,				// minor groove int (dimensions * 10000) e.g 10.784 ==> 107840
	eSSmajorgroove,				// inferenced major groove from twist + rise (dimensions * 10000) e.g 18.781 ==> 187810
	eSStwist,					// twist int(angle * 10000) e.g 	37.6262 ==> 376262
	eSSroll,					// roll int(angle * 10000) e.g 	2.4139 ==> 24139
	eSStilt,					// tilt int(tilt * 10000) e.g 	-0.022 ==> -220
	eSSrise,					// rise int(rise * 10000) e.g 	3.1409 ==> 31409
	eSSslide,					// slide int(slide * 10000) e.g 	-0.0968 ==> -968	
	eSSshift,					// shift int(shift * 10000) e.g 	0.0645 ==> 645
	eSSrmsd,					// rmsd int(rmsd * 10000) e.g 	0.3078 ==> 3078
	eSSORChidVal,				// hydroxyl radical cleavage value from ORChid dataset
	eSSNumStatParams			// placeholder which also defines the enumeration range
} teOctStructStats;

typedef struct TAG_sOctStructParam {
	union {
		struct {
		int energy;					// minimal energy int(energy * 10000) e.g. -408.2632 ==> -4082632
		int minorgroove;			// groove int(dimensions * 10000) e.g 10.784 ==> 107840
		int majorgroove;			// inferenced major groove from twist + rise (dimensions * 10000) e.g 18.781 ==> 187810
		int	twist;					// twist int(angle * 10000) e.g 	37.6262 ==> 376262
		int roll;					// roll int(angle * 10000) e.g 	2.4139 ==> 24139
		int tilt;					// tilt int(tilt * 10000) e.g 	-0.022 ==> -220
		int rise;					// rise int(rise * 10000) e.g 	3.1409 ==> 31409
		int slide;					// slide int(slide * 10000) e.g 	-0.0968 ==> -968	
		int shift;					// shift int(shift * 10000) e.g 	0.0645 ==> 645
		int rmsd;					// rmsd int(rmsd * 10000) e.g 	0.3078 ==> 3078
		int orchid;					// hydroxyl radical cleavage value from ORChid dataset
		} Param;
		int Params[eSSNumStatParams];
		};
} tsOctStructParam;

const int cNumParamOctamers = 4 * 4 * 4 * 4 * 4 * 4 * 4 * 4;	// number of structural parameter entries for octamer (4^8)
const int cOctStructParamAllocSize = sizeof(tsOctStructParam) * (1 + cNumParamOctamers);

const int cDiStructParamAllocSize = sizeof(tsDiStructParam) * (1 + 124);

#pragma pack()

class CConformation : public CErrorCodes
{
protected:
	tsOctStructParam *m_pOctStructParams;		// holds all octamer structural base stacking parameters
	char m_szOctStructParamFile[_MAX_PATH];		// file from which octamer structural parameters have been loaded

	int m_NumDiStructParams;					// number of dimer structural base stacking parameters loaded
	tsDiStructParam *m_pDiStructParams;			// holds all dimer structural base stacking parameters
	char m_szDiStructParamFile[_MAX_PATH];		// file from which dimer structural parameters have been loaded

public:
	CConformation(void);
	~CConformation(void);
	bool StructParamsLoaded(void);
	teBSFrsltCodes LoadStructOctamersParams(char *pszStructParamsFile); // load structural parameters + ORChid from file
	teBSFrsltCodes LoadStructDimersParams(char *pszStructParamsFile);	// load structural parameters from file

	int StructParamIdx(etSeqBase *pOctamer);		// sequence
	teBSFrsltCodes GetSequenceConformation(teOctStructStats Param,	// which structural parameter value to return
				 unsigned int iStartOfs,			// initial starting offset (0..n) in pSeq
				  unsigned int iNumSteps,		  // number of steps (0 for all) to process starting at pSeq[iStartPsn]|pSeq[iStartPsn+1]
  				  unsigned int SeqLen,			  // total length of sequence
				  etSeqBase *pSeq,				  // sequence to be processed
 				  int *pRetConfValue,			  // where to return conformation
				  int UndefBaseValue);

	int StructValue(teOctStructStats Param,	// which structural parameter value to return
			unsigned int Step,				// which step in sequence to return structural value for
			unsigned int SeqLen,			// total length of sequence
			etSeqBase *pSeq,				// sequence to be processed
			int UndefBaseValue);			// value to return for undefined or indeterminate ('N') bases 

	
	int PsudeoRandomise();			// reproducibly pseudo-randomise the conformation characteristics

};
