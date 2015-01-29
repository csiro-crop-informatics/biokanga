#pragma once
#include "./commdefs.h"

#pragma pack(1)
typedef struct TAG_sStepValues {
	int StepPsn;
	int StepValue;
} tsStepValues;

#pragma pack()

class CCSV2BED : public CErrorCodes
{
	int m_NumStepsAllocd;
	tsStepValues *m_pStepValues;
	int m_NumSteps;

public:
	CCSV2BED(void);
	~CCSV2BED(void);

	int	process(char *pszCSV,			// csv file containing structual values
				  char *pszBEDorWIG,	// name of BED or WIG file to create/append
				  char *pszTrackName,	// track name 
				  char *pszDescr,		// description
				  char *pszChrom,		// chromosome 
				  char *pValName,		// entry name
				  unsigned int NumClasses,	// number of classes into which scored values are to be segmented
				  unsigned int Windowsize,	// window size to use for data mean smoothing (1 if no smoothing)
				  unsigned int EntryFld,	// which field in the csv file contains the entry identifier (1..n)
				  unsigned int EntryID,		// entry identifer to match
				  unsigned int EntryOfsFld,	// which field in the csv file contains the entry offset base to left of step
				  unsigned int StartOfs,	// starting offset (base to left of last step)
				  unsigned int EndOfs,		// ending offset (base to left of last step)
  				  int Remap,				// remap offsets by Remap delta
				  unsigned int ValFld,		// which field contains the value of interest
  				  bool bAppend,			// true to append on to existing file otherwise (false) truncate existing file
				  bool bWIG);			// true if to output as WIG file, default is as BED

	int AddStepValue(int StepPsn,long StepValue);
	static int CompareStepPsns( const void *arg1, const void *arg2 );

	bool SmoothData(int WindowSize);		// smooth data using specified window size

	int Output2BEDorWIG(char *pszBED,			// BED or WIG file to create/append
					 char *pszTrackName,	// track name="???"
					 char *pszDescr,		// description="???"
					 char *pszChrom,		// name of chromosome
					 char *pValName,		// name to output as BED or WIG
					 int NumClasses,		// number of classes to bin mean values into
					int Windowsize,			// window size to use for data mean smoothing (1 if no smoothing)
					 bool bAppend,			// true to append on to existing file
					 bool bWIG);			// true if to output as wiggle file, default is as BED

	char *Fmt2FixedDec(int Value,char *pszRet);

};

