#pragma once

#include "../libbiokanga/commhdrs.h"

const int cExpNumInputFields = 4;		// expecting input transcriptome input file to contain fileds <FeatID>,<FeatName>,<FeatLen>,<FeatInsts>
const int cTransRegionAllocNum = 20000;	// alloc for tsTransRegions in this many increments
const int cTransSeqAllocSize = 0x3fffffff; // allocate to hold transcripts in this sized increments

const int cDfltTargFragLen = 300;			// default targeted maximum fragment length
const int cDfltFragLenBins = 10;			// default number of fragment length bins
const int cDfltMinFragLen = 125;			// default minimum length fragment
const int cDfltMaxFragLen = 175;			// default maximum length fragment
const int cDfltStericOcc5Len = 30;			// default 5' fragmentation steric occlusion length
const int cDfltStericOcc3Len = 30;			// default 3' fragmentation steric occlusion length

const UINT8 cFragStartFlag = 0x01;				// flags 5' start of a fragment
const UINT8 cFragEndFlag = 0x02;				// flags 3' end of a fragment
const UINT8 cFragStericOcc5Flag = 0x04;			// flags 5' steric occlusion	
const UINT8 cFragStericOcc3Flag = 0x08;			// flags 3' steric occlusion	



#pragma pack(1)
// Each transcribed region will have a number of transcripts transcribed from that region
typedef struct TAG_sTransRegion {
	int ElIdx;			// element identifier corresponding to nTh element in input source file
	int NumInsts;		// number of instances of this transcript
	int TransLen;		// transcript length
	size_t TransSeqOfs;	// offset in m_pTransSeqs of where the first transcipt starts for this transcribed region
	} tsTransRegion;

#pragma pack()

class CTranscriptome
{
	char m_szTransFile[_MAX_PATH];			// input transcriptome file
	char m_szFragsFile[_MAX_PATH];			// output fragment counts file
	char m_szDistFile[_MAX_PATH];			// output fragment length distributions file
	int m_hTransFile;						// input transcriptome file handle
	int m_hFragsFile;						// output fragment counts handle
	int m_hDistFile;						// output fragment length distribution handle
	CCSVFile *m_pCSV;							// for parsing the input CSV transcriptome file 

	int m_StericOcc5Len;			// 5' fragmentation steric occlusion length
	int m_StericOcc3Len;			// 3' fragmentation steric occlusion length
	int m_MinFragLen;				// count fragments with >= this length and
	int m_MaxFragLen;				// <= this length
	int m_NumFragLenBins;			// bin fragment lengths into this many bins so as to obtain a distribution
	int m_TargFragLen;				// fragment until mean fragmentation length is <= this				

	int m_NumRegionsParsed;			// number of transcribed regions parsed from file
	int m_NumTransRegions;			// number of transcribed regions
	int m_NumAllocdTransRegions;	// currently allocated to hold this many transcribed regions
	size_t m_MemAllocdTransRegions; // memory size currently allocated to hold m_NumAllocdTransRegions transcribed regions
	tsTransRegion *m_pTransRegions;	// all transcribed regions


	int m_NumTranscripts;			// current total number of transcripts
	size_t m_CurTotTransLen;		// current total concatenated transcript length
	size_t m_MemAllocdTrnscripts;	// memory size currently allocated to hold 
	UINT8 *m_pTranscripts;			// to hold concatenated transcripts pTranscripts

	int AddTranscripts(int TranRegionID,	// uniquely identifies transcribed region 
			int TransLen,				// transcripts are of this length 
			int NumInsts);				// and there are this many instances

public:
	CTranscriptome(void);
	~CTranscriptome(void);
	
	int Reset(bool bSync = true);	// if bSync true then force file sync to disk before closing, then release allocated memory etc

	int
		Fragmentate(char *pszTransFile,		// transcripts and their abundance are from this CSV file
                     char *pszFragsFile,	// fragmentation counts to this file
					 char *pszDistFile,		// output fragment length distributions (CSV) file
					 int TargFragLen,				// fragment until mean fragmentation length is <= this				
					 int NumFragLenBins = cDfltFragLenBins, // bin fragment lengths into this many bins so as to obtain a distribution
					 int MinFragLen = cDfltMinFragLen,		// count fragments with >= this length and
					 int MaxFragLen = cDfltMaxFragLen,		// <= this length
			 		 int StericOcc5Len = cDfltStericOcc5Len,		// 5' steric occlusion region length
					 int StericOcc3Len = cDfltStericOcc3Len);	// 3' steric occlusion region length


	

};

