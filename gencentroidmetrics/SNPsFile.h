#pragma once

// SNP class - as currently only single nucleotide polymorphisms are of interest then this
// class is not complete and contains limited functionality

const int cAllocPMChunk = 1000000;		// allocation chunck size
const int cMaxPMChromNameLen= 25;		// max chromosome name length
const int cMaxPMNumChroms   = 200;		// maximum number of chromosomes supported

#pragma pack(1)
// struct used to map from text to enum, and back from enum to text
typedef struct {
		const char *pszName;
		int EnumValue;
	} tsSNPEnumValue;

// molecular type
typedef enum eSNPMoleType {
	eSNPMTUnknown = 0,
	eSNPNTGenomic,
	eSNPMTcDNA,
	eSNPMTmito,
	eSNPMTchloro
	} teSNPMoleType;

// Class
typedef enum eSNPClass {
	eSNPCunknown = 0,
	eSNPCsnp,
	eSNPCindel,
	eSNPChet,
	eSNPCmicrosat,
	eSNPCnamed,
	eSNPCnovariation,
	eSNPCmixed,
	eSNPCmnp
} teSNPClass;

// actually a set rather than an enumeration
typedef enum eSNPValid {
	eSNPVunknown = 0,
	eSNPVotherpop,
	eSNPVbyfrequency,
	eSNPVbycluster,
	eSNPVby2hit2allele,
	eSNPVbyhapmap,
	eSNPVgenotype
} teSNPValid;


typedef enum eSNPFunc {
	eSNPFunknown = 0,
	eSNPFlocusregion,
	eSNPFcoding,
	eSNPFcodingsynon,
	eSNPFcodingnonsynon,
	eSNPFmrnautr,
	eSNPFintron,
	eSNPFsplicesite,
	eSNPFreference,
	eSNPFexception
} teSNPFunc;

typedef enum eSNPSource {
	eSNPSunknown = 0,
	eSNPSdbSnp,
	eSNPSAffy10K,
	eSNPSAffy10Kv2,
	eSNPSAffy50K_HindIII,
	eSNPSAffy50K_XbaI
} teSNPSource;

typedef enum eSNPLocType {
	eSNPLunknown = 0,
	eSNPLrange,
	eSNPLexact,
	eSNPLbetween
} teSNPLocType;

typedef struct TAG_sSNPs {
	int ChromID;				// on which chromosome is this SNP is located
	int StartOfs;				// locii start (0..n) on chromosome
	char Strand;				// on which strand
	unsigned char PMs;			// bitmap of observed polymorphisms
} tsSNPs;


typedef struct TAG_sPMChrom {
	int ChromID;				// chromosome identifier
	int FirstSNPID;				// identifier of 1st SNP on this chromosome
	char szName[cMaxPMChromNameLen+1];	// chromosome name
} tsPMChrom;

#pragma pack()

class CSNPsFile
{
	int m_hFile;						// opened ascii SNPs file currently being parsed
	static tsSNPEnumValue m_SNPMoleTypeEnums[];
	static tsSNPEnumValue m_SNPClassEnums[];
	static tsSNPEnumValue m_SNPValidEnums[];
	static tsSNPEnumValue m_SNPFuncEnums[];
	static tsSNPEnumValue m_SNPSourceEnums[];
	static tsSNPEnumValue m_SNPLocTypeEnums[];

	int m_NumPMs;				// current number of polymorphisms in m_pPMs
	int m_AllocdPMs;			// how many instances of tsSNPs m_pPMs can hold before realloc required
	tsSNPs *m_pPMs;				// pts to array of tsSNPs

	int m_NumChroms;			// current number of chromosomes in m_Chroms
	tsPMChrom m_Chroms[cMaxPMNumChroms];

	int AddChrom(char *pszChrom);
	static int SortPMs(const void *arg1, const void *arg2);  // used to sort by ChromID->Start

public:
	CSNPsFile(void);
	~CSNPsFile(void);
	void Reset(void);
	teBSFrsltCodes ProcessSNPsFile(char *pszFileName);	// file to process
	teBSFrsltCodes Parse(void);	// parse the opened SNPs ascii file into tsSNPs 
	int								// returns count of SNPs added, 0 if SNP filtered out, or -1 if internal errors
	AddSNP(char *pszChrom,		    // on which chromosome
				int StartOfs,		// start offset (0..n)
				char Strand,		// which strand ('?','+','-')
				char *pszObserved,	// observed mutations separated by '/'
				char *pszMoleType,	// enum: Sample used to find this variant
				char *pszClass,		// enum: Variant classification
				char *pszLocType);	// enum: Describes how a segment of the reference assembly must be altered to represent the variant SNP allele

    int GetNextSNP(int CurSNPID,int ChromID=0); // returns next SNP identifier, identifiers are ordered by chromosome then offset
	int GetChromID(int SNPID);		// returns chromosome identifier on which SNPID is located
	char *LocateChrom(int ChromID);	// returns chromosome name or NULL if unable to locate 
	int LocateChrom(char *pszChrom); // returns chromosome identifier or -1 if unable to locate
	int GetOffset(int SNPID);		// returns offset on chromosome at which SNPID is located
	char GetStrand(int SNPID);		// returns strand for this SNP

	bool IsPM(int SNPID,teSeqBases Base); // returns true if Base is a 	polymorphism
	int LocateSNP(int ChromID,int Offset); // returns SNPID of any polymorphism at this offset 
	
	teSNPMoleType MolType2Enum(char *pszTxt);
	teSNPClass Class2Enum(char *pszTxt);
	teSNPValid Valid2Enum(char *pszTxt);
	teSNPFunc Funct2Enum(char *pszTxt);
	teSNPSource Source2Enum(char *pszTxt);
	eSNPLocType LocType2Enum(char *pszTxt);

	char *MolTypeEnum2Txt(teSNPMoleType Enum);
	char *ClassEnum2Txt(teSNPClass Enum);
	char *ValidEnum2Txt(teSNPValid Enum);
	char *FunctEnum2Txt(teSNPFunc Enum);
	char *SourceEnum2Txt(teSNPSource Enum);
	char *LocTypeEnum2Txt(eSNPLocType Enum);
};
