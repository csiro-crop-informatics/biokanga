#pragma once
#include "./commdefs.h"

// Reference: http://mblab.wustl.edu/GTF22.html
//GTF stands for Gene transfer format. It borrows from GFF, but has additional structure that warrants a separate definition and format name. 
//Structure is as GFF, so the fields are: 
//<seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes] [comments] 

//Here is a simple example with 3 translated exons. Order of rows is not important. 

//381 Twinscan  CDS          380   401   .   +   0  gene_id "001"; transcript_id "001.1";
//381 Twinscan  CDS          501   650   .   +   2  gene_id "001"; transcript_id "001.1";
//381 Twinscan  CDS          700   707   .   +   2  gene_id "001"; transcript_id "001.1";
//381 Twinscan  start_codon  380   382   .   +   0  gene_id "001"; transcript_id "001.1";
//381 Twinscan  stop_codon   708   710   .   +   0  gene_id "001"; transcript_id "001.1";
//The whitespace in this example is provided only for readability. In GTF, fields must be separated by a single TAB and no white space


const int cGTFMaxAttribs	= 50;	// max number of attributes
const int cGTFMaxValLen		= 1024;		// maximum length (chars) of any single attribute value supported
const int cGTFMaxLineLen  = (1024 + (cGTFMaxAttribs * cGTFMaxValLen));	// max number of characters which can be buffered

typedef enum TAG_etGTFFeatType {
	eGGTFany = 0,		// any GTF record
	eGGTFCDS,			// coding sequence
	eGGTFStartCodon,	// start codon
	eGGTFStopCodon,		// stop codon
	eGGTF5UTR,			// 5' UTR
	eGGTF3UTR,			// 3' UTR
	eGGTFinter,			// intergenic region
	eGGTFinter_CNS,		// intergenic conserved noncoding sequence region
	eGGTFintron_CNS,	// conserved noncoding sequence region within an intron
	eGGTFexon,			// exon
	eGGTFplaceholder
} etGGTFFeatType;

#pragma pack(1)
typedef struct TAG_sGTFFields {
	char szRawLine[cGTFMaxLineLen];					// raw record line
	char szSeqName[cMaxDatasetSpeciesChrom];		// sequence name
	char szSeqSource[cMaxDatasetSpeciesChrom];		// sequence source
	char szFeature[cMaxDatasetSpeciesChrom];		// feature type
	UINT32 Start;									// feature start locus (1..End)
	UINT32 End;										// feature end locus (inclusive and >= start)
	int ScoreOfs;									// offset into szRawLine at score starts
	bool bDfltScore;								// true if no score specified and default (0.0) is used
	float Score;									// if none then 0.0
	bool bDfltStrand;								// true if no strand specified and default '?' is used
	char Strand;									// '+', '-' or '?'
	bool bDfltFrame;								// true if no frame specified and default (0) is used
	int Frame;										// 0..2
	int AttribOfs;									// offset into szRawLine at which attributes start
	} tsGTFFields;

#pragma pack()


class CGTFFile  : public CErrorCodes
{
	FILE *m_pGTFStream;		// GTF file opened as a stream
	char m_szInFile[_MAX_PATH+1];	// input gff file

	bool m_bEOF;			// set TRUE when last record has been parsed

	int m_CurLineNum;		// current line number
	int m_CurLineLen;		// current line length
	int m_CurNumFields;		// number of fields parsed from current line
	int m_CurMaxFieldLen;	// maximum field length of any field in current line
	tsGTFFields m_Fields;	// parsed fields

	void Reset(void);
	int ParseLine(void);	// read next line and parse
	char *TrimWhitespace(char *pTxt); // inplace trim leading/trailing whitespace

public:
	CGTFFile(void);
	~CGTFFile(void);
	int Open(char *pszFileName);	// file to open
	int Close(void);

	int GetCurFields(void);				// get current number of attributes
	int GetCurFieldLen(void);			// get current maximum length of any attribute value
	int GetCurLineLen(void);			// get current line length

	int NextRecordOfType(etGGTFFeatType GGTFFeatType); // moves to next record of specified type

	int NextLine(void);					// move to next line in GTF file, skips blank lines and comment lines starting with comment char
	int GetLineNumber(void);			// get current line number
	char *GetRecord(void);				// returns complete current record
	tsGTFFields *GetFields(void);		// returns all parsed out fields for current record
	char *GetSeqName(void);
	char *GetSource(void);
	char *GetFeature(void);
	int GetStart(void);
	int GetEnd(void);
	float GetScore(void);
	char GetStrand(void);
	char GetFrame(void);
	char *GetNoteValue(char *pszAttribute);

};


